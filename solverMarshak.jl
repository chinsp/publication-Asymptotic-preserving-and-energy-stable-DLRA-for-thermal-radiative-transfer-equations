__precompile__
include("quadrature.jl")

using ProgressMeter
using ProgressBars
using LinearAlgebra
using FastGaussQuadrature, LegendrePolynomials
using LinearSolve
using LowRankApprox
using PyCall

struct solverMarshak
    # Spatial grid of cell vertices
    x::Array{Float64,1};
    xMid::Array{Float64}

    # Solver settings
    settings::Settings;

    # Quadrature for flux matrix
    w::Array{Float64,1};
    # v::Array{Float64,2};

    # Pn discretisation
    AFull::Array{Float64,2};
    A::Array{Float64,2};
    absA::Array{Float64};
    Ap::Array{Float64,2};
    Am::Array{Float64,2};
    gamma::Array{Float64,1};

    # Stencil matrices for spatial discretisation
    Dx::Array{Float64,2};
    Dxx::Array{Float64,2};
    Dp::Array{Float64,2};
    Dm::Array{Float64,2};
    Do::Array{Float64,2};
    deltao::Array{Float64,2};
    b::Array{Float64,1};

    # Physical parameters
    sigmaA::Float64;
    sigmaS::Float64;

    # Temporal step size constant
    Const::Float64;
    
    # Constructor
    function solverMarshak(settings)
        x = settings.x;
        xMid = settings.xMid;

        Nx = settings.Nx;
        NxC = settings.NxC;

        # Setting up the flux matrix
        N = settings.N;
        gamma = zeros(Float64,N+1);

        for i = 1:N+1
            gamma[i] = 2/(2*(i-1) + 1);
        end

        a_norm = zeros(Float64,N+1);
        for i = 1:N+1
            a_norm[i] = i/sqrt((2*(i-1)+1)*(2*(i-1)+3));
        end

        AFull = zeros(Float64,N+1,N+1);
        AFull = Tridiagonal(a_norm[1:end-1],zeros(N+1),a_norm[1:end-1]);

        # A = T*diagm(mu)*transpose(T);
        A = AFull[2:end,2:end];

        TFull = zeros(N+1,N+1);
        mu,w = gausslegendre(N+1);
        for k = 1:N+1
            P = collectPl(mu[k],lmax = N);
            for i = 1:N+1
                TFull[i,k] = P[i-1]*sqrt(w[k])/sqrt(gamma[i]);
            end
        end
        T = TFull[2:end,:];
        absA = T*abs.(diagm(mu))*transpose(T);
        Ap = 0.5 .* (A + absA);
        Am = 0.5 .* (A - absA);

        # Setup temporal discretisation
        BetaN = maximum(w) * (N+1);
        mu1 = deleteat!(mu,findall(x->x==0,mu));
        dt_array = zeros(Float64,length(mu1));

        if settings.cflType == "parabolic"
            for k = eachindex(mu1)
                dt_array[k] = settings.dt /1/5/BetaN/settings.c/ mu1[k]^2;
            end
            _,k_min = findmin(dt_array);
            musq_max = mu1[k_min]^2;
        elseif settings.cflType == "hyperbolic"
            for k = eachindex(mu1)
                dt_array[k] = settings.dt /1/5/BetaN/settings.c/ abs(mu1[k]);
            end
            _,k_min = findmin(dt_array);
            musq_max = abs(mu1[k_min]);
        elseif settings.cflType == "mixed"
            for k = eachindex(mu1)
                dt_array[k] = (settings.cfl1*settings.sigmaA*settings.dx^2/ mu1[k]^2 + settings.cfl2*2*settings.epsilon*settings.dx/ abs(mu1[k]));
            end
            dt_min,k_min = findmin(dt_array);
            settings.dt = dt_min;
            musq_max = 1.0;
        end
        # println("Max. weight: ",w[k_min]," , Max. quadrature point: ",mu1[k_min])
        
        Const = 1/5/BetaN/settings.c/musq_max;

        # Stencil matrices for spatial discretisation
        dx = settings.dx;

        Dp = zeros(NxC,NxC);
        Dm = zeros(NxC,NxC);
        
        for i = 1:NxC-1
            Dp[i,i] = -1/dx;
            Dp[i,i+1] = 1/dx;
        end

        for i = 2:NxC
            Dm[i,i] = 1/dx;
            Dm[i-1,i] = -1/dx;
        end
        # set up spatial stencil matrices
        Dx = Tridiagonal(-ones(NxC-1)./dx/2.0,zeros(NxC),ones(NxC-1)./dx/2.0) # central difference matrix
        Dxx = Tridiagonal(ones(NxC-1)./dx/2.0,-ones(NxC)./dx,ones(NxC-1)./dx/2.0) # stabilization stencil matrix

        # Stencil matrices for limiting diffusion equation
        deltao = zeros(NxC, Nx);
        Do = zeros(Nx, NxC);
        for i = 2:Nx
            deltao[i,i] = 1/dx;
            deltao[i,i-1] = -1/dx;
        end
        deltao[1,1], deltao[end,end] = 1/dx, -1/dx;

        for i = 1:Nx
            Do[i,i] = -1/dx;
            Do[i,i+1] = 1/dx;
        end
        # Do[1,:],Do[end,:] = zeros(NxC),zeros(NxC);

        b = zeros(N);
        b[1] = sqrt(gamma[2]);
        
        new(x,xMid,settings,w,AFull,A,absA,Ap,Am,gamma,Dx,Dxx,Dp,Dm,Do,deltao,b,settings.sigmaA,settings.sigmaS,Const)
    end
end

# py"""
# import numpy
# def qr(A):
#     return numpy.linalg.qr(A)
# """

function setupIC(obj::solverMarshak)
    h0 = zeros(obj.settings.Nx);
    g0 = zeros(obj.settings.NxC,obj.settings.N);
    T0 = zeros(obj.settings.Nx);
    h0,g0,T0 = IC(obj.settings);
    return h0,g0,T0;
end

function BCT(obj::solverMarshak,T::Array)
    if obj.settings.problem == "1DLinearTestcase"
        T[1],T[end] = 0.0,0.0;
    end
    return T;
end

function BCg(obj::solverMarshak,g::Array)
    if obj.settings.problem == "1DLinearTestcase"
        g[1,:] = zeros(obj.settings.N);
        g[end,:] = zeros(obj.settings.N);
    end
    return g;
end

function BCh(obj::solverMarshak,h::Array,T::Array)
    aRad = obj.settings.aRad;
    c = obj.settings.c;
    epsilon = obj.settings.epsilon;
    if obj.settings.problem == "1DLinearTestcase"
        h[1],h[end] = -1/epsilon^2 * aRad * c * T[1],-1/epsilon^2 * aRad * c * T[end];
    elseif obj.settings.problem == "1DAbsorberTestcase"
        # h[1],h[end] = -1/epsilon^2 * aRad * c * T[1],-1/epsilon^2 * aRad * c * T[end];
    end
    return h;
end

function ComputeEnergy(obj::solverMarshak,h::Array,g::Array, T::Array)
    aRad = obj.settings.aRad;
    epsilon = obj.settings.epsilon;
    # epsi_array_Nx = obj.settings.epsi_array_Nx;
    # epsi_array_NxC = obj.settings.epsi_array_NxC;
    c = obj.settings.c;
    gamma0 = sqrt(2);
    c_nu = obj.settings.c_nu;
    dx = obj.settings.dx;
    energy = 0.0;
    e1,e2,e3 = 0.0,0.0,0.0;
    Nx, NxC = obj.settings.Nx, obj.settings.NxC;
    for i = 1:Nx
        e1 = e1 + (aRad * T[i] + epsilon^2 / c * h[i])^2 * dx;
        e2 = e2 + (sqrt(aRad * c_nu / 2) * T[i])^2 * dx;
    end
    for i = 1:NxC
        e3 = e3 + (epsilon/gamma0/c)^2 * transpose(g[i,:]) * g[i,:] * dx;
    end
    energy = e1 + e2 + e3;
    return energy;
end

function DTPlanck(obj::solverMarshak, T::Array,type::String)
    BT = zeros(size(T));
    if type == "NL"
        BT = 4 .* T .^3;
    elseif type == "Lin"
        BT = ones(size(BT));
    else
        println("Give valid input for the derivative of the frequency integrated Planckian.");
    end
    return BT;
end

function AbsorpMatrix(obj::solverMarshak)
    INx = I(obj.settings.Nx);
    INxC = I(obj.settings.NxC);
    if obj.settings.problem == "1DAbsorberTestcase"
        SigmaA = obj.sigmaA*diagm(ones(obj.settings.NxC));
        SigmaAf = obj.sigmaA*diagm(ones(obj.settings.Nx));
        alim = obj.settings.alim;
        ## Implementation for an absorber in the middle of the domain
        for j = 1:obj.settings.NxC
            if obj.xMid[j] >= -alim && obj.xMid[j] <= alim
                SigmaA[j,j] = 5.0;
            end
        end
        for j = 1:obj.settings.Nx
            if obj.x[j] >= -alim && obj.x[j] <= alim
                SigmaAf[j,j] = 5.0;
            end
        end
    else
        SigmaA = obj.sigmaA*INxC;
        SigmaAf = obj.sigmaA*INx;
    end
    return SigmaA, SigmaAf;
end

function PrintSolverInformation(obj::solverMarshak,solverName::String, dt::Float64)
    println("=======================================================");
    println("Solver: ", solverName);
    println("Problem: ", obj.settings.problem);
    println("-------------------------------------------------------");
    println("Discretisation parameters");
    println("Spatial resolution: Δx = ", obj.settings.dx);
    println("Mean free path: ϵ = ", obj.settings.epsilon);
    println("Step size: Δt = ", dt);
    println("-------------------------------------------------------");
    println("Settings for Dynamical low-rank integrator")
    println("Rank: r = ", obj.settings.r);
    println("Adaptive BUG tolerance: ϑ = ", obj.settings.epsAdapt);
    println("=======================================================");
end

function solveRosselandLimit(obj::solverMarshak)
    t = 0.0;
    dt = obj.Const *obj.settings.dt;
    Tend = obj.settings.Tend;
    Nx = obj.settings.Nx;
    NxC = obj.settings.NxC;
    c_nu = obj.settings.c_nu;
    aRad = obj.settings.aRad;
    c = obj.settings.c; # Speed of particles

    _, _, T = setupIC(obj);

    ## Initialising stencil matrices
    deltao = obj.deltao;
    Do = obj.Do;
    
    # Linearisation of the problem
    LinTyp = obj.settings.LinTyp;
    
    Nt = round(Tend/dt);
    
    psi_f = diagm(DTPlanck(obj,T,LinTyp)); # For the non-linear system
    psi_aux = zeros(Nx-1);
    alpha1, alpha2 = 0.5,0.5;
    for i = 1:Nx-1
        psi_aux[i] = alpha1 *psi_f[i,i] + alpha2 * psi_f[i+1,i+1];
    end
    psi_avg = diagm([psi_f[1,1];psi_aux;psi_f[end,end]]);
    
    SigmaA,_ = AbsorpMatrix(obj);
    SigmaAinv = inv(SigmaA);

    PrintSolverInformation(obj,"Explicit Rosseland diffusion limit",dt)
    println("Running explicit solver for limiting Rosseland diffusion equation")

    for k = ProgressBar(1:Nt)
        # Update temperature
        M = inv(c_nu*I(Nx) + 2*aRad*psi_f);
        T .= T .+ (2*aRad*c*dt/3) .* M * Do * SigmaAinv * psi_avg * deltao * T;
        T .= BCT(obj,T);
        
        ## Update psi_f and psi_avg
        psi_f = diagm(DTPlanck(obj,T,LinTyp)); # For the non-linear system
        psi_aux = zeros(Nx-1);
        for i = 1:Nx-1
            psi_aux[i] = alpha1 *psi_f[i,i] + alpha2 * psi_f[i+1,i+1];
        end
        psi_avg = diagm([psi_f[1,1];psi_aux;psi_f[end,end]]);
        
        t = t + dt;
    end
    return t, T
end

py"""
import numpy
def qr(A):
    return numpy.linalg.qr(A)
"""

function solveFullMacroMicro(obj::solverMarshak)
    t = 0.0;
    dt = obj.Const * obj.settings.dt;
    dx = obj.settings.dx;
    Tend = obj.settings.Tend;
    Nx = obj.settings.Nx;
    NxC = obj.settings.NxC;
    epsilon = obj.settings.epsilon;

    ## Physical constants
    c_nu = obj.settings.c_nu;
    aRad = obj.settings.aRad;
    c = obj.settings.c; # Speed of particles
    kappa = 2/c_nu;
    gamma1 = sqrt(obj.gamma[2]);

    ## Spatial stencil matrices
    deltao = obj.deltao;
    Do = obj.Do;
    Dx = obj.Dx;
    Dxx = obj.Dxx;

    ## Flux matrix and stabilisation matrices
    A = obj.A;
    absA = obj.absA;
    b = obj.b;


    #Identity matrices required 
    INx = I(Nx);
    INxC = I(NxC);

    SigmaA,SigmaAf = AbsorpMatrix(obj);


    e1 = zeros(Float64,obj.settings.N);
    e1[1] = 1;

    h, g, T = setupIC(obj);
    mass_0 = sum((aRad*c*T + epsilon^2 .* h) + c_nu*c/2 .* T) * dx;

    # Linearisation of the problem
    LinTyp = obj.settings.LinTyp;

    psi_f = diagm(DTPlanck(obj,T,LinTyp)); # For the non-linear system
    psi_aux = zeros(Nx-1);
    for i = 1:Nx-1
        psi_aux[i] = (psi_f[i,i] + psi_f[i+1,i+1])/2;
    end
    psi_avg = diagm([psi_f[1,1];psi_aux;psi_f[end,end]]);

    ## Initialising matrices used for computation
    Beta = inv(1 / c .* INxC + dt*SigmaA ./epsilon^2);
    alpha_n = diagm(zeros(Nx));
    Nt = Int(round(Tend/dt));
    
    ## Initiate array to save energy
    energy = zeros(Float64, Nt+1);
    energy[1] = ComputeEnergy(obj,h,g,T);

    ## Checking for mass conservation
    mass = zeros(Float64,Nt);

    PrintSolverInformation(obj,"Full modal macro-micro explicit solver",dt)
    println("Running explicit solver for full modal macro-micro solver:")
    for k = ProgressBar(1:Nt)
        ## Micro Update
        g .= 1 /c .* g .- dt .* Dx * g * transpose(A) ./ epsilon .+ dt .* Dxx * g * transpose(absA) ./ epsilon .- aRad * c* dt / epsilon^2 .* psi_avg * deltao * T * transpose(b) .- dt .* deltao * h * transpose(b);
        g .= Beta * g;
        # Boundary condition
        g .= c.*BCg(obj,g);

        ## Meso Update
        alpha_n .= inv(1 / c * INx + dt*aRad*kappa/epsilon^2 * SigmaAf * psi_f + dt*SigmaAf ./epsilon^2);
        h .= 1 / c .* h .- gamma1 * dt  /2 /epsilon^2 .* Do * g * e1;
        h .=  alpha_n * h;
        h .= c.*BCh(obj,h,T);
        
        ## Macro update
        T .= T + dt * kappa .* SigmaAf * h;
        # Boundary conditions
        T .= BCT(obj,T);

        ## Update gradient of Planckian
        # Update gradient
        psi_f .= diagm(DTPlanck(obj,T,LinTyp)); # For the non-linear system
        psi_aux = zeros(Nx-1);
        for i = 1:Nx-1
            psi_aux[i] = (psi_f[i,i] + psi_f[i+1,i+1])/2;
        end
        psi_avg = diagm([psi_f[1,1];psi_aux;psi_f[end,end]]);
        
        ## Compute the energy in the system
        energy[k+1] = ComputeEnergy(obj,h,g,T);

        ## Compute mass of the system 
        mass_n = sum((aRad*c*T + epsilon^2 .* h) + c_nu*c/2 .* T) * dx;
        mass[k] = abs(mass_0 - mass_n)/abs(mass_0);
        # mass_0 = mass_n;
 
        t = t + dt;
    end

    return t, h, g, T, energy, mass;
end


function solveBUGintegrator(obj::solverMarshak)
    t = 0.0;
    dt = obj.Const * obj.settings.dt;
    Tend = obj.settings.Tend;
    Nx = obj.settings.Nx;
    NxC = obj.settings.NxC;
    dx = obj.settings.dx;
    epsilon = obj.settings.epsilon;
    r = obj.settings.r;

    ## Physical constants
    c_nu = obj.settings.c_nu;
    aRad = obj.settings.aRad;
    c = obj.settings.c; # Speed of particles
    kappa = 2/c_nu;
    gamma1 = sqrt(obj.gamma[2]);

    ## Spatial stencil matrices
    deltao = obj.deltao;
    Do = obj.Do;
    Dx = obj.Dx;
    Dxx = obj.Dxx;

    ## Flux matrix and stabilisation matrices
    A = obj.A;
    absA = obj.absA;
    b = obj.b;

    #Identity matrices required 
    INx = I(Nx);
    INxC = I(NxC);

    SigmaA,SigmaAf = AbsorpMatrix(obj);

    e1 = zeros(Float64,obj.settings.N);
    e1[1] = 1;

    h, g, T = setupIC(obj);
    mass_0 = sum((aRad*c*T + epsilon^2 .* h) + c_nu*c/2 .* T) * dx;

    X0,s,V0 = svd(g);
    X = X0[:,1:r];
    V = V0[:,1:r];
    S = diagm(s[1:r]);
    K = zeros(size(X));
    L = zeros(size(transpose(V)));

    # Linearisation of the problem
    LinTyp = obj.settings.LinTyp;

    psi_f = diagm(DTPlanck(obj,T,LinTyp)); # For the non-linear system
    # psi_f = diagm(ones(size(T))); # For the linearised system
    psi_aux = zeros(Nx-1);
    for i = 1:Nx-1
        psi_aux[i] = (psi_f[i,i] + psi_f[i+1,i+1])/2;
    end
    psi_avg = diagm([psi_f[1,1];psi_aux;psi_f[end,end]]);

    ## Initialising matrices used for computation
    Beta_K = 1 / c .* INxC + dt*SigmaA ./epsilon^2;
    alpha_n = diagm(zeros(Nx));
    Nt = Int(round(Tend/dt));

    ## Initiate array to save energy
    energy = zeros(Float64, Nt+1);
    energy[1] = ComputeEnergy(obj,h,g,T);
    
    ## Checking for mass conservation
    mass = zeros(Float64,Nt);
    abc = zeros(Float64,Nt);

    PrintSolverInformation(obj,"BUG solver",dt)
    println("Running modal macro-micro BUG solver:")
    for k = ProgressBar(1:Nt)
        ## Micro Update
        # K-step 
        K = X * S;
        VAV = transpose(V) * transpose(A) * V;
        VabsAV = transpose(V) * transpose(absA) * V;
        K = 1 / c .* K .- dt / epsilon .* Dx * K * VAV .+ dt / epsilon .* Dxx * K * VabsAV .- dt * aRad * c /epsilon^2 .* psi_avg * deltao *  T * transpose(b) * V .- dt .* deltao * h * transpose(b) * V;
        K = Beta_K \ K;
        # BC condition
        # K[1,:],K[end,:] = zeros(Float64,r),zeros(Float64,r);
        X1,_ = qr(K);
        X1 = Matrix(X1);
        X1[1,:],X1[end,:] = zeros(Float64,r),zeros(Float64,r);
        M_BUG = transpose(X1) * X;

        # L-step
        L = S * transpose(V);
        Beta_L = (1 / c .* I(r) + dt/epsilon^2 .* transpose(X) * SigmaA * X);
        
        L = 1 / c .* L .- dt / epsilon .* transpose(X) * Dx * X * L * transpose(A) .+ dt / epsilon .* transpose(X) * Dxx * X * L * transpose(absA) .- dt * aRad * c /epsilon^2 .* transpose(X) * psi_avg * deltao * T * transpose(b) .- dt .* transpose(X) * deltao * h * transpose(b);
        for k = 1:N
            L[:,k] = Beta_L \ L[:,k];
        end
        # L = Beta_L * L;
        V1,_ = qr(transpose(L));
        V1 = Matrix(V1);
        N_BUG = transpose(V1) * V;
        # Update X,V
        X = X1;
        V = V1;

        #S-step
        Beta_S = (1 / c .* I(r) + dt/epsilon^2 .* transpose(X) * SigmaA * X);
        S = M_BUG * S * transpose(N_BUG);
        S = 1 / c .* S .- dt / epsilon .* transpose(X) * Dx * X * S * transpose(V) * transpose(A) * V .+ dt / epsilon .* transpose(X) * Dxx * X * S * transpose(V) * transpose(absA) * V .- dt * aRad * c /epsilon^2 .* transpose(X) * psi_avg * deltao * T * transpose(b) * V .- dt .* transpose(X) * deltao * h * transpose(b) * V;
        for k = 1:r
            S[:,k] = Beta_S \ S[:,k];
        end
        # S = Beta_S * S;
  

        ## Meso Update
        alpha_n = (1 / c * INx + dt*aRad*kappa/epsilon^2 * SigmaAf * psi_f + dt*SigmaAf./epsilon^2);
        h = 1 / c .* h .- gamma1 * dt  /2 /epsilon^2 .* Do * X * S * transpose(V) * e1;
        h =  alpha_n \ h;
        h = BCh(obj,h,T);
        
        ## Macro update
        T = T .+ dt * kappa .* SigmaAf * h;
        # Boundary conditions
        T = BCT(obj,T);

        ## Update gradient of Planckian
        psi_f = diagm(DTPlanck(obj,T,LinTyp)); # For the non-linear system
        psi_aux = zeros(Nx-1);
        for i = 1:Nx-1
            psi_aux[i] = (psi_f[i,i] + psi_f[i+1,i+1])/2;
        end
        psi_avg = diagm([psi_f[1,1];psi_aux;psi_f[end,end]]);

        ## Compute the energy in the system
        energy[k+1] = ComputeEnergy(obj,h,X*S*transpose(V),T);

        ## Compute mass of the system 
        mass_n = sum((aRad*c*T .+ epsilon^2 .* h) .+ c_nu*c/2 .* T) * dx;
        mass[k] = abs(mass_0 - mass_n)/abs(mass_0);
        # mass_0 = mass_n;
        abc[k] = sum(Do * X * S * transpose(V) * e1);
        t = t + dt;
    end
    return t, h, X*S*transpose(V), T, energy, mass,abc;
end


function solveBUG_rankadaptive(obj::solverMarshak)
    t = 0.0;
    dt = obj.Const * obj.settings.dt;
    dx = obj.settings.dx;
    Tend = obj.settings.Tend;
    Nx = obj.settings.Nx;
    NxC = obj.settings.NxC;
    epsilon = obj.settings.epsilon;
    r = obj.settings.r;
    ## Physical constants
    c_nu = obj.settings.c_nu;
    aRad = obj.settings.aRad;
    c = obj.settings.c; # Speed of particles
    kappa = 2/c_nu;
    gamma1 = sqrt(obj.gamma[2]);

    ## Spatial stencil matrices
    deltao = obj.deltao;
    Do = obj.Do;
    Dx = obj.Dx;
    Dxx = obj.Dxx;

    ## Flux matrix and stabilisation matrices
    A = obj.A;
    absA = obj.absA;
    b = obj.b;

    #Identity matrices required 
    INx = I(Nx);
    INxC = I(NxC);

    SigmaA,SigmaAf = AbsorpMatrix(obj);

    e1 = zeros(Float64,obj.settings.N);
    e1[1] = 1;

    h, g, T = setupIC(obj);
    mass_0 = sum((aRad*c*T + epsilon^2 .* h) + c_nu*c/2 .* T) * dx;

    X0,s,V0 = svd(g);
    X0 = Matrix(X0);
    X = X0[:,1:r];
    V1 = Matrix(V0);
    V = V1[:,1:r];
    S = diagm(s[1:r]);
    K = zeros(size(X));
    L = zeros(size(transpose(V)));

    #  for the non-linear case psi_f and psi_avg need to be recomputed after every iteration
    # Linearisation of the problem
    LinTyp = obj.settings.LinTyp;

    psi_f = diagm(DTPlanck(obj,T,LinTyp)); # For the non-linear system
    psi_aux = zeros(Nx-1);
    for i = 1:Nx-1
        psi_aux[i] = (psi_f[i,i] + psi_f[i+1,i+1])/2;
    end
    psi_avg = diagm([psi_f[1,1];psi_aux;psi_f[end,end]]);


    ## Initialising matrices used for computation
    Beta_K = inv(1 / c .* INxC + dt*SigmaA ./epsilon^2);
    alpha_n = diagm(zeros(Nx));
    Nt = Int(round(Tend/dt));

    ## Initiate array to save energy and ranks
    energy = zeros(Float64, Nt+1);
    energy[1] = ComputeEnergy(obj,h,g,T);
    ranks = zeros(Int64,Nt+1);
    ranks[1] = r;

    
    ## Checking for mass conservation
    mass = zeros(Float64,Nt);

    PrintSolverInformation(obj,"rank-adaptive BUG solver",dt)
    println("Running asymptotic-preserving rank-adaptive BUG solver for modal macro-micro equations:")
    for k = ProgressBar(1:Nt)
        ## Micro Update
        # K-step 
        K = X * S;
        K .= 1 / c .* K .- dt / epsilon .* Dx * K * transpose(V) * transpose(A) * V .+ dt / epsilon .* Dxx * K * transpose(V) * transpose(absA) * V - dt * aRad * c /epsilon^2 .* psi_avg * deltao *  T * transpose(b) * V - dt .* deltao * h * transpose(b) * V;
        K .= Beta_K * K;
        # Boundary condition
        # K[1,:],K[end,:] = zeros(Float64,r),zeros(Float64,r);
        X1,STmp = qr([aRad*c.*inv(SigmaA)*psi_avg*deltao*T K X]);
        X1 = Matrix(X1);
        X1[1,:], X1[end,:] = zeros(Float64,2*r+1),zeros(Float64,2*r+1);
        M_BUG = transpose(X1) * X;

        # L-step
        L = S * transpose(V);
        Beta_L = inv(1 / c .* I(r) + dt/epsilon^2 .* transpose(X) * SigmaA * X);
        L .= 1 / c .* L .- dt / epsilon .* transpose(X) * Dx * X * L * transpose(A) .+ dt / epsilon .* transpose(X) * Dxx * X * L * transpose(absA) - dt * aRad * c /epsilon^2 .* transpose(X) * psi_avg * deltao * T * transpose(b) - dt .* transpose(X) * deltao * h * transpose(b);
        L .= Beta_L * L;
        V1,STmp = qr([b transpose(L) V]);
        V1 = Matrix(V1);
        N_BUG = transpose(V1) * V;

        #S-step
        Beta_S = inv(1 / c .* I(2*r + 1) + dt/epsilon^2 .* transpose(X1) * SigmaA * X1);
        S = M_BUG * S * transpose(N_BUG);
        S1 = 1 / c .* S .- dt / epsilon .* transpose(X1) * Dx * X1 * S * transpose(V1) * transpose(A) * V1 .+ dt / epsilon .* transpose(X1) * Dxx * X1 * S * transpose(V1) * transpose(absA) * V1 .- dt * aRad * c /epsilon^2 .* transpose(X1) * psi_avg * deltao * T * transpose(b) * V1 .- dt .* transpose(X1) * deltao * h * transpose(b) * V1
        S1 .= Beta_S * S1;

        #Conservatie truncation
        Khat = X1 * S1;
        m = 1; # Number of basis vectors to left unchanged 
        Khat_ap, Khat_rem = Khat[:,1:m], Khat[:,m+1:end]; # Splitting Khat into basis required for Ap and remaining vectors
        Vap, Vrem = V1[:,1:m], V1[:,m+1:end]; # Splitting Khat into basis required for Ap and remaining vectors

        Xhrem, Shrem = qr(Khat_rem);
        Xhrem = Matrix(Xhrem);

        U, sigma, W = svd(Shrem);
        U = Matrix(U);
        W = Matrix(W);

        rmax = -1
        tmp = 0.0;
        tol = obj.settings.epsAdapt * norm(sigma);

        # Truncating the rank
        for i = 1:2*r+1-m
            tmp = sqrt(sum(sigma[i:2r+1-m].^2));
            if tmp < tol
                rmax = i
                # println(rmax,tmp,tol);
                break
            end
        end
        rmaxTotal = Int(round(max(obj.settings.NxC, obj.settings.N)/2));

        rmax = min(rmax,rmaxTotal);
        rmax = max(rmax,1);


        if rmax == -1
            rmax = rmaxTotal;
        end

        Uhat = U[:,1:rmax];
        What = W[:,1:rmax];
        sigma_hat = diagm(sigma[1:rmax]);

        W1 = Vrem * What;
        Xrem = Xhrem * Uhat;

        V = [Vap W1];
        Xap, Sap = qr(Khat_ap);
        Xap = Matrix(Xap);
        X, R2 = qr([Xap Xrem]);
        X = Matrix(X);
        X[1,:], X[end,:] = zeros(Float64,rmax+m),zeros(Float64,rmax+m);
        S = zeros(rmax+m,rmax+m);
        S[1:m,1:m] = Sap;
        S[m+1:end,m+1:end] = sigma_hat;
        S .= R2 * S;
        r = rmax+m;
        ranks[k+1] = r;

        ## Meso Update
        alpha_n = (1 / c * INx + dt*aRad*kappa/epsilon^2 * SigmaAf * psi_f + dt*SigmaAf ./epsilon^2);
        h .= 1 / c .* h .- gamma1 * dt  /2 /epsilon^2 .* Do * X * S * transpose(V) * e1;
        h .=  inv(alpha_n) * h;
        h .= BCh(obj,h,T);
        
        ## Macro update
        T .= T + dt * kappa .* SigmaAf * h;
        # Boundary conditions
        T .= BCT(obj,T);

        ## Update gradient of Planckian
        # Update gradient
        psi_f .= diagm(DTPlanck(obj,T,LinTyp)); # For linear problem set to "L"
        psi_aux = zeros(Nx-1);
        for i = 1:Nx-1
            psi_aux[i] = (psi_f[i,i] + psi_f[i+1,i+1])/2;
        end
        psi_avg = diagm([psi_f[1,1];psi_aux;psi_f[end,end]]);

        ## Compute the energy in the system
        energy[k+1] = ComputeEnergy(obj,h,X*S*transpose(V),T);

        ## Compute mass of the system 
        mass_n = sum((aRad*c*T + epsilon^2 .* h) + c_nu*c/2 .* T) * dx;
        mass[k] = abs(mass_0 - mass_n)/abs(mass_0);
        # mass_0 = mass_n;
    
        t = t + dt;
    end
    return t, h, X*S*transpose(V), T, energy, ranks, mass;
end

