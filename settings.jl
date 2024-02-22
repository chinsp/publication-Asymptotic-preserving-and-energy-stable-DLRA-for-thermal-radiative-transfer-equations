__precompile__

mutable struct Settings
        ## Settings of the staggered grids
    # Number of spatial vertices
    Nx::Int64;
    # Number of cell centres
    NxC::Int64;
    # Start and end point of spatial domain
    a::Float64;
    b::Float64;
    # Grid cell width
    dx::Float64;
    
    ## Settings of temporal domain
    # End time
    Tend::Float64;
    # Time increment width
    dt::Float64;
    # CFL number 
    cfl1::Float64; # For parabolic CFL condition
    cfl2::Float64; # For hyperbolic CFL condition
    # CFL condition type
    cflType::String;

    ## Settings for angular approximation
    # Number of quadrature points
    N::Int64;
    
    ## Spatial grid
    x
    xMid

    ## Problem Settings
    problem::String;

    ## Rescaled collision length
    epsilon::Float64;    
    epsi_array_Nx::Array;
    epsi_array_NxC::Array;

    ## Initial conditions
    ICType::String;
    BCType::String;

    ## Physical parameters
    sigmaS::Float64; ## Try to change this to get a non-constant value for the scattering coefficients
    sigmaA::Float64;

    ## Physical parameters for the Marshak Wave
    c::Float64   # Speed of light
    StefBoltz::Float64   # Stefan Boltzmann constant
    aRad::Float64   # Radition constant
    c_nu::Float64   # Specific heat

    ## Dynamcial low-rank approximation
    r::Int; # rank of the system
    epsAdapt::Float64; # Tolerance for rank adaptive BUG integrator

    # Linearisation of the problem
    LinTyp::String;

    # If the problem impliments and absorber
    alim::Float64; # The position of the absorber in the domain

    function Settings(Nx::Int=1001,N::Int=500,epsilon::Float64=1.0,cflType::String="hyperbolic")
        # Setup spatial grid
        NxC = Nx + 1;
        a = -10 #0.0; # Starting point for the spatial interval
        b = 10 #0.002; # End point for the spatial interval

        x = collect(range(a,stop = b,length = Nx));
        dx = x[2] - x[1];
        xMid = [x[1]-dx;x];
        xMid = xMid .+ dx/2
        # xMid = x .+ dx/2;
        # xMid = xMid[1:(end-1)];

        # println("Number of points for spatial discretisation of macro and meso variables = ",Nx);
        # println("Number of points for spatial discretisation of micro variable = ",NxC);

        # Problem 
        problem = "1DLinearTestcase" #  1DLinearTestcase, 1DAbsorberTestcase

        # Scattering and absorption coefficients
        sigmaA = 1.0 / 0.926 / 1e-6 / 100.0;
        sigmaS = 0.0;

        # Initial and Boundary condition
        ICType = "LS";
        BCType = "exact";

        # Linearisation of the problem
        LinTyp = "NL"; # The options are NL for non-linear  and Lin for linear

        # Defining the constants related to the simulation
        c = 299792458.0 * 100.0;                # speed of light in [cm/s]
        StefBoltz = 5.6704 * 1e-5;                # Stefan Boltzmann constant in [erg/cm^2/s/K^4]
        aRad = 4 * StefBoltz / c;     # Radiation constant
        density = 0.01;
        c_nu = density * 0.831 * 1e7;    # heat capacity: [kJ/(kg K)] = [1000 m^2 / s^2 / K] therefore density * 0.831 * 1e7
        # c_nu = 0.831 * 1e9;
        # Setup temporal discretisation
        Tend = 0.05 * 1e-10;
        cfl1 = 1.0; # CFL condition parabolic
        cfl2 = 1.0; # CFL condition hyperbolic
        alim = 0.0;
        epsi_array_Nx = zeros(Nx);
        epsi_array_NxC = zeros(NxC);
        if problem == "1DLinearTestcase"
            Tend = 1.5;
            sigmaA = 0.5;
            aRad = 1.0;
            c = 1.0;
            c_nu = 1.0;
            LinTyp = "Lin";
            alim = 0.25;
        elseif problem == "1DAbsorberTestcase"
            Tend = 1.5;
            sigmaA = 0.5;
            aRad = 1.0;
            c = 1.0;
            c_nu = 1.0;
            LinTyp = "Lin";
            alim = 0.25;
        end

        if cflType == "parabolic"
            dt = cfl1*(sigmaA*dx^2);
        elseif cflType == "hyperbolic"
            dt = cfl2*2*epsilon*dx;
        elseif cflType == "mixed"
            dt = cfl1*sigmaA*dx^2 + cfl2*2*epsilon*dx; # The step size for mixed epsilon is not set here
        else
            println("Please enter valid type for CFL condition")
        end

        # Settings for BUG integrator
        r = 5;
        epsAdapt = 0.05; # Tolerance for rank adaptive integrator

        new(Nx,NxC,a,b,dx,Tend,dt,cfl1,cfl2,cflType,N,x,xMid,problem,epsilon,epsi_array_Nx,epsi_array_NxC,ICType,BCType,sigmaS,sigmaA,c,StefBoltz,aRad,c_nu,r,epsAdapt,LinTyp,alim);
    end 
end

function IC(obj::Settings)
    f0 = zeros(obj.NxC,obj.N);
    f0_mu = zeros(obj.Nx);
    T0 = zeros(obj.Nx);
    x = obj.x;
    if obj.problem == "1DLinearTestcase"
        T0 = zeros(obj.Nx);
        for j = 1:obj.Nx
            if obj.x[j] >= -0.5 && obj.x[j] <= 0.5
                T0[j] = 100.0/obj.aRad/obj.sigmaA;
            end
        end
        # T0[1] =  50.0;
        h0 = f0_mu ; #-1/epsilon^2 .* obj.aRad * obj.c * T0;
        # h0 = f0_mu
        g0 = f0;
    elseif obj.problem == "1DAbsorberTestcase"
        T0 = zeros(obj.Nx);
        for j = 1:obj.Nx
            if obj.x[j] >= -0.5 && obj.x[j] <= 0.5
                T0[j] = 100.0/obj.aRad/obj.sigmaA;
                if obj.x[j] >= -obj.alim && obj.x[j] <= obj.alim
                    T0[j] = 100.0/obj.aRad/5.0;
                end
            end
        end
        h0 = f0_mu ; #-1/epsilon^2 .* obj.aRad * obj.c * T0;
        # h0 = f0_mu
        g0 = f0;
    else
        println("Initial condition not coded yet")
        
    end
    return h0,g0,T0;
end