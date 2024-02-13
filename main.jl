include("settings.jl")
include("solverMarshak.jl")

using Base.Threads
using PyPlot
using DelimitedFiles
using BenchmarkTools
using LaTeXStrings

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 30

close("all")
epsilon = 1.0;

# Nx = 501;
# N = 5;
# s = Settings(Nx,N,epsilon,"hyperbolic");
# solver = solverMarshak(s);
# t_ex5, h_ex5, g_ex5, T_ex5, energy_ex5 = solveFullMacroMicroExplicit(solver);


# Nx = 501;
# N = 15;
# s = Settings(Nx,N,epsilon,"hyperbolic");
# solver = solverMarshak(s);
# t_ex15, h_ex15, g_ex15, T_ex15, energy_ex15 = solveFullMacroMicroExplicit(solver);

Nx = 501;
N = 100;
s = Settings(Nx,N,epsilon,"hyperbolic");
s = Settings(Nx,N,epsilon,"hyperbolic");
solver = solverMarshak(s);

# t, T = solveRosselandLimitExplicit(solver);

# t_ex, h_ex, g_ex, T_ex, energy_ex = solveFullMacroMicro_MS(solver);
t_ex, h_ex, g_ex, T_ex, energy_ex, mass_ex = solveFullMacroMicroExplicit(solver);

t_exBUG, _, _, T_exBUG, _, _ = solveBUGintegratorExplicit(solver);

# fig1,ax1 = subplots(figsize=(15,12), dpi = 200);
# fig2,ax2 = subplots(figsize=(15,12), dpi = 200);
# rank = [1,2,4,5,8,10,15];
# error = zeros(length(rank));
# i = 1;
# for r in rank
#     s.r = r;
#     r_BUG = s.r;
#     # t_exBUG, h_exBUG, g_exBUG, T_exBUG, energy_exBUG = solveBUGintegrator_MS(solver);
#     t_exBUG, _, _, T_exBUG, _, _ = solveBUGintegratorExplicit(solver);
#     ax1.plot(solver.x,T_exBUG,"--",label = string("r = ", r_BUG));
#     error[i] = LinearAlgebra.norm(T_ex-T_exBUG)/LinearAlgebra.norm(T_ex);
#     global i = i+1;
# end
# ax1.plot(solver.x,T_ex,"--",label = "Full");
# ax1.set_xlim([0,3]);
# ax1.set_xticks(range(0,3, step=0.5));
# ax1.legend()


# ax2.plot(rank,error)

# s.r = 15;
# r_BUGH = s.r;
# t_exBUGH, h_exBUGH, g_exBUGH, T_exBUGH, energy_exBUGH = solveBUGintegratorExplicit(solver);

# s.r = 1;
# t_exraBUG, h_exraBUG, g_exraBUG, T_exraBUG, energy_exraBUG, ranks_exraBUG, mass_exraBUG = solveExplicitBUG_rankadaptive(solver);


# x1 = collect(range(s.a,s.b,size(T)[1]));

fig1, ax1 = subplots(figsize=(15, 12), dpi=200);
# ax1.plot(x1, T,color = "black", "--", label = "Rosseland limit");
ax1.plot(solver.x, T_ex ,"--",color = "red", label = string(L"P_{100}"));
ax1.plot(solver.x, T_exBUG ,"--",color = "green", label = string(L"P_{5}"));
# ax1.plot(solver.x, T_ex15 ,"--",color = "brown", label = string(L"P_{15}"));

# ax1.plot(-s.alim.*ones(100),collect(range(-10,1000,100)),"--",color = "black");
# ax1.plot(s.alim.*ones(100),collect(range(-10,1000,100)),"--",color = "black");
# ax1.set_ylim([minimum(T_ex)-1.0,maximum(T_ex)+5.0]);
# ax1.set_xlim([0,3]);
# ax1.set_xticks(range(0,3, step=0.5));
# ax1.set_ylabel("Temperature");
# ax1.set_xlabel(L"x");
# ax1.legend();
# fig1.canvas.draw();
# fig1
# savefig("Temperature_moment.pdf")

fig1, ax1 = subplots(figsize=(15, 12), dpi=200);
# ax1.plot(x1, T,color = "black", "--", label = "Rosseland limit");
ax1.plot(solver.x, T_ex ,"--",color = "red", label = string(L"P_{100}"));
ax1.plot(solver.x, T_exBUG ,"--",color = "orange", label = string(L"BUG_{5}"));
# ax1.plot(solver.x, T_exBUGH,"--",color = "blue", label = string(L"BUG_{15}"));

# ax1.plot(solver.x, T_exraBUG ,"--",color = "purple", label = string("raBUG"));
# ax1.plot(solver.x, T_exraBUG ,"--",color = "purple", label = string("raBUG"));

# ax1.plot(-s.alim.*ones(100),collect(range(-10,1000,100)),"--",color = "black");
# ax1.plot(s.alim.*ones(100),collect(range(-10,1000,100)),"--",color = "black");
# ax1.set_ylim([minimum(T_ex)-1.0,maximum(T_ex)+5.0]);
ax1.set_xlim([0,3]);
ax1.set_xticks(range(0,3, step=0.5));
ax1.set_ylabel("Temperature");
ax1.set_xlabel(L"x");
ax1.legend();
fig1.canvas.draw();
fig1
# savefig("Temperature_moment.pdf")

# h0,g0,T0 = IC(s);

# fig1, ax1 = subplots(figsize=(15, 12), dpi=200);
# # ax1.plot(x1, T,color = "black", "--", label = "Rosseland limit");
# # ax1.plot(solver.x, T0 ,"--",color = "gray", label = string(L"t = 0"));
# ax1.plot(solver.x, T_ex ,"--",color = "red", label = string(L"P_{100}"));

# ax1.plot(solver.x, T_exBUG ,"--",color = "orange", label = string(L"BUG_{5}"));
# # # ax1.plot(solver.x, T_exBUGH,"--",color = "blue", label = string(L"BUG_{15}"));

# ax1.plot(solver.x, T_exraBUG ,"--",color = "purple", label = string("raBUG"));

# # ax1.plot(-1.0.*ones(100),collect(range(-10,1000,100)),"--",color = "black");
# # ax1.plot(0.5.*ones(100),collect(range(-10,1000,100)),"--",color = "black");
# # ax1.set_ylim([minimum(T0)-1.0,maximum(T0)+10.0]);
# ax1.set_xlim([0,3]);
# ax1.set_xticks(range(0,3, step=0.5));
# ax1.set_ylabel("Temperature");
# ax1.set_xlabel(L"x");
# # ax1.set_title(string(L"\sigma_{a} = ", s.sigmaA))
# ax1.legend();
# fig1.canvas.draw();
# fig1
# savefig("Temperature_lowrank.pdf")

# fig2, ax2 = subplots(figsize=(15, 12), dpi=200);
# ax2.plot(solver.x, (s.aRad * s.c .* T_ex + epsilon^2 .* h_ex),"--",color = "red", label = L"P_{100}");
# ax2.plot(solver.x, (s.aRad * s.c .* T_ex5 + epsilon^2 .* h_ex5),"--",color = "green", label = L"P_{5}");
# ax2.plot(solver.x, (s.aRad * s.c .* T_ex15 + epsilon^2 .* h_ex15),"--",color = "brown", label = L"P_{15}");
# ax2.plot(-s.alim.*ones(100),collect(range(-100,200,100)),"--",color = "black");
# ax2.plot(s.alim.*ones(100),collect(range(-100,200,100)),"--",color = "black");
# ax2.set_ylim([minimum(s.aRad * s.c .* T_ex + epsilon^2 .* h_ex)-1.0,maximum(s.aRad * s.c .* T_ex + epsilon^2 .* h_ex)+5.0]);
# ax2.set_xlim([0,3]);
# ax2.set_xticks(range(0,3, step=0.5));
# ax2.set_ylabel(L"Scalar flux, $\phi(t,x) = \frac{1}{2}\int_{-1}^{1}f(t,x,\mu) \mathrm{d}\mu$");
# ax2.set_xlabel(L"x");
# ax2.legend();
# fig2.canvas.draw();
# fig2
# savefig("Scalar_flux_moment.pdf")

# fig2, ax2 = subplots(figsize=(15, 12), dpi=200);
# ax2.plot(solver.x, s.aRad * s.c .* T0 + epsilon^2 .* h0 ,"--",color = "gray", label = string(L"t = 0"));
# ax2.plot(solver.x, (s.aRad * s.c .* T_ex + epsilon^2 .* h_ex),"--",color = "red", label = L"P_{100}");
# ax2.plot(solver.x, (s.aRad * s.c .* T_exBUG + epsilon^2 .* h_exBUG),"--",color = "orange", label = L"BUG_{5}");
# # ax2.plot(solver.x, (s.aRad * s.c .* T_exBUGH + epsilon^2 .* h_exBUGH),"--",color = "blue", label = L"BUG_{15}");
# ax2.plot(solver.x, (s.aRad * s.c .* T_exraBUG + epsilon^2 .* h_exraBUG),"--",color = "purple", label = "raBUG");
# ax2.plot(-1.0.*ones(100),collect(range(-100,1000,100)),"--",color = "black");
# ax2.plot(0.5.*ones(100),collect(range(-100,1000,100)),"--",color = "black");
# ax2.set_ylim([minimum(s.aRad * s.c .* T0 + epsilon^2 .* h0)-1.0,maximum(s.aRad * s.c .* T0 + epsilon^2 .* h0)+10.0]);
# ax2.set_xlim([-3,6]);
# ax2.set_xticks(range(-3,6, step=1.0));
# ax2.set_ylabel(L"Scalar flux, $\phi(t,x) = \frac{1}{2}\int_{-1}^{1}f(t,x,\mu) \mathrm{d}\mu$");
# ax2.set_xlabel(L"x");
# ax2.legend();
# # ax2.set_title(string(L"\sigma_{a} = ", s.sigmaA))
# fig2.canvas.draw();
# fig2
# savefig("Scalar_flux_lowrank.pdf")


# t = collect(range(0,s.Tend,size(energy_ex)[1]));
# t1 = collect(range(0,s.Tend,size(energy_ex5)[1]));
# t2 = collect(range(0,s.Tend,size(energy_ex15)[1]));

# fig3, ax3 = subplots(figsize=(15, 12), dpi=200);
# start = 1;
# ax3.semilogy(t[start:end],energy_ex[start:end],"--",color = "red", label = L"P_{100}");
# # ax3.plot(t1[start:end],energy_ex5[start:end],"--",color = "green", label = L"P_{5}");
# # ax3.plot(t2[start:end],energy_ex15[start:end],"--",color = "brown", label = L"P_{15}");

# # ax3.semilogy(t[start:end],energy_exBUG[start:end],"--",color = "orange", label = L"BUG_{5}");
# # ax3.plot(t[start:end],energy_exBUGH[start:end],"--",color = "blue", label = L"BUG_{15}");

# # ax3.semilogy(t[start:end],energy_exraBUG[start:end],"--",color = "purple",label = "raBUG");

# ax3.legend();
# ax3.set_xlim([0,s.Tend]);
# ax3.set_ylabel(L"e");
# ax3.set_xlabel(L"t");
# fig3.canvas.draw();
# fig3
# savefig("energy.pdf")


# fig,ax = subplots(figsize = (15,12),dpi=200);
# ax.semilogy(t[2:end],mass_ex,"--", color = "red", label = "full");
# ax.semilogy(t[2:end],mass_exBUG,"--", color = "orange", label = "BUG");
# ax.semilogy(t[2:end],mass_exraBUG,"--", color = "purple", label = "raBUG");
# ax.set_xlim([0,s.Tend]);
# ax.set_ylabel("mass");
# ax.set_xlabel(L"t");
# ax.legend();
# fig.canvas.draw();
# fig


# fig4, ax4 = subplots(figsize=(15, 12), dpi=200);
# ax4.plot(t[1:end],ranks_exraBUG,"-",color = "purple", label = "rank");
# # ax4.legend();
# ax4.set_xlim([0,s.Tend]);
# ax4.set_ylim([0,maximum(ranks_exraBUG)+1]);
# ax4.set_yticks(range(1,maximum(ranks_exraBUG)+1,step=2))
# ax4.set_ylabel("rank");
# ax4.set_xlabel(L"t");
# fig4.canvas.draw();
# fig4
# savefig("raBUG_ranks.pdf")