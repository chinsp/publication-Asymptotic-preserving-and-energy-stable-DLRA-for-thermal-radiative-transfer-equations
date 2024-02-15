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

#########################################################################
#                  Computing the Rosseland limit                        #
#########################################################################
Nx = 201;
N = 100;
epsilon = 1e-5;
regime = "parabolic"
s = Settings(Nx,N,epsilon,regime);
solver = solverMarshak(s);
t,T = solveRosselandLimit(solver);
#########################################################################

epsilon = 1.0;
regime = "hyperbolic"
Nx = 501;

N = 5;
s = Settings(Nx,N,epsilon,regime);
solver = solverMarshak(s);
t_5, h_5, g_5, T_5, energy_5 = solveFullMacroMicro(solver);


N = 15;
s = Settings(Nx,N,epsilon,regime);
solver = solverMarshak(s);
t_15, h_15, g_15, T_15, energy_15 = solveFullMacroMicro(solver);


#########################################################################
#            The numerical experiments for epsilon = 1.0                #
#########################################################################
N = 100;

s = Settings(Nx,N,epsilon,regime);
solver = solverMarshak(s);

t_Full, h_Full, g_Full, T_Full, energy_Full, Resmass_Full = solveFullMacroMicro(solver);

t_BUG, h_BUG, g_BUG, T_BUG, energy_BUG, Resmass_BUG = solveBUGintegrator(solver);

s.r = 15;
r_BUGH = s.r;
t_BUGH, h_BUGH, g_BUGH, T_BUGH, energy_BUGH, Resmass_BUGH = solveBUGintegrator(solver);

s.r = 1;
t_raBUG, h_raBUG, g_raBUG, T_raBUG, energy_raBUG, ranks_raBUG, Resmass_raBUG = solveBUG_rankadaptive(solver);


x1 = collect(range(s.a,s.b,size(T)[1]));

fig1, ax1 = subplots(figsize=(15, 12), dpi=200);
ax1.plot(x1, T,color = "black", "--", label = "Rosseland limit");
ax1.plot(solver.x, T_Full ,"--",color = "red", label = string(L"P_{100}"));
ax1.plot(solver.x, T_5 ,"--",color = "green", label = string(L"P_{5}"));
ax1.plot(solver.x, T_15 ,"--",color = "brown", label = string(L"P_{15}"));

ax1.set_xlim([0,3]);
ax1.set_xticks(range(0,3, step=0.5));
ax1.set_ylabel("Temperature");
ax1.set_xlabel(L"x");
ax1.legend();
fig1.canvas.draw();
fig1
# savefig("Temperature_moment_hyperbolic.pdf")

fig2, ax2 = subplots(figsize=(15, 12), dpi=200);
ax2.plot(x1, T,color = "black", "--", label = "Rosseland limit");
ax2.plot(solver.x, T_Full ,"--",color = "red", label = string(L"P_{100}"));
ax2.plot(solver.x, T_BUG ,"--",color = "orange", label = string(L"BUG_{5}"));
ax2.plot(solver.x, T_BUGH,"--",color = "blue", label = string(L"BUG_{15}"));

ax2.plot(solver.x, T_raBUG ,"--",color = "purple", label = string("raBUG"));

ax2.set_ylim([minimum(T_Full)-1.0,maximum(T_Full)+5.0]);
ax2.set_xlim([0,3]);
ax2.set_xticks(range(0,3, step=0.5));
ax2.set_ylabel("Temperature");
ax2.set_xlabel(L"x");
ax2.legend();
fig2.canvas.draw();
fig2
# savefig("Temperature_low_rank_hyperbolic.pdf")


fig3, ax3 = subplots(figsize=(15, 12), dpi=200);
ax3.plot(solver.x, (s.aRad * s.c .* T_Full + epsilon^2 .* h_Full),"--",color = "red", label = L"P_{100}");
ax3.plot(solver.x, (s.aRad * s.c .* T_5 + epsilon^2 .* h_5),"--",color = "green", label = L"P_{5}");
ax3.plot(solver.x, (s.aRad * s.c .* T_15 + epsilon^2 .* h_15),"--",color = "brown", label = L"P_{15}");

ax3.set_ylim([minimum(s.aRad * s.c .* T_Full + epsilon^2 .* h_Full)-1.0,maximum(s.aRad * s.c .* T_Full + epsilon^2 .* h_Full)+5.0]);
ax3.set_xlim([0,3]);
ax3.set_xticks(range(0,3, step=0.5));
ax3.set_ylabel(L"Scalar flux, $\phi(t,x) = \frac{1}{2}\int_{-1}^{1}f(t,x,\mu) \mathrm{d}\mu$");
ax3.set_xlabel(L"x");
ax3.legend();
fig3.canvas.draw();
fig3
# savefig("Scalar_flux_moment_hyperbolic.pdf")

fig4, ax4 = subplots(figsize=(15, 12), dpi=200);
ax4.plot(solver.x, (s.aRad * s.c .* T_Full + epsilon^2 .* h_Full),"--",color = "red", label = L"P_{100}");

ax4.plot(solver.x, (s.aRad * s.c .* T_BUG + epsilon^2 .* h_BUG),"--",color = "orange", label = L"BUG_{5}");
ax4.plot(solver.x, (s.aRad * s.c .* T_BUGH + epsilon^2 .* h_BUGH),"--",color = "blue", label = L"BUG_{15}");

ax4.plot(solver.x, (s.aRad * s.c .* T_raBUG + epsilon^2 .* h_raBUG),"--",color = "purple", label = "raBUG");

ax4.set_ylim([minimum(s.aRad * s.c .* T_Full + epsilon^2 .* h_Full)-1.0,maximum(s.aRad * s.c .* T_Full + epsilon^2 .* h_Full)+5.0]);
ax4.set_xlim([0,3]);
ax4.set_xticks(range(0,3, step=1.0));
ax4.set_ylabel(L"Scalar flux, $\phi(t,x) = \frac{1}{2}\int_{-1}^{1}f(t,x,\mu) \mathrm{d}\mu$");
ax4.set_xlabel(L"x");
ax4.legend();
fig4.canvas.draw();
fig4
# savefig("Scalar_flux_lowrank_hyperbolic.pdf")


t = collect(range(0,s.Tend,size(energy_Full)[1]));
t1 = collect(range(0,s.Tend,size(energy_5)[1]));
t2 = collect(range(0,s.Tend,size(energy_15)[1]));

fig5, ax5 = subplots(figsize=(15, 12), dpi=200);
start = 1;
ax5.semilogy(t[start:end],energy_Full[start:end],"--",color = "red", label = L"P_{100}");
ax5.semilogy(t1[start:end],energy_5[start:end],"--",color = "green", label = L"P_{5}");
ax5.semilogy(t2[start:end],energy_15[start:end],"--",color = "brown", label = L"P_{15}");

ax5.semilogy(t[start:end],energy_BUG[start:end],"--",color = "orange", label = L"BUG_{5}");
ax5.semilogy(t[start:end],energy_BUGH[start:end],"--",color = "blue", label = L"BUG_{15}");

ax5.semilogy(t[start:end],energy_raBUG[start:end],"--",color = "purple",label = "raBUG");

ax5.legend();
ax5.set_xlim([0,s.Tend]);
ax5.set_ylabel(L"energy");
ax5.set_xlabel(L"t");
fig5.canvas.draw();
fig5
# savefig("energy_hyperbolic.pdf")


fig6,ax6 = subplots(figsize = (15,12),dpi=200);
ax6.semilogy(t[2:end],Resmass_Full,"--", color = "red", label = L"P_{100}");
ax6.semilogy(t[2:end],Resmass_BUG,"--", color = "orange", label = L"BUG_{5}");
ax6.semilogy(t[2:end],Resmass_raBUG,"--", color = "purple", label = "raBUG");
ax6.set_xlim([0,s.Tend]);
ax6.set_ylabel(L"Relative residual mass, $ \frac{|m^{0} - m^{n}|}{|m^{0}|} $");
ax6.set_xlabel(L"t");
ax6.legend();
fig6.canvas.draw();
fig6
# savefig("Rel_residual_mass_hyperbolic.pdf")


fig7, ax7 = subplots(figsize=(15, 12), dpi=200);
ax7.plot(t[1:end],ranks_raBUG,"-",color = "purple", label = "rank");

ax7.set_xlim([0,s.Tend]);
ax7.set_ylim([0,maximum(ranks_raBUG)+1]);
ax7.set_yticks(range(1,maximum(ranks_raBUG)+1,step=2))
ax7.set_ylabel("rank");
ax7.set_xlabel(L"t");
# ax4.legend();
fig7.canvas.draw();
fig7
# savefig("raBUG_ranks_hyperbolic.pdf")


#########################################################################
#            The numerical experiments for epsilon = 1e-5               #
#########################################################################
Nx = 201;
N = 100;
epsilon = 1e-5;
regime = "parabolic"
s = Settings(Nx,N,epsilon,regime);
solver = solverMarshak(s);

t_Full, h_Full, g_Full, T_Full, energy_Full, Resmass_Full = solveFullMacroMicro(solver);

s.r = 1;
t_BUG, h_BUG, g_BUG, T_BUG, energy_BUG, Resmass_BUG = solveBUGintegrator(solver);

s.r = 1;
t_raBUG, h_raBUG, g_raBUG, T_raBUG, energy_raBUG, ranks_raBUG, Resmass_raBUG = solveBUG_rankadaptive(solver);




fig8, ax8 = subplots(figsize=(15, 12), dpi=200);
ax8.plot(x1, s.aRad * s.c * T,color = "black", "--", label = "Rosseland limit");
ax8.plot(solver.x, (s.aRad * s.c .* T_Full + epsilon^2 .* h_Full),"--",color = "red", label = L"P_{100}");
ax8.plot(solver.x, (s.aRad * s.c .* T_BUG + epsilon^2 .* h_BUG),"--",color = "orange", label = L"BUG_{1}");
ax8.plot(solver.x, (s.aRad * s.c .* T_raBUG + epsilon^2 .* h_raBUG),"--",color = "purple", label = "raBUG");
ax8.set_ylim([minimum(s.aRad * s.c .* T_Full + epsilon^2 .* h_Full)-1.0,maximum(s.aRad * s.c .* T_Full + epsilon^2 .* h_Full)+5.0]);
ax8.set_xlim([0,3]);
ax8.set_xticks(range(0,3, step=1.0));
ax8.set_ylabel(L"Scalar flux, $\phi(t,x) = \frac{1}{2}\int_{-1}^{1}f(t,x,\mu) \mathrm{d}\mu$");
ax8.set_xlabel(L"x");
ax8.legend();
fig8.canvas.draw();
fig8
# savefig("Scalar_flux_diffusive.pdf")

fig9, ax9 = subplots(figsize=(15, 12), dpi=200);
ax9.plot(x1, T,color = "black", "--", label = "Rosseland limit");
ax9.plot(solver.x, T_Full ,"--",color = "red", label = string(L"P_{100}"));
ax9.plot(solver.x, T_BUG ,"--",color = "orange", label = string(L"BUG_{1}"));
ax9.plot(solver.x, T_raBUG ,"--",color = "purple", label = string("raBUG"));
ax9.set_ylim([minimum(T_Full)-1.0,maximum(T_Full)+5.0]);
ax9.set_xlim([0,3]);
ax9.set_xticks(range(0,3, step=0.5));
ax9.set_ylabel("Temperature");
ax9.set_xlabel(L"x");
ax9.legend();
fig9.canvas.draw();
fig9
# savefig("Temperature_diffusive.pdf")

t = collect(range(0,s.Tend,size(energy_Full)[1]));

fig10, ax10 = subplots(figsize=(15, 12), dpi=200);
start = 1;
ax10.semilogy(t[start:end],energy_Full[start:end],"--",color = "red", label = L"P_{100}",linewidth = 4);
ax10.semilogy(t[start:end],energy_BUG[start:end],"--",color = "orange", label = L"BUG_{1}",linewidth = 4);
ax10.semilogy(t[start:end],energy_raBUG[start:end],"--",color = "purple",label = "raBUG",linewidth = 4);
ax10.legend();
ax10.set_xlim([0,s.Tend]);
ax10.set_ylabel(L"energy");
ax10.set_xlabel(L"t");
fig10.canvas.draw();
fig10
# savefig("energy_diffusive.pdf")


fig11,ax11 = subplots(figsize = (15,12),dpi=200);
ax11.semilogy(t[2:end],Resmass_Full,"--", color = "red", label = L"P_{100}",linewidth = 4);
ax11.semilogy(t[2:end],Resmass_BUG,"--", color = "orange", label = L"BUG_{1}",linewidth = 4);
ax11.semilogy(t[2:end],Resmass_raBUG,"--", color = "purple", label = "raBUG",linewidth = 4);
ax11.set_xlim([0,s.Tend]);
ax11.set_ylabel(L"Relative residual mass, $ \frac{|m^{0} - m^{n}|}{|m^{0}|} $");
ax11.set_xlabel(L"t");
ax11.legend();
fig11.canvas.draw();
fig11
# savefig("Rel_residual_mass_diffusive.pdf")

fig12, ax12 = subplots(figsize=(15, 12), dpi=200);
ax12.plot(t[1:end],ranks_raBUG,"-",color = "purple", label = "rank");
ax12.set_xlim([0,s.Tend]);
ax12.set_ylim([0,maximum(ranks_raBUG)+1]);
ax12.set_yticks(range(1,maximum(ranks_raBUG)+1,step=2))
ax12.set_ylabel("rank");
ax12.set_xlabel(L"t");
fig12.canvas.draw();
fig12
# savefig("raBUG_ranks_diffusive.pdf")

# h0,g0,T0 = IC(s);
# fig,ax = subplots(figsize = (15,12), dpi = 200)
# ax.plot(x1, s.aRad * s.c * T,color = "red", "--", label = "Diffusive regime",linewidth = 6);
# ax.plot(solver.x, (s.aRad * s.c .* T_Full + h_Full),"--",color = "purple", label = "Hyperbolic regime",linewidth = 6);
# ax.plot(solver.x, (s.aRad * s.c .* T0 .+ h0),"--",color = "black", label = "Initial distribution",linewidth = 6);
# ax.set_ylim([minimum(T0)-1.0,maximum(T0)+5.0]);
# ax.set_xlim([0,3]);
# ax.set_xticks(range(0,3, step=0.5));
# ax.set_ylabel("Particle distribution");
# ax.set_xlabel(L"x");
# ax.legend();
# fig.canvas.draw();
# fig
# savefig("presentation_plot.pdf")