__precompile__

using PyPlot
using LaTeXStrings


rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 30


function plothyperbolic(s,solver,)
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
    # savefig("Temperature_moment.pdf")

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
    # savefig("Temperature_low_rank.pdf")


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
    # savefig("Scalar_flux_moment.pdf")

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
    fig2.canvas.draw();
    fig4
    # savefig("Scalar_flux_lowrank.pdf")


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
    # savefig("energy.pdf")


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
    # savefig("raBUG_ranks.pdf")


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
    # savefig("raBUG_ranks.pdf")
end