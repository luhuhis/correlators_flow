#!/usr/local/bin/python
import numpy
import matplotlib
import plotlib as pl

fig, ax, plots = pl.create_figure(xlims=[0,0.5], ylims=[0.75,4.25], xlabel=r'$\tau T$', ylabel=r'$\displaystyle \frac{G_{\tau_F}(\tau)}{G_{\tau_F = 0}^{\mathrm{ norm } }(\tau)}$')

pl.titlestyle.update(dict(y=0.96))
pl.legendstyle.update(dict(loc="upper left", framealpha=0.5, bbox_to_anchor=(0,0.65), handler_map={matplotlib.container.ErrorbarContainer: matplotlib.legend_handler.HandlerErrorbar(xerr_size=0.4)}))
ax.xaxis.set_label_coords(0.99,0.01)
ax.yaxis.set_label_coords(0.01,0.97)

#for i in range(0,len(flow_radius)):
#for i in (0,):
for i in (0,5,10,15,20,22):
    ax.set_title(r'$\sqrt{8\tau_F} T=$ '+'{0:.3f}'.format(pl.flow_radius[i]), **pl.titlestyle)
    plots.append(ax.errorbar(pl.tauT_16, pl.EE_16[i,:], pl.EE_err_16[i,:], label=r'$N_\tau = 16$', **pl.plotstyle_points, color=matplotlib.cm.gnuplot(0.03), zorder=-3))
    plots.append(ax.errorbar(pl.tauT_20, pl.EE_20[i,:], pl.EE_err_20[i,:], label=r'$N_\tau = 20$', **pl.plotstyle_points, color=matplotlib.cm.gnuplot(0.40), zorder=-2))
    plots.append(ax.errorbar(pl.tauT_24, pl.EE_24[i,:], pl.EE_err_24[i,:], label=r'$N_\tau = 24$', **pl.plotstyle_points, color=matplotlib.cm.gnuplot(0.65), zorder=-1))
    plots.append(ax.errorbar(pl.tauT_30, pl.EE_30[i,:], pl.EE_err_30[i,:], label=r'$N_\tau = 30$', **pl.plotstyle_points, color=matplotlib.cm.gnuplot(0.90), zorder=0))
    try:
        EE_cont = numpy.loadtxt(pl.inputfolder+"/continuum_limit/EE_"+'{0:.4f}'.format(pl.flow_radius[i])+"_cont.txt", skiprows=pl.n_skiprows)
        for j in range(0,len(EE_cont)):
            if EE_cont[j,0] < pl.flow_radius[i]/numpy.sqrt(8*0.014)+pl.offset:
                EE_cont[j,:] = None
        plots.append(ax.fill_between(EE_cont[:,0], EE_cont[:,1]-EE_cont[:,2], EE_cont[:,1]+EE_cont[:,2], facecolor='grey', zorder=-6))
        ax.legend(handles=plots, labels=[r'$N_\tau = 16$', r'$N_\tau = 20$', r'$N_\tau = 24$', r'$N_\tau = 30$', r'cont.'], title="", **pl.legendstyle)
        #plot a minimum linewidth of cont extr for visibility
        ax.errorbar(EE_cont[:,0], EE_cont[:,1], color='grey', lw=1.5, fmt='-', zorder=-5, solid_capstyle='butt')
    except:
        ax.legend(handles=plots, labels=[r'$N_\tau = 16$', r'$N_\tau = 20$', r'$N_\tau = 24$', r'$N_\tau = 30$'], title="", **pl.legendstyle)
        #pass
    if i != 0:
        ax.axvline(x=(pl.flow_radius[i]/numpy.sqrt(8*0.014)+pl.offset), **pl.verticallinestyle)
    fig.savefig(pl.outputfolder+"/single_flow/EE_flow_"+'{0:.4f}'.format(pl.flow_radius[i])+".pdf", **pl.figurestyle) 
    print("saved lattice effects plot", pl.outputfolder+"/single_flow/EE_flow_"+'{0:.4f}'.format(pl.flow_radius[i])+".pdf") 
    ax.lines.clear() ; ax.collections.clear() ; plots.clear();
    ax.set_title(r'')

print("done with lattice effects plot") 
