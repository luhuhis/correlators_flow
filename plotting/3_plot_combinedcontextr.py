#!/usr/local/bin/python
import numpy
import matplotlib
import lib_plot as pl

fig, ax, plots = pl.create_figure(xlims=[0.15,0.5], ylims=[2,4], xlabel=r'$\tau T$', ylabel=r'$\displaystyle \frac{G_{\tau_F}^\mathrm{cont } (\tau)}{G^\mathrm{free }_{\tau_F=0} (\tau)}$')
pl.legendstyle.update(dict(loc='center right', bbox_to_anchor=(1,0.4) ))

flowstart=10
flowend=21


#for i in range(0,len(pl.flow_radius)):
    #pl.flow_radius[i] = round(pl.flow_radius[i],4)
flow_index_range = range(flowstart,flowend)
for i in flow_index_range:
    EE_cont = numpy.loadtxt(pl.inputfolder+"/continuum_limit/EE_"+'{0:.4f}'.format(pl.flow_radius[i])+"_cont.txt", skiprows=pl.n_skiprows)
    for j in range(0,len(EE_cont)):
        if EE_cont[j,0] < pl.flow_radius[i]/numpy.sqrt(8*0.014)+pl.offset:
            EE_cont[j,:] = None
    plots.append(ax.fill_between(EE_cont[:,0], EE_cont[:,1]-EE_cont[:,2], EE_cont[:,1]+EE_cont[:,2], label='{0:.3f}'.format(pl.flow_radius[i]), facecolor=pl.get_color(pl.flow_radius, i, flowstart, flowend), zorder=-2*len(flow_index_range)+i))
ax.legend(handles=plots, labels=['{0:.3f}'.format(i) for i in pl.flow_radius[flowstart:flowend]], title=r'$\sqrt{8\tau_F}T$', **pl.legendstyle)

#plot zoomed plot
zoom_plots=[]
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
axins = zoomed_inset_axes(parent_axes=ax, zoom=2.5, loc='lower center')
axins.set_xlim(0.36,0.39)
axins.set_ylim(3.2, 3.45)
axins.xaxis.set_visible(False)
axins.yaxis.set_visible(False)
#plot zoom
for i in flow_index_range:
    EE_cont = numpy.loadtxt(pl.inputfolder+"/continuum_limit/EE_"+'{0:.4f}'.format(pl.flow_radius[i])+"_cont.txt", skiprows=pl.n_skiprows)
    for j in range(0,len(EE_cont)):
        if EE_cont[j,0] < pl.flow_radius[i]/numpy.sqrt(8*0.014)+pl.offset:
            EE_cont[j,:] = None
    zoom_plots.append(axins.fill_between(EE_cont[:,0], EE_cont[:,1]-EE_cont[:,2], EE_cont[:,1]+EE_cont[:,2], label='{0:.3f}'.format(pl.flow_radius[i]), facecolor=pl.get_color(pl.flow_radius, i, flowstart, flowend), zorder=-2*len(flow_index_range)+i))
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
mark_inset(parent_axes=ax, inset_axes=axins, linewidth=1, alpha=0.5, loc1=2, loc2=4, zorder=1000)

fig.savefig(pl.outputfolder+"/EE_aeq0_tneq0.pdf", **pl.figurestyle) 
print("saved continuum extr plot", pl.outputfolder+"/EE_aeq0_tneq0.pdf") 
