#!/usr/local/bin/python
import numpy
import re
import sys
import matplotlib
import lib_plot as pl
import lib_process_data as lpd

def skip_large_errors(EE, EE_err, boundary):
    for i in range(0,EE.shape[0]):
        for j in range(0,EE.shape[1]):
            if EE_err[i,j] > boundary :
                EE[i,j] = None

qcdtype, conftype, beta, ns, nt, nt_half = lpd.read_args()

 
"""load data"""
inputfolder=pl.get_inputfolder(qcdtype)+conftype+"/"
outputfolder=pl.get_outputfolder(qcdtype)+conftype+"/"
lpd.create_folder(outputfolder) 
flow_radius = numpy.loadtxt(inputfolder+"flowradii_"+conftype+".dat")
EE        = numpy.loadtxt(inputfolder+"EE_"+conftype+".dat")
EE_err = numpy.loadtxt(inputfolder+"EE_err_"+conftype+".dat")
for i in range(len(flow_radius)):
#for i in range(150):
    for j in range(nt_half):
        EE[i,j] = EE[i,j] * lpd.improve_corr_factor(j, nt, i)
        EE_err[i,j] = EE_err[i,j] * lpd.improve_corr_factor(j, nt, i)
tauT = lpd.tauT_imp[nt]
 
 
#"""PLOT: x-axis tauT, flowtimes legend"""
#fig, ax, plots = pl.create_figure(xlims=[0,0.5], xlabel=r'$\tau T$', ylabel=r'$\displaystyle \frac{G_{\tau_F} (\tau)}{G^\mathrm{ norm } (\tau)}$')
##pl.legendstyle.update(dict(loc="upper right"))
##ax.yaxis.set_label_coords(0.01,0.99)

##skip_large_errors(EE, EE_err, EE_backup, EE_err_backup, 10000) 
#flowstart=0; 
#if qcdtype == "quenched":
    #flowend=16
    ##ax.set_ylim([0.7,4])
    #ax.set_ylim([-1,4])
#if qcdtype == "hisq":
    #flowstart=5
    #flowend=21
    #ax.set_ylim([0,12])
      
#flow_selection = range(0,150,10)

#for i in flow_selection:
    #plots.append(ax.fill_between(list(tauT), EE[i,:]-EE_err[i,:], EE[i,:]+EE_err[i,:], facecolor=pl.get_color(flow_radius, i, flow_selection[0], flow_selection[-1]), linewidth=pl.mylinewidth))
    #ax.errorbar(list(tauT), EE[i,:], color=pl.get_color(flow_radius, i, flow_selection[0], flow_selection[-1]), **pl.plotstyle_lines)

#leg = ax.legend(handles=plots, labels=['{0:.3f}'.format(j) for i,j in enumerate(flow_radius) if i in flow_selection], title=r"$ \sqrt{8\tau_F}T$", **pl.legendstyle)
#matplotlib.pyplot.tight_layout(0)
#fig.savefig(outputfolder+conftype+"_EE_fill.pdf") 
#print("saved correlator plot", outputfolder+conftype+"_EE_fill.pdf")
#ax.lines.clear() ; ax.collections.clear() ; plots.clear()


"""PLOT: x-axis flowtimes, tautT legend"""
if qcdtype == "hisq":
    ylims=([0,12])
else:
    ylims=([1,4])
fig, ax, plots = pl.create_figure(xlims=[-0.0001,0.0028125+0.00015], ylims=ylims, 
                                  xlabelpos=[0.94,0.05],
                                  ylabelpos=[0.01,0.97],
                                  xlabel=r'$ \tau_F T^2 $', 
                                  ylabel=r'$\displaystyle\frac{G^\mathrm{latt }_{\tau}(\tau_F)}{G_{\tau,{\tau}_{ F }=0}^{\substack{ \text{\tiny  norm} \\[-0.4ex] \text{\tiny latt } } } }$')
pl.plotstyle_add_point.update(dict(fmt=',', capsize=0.6))


"""calculate flow limits"""
flow_limit=numpy.zeros(len(tauT))
interpolation_value=numpy.zeros(len(tauT))
for i in range(len(tauT)):
    flow_limit[i] = lpd.upper_flow_limit(tauT[i])
    index = (numpy.abs(flow_radius[:] - flow_limit[i])).argmin()
    offset = 1 if flow_limit[i]-flow_radius[index] > 0 else -1
    offset2 = -1 if offset == -1 else 0
    interpolation_value[i] = EE[index+offset2,i]-abs(EE[index,i]-EE[index+offset,i])/abs(flow_radius[index]-flow_radius[index+offset])*abs(flow_limit[i]-flow_radius[index+offset2])

max_flow_index = 134 #this is the max because of the coarsest lattice (too few data points to put a spline through at this flow time)
for i in range(nt_half):
    flow_extr_filter_high = 0
    flow_extr_filter_low = -1
    flow_extr_filter_low_grey = -1
    rel_err = numpy.fabs(EE_err[:,i]/EE[:,i])
    for j,val in enumerate(rel_err):
        if flow_extr_filter_high == 0:
            if flow_limit[i] < flow_radius[j]:
                flow_extr_filter_high = j
        if flow_extr_filter_low == -1:
            if val <= 0.0075 *tauT[i]:
                flow_extr_filter_low = max(j,50)
        if flow_extr_filter_low_grey == -1:
            if val <= 0.05:
                flow_extr_filter_low_grey = j
    flow_extr_filter_high = min(flow_extr_filter_high,max_flow_index)
    
    xdata = flow_radius[flow_extr_filter_low:flow_extr_filter_high]**2/8
    ydata = EE[flow_extr_filter_low:flow_extr_filter_high,i]
    edata = EE_err[flow_extr_filter_low:flow_extr_filter_high,i]
    ax.errorbar(xdata, ydata, edata, color=pl.get_color(tauT, i, 0, nt_half), zorder=(-i), **pl.plotstyle_add_point)
    
    xdata = flow_radius[flow_extr_filter_low_grey:]**2/8
    ydata = EE[flow_extr_filter_low_grey:,i]
    edata = EE_err[flow_extr_filter_low_grey:,i]
    #ax.errorbar(xdata, ydata, color='lightgrey', zorder=(-100*i), fmt='-', lw=0.5, markersize=0.5)  
    ax.fill_between(xdata, ydata-edata, ydata+edata, linewidth=0.5, color=(0.9,0.9,0.9,1), zorder=(-100*i))  #
    ax.errorbar(xdata, ydata, linewidth=0.5, color=(0.7,0.7,0.7,1), zorder=(-100*i+1))  #
    ax.axvline(x=flow_radius[50]**2/8, **pl.verticallinestyle)
    plots.append(ax.errorbar(flow_limit[i]**2/8, interpolation_value[i], marker=pl.markers[i], fillstyle='none', markersize=4, mew=0.25, color=pl.get_color(tauT, i), label='{0:.3f}'.format(tauT[i]), capsize=0, lw=0 ))
    ax.errorbar(flow_limit[i]**2/8, interpolation_value[i], marker='|', fillstyle='none', markersize=4, mew=0.25, color=pl.get_color(tauT, i, 0, nt_half), capsize=0 )
    
ax2=ax.twiny()
new_tick_locations = numpy.array([0,0.05, 0.08, 0.1, 0.12, 0.13, 0.14, 0.15])**2/8
#new_tick_locations = numpy.array([0,0.05**2/8, 0.08**2/8, 0.1**2/8, 0.12**2/8, 0.13**2/8, 0.14**2/8, 0.15**2/8])
def tick_function(X):
    V = numpy.sqrt(X*8)
    return ["%.1f" % z if z == 0 else "%.2f" % z for z in V ]
ax2.tick_params(direction='in', pad=0, width=0.5)
ax2.set_xlim(ax.get_xlim())
ax2.set_xticks(new_tick_locations)
ax2.set_xticklabels(tick_function(new_tick_locations))
ax2.set_xlabel(r'$\sqrt{8\tau_F}T$')#, horizontalalignment='right', verticalalignment='top', bbox=pl.labelboxstyle, zorder=99)
ax2.xaxis.set_label_coords(0.92,0.92)


pl.legendstyle.update(dict(loc="center right", bbox_to_anchor=(1,0.4), columnspacing=0.3, handletextpad=0.2, handlelength=0.5, handler_map={matplotlib.container.ErrorbarContainer: matplotlib.legend_handler.HandlerErrorbar(xerr_size=0)}))
ax.legend(handles=plots)    
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::-1], labels[::-1], ncol=2, title=r"$\tau T=$", **pl.legendstyle)
matplotlib.pyplot.tight_layout(0)
### change first xtick label
fig.canvas.draw()
ticks=ax.get_xticks().tolist()
ticks[0]='0.0'
ax.set_xticklabels(ticks)
fig.savefig(outputfolder+"/"+conftype+"_EE_flowlim_fill.pdf") 
print("saved flow effect plot", outputfolder+"/"+conftype+"_EE_flowlim_fill.pdf")
