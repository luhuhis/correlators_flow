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
inputfolder=pl.inputfolder+conftype+"/"
outputfolder=pl.outputfolder+conftype+"/"
lpd.create_folder(outputfolder) 
flow_radius = numpy.loadtxt(inputfolder+"flowradii_"+conftype+".dat")
EE        = numpy.loadtxt(inputfolder+"EE_"+conftype+".dat")
EE_backup = numpy.loadtxt(inputfolder+"EE_"+conftype+".dat")
EE_err = numpy.loadtxt(inputfolder+"EE_err_"+conftype+".dat")
EE_err_backup = numpy.loadtxt(inputfolder+"EE_err_"+conftype+".dat")
n_datafiles, n_streams=[int(i) for i in numpy.loadtxt(inputfolder+"n_datafiles_"+conftype+".dat")]
tauT = lpd.tauT_imp[nt]
 
 
"""PLOT: x-axis tauT, flowtimes legend"""
fig, ax, plots = pl.create_figure(xlims=[0,0.5], xlabel=r'$\tau T$', ylabel=r'$\displaystyle \frac{G_{\tau_F} (\tau)}{G^\mathrm{ norm } (\tau)}$')
#pl.legendstyle.update(dict(loc="upper right"))
#ax.yaxis.set_label_coords(0.01,0.99)

#skip_large_errors(EE, EE_err, EE_backup, EE_err_backup, 10000)
start=0; end=int(nt/2); 
flowstart=0; 
if qcdtype == "quenched":
    flowend=16
    #ax.set_ylim([0.7,4])
    ax.set_ylim([-1,4])
if qcdtype == "hisq":
    flowstart=5
    flowend=21
    ax.set_ylim([0,12])
      
flow_selection = range(0,150,1)

for i in flow_selection:
    plots.append(ax.fill_between(list(tauT), EE[i,:]-EE_err[i,:], EE[i,:]+EE_err[i,:], facecolor=pl.get_color(flow_radius, i, flow_selection[0], flow_selection[-1]), linewidth=pl.mylinewidth))
    ax.errorbar(list(tauT), EE[i,:], color=pl.get_color(flow_radius, i, flow_selection[0], flow_selection[-1]), **pl.plotstyle_lines)

leg = ax.legend(handles=plots, labels=['{0:.3f}'.format(j) for i,j in enumerate(flow_radius) if i in flow_selection], title=r"$ \sqrt{8\tau_F}T$", **pl.legendstyle)
fig.savefig(outputfolder+conftype+"_EE_fill.pdf", **pl.figurestyle) 
print("saved correlator plot", outputfolder+conftype+"_EE_fill.pdf")
ax.lines.clear() ; ax.collections.clear() ; plots.clear()




"""PLOT: x-axis flowtimes, tautT legend"""
#pl.legendstyle.update(dict(loc="center right"))
#ax.yaxis.set_label_coords(0.01,0.98)

"""calculate flow limits"""
flow_limit=numpy.zeros(len(tauT))
interpolation_value=numpy.zeros(len(tauT))
for i in range(0,len(tauT)):
    flow_limit[i] = tauT[i]*numpy.sqrt(8*0.014)
    index = (numpy.abs(flow_radius[:] - flow_limit[i])).argmin()
    offset = 1 if flow_limit[i]-flow_radius[index] > 0 else -1
    offset2 = -1 if offset == -1 else 0
    interpolation_value[i] = EE[index+offset2,i]-abs(EE[index,i]-EE[index+offset,i])/abs(flow_radius[index]-flow_radius[index+offset])*abs(flow_limit[i]-flow_radius[index+offset2])
flow_limit_lo=numpy.zeros(len(tauT))
interpolation_value_lo=numpy.zeros(len(tauT))
for i in range(0,len(tauT)):
    flow_limit_lo[i] = tauT[i]*numpy.sqrt(8*0.014)/3
    index = (numpy.abs(flow_radius[:] - flow_limit_lo[i])).argmin()
    offset = 1 if flow_limit_lo[i]-flow_radius[index] > 0 else -1
    offset2 = -1 if offset == -1 else 0
    interpolation_value_lo[i] = EE[index+offset2,i]-abs(EE[index,i]-EE[index+offset,i])/abs(flow_radius[index]-flow_radius[index+offset])*abs(flow_limit_lo[i]-flow_radius[index+offset2])


ax.set_xlim([0,0.003])
if qcdtype == "quenched":
    ax.set_ylim([0.5,4])
if qcdtype == "hisq":
    ax.set_ylim([0,12])
ax.set_xlabel(r'$ \tau_F T^2 $')
ax.set_ylabel(r'$\displaystyle\frac{G_{\tau}(\tau_F)}{G_{\tau,{\tau}_{ F }=0}^{\mathrm{ norm } } }$')

for i in range(start,end):
    plots.append(ax.fill_between(flow_radius**2/8, EE[:,i]-EE_err[:,i], 
                                 EE[:,i]+EE_err[:,i], facecolor=pl.get_color(tauT, i, start, end), zorder=(-(20+i)), linewidth=pl.mylinewidth))
    ax.errorbar(flow_radius**2/8, EE[:,i], color=pl.get_color(tauT, i, start, end), zorder=(-i), **pl.plotstyle_add_point)
    ax.errorbar(flow_limit[i]**2/8, interpolation_value[i], color=pl.get_color(tauT, i, start, end), **pl.flowlimitplotstyle)
    
ax2=ax.twiny()
new_tick_locations = numpy.array([0,0.05**2/8,0.075**2/8,0.1**2/8,0.125**2/8,0.15**2/8])
def tick_function(X):
    V = numpy.sqrt(X*8)
    return ["%.3f" % z for z in V]

ax2.set_xlim(ax.get_xlim())
ax2.set_xticks(new_tick_locations)
ax2.set_xticklabels(tick_function(new_tick_locations))
ax2.set_xlabel(r'$\sqrt{8\tau_F}T$')#, horizontalalignment='right', verticalalignment='top', bbox=pl.labelboxstyle, zorder=99)
#ax2.xaxis.set_label_coords(0.99,0.98)
    
ax.legend(handles=plots, labels=['{0:.3f}'.format(j) for j in tauT[start:end]], title=r"$\tau T$", **pl.legendstyle)
fig.savefig(outputfolder+"/"+conftype+"_EE_flowlim_fill.pdf", **pl.figurestyle) 
print("saved flow effect plot", outputfolder+"/"+conftype+"_EE_flowlim_fill.pdf")
ax.lines.clear() ; ax.collections.clear() ; plots.clear()
