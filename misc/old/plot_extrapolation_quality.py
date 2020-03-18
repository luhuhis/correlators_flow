#!/usr/local/bin/python
import math
import numpy
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from matplotlib.legend_handler import HandlerErrorbar
from os import listdir
import sys
from latqcdtools import fitting as fit

def get_color(myarray, i, start, end):
    return cm.gnuplot((myarray[i]-myarray[start])/(myarray[end-1]-myarray[start]))

qcdtype="quenched"
inputfolder="../data_merged/"+qcdtype+"/"


flow_radius = numpy.loadtxt(inputfolder+"s064t16_b0687361/flowradius_s064t16_b0687361.dat")


figurestyle        = dict(bbox_inches="tight")
plotstyle_points   = dict(fmt='D', linewidth=1, markersize=4, capsize=2, mew=0.5, fillstyle='none')
labelboxstyle      = dict(boxstyle="square", fc="w", ec='none', alpha=0.7, pad=0.15, zorder=999999)
legendstyle        = dict(loc="center right", frameon=True, framealpha=0.8, edgecolor='none', fancybox=False, facecolor="w", title="", labelspacing=0.1, borderpad=0.1, handletextpad=0.4)#,
#legendstyle        = dict(loc="center left", frameon=True, framealpha=0.8, edgecolor='none', fancybox=False, facecolor="w", title="", labelspacing=0.1, borderpad=0.1, handletextpad=0.4, prop={'size': 9})
                          #, handlelength=1)#, prop={'size': 9})
figurestyle        = dict(bbox_inches="tight", pad_inches=0)
plotstyle_fill     = dict(linewidth=0.5)
xlabelstyle        = dict(horizontalalignment='right', verticalalignment='bottom',  bbox=labelboxstyle, zorder=999998)
ylabelstyle        = dict(horizontalalignment='left', verticalalignment='top', rotation=0, bbox=labelboxstyle, 
                          zorder=999999)
titlestyle         = dict(x=0.5, y=0.95, bbox=labelboxstyle, verticalalignment='top', zorder=999999)


plt.rc('text', usetex=True)
plt.rc('text.latex')
plt.rc('font', family='serif', size='14')

fig = plt.figure() 
ax = fig.add_subplot(1,1,1) 
ax.xaxis.set_label_coords(0.99,0.01)
ax.yaxis.set_label_coords(0.01,0.97)

ax.set_xlim([-0.0001,1/15**2+0.001])
#ax.set_ylim([1.9,4])
ax.set_ylim([2.9,3.8])

plots = []
ax.set_ylabel(r'$\displaystyle \frac{G_{\tau,\tau_F} }{G^{\mathrm{free}}_{\tau,\tau_F=0} }$', **ylabelstyle)
ax.set_xlabel(r'$N_\tau^{-2}$', **xlabelstyle)

start=6
end=14
#be careful with the 

flowstart=20
flowend=21

for i in range(flowstart,flowend):
    for tauT_index in range(start,end):
        EE_16 = numpy.loadtxt(inputfolder+"continuum_limit/EE_"+'{0:.4f}'.format(flow_radius[i])+"_interpolation_16.txt", max_rows=14)
        EE_20 = numpy.loadtxt(inputfolder+"continuum_limit/EE_"+'{0:.4f}'.format(flow_radius[i])+"_interpolation_20.txt", max_rows=14)
        EE_24 = numpy.loadtxt(inputfolder+"continuum_limit/EE_"+'{0:.4f}'.format(flow_radius[i])+"_interpolation_24.txt", max_rows=14)
        EE_30 = numpy.loadtxt(inputfolder+"continuum_limit/EE_"+'{0:.4f}'.format(flow_radius[i])+"_interpolation_30.txt", max_rows=14)
        EE_cont = numpy.loadtxt(inputfolder+"continuum_limit/EE_"+'{0:.4f}'.format(flow_radius[i])+"_cont.txt")
        ax.set_title(r'$\sqrt{8\tau_F}T=$'+'{0:.2f}'.format(flow_radius[i]),**titlestyle) #r'$ \tau T = $' + '{0:.2f}'.format(EE_cont[tauT_index,1]) +
        
        data = []
        Nts = [16, 20, 24, 30]
        data.append([0, EE_cont[tauT_index,1], EE_cont[tauT_index,2]])
        data.append([1/16**2, EE_16[tauT_index,1], EE_16[tauT_index,2]])
        data.append([1/20**2, EE_20[tauT_index,1], EE_20[tauT_index,2]])
        data.append([1/24**2, EE_24[tauT_index,1], EE_24[tauT_index,2]])
        data.append([1/30**2, EE_30[tauT_index,1], EE_30[tauT_index,2]])
        
        _xdata = [j[0] for j in data[:]]
        _ydata = [j[1] for j in data[:]]
        _edata = [j[2] for j in data[:]]
        
        plots.append(ax.errorbar(_xdata, _ydata, _edata, **plotstyle_points, color=get_color(EE_cont[:,0], 
                                 tauT_index,start,end), zorder=-tauT_index-1000, label='{0:.3f}'.format(EE_cont[tauT_index,0])))
        
        #fit
        func = lambda x,a,b:a*x+b
        fitter = fit.Fitter(func, xdata = _xdata, ydata = _ydata, edata = _edata)
        res, res_err, chi_dof = fitter.do_fit(start_params = [0, 3])
        x = numpy.linspace(0,0.1,1000)
        ax.errorbar(x, func(x, *res), color='grey', alpha=0.8, fmt='--', lw=0.75, zorder=-100)
    ax.legend(handles=plots, loc="center right", frameon=True, framealpha=0.8, edgecolor='none', fancybox=False, facecolor="w", labelspacing=0.25, borderpad=0, handletextpad=0, title=r'$\tau T$')
    fig.savefig(inputfolder+"continuum_limit/EE_rF"+'{0:.4f}'.format(flow_radius[i])+"_interpolation.pdf", **figurestyle) 
    #"_tauT"+'{0:.4f}'.format(EE_cont[tauT_index,1])+
    ax.lines.clear() ; ax.collections.clear() ; plots.clear();
    print("done "+'{0:.4f}'.format(flow_radius[i]))
