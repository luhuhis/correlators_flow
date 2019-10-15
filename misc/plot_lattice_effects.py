#!/usr/local/bin/python
import math
import numpy
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from os import listdir
import sys
from latqcdtools import scales_quenched

outputfoldermisc="../plots/quenched/misc/"

inputfolder="../data_merged/quenched/"

EE_16 = numpy.loadtxt(inputfolder+"/s064t16_b0687361/single_flow/EE_0.075_Nt16.dat")
EE_20 = numpy.loadtxt(inputfolder+"/s080t20_b0703500/single_flow/EE_0.075_Nt20.dat")
EE_24 = numpy.loadtxt(inputfolder+"/s096t24_b0719200/single_flow/EE_0.075_Nt24.dat")
EE_30 = numpy.loadtxt(inputfolder+"/s120t30_b0739400/single_flow/EE_0.075_Nt30.dat")
EE_cont = numpy.loadtxt(inputfolder+"/continuum_limit/EE_0.075_cont.txt")
EE_cont_2015 = numpy.loadtxt(inputfolder+"/multi-level_2015/cont_thomas.dat")
EE_cont_2015_new = numpy.loadtxt(inputfolder+"/multi-level_2015/EE_2015_cont.txt")

plots = []
plotstyle_fill     = dict(linewidth=0.5)
plotstyle_points   = dict(fmt='D-', linewidth=1, markersize=4, capsize=2, mew=0.5, fillstyle='none')
legendstyle        = dict(bbox_to_anchor=(1,0.1), loc="lower right", frameon=True, framealpha=0, title="", prop={'size': 9})
figurestyle        = dict(bbox_inches="tight")
labelboxstyle      = dict(boxstyle="round", fc="w", ec='none', alpha=0.7)


plt.rc('text', usetex=True)
plt.rc('text.latex')
plt.rc('font', family='serif')

fig = plt.figure() 
ax = fig.add_subplot(1,1,1) 
ax.set_title(r'pure SU(3) $|$ $T\approx 1.5 T_C$ $|$ imp. dist.')
ax.xaxis.set_label_coords(0.975,0.050)
ax.yaxis.set_label_coords(0.025,0.975)
ax.set_ylabel(r'$\displaystyle \frac{G_{EE}(\tau T)}{G_{EE}^{\,\mathrm{free}}(\tau T)}$', 
              bbox=labelboxstyle, horizontalalignment='left', verticalalignment='top', rotation=0)

ax.set_xlim([0.25,0.5]); 
ax.set_ylim([2.3,4])
ax.set_xlabel(r'$\tau T$', horizontalalignment='right', bbox=labelboxstyle)
legendstyle["title"]=r'Flow $\scriptstyle \sqrt{\tilde{\tau}_F}T=0.075$'

plots.append(ax.errorbar(EE_16[:,0], EE_16[:,1], EE_16[:,2], **plotstyle_points))
plots.append(ax.errorbar(EE_20[:,0], EE_20[:,1], EE_20[:,2], **plotstyle_points))
plots.append(ax.errorbar(EE_24[:,0], EE_24[:,1], EE_24[:,2], **plotstyle_points))
plots.append(ax.errorbar(EE_30[:,0], EE_30[:,1], EE_30[:,2], **plotstyle_points))
plots.append(ax.fill_between(EE_cont[:,0], EE_cont[:,1]-EE_cont[:,2], EE_cont[:,1]+EE_cont[:,2], facecolor='black', alpha=0.3, **plotstyle_fill))
plots.append(ax.plot([],marker="", ls="")[0])
#plots.append(ax.fill_between(EE_cont_2015[:,0], EE_cont_2015[:,1]-EE_cont_2015[:,2], EE_cont_2015[:,1]+EE_cont_2015[:,2], facecolor='black', alpha=0.5, **plotstyle_fill))
plots.append(ax.fill_between(EE_cont_2015_new[:,0], EE_cont_2015_new[:,1]-EE_cont_2015_new[:,2], EE_cont_2015_new[:,1]+EE_cont_2015_new[:,2], facecolor='magenta', alpha=0.75, **plotstyle_fill))
   
leg = ax.legend(handles=plots, labels=list((r'$N_\tau=16$',r'$N_\tau=20$',r'$N_\tau=24$',r'$N_\tau=30$', "cont.", "Multi-Level (2015)", "cont.", "cont. \n"+r'new $Z_\mathrm{pert}$')), **legendstyle)
leg._legend_handle_box.get_children()[0].get_children()[5].get_children()[0].set_width(0)
fig.savefig(outputfoldermisc+"EE_lat_effects_flow0075.pdf", **figurestyle) 
