#!/usr/local/bin/python
import math
import numpy
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from matplotlib.legend_handler import HandlerErrorbar
from os import listdir
import sys
from latqcdtools import scales_quenched

#a in 1/GeV
a_quenched = 0.0176
a_hisq = 0.0716

qcdtype="quenched"
outputfolder="../plots/"+qcdtype+"/"
inputfolder="../data_merged/"+qcdtype+"/"

EE_quenched = numpy.loadtxt(inputfolder+"/s096t32_b0719200/EE_s096t32_b0719200.dat")
EE_err_quenched = numpy.loadtxt(inputfolder+"/s096t32_b0719200/EE_err_s096t32_b0719200.dat")
flow_radius_quenched = list(numpy.around(numpy.loadtxt(inputfolder+"/s096t32_b0719200/flowradius_s096t32_b0719200.dat"), 6)) #need to change this to flow -times-
flow_times_quenched = [(32*i)**2 /8 for i in flow_radius_quenched] #they are the same for all these lattices
tauT_32 = list(numpy.around(numpy.loadtxt(inputfolder+"/s096t32_b0719200/tauT_imps096t32_b0719200.dat"), 4))
tau_32 = [i* 32 * a_quenched for i in tauT_32]

qcdtype="hisq"
outputfolder="../plots/"+qcdtype+"/"
inputfolder="../data_merged/"+qcdtype+"/"

EE_hisq = numpy.loadtxt(inputfolder+"/s064t16_b07188_m00113_m0306/EE_s064t16_b07188_m00113_m0306.dat")
EE_err_hisq = numpy.loadtxt(inputfolder+"/s064t16_b07188_m00113_m0306/EE_err_s064t16_b07188_m00113_m0306.dat")
flow_radius_hisq = list(numpy.around(numpy.loadtxt(inputfolder+"/s064t16_b07188_m00113_m0306/flowradius_s064t16_b07188_m00113_m0306.dat"), 6)) #need to change this to flow -times-
flow_times_hisq = [(32*i)**2 /8 for i in flow_radius_hisq] #they are the same for all these lattices
tauT_16 = list(numpy.around(numpy.loadtxt(inputfolder+"/s064t16_b07188_m00113_m0306/tauT_imps064t16_b07188_m00113_m0306.dat"), 4))
tau_16 = [i* 16 * a_hisq for i in tauT_16 ]

print(flow_radius_quenched)
print(flow_radius_hisq)


#-------------------
#set plot parameters
#-------------------
figurestyle        = dict(bbox_inches="tight")
plotstyle_points   = dict(fmt='D-', linewidth=1, markersize=4, capsize=2, mew=0.5, fillstyle='none')
legendstyle        = dict(loc="center left", frameon=True, framealpha=0.8, edgecolor='none', fancybox=False, facecolor="w", title="", labelspacing=0.1, borderpad=0.1, handletextpad=0.4)
figurestyle        = dict(bbox_inches="tight", pad_inches=0)
labelboxstyle      = dict(boxstyle="square", fc="w", ec='none', alpha=0.7, pad=0.15, zorder=999999)
plotstyle_fill     = dict(linewidth=0.5)
xlabelstyle        = dict(horizontalalignment='right', verticalalignment='bottom',  bbox=labelboxstyle, zorder=999998)
ylabelstyle        = dict(horizontalalignment='left', verticalalignment='top', rotation=0, bbox=labelboxstyle, zorder=999999)



plt.rc('text', usetex=True)
plt.rc('text.latex')
plt.rc('font', family='serif', size=9)
fig = plt.figure() 
ax = fig.add_subplot(1,1,1) 
ax.xaxis.set_label_coords(0.975,0.025)
ax.yaxis.set_label_coords(0.01,0.97)


ax.set_xlabel(r'$\tau [\mathrm{fm}]$', **xlabelstyle)
ax.set_ylabel(r'$\displaystyle \frac{G_{r_F}(\tau T)}{G_\mathrm{norm}(\tau T)}$', **ylabelstyle)


ax.set_xlim([0,0.6])
ax.set_ylim([0,12])
aspect=1/2
ax.set_aspect(1.0/ax.get_data_ratio()*aspect)
plots = []

offset=0.02
flowstart=10
flowend=21

ax.set_title(r'$r_F =0.075$', x=0.5, y=0.9, bbox=labelboxstyle, zorder=999999)

#flow radii that are the same:
#hisq 15
#quenched 20
#or
#hisq 22
#quenched 23

#PLOT LATTICE EFFECTS FOR ANIMATION
plots.append(ax.errorbar(tau_16, EE_hisq[15,:], EE_err_hisq[15,:], label=r'$1.10 \,T_C$, HISQ', **plotstyle_points, color=cm.gnuplot(0.80), zorder=-4))
plots.append(ax.errorbar(tau_32, EE_quenched[20,:], EE_err_quenched[20,:], label=r'$1.11 \,T_C$, quenched', **plotstyle_points, color=cm.gnuplot(0.4), zorder=-6))
b = ax.axvline(x=(flow_radius_hisq[15]/numpy.sqrt(8*0.014)+offset)*16*a_hisq, ymin=0, ymax=1, color=cm.gnuplot(0.8), alpha=0.8, zorder=-5, dashes=(4,4), lw=0.5)
a = ax.axvline(x=(flow_radius_quenched[20]/numpy.sqrt(8*0.014)+offset)*32*a_quenched, ymin=0, ymax=1, color=cm.gnuplot(0.4), alpha=0.8, zorder=-6, dashes=(4,4), lw=0.5)
ax.legend(handles=plots, **legendstyle, handler_map={type(plots[0]): HandlerErrorbar(xerr_size=0.4)}, handlelength=1.5)
fig.savefig(outputfolder+"/EE_hisq_quenched_1.1Tc.pdf", **figurestyle) 
