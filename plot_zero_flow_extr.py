#!/usr/local/bin/python
import math
import numpy
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from os import listdir
import sys
from latqcdtools import fitting as fit

def get_color(myarray, i, start, end):
    return cm.gnuplot((myarray[i]-myarray[start])/(myarray[end-1]-myarray[start]))

qcdtype="quenched"
outputfolder="../plots/"+qcdtype+"/"
inputfolder="../data_merged/"+qcdtype+"/"

EE_16 = numpy.loadtxt(inputfolder+"/s064t16_b0687361/EE_s064t16_b0687361_zero_flow_extr.dat")
EE_20 = numpy.loadtxt(inputfolder+"/s080t20_b0703500/EE_s080t20_b0703500_zero_flow_extr.dat")
EE_24 = numpy.loadtxt(inputfolder+"/s096t24_b0719200/EE_s096t24_b0719200_zero_flow_extr.dat")
EE_30 = numpy.loadtxt(inputfolder+"/s120t30_b0739400/EE_s120t30_b0739400_zero_flow_extr.dat")
EE_cont_2015_new = numpy.loadtxt(inputfolder+"/multi-level_2015/EE_2015_cont.txt")

tauT_16 = list(numpy.around(numpy.loadtxt(inputfolder+"/s064t16_b0687361/tauT_imps064t16_b0687361.dat"), 4))
tauT_20 = list(numpy.around(numpy.loadtxt(inputfolder+"/s080t20_b0703500/tauT_imps080t20_b0703500.dat"), 4))
tauT_24 = list(numpy.around(numpy.loadtxt(inputfolder+"/s096t24_b0719200/tauT_imps096t24_b0719200.dat"), 4))
tauT_30 = list(numpy.around(numpy.loadtxt(inputfolder+"/s120t30_b0739400/tauT_imps120t30_b0739400.dat"), 4))

#-------------------
#set plot parameters
#-------------------
figurestyle        = dict(bbox_inches="tight")
plotstyle_points   = dict(fmt='D-', linewidth=1, markersize=4, capsize=2, mew=0.5, fillstyle='none')
legendstyle        = dict(loc="center right", frameon=True, framealpha=0, title=r'$N_\tau$', prop={'size': 9})
figurestyle        = dict(bbox_inches="tight")
labelboxstyle      = dict(boxstyle="round", fc="w", ec='none', alpha=0.7)
plotstyle_fill     = dict(linewidth=0.5)
xlabelstyle        = dict(horizontalalignment='right', verticalalignment='bottom',  bbox=labelboxstyle, zorder=999999)
ylabelstyle        = dict(horizontalalignment='left', verticalalignment='top', rotation=0, bbox=labelboxstyle, zorder=999999)

plt.rc('text', usetex=True)
plt.rc('text.latex')
plt.rc('font', family='serif')
fig = plt.figure() 
ax = fig.add_subplot(1,1,1) 
ax.set_title(r'')
ax.xaxis.set_label_coords(0.975,0.050)
ax.yaxis.set_label_coords(0.025,0.975)
ax.set_xlabel(r'$\tau T$', **xlabelstyle)
ax.set_ylabel(r'$\displaystyle \frac{G_{r_F \rightarrow 0}(\tau T)}{G_\mathrm{norm}(\tau T)}$', **ylabelstyle)
#legendstyle["title"]=r"$N_\tau$"


ax.set_xlim([0.1,0.49]); 
ax.set_ylim([2,3.99])

ax.errorbar(tauT_16[3:], EE_16[:,0], EE_16[:,1], label=r'$16$', **plotstyle_points, color=cm.gnuplot(0.03))
ax.errorbar(tauT_20[3:], EE_20[:,0], EE_20[:,1], label=r'$20$', **plotstyle_points, color=cm.gnuplot(0.40))
ax.errorbar(tauT_24[3:], EE_24[:,0], EE_24[:,1], label=r'$24$', **plotstyle_points, color=cm.gnuplot(0.65))
ax.errorbar(tauT_30[3:], EE_30[:,0], EE_30[:,1], label=r'$30$', **plotstyle_points, color=cm.gnuplot(0.90))
#ax.fill_between(EE_cont_2015_new[:,0], EE_cont_2015_new[:,1]-EE_cont_2015_new[:,2], EE_cont_2015_new[:,1]+EE_cont_2015_new[:,2], facecolor='magenta', alpha=0.75, **plotstyle_fill, label="2015 ML \n cont. extr. \n new Z")
ax.legend(**legendstyle)
fig.savefig(outputfolder+"EE_flow0_extr.pdf", **figurestyle) 
