#!/usr/local/bin/python3

#input: 
inputfolder="../../data_merged/quenched/pert_corr/"
#   inputfolder+"EE_latt_flow_"+str(Ntau)+".dat"
#   inputfolder+"EE_cont_flow_"+str(Ntau)+".dat"

#output:
outputfolder="../../plots/quenched/pert_corr/"
#   outputfolder+"/EE_pert_contvslatt_flow.pdf"

import lib_process_data as lpd
lpd.create_folder(outputfolder)

import math
import numpy
import matplotlib

xdata,I2,I3 = numpy.loadtxt(inputfolder+"EE_QED_LPT.dat", unpack=True)

fig, ax, plots = lpd.create_figure(xlims=[0,0.7584], 
                                  ylims=[-1.62,0.0],
                                  xlabel=r'$\tau_\mathrm{F}/a^2$', 
                                  xlabelpos=(0.93,0.06),
                                  ylabel=r'')

ax.errorbar(xdata, I2, lw=0.75, label=r'$-a^4\frac{3}{2}\mathcal{I}_2(\tau_\mathrm{F})$')
ax.errorbar(xdata, I3, lw=0.75, fmt='--', label=r'$a^4\frac{1}{4}\mathcal{I}_3(\tau_\mathrm{F})$')
ax.axhline(**lpd.horizontallinestyle)

lpd.legendstyle.update(dict(loc="center right", bbox_to_anchor=(1,0.5)))
ax.legend(**lpd.legendstyle)
fig.canvas.draw()
matplotlib.pyplot.tight_layout(0)
fig.savefig(outputfolder+"/EE_QED_LPT.pdf") 
print("saved QED LPT plot", outputfolder+"/EE_QED_LPT.pdf")
