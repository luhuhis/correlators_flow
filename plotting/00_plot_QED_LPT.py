#!/usr/bin/env python3
import lib_process_data as lpd

import numpy
import matplotlib

inputfolder="../../data/merged/"
outputfolder="../../plots/"
lpd.create_folder(outputfolder)


xdata,I2,I3 = numpy.loadtxt(inputfolder+"EE_QED_LPT.dat", unpack=True)

#xlims=[0,0.7584],                                  ylims=[-1.62,0.0],
fig, ax, plots = lpd.create_figure(xlabel=r'$\sqrt{8\tau_\mathrm{F}}/a$',
                                  xlabelpos=(0.93,0.06),
                                  ylabel=r'', ylabelpos=(0, 0.95), figsize=(3.5, 3))

ax.errorbar(numpy.sqrt(xdata*8), I2/I2[0]*100, lw=0.75, label=r'operator mixing term')  # $-a^4\frac{3}{2}\mathcal{I}_2(\tau_\mathrm{F})$
ax.errorbar(numpy.sqrt(xdata*8), I3/I3[-1]*100, lw=0.75, label=r'flow time renorm.\ term')  # $a^4\frac{1}{4}\mathcal{I}_3(\tau_\mathrm{F})$
ax.axhline(y=100, **lpd.horizontallinestyle)
ax.axhline(y=0, **lpd.horizontallinestyle)

import matplotlib.ticker as mtick
ax.yaxis.set_major_formatter(mtick.PercentFormatter())

lpd.legendstyle.update(dict(loc="center right", bbox_to_anchor=(1, 0.5)))
ax.legend(**lpd.legendstyle, title="QED: lattice artifacts at NLO", fontsize=8, title_fontsize=8)
fig.savefig(outputfolder+"/EE_QED_LPT.pdf")
print("saved QED LPT plot", outputfolder+"/EE_QED_LPT.pdf")
