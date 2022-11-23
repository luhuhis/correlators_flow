#!/usr/bin/env python3

import lib_process_data as lpd
import sys
import numpy
import matplotlib.pyplot as plt

# input:
inputfolder = "../../data/merged/quenched_pertLO_wilsonFlow/EE/"
#   inputfolder+"EE_latt_flow_"+str(Ntau)+".dat"
#   inputfolder+"EE_cont_flow_"+str(Ntau)+".dat"

# output:
outputfolder = "./"
#   outputfolder+"/EE_pert_contvslatt_flow.pdf"

lpd.create_folder(outputfolder)

try:
    Ntau = int(sys.argv[1])
except IndexError:
    exit("Invalid Arguments. Usage: script.py <Ntau>, e.g. script.py 20")


# -------------------
# set plot parameters
# -------------------
figurestyle = dict(bbox_inches="tight")
plotstyle_points = dict(linewidth=1, markersize=6, capsize=2, mew=1, fillstyle='none')
labelboxstyle = dict(boxstyle="square", fc="w", ec='none', alpha=0.8, pad=0.15, zorder=999999)
legendstyle = dict(loc="upper right", frameon=True, framealpha=0.8, edgecolor='none', fancybox=False, facecolor="w", labelspacing=0.1, borderpad=0.1,
                   handletextpad=0.4)  # , handlelength=1)#, prop={'size': 9})
figurestyle = dict(bbox_inches="tight", pad_inches=0.03)
plotstyle_fill = dict(linewidth=0.5)
xlabelstyle = dict(horizontalalignment='right', verticalalignment='bottom', bbox=labelboxstyle, zorder=999998)
ylabelstyle = dict(horizontalalignment='left', verticalalignment='top', rotation=0, bbox=labelboxstyle, zorder=999999)
titlestyle = dict(x=0.5, y=0.9, bbox=labelboxstyle, verticalalignment='top', zorder=999999)
flowlimitplotstyle = dict(fmt='|', mew=1, markersize=25)

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=10)

# fig = matplotlib.pyplot.figure(figsize=figsize, constrained_layout=True)
fig = plt.figure(figsize=(3, 2.5), constrained_layout=True)
ax = fig.add_subplot(1, 1, 1)
ax.set_yscale('log')
ax.minorticks_off()
ax.tick_params(direction='in', pad=1.5)

# copy these from the corresponding mathematica notebook!
flowresolution = 1 / 1000
maxflowradius = 150 / 1000
contCorrResolution = 10

nflowtimes = len(numpy.arange(0, maxflowradius + flowresolution / 2, flowresolution))

ax.xaxis.set_label_coords(0.98, 0.03)
ax.yaxis.set_label_coords(0.2, 0.98)
ax.set_ylabel(r'$\displaystyle \frac{G^\mathrm{LO}}{(g^2 C_F)T^4}$', **ylabelstyle)
ax.set_xlabel(r'$\tau T$', **xlabelstyle)
ax.set_ylim([0.5, 50000])
ax.set_xlim([0, 0.52])
plots = []

EE_latt = numpy.loadtxt(inputfolder + "EE_latt_flow_" + str(Ntau) + ".dat")
EE_cont = numpy.loadtxt(inputfolder + "EE_cont_flow_" + str(Ntau) + ".dat")
ntauTcont = len(numpy.arange(1 / Ntau / contCorrResolution, Ntau / 2 + 1 / Ntau / contCorrResolution / 2, 1 / Ntau / contCorrResolution))

for i, fmtlat, fmtcont, color in zip((0, 50, 100), ('x', 'x', '+'), ('-', '--', ':'), ('k', 'C0', 'C1')):
    i0 = int(i * (Ntau / 2))
    i1 = int((i + 1) * (Ntau / 2))
    j0 = int(i * ntauTcont)
    j1 = int((i + 1) * ntauTcont)
    flow_radius = EE_latt[i0][1]

    flowstr = '{:.2f}'.format(flow_radius)
    plots.append(
        ax.errorbar([row[0] for row in EE_cont[j0:j1]], [row[2] for row in EE_cont[j0:j1]], fmt=fmtcont, color=color, **plotstyle_points, label=flowstr,
                    zorder=-1))
    if i == 10:
        plots.append(
            ax.errorbar([row[0] for row in EE_latt[i0:i1]], [row[2] for row in EE_latt[i0:i1]], fmt=fmtlat, color=color, **plotstyle_points, label=flowstr,
                        zorder=-2))
    xval = numpy.abs(flow_radius / numpy.sqrt(8 * 0.014) - EE_cont[j0:j1, 0]).argmin()

    if i != 0:
        ax.errorbar(EE_cont[j0 + xval][0], EE_cont[j0 + xval][2], **flowlimitplotstyle, color=color, alpha=1, zorder=-1000)
ax.legend(handles=plots, title=r'$\sqrt{8\tau_\mathrm{F}} T=$', **legendstyle)
# plt.tight_layout(0)
fig.savefig(outputfolder + "/EE_pert_contvslatt_flow.pdf")  # , **figurestyle)
print("saved pert corr plot", outputfolder + "/EE_pert_contvslatt_flow.pdf")
ax.lines.clear()
ax.collections.clear()
plots.clear()

# for i in range(0,nflowtimes):
# i0=int(i*(Ntau/2))
# i1=int((i+1)*(Ntau/2))
# j0=int(i*ntauTcont)
# j1=int((i+1)*ntauTcont)
# flow_radius=EE_latt[i0][1]

# ax.set_title(r'$\sqrt{8\tau_\mathrm{F}} T=$ '+'{0:.3f}'.format(flow_radius), x=0.5, y=0.96, bbox=labelboxstyle, verticalalignment='top', zorder=999999)
# plots.append(ax.errorbar([row[0] for row in EE_cont[j0:j1]], [row[2] for row in EE_cont[j0:j1]], fmt='-', **plotstyle_points, color="orange", label="LO cont.", zorder=-1))
# plots.append(ax.errorbar([row[0] for row in EE_latt[i0:i1]], [row[2] for row in EE_latt[i0:i1]], fmt='x', **plotstyle_points, color="blue", label=r'LO latt. $N_\tau=$'+str(Ntau), zorder=-2))
# if i != 0:
# ax.axvline(x=(flow_radius/numpy.sqrt(8*0.014)), ymin=0, ymax=1, color='grey', alpha=0.8, dashes=(4,4), lw=0.5, zorder=-10000)
# ax.legend(handles=plots, **legendstyle)
##fig.savefig(outputfolder+"Ntau_"+str(Ntau)+"/EE_free_contvslatt_flow_"+'{0:.4f}'.format(flow_radius)+".png", dpi=300, **figurestyle)
# fig.savefig(outputfolder+"Ntau_"+str(Ntau)+"/EE_free_contvslatt_flow_"+str(i+1)+".png", dpi=300, **figurestyle)
# ax.lines.clear() ; ax.collections.clear() ; plots.clear()

# plot relative
# for i in range(nflowtimes):
# i0=int(i*(Ntau/2))
# i1=int((i+1)*(Ntau/2))
# j0=int(i*ntauTcont)
# j1=int((i+1)*ntauTcont)
# flow_radius=EE_latt[i0][1]

# ax.set_yscale('linear')
# ax.set_ylim([-1,4])
# ax.set_xlim([-0.01,0.51])

# ax.set_title(r'$\sqrt{8\tau_\mathrm{F}} T=$ '+'{0:.3f}'.format(flow_radius), x=0.5, y=0.96, bbox=labelboxstyle, verticalalignment='top', zorder=999999)
# xdata_lat = [row[0] for row in EE_latt[i0:i1]]
# ydata_lat = [row[2] for row in EE_latt[i0:i1]]
# xdata_cont = [row[0] for row in EE_cont[j0:j1]]
# ydata_cont = [row[2] for row in EE_cont[j0:j1]]
# for k,tauT in enumerate(xdata_lat):
# cont_index = xdata_cont.index(tauT)
# ydata_lat[k] /= ydata_cont[cont_index]

# plots.append(ax.errorbar(xdata_lat, ydata_lat, fmt='x-', **plotstyle_points, color="orange", label="LO: ratio latt/cont.", zorder=-1))
# if i != 0:
# ax.axvline(x=(flow_radius/numpy.sqrt(8*0.014)), ymin=0, ymax=1, color='grey', alpha=0.8, dashes=(4,4), lw=0.5, zorder=-10000)
# ax.legend(handles=plots, **legendstyle)
# print(i)
# fig.savefig(outputfolder+"/relative/Ntau_"+str(Ntau)+"/EE_free_contvslatt_flow_"+str(i+1)+".pdf", **figurestyle)
# ax.lines.clear() ; ax.collections.clear() ; plots.clear()
