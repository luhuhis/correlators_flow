inputfolder = "/work/home/altenkort/2piTD/"
outputfolder = "/work/home/altenkort/2piTD/"

import lib_process_data as lpd
import numpy as np

lpd.create_folder(outputfolder)

xlims = [0.9, 3.2]
ylims = [0, 20]
xlabel = r'$T/T_\mathrm{c}$'
ylabel = r'$2\pi TD $'
xlabelpos = None
ylabelpos = [0.02, 0.97]
figsize = ((3 + 3 / 8), (3 + 3 / 8) - 1 / 2.54)

plotstyle = dict(fmt='.', linewidth=1, markersize=0, capsize=2, mew=1, fillstyle='none')

fig, ax, plots = lpd.create_figure(xlims=xlims, ylims=ylims, xlabel=xlabel, xlabelpos=xlabelpos, ylabel=ylabel, ylabelpos=ylabelpos, UseTex=True,
                                   figsize=figsize)

pqcd = 8.4
adscft = 1.0

plots.append(ax.errorbar(xlims, [pqcd, pqcd], linewidth=0.75, color='black', fmt='--', label=r'pQCD (NLO, $\alpha\approx 0.2)$'))
plots.append(ax.errorbar(xlims, [adscft, adscft], linewidth=0.75, fmt=':', color='black', label='AdS/CFT'))
plots.append(ax.errorbar(0, 0, label=' ', markersize=0, alpha=0, lw=0))

# x y yerr
AL_bottom = [[1.3, 1.13, 0.91], [1.5, 0.345, 0.165], [2.25, 1.995, 1.665]]
AL_charm = [[1.5, 1.785, 0.945], [2.26, 0.855, 0.195]]

brambilla = [[1.1, 4.45, 2.13], [1.52, 6.52, 3.07], [3, 12.83, 7.12]]

hisq = [ [196.0, 8.26, 1.94, 2.43], [220.0, 7.77, 2.16, 2.22], [251.0, 6.82, 1.12, 1.21], [352.0, 5.32, 0.26, 0.33]]

hisqplot =[]
Tc = 180

for arr in hisq:
    x = arr[0]/Tc
    max_err = np.fmax(arr[2], arr[3])
    y = (4*np.pi/(arr[1]-max_err) + 4*np.pi/(arr[1]+max_err))/2
    err = y - 4*np.pi/(arr[1]+max_err)
    hisqplot.append([x,y,err])

altenkort = [1.5, 4.42, 1.02]

kaczmarek = [1.48, 5.31, 1.6]


# LQCD hadronic corrs
# plots.append(ax.errorbar([0,],[0,], fmt='.', linewidth=0, markersize=0, label=r'finite $M$, quenched'))
# plots.append(ax.errorbar([x[0] for x in AL_bottom], [y[1] for y in AL_bottom], [yerr[2] for yerr in AL_bottom], color='tab:purple', **plotstyle, label=r'bottom, Lorenz et al. \textquotesingle 21' ))
# plots.append(ax.errorbar([x[0] for x in AL_charm], [y[1] for y in AL_charm], [yerr[2] for yerr in AL_charm], color='tab:pink', **plotstyle, label=r'charm, Lorenz et al. \textquotesingle 21' ))
# plots.append(ax.errorbar(0, 0, label=' ', markersize=0, alpha=0, lw=0))

# LQCD HQ limit
plots.append(ax.errorbar([0,],[0,], fmt='.', linewidth=0, markersize=0, label=r'$M\rightarrow \infty$, quenched'))
plots.append(ax.errorbar(kaczmarek[0], kaczmarek[1], kaczmarek[2], **plotstyle, color='C0', label=r'Francis et al. \textquotesingle 15'))
plots.append(ax.errorbar([x[0] for x in brambilla], [y[1] for y in brambilla], [yerr[2] for yerr in brambilla], color='C1', **plotstyle,
                         label=r'Brambilla et al. \textquotesingle 20'))
plotstyle.update(linewidth=1.5, capsize=3, zorder=10)
plots.append(ax.errorbar(altenkort[0], altenkort[1], altenkort[2], **plotstyle, color='C2', label=r'Altenkort et al. \textquotesingle 21'))
# ax.axhline(y=1, **lpd.horizontallinestyle)

plots.append(ax.errorbar(0, 0, label=' ', markersize=0, alpha=0, lw=0))
plots.append(ax.errorbar([0,],[0,], fmt='.', linewidth=0, markersize=0, label="$M\\rightarrow \\infty$, 2+1 flavor"))
plots.append(ax.errorbar([x[0] for x in hisqplot], [y[1] for y in hisqplot], [yerr[2] for yerr in hisqplot], color='tab:gray', **plotstyle,
                         label="PRELIMINARY \n (systematic errors missing)"))

lpd.legendstyle.update(loc="upper left", bbox_to_anchor=(0.3, 0.99), fontsize=8, framealpha=0.99)





leg = ax.legend(handles=plots, **lpd.legendstyle)

ax.text(0.61, 0.28, 'HotQCD preliminary', transform=ax.transAxes,
        fontsize=14, color='gray', alpha=0.5,
        ha='center', va='center', rotation='20', zorder=-1000000)

filename = outputfolder + "/2piTD.pdf"
fig.savefig(filename)
print("saved correlator plot", filename)
