
inputfolder="/work/home/altenkort/2piTD/"
outputfolder="/work/home/altenkort/2piTD/"

import matplotlib
import lib_process_data as lpd
lpd.create_folder(outputfolder) 


xlims = [1,2.5]
ylims = [0,10]
xlabel = r'$T/T_\mathrm{c}$'
ylabel = r'$2\pi TD $'
xlabelpos = None
ylabelpos = [0.02,0.97]
figsize = ((3+3/8),(3+3/8)-1/2.54)

plotstyle = dict(fmt='.', linewidth=1, markersize=0, capsize=2, mew=1, fillstyle='none')

fig, ax, plots = lpd.create_figure(xlims=xlims, ylims=ylims, xlabel=xlabel, xlabelpos=xlabelpos, ylabel=ylabel, ylabelpos=ylabelpos, UseTex=True, figsize=figsize)


pqcd=8.4
adscft=1.0

#plots.append(ax.errorbar(xlims, [pqcd,pqcd], linewidth=0.75, color='black', fmt='--', label=r'pQCD (NLO, $\alpha\approx 0.2)$'))
#plots.append(ax.errorbar(xlims, [adscft,adscft], linewidth=0.75, fmt='-.', color='black', label='AdS/CFT'))
plots.append(ax.errorbar([0,],[0,], fmt='.', linewidth=0, markersize=0, label=r'LQCD:'))

# x y yerr
AL_bottom = [[1.3, 1.13, 0.91] , [1.5, 0.345, 0.165], [2.25, 1.995, 1.665]]
AL_charm = [[1.5, 1.785, 0.945],[2.26, 0.855, 0.195]]

brambilla = [[1.1, 4.45, 2.13],[1.52, 6.52, 3.07],[3, 12.83, 7.12]]


altenkort = [1.5,4.42,1.02]

kaczmarek = [1.48,5.31,1.6]

#LQCD hadronic corrs
#plots.append(ax.errorbar([x[0] for x in AL_bottom], [y[1] for y in AL_bottom], [yerr[2] for yerr in AL_bottom], color='darkgreen', **plotstyle, label=r'bottom, Lorenz et al. \textquotesingle 21' ))
#plots.append(ax.errorbar([x[0] for x in AL_charm], [y[1] for y in AL_charm], [yerr[2] for yerr in AL_charm], color='lightgreen', **plotstyle, label=r'charm, Lorenz et al. \textquotesingle 21' ))

#LQCD HQ limit
plots.append(ax.errorbar(kaczmarek[0], kaczmarek[1], kaczmarek[2], **plotstyle, color='lightblue', label=r'$\kappa$, Francis et al. \textquotesingle 15' ))
plots.append(ax.errorbar([x[0] for x in brambilla], [y[1] for y in brambilla], [yerr[2] for yerr in brambilla], color='#2B6BB5', **plotstyle,label=r'$\kappa$, Brambilla et al. \textquotesingle 20' ))
plotstyle.update(linewidth=1.5, capsize=3, zorder=10)
plots.append(ax.errorbar(altenkort[0], altenkort[1], altenkort[2], **plotstyle, color='darkblue', label=r'$\kappa$, Altenkort et al. \textquotesingle 21'))

lpd.legendstyle.update(loc="center right", bbox_to_anchor=(1,0.6))

leg = ax.legend(handles=plots, **lpd.legendstyle)
matplotlib.pyplot.tight_layout(0)
filename = outputfolder+"/2piTD.pdf"
fig.savefig(filename) 
print("saved correlator plot", filename)
ax.lines.clear() ; ax.collections.clear() ; plots.clear()
