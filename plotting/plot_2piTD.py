#!/usr/bin/env python3

outputfolder = "/home/altenkort/work/correlators_flow/plots/"

import lib_process_data as lpd
import numpy as np
# import argparse


def main():

    # parser = argparse.ArgumentParser()
    # parser.add_argument('--test')
    # args = parser.parse_args()

    lpd.create_folder(outputfolder)

    # xlims = [0.9, 3.2]
    xlims = [0.77, 2.2]
    ylims = [0, 14]
    xlabel = r'$T/T_\mathrm{c}$'
    ylabel = r'$2\pi TD $'



    fig, ax, plots = lpd.create_figure(figsize="wide", xlims=xlims, ylims=ylims, xlabel=xlabel, ylabel=ylabel)

    pqcd = 8.4
    adscft = 0.9

    plots.append(ax.errorbar([1.5,2.2], [pqcd, pqcd], alpha=0.5, zorder=-10000, color='black', fmt='--', label=r'pQCD (NLO, $\alpha\approx 0.2)$'))
    plots.append(ax.errorbar(xlims, [adscft, adscft], alpha=0.5, zorder=-10000, fmt=':', color='black', label='AdS/CFT estimate'))
    plots.append(ax.errorbar(0, 0, label=' ', markersize=0, alpha=0, lw=0))

    # x y yerr
    AL_bottom = [[1.3, 1.13, 0.91], [1.5, 0.345, 0.165], [2.25, 1.995, 1.665]]
    AL_charm = [[1.5, 1.785, 0.945], [2.26, 0.855, 0.195]]

    brambilla = [[1.1, 4.45, 2.13], [1.52, 6.52, 3.07], [3, 12.83, 7.12]]

    hisq = [[196.0, 8.26, 1.94, 2.43], [220.0, 7.77, 2.16, 2.22], [251.0, 6.82, 1.12, 1.21], [352.0, 5.32, 0.26, 0.33]]

    hisq_new_low = [[220.0, 6.86, 1.00, 1.07], [251.0, 4.88, 0.46, 0.51], [295, 3.64, 0.29, 0.43], [352.0, 2.69, 0.07, 0.09]]
    hisq_new_high = [[220.0, 11.78, 0.76, 0.74], [251.0, 9.44, 0.43, 0.45], [295, 7.62, 0.32, 0.43],  [352.0, 6.17, 0.08, 0.12]]

    hisqplot =[]
    Tc = 180  # ask dibyendu

    for arr_low, arr_high in zip(hisq_new_low, hisq_new_high):
        x = arr_low[0]/Tc
        max_err_stat = np.amax((arr_low[2], arr_low[3], arr_high[2], arr_high[3]))
        err_sys = np.fabs(arr_low[1]-arr_high[1])/2
        err_tot = np.sqrt(max_err_stat**2+err_sys**2)
        y_avg = (arr_low[1]+arr_high[1])/2
        print(y_avg, err_tot, err_tot/y_avg)
        y = (4*np.pi/(y_avg-err_tot) + 4*np.pi/(y_avg+err_tot))/2
        err = y - 4*np.pi/(y_avg+err_tot)
        print(y, err, err/y)
        hisqplot.append([x, y, err])

    altenkort = [1.5, 4.42, 1.02]

    kaczmarek = [1.48, 5.31, 1.6]


    # LQCD hadronic corrs
    # plots.append(ax.errorbar([0,],[0,], fmt='.', linewidth=0, markersize=0, label=r'finite $M$, quenched'))
    # plots.append(ax.errorbar([x[0] for x in AL_bottom], [y[1] for y in AL_bottom], [yerr[2] for yerr in AL_bottom], color='tab:purple', **plotstyle, label=r'bottom, Lorenz et al. \textquotesingle 21' ))
    # plots.append(ax.errorbar([x[0] for x in AL_charm], [y[1] for y in AL_charm], [yerr[2] for yerr in AL_charm], color='tab:pink', **plotstyle, label=r'charm, Lorenz et al. \textquotesingle 21' ))
    # plots.append(ax.errorbar(0, 0, label=' ', markersize=0, alpha=0, lw=0))

    # LQCD HQ limit
    plotstyle = dict(fmt='.', markersize=0, fillstyle='none')
    plots.append(ax.errorbar([0,],[0,], fmt='.', linewidth=0, markersize=0, label=r'quenched'))  # $M\rightarrow \infty$
    plots.append(ax.errorbar(kaczmarek[0], kaczmarek[1], kaczmarek[2], **plotstyle, color='tab:gray', label=r'Francis et al. \textquotesingle 15'))
    plots.append(ax.errorbar([x[0] for x in brambilla], [y[1] for y in brambilla], [yerr[2] for yerr in brambilla], color='C0', **plotstyle,
                             label=r'Brambilla et al. \textquotesingle 20'))
    # plotstyle.update(linewidth=1.5, capsize=3, zorder=10)
    plots.append(ax.errorbar(altenkort[0], altenkort[1], altenkort[2], **plotstyle, color='C2', label=r'Altenkort et al. \textquotesingle 21'))
    # ax.axhline(y=1, **lpd.horizontallinestyle)

    plots.append(ax.errorbar(0, 0, label=' ', markersize=0, alpha=0, lw=0))
    plots.append(ax.errorbar([0,],[0,], fmt='.', linewidth=0, markersize=0, label="2+1 flavor"))  # M\\rightarrow \\infty$,
    print(np.column_stack((np.asarray([x[0] for x in hisqplot]),np.asarray([y[1] for y in hisqplot]), np.asarray([yerr[2] for yerr in hisqplot]))))
    plots.append(ax.errorbar([x[0] for x in hisqplot], [y[1] for y in hisqplot], [yerr[2] for yerr in hisqplot], color='C1', **plotstyle,
                             label="HotQCD \\textquotesingle 22"))  # (systematic errors \nunderestimated) \n (finite $a$, $\\tau_F$)

    lpd.legendstyle.update(loc="center left", bbox_to_anchor=(1, 0.5), framealpha=0)

    leg = ax.legend(handles=plots, **lpd.legendstyle, **lpd.leg_err_size(x=1, y=0.3))

    ax.tick_params(axis='y', direction='out')

    # ax.text(0.78, 0.45, 'HotQCD preliminary', transform=ax.transAxes,
    #         fontsize=9, color='C1', alpha=0.5,
    #         ha='center', va='center', rotation='0', zorder=-1000000)

    filename = outputfolder + "/2piTD.pdf"
    fig.savefig(filename)
    print("saved correlator plot", filename)

    return


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()