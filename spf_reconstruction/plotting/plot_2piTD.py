#!/usr/bin/env python3

outputfolder = "/home/altenkort/work/correlators_flow/plots/"

import lib_process_data as lpd
import numpy as np
import argparse
import rundec


def D2piT(mu):
    crd = rundec.CRunDec()
    Alphas = crd.AlphasLam(0.339, mu, 3, 5)
    g2 = 4. * np.pi * Alphas
    g = np.sqrt(g2)
    kappa_NLO_by_T3 = (16*np.pi)/3 * Alphas**2 * (np.log(1/g) + 0.07428 + 1.9026 * g)
    D2piT_NLO = 4 * np.pi / kappa_NLO_by_T3
    return D2piT_NLO


def mean(low, high):
    return (low+high)/2


def err(low, high):
    return high - (low+high)/2


def mean_and_err(low, high):
    return mean(low,high), err(low,high)


def mean_and_err_kappaByT3_to_2piTD(low, high):
    D2piT_low = 4 * np.pi / high
    D2piT_high = 4 * np.pi / low
    return mean(D2piT_low, D2piT_high), err(D2piT_low, D2piT_high)


def thiscolor(i):
    return lpd.get_discrete_color(i)


def plot_EE_quenched_literature(plots, ax, i):

    plot = ax.errorbar(0, 0, fmt='.', markersize=0, label=r'Quenched QCD')  # via $\kappa_E$
    plots.append(plot)

    ms=0
    myfmt='.'

    #  2015, A. Francis, O. Kaczmarek, M. Laine, T. Neuhaus, and H. Ohno, Phys. Rev. D 92, 116003
    plot = ax.errorbar(1.46, *mean_and_err_kappaByT3_to_2piTD(2.6-0.8,2.6+0.8),
                color=thiscolor(3-i), fmt=myfmt, markersize=ms, label=r'Francis \textquotesingle 15 (ML)', zorder=5)
    i += 1
    plots.append(plot)

    #  2020, Nora Brambilla, Viljami Leino, Peter Petreczky, and Antonio Vairo, Phys. Rev. D 102, 074503
    color = thiscolor(3-i)
    plot = ax.errorbar(1.1, *mean_and_err_kappaByT3_to_2piTD(1.91,5.40), markersize=ms, fmt=myfmt, color=color, label=r'TUMQCD \textquotesingle 20 (ML)')
    ax.errorbar(1.48, *mean_and_err_kappaByT3_to_2piTD(1.31,3.64), markersize=ms, fmt=myfmt, color=color)
    ax.errorbar(3,   *mean_and_err_kappaByT3_to_2piTD(0.63,2.20), markersize=ms, fmt=myfmt, color=color, zorder=-10)
    i += 1
    plots.append(plot)

    # 2021 Altenkort
    plot = ax.errorbar(1.5, *mean_and_err_kappaByT3_to_2piTD(1.5, 3.2),
                markersize=ms, fmt=myfmt, color='k', label=r'\textbf{this work}* (flow)', zorder=10)
    plots.append(plot)

    #  2022, Nora Brambilla, Viljami Leino, Julian Mayer-Steudte, Peter Petreczky, 	arXiv:2206.02861
    color = thiscolor(3-i)
    plot = ax.errorbar(1.52, *mean_and_err_kappaByT3_to_2piTD(1.70,3.12), zorder=5,  markersize=ms, fmt=myfmt, color=color, label=r'TUMQCD \textquotesingle 22 (flow)')
    i += 1
    plots.append(plot)

    # 2022, Debasish Banerjee, Rajiv Gavai, Saumen Datta, Pushan Majumdar, arXiv:2206.15471v1
    color = thiscolor(3-i)
    plot = ax.errorbar(1.2, *mean_and_err_kappaByT3_to_2piTD(2.1,3.5), markersize=ms, fmt=myfmt, color=color,
                label=r'Banerjee \textquotesingle 22 (ML)')
    ax.errorbar(1.54, *mean_and_err_kappaByT3_to_2piTD(1.5,2.8), markersize=ms, fmt=myfmt, color=color)
    ax.errorbar(2.0, *mean_and_err_kappaByT3_to_2piTD(1.0,2.3), markersize=ms, fmt=myfmt, color=color)
    ax.errorbar(2.5, *mean_and_err_kappaByT3_to_2piTD(0.9,2.1), markersize=ms, fmt=myfmt, color=color)
    ax.errorbar(3.0, *mean_and_err_kappaByT3_to_2piTD(0.8,1.8), markersize=ms, fmt=myfmt, color=color)
    ax.errorbar(3.5, *mean_and_err_kappaByT3_to_2piTD(0.75, 1.5), markersize=ms, fmt=myfmt, color=color)
    i += 1
    plots.append(plot)

    return i


def plot_PT_and_AdsCFT(plots, ax):

    adscft = 0.9
    plot = ax.errorbar([0, 1.5], [adscft, adscft], zorder=-10000, lw=0.75, fmt='--', color='black', label='AdS/CFT estimate')
    plots.append(plot)

    Tc = 0.156  # GeV
    xpoints = np.linspace(1.5, 3.2, 100)
    pqcd_2piT = []
    pqcd_4piT = []
    for x in xpoints:
        pqcd_2piT.append(D2piT(2 * np.pi * x * Tc))
        pqcd_4piT.append(D2piT(4 * np.pi * x * Tc))

    plot = ax.fill_between(xpoints, pqcd_2piT, pqcd_4piT, zorder=-10000, color='palegreen', ec=None,
                    label=r'\begin{flushleft}pQCD NLO, \newline $\mu \in [2\pi T, 4\pi T]$\end{flushleft}')  # via $\kappa_E$
    plots.append(plot)

    # ax.errorbar(xpoints, pqcd_4piT, zorder=-100, lw=0.5, fmt='-',  color='palegreen')
    # ax.errorbar(xpoints, pqcd_2piT, zorder=-100, lw=0.5, fmt='-', color='palegreen')





def parse_args():
    parser = argparse.ArgumentParser()
    # parser.add_argument('--add_leg_titles', action="store_true")
    args = parser.parse_args()
    return args


def plot_T_matrix(plots, ax):
    Tc = 156
    T = np.asarray([194, 258, 320, 400])
    D2piT_chi1 = [2.309325630638798, 3.559583919794666, 5.49821428630848, 9.944405023982574]
    D2piT_chi06 = [1.654568323280973, 2.5372613158336925, 3.792920321249195, 7.078861688212292]

    plot = ax.fill_between(T/Tc, D2piT_chi06, D2piT_chi1, zorder=-100, color='lemonchiffon', ec=None,
                    label=r'T-matrix')
    plots.append(plot)


def plot_Bayesian(plots, ax):
    Ts = np.array([1.00, 1.50, 2.00, 2.48])
    low = [1.94, 3.31, 4.81, 6.38]
    high = [2.75, 4.63, 7.06, 10.06]
    plot = ax.fill_between(Ts, low, high, zorder=-500, color='pink', ec=None,
                    label=r'Bayesian')
    plots.append(plot)

    # ax.errorbar(Ts, low, zorder=-100, lw=0.5, fmt='-',  color='pink')
    # ax.errorbar(Ts, high, zorder=-50, lw=0.5, fmt='-', color='pink')


def plot_hisq(plots, ax):
    plot = ax.errorbar(0, 0, fmt='.', markersize=0, label=r'2+1-flavor QCD') # via $\kappa_E$
    plots.append(plot)
    Tc = 180

    data = [[195, 8.5, 13.4],
            [220, 6.0, 10.8],
            [251, 4.7, 9.0],
            [293, 3.9, 7.7]]
    for i, dat in enumerate(data):
        T = dat[0]
        kappa_low = dat[1]
        kappa_high = dat[2]
        if i == 0:
            label = r'\textbf{this work} (flow)'
        else:
            label = None

        plot = ax.errorbar(T/Tc, *mean_and_err_kappaByT3_to_2piTD(kappa_low, kappa_high), markersize=0, fmt='.', color='m', label=label)

        if i == 0:
            plots.append(plot)


def add_empty_line(plots, ax):
    plot = ax.errorbar(0, 0, label=' ', markersize=0, alpha=0, lw=0)
    plots.append(plot)


def plot_ALICE(plots, ax):
    plot = ax.errorbar(1, *mean_and_err(1.5, 4.5), fmt='.', ms=0, color='mediumaquamarine', label='ALICE')
    plots.append(plot)


def plot_fig(args):

    plots = []

    # xlims = [0.9, 3.2]
    xlims = [0.9, 2.6]
    ylims = [0, 14.25]
    xlabel = r'$T/T_\mathrm{c}$'
    ylabel = r'$2\pi TD  $'

    fig, ax, axtwiny = lpd.create_figure(figsize="wide", xlims=xlims, ylims=ylims, xlabel=xlabel, ylabel=ylabel)
    ax.set_yticks(np.arange(0, ylims[1], 2.5))
    axtwiny.set_yticks(np.arange(0, ylims[1], 2.5))
    ax.set_xticks(np.arange(1, xlims[1], 0.5))

    plot_ALICE(plots, ax)
    plot_T_matrix(plots, ax)
    plot_Bayesian(plots, ax)

    add_empty_line(plots, ax)
    plot_PT_and_AdsCFT(plots, ax)
    add_empty_line(plots, ax)

    plot_EE_quenched_literature(plots, ax, 0)
    add_empty_line(plots, ax)

    plot_hisq(plots, ax)




    # LQCD hadronic corrs
    # AL_bottom = [[1.3, 1.13, 0.91], [1.5, 0.345, 0.165], [2.25, 1.995, 1.665]]
    # AL_charm = [[1.5, 1.785, 0.945], [2.26, 0.855, 0.195]]
    # ax.errorbar([0,],[0,], fmt='.', linewidth=0, markersize=0, label=r'finite $M$, quenched')
    # ax.errorbar([x[0] for x in AL_bottom], [y[1] for y in AL_bottom], [yerr[2] for yerr in AL_bottom],
    #             fmt='.', markersize=0, color='tab:purple', label=r'bottom, Lorenz et al. \textquotesingle 21' )
    # ax.errorbar([x[0] for x in AL_charm], [y[1] for y in AL_charm], [yerr[2] for yerr in AL_charm],
    #             fmt='.', markersize=0, color='tab:pink', label=r'charm, Lorenz et al. \textquotesingle 21' )
    # ax.errorbar(0, 0, label=' ', markersize=0, alpha=0, lw=0)

    handles, labels = [], []
    for plot in plots:
        labels.append(plot.get_label())

    ax.legend(plots, labels, loc="center left", bbox_to_anchor=(1, 0.5), framealpha=0, **lpd.leg_err_size(x=1, y=0.3), fontsize=9, handlelength=1)

    return fig


def main():

    args = parse_args()

    lpd.create_folder(outputfolder)

    fig = plot_fig(args)

    filename = outputfolder + "/2piTD.pdf"
    fig.savefig(filename)
    print("saved correlator plot", filename)

    return


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()