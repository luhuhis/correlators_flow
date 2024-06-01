#!/usr/bin/env python3

import lib_process_data as lpd
import numpy
import argparse
import pandas as pd


def plot_fig_tauT_on_xaxis(args):
    fig, ax, axtwiny = lpd.create_figure(ylabelpos=(0.15, 0.98), ylims=[0.5, 50000], xlims=[0, 0.52], ylabel=r'$\displaystyle \frac{G^\mathrm{LO}}{(g^2 C_F)T^4}$', xlabel=r'$\tau T$')
    ax.set_yscale('log')
    ax.minorticks_off()
    axtwiny.set_yscale('log')
    axtwiny.minorticks_off()


    # copy these from the corresponding mathematica notebook!
    contCorrResolution = 10

    ax.set_xticks((0.0, 0.1, 0.2, 0.3, 0.4, 0.5))
    plots = []

    EE_latt = numpy.loadtxt(args.inputfolder + "EE_latt_flow_" + str(args.Ntau) + ".dat")
    EE_cont = numpy.loadtxt(args.inputfolder + "EE_cont_flow_" + str(args.Ntau) + ".dat")
    ntauTcont = len(numpy.arange(1 / args.Ntau / contCorrResolution, args.Ntau / 2 + 1 / args.Ntau / contCorrResolution / 2, 1 / args.Ntau / contCorrResolution))

    plotstyle_points = dict(linewidth=1, markersize=6, fillstyle='none')
    for i, fmtlat, fmtcont, color in zip((0, 50, 100), ('x', '+', 'o'), ('-', '--', ':'), ('k', 'C0', 'C1')):
        i0 = int(i * (args.Ntau / 2))
        i1 = int((i + 1) * (args.Ntau / 2))
        j0 = int(i * ntauTcont)
        j1 = int((i + 1) * ntauTcont)
        flow_radius = EE_latt[i0][1]

        flowstr = '{:.2f}'.format(flow_radius)

        plots.append(
            ax.errorbar([row[0] for row in EE_cont[j0:j1]], [row[2] for row in EE_cont[j0:j1]], fmt=fmtcont, color=color, **plotstyle_points, label=flowstr,
                        zorder=-3)
            # ax.errorbar([row[0] for row in EE_cont[j0:j1]], [Gcont(row[0], flow_radius) for row in EE_cont[j0:j1]], fmt=fmtcont, color=color, **plotstyle_points, label=flowstr,
            #             zorder=-3)
        )
        plots.append(
                ax.errorbar([row[0] for row in EE_latt[i0:i1]], [row[2] for row in EE_latt[i0:i1]], fmt=fmtlat, color=color, **plotstyle_points, label=flowstr,
                            zorder=-2))
        xval = numpy.abs(flow_radius / numpy.sqrt(8 * 0.014) - EE_cont[j0:j1, 0]).argmin()

        if i != 0:
            ax.errorbar(EE_cont[j0 + xval][0], EE_cont[j0 + xval][2], fmt='|', mew=1, markersize=25, color=color, alpha=1, zorder=-1)
    ax.legend(handles=plots, title=r'$\sqrt{8\tau_\mathrm{F}} T$', loc="upper right", frameon=True, framealpha=0.8, edgecolor='none', fancybox=False, facecolor="w", labelspacing=0.1, borderpad=0.1,
                       handletextpad=0.4, handlelength=1.25)
    lpd.create_folder(args.outputfolder)
    file = args.outputfolder + "/EE_pert_contvslatt_flow.pdf"
    fig.savefig(file)
    print("saved pert corr plot", file)

def plot_fig_flow_on_xaxis(args):

    fig, ax, axtwiny = lpd.create_figure(ylabelpos=(0.17, 0.98), xlabelpos=(0.96,0.01), ylims=[0.5, 50000], xlims=[0, 0.15], ylabel=R'$\displaystyle \frac{G^\mathrm{LO}}{(g^2 C_F)T^4}$', xlabel=R'$\sqrt{8\tau_\mathrm{F}} T$', ylabelbox=None)
    ax.set_yscale('log')
    ax.minorticks_off()
    axtwiny.set_yscale('log')
    axtwiny.minorticks_off()

    plots = []
    zorder = 0
    plot_tauTs = ["0.041667", "0.125", "0.25", "0.375", "0.5"]
    colors = [f"C{i}" for i in range(len(plot_tauTs))]

    for prefix, fmt, subtitle in zip(["EE_cont_flow_", "EE_latt_flow_"], ["-", "--"], ["cont.", "latt."]):

        plots.append(ax.errorbar(0,0, label=subtitle, ms=0, lw=0))

        # Load the data from the text file
        file_path = f"{args.inputfolder}{prefix}{args.Ntau}.dat"
        data = pd.read_csv(file_path, delim_whitespace=True, header=None, names=['x', 'y', 'value'])

        # Pivot the data to get it in the right format for plotting
        pivot_data = data.pivot(index='x', columns='y', values='value')

        i = 0
        for key in pivot_data.index:
            if str(key) in plot_tauTs:
                plots.append(ax.errorbar(pivot_data.columns, pivot_data.loc[key], label=f'{round(key,3)}', fmt=fmt, zorder=zorder, color=colors[i]))
                zorder -= 1
                i += 1
        


    ax.legend(handles=plots, title=R'$\tau T$', loc="upper right", frameon=False, edgecolor='none', fancybox=False, facecolor="w", labelspacing=0.1, borderpad=0.1,
                       handletextpad=0.4, handlelength=1, ncol=2, columnspacing=0.5)
    lpd.create_folder(args.outputfolder)
    file = args.outputfolder + "/EE_pert_contvslatt_flow_dep.pdf"
    fig.savefig(file)
    print("saved pert corr plot", file)


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--inputfolder', type=str)
    parser.add_argument('--outputfolder', type=str)
    parser.add_argument('--Ntau', required=True, type=int)
    # parser.add_argument('--swap_xaxis_and_leg', action="store_true", type=bool, help="Put flow time on the x-axis and tauT in the legend.")
    args = parser.parse_args()

    # plot_fig_tauT_on_xaxis(args)
    plot_fig_flow_on_xaxis(args)


    return


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()