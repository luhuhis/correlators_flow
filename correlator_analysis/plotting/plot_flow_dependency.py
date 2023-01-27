#!/usr/bin/env python3
import lib_process_data as lpd
import numpy
import argparse
from _2_plot_lateffects import apply_tree_level_imp


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--qcdtype', required=True, type=str)
    parser.add_argument('--conftype', type=str)
    parser.add_argument('--corr', choices=["EE", "BB"], required=True)
    parser.add_argument('--basepath', type=str, default="../../../../data/merged/")
    parser.add_argument('--basepath_plot', type=str, default="../../../../plots/")
    parser.add_argument('--ylims', type=float, nargs=2, default=None)
    parser.add_argument('--xlims', type=float, nargs=2, default=None)
    parser.add_argument('--outputfolder', type=str, default=None)
    parser.add_argument('--ticklocations', nargs='*', type=float)
    parser.add_argument('--leg_loc', type=str)
    parser.add_argument('--leg_pos', nargs=2, type=float)
    parser.add_argument('--leg_ncol', type=int, default=2)
    parser.add_argument('--suffix', type=str, default="")
    parser.add_argument('--leg_lw', default=0.5, type=float)
    parser.add_argument('--leg_pad', default=0, type=float)
    args = parser.parse_args()

    if args.outputfolder is None:
        args.outputfolder = lpd.get_plot_path(args.qcdtype, args.corr, args.conftype, args.basepath_plot)

    return args


def load_data(args):

    fermions, temp, flowtype, gaugeaction, flowaction = lpd.parse_qcdtype(args.qcdtype)
    inputfolder = lpd.get_merged_data_path(args.qcdtype, args.corr, args.conftype, args.basepath)
    _, _, nt, _ = lpd.parse_conftype(args.conftype)
    flowtimesBy_a2 = numpy.loadtxt(inputfolder+"flowtimes_"+args.conftype+".dat")
    flowtimesT2 = flowtimesBy_a2 / nt**2
    XX = numpy.loadtxt(inputfolder + "/" + args.corr + "_" + args.conftype + ".dat")
    XX = apply_tree_level_imp(XX, nt, [0 for _ in range(len(flowtimesBy_a2))], args.corr, flowaction, gaugeaction)
    XX_err = numpy.loadtxt(inputfolder + "/" + args.corr + "_err_" + args.conftype + ".dat")
    XX_err = numpy.fabs(apply_tree_level_imp(XX_err, nt, [0 for _ in range(len(flowtimesBy_a2))], args.corr, flowaction, gaugeaction))

    tauT = lpd.get_tauTs(nt)
    return tauT, flowtimesT2, XX, XX_err


def plot(args, tauT, flowtimesT2, XX, XX_err):

    XX = numpy.swapaxes(XX, 0, 1)
    XX_err = numpy.swapaxes(XX_err, 0, 1)

    fig, ax, plots = lpd.create_figure(xlims=args.xlims, ylims=args.ylims, xlabel=r'$8\tau_\mathrm{F} / \tau^2$', ylabel=r'$\displaystyle \frac{G}{G^{\mathrm{norm}}}$')

    ax.yaxis.set_label_coords(0.07, 0.97)
    labels = []
    for i in range(len(XX)):
        color = lpd.get_color(range(len(XX)), i, 0, len(XX)-1)
        handle1 = ax.errorbar(8*flowtimesT2/tauT[i]**2, XX[i], color=color, fmt='o', fillstyle='full', markersize=0.5, zorder=-i)
        # ax.errorbar(8 * flowtimesT2 / tauT[i] ** 2, XX[i], color=color, fmt='-', markersize=0, lw=0.25, zorder=-i)
        handle2 = ax.fill_between(8*flowtimesT2/tauT[i]**2, XX[i]-XX_err[i], XX[i]+XX_err[i], facecolor=color, alpha=0.5, zorder=-100-i)
        plots.append((handle1, handle2))
        labels.append(lpd.format_float(tauT[i]))

    # second x-axis for flow radius
    ax2 = ax.twiny()
    new_tick_locations = numpy.array(args.ticklocations)**2
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(new_tick_locations)
    ax2.set_xticklabels(["%.1f" % z if z == 0 else "%.2f" % z for z in numpy.sqrt(new_tick_locations)])
    ax2.set_xlabel(r'$' + r'{\sqrt{8\tau_\mathrm{F}}}/{\tau}$', horizontalalignment='right', verticalalignment='top', zorder=999999,
                   bbox=lpd.labelboxstyle)
    ax2.xaxis.set_label_coords(0.99, 0.97)
    ax2.tick_params(direction='in', pad=0, width=0.5)

    # ax.legend(plots, labels)
    # handles, labels = ax.get_legend_handles_labels()
    leg = ax.legend(plots[::-1], labels[::-1], title=r'$\tau T$', loc=args.leg_loc, borderaxespad=args.leg_pad, handletextpad=0.2, columnspacing=0.5,
                    labelspacing=0.1, bbox_to_anchor=args.leg_pos, **lpd.leg_err_size(1, 0.3), framealpha=0.9, edgecolor="k", ncol=args.leg_ncol, fontsize=10)
    leg.get_frame().set_linewidth(args.leg_lw)

    file = args.outputfolder+"/"+args.corr+"_"+args.conftype+"_flow_dep"+args.suffix+".pdf"
    print("saving", file)
    fig.savefig(file)


def main():

    args = get_args()
    tauT, flowtimesT2, XX, XX_err = load_data(args)
    plot(args, tauT, flowtimesT2, XX, XX_err)

    return


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
