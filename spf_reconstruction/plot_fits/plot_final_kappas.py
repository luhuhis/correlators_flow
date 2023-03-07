#!/usr/bin/env python3
import lib_process_data as lpd
import numpy
import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--corr', choices=["EE", "BB"])
    parser.add_argument('--input_kappa_files', nargs='*')
    parser.add_argument('--labels', nargs='*')
    parser.add_argument('--fillstyles', nargs='*')
    parser.add_argument('--fmts', nargs='*')
    parser.add_argument('--colors', nargs='*')
    parser.add_argument('--zorders', nargs='*', type=int)
    parser.add_argument('--markersize', type=int)
    parser.add_argument('--outputpath')
    parser.add_argument('--suffix')
    parser.add_argument('--temps_in_GeV', type=float, nargs="*")
    parser.add_argument('--Tc_in_GeV', type=float, required=True)
    parser.add_argument('--leg_hide', action="store_true")
    parser.add_argument('--leg_ncol', default=1, type=int)
    parser.add_argument('--xlims', nargs=2, type=float)
    parser.add_argument('--ylims', nargs=2, type=float)
    parser.add_argument('--plot_EE_quenched_lit', action="store_true")
    parser.add_argument('--plot_BB_quenched_lit', action="store_true")
    parser.add_argument('--plot_analytical_results', action="store_true")
    parser.add_argument('--no_subscript', action="store_true")
    args = parser.parse_args()
    return args


def load_data(args):
    data = []
    for file in args.input_kappa_files:
        data.append(numpy.loadtxt(file))
    return data


def mean(low, high):
    return (low+high)/2


def err(low, high):
    return high - (low+high)/2


def mean_and_err(low, high):
    return mean(low,high), err(low,high)


def thiscolor(i):
    return lpd.get_discrete_color(i)
    # return lpd.lighten_color(lpd.get_discrete_color(i), 0.75)


def plot_EE_quenched_literature(args, ax, i):

    ms=0
    myfmt='.'

    #  2015, A. Francis, O. Kaczmarek, M. Laine, T. Neuhaus, and H. Ohno, Phys. Rev. D 92, 116003
    ax.errorbar(1.46, 2.6, 0.8, color=thiscolor(3-i), fmt=myfmt, markersize=ms, label=r'Francis \textquotesingle 15 (ML)')
    i += 1

    #  2020, Nora Brambilla, Viljami Leino, Peter Petreczky, and Antonio Vairo, Phys. Rev. D 102, 074503
    color = thiscolor(3-i)
    ax.errorbar(1.1, *mean_and_err(1.91,5.40), markersize=ms, fmt=myfmt, color=color, label=r'TUMQCD \textquotesingle 20 (ML)')
    ax.errorbar(1.48, *mean_and_err(1.31,3.64), markersize=ms, fmt=myfmt, color=color)
    ax.errorbar(3,   *mean_and_err(0.63,2.20), markersize=ms, fmt=myfmt, color=color, zorder=-10)
    i += 1

    #  2022, Nora Brambilla, Viljami Leino, Julian Mayer-Steudte, Peter Petreczky, 	arXiv:2206.02861
    color = thiscolor(3-i)
    ax.errorbar(1.5, *mean_and_err(1.70,3.12), markersize=ms, fmt=myfmt, color=color, label=r'TUMQCD \textquotesingle 22 (flow)')
    i += 1

    # 2022, Debasish Banerjee, Rajiv Gavai, Saumen Datta, Pushan Majumdar, arXiv:2206.15471v1
    color = thiscolor(3-i)
    ax.errorbar(1.2, *mean_and_err(2.1,3.5), markersize=ms, fmt=myfmt, color=color, label=r'Banerjee \textquotesingle 22 (ML)')
    ax.errorbar(1.52, *mean_and_err(1.5,2.8), markersize=ms, fmt=myfmt, color=color)
    ax.errorbar(2.0, *mean_and_err(1.0,2.3), markersize=ms, fmt=myfmt, color=color)
    ax.errorbar(2.5, *mean_and_err(0.9,2.1), markersize=ms, fmt=myfmt, color=color)
    ax.errorbar(3.0, *mean_and_err(0.8,1.8), markersize=ms, fmt=myfmt, color=color)
    ax.errorbar(3.5, *mean_and_err(0.75, 1.5), markersize=ms, fmt=myfmt, color=color)
    i += 1

    return i


def plot_BB_quenched_literature(args, ax, i):
    ms = 0
    myfmt = '.'

    #  2022, D. Banerjee, S. Datta & M. Laine, JHEP08(2022)128
    ax.errorbar(1.2, *mean_and_err(1,2.6), color=thiscolor(i), fmt=myfmt, markersize=ms, label=r'Banerjee \textquotesingle 22 (ML)')
    ax.errorbar(1.47, *mean_and_err(1,2.1), color=thiscolor(i), fmt=myfmt, markersize=ms)
    ax.errorbar(2.0, *mean_and_err(0.6,1.8), color=thiscolor(i), fmt=myfmt, markersize=ms)
    i += 1

    # 2022, Nora Brambilla, Viljami Leino, Julian Mayer-Steudte, Peter Petreczky, 	arXiv:2206.02861
    ax.errorbar(1.5, *mean_and_err(1.23,2.74), color=thiscolor(i), fmt=myfmt, markersize=ms, label=r'TUMQCD \textquotesingle 22 (flow*)')
    i += 1

    return i


def plot_analytical_results(args, ax):
    xlims = ax.get_xlim()
    ax.errorbar([0, 100], [4 * numpy.pi / 0.9, 4 * numpy.pi / 0.9], alpha=0.5, zorder=-10000, color='grey', fmt=':')
    ax.errorbar([1.5, 100], [4 * numpy.pi / 6, 4 * numpy.pi / 6], alpha=0.5, zorder=-10000, color='grey', fmt='--')
    ax.text(xlims[1], 4 * numpy.pi / 0.9, r'AdS/CFT estimate', ha='right', va='bottom')
    ax.text(xlims[1], 4 * numpy.pi / 6, r'pQCD (NLO, $\alpha\approx 0.5)$', ha='right', va='bottom')


def plot_data(args, ax, data, offset, temps):
    for i, [kappa, kappa_err] in enumerate(data):
        if i == 0:
            color = 'k'
        elif args.plot_EE_quenched_lit or args.plot_BB_quenched_lit:
            color = lpd.get_discrete_color(offset)
            offset += 1
        ax.errorbar(temps[i]/args.Tc_in_GeV, kappa, kappa_err,
                    fmt='.' if args.fmts is None else args.fmts[i],
                    fillstyle=None if args.fillstyles is None else args.fillstyles[i],
                    markersize=0 if args.markersize is None else args.markersize,
                    color=color if args.colors is None else args.colors[i],
                    zorder=None if args.zorders is None else args.zorders[i],
                    label=args.labels[i])

    return offset


def do_plot(args, data, temps):

    if args.no_subscript:
        corr = ""
    else:
        corr = args.corr

    fig, ax, plots = lpd.create_figure(xlims=args.xlims, ylims=args.ylims, xlabel=r'$T/T_c$', ylabel=r'$\displaystyle \frac{\kappa'+lpd.get_corr_subscript(corr)+r'}{T^3}$')

    offset = 0
    if args.plot_EE_quenched_lit:
        offset = plot_EE_quenched_literature(args, ax, offset)

    if args.plot_BB_quenched_lit:
        offset = plot_BB_quenched_literature(args, ax, offset)

    offset = plot_data(args, ax, data, offset, temps)

    if args.plot_analytical_results:
        plot_analytical_results(args, ax)

    if not args.leg_hide:
        leg = ax.legend(**lpd.leg_err_size(1, 0.5), ncol=args.leg_ncol, columnspacing=1, labelspacing=0.5)
        for t in leg.texts:
            t.set_multialignment('left')

    filename = args.outputpath + "/kappa_"+args.suffix+".pdf"
    fig.savefig(filename)
    print("saved correlator plot", filename)


def main():

    args = parse_args()

    kappa = load_data(args)

    do_plot(args, kappa, args.temps_in_GeV)

    return


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
