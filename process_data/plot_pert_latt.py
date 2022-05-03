#!/usr/local/bin/python3.7 -u
import lib_process_data as lpd
import numpy
import argparse


def Gnorm(tauT):
    norm = numpy.pi ** 2 * (numpy.cos(numpy.pi * tauT) ** 2 / numpy.sin(numpy.pi * tauT) ** 4 + 1 / (3 * numpy.sin(numpy.pi * tauT) ** 2))
    return norm


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--Nt', type=int, required=True)
    parser.add_argument('--corr', choices=["EE", "BB"], type=str, required=True)
    parser.add_argument('--flowtime_file', required=True, type=str)
    parser.add_argument('--outputpath', required=True, type=str)
    parser.add_argument('--inputpath', required=True, type=str)
    parser.add_argument('--tau', type=int, required=True)
    parser.add_argument('--xlims', help="xlim. default: 0 to largest flow time", type=float, nargs=2)
    parser.add_argument('--ylims', type=float, nargs=2)
    args = parser.parse_args()

    flowtimes = numpy.loadtxt(args.flowtime_file)

    if not args.xlims:
        args.xlims = (0, flowtimes[-1]*1.02)
    if not args.ylims:
        args.ylims = (0, 1.1)

    fig, ax, _ = lpd.create_figure(figsize=(16/9*(3 + 3 / 8), (3 + 3 / 8 - 1 / 2.54)), xlabelpos=(0.93, 0.15),
                                   ylabelpos=(0.03, 0.95), xlims=args.xlims, ylims=args.ylims)

    ax.set_ylabel("$G^\\mathrm{latt}/G^\\mathrm{cont}$", fontsize=12)
    ax.set_xlabel("$\\sqrt{8\\tau_F}T$", fontsize=12)

    gauge_actions = ("Wilson", "LW")
    flow_actions = ("Wilson", "Zeuthen", "LW")

    tauT = args.tau/args.Nt

    counter = 0
    fmts = ('-', '--', ':', '-', '--', ':')

    for gauge_action in gauge_actions:
        for flow_action in flow_actions:
            try:
                corr = numpy.loadtxt(args.inputpath+"/"+args.corr+"_pert_latt_" + flow_action + "_" + gauge_action + "_Nt" + str(args.Nt) + ".dat")
                # for i, tauT in enumerate(tauTs):
                i = args.tau - 1
                ax.errorbar(numpy.sqrt(8*flowtimes)/args.Nt, numpy.fabs(corr[:, i]/Gnorm(tauT)), label=gauge_action+", "+flow_action, lw=1, fmt=fmts[counter])
                counter += 1
            except OSError:
                pass

    ax.set_title(args.corr + ", $\\tau T = "+'{0:.3f}'.format(tauT)+"$")
    ax.axhline(y=1, **lpd.horizontallinestyle)
    ax.legend(title="$S_\\mathrm{gauge}, S_\\mathrm{flow}$", **lpd.legendstyle)

    fig.savefig(args.outputpath + "/pert_latt_comparison_"+ args.corr + "_" + "Nt" + str(args.Nt) + "_" + '{0:.3f}'.format(tauT) + ".pdf")


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
