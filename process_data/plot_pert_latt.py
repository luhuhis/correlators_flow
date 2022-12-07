#!/usr/bin/env python3
import lib_process_data as lpd
import numpy
import argparse


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
    parser.add_argument('--exp', default=False, action="store_true")
    args = parser.parse_args()

    flowtimes = numpy.loadtxt(args.flowtime_file)

    if not args.xlims:
        args.xlims = (0, 1)
    if not args.ylims:
        if args.exp:
            args.ylims = (0.01, 1000000)
        else:
            args.ylims = (0.8, 1.5)

    ylabel = r'$\displaystyle \frac{{G^\mathrm{latt}}}{{G^\mathrm{cont}}}$'
    if args.exp:
        ylabel = r'$\displaystyle \frac{\exp{G^\mathrm{latt}}}{\exp{G^\mathrm{cont}}}$'

    fig, ax, _ = lpd.create_figure(xlims=args.xlims, ylims=args.ylims, ylabel=ylabel, xlabel="$\\sqrt{8\\tau_F}/\\tau$")

    gauge_actions = ("Wilson", "LW")
    flow_actions = ("Wilson", "Zeuthen") #, "LW")

    tauT = args.tau/args.Nt

    counter = 0
    fmts = ('-', '--', ':', '-', '--', ':')
    if args.exp:
        ax.set_yscale('log')
    for gauge_action in gauge_actions:
        for flow_action in flow_actions:
            file = args.inputpath+"/"+args.corr+"_pert_latt_" + flow_action + "_" + gauge_action + "_Nt" + str(args.Nt) + ".dat"
            try:
                corr = numpy.loadtxt(file)
                i = args.tau - 1
                flowradiitimesT = numpy.sqrt(8*flowtimes)/args.Nt
                x = numpy.sqrt(8*flowtimes)/args.tau
                if args.exp:
                    y = numpy.exp(corr[:, i]) / numpy.exp(lpd.EE_cont_LO(tauT, flowradiitimesT))
                else:
                    y = (corr[:, i])/(lpd.EE_cont_LO(tauT, flowradiitimesT))
                ax.errorbar(x, y, label=gauge_action+", "+flow_action, fmt=fmts[counter])
                counter += 1
            except OSError:
                print("did not find", file)
                pass

    ax.set_title(args.corr + ", $\\tau/a = "+str(args.tau)+"$")
    ax.axhline(y=1, **lpd.horizontallinestyle)
    ax.axvline(x=0.33, **lpd.verticallinestyle)
    ax.legend(title="$S_\\mathrm{gauge}, S_\\mathrm{flow}$", **lpd.legendstyle)
    add_suffix = ""
    if args.exp:
        add_suffix = "_exp"
    fig.savefig(args.outputpath + "/pert_latt_comparison_"+ args.corr + "_" + "Nt" + str(args.Nt) + "_" + str(args.tau) + add_suffix + ".pdf")


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
