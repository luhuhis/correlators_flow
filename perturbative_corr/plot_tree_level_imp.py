#!/usr/bin/env python3
import lib_process_data as lpd
import numpy
import argparse
import matplotlib.pyplot
import scipy.optimize
import scipy.interpolate


def fit_ansatz(x, b, d):
    return b/x**2 + d/x**4 + 1


def chisq_dof(fitparams, ydata, xdata, edata, extrapolation_ansatz):
    ndata = len(ydata)
    chisq = 0
    for i in range(ndata):
        chisq += (extrapolation_ansatz(xdata[i], *fitparams) - ydata[i])**2 / edata[i]**2

    nfitparams = len(fitparams)
    chisqdof = chisq / (ndata-nfitparams)
    return chisqdof


def fit_sample(ydata, xdata, edata, extrapolation_ansatz):
    start_params = [0, 0]

    fitparams, _ = scipy.optimize.curve_fit(extrapolation_ansatz, xdata, ydata, p0=start_params, sigma=edata)
    chisqdof = chisq_dof(fitparams, ydata, xdata, edata, extrapolation_ansatz)
    return [*fitparams, chisqdof]


def combined_fit_ansatz(x, tauT, a, b0, c0): #, b1, c0, c1, d0, d1):
    return a + (b0/tauT**2+c0/tauT**4)*x #+ (d0/tauT**4)*x**2 #+ # +c1/tauT**8  (d0/tauT**6+d1/tauT**12)*x**3  +b1/tauT**4
    # return 1 + b0*(x/tauT) + c0*(x/tauT)**2+c1 * (x/tauT)**3


def combined_chisqdof(fitparams, relflow_corrs, Nts, tauTs):
    ndata = relflow_corrs.size
    ntauT = relflow_corrs.shape[1]
    chisq = 0

    combined_fitparams = fitparams[ntauT:]

    Nts = numpy.asarray(Nts)
    for i in range(ntauT):
        chisq += numpy.sum((combined_fit_ansatz(1/Nts**2, tauTs[i], fitparams[i], *combined_fitparams) - relflow_corrs[:,i])**2)

    nfitparams = len(fitparams)
    chisqdof = chisq / (ndata-nfitparams)
    return chisqdof


def perform_combined_fit(relflowcorrs, Nts, tauTs):
    ntauT = len(tauTs)
    fitparams = scipy.optimize.minimize(combined_chisqdof, x0=numpy.asarray([*[1 for _ in range(ntauT)], 1, 1]), args=(relflowcorrs, Nts, tauTs))
    fitparams = fitparams.x

    chisqdof = combined_chisqdof(fitparams, relflowcorrs, Nts, tauTs)

    return [*fitparams, chisqdof]


def plot_extr(args):
    ylabel = r'$\displaystyle \frac{{G^\mathrm{latt}_\mathrm{LO}}}{{G^\mathrm{cont}_\mathrm{LO}}}$'

    xlims = (0.24, 0.51)
    ylims = (0.995, 1.02)

    gauge_action = "LW"
    flow_action = "Zeuthen"

    flowtimes = numpy.loadtxt(args.flowtime_file)
    nflow = len(flowtimes)
    Nts = numpy.asarray([40, 48, 56, 64])

    ref_Nt = 20
    mindtauTidx = 4 # int(ref_Nt / 4 - 1)
    tauTs = lpd.get_tauTs(ref_Nt)[mindtauTidx:]
    # tauTs = lpd.get_tauTs(Nts[-1])

    corrs = []
    for Nt in Nts:
        file = args.inputpath + "/" + args.corr + "_pert_latt_" + flow_action + "_" + gauge_action + "_Nt" + str(Nt) + ".dat"
        corr = numpy.loadtxt(file)
        print(corr.shape)
        tmp = []
        nt_quart = int(Nt/4)
        minidx = nt_quart-1
        # interpolate to correct tauTs
        for f in range(nflow):
            validtauTs = lpd.get_tauTs(Nt)[minidx:]
            y = corr[f, minidx:]/lpd.EE_cont_LO(validtauTs, numpy.sqrt(8*flowtimes[f])/Nt)
            spline = scipy.interpolate.CubicSpline(validtauTs, y, bc_type=((2, 0.0), (1, 0.0)))  # true interpolating spline
            tmp2 = []
            for tauT in tauTs:
                tmp2.append(spline(tauT))
            tmp.append(numpy.asarray(tmp2))
        corrs.append(numpy.asarray(tmp))

    corrs = numpy.asarray(corrs)
    print(corrs.shape)

    # absolute flow time
    fig, ax, _ = lpd.create_figure(xlims=xlims, ylims=ylims, ylabel=ylabel, xlabel=r'$\tau T$')
    flowradius_T = 0.1
    for i, Nt in enumerate(Nts):
        flowradii_T = numpy.sqrt(8 * flowtimes) / Nt
        flowidx = numpy.fabs(flowradii_T - flowradius_T).argmin()
        ax.errorbar(tauTs, corrs[i][flowidx, :], fmt='x')

    ax.text(0.99, 0.99, r'$\sqrt{8\tau_\mathrm{F}}T=' + lpd.format_float(flowradius_T) + r'$', ha='right', va='top', transform=ax.transAxes)
    ax.legend(title=r'$\tau T$', handlelength=1, loc="upper right", bbox_to_anchor=(1, 0.95), alignment="left")
    fig.savefig(args.outputpath + "/pert_latt_comparison_" + args.corr + "_latt_effects.pdf")
    matplotlib.pyplot.close(fig)

    # relative flow time!
    flowradiusbytauT = 0.3
    relflowcorrs = []
    for i, Nt in enumerate(Nts):
        flowradii_T = numpy.sqrt(8 * flowtimes) / Nt
        Ntcorr = []
        for j, tauT in enumerate(tauTs):
            spline = scipy.interpolate.CubicSpline(flowradii_T, corrs[i][:, j], bc_type=((2, 0.0), (1, 0.0)))  # true interpolating spline
            # flowidx = numpy.fabs(flowradii_T - flowradiusbytauT * tauT).argmin()
            # Ntcorr.append(corrs[i][flowidx, j])
            Ntcorr.append(spline(flowradiusbytauT * tauT))
        relflowcorrs.append(numpy.asarray(Ntcorr))

    relflowcorrs = numpy.asarray(relflowcorrs)

    fig, ax, _ = lpd.create_figure(xlims=xlims, ylims=ylims, ylabel=ylabel, xlabel=r'$\tau T$')
    for i, Nt in enumerate(Nts):
        ax.errorbar(tauTs, relflowcorrs[i], fmt='x')
    ax.text(0.99, 0.99, r'$\sqrt{8\tau_\mathrm{F}}/\tau=' + lpd.format_float(flowradiusbytauT) + r'$', ha='right', va='top', transform=ax.transAxes)
    ax.legend(title=r'$\tau T$', handlelength=1, loc="upper right", bbox_to_anchor=(1, 0.95), alignment="left")
    fig.savefig(args.outputpath + "/pert_latt_comparison_" + args.corr + "_latt_effects_relflow.pdf")
    matplotlib.pyplot.close(fig)

    colors = ["C"+str(i) for i in range(10)]

    # cont extr
    x = 1/Nts**2
    xlims = (0, 1/Nts[0]**2 * 1.05)
    fig, ax, _ = lpd.create_figure(xlims=xlims, ylims=ylims, ylabel=ylabel, xlabel=r'$1/N_\tau^2$')
    # ax.set_yscale('log')

    for i, tauT in enumerate(tauTs):
        ax.errorbar(x, relflowcorrs[:, i], fmt='x', label=lpd.format_float(tauT), color=colors[i])

    fitparams = perform_combined_fit(relflowcorrs, Nts, tauTs)
    print(fitparams)
    xpoints = numpy.linspace(0, 1/Nts[0]**2, 100)
    for i, tauT in enumerate(tauTs):
        print(lpd.format_float(tauT), "-> chisqdof=", fitparams[-1])
        ax.errorbar(xpoints, combined_fit_ansatz(xpoints, tauT, fitparams[i], *fitparams[len(tauTs):-1]), fmt='-', color=colors[i])

    ax.text(0.99, 0.99, r'$\sqrt{8\tau_\mathrm{F}}/\tau=' + lpd.format_float(flowradiusbytauT) + r'$', ha='right', va='top', transform=ax.transAxes)
    ax.legend(title=r'$\tau T$', handlelength=1, loc="center left", bbox_to_anchor=(0, 0.5), alignment="left")
    fig.savefig(args.outputpath + "/pert_latt_comparison_" + args.corr + "_extr.pdf")
    matplotlib.pyplot.close(fig)




def plot_tau(args):
    ylabel = r'$\displaystyle \frac{{G^\mathrm{latt}_\mathrm{LO}}}{{G^\mathrm{cont}_\mathrm{LO}}}$'

    xlims = (0, args.Nt/2+1)
    ylims = (0.95, 1.45)

    fig, ax, _ = lpd.create_figure(xlims=xlims, ylims=ylims, ylabel=ylabel, xlabel=r'$\tau/a$')

    gauge_actions = ("Wilson", "LW")
    flow_actions = ("Wilson", "Zeuthen")  # , "LW")

    nt_half = int(args.Nt/2)

    flowtimes = numpy.loadtxt(args.flowtime_file)
    flowradiitimesT = numpy.sqrt(8 * flowtimes) / args.Nt

    flowradiusbytauT = 0.3

    start_idx = 4

    counter = 0
    fmts = ('x', '+', 'o', 'o', '--', ':')
    colors = ('C1', 'C2', 'C3', 'C4')
    for gauge_action in gauge_actions:
        for flow_action in flow_actions:
            if not (gauge_action == "LW" and flow_action == "Wilson"):
                file = args.inputpath + "/" + args.corr + "_pert_latt_" + flow_action + "_" + gauge_action + "_Nt" + str(args.Nt) + ".dat"
                try:
                    corr = numpy.loadtxt(file)
                    x = range(1,nt_half+1)
                    y = []
                    for i in x:
                        tauT = float(i)/args.Nt
                        flowidx = numpy.fabs(flowradiitimesT-flowradiusbytauT*tauT).argmin()
                        y.append(corr[flowidx, i-1] / (lpd.EE_cont_LO(float(i)/args.Nt, flowradiitimesT[flowidx])))
                    if gauge_action == "LW":
                        gauge_action = "Symanzik 2x1"  # rename
                    ax.errorbar(x[start_idx:], y[start_idx:], label=gauge_action + ", " + flow_action, fmt=fmts[counter], fillstyle='none', color=colors[counter])

                    x = numpy.asarray(x)
                    y = numpy.asarray(y)
                    results = fit_sample(y[start_idx:], x[start_idx:], [1 for _ in range(len(x[start_idx:]))], fit_ansatz)
                    print(results)

                    xplot = numpy.linspace(5, nt_half, 100)
                    ax.errorbar(xplot, fit_ansatz(xplot, *results[:-1]), fmt='-', color=colors[counter])

                    counter += 1
                except OSError:
                    print("did not find", file)
                    pass

    ax.text(0.99, 0.99, r'$\sqrt{8\tau_\mathrm{F}}/\tau=' + lpd.format_float(flowradiusbytauT) + r'$', ha='right', va='top', transform=ax.transAxes)
    ax.axhline(y=1, **lpd.horizontallinestyle)
    # ax.axvline(x=0.33, **lpd.verticallinestyle)
    ax.legend(title="$S_\\mathrm{gauge}, S_\\mathrm{flow}$", handlelength=1, loc="upper right", bbox_to_anchor=(1, 0.95), alignment="left")
    fig.savefig(args.outputpath + "/pert_latt_comparison_" + args.corr + "_" + "Nt" + str(args.Nt) + ".pdf")
    matplotlib.pyplot.close(fig)


def plot_flow(args):

    if not args.xlims:
        args.xlims = (0, 0.33)
    if not args.ylims:
            args.ylims = (0.9, 1.55)

    flowtimes = numpy.loadtxt(args.flowtime_file)

    ylabel = r'$\displaystyle \frac{{G^\mathrm{latt}_\mathrm{LO}}}{{G^\mathrm{cont}_\mathrm{LO}}}$'

    fig, ax, _ = lpd.create_figure(xlims=args.xlims, ylims=args.ylims, ylabel=ylabel, xlabel="$\\sqrt{8\\tau_F}/\\tau$")

    ax.errorbar(0, 0, label="$S_\\mathrm{gauge}, S_\\mathrm{flow}$", markersize=0, alpha=0, lw=0)

    gauge_actions = ("Wilson", "LW")
    flow_actions = ("Wilson", "Zeuthen")  # "LW")

    tauT = args.tau/args.Nt

    counter = 0
    fmts = ('-', '--', ':', '-', '--', ':')
    for gauge_action in gauge_actions:
        for flow_action in flow_actions:
            if not (gauge_action == "LW" and flow_action == "Wilson"):
                file = args.inputpath+"/"+args.corr+"_pert_latt_" + flow_action + "_" + gauge_action + "_Nt" + str(args.Nt) + ".dat"
                try:
                    corr = numpy.loadtxt(file)
                    i = args.tau - 1
                    flowradiitimesT = numpy.sqrt(8*flowtimes)/args.Nt
                    x = numpy.sqrt(8*flowtimes)/args.tau
                    y = (corr[:, i])/(lpd.EE_cont_LO(tauT, flowradiitimesT))  # flowradiitimesT
                    if gauge_action == "LW":
                        gauge_action = "Symanzik 2x1"  # rename
                    ax.errorbar(x, y, label=gauge_action+", "+flow_action, fmt=fmts[counter])
                    counter += 1
                except OSError:
                    print("did not find", file)
                    pass

    ax.text(0.99, 0.99, r'$\tau/a='+str(args.tau)+"$", ha='right', va='top', transform=ax.transAxes)
    ax.axhline(y=1, **lpd.horizontallinestyle)
    # ax.axvline(x=0.33, **lpd.verticallinestyle)
    ax.legend(handlelength=1, loc="upper right", bbox_to_anchor=(1, 0.95), alignment="left")
    add_suffix = ""
    fig.savefig(args.outputpath + "/pert_latt_comparison_" + args.corr + "_" + "Nt" + str(args.Nt) + "_" + str(args.tau) + add_suffix + ".pdf")


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

    plot_extr(args)
    # plot_tau(args)
    # plot_flow(args)


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
