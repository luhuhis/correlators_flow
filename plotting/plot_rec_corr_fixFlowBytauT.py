#!/usr/local/bin/python3.7m -u

import lib_process_data as lpd
import numpy
import scipy.integrate
import scipy.interpolate
import argparse


def Gnorm(tauT):
    norm = numpy.pi ** 2 * (numpy.cos(numpy.pi * tauT) ** 2 / numpy.sin(numpy.pi * tauT) ** 4 + 1 / (3 * numpy.sin(numpy.pi * tauT) ** 2))
    return norm


# ====================


def model3(OmegaByT, kappa, PhiuvByT3, UV_prefactor):
    return numpy.maximum(PhiIR(OmegaByT, kappa), UV_prefactor * PhiuvByT3(OmegaByT))


def PhiIR(OmegaByT, kappa):
    return kappa * OmegaByT / 2


def Kernel(OmegaByT, tauT):
    return numpy.cosh(OmegaByT / 2 - OmegaByT * tauT) / numpy.sinh(OmegaByT / 2)


def Integrand(OmegaByT, tauT, spf, *args):
    return 1. / numpy.pi * Kernel(OmegaByT, tauT) * spf(OmegaByT, *args)


def main():
    parser = argparse.ArgumentParser()

    # TODO: make this script read in any correlator data to find the constant flowradiusBytauT.
    #  plain lattice correlator data has a different layout (all flow times in one file) so need to adjust for that
    parser.add_argument('--fitparam_files', nargs='*', help='plot model 3 spf for these fit params', required=True)
    parser.add_argument('--use_tex', action="store_true")
    parser.add_argument('--npoints', help='how many tauT to plot between 0 and 0.5 for the model correlators', default=100, type=int)
    parser.add_argument('--ylims', help='custom ylims', nargs=2, type=float, default=(-0.03, 4))
    parser.add_argument('--xlims', help='custom xlims', nargs=2, type=float, default=(-0.01, 0.51))
    parser.add_argument('--flowradiusBytauT', nargs='*', type=float, help='which sqrt(8tauF)/tauT to plot', required=True)
    parser.add_argument('--min_tauT', default=0.2, type=float, help='from which tauT on the reconstructed corr should be plotted')

    args = parser.parse_args()

    basepath = "/work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/EE/"
    basepath_plot = "/work/home/altenkort/work/correlators_flow/plots/quenched_1.50Tc_zeuthenFlow/EE/"

    # TODO abstract this or make an argument out of this
    XX = []
    XX.append(numpy.loadtxt(basepath+"EE_flow_extr.txt", unpack=True))

    # load PhiUV data
    PhiUVa = numpy.loadtxt("/work/data/htshu/ee_spf/PHIUV_a.dat", unpack=True)[0:3]
    PhiUVaByT3 = scipy.interpolate.InterpolatedUnivariateSpline(PhiUVa[0], PhiUVa[1], k=3)

    MinOmegaByT = 0
    MaxOmegaByT = PhiUVa[0][-1]

    # calc reconstructed corrs for different parameters
    xpoints = numpy.linspace(args.min_tauT, 0.5, args.npoints)
    model3_means = []
    model3_uppers = []
    model3_lowers = []
    fitparams = []
    for file in args.fitparam_files:
        params, err_left, err_right = numpy.loadtxt(file, unpack=True)
        fitparams.append(numpy.column_stack((params, err_left, err_right)))
        model3_mean = []
        model3_lower = []
        model3_upper = []
        for flowindex in xpoints:
            model3_mean.append(scipy.integrate.quad(
                lambda OmegaByT: Integrand(OmegaByT, flowindex, model3, params[0], PhiUVaByT3, params[1]), MinOmegaByT, MaxOmegaByT)[0] / Gnorm(flowindex))
            model3_lower.append(scipy.integrate.quad(
                lambda OmegaByT: Integrand(OmegaByT, flowindex, model3, params[0]-err_left[0], PhiUVaByT3, params[1]-err_left[1]), MinOmegaByT, MaxOmegaByT)[0] / Gnorm(flowindex))
            model3_upper.append(scipy.integrate.quad(
                lambda OmegaByT: Integrand(OmegaByT, flowindex, model3, params[0]+err_right[0], PhiUVaByT3, params[1]+err_right[1]), MinOmegaByT, MaxOmegaByT)[0] / Gnorm(flowindex))
        model3_means.append(model3_mean)
        model3_lowers.append(model3_lower)
        model3_uppers.append(model3_upper)

    # ==========================================================================================================================================================

    EE_cont_arr = []
    flowradii = numpy.loadtxt("/work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/EE/s144t36_b0754400/flowradii_s144t36_b0754400.dat")
    min_idx = 0
    max_idx = 0
    stop = False
    for flowradius in flowradii:
        flowradius_str = r'{0:.4f}'.format(flowradius)
        file = "/work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/EE/cont_extr/EE_" + flowradius_str + "_cont.txt"
        try:
            EE_cont_arr.append(numpy.loadtxt(file, unpack=True))
            stop = True
            max_idx += 1
        except OSError:
            # print("could not find ", file)
            if not stop:
                min_idx += 1
    max_idx += min_idx
    flowradii = flowradii[min_idx:max_idx]  # filter which flowradii got actually read in
    possible_tauTs = EE_cont_arr[0][0]  # doesn't matter which one we choose here, they are all the same

    # interpolate between flowtimes for each tauT. for that we first need to remove the nan's
    cont_int = []
    cont_err_int = []
    for i in range(len(possible_tauTs)):
        ydata = numpy.asarray([EE_cont[1][i] for EE_cont in EE_cont_arr])
        edata = numpy.asarray([EE_cont[2][i] for EE_cont in EE_cont_arr])
        mask = ~numpy.isnan(ydata)
        xdata = flowradii[mask]
        ydata = ydata[mask]
        edata = edata[mask]
        if len(xdata) >= 3:
            cont_int.append(scipy.interpolate.InterpolatedUnivariateSpline(xdata, ydata, k=3, ext=2, check_finite=True))
            cont_err_int.append(scipy.interpolate.InterpolatedUnivariateSpline(xdata, edata, k=3, ext=2, check_finite=True))
        else:
            cont_int.append(None)
            cont_err_int.append(None)

    fig2, ax2, _ = lpd.create_figure(xlims=args.xlims, ylims=args.ylims, xlabel=r'$\tau T$', xlabelpos=(0.95, 0.05), ylabelpos=(0.08, 0.9),
                                     UseTex=args.use_tex,
                                     ylabel=r'$\frac{G^*}{ G^\mathrm{norm} }$', figsize=(1.5 * (3 + 3 / 8), 1.5 * (3 + 3 / 8) / 16 * 9))

    # plot cont+flow extr. corr
    ax2.errorbar(XX[-1][0], XX[-1][1], XX[-1][2], zorder=0, **lpd.chmap(lpd.plotstyle_add_point, color='green', fmt='-', label=r'$G^{\tau_F\rightarrow 0}_\mathrm{cont}$'))

    # find frankenstein correlator at various fixed flowradiusBytauT
    frankenstein_arr = []
    for f, flowradiusBytauT in enumerate(args.flowradiusBytauT):
        frankenstein_edata = []
        frankenstein_ydata = []
        frankenstein_xdata = []
        for i, tauT in enumerate(possible_tauTs):
            wanted_flowradius = flowradiusBytauT * tauT
            if cont_int[i] is not None:
                try:
                    frankenstein_edata.append(cont_err_int[i](wanted_flowradius))
                    frankenstein_ydata.append(cont_int[i](wanted_flowradius))
                    frankenstein_xdata.append(tauT)
                except ValueError:
                    pass
        frankenstein_arr.append(numpy.column_stack((frankenstein_xdata, frankenstein_ydata, frankenstein_edata)))
        numpy.savetxt("/work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/EE/frankenstein_"+r'{0:.3f}'.format(flowradiusBytauT)+".txt",
                      numpy.column_stack((frankenstein_xdata, frankenstein_ydata, frankenstein_edata)),
                      header="tauT      G/Gnorm_sqrt(8tauF)/tau="+r'{0:.3f}'.format(flowradiusBytauT)+"       err")
        ax2.errorbar(frankenstein_xdata, frankenstein_ydata, frankenstein_edata, zorder=-len(args.flowradiusBytauT)+f,
                     **lpd.chmap(lpd.plotstyle_add_point, markersize=0, fmt='x-', color=lpd.get_color(args.flowradiusBytauT, f),
                                 label=r'$G_{\sqrt{8\tau_F}/\tau='+r'{0:.3f}'.format(flowradiusBytauT)+r'}}$'))  #

    # plot reconstructed correlators with "error bands"
    colors = []
    for i in range(len(model3_means)):
        colors.append(lpd.get_color(range(len(model3_means)), i))
    for model3_mean, model3_upper, model3_lower, fitparam, color in zip(model3_means, model3_uppers, model3_lowers, fitparams, colors):
        ax2.fill_between(xpoints, model3_lower, model3_upper, linewidth=0,  color=color, alpha=0.3, zorder=-1000)
        ax2.errorbar(xpoints, model3_mean, color=color, alpha=0.3,
                     label="$G^{\mathrm{max(IR,c UV_{a})}}$ \n $\kappa/T^3=" + "{0:.3f}".format(fitparam[0][0])
                           + "\pm " + "{0:.3f}".format(numpy.fmax(fitparam[0][1], fitparam[0][2]))
                           +"$\n $c=" + "{0:.3f}".format(fitparam[1][0]) + "\pm " + "{0:.3f}".format(numpy.fmax(fitparam[1][1], fitparam[1][2])) + "$",
                     **lpd.chmap(lpd.plotstyle_add_point, fmt='--', zorder=-100))

    lpd.legendstyle.update(dict(loc="lower right", bbox_to_anchor=(1.01, 0.07), labelspacing=1))
    ax2.legend(title=r'', **lpd.chmap(lpd.legendstyle, bbox_to_anchor=(1, 0.5), loc='center left', ncol=1, fontsize=6))
    ax2.set_title(r'quenched, $1.5 T_c$, IR=$\kappa \omega /2T$', x=0.5, y=1, fontsize=8)

    fig2.savefig(basepath_plot + "EE_frankenstein_comparison.pdf")

    # ==========================================================================================================================================================

    # plot normalized to UV part

    # fig3, ax3, _ = lpd.create_figure(xlims=args.xlims, xlabel=r'$\tau T$', xlabelpos=(0.95, 0.05), ylabelpos=(0.08, 0.9),
    #                                  UseTex=args.use_tex,
    #                                  ylabel=None, figsize=(1.5 * (3 + 3 / 8), 1.5 * (3 + 3 / 8) / 16 * 9))
    #
    # for i, flowradiusBytauT in enumerate(args.flowradiusBytauT):
    #     # spfs
    #     PhiUVa_corr = []
    #     PhiIRacorr = []
    #     model3acorr = []
    #     for flowindex in frankenstein_arr[i][:, 0]:
    #         PhiUVa_corr.append(args.c_a*scipy.integrate.quad(lambda OmegaByT: Integrand(OmegaByT, flowindex, PhiUVaByT3), 0, MaxOmegaByT)[0] / Gnorm(flowindex))
    #         PhiIRacorr.append(scipy.integrate.quad(lambda OmegaByT: Integrand(OmegaByT, flowindex, PhiIR, args.kappa_1), 0, MaxOmegaByT)[0] / Gnorm(flowindex))
    #         model3acorr.append(scipy.integrate.quad(lambda OmegaByT: Integrand(OmegaByT, flowindex, model3, args.kappa_1, PhiUVaByT3, args.c_a), 0, MaxOmegaByT)[0] / Gnorm(flowindex))
    #
    #     ax3.errorbar(frankenstein_arr[i][:, 0], frankenstein_arr[i][:, 1]/PhiUVa_corr, frankenstein_arr[i][:, 2]/PhiUVa_corr,
    #                  **lpd.chmap(lpd.plotstyle_add_point, fmt='-', label=r'$\frac{G_{\sqrt{8\tau_F}/\tau='+r'{0:.3f}'.format(flowradiusBytauT)+r'}}{ c\Phi_\mathrm{UV}^a}$'))
    #
    # lpd.legendstyle.update(dict(loc="lower right", bbox_to_anchor=(1.01, 0.07)))
    # ax3.legend(**lpd.chmap(lpd.legendstyle, bbox_to_anchor=(1, 0.5), loc='center left'))
    # ax3.set_title(r'quenched, $1.5 T_c$', x=0.5, y=1, fontsize=8)
    #
    # fig3.savefig(basepath_plot + "EE_frankenstein_comparison_ratio.pdf")


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
