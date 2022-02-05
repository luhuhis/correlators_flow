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
    parser.add_argument('--kappa_a', required=True, type=float, default=1)
    parser.add_argument('--c_a', help='prefactor of UV spf in model3b', required=True, type=float, default=1)
    parser.add_argument('--kappa_b', required=True, type=float, default=1)
    parser.add_argument('--c_b', help='prefactor of UV spf in model3a', required=True, type=float, default=1)
    parser.add_argument('--use_tex', action="store_true")
    parser.add_argument('--npoints', help='how many tauT to plot between 0 and 0.5 for the model correlators', default=100, type=int)
    parser.add_argument('--ylims', help='custom ylims', nargs=2, type=float, default=(-0.03, 4))
    parser.add_argument('--xlims', help='custom xlims', nargs=2, type=float, default=(-0.01, 0.51))

    args = parser.parse_args()

    basepath = "/work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/EE/"

    conftypes = ("s080t20_b0703500", "s096t24_b0719200", "s120t30_b0739400", "s144t36_b0754400")
    labels = ("20", "24", "30", "36", "cont")

    XX = []
    for conftype in conftypes:
        tmp = numpy.loadtxt(basepath+conftype+"/EE_flow_extr.txt", unpack=True)
        XX.append(tmp)
    XX.append(numpy.loadtxt(basepath+"EE_flow_extr.txt", unpack=True))

    # load PhiUV data
    PhiUVa = numpy.loadtxt("/work/data/htshu/ee_spf/PHIUV_a.dat", unpack=True)[0:3]
    PhiUVaByT3 = scipy.interpolate.InterpolatedUnivariateSpline(PhiUVa[0], PhiUVa[1], k=3) #, ext=2)
    PhiUVb = numpy.loadtxt("/work/data/htshu/ee_spf/PHIUV_b.dat", unpack=True)[0:3]
    PhiUVbByT3 = scipy.interpolate.InterpolatedUnivariateSpline(PhiUVb[0], PhiUVb[1], k=3) #, ext=2)

    MaxOmegaByT = PhiUVa[0][-1]

    # spfs
    PhiUVa_corr = []
    PhiUVb_corr = []
    PhiIRacorr = []
    PhiIRbcorr = []
    model3acorr = []
    model3bcorr = []
    xpoints = numpy.linspace(0, 0.5, args.npoints)
    for flowindex in xpoints:
        PhiUVa_corr.append(scipy.integrate.quad(lambda OmegaByT: Integrand(OmegaByT, flowindex, PhiUVaByT3), 0, MaxOmegaByT)[0] / Gnorm(flowindex))
        PhiUVb_corr.append(scipy.integrate.quad(lambda OmegaByT: Integrand(OmegaByT, flowindex, PhiUVbByT3), 0, MaxOmegaByT)[0] / Gnorm(flowindex))
        PhiIRacorr.append(scipy.integrate.quad(lambda OmegaByT: Integrand(OmegaByT, flowindex, PhiIR, args.kappa_a), 0, MaxOmegaByT)[0] / Gnorm(flowindex))
        PhiIRbcorr.append(scipy.integrate.quad(lambda OmegaByT: Integrand(OmegaByT, flowindex, PhiIR, args.kappa_b), 0, MaxOmegaByT)[0] / Gnorm(flowindex))
        model3acorr.append(scipy.integrate.quad(lambda OmegaByT: Integrand(OmegaByT, flowindex, model3, args.kappa_a, PhiUVaByT3, args.c_a), 0, MaxOmegaByT)[0] / Gnorm(flowindex))
        model3bcorr.append(scipy.integrate.quad(lambda OmegaByT: Integrand(OmegaByT, flowindex, model3, args.kappa_b, PhiUVbByT3, args.c_b), 0, MaxOmegaByT)[0] / Gnorm(flowindex))

    fig, ax, plots = lpd.create_figure(xlims=[-0.01, 0.51], ylims=(-0.03, 4), xlabel=r'$\tau T$', xlabelpos=(0.95, 0.05), ylabelpos=(0.08, 0.9), UseTex=args.use_tex,
                                       ylabel=r'$\frac{G^*}{ G^\mathrm{norm} }$', figsize=(1.5 * (3 + 3 / 8), 1.5 * (3 + 3 / 8) / 16 * 9))

    for i, data in enumerate(XX):
        ax.errorbar(data[0], data[1], data[2], **lpd.chmap(lpd.plotstyle_add_point, fmt='-', label=r'$G^{\tau_F\rightarrow 0}_{N_\tau='+labels[i]+r'}$'))
    ax.errorbar(xpoints, model3bcorr, **lpd.chmap(lpd.plotstyle_add_point, fmt='--', zorder=-10, label=r'$G_\mathrm{\kappa/T^3=' + str(args.kappa_b) + r', c=' + str(
        args.c_b) + r'}^{\mathrm{max(IR,c UV_{b})}}$'))
    ax.errorbar(xpoints, model3acorr, **lpd.chmap(lpd.plotstyle_add_point, fmt='--', zorder=-11, label=r'$G_\mathrm{\kappa/T^3=' + str(args.kappa_a) + r', c=' + str(
        args.c_a) + r'}^{\mathrm{max(IR,c UV_{a})}}$'))
    ax.errorbar(xpoints, PhiUVa_corr, **lpd.chmap(lpd.plotstyle_add_point, fmt='-', label=r'$G_{\mathrm{UV_a}}$'))
    ax.errorbar(xpoints, PhiUVb_corr, **lpd.chmap(lpd.plotstyle_add_point, fmt='-', label=r'$G_{\mathrm{UV_b}}$'))
    ax.errorbar(xpoints, PhiIRbcorr, **lpd.chmap(lpd.plotstyle_add_point, fmt='-.', label=r'$G_{\Phi_\mathrm{IR_b}}$'))
    ax.errorbar(xpoints, PhiIRacorr, **lpd.chmap(lpd.plotstyle_add_point, fmt='-.', label=r'$G_{\Phi_\mathrm{IR_a}}$'))

    lpd.legendstyle.update(dict(loc="lower right", bbox_to_anchor=(1.01, 0.07)))
    ax.legend(**lpd.chmap(lpd.legendstyle, bbox_to_anchor=(1,0.5), loc='center left'))
    ax.set_title(r'quenched, $1.5 T_c$, IR=$\kappa \omega /2T$', x=0.5, y=1, fontsize=8)

    basepath_plot = "/work/home/altenkort/work/correlators_flow/plots/quenched_1.50Tc_zeuthenFlow/EE/"
    fig.savefig(basepath_plot+"EE_flow_extr_comparison.pdf")

    # ==========================================================================================================================================================

    # /work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/EE/cont_extr/
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
    max_tauT = 0.5
    min_flowradius = 0.05

    flowradiusBytauT_arr = numpy.linspace(min_flowradius/max_tauT, 0.3, 21)
    frankenstein_arr = []
    for f, flowradiusBytauT in enumerate(flowradiusBytauT_arr):
        frankenstein_edata = []
        frankenstein_ydata = []
        frankenstein_xdata = []
        for i, tauT in enumerate(possible_tauTs):
            wanted_flowradius = flowradiusBytauT * tauT
            # flowindex = (numpy.fabs(flowradii - wanted_flowradius)).argmin()
            if cont_int[i] is not None:
                try:
                    frankenstein_edata.append(cont_err_int[i](wanted_flowradius))
                    frankenstein_ydata.append(cont_int[i](wanted_flowradius))
                    frankenstein_xdata.append(tauT)
                except ValueError:
                    pass
        print(flowradiusBytauT, len(frankenstein_xdata))
        frankenstein_arr.append(numpy.column_stack((frankenstein_xdata, frankenstein_ydata, frankenstein_edata)))
        numpy.savetxt("/work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/EE/frankenstein_"+r'{0:.3f}'.format(flowradiusBytauT)+".txt",
                      numpy.column_stack((frankenstein_xdata, frankenstein_ydata, frankenstein_edata)),
                      header="tauT      G/Gnorm_sqrt(8tauF)/tau="+r'{0:.3f}'.format(flowradiusBytauT)+"       err")
        ax2.errorbar(frankenstein_xdata, frankenstein_ydata, frankenstein_edata, zorder=-len(flowradiusBytauT_arr)+f,
                     **lpd.chmap(lpd.plotstyle_add_point, markersize=0, fmt='x-', color=lpd.get_color(flowradiusBytauT_arr, f),
                                 label=r'$'+r'{0:.3f}'.format(flowradiusBytauT)+r'}$'))  # G_{\sqrt{8\tau_F}/\tau='

    # ax2.errorbar(xpoints, model3bcorr, **lpd.chmap(lpd.plotstyle_add_point, fmt='--', zorder=-10, label=r'$G_\mathrm{\kappa/T^3=' + str(args.kappa_b) + r', c=' + str(
    #     args.c_b) + r'}^{\mathrm{max(IR,c UV_{b})}}$'))
    # ax2.errorbar(xpoints, model3acorr, **lpd.chmap(lpd.plotstyle_add_point, fmt='--', zorder=-11, label=r'$G_\mathrm{\kappa/T^3=' + str(args.kappa_a) + r', c=' + str(
    #     args.c_a) + r'}^{\mathrm{max(IR,c UV_{a})}}$'))
    # ax2.errorbar(xpoints, PhiUVa_corr, **lpd.chmap(lpd.plotstyle_add_point, fmt='-', label=r'$G_{\mathrm{UV_a}}$'))
    # ax2.errorbar(xpoints, PhiUVb_corr, **lpd.chmap(lpd.plotstyle_add_point, fmt='-', label=r'$G_{\mathrm{UV_b}}$'))
    # ax2.errorbar(xpoints, PhiIRbcorr, **lpd.chmap(lpd.plotstyle_add_point, fmt='-.', label=r'$G_{\Phi_\mathrm{IR_b}}$'))
    # ax2.errorbar(xpoints, PhiIRacorr, **lpd.chmap(lpd.plotstyle_add_point, fmt='-.', label=r'$G_{\Phi_\mathrm{IR_a}}$'))

    lpd.legendstyle.update(dict(loc="lower right", bbox_to_anchor=(1.01, 0.07))) # labelspacing=1
    ax2.legend(**lpd.chmap(lpd.legendstyle, bbox_to_anchor=(1, 0.5), loc='center left', ncol=2, fontsize=3))
    ax2.set_title(r'quenched, $1.5 T_c$, IR=$\kappa \omega /2T$', x=0.5, y=1, fontsize=8)

    fig2.savefig(basepath_plot + "EE_frankenstein_comparison.pdf")

    # ==========================================================================================================================================================

    if False:
        fig3, ax3, _ = lpd.create_figure(xlims=args.xlims, xlabel=r'$\tau T$', xlabelpos=(0.95, 0.05), ylabelpos=(0.08, 0.9),
                                         UseTex=args.use_tex,
                                         ylabel=None, figsize=(1.5 * (3 + 3 / 8), 1.5 * (3 + 3 / 8) / 16 * 9))

        for i, flowradiusBytauT in enumerate(flowradiusBytauT_arr):
            # spfs
            PhiUVa_corr = []
            PhiUVb_corr = []
            PhiIRacorr = []
            PhiIRbcorr = []
            model3acorr = []
            model3bcorr = []
            for flowindex in frankenstein_arr[i][:, 0]:
                PhiUVa_corr.append(args.c_a*scipy.integrate.quad(lambda OmegaByT: Integrand(OmegaByT, flowindex, PhiUVaByT3), 0, MaxOmegaByT)[0] / Gnorm(flowindex))
                PhiUVb_corr.append(args.c_b*scipy.integrate.quad(lambda OmegaByT: Integrand(OmegaByT, flowindex, PhiUVbByT3), 0, MaxOmegaByT)[0] / Gnorm(flowindex))
                PhiIRacorr.append(scipy.integrate.quad(lambda OmegaByT: Integrand(OmegaByT, flowindex, PhiIR, args.kappa_a), 0, MaxOmegaByT)[0] / Gnorm(flowindex))
                PhiIRbcorr.append(scipy.integrate.quad(lambda OmegaByT: Integrand(OmegaByT, flowindex, PhiIR, args.kappa_b), 0, MaxOmegaByT)[0] / Gnorm(flowindex))
                model3acorr.append(scipy.integrate.quad(lambda OmegaByT: Integrand(OmegaByT, flowindex, model3, args.kappa_a, PhiUVaByT3, args.c_a), 0, MaxOmegaByT)[0] / Gnorm(flowindex))
                model3bcorr.append(scipy.integrate.quad(lambda OmegaByT: Integrand(OmegaByT, flowindex, model3, args.kappa_b, PhiUVbByT3, args.c_b), 0, MaxOmegaByT)[0] / Gnorm(flowindex))

            ax3.errorbar(frankenstein_arr[i][:, 0], frankenstein_arr[i][:, 1]/PhiUVa_corr, frankenstein_arr[i][:, 2]/PhiUVa_corr,
                         **lpd.chmap(lpd.plotstyle_add_point, fmt='-', label=r'$\frac{G_{\sqrt{8\tau_F}/\tau='+r'{0:.3f}'.format(flowradiusBytauT)+r'}}{ c\Phi_\mathrm{UV}^a}$'))

            ax3.errorbar(frankenstein_arr[i][:, 0], frankenstein_arr[i][:, 1] / PhiUVb_corr, frankenstein_arr[i][:, 2] / PhiUVb_corr,
                         **lpd.chmap(lpd.plotstyle_add_point, fmt='-',
                                     label=r'$\frac{G_{\sqrt{8\tau_F}/\tau=' + r'{0:.3f}'.format(flowradiusBytauT) + r'}}{c \Phi_\mathrm{UV}^b}$'))

        lpd.legendstyle.update(dict(loc="lower right", bbox_to_anchor=(1.01, 0.07)))
        ax3.legend(**lpd.chmap(lpd.legendstyle, bbox_to_anchor=(1, 0.5), loc='center left'))
        ax3.set_title(r'quenched, $1.5 T_c$', x=0.5, y=1, fontsize=8)

        fig3.savefig(basepath_plot + "EE_frankenstein_comparison_ratio.pdf")


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
