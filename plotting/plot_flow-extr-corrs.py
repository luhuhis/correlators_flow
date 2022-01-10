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
    PhiUVaByT3 = scipy.interpolate.InterpolatedUnivariateSpline(PhiUVa[0], PhiUVa[1], k=3)
    PhiUVb = numpy.loadtxt("/work/data/htshu/ee_spf/PHIUV_b.dat", unpack=True)[0:3]
    PhiUVbByT3 = scipy.interpolate.InterpolatedUnivariateSpline(PhiUVb[0], PhiUVb[1], k=3)

    MaxOmegaByT = PhiUVa[0][-1]

    # spf = PhiUVa
    PhiUVa_corr = []
    PhiUVb_corr = []
    PhiIRacorr = []
    PhiIRbcorr = []
    model3acorr = []
    model3bcorr = []
    xpoints = numpy.linspace(0, 0.5, args.npoints)
    for x in xpoints:
        PhiUVa_corr.append(scipy.integrate.quad(lambda OmegaByT: Integrand(OmegaByT, x, PhiUVaByT3), 0, MaxOmegaByT)[0] / Gnorm(x))
        PhiUVb_corr.append(scipy.integrate.quad(lambda OmegaByT: Integrand(OmegaByT, x, PhiUVbByT3), 0, MaxOmegaByT)[0] / Gnorm(x))
        PhiIRacorr.append(scipy.integrate.quad(lambda OmegaByT: Integrand(OmegaByT, x, PhiIR, args.kappa_a), 0, MaxOmegaByT)[0] / Gnorm(x))
        PhiIRbcorr.append(scipy.integrate.quad(lambda OmegaByT: Integrand(OmegaByT, x, PhiIR, args.kappa_b), 0, MaxOmegaByT)[0] / Gnorm(x))
        model3acorr.append(scipy.integrate.quad(lambda OmegaByT: Integrand(OmegaByT, x, model3, args.kappa_a, PhiUVaByT3, args.c_a), 0, MaxOmegaByT)[0] / Gnorm(x))
        model3bcorr.append(scipy.integrate.quad(lambda OmegaByT: Integrand(OmegaByT, x, model3, args.kappa_b, PhiUVbByT3, args.c_b), 0, MaxOmegaByT)[0] / Gnorm(x))

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


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
