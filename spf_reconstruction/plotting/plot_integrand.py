#!/usr/bin/env python3

import lib_process_data as lpd
import numpy
import scipy.interpolate
import scipy.integrate
import argparse
import matplotlib.pyplot
from spf_reconstruction.model_fitting.EE_UV_spf import get_spf, add_args


def Kernel(OmegaByT, tauT):
    return numpy.cosh(OmegaByT / 2 - OmegaByT * tauT) / numpy.sinh(OmegaByT / 2)


def Gnorm(tauT):
    return numpy.pi ** 2 * (numpy.cos(numpy.pi * tauT) ** 2 / numpy.sin(numpy.pi * tauT) ** 4 + 1 / (3 * numpy.sin(numpy.pi * tauT) ** 2))


def PHIUVLOByT3(OmegaByT):
    return OmegaByT**3/(6*numpy.pi)


def PhiIR(OmegaByT, kappaByT3):
    return kappaByT3 / 2 * OmegaByT


def integrand(args):
    fig, ax, _ = lpd.create_figure(ylims=[0, 7.9], xlims=[0, 45],
                                   xlabel="$\\omega/T$", ylabel=r'$ \frac{1}{\pi} \displaystyle \frac{\rho_{\mathrm{smax}_\mathrm{LO}}}{T^3}  K(\omega,\tau) $')

    # ylims=[0.001, 100], xlims=[0.1, 100]
    # if not args.PathPhiUV:
    #     thisPhiUVByT3 = PHIUVLOByT3
    # else:
    OmegaByT_arr, g2_arr, LO_SPF, NLO_SPF = get_spf(0, args.max_type, "eff", 0.472, 1, args.Npoints, args.Nloop)
    thisPhiUVByT3 = scipy.interpolate.InterpolatedUnivariateSpline(OmegaByT_arr, LO_SPF, k=3, ext=2)
        # PhiUV = numpy.loadtxt(args.PathPhiUV)
        # PhiUV = PhiUV[:, 0:2]
        # thisPhiUVByT3 = scipy.interpolate.InterpolatedUnivariateSpline(PhiUV[:, 0], PhiUV[:, 1], k=3, ext=2)

    OmegaByT = numpy.logspace(-2, 3, 10000)

    tauTs = (0.2, 0.3, 0.4, 0.5)

    for i, tauT in enumerate(tauTs):
        y = numpy.sqrt(PhiIR(OmegaByT, 3)**2 + thisPhiUVByT3(OmegaByT)**2) * Kernel(OmegaByT, tauT) / numpy.pi  # / OmegaByT
        y_UV = thisPhiUVByT3(OmegaByT) * Kernel(OmegaByT, tauT) / numpy.pi  # / OmegaByT
        ax.errorbar(OmegaByT, y, fmt='-', label=str(tauT), lw=1, color=lpd.get_discrete_color(i), zorder=-10)
        ax.errorbar(OmegaByT, y_UV, fmt='--', lw=0.5, zorder=-100, color=lpd.get_discrete_color(i))
    # for maxomega in (0.1, 1, 10, 100):
    #     corr = []
    #     for tauT in tauTs:
    #         y = Kernel(OmegaByT, tauT)
    #
    #         spline = scipy.interpolate.InterpolatedUnivariateSpline(OmegaByT, y*OmegaByT, k=3, ext=2, check_finite=True)
    #         corr.append(scipy.integrate.quad(spline, OmegaByT[0], maxomega)[0]/Gnorm(tauT))  #
    #         # ax.errorbar(OmegaByT, OmegaByT * y, fmt='-', label=str(tauT), lw=0.5)
    #
    #     ax.errorbar(tauTs, corr, fmt='-', label=str(maxomega), lw=0.5)
    #     # ax.errorbar(tauTs, Gnorm(tauTs), fmt='-', label="corr2", lw=0.5)

    # \\frac{\\tau/a}{N_\\tau}=

    ax.legend(title="$\\tau T$", loc='upper right', bbox_to_anchor=(1, 1), handlelength=0.9)

    # ax.axvline(x=numpy.pi, **lpd.verticallinestyle)

    # ax.set_yscale('log')
    # ax.set_xscale('log')

    # save figure
    file = args.outputpath + "/integrand.pdf"
    print("saving ", file)
    fig.savefig(file)
    matplotlib.pyplot.close(fig)


def model_corrs(args):

    tauTs = numpy.linspace(0.015, 0.5, 100)
    Gmodel_norm = []
    Gmodel_LO_UV = []
    Gmodel_smax_LO_UV_kappa1 = []
    Gmodel_smax_LO_UV_kappa2 = []
    Gmodel_smax_LO_UV_kappa3 = []
    OmegaByT_arr, g2_arr, LO_SPF, NLO_SPF = get_spf(0, args.max_type, "eff", 0.472, 1, args.Npoints, args.Nloop)

    for tauT in tauTs:
        tmp = scipy.integrate.quad(lambda OmegaByT: PHIUVLOByT3(OmegaByT) * Kernel(OmegaByT, tauT)/numpy.pi, 0, 1000)[0]
        Gmodel_norm.append(tmp)

        PhiUVByT3_interpolation = scipy.interpolate.InterpolatedUnivariateSpline(OmegaByT_arr, LO_SPF, k=3, ext=2)
        tmp = scipy.integrate.quad(lambda OmegaByT: PhiUVByT3_interpolation(OmegaByT) * Kernel(OmegaByT, tauT)/numpy.pi, 0, 1000)[0]
        Gmodel_LO_UV.append(tmp)

        tmp = scipy.integrate.quad(lambda OmegaByT: numpy.sqrt(PhiIR(OmegaByT, 1)**2 + PhiUVByT3_interpolation(OmegaByT)**2) * Kernel(OmegaByT, tauT)/numpy.pi, 0, 1000)[0]
        Gmodel_smax_LO_UV_kappa1.append(tmp)

        tmp = \
        scipy.integrate.quad(lambda OmegaByT: numpy.sqrt(PhiIR(OmegaByT, 2) ** 2 + PhiUVByT3_interpolation(OmegaByT) ** 2) * Kernel(OmegaByT, tauT) / numpy.pi,
                             0, 1000)[0]
        Gmodel_smax_LO_UV_kappa2.append(tmp)

        tmp = \
        scipy.integrate.quad(lambda OmegaByT: numpy.sqrt(PhiIR(OmegaByT, 3) ** 2 + PhiUVByT3_interpolation(OmegaByT) ** 2) * Kernel(OmegaByT, tauT) / numpy.pi,
                             0, 1000)[0]
        Gmodel_smax_LO_UV_kappa3.append(tmp)



    fig, ax, _ = lpd.create_figure(ylims=[1, 3.6], xlims=[0, 0.5],
                                   xlabel=r'$\tau T$',
                                   ylabel=r'$\displaystyle \frac{G^{\mathrm{smax}_\mathrm{LO}}}{G^\mathrm{norm}}$')
    ax.set_xticks((0.0, 0.1, 0.2, 0.3, 0.4, 0.5))


    Gnorms = Gnorm(tauTs)

    # ax.errorbar(tauTs, Gmodel_norm/Gnorms, label=r'$\frac{\phi_\mathrm{UV}^\mathrm{LO}}{C_\mathrm{F} g^2}$', lw=0.5, fmt='--', color='k')
    ax.errorbar(tauTs, Gmodel_LO_UV / Gnorms, fmt='-', label=r'$0$', lw=1, color='C0')
    ax.errorbar(tauTs, Gmodel_smax_LO_UV_kappa1 / Gnorms, fmt='-', label=r'$1$', lw=1, color='C1')
    ax.errorbar(tauTs, Gmodel_smax_LO_UV_kappa2 / Gnorms, fmt='-', label=r'$2$', lw=1, color='C2')
    ax.errorbar(tauTs, Gmodel_smax_LO_UV_kappa3 / Gnorms, fmt='-', label=r'$3$', lw=1, color='C3')

    # ax.errorbar(0, 0, label=r'$\kappa_E/T^3$', markersize=0, alpha=0, lw=0)

    # reverse ordering of legend
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1], title=r'$\kappa_E/T^3$', loc='lower right', bbox_to_anchor=(1, 0.1), handlelength=0.9, framealpha=0)

    # save figure
    file = args.outputpath + "/model_corrs.pdf"
    print("saving ", file)
    fig.savefig(file)
    matplotlib.pyplot.close(fig)

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--outputpath', type=str, help='where to save the plot', required=True)
    parser.add_argument('--PathPhiUV', type=str)
    add_args(parser)
    args = parser.parse_args()

    integrand(args)

    model_corrs(args)




if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()



