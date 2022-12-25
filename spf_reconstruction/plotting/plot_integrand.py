#!/usr/bin/env python3

import lib_process_data as lpd
import numpy
import scipy.interpolate
import scipy.integrate
import argparse


def Kernel(OmegaByT, tauT):
    return numpy.cosh(OmegaByT / 2 - OmegaByT * tauT) / numpy.sinh(OmegaByT / 2)


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--outputpath', type=str, help='where to save the plot', required=True)
    parser.add_argument('--PathPhiUV', type=str, required=True)
    args = parser.parse_args()

    fig, ax, _ = lpd.create_figure(ylims=[0.1, 500], xlims=[0, 45],
                                   xlabel="$\\omega/T$", ylabel="$\\rho(\\omega) K(\\omega, \\tau) / T^2$")


    #ylims=[0.001, 100], xlims=[0.1, 100]
    PhiUV = numpy.loadtxt(args.PathPhiUV)
    PhiUV = PhiUV[:, 0:2]
    PhiuvByT3 = scipy.interpolate.InterpolatedUnivariateSpline(PhiUV[:, 0], PhiUV[:, 1], k=3, ext=2)

    OmegaByT = numpy.logspace(-2, 3, 1000)

    tauTs = (0.1, 0.2, 0.3, 0.4, 0.5)

    for tauT in tauTs:
        y = numpy.sqrt((3*OmegaByT)**2+PhiuvByT3(OmegaByT)**2) * Kernel(OmegaByT, tauT)# / OmegaByT
        ax.errorbar(OmegaByT, y, fmt='-', label=str(tauT), lw=1)
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

    #\\frac{\\tau/a}{N_\\tau}=

    ax.legend(title="$\\tau T$", loc='upper right', bbox_to_anchor=(1, 0.84), handlelength=1)

    # ax.axvline(x=numpy.pi, **lpd.verticallinestyle)

    ax.set_yscale('log')
    # ax.set_xscale('log')

    # save figure
    file = args.outputpath + "/integrand.pdf"
    print("saving ", file)
    fig.savefig(file)


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()



