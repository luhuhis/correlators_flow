#!/usr/bin/env python3
import numpy
from scipy import integrate
import argparse
import lib_process_data as lpd


def find_index(array, value):
    array = numpy.asarray(array)
    return (numpy.abs(array - value)).argmin()


def plot_integrand(args, tauFT2, integrand, integrand_err, reference_tauF_T2):
    fig, ax, _ = lpd.create_figure(xlabel=r'$\tau_F T^2$', ylabel=r'$ \frac{1}{\tau_FT^2} \left(\frac{-3g^2(\tau_F)}{8\pi^2}\right), \quad T=1.5T_c$')
    ax.set_ylim((-250, 0))
    ax.errorbar(tauFT2, integrand, fmt='-')
    ax.axvline(x=reference_tauF_T2, **lpd.verticallinestyle)
    ax.fill_between([convert_sqrt8taufT_to_tfT2(0.25*0.25), convert_sqrt8taufT_to_tfT2(0.5*0.3)],
                     [-250, -250], [400, 400], facecolor='grey', alpha=0.25, zorder=-1000)
    file = args.outputpath_plot+"/integrand.pdf"
    print("save", file)
    fig.savefig(file)


def convert_sqrt8taufT_to_tfT2(sqrt8taufT):
    return sqrt8taufT**2/8


def plot_Zf2(args, taufT2, Z2, reference_tauF_T2):
    fig2, ax2, plots2 = lpd.create_figure(xlabel=r'$\tau_F T^2$', ylabel=r'$Z_f^2$')

    min_idx = numpy.abs(numpy.asarray(taufT2)-convert_sqrt8taufT_to_tfT2(0.25*0.25)).argmin()-1
    max_idx = numpy.abs(numpy.asarray(taufT2)-convert_sqrt8taufT_to_tfT2(0.5*0.3)).argmin()+1

    ax2.errorbar(taufT2[min_idx:max_idx], Z2[min_idx:max_idx], fmt='-')
    ax2.axvline(x=reference_tauF_T2, **lpd.verticallinestyle)
    # ax2.fill_between([convert_sqrt8taufT_to_tfT2(0.25*0.25), convert_sqrt8taufT_to_tfT2(0.5*0.3)],
    #                  [-1, -1], [1.25, 1.25], facecolor='grey', alpha=0.25, zorder=-1000)
    ax2.axhline(y=1, color='k', zorder=-10000, lw=0.5, dashes=(6,2))
    ax2.set_ylim((0.85,  1.35))
    ax2.set_xlim(0, 0.0045)
    file = args.outputpath_plot+"Zf2.pdf"
    print("saving", file)
    fig2.savefig(file)


def calc_Zf(args, tauFT2, integrand, reference_tauF_T2):

    lowindex = find_index(tauFT2, reference_tauF_T2)

    Z2 = []
    x = []
    for tfT2 in tauFT2:
        highindex = find_index(tauFT2, tfT2)

        if highindex < lowindex:
            this_y = integrand[highindex:lowindex+1]
            this_x = tauFT2[highindex:lowindex+1]
            sign = -1
        else:
            this_y = integrand[lowindex:highindex+1]
            this_x = tauFT2[lowindex:highindex+1]
            sign = 1
        # TODO: DO BOOTSTRAP FOR THE ERRORS ?
        integral = integrate.trapz(this_y, this_x)
        Z2.append(numpy.exp(sign * integral))
        x.append(tfT2)

    file = args.outputpath_data + "/Z2_cont.dat"
    print("saving", file)
    numpy.savetxt(file, numpy.stack((x, Z2), axis=-1), header="tf T^2        Z^2")

    return x, Z2


def load_data(args):
    tauF_by_r0sq, g2, g2_err, _, _, _, _ = numpy.loadtxt(args.g2_file, unpack=True)
    tauFT2 = tauF_by_r0sq*(0.7457*args.T_by_Tc)**2
    integrand = (-3*g2)/(8*numpy.pi**2) / tauFT2
    integrand_err = numpy.fabs((-3*g2_err)/(8*numpy.pi**2) / tauFT2)

    return tauFT2, g2, g2_err, integrand, integrand_err


def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('--T_by_Tc', type=float, required=True)
    parser.add_argument('--g2_file', type=str, required=True, help="path to file containing g^2 in flow scheme")
    parser.add_argument('--outputpath_plot', default="/work/home/altenkort/work/correlators_flow/plots/quenched_1.50Tc_zeuthenFlow/coupling/")
    parser.add_argument('--outputpath_data', default="/work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/coupling/")
    args = parser.parse_args()

    return args


def main():
    reference_tauF_T2 = 0.002719

    args = parse_args()
    tauFT2, g2, g2_err, integrand, integrand_err = load_data(args)
    plot_integrand(args, tauFT2, integrand, integrand_err, reference_tauF_T2)
    tauFT2, Zf2 = calc_Zf(args, tauFT2, integrand, reference_tauF_T2)
    plot_Zf2(args, tauFT2, Zf2, reference_tauF_T2)


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
