#!/usr/bin/env python3
import numpy
from scipy import integrate
import argparse
import lib_process_data as lpd

fmts = ['-', '-', '--', ':']
fontsize = 8

def find_index(array, value):
    array = numpy.asarray(array)
    return (numpy.abs(array - value)).argmin()


def plot_integrand(args, data, reference_muF_by_T):
    fig, ax, _ = lpd.create_figure(xlabel=r'$\mu_\mathrm{F}/T$', ylabel=r'$ \frac{T}{\mu_\mathrm{F}} \displaystyle\gamma_0 g^2$')
    ax.set_ylim((-0.06, 0.005))
    ax.set_xlim((0, 21))

    counter = 0
    for _, value in data.items():
        muF_by_T, integrand, plotlabel = value
        ax.errorbar(muF_by_T, integrand, fmt=fmts[counter % 4], label=plotlabel, zorder=counter)
        counter += 1
    ax.axvline(x=reference_muF_by_T, **lpd.verticallinestyle)
    ax.fill_between([6.66, 16], [-250, -250], [400, 400], facecolor='grey', alpha=0.25, zorder=-1000)
    ax.legend(title=r'$g^2$', fontsize=fontsize, title_fontsize=fontsize, framealpha=0)
    file = args.outputpath_plot+"/integrand.pdf"
    fig.savefig(file)
    print("saved", file)


def convert_sqrt8taufT_to_tfT2(sqrt8taufT):
    return sqrt8taufT**2/8


def plot_Zf2(args, data, reference_tauF_T2):
    fig2, ax2, plots2 = lpd.create_figure(xlabel=r'$\mu_\mathrm{F}/T$', ylabel=r'$Z_f^2$')

    # min_idx = numpy.abs(numpy.asarray(taufT2)-convert_sqrt8taufT_to_tfT2(0.25*0.25)).argmin()-1
    # max_idx = numpy.abs(numpy.asarray(taufT2)-convert_sqrt8taufT_to_tfT2(0.5*0.3)).argmin()+1

    counter = 0
    for _, value in data.items():
        muF_by_T, Z2, plotlabel = value
        ax2.errorbar(muF_by_T, Z2, label=plotlabel, fmt=fmts[counter % 4], zorder=counter)
        counter += 1

    ax2.axvline(x=reference_tauF_T2, **lpd.verticallinestyle)
    ax2.fill_between([6.66, 16],
                     [-1, -1], [1.25, 1.25], facecolor='grey', alpha=0.25, zorder=-1000)
    ax2.axhline(y=1, **lpd.horizontallinestyle)
    ax2.set_ylim((0.95,  1.225))
    ax2.set_xlim(0, 20)
    ax2.legend(title=r'$g^2$', fontsize=fontsize, title_fontsize=fontsize, framealpha=0)
    file = args.outputpath_plot+"Zf2.pdf"
    print("saving", file)
    fig2.savefig(file)


def calc_Zf(args, muF_by_Ts, integrand, reference_muF_by_T, label):

    lowindex = find_index(muF_by_Ts, reference_muF_by_T)

    Z2 = []
    x = []
    for muF_by_T in muF_by_Ts:
        highindex = find_index(muF_by_Ts, muF_by_T)

        if highindex < lowindex:
            this_y = integrand[highindex:lowindex+1]
            this_x = muF_by_Ts[highindex:lowindex+1]
            sign = -1
        else:
            this_y = integrand[lowindex:highindex+1]
            this_x = muF_by_Ts[lowindex:highindex+1]
            sign = 1

        integral = integrate.trapz(this_y, this_x)
        Z2.append(numpy.exp(sign * integral))
        x.append(muF_by_T)

    muFbyT = numpy.asarray(x)
    taufT2 = 1/(muFbyT**2 / 8)

    lpd.save_columns_to_file(args.outputpath_data + "/Z2_muFByT_"+label+".dat", (x, Z2), ["mu_F/T", "Z_f"])
    lpd.save_columns_to_file(args.outputpath_data + "/Z2_taufT2_" + label + ".dat", (taufT2, Z2), ["tau_F T^2", "Z_f"])

    return x, Z2


T_in_GeV = 0.472


def load_data(file):
    muF_by_T, g2 = numpy.loadtxt(file, unpack=True)[:2]
    integrand = (-3*g2)/(8*numpy.pi**2) / muF_by_T
    return muF_by_T, g2, integrand


def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('--T_by_Tc', type=float, required=True)
    parser.add_argument('--g2_files', type=str, nargs='*', required=True, help="path to files containing g^2 in flow scheme")
    parser.add_argument('--filelabels', type=str, nargs='*', required=True, help="file labels")
    parser.add_argument('--plotlabels', type=str, nargs='*', required=True, help="plot labels")
    parser.add_argument('--outputpath_plot', default="/work/home/altenkort/work/correlators_flow/plots/quenched_1.50Tc_zeuthenFlow/coupling/")
    parser.add_argument('--outputpath_data', default="/work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/coupling/")
    args = parser.parse_args()

    return args


def main():
    reference_muF_by_T = 19.179

    args = parse_args()
    data_Zf = {}
    data_integrand = {}

    for file, filelabel, plotlabel in zip(args.g2_files, args.filelabels, args.plotlabels):
        muF_by_T, g2, integrand = load_data(file)
        data_integrand[filelabel] = (muF_by_T, integrand, plotlabel)
        muF_by_T, Zf2 = calc_Zf(args, muF_by_T, integrand, reference_muF_by_T, filelabel)
        data_Zf[filelabel] = (muF_by_T, Zf2, plotlabel)

    plot_integrand(args, data_integrand, reference_muF_by_T)
    plot_Zf2(args, data_Zf, reference_muF_by_T)


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
