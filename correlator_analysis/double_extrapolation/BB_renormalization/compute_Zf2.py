#!/usr/bin/env python3
import numpy
from scipy import integrate, interpolate
import argparse
import lib_process_data as lpd

# Type hints
from nptyping import NDArray, Float64
from typing import Literal as Shape



fmts = ['-', '-', '--', ':']
fontsize = 8
gamma_0 = 3/(8 * numpy.pi**2)


def find_index(array, value):
    array = numpy.asarray(array)
    return (numpy.abs(array - value)).argmin()


def plot_integrand(args, data):  #, reference_muF_by_T, flow_extr_window):
    fig, ax, _ = lpd.create_figure(xlabel=r'$\mu_\mathrm{F}/T$', ylabel=r'$ \frac{T}{\mu_\mathrm{F}} \displaystyle\gamma_0 g^2$')
    ax.set_ylim((-0.06, 0.005))
    ax.set_xlim((0, 21))

    counter = 0
    ax.errorbar(data.mu_by_T, data.integrand, fmt=fmts[counter % 4], zorder=counter)
    counter += 1
    # ax.axvline(x=reference_muF_by_T, **lpd.verticallinestyle)
    # ax.fill_between(flow_extr_window, [-250, -250], [400, 400], facecolor='grey', alpha=0.25, zorder=-1000)
    ax.legend(title=r'$g^2$', fontsize=fontsize, title_fontsize=fontsize, framealpha=0)
    file = args.outputpath_plot+"/integrand.pdf"
    fig.savefig(file)
    print("saved", file)


def convert_sqrt8taufT_to_tfT2(sqrt8taufT):
    return sqrt8taufT**2/8


# def plot_Zf2(args, data, reference_tauF_T2, flow_extr_window):
#     fig2, ax2, plots2 = lpd.create_figure(xlabel=r'$\mu_\mathrm{F}/T$', ylabel=r'$Z_f^2$')
#
#     counter = 0
#     for _, value in data.items():
#         muF_by_T, Z2, plotlabel = value
#         ax2.errorbar(muF_by_T, Z2, label=plotlabel, fmt=fmts[counter % 4], zorder=counter)
#         counter += 1
#
#     ax2.axvline(x=reference_tauF_T2, **lpd.verticallinestyle)
#     ax2.fill_between(flow_extr_window,
#                      [-1, -1], [2, 2], facecolor='grey', alpha=0.25, zorder=-1000)
#     ax2.axhline(y=1, **lpd.horizontallinestyle)
#     ax2.set_ylim((0.7, 1.3))
#     ax2.set_xlim(0, 50)
#     ax2.legend(title=r'$g^2$', fontsize=fontsize, title_fontsize=fontsize, framealpha=0)
#     if args.eight:
#         eight_suffix = "_8"
#     else:
#         eight_suffix = ""
#     file = args.outputpath_plot+"Zf2"+eight_suffix+".pdf"
#     print("saving", file)
#     fig2.savefig(file)


def compute_Z_run(coupling_container, scale_choice):

    def integrand(muBar_By_T):
        return 2*gamma_0*coupling_container.g2_spline(muBar_By_T)

    def compute_Z_run(muF_by_T):
        muBarUV_By_T = scale_choice.muBarUV_by_muF * muF_by_T
        muBarIR_By_T = scale_choice.muBarIR_by_T
        this_Z_run = integrate.quad(integrand, muBarUV_By_T, muBarIR_By_T)[0]
        return this_Z_run

    integral = numpy.asarray([compute_Z_run(muF_by_T) for muF_by_T in coupling_container.mu_by_T])
    Z_run = numpy.exp(integral)

    return Z_run


def calc_Zf2(coupling_container, scale_choice):

    ref_index = find_index(coupling_container.mu_by_T, scale_choice.muBarIR_by_T)

    Z2 = []
    for muF_by_T in coupling_container.mu_by_T:
        target_index = find_index(coupling_container.mu_by_T, muF_by_T)

        if ref_index < target_index:
            this_y = coupling_container.integrand[ref_index:target_index+1]
            this_x = coupling_container.mu_by_T[ref_index:target_index+1]
            sign = 1
        else:
            this_y = coupling_container.integrand[target_index:ref_index+1]
            this_x = coupling_container.mu_by_T[target_index:ref_index+1]
            sign = -1  # account for reversed integral bounds

        integral = integrate.trapz(this_y, this_x)

        if integral < 0:
            print("Warn: integral should be positive", muF_by_T, ref_index, target_index, integral)
        this_Z2 = numpy.exp(sign * integral)
        Z2.append(this_Z2)

    return numpy.asarray(Z2)


@lpd.typed_frozen_data
class CouplingContainer:
    mu_by_T: NDArray[Shape["*"], Float64]
    g2: NDArray[Shape["*"], Float64]
    g2_spline: interpolate._fitpack2.InterpolatedUnivariateSpline
    integrand: NDArray[Shape["*"], Float64]


def load_data(args: argparse.Namespace) -> CouplingContainer:
    muF_by_T, g2 = numpy.loadtxt(args.g2_file, unpack=True)[:2]
    g2_spline = interpolate.InterpolatedUnivariateSpline(muF_by_T, g2, k=3, ext=2)
    integrand = 2 * gamma_0 * g2 / muF_by_T
    return CouplingContainer(muF_by_T, g2, g2_spline, integrand)


def compute_Zk(coupling_container, scale_choice):
    muF_by_T = coupling_container.mu_by_T
    muBarUV_by_T = scale_choice.muBarUV_by_muF * muF_by_T
    g2_MSBar = coupling_container.g2_spline(muF_by_T)  # TODO muBarUV_by_T ?
    inner_bracket = numpy.log(muBarUV_by_T ** 2 / coupling_container.mu_by_T ** 2) - 2 * numpy.log(2) - numpy.euler_gamma - 8 / 3
    Zk = 1/(1 + gamma_0 * g2_MSBar * inner_bracket)
    return Zk



def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('--g2_file', type=str, required=True, help="path to file containing g^2 in flow scheme")
    parser.add_argument('--outputpath_plot', default="/work/home/altenkort/work/correlators_flow/plots/quenched_1.50Tc_zeuthenFlow/coupling/")
    parser.add_argument('--outputpath_data', default="/work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/coupling/")

    args = parser.parse_args()

    return args


@lpd.typed_frozen_data
class ScaleChoice:
    muRef_by_T: float
    muBarUV_by_muF: float
    muBarIR_by_T: float
    choice_label: str



def get_scale_choices():
    scale_choices = []

    muBarIR_by_T_choices = [2 * numpy.pi, 4 * numpy.pi * numpy.exp(1 - numpy.euler_gamma)]
    muBarUV_by_muF_choices = [1., numpy.sqrt(4 * numpy.exp(numpy.euler_gamma + 8 / 3))]
    order_string = ["LO", "NLO"]
    muRef_by_T_choices = [8.,]

    for UV in range(2):  # UV
        for IR in range(2):
            for ref in range(1):
                choice_label = "ref=" + str(muRef_by_T_choices[ref])+"_UV="+order_string[UV] + "_IR=" + order_string[IR]
                muBarUV_by_muF_choice = muBarUV_by_muF_choices[UV]
                muBarIR_by_T_choice = muBarIR_by_T_choices[IR]
                muRef_by_T = muRef_by_T_choices[ref]

                scale_choices.append(ScaleChoice(muRef_by_T, muBarUV_by_muF_choice, muBarIR_by_T_choice, choice_label))
    return scale_choices


def compute_Z_phys(coupling_container, scale_choice):
    inner_bracket = numpy.log(scale_choice.muBarIR_by_T**2 / (numpy.pi*4)**2 - 2 + 2*numpy.euler_gamma)
    Z_phys = 1 + gamma_0 * coupling_container.g2_spline(scale_choice.muBarIR_by_T) * inner_bracket
    return Z_phys


def main():

    args = parse_args()
    coupling_container = load_data(args)
    scale_choices = get_scale_choices()

    print(coupling_container.mu_by_T[0], coupling_container.mu_by_T[-1])

    plot_integrand(args, coupling_container)  # , reference_muF_by_T, flow_extr_window)

    for scale_choice in scale_choices:
        Z_K = compute_Zk(coupling_container, scale_choice)
        Z_run = compute_Z_run(coupling_container, scale_choice)
        Z_phys = compute_Z_phys(coupling_container, scale_choice)

        print([*scale_choice])

        print(lpd.format_float(Z_K[0]),
              lpd.format_float(Z_K[-1]),
              lpd.format_float(Z_run[0]),
              lpd.format_float(Z_run[-1]),
              lpd.format_float(Z_phys))

        Z_total = Z_K * Z_run * Z_phys

        print(lpd.format_float(Z_total[0]), lpd.format_float(Z_total[-1]))

        # taufT2 = (1 / coupling_container.mu_by_T) ** 2 / 8
        #
        # lpd.save_columns_to_file(args.outputpath_data + "/Z2_muFByT_" + scale_choice.choice_label + ".dat", (coupling_container.mu_by_T, Z2), ["mu_F/T", "Z_f"])
        # lpd.save_columns_to_file(args.outputpath_data + "/Z2_taufT2_" + scale_choice.choice_label + ".dat", (taufT2, Z2), ["tau_F T^2", "Z_f"])



    # plot_Zf2(args, data_Zf, reference_muF_by_T, flow_extr_window)


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
