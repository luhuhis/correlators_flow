#!/usr/bin/env python3
import scipy.interpolate

import lib_process_data as lpd
import numpy as np
import argparse
import compute_UV_spf
import matplotlib
import spf_reconstruct

from nptyping import NDArray, Float64
from typing import Literal as Shape

from correlator_analysis.double_extrapolation.BB_renormalization.compute_Zf2 import get_scale_choices

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--outputpath_plot', type=str)
    parser.add_argument('--inputfolder', type=str)
    # compute_UV_spf.add_args(parser)
    args = parser.parse_args()
    return args


@lpd.typed_frozen_data
class InputParams:
    corr: str
    order: str
    min_scale: str
    omega_prefactor: str
    max_type: str
    plotlabel: str
    fmt: str


@lpd.typed_frozen_data
class UVSpfData:
    OmegaByT_arr: NDArray[Shape["*"], Float64]
    g2_arr: NDArray[Shape["*"], Float64]
    PhiUVByT3: NDArray[Shape["*"], Float64]
    corr: NDArray[Shape["*"], Float64]
    input_params: InputParams


def twoRhoByomegaTsq(omegaByT, rhoByT3):
    return 2 * rhoByT3 / omegaByT


class SpfPlotter:
    def __init__(self, args, uv_spf_data_list):
        self.uv_spf_data_list = uv_spf_data_list
        self.outputpath_plot = args.outputpath_plot
        self.inputfolder = args.inputfolder

    def __setup_plot(self):
        self.fig, self.ax1, _ = lpd.create_figure(xlabel=r'$\omega/T$',
                                                  ylabel=r'$ \frac{2 \rho_\mathrm{UV}}{\omega T^2}$',
                                                  subplot=121, figsize=(14, 7))
        self.ax1.set_ylim((0.01, 1000))
        self.ax1.set_xlim((0.1, 100))
        self.ax1.set_yscale('log')
        self.ax1.set_xscale('log')

        self.ax2 = self.fig.add_subplot(122)
        ax = lpd.apply_ax_settings(self.ax2, (0, 0.51), (2, 4.5), r'$\tau T$', r'$G^\text{UV}\frac{G_B(0.25)}{G^\text{UV}(0.25)}$')

    def _plot_BB(self):
        scale_choices = get_scale_choices()
        for scale_choice in scale_choices:
            if scale_choice.choice_label == "ref4.0_UVLO_IRLO":
                BB = np.loadtxt(self.inputfolder + "/BB/BB_flow_extr_relflow_" + scale_choice.choice_label + ".txt", unpack=True)
                self.ax2.errorbar(BB[0], BB[1], BB[2])
                self.norm_value = BB[1][8]

    def _plot(self):
        self.__setup_plot()

        # TODO also plot our BB correlator data here

        title = r'Order, $\bar{\mu}_\omega$'
        self.ax1.errorbar(1, 1, fmt='.', markersize=0, label=title)
        self.ax2.errorbar(1, 1, fmt='.', markersize=0, label=title)
        self.ax1.set_prop_cycle(None)
        self.ax2.set_prop_cycle(None)

        self._plot_BB()

        for uv_spf_data in self.uv_spf_data_list:
            self.ax1.errorbar(uv_spf_data.OmegaByT_arr, twoRhoByomegaTsq(uv_spf_data.OmegaByT_arr, uv_spf_data.PhiUVByT3),
                              fmt=uv_spf_data.input_params.fmt, label=uv_spf_data.input_params.plotlabel)

            self.ax2.errorbar(np.arange(1/36, 0.501, 1/36), uv_spf_data.corr/uv_spf_data.corr[8]*self.norm_value, label=uv_spf_data.input_params.plotlabel, fmt=uv_spf_data.input_params.fmt)

        self.ax1.axvline(x=2*np.pi, **lpd.verticallinestyle)
        self.ax1.axvline(x=19.18, **lpd.verticallinestyle)

        self.__finalize_plot()

    def __finalize_plot(self):
        self.ax1.legend(loc="center left", bbox_to_anchor=(0, 0.7), fontsize=6, title_fontsize=6, framealpha=0, handlelength=1)
        self.ax2.legend(loc="upper right", bbox_to_anchor=(1, 0.5), fontsize=6, title_fontsize=6, framealpha=0, handlelength=1)
        file = self.outputpath_plot+"/UV_spf.pdf"
        self.fig.savefig(file)
        print("saved", file)
        matplotlib.pyplot.close(self.fig)

    @classmethod
    def plot(cls, *args):
        # This method only exists such that you don't have to manually instantiate the class.
        instance = cls(*args)
        instance._plot()


def compute_correlator(OmegaByT_arr, SpfByT3_arr):
    def calc_corr(OmegaByT_arr, SpfByT3_arr, tauT):
        integrand = scipy.interpolate.InterpolatedUnivariateSpline(OmegaByT_arr,
                                                                   [Integrand(OmegaByT, SpfByT3, tauT)
                                                                    for OmegaByT, SpfByT3 in
                                                                    zip(OmegaByT_arr, SpfByT3_arr)],
                                                                   k=3, ext=2)
        corr = scipy.integrate.quad(integrand, OmegaByT_arr[0], OmegaByT_arr[-1])[0] / spf_reconstruct.Gnorm(tauT)
        return corr

    def Integrand(OmegaByT, SpfByT3, tauT):
        return 1. / np.pi * spf_reconstruct.Kernel(OmegaByT, tauT) * SpfByT3

    Nt = 36

    correlator = np.asarray([calc_corr(OmegaByT_arr, SpfByT3_arr, tauT) for tauT in np.arange(1/Nt, 0.501, 1/Nt)])

    return correlator


def get_spf_and_corr(Nf, T_in_GeV, Npoints, Nloop, input_params):
    OmegaByT_arr, g2_arr, SpfByT3_arr = compute_UV_spf.get_spf(Nf, input_params.max_type, input_params.min_scale, T_in_GeV, input_params.omega_prefactor, Npoints, Nloop, input_params.order, input_params.corr)

    corr = compute_correlator(OmegaByT_arr, SpfByT3_arr)

    return OmegaByT_arr, g2_arr, SpfByT3_arr, corr, input_params


def parallelwrapper(input_param, Nf, T_in_GeV, Npoints, Nloop):
    return UVSpfData(*get_spf_and_corr(Nf, T_in_GeV, Npoints, Nloop, input_param))


def main():

    args = parse_args()

    Nf = 0
    T_in_GeV = 0.472
    Npoints = 1000
    Nloop = 5
    corr = "BB"

    input_params = [InputParams(corr, "LO", "2piT",       "1", "hard", r'LO,\ \ \  $\sqrt{\omega^2 + (2 \pi T)^2}$',    '-'),
                    InputParams(corr, "LO", "mu_IR_NLO",  "1", "hard", r'LO,\ \ \   $\sqrt{\omega^2 + (19.18 T)^2}$',    '-'),
                    InputParams(corr, "NLO", "2piT",      "2", "hard", r'NLO,   $\sqrt{(2\omega)^2 + (2 \pi T)^2}$', '-'),
                    InputParams(corr, "NLO", "mu_IR_NLO", "2", "hard", r'NLO,   $\sqrt{(2\omega)^2 + (19.18 T)^2}$', '-'),
                    InputParams(corr, "LO", "2piT",       "1", "smooth", r'LO,\ \ \ $\mathrm{max}(\omega, 2 \pi T)$',    '--'),
                    InputParams(corr, "LO", "mu_IR_NLO",  "1", "smooth", r'LO,\ \ \ $\mathrm{max}(\omega, 19.18 T)$',    '--'),
                    InputParams(corr, "NLO", "2piT",      "2", "smooth", r'NLO, $\mathrm{max}(2\omega, 2 \pi T)$',   '--'),
                    InputParams(corr, "NLO", "mu_IR_NLO", "2", "smooth", r'NLO, $\mathrm{max}(2\omega, 19.18 T)$',   '--'),
                    ]

    uv_spf_data = lpd.parallel_function_eval(parallelwrapper, input_params, 10, Nf, T_in_GeV, Npoints, Nloop)

    SpfPlotter.plot(args, uv_spf_data)

    return


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
