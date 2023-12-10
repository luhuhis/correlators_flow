#!/usr/bin/env python3
import numpy as np
from scipy import integrate, interpolate
import argparse
import lib_process_data as lpd
import matplotlib.pyplot
from correlator_analysis.double_extrapolation.BB_renormalization.extrapolate_coupling import ScaleFunctions

# Type hints
from nptyping import NDArray, Float64
from typing import Literal as Shape

gamma_0 = 3 / (8 * np.pi ** 2)


class ZIndividualPlotter:
    def __init__(self, args, Z_container, fontsize, scale_choice):
        self.Z_container = Z_container
        self.outputpath_plot = args.outputpath_plot
        self.fontsize = fontsize
        self.scale_choice = scale_choice

    def __setup_plot(self):
        self.fig, self.ax, _ = lpd.create_figure(xlabel=r'$\mu_\mathrm{F}/T$',
                                                 ylabel=r'$ Z$')
        self.ax.set_ylim((0.5, 2))
        self.ax.set_xlim((0, 21))

    def __plot_Z(self):
        self.ax.errorbar(self.Z_container.mu_By_T, self.Z_container.Z_K, fmt='--', label="$Z_K$")
        self.ax.errorbar(self.Z_container.mu_By_T, self.Z_container.Z_run, fmt=':', label=r'$Z_\text{run}$')
        # Retrieve the current x-axis limits
        xlim = self.ax.get_xlim()
        axis_range = xlim[1] - xlim[0]

        # Calculate fractional positions for xmin and xmax
        xmin_fraction = (self.Z_container.mu_By_T[0] - xlim[0]) / axis_range
        xmax_fraction = (self.Z_container.mu_By_T[-1] - xlim[0]) / axis_range

        self.ax.axhline(y=self.Z_container.Z_phys, xmin=xmin_fraction, xmax=xmax_fraction,
                        dashes=(1, 1), color='C4', label=r'$Z_\text{phys}$')
        self.ax.errorbar(self.Z_container.mu_By_T, self.Z_container.Z_total, fmt='-', label=r'$Z$', zorder=-1)

    def _plot(self):
        self.__setup_plot()
        self.__plot_Z()
        self.__finalize_plot()

    def __finalize_plot(self):
        self.ax.legend(fontsize=self.fontsize, title_fontsize=self.fontsize, framealpha=0, handlelength=4)
        file = self.outputpath_plot+"/Z"+self.scale_choice.choice_label+".pdf"
        self.fig.savefig(file)
        print("saved", file)
        matplotlib.pyplot.close(self.fig)

    @classmethod
    def plot(cls, *args):
        # This method only exists such that you don't have to manually instantiate the class.
        instance = cls(*args)
        instance._plot()

global ylims
ylims = (0.5, 1.65)

class ZTotalPlotter:
    def __init__(self, args, Z_containers, fontsize, scale_choices, output_suffix=""):
        self.Z_containers = Z_containers
        self.outputpath_plot = args.outputpath_plot
        self.fontsize = fontsize
        self.scale_choices = scale_choices
        self.output_suffix = output_suffix

    def _setup_plot(self):
        self.fig, self.ax, _ = lpd.create_figure(xlabel=r'$\mu_\mathrm{F}/T$',
                                                 ylabel=r'$ Z$')
        global ylims
        self.ax.set_ylim(ylims)
        self.ax.set_xlim((0, 17))

    def _plot_grey_flow_extr_band(self):
        self.ax.fill_between([6.66, 16],
                             [-1, -1], [100, 100], facecolor='k', alpha=0.15, zorder=-1000)

    def _plot_Z(self):
        for Z_container, scale_choice in zip(self.Z_containers, self.scale_choices):
            self.ax.errorbar(Z_container.mu_By_T, Z_container.Z_total, label=scale_choice.choice_label_for_plot)

    def _plot(self):
        self._setup_plot()
        self._plot_grey_flow_extr_band()
        self._plot_Z()
        self._finalize_plot()

    def _finalize_plot(self):
        self.ax.legend(fontsize=self.fontsize, title_fontsize=self.fontsize, framealpha=0, handlelength=3)
        file = self.outputpath_plot + "/Z_total" + self.output_suffix + ".pdf"
        self.fig.savefig(file)
        print("saved", file)
        matplotlib.pyplot.close(self.fig)

    @classmethod
    def plot(cls, *args):
        # This method only exists such that you don't have to manually instantiate the class.
        instance = cls(*args)
        instance._plot()


class ZTotalPlotterFlowtime(ZTotalPlotter):
    def _setup_plot(self):
        self.fig, self.ax, _ = lpd.create_figure(xlabel=r'$\sqrt{8 \tau_\text{F}}T$',
                                                 ylabel=r'$ Z$')
        global ylims
        self.ax.set_ylim(ylims)
        self.ax.set_xlim((0, 0.2))

    def _plot_grey_flow_extr_band(self):
        self.ax.fill_between([1/16, 1/6.66],
                             [-1, -1], [100, 100], facecolor='k', alpha=0.15, zorder=-1000)

    def _plot_Z(self):
        for Z_container, scale_choice in zip(self.Z_containers, self.scale_choices):
            self.ax.errorbar(np.flip(1/Z_container.mu_By_T), np.flip(Z_container.Z_total),
                             label=scale_choice.choice_label_for_plot)


class IntegrandPlotter:
    def __init__(self, args, data, fontsize, fmts):
        self.data: CouplingContainer = data
        self.outputpath_plot = args.outputpath_plot
        self.fontsize = fontsize
        self.fmts = fmts

    def __setup_plot(self):
        self.fig, self.ax, _ = lpd.create_figure(xlabel=r'$\mu_\mathrm{F}/T$',
                                                 ylabel=r'$ \frac{T}{\mu_\mathrm{F}} \displaystyle\gamma_0 g^2$')
        # self.ax.set_ylim((-0.06, 0.005))
        # self.ax.set_xlim((0, 21))

    def __plot_integrand(self):
        counter = 0
        self.ax.errorbar(self.data.mu_by_T, self.data.integrand_spline(self.data.mu_by_T), fmt=self.fmts[counter % 4],
                         zorder=counter)
        counter += 1
        # ax.axvline(x=reference_muF_by_T, **lpd.verticallinestyle)
        # ax.fill_between(flow_extr_window, [-250, -250], [400, 400], facecolor='grey', alpha=0.25, zorder=-1000)

    def __finalize_plot(self):
        file = self.outputpath_plot+"/integrand.pdf"
        self.fig.savefig(file)
        print("saved", file)
        matplotlib.pyplot.close(self.fig)

    def _plot(self):
        self.__setup_plot()
        self.__plot_integrand()
        self.__finalize_plot()

    @classmethod
    def plot(cls, *args):
        # This method only exists such that you don't have to manually instantiate the class.
        instance = cls(*args)
        instance._plot()


@lpd.typed_frozen_data
class CouplingContainer:
    mu_by_T: NDArray[Shape["*"], Float64]
    g2: NDArray[Shape["*"], Float64]
    g2_spline: interpolate._fitpack2.InterpolatedUnivariateSpline
    integrand_spline: interpolate._fitpack2.InterpolatedUnivariateSpline
    input_mu_by_T: NDArray[Shape["*"], Float64]


@lpd.typed_frozen_data
class ZContainer:
    mu_By_T: NDArray[Shape["*"], Float64]
    Z_K: NDArray[Shape["*"], Float64]
    Z_phys: float
    Z_run: NDArray[Shape["*"], Float64]
    Z_total: NDArray[Shape["*"], Float64]


@lpd.typed_frozen_data
class ScaleChoice:
    muRef_by_T: float
    muBarUV_by_muF: float
    muBarIR_by_T: float
    choice_label: str
    choice_label_for_plot: str
    choice_label_for_plot_short: str


def load_data(args: argparse.Namespace) -> CouplingContainer:
    muF_by_T, g2 = np.loadtxt(args.g2_file, unpack=True)[:2]
    g2_spline = interpolate.InterpolatedUnivariateSpline(muF_by_T, g2, k=3, ext=2)
    integrand_spline = interpolate.InterpolatedUnivariateSpline(muF_by_T, 2*gamma_0 * g2 / muF_by_T, k=3, ext=2)  # The factor 2 here is essential as we write the integral measure as dmu and not dmu^2
    input_mu_by_T = ScaleFunctions.get_muF_by_T_where_we_measure()
    return CouplingContainer(muF_by_T, g2, g2_spline, integrand_spline, input_mu_by_T)


class ZFactorComputer:
    def __init__(self, coupling_container, scale_choice):
        self.coupling_container = coupling_container
        self.scale_choice = scale_choice

    def compute_Zk(self):
        muF_by_T = self.coupling_container.input_mu_by_T
        muBarUV_by_T = self.scale_choice.muBarUV_by_muF * muF_by_T
        g2_MSBar = self.coupling_container.g2_spline(muBarUV_by_T)
        inner_bracket = np.log(muBarUV_by_T ** 2 / (4*self.coupling_container.input_mu_by_T ** 2)) + np.euler_gamma
        Zk = np.exp(- gamma_0 * g2_MSBar * inner_bracket)
        return Zk

    def compute_Z_run(self):

        def compute_Z_run(muF_by_T):
            muBarUV_By_T = self.scale_choice.muBarUV_by_muF * muF_by_T
            muBarIR_By_T = self.scale_choice.muBarIR_by_T
            # print("integrate", self.scale_choice.choice_label, lpd.format_float(muBarIR_By_T), lpd.format_float(muBarUV_By_T))
            this_Z_run = integrate.quad(self.coupling_container.integrand_spline, muBarIR_By_T, muBarUV_By_T, limit=300,
                                        epsabs=1e-5, epsrel=1e-5)[0]
            return this_Z_run

        integral = np.asarray([compute_Z_run(muF_by_T) for muF_by_T in self.coupling_container.input_mu_by_T])
        Z_run = np.exp(integral)

        return Z_run

    def compute_Z_phys(self):
        inner_bracket = np.log(self.scale_choice.muBarIR_by_T**2 / ((np.pi*4)**2)) - 2 + 2*np.euler_gamma
        Z_phys = np.exp(gamma_0 * self.coupling_container.g2_spline(self.scale_choice.muBarIR_by_T) * inner_bracket)
        return Z_phys

    @classmethod
    def compute(cls, *args):
        # This method only exists such that you don't have to manually instantiate the class.
        instance = cls(*args)
        Z_k = instance.compute_Zk()
        Z_phys = instance.compute_Z_phys()
        Z_run = instance.compute_Z_run()
        Z_total = Z_k * Z_phys * Z_run
        return ZContainer(instance.coupling_container.input_mu_by_T, Z_k, Z_phys, Z_run, Z_total)


def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('--g2_file', type=str, required=True, help="path to file containing g^2 in flow scheme")
    parser.add_argument('--outputpath_plot', default="/work/home/altenkort/work/correlators_flow/plots/quenched_1.50Tc_zeuthenFlow/coupling/")
    parser.add_argument('--outputpath_data', default="/work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/coupling/")

    args = parser.parse_args()

    return args


def get_scale_choices():
    scale_choices = []

    muBarIR_by_T_choices = [2 * np.pi, 4 * np.pi * np.exp(1 - np.euler_gamma)]
    muBarUV_by_muF_choices = [1., np.sqrt(4 * np.exp(-np.euler_gamma))]
    order_string = ["LO", "NLO"]
    muRef_by_T_choices = [4.,]

    for UV in range(len(muBarUV_by_muF_choices)):  # UV
        for IR in range(len(muBarIR_by_T_choices)):
            for ref in range(1):
                choice_label = "ref" + str(muRef_by_T_choices[ref])+"_UV"+order_string[UV] + "_IR" + order_string[IR]
                muBarUV_by_muF_choice = muBarUV_by_muF_choices[UV]
                muBarIR_by_T_choice = muBarIR_by_T_choices[IR]
                muRef_by_T = muRef_by_T_choices[ref]
                choice_label_for_plot = (r'$\bar{\mu}_\text{UV}/\mu_\text{F}='
                                         + lpd.format_float_latex(muBarUV_by_muF_choice, 2, 5)
                                         + r',\ \bar{\mu}_\text{IR}/T='
                                         + lpd.format_float_latex(muBarIR_by_T_choice, 2, 5)+r'$')
                # r'$\mu_\text{ref}/T='+lpd.format_float(muRef_by_T,1)
                choice_label_for_plot_short = lpd.format_float_latex(muBarUV_by_muF_choice, 2, 5) + r',\ ' + lpd.format_float_latex(muBarIR_by_T_choice, 2, 5)
                scale_choices.append(ScaleChoice(muRef_by_T, muBarUV_by_muF_choice, muBarIR_by_T_choice, choice_label, choice_label_for_plot, choice_label_for_plot_short))
    return scale_choices


def save_Z_to_file(args, Z_container, scale_choice):
    outputfolder = args.outputpath_data
    file = "Z_match_"+scale_choice.choice_label+".dat"
    filename = outputfolder+file
    sqrt8taufT = np.flip(1/Z_container.mu_By_T)
    Z = np.flip(Z_container.Z_total)
    lpd.save_columns_to_file(filename, (sqrt8taufT, Z), ("sqrt(8tauF)T", "Z"))


def main():
    args = parse_args()
    coupling_container = load_data(args)
    scale_choices = get_scale_choices()

    fontsize = 8
    fmts = ['-', '-', '--', ':']

    IntegrandPlotter.plot(args, coupling_container, fontsize, fmts)

    Z_containers = [ZFactorComputer.compute(coupling_container, scale_choice) for scale_choice in scale_choices]
    for Z_container, scale_choice in zip(Z_containers, scale_choices):
        ZIndividualPlotter.plot(args, Z_container, fontsize, scale_choice)
        save_Z_to_file(args, Z_container, scale_choice)

    ZTotalPlotter.plot(args, Z_containers, fontsize, scale_choices)
    ZTotalPlotterFlowtime.plot(args, Z_containers, fontsize, scale_choices, "_flowtime")


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
