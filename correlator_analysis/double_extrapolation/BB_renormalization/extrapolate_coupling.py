#!/usr/bin/env python3
import lib_process_data as lpd
import numpy as np
from latqcdtools.statistics import bootstr
import argparse
import scipy.optimize
import scipy.interpolate
import latqcdtools.physics.referenceScales as rS
import rundec
import os

# Type hints
from nptyping import NDArray, Float64
from typing import Literal as Shape
from beartype.typing import List


class ScaleFunctions:
    @classmethod
    def convert_sqrt8taufTByTau_to_muFByT(cls, sqrt8taufTByTau, tauT):
        muFByT = 1 / (sqrt8taufTByTau * tauT)
        return muFByT

    @classmethod
    def convert_taufT2_to_muFByT(cls, taufT2):
        muFByT = 1 / np.sqrt(8 * taufT2)
        return muFByT

    @classmethod
    def get_muF_by_T_where_we_measure(cls):
        muF_by_T_where_we_measure = np.sort(np.unique(np.concatenate([cls.convert_sqrt8taufTByTau_to_muFByT(lpd.get_relflow_range(), tauT)
                                                                      for tauT in np.arange(0.25, 0.501, 1 / 36)])))
        return muF_by_T_where_we_measure

    @classmethod
    def get_min_and_max_indices(cls, muF_by_T, min_muF_by_T, max_muF_by_T):
        min_idx = np.where(muF_by_T <= min_muF_by_T)[0][-1]
        max_idx = np.where(muF_by_T >= max_muF_by_T)[0][0] + 1
        return min_idx, max_idx


@lpd.typed_frozen_data
class GeneralParameters:
    Ntexp: int  # We extrapolate linear in a**Ntexp
    factor: int  # Scale factor to make cont. extr. converge
    T_in_GeV: float
    muUV_by_T_LO: float  # The "physical" scale we run to
    muUV_by_T_NLO: float
    min_muF_by_T_in_flow_extr: float  # The minimum flow scale inside the flow extrapolation at which we measure the correlator on the lattice.
    max_muF_by_T_in_flow_extr: float  # The maximum "
    max_muBar_by_T: float  # Largest scale we ever want to run to.
    muF_by_T_where_we_measure: NDArray[Shape["*"], Float64]


@lpd.typed_frozen_data
class RawDataContainer:
    muF_by_T_arr: NDArray[Shape["*, *"], Float64]
    g2_arr: NDArray[Shape["*, *"], Float64]
    g2_err_arr: NDArray[Shape["*, *"], Float64]
    g2_ints: List[scipy.interpolate._fitpack2.InterpolatedUnivariateSpline]
    g2_err_ints: List[scipy.interpolate._fitpack2.InterpolatedUnivariateSpline]
    target_muF_by_T_for_cont_extr: NDArray[Shape["*"], Float64]


@lpd.typed_frozen_data
class ContDataContainer:
    muF_by_T_cont: NDArray[Shape["*"], Float64]
    g2_latt_from_int: NDArray[Shape["* nflow, * nNts"], Float64]
    g2_latt_from_int_err: NDArray[Shape["* nflow, * nNts"], Float64]
    g2_cont: NDArray[Shape["* nflow, [g2, g2err]"], Float64]
    g2_extr_slope: NDArray[Shape["* nflow, [slope, slope_err]"], Float64]
    g2_extr_chisqdof: NDArray[Shape["* nflow, [chisqdof, chisqdof_err]"], Float64]
    g2_cont_int: scipy.interpolate._fitpack2.InterpolatedUnivariateSpline


@lpd.typed_frozen_data
class MSBarDataContainer:
    mubar_by_T: NDArray[Shape["*"], Float64]
    g2_MSBAR_spline: scipy.interpolate._fitpack2.InterpolatedUnivariateSpline


@lpd.typed_frozen_data
class PertRunMSBarDataContainer:
    mubar_by_T: NDArray[Shape["*"], Float64]
    g2_MSBAR_pert_run: scipy.interpolate._fitpack2.InterpolatedUnivariateSpline
    reference_muF_by_T: float
    nloop: int


def plot_g2_vs_1overNt2(args, general_params, cont_data_container):
    # plot continuum extrapolation at different flow times
    xlabel = r'$N_\tau^{-' + str(general_params.Ntexp) + r'}$'
    fig3, ax3, _ = lpd.create_figure(xlabel=xlabel, ylabel=r'$ \displaystyle \alpha_\text{flow}$', ylims=(0.95/(np.pi*4), 3.55/(np.pi*4)))

    xpoints = np.linspace(0, 1 / args.Nts[-1] ** general_params.Ntexp, 10)

    plot_muF_by_T = np.geomspace(4, general_params.max_muF_by_T_in_flow_extr, 10)

    for j, muF_by_T in enumerate(plot_muF_by_T):
        i = (np.fabs(cont_data_container.muF_by_T_cont - muF_by_T)).argmin()
        color = lpd.get_color(plot_muF_by_T, j)
        ax3.errorbar(np.insert(1 / np.asarray(args.Nts) ** general_params.Ntexp, 0, 0),
                     np.insert(cont_data_container.g2_latt_from_int[i]/(np.pi*4), 0, cont_data_container.g2_cont[i, 0]/(np.pi*4)),
                     np.insert(cont_data_container.g2_latt_from_int_err[i]/(np.pi*4), 0, cont_data_container.g2_cont[i, 1]/(np.pi*4)),
                     color=color,
                     fmt='|', label=lpd.format_float(cont_data_container.muF_by_T_cont[i], 2), zorder=-i)
        ax3.errorbar(xpoints, extrapolation_ansatz(xpoints, cont_data_container.g2_extr_slope[i, 0] * general_params.factor/(np.pi*4),
                                                   cont_data_container.g2_cont[i, 0]/(np.pi*4)),
                     **lpd.fitlinestyle, color=color)
        # counter += 1
    ax3.legend(title=r'$\displaystyle \mu_\mathrm{F}/ T$', loc="center left", bbox_to_anchor=(1, 0.5),
               **lpd.leg_err_size())

    thisxlims = ax3.get_xlim()
    ax3.set_xlim((thisxlims[0], thisxlims[1] * 1.025))

    ax3.figure.canvas.draw()
    offset = ax3.xaxis.get_major_formatter().get_offset()
    ax3.xaxis.offsetText.set_visible(False)
    ax3.xaxis.set_label_text(xlabel + " " + offset)

    file = args.outputpath_plot + "g2" + "_cont_extr.pdf"
    lpd.create_folder(os.path.dirname(file))
    fig3.savefig(file)
    print("saved", file)


class Plotter_g2_vs_mu:
    def __init__(self, args, general_params, raw_data_container, cont_data_container, msbar_data_container,
                 list_of_pertRun_MSBar_data_container, file_suffix=""):
        self.args = args
        self.general_params: GeneralParameters = general_params
        self.raw_data_container: RawDataContainer = raw_data_container
        self.cont_data_container: ContDataContainer = cont_data_container
        self.msbar_data_container: MSBarDataContainer = msbar_data_container
        self.list_of_pertRun_MSBar_data_container: List[PertRunMSBarDataContainer] = list_of_pertRun_MSBar_data_container
        self.linewidth: float = 1
        self.file_suffix: str = file_suffix

    def __setup_plot(self):
        self.fig, self.ax, _ = lpd.create_figure(xlabel=r'$\displaystyle \mu_{\mathrm{F}} / T$', ylabel=None, ylims=(0, 4.5/(np.pi*4)), xlims=(1, 200))
        self.ax.set_xscale('log')

    def __plot_pert_g2(self, mus, suffix="pert"):
        g2s = np.asarray([lpd.get_g2_pert(mu * self.general_params.T_in_GeV, Nf=0, Nloop=3) for mu in mus])
        self.ax.errorbar(mus, g2s/(np.pi*4), fmt='-.', lw=self.linewidth, markersize=0,
                         label=r'pert. ($\mu = \mu_\mathrm{F} \equiv 1/\sqrt{8\tau_F} $)', zorder=-1)
        lpd.save_columns_to_file(self.args.outputpath_data + "g2_muF_" + suffix + ".txt", (mus, g2s), ["mu_F/T", "g2s"])

    def __plot_lattice_data(self):
        for i, Nt in enumerate(self.args.Nts):
            self.ax.errorbar(self.raw_data_container.muF_by_T_arr[i], self.raw_data_container.g2_arr[i]/(np.pi*4), fmt=':',
                             lw=self.linewidth*0.75, markersize=0, label=r'$\alpha_\text{flow}$ ($N_\tau = ' + str(Nt) + r'$)', zorder=-i - 10)

    def __plot_grey_flow_extr_band(self):
        self.ax.axvline(self.general_params.muUV_by_T_NLO, alpha=0.5, dashes=(2, 2), zorder=-10000, lw=self.linewidth/2, color='k')
        self.ax.axvline(self.general_params.muUV_by_T_LO, alpha=0.5, dashes=(2, 2), zorder=-10000, lw=self.linewidth/2, color='k')
        self.ax.fill_between([self.general_params.min_muF_by_T_in_flow_extr, self.general_params.max_muF_by_T_in_flow_extr],
                             [-1, -1], [100, 100], facecolor='k', alpha=0.15, zorder=-1000)

    def __plot_cont(self):
        cont_data = self.cont_data_container.g2_cont[:, 0]
        self.ax.errorbar(self.cont_data_container.muF_by_T_cont, cont_data/(np.pi*4), fmt='-', markersize=0, lw=self.linewidth,
                         label=r'$\alpha_\text{flow}$ (cont.\ extr.)', zorder=-1)

    def __plot_pert_run_from_converted_g2MSbar(self):
        for data in self.list_of_pertRun_MSBar_data_container:
            self.ax.errorbar(data.mubar_by_T, data.g2_MSBAR_pert_run(data.mubar_by_T)/(np.pi*4), fmt='-', markersize=0,
                             lw=self.linewidth, zorder=(int(data.mubar_by_T[0])) + 20,
                             label=r'\begin{flushleft}$\alpha_s$ via Eq.\,'+self.args.eq_no+r' at \newline  ' + r' ${\mu_\text{F}}=\bar{\mu}_\mathrm{ref}'
                                   r'$, then \newline run via ' + str(
                             data.nloop) + r'-loop $\beta$-fct. \end{flushleft}')
            self.ax.axvline(x=data.reference_muF_by_T, alpha=1, dashes=(1, 1), zorder=-10000, lw=self.linewidth/2, color='k')

    def __plot_MSBAR(self):
        self.ax.errorbar(self.msbar_data_container.mubar_by_T,
                    self.msbar_data_container.g2_MSBAR_spline(self.msbar_data_container.mubar_by_T)/(np.pi*4),
                    fmt='-', markersize=0, lw=self.linewidth, label=r'$\alpha_s$ via Eq.\,'+self.args.eq_no)

    def __finalize_plot(self):
        self.ax.legend(**lpd.leg_err_size(), loc="upper right", bbox_to_anchor=(1.01, 1.01), handlelength=1, title_fontsize=9, fontsize=10, framealpha=0)
        file = self.args.outputpath_plot + "g2" + self.file_suffix + ".pdf"
        print("saving", file)
        lpd.create_folder(os.path.dirname(file))
        self.fig.savefig(file)

    def _plot(self):
        self.__setup_plot()
        if False:
            self.__plot_grey_flow_extr_band()
        self.__plot_pert_run_from_converted_g2MSbar()
        self.__plot_MSBAR()
        self.__plot_cont()
        if False:
            self.__plot_lattice_data()

        # self.__plot_pert_g2(ax, np.geomspace(2, muB_by_T, 200))
        self.__finalize_plot()

    @classmethod
    def plot(cls, *args, **kwargs):
        # This method only exists such that you don't have to manually instantiate the class.
        instance = cls(*args, **kwargs)
        return instance._plot()


def extrapolation_ansatz(x, m, b):
    return m * x + b


class ContinuumExtrapolator:
    def __init__(self, args, general_parameters, raw_data_container):
        self.calc_cont = args.calc_cont
        self.output_file = args.outputpath_data + "g2_muF" + "_cont_extr.txt"

        # raw data and corresponding parameters
        self.general_parameters = general_parameters
        self.raw_data_container = raw_data_container
        self.nflow = len(raw_data_container.target_muF_by_T_for_cont_extr)
        self.Nts = args.Nts
        self.nNts = len(args.Nts)

        # some empty numpy arrays
        self.g2_latt_from_int = np.empty((self.nflow, self.nNts))
        self.g2_latt_from_int_err = np.empty((self.nflow, self.nNts))
        self.cont = np.empty((self.nflow, 2))
        self.slope = np.empty((self.nflow, 2))
        self.chisqdof = np.empty((self.nflow, 2))

        # start reading here
        self.__get_cont_data()

    def __get_interpolated_lattice_data(self):
        for i, muF_by_T in enumerate(self.raw_data_container.target_muF_by_T_for_cont_extr):
            self.g2_latt_from_int[i] = [spline(muF_by_T) for spline in self.raw_data_container.g2_ints]
            self.g2_latt_from_int_err[i] = [spline(muF_by_T) for spline in self.raw_data_container.g2_err_ints]
        print("max relative error in % in the lattice data:",
              np.amax(self.g2_latt_from_int_err / self.g2_latt_from_int) * 100)

    def __save_cont_extr_to_disk(self):
        columns = (self.raw_data_container.target_muF_by_T_for_cont_extr, self.cont, self.slope, self.chisqdof)
        labels = ["mu_F/T", "g2", "err", "slope", "err", "chisqdof", "err"]

        lpd.save_columns_to_file(self.output_file, columns, labels)

    def __load_cont_extr_from_disk(self):
        try:
            _, b, berr, m, merr, chisqdof, chisqdof_err = np.loadtxt(self.output_file, unpack=True)
            print("load", self.output_file)
        except OSError:
            print("Error: could not find ", self.output_file)
            exit(1)
        for i in range(self.nflow):
            self.cont[i, 0] = b[i]
            self.cont[i, 1] = berr[i]
            self.slope[i, 0] = m[i]
            self.slope[i, 1] = merr[i]
            # here we could also load the chisqdof, but we don't do anything with it later, so we don't load it.

    @classmethod
    def _fit_sample(cls, ydata, xdata, edata):
        def chisqdof_func(fitparams_, xdata_, ydata_, edata_):
            return (np.sum(((ydata_ - extrapolation_ansatz(xdata_, *fitparams_)) / edata_) ** 2)
                    / (len(ydata_) - len(fitparams_)))
        fitparams = scipy.optimize.minimize(chisqdof_func, x0=np.asarray([-1, 2]), args=(xdata, ydata, edata))
        fitparams = fitparams.x
        this_chisqdof = chisqdof_func(fitparams, xdata, ydata, edata)
        return [*fitparams, this_chisqdof]

    def __do_cont_extr(self):
        for i in range(self.nflow):
            fitparams, fitparams_err = bootstr.bootstr_from_gauss(
                ContinuumExtrapolator._fit_sample, data=self.g2_latt_from_int[i],
                data_std_dev=self.g2_latt_from_int_err[i],
                numb_samples=100, sample_size=1, return_sample=False,
                args=[1 / np.asarray(self.Nts) ** self.general_parameters.Ntexp * self.general_parameters.factor, self.g2_latt_from_int_err[i]], parallelize=True,
                nproc=20)  # TODO make nproc an input parameter
            self.chisqdof[i, 0] = fitparams[2]
            self.chisqdof[i, 1] = fitparams_err[2]
            self.cont[i, 0] = fitparams[1]
            self.cont[i, 1] = fitparams_err[1]
            self.slope[i, 0] = fitparams[0]
            self.slope[i, 1] = fitparams_err[0]

    def __get_cont_data(self):
        self.__get_interpolated_lattice_data()
        if self.calc_cont:
            self.__do_cont_extr()
            self.__save_cont_extr_to_disk()
        else:
            self.__load_cont_extr_from_disk()

        # don't save/load this to/from disk. just recompute it every time.
        g2_cont_spline = scipy.interpolate.InterpolatedUnivariateSpline(
            self.raw_data_container.target_muF_by_T_for_cont_extr,
            self.cont[:, 0], k=3, ext=2)

        self.cont_data = ContDataContainer(self.raw_data_container.target_muF_by_T_for_cont_extr,
                                           self.g2_latt_from_int,
                                           self.g2_latt_from_int_err,
                                           self.cont,
                                           self.slope,
                                           self.chisqdof,
                                           g2_cont_spline)

    @classmethod
    def extrapolate(cls, *args):
        # This method only exists such that you don't have to manually instantiate the class.
        instance = cls(*args)
        return instance.cont_data


def flow_coupling(tauF2E: float, a2bytauF: float) -> float:
    # some constants
    prefactor = np.pi ** 2 / 8
    c_0 = 3 / 128
    # c_2 = -1/256
    # c_4 = 343/327680
    c_2 = 0
    c_4 = 0

    this_flow_coupling = prefactor * tauF2E / (c_0 + c_2 * a2bytauF + c_4 * a2bytauF ** 2)

    return this_flow_coupling


def load_data_and_interpolate(args: argparse.Namespace, general_params: GeneralParameters) -> RawDataContainer:
    def pseudo_nt_for_1point5_Tc(beta):
        r0Tc = 0.7457
        TbyTc = 1.5
        r0_T = r0Tc * TbyTc

        # 1/(a T). a is lattice spacing. this is a "pseudo" Nt
        # (= the Nt, that we should have chosen to get T=1.5Tc while keeping beta fixed.)
        one_over_aT = rS.r0_div_a(beta, 2017) / r0_T
        return one_over_aT

    # declare a few arrays
    muF_by_T_arr = []
    g2_arr = []
    g2_err_arr = []
    g2_ints = []
    g2_err_ints = []

    order = 3

    largest_min_muF_by_T = 0

    # load data and save various "versions" of it
    for i in range(len(args.input_files)):
        tauFbya2, tauF2E, tauF2E_err = np.loadtxt(args.input_basepath + "/" + args.input_files[i], unpack=True)

        # The continuum extrapolation happens at fixed physical scale, so we need to "align" the lattice measurements as they
        # have different lattice spacings. The most intuitive scale is the "target" temperature T=1.500Tc.
        # But, due to slight mistuning of the beta values, not all lattices are at this exact same temperature.
        # However, using r0Tc and the r0-scale, we can compute a non-integer "pseudo" Nt that we can use to convert
        # from the individual lattice spacing scales to a common 1.500Tc=1/(a Nt_pseudo) scale.
        tauFbyScale = tauFbya2 / pseudo_nt_for_1point5_Tc(args.betas[i]) ** 2
        muF_by_T = np.flip(ScaleFunctions.convert_taufT2_to_muFByT(tauFbyScale))
        muF_by_T_arr.append(muF_by_T)

        largest_min_muF_by_T = max(largest_min_muF_by_T, np.min(muF_by_T))

        # coupling
        thisg2 = np.flip(flow_coupling(tauF2E, 1 / tauFbya2))
        thisg2_err = np.flip(np.fabs(flow_coupling(tauF2E_err, 1 / tauFbya2)))
        g2_arr.append(thisg2)
        g2_err_arr.append(thisg2_err)

    for i in range(len(args.input_files)):
        muF_by_T = muF_by_T_arr[i]
        thisg2 = g2_arr[i]
        thisg2_err = g2_err_arr[i]
        g2_ints.append(scipy.interpolate.InterpolatedUnivariateSpline(muF_by_T, thisg2, k=order, ext=2))
        g2_err_ints.append(scipy.interpolate.InterpolatedUnivariateSpline(muF_by_T, thisg2_err, k=order, ext=2))

    target_muF_by_T_for_cont_extr = np.sort(np.unique(np.concatenate((
        np.geomspace(largest_min_muF_by_T, general_params.muUV_by_T_NLO, 500),
        general_params.muF_by_T_where_we_measure,
        [general_params.muUV_by_T_LO, general_params.muUV_by_T_NLO]))))

    return RawDataContainer(np.asarray(muF_by_T_arr), np.asarray(g2_arr), np.asarray(g2_err_arr), g2_ints, g2_err_ints,
                            target_muF_by_T_for_cont_extr)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('--calc_cont',
                        help='calc continuum extrapolation and save to file instead of reading it from the file',
                        action="store_true")
    parser.add_argument('--input_basepath',
                        default="/work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/coupling/")
    parser.add_argument('--input_files', nargs='*', type=str, required=True)
    parser.add_argument('--outputpath_plot',
                        default="/work/home/altenkort/work/correlators_flow/plots/quenched_1.50Tc_zeuthenFlow/coupling/")
    parser.add_argument('--outputpath_data',
                        default="/work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/coupling/")
    parser.add_argument('--Nts', nargs='*', type=int)
    parser.add_argument('--betas', nargs='*', type=float)
    parser.add_argument('--eq_no', type=str, default="20", help="Equation number for matching in alphas vs mu plot.")
    args = parser.parse_args()

    if len(args.betas) != len(args.input_files) or len(args.Nts) != len(args.input_files):
        print("ERROR: --betas and --input_files must have same number of arguments.")
        exit(1)

    return args


def convert_g2_flow_to_MSBAR(g2_flow: float, nf: int = 0) -> float:
    # Harlander, R.V., Neumann, T. The perturbative QCD gradient flow to three loops. J. High Energ. Phys. 2016, 161 (2016).
    # https://doi.org/10.1007/JHEP06(2016)161

    k2 = -0.982 - 0.070 * nf + 0.002 * nf ** 2
    k1 = 1.098 + 0.008 * nf

    alpha_flow = g2_flow / (4 * np.pi)

    coefficients = [k2, k1, 1, -alpha_flow]

    solutions = np.roots(coefficients)

    # Filter the real solutions
    real_solution = [sol * 4 * np.pi for sol in solutions if np.isreal(sol) and sol >= 0 and sol * 4 * np.pi <= 5][0]

    return real_solution


def calc_g2_MSBAR(cont_data_container: ContDataContainer, target_muF_by_T_for_cont_extr) -> MSBarDataContainer:
    g2_MSBAR = np.asarray([convert_g2_flow_to_MSBAR(cont_data_container.g2_cont_int(mu_by_T)) for mu_by_T in target_muF_by_T_for_cont_extr])
    g2_MSBAR_spline = scipy.interpolate.InterpolatedUnivariateSpline(target_muF_by_T_for_cont_extr, g2_MSBAR, k=3, ext=2)
    return MSBarDataContainer(cont_data_container.muF_by_T_cont, g2_MSBAR_spline)


class BetaFunction:
    def __init__(self, args, MSBar_data_container, general_params):
        self.MSBar_data_container = MSBar_data_container
        self.crd = rundec.CRunDec()
        self.nf: int = 0
        self.nloop: int = 5
        self.outputpath_data = args.outputpath_data
        self.list_of_data_containers = []
        self.min_muF_by_T_where_we_measure = np.min(general_params.muF_by_T_where_we_measure)

        self.calc_pert_run_of_g2_MSBAR()

    def calc_mu_and_g2(self, mu, g2_0, mu_0):
        alpha_s0 = g2_0 / (4 * np.pi)
        Alphas = self.crd.AlphasExact(alpha_s0, mu_0, mu, self.nf, self.nloop)
        g2 = 4. * np.pi * Alphas
        return g2

    def calc_pert_run_of_g2_MSBAR(self):
        muF_by_T_targets_extra = np.geomspace(np.max(self.MSBar_data_container.mubar_by_T), 220, 200)
        target_muF_by_T = (
            np.sort(
                np.unique(
                    np.concatenate((self.MSBar_data_container.mubar_by_T, muF_by_T_targets_extra)))))

        ref_mu_by_T_choices = [2*np.pi, ]  # TODO make this an input, read from somewhere else not hard code
        for reference_muF_by_T in ref_mu_by_T_choices:
            mu_0 = reference_muF_by_T
            g2_0 = self.MSBar_data_container.g2_MSBAR_spline(mu_0)

            these_target_mus = np.asarray([mu for mu in target_muF_by_T if mu >= self.min_muF_by_T_where_we_measure])

            g2s = np.asarray([self.calc_mu_and_g2(mu, g2_0, mu_0) for mu in these_target_mus])
            g2_spline = scipy.interpolate.InterpolatedUnivariateSpline(these_target_mus, g2s, k=3, ext=2)
            self.list_of_data_containers.append(
                PertRunMSBarDataContainer(these_target_mus, g2_spline, reference_muF_by_T, self.nloop))

            lpd.save_columns_to_file(
                self.outputpath_data + "g2_MSBAR_runFromMu_" + lpd.format_float(reference_muF_by_T, 2) + ".txt",
                (these_target_mus, g2s), ["mu/T", "g^2"])

    @classmethod
    def run(cls, *args):
        # This method only exists such that you don't have to manually instantiate the class.
        instance = cls(*args)
        return instance.list_of_data_containers


def get_general_parameters():
    NTEXP = 2
    max_muF_by_T_in_flow_extr = ScaleFunctions.convert_sqrt8taufTByTau_to_muFByT(sqrt8taufTByTau=0.25, tauT=0.25)
    muUV_by_T_NLO = 4*np.pi*np.exp(1-np.euler_gamma)       # TODO rename these according to the paper
    muBar_by_muF_NLO = np.sqrt(4*np.exp(-np.euler_gamma))
    general_params = GeneralParameters(
        Ntexp=NTEXP,
        factor=96 ** NTEXP,
        T_in_GeV=0.472,
        muUV_by_T_LO=2*np.pi,
        muUV_by_T_NLO=muUV_by_T_NLO,
        min_muF_by_T_in_flow_extr=ScaleFunctions.convert_sqrt8taufTByTau_to_muFByT(sqrt8taufTByTau=0.3, tauT=0.5),
        max_muF_by_T_in_flow_extr=max_muF_by_T_in_flow_extr,
        max_muBar_by_T=np.fmax(muBar_by_muF_NLO*max_muF_by_T_in_flow_extr, muUV_by_T_NLO),
        muF_by_T_where_we_measure=ScaleFunctions.get_muF_by_T_where_we_measure())

    return general_params


def main():
    # this script load tf^2 <E> from some files, then calculates the tree-level improved coupling g^2 and does a continuum extrapolation of g^2

    args = parse_args()

    general_params = get_general_parameters()
    raw_data = load_data_and_interpolate(args, general_params)
    cont_data = ContinuumExtrapolator.extrapolate(args, general_params, raw_data)
    MSBar_data = calc_g2_MSBAR(cont_data, raw_data.target_muF_by_T_for_cont_extr)
    list_of_pertRun_MSBar_data = BetaFunction.run(args, MSBar_data, general_params)

    Plotter_g2_vs_mu.plot(args, general_params, raw_data, cont_data, MSBar_data,
                          list_of_pertRun_MSBar_data)
    plot_g2_vs_1overNt2(args, general_params, cont_data)


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
