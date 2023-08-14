#!/usr/bin/env python3
import lib_process_data as lpd
import numpy
from latqcdtools.statistics import bootstr
import argparse
import scipy.optimize
import scipy.interpolate
import latqcdtools.physics.referenceScales as rS

# Type hints
from nptyping import NDArray, Float64
from typing import Literal as Shape
from beartype.typing import List


Ntexp = 2
factor = 96**Ntexp


def convert_sqrt8taufTByTau_to_muFByT(sqrt8taufTByTau, tauT):
    muFByT = 1/(sqrt8taufTByTau*tauT)
    return muFByT


def convert_taufT2_to_muFByT(taufT2):
    muFByT = 1/numpy.sqrt(8*taufT2)
    return muFByT


def get_necessary_muF_by_T():
    necessary_muF_by_T = numpy.asarray([])

    for tauT in numpy.arange(0.25, 0.51, 1 / 36):
        tmp = convert_sqrt8taufTByTau_to_muFByT(lpd.get_relflow_range(), tauT)
        necessary_muF_by_T = numpy.concatenate((necessary_muF_by_T, tmp))

    return necessary_muF_by_T


muB_by_T = 19.179
min_muF_by_T_in_flow_extr = convert_sqrt8taufTByTau_to_muFByT(sqrt8taufTByTau=0.3, tauT=0.5)
max_muF_by_T_in_flow_extr = convert_sqrt8taufTByTau_to_muFByT(sqrt8taufTByTau=0.25, tauT=0.25)

print(min_muF_by_T_in_flow_extr, max_muF_by_T_in_flow_extr)

necessary_muF_by_T = get_necessary_muF_by_T()
max_muF_by_T_for_running = numpy.max([*necessary_muF_by_T, muB_by_T])
min_muF_by_T_for_running = numpy.min([*necessary_muF_by_T, muB_by_T])


T_in_GeV = 0.472


def get_min_and_max_indices(muF_by_T, min_muF_by_T, max_muF_by_T):
    # print(numpy.where(muF_by_T > min_muF_by_T))
    min_idx = numpy.where(muF_by_T <= min_muF_by_T)[0][-1]
    max_idx = numpy.where(muF_by_T >= max_muF_by_T)[0][0]+1
    return min_idx, max_idx


def get_ylabel():
    # return r'$ g^{2 \text{(flow)}} = \frac{\pi^2}{8} \tau_F^2 \langle E \rangle / \frac{3}{128}$'
    # return r'$ \scriptstyle g^{2 \text{(flow)}} = \frac{\pi^2\tau_F^2 \langle E\rangle}{8}  / \left( \frac{3}{128} + \scriptstyle \underset{n=1,2}{\sum} \scriptscriptstyle  c_{2n} \frac{a^{2n}}{\tau_\mathrm{F}^{n}} \right)$'
    return r'$ \displaystyle g^{2}$'


def plot_g2_vs_1overNt2(args, cont_data_container):

    # plot continuum extrapolation at different flow times
    xlabel = r'$N_\tau^{-'+str(Ntexp)+r'}$'
    fig3, ax3, _ = lpd.create_figure(xlabel=xlabel, ylabel=get_ylabel())

    xpoints = numpy.linspace(0, 1 / args.Nts[-1] ** Ntexp, 10)

    plot_muF_by_T = numpy.geomspace(min_muF_by_T_in_flow_extr, max_muF_by_T_for_running, 10)

    for j, muF_by_T in enumerate(plot_muF_by_T):
        i = (numpy.fabs(cont_data_container.muF_by_T_cont - muF_by_T)).argmin()
        color = lpd.get_color(plot_muF_by_T, j)
        ax3.errorbar(numpy.insert(1 / numpy.asarray(args.Nts)**Ntexp, 0, 0),
                     numpy.insert(cont_data_container.g2_latt_from_int[i], 0, cont_data_container.g2_cont[i, 0]),
                     numpy.insert(cont_data_container.g2_latt_from_int_err[i], 0, cont_data_container.g2_cont[i, 1]),
                     color=color,
                     fmt='|', label=lpd.format_float(cont_data_container.muF_by_T_cont[i], 2), zorder=-i)
        ax3.errorbar(xpoints, extrapolation_ansatz(xpoints, cont_data_container.g2_extr_slope[i, 0] * factor, cont_data_container.g2_cont[i, 0]),
                     **lpd.fitlinestyle, color=color)
        # counter += 1
    ax3.legend(title=r'$\displaystyle \mu_\mathrm{F}/ T$', loc="center left", bbox_to_anchor=(1, 0.5), **lpd.leg_err_size())

    # thisylims = ax3.get_ylim()
    # ax3.set_ylim((0, thisylims[1]*1.1))
    ax3.set_ylim(0.7, 2.8)
    thisxlims = ax3.get_xlim()
    ax3.set_xlim((thisxlims[0], thisxlims[1] * 1.025))

    ax3.figure.canvas.draw()
    offset = ax3.xaxis.get_major_formatter().get_offset()
    print(offset)
    ax3.xaxis.offsetText.set_visible(False)
    ax3.xaxis.set_label_text(xlabel + " " + offset)

    file = args.outputpath_plot + "g2" + "_cont_extr.pdf"
    print("saving", file)
    fig3.savefig(file)


class Plotter_g2_vs_mu:
    def __init__(self, args, raw_data_container, cont_data_container, msbar_data_container):
        self.args = args
        self.raw_data_container: RawDataContainer = raw_data_container
        self.cont_data_container: ContDataContainer = cont_data_container
        self.msbar_data_container: MSBarDataContainer = msbar_data_container
        self.linewidth = 1

    def __setup_plot(self):
        fig, ax, _ = lpd.create_figure(xlabel=r'$\displaystyle \mu_{\mathrm{F}} / T$', ylabel=get_ylabel())
        ax.set_ylim(0, 4.5)
        ax.set_xlim(1, 100)
        ax.set_xscale('log')
        return fig, ax

    def __finalize_plot(self, fig, ax, file_suffix):
        ax.axvline(muB_by_T, alpha=1, dashes=(2, 2), zorder=-10000, lw=self.linewidth, color='k')

        ax.legend(**lpd.leg_err_size(), loc="lower left", bbox_to_anchor=(0, 0), handlelength=1, fontsize=8, framealpha=0)

        file = self.args.outputpath_plot + "g2" + file_suffix + ".pdf"
        print("saving", file)
        fig.savefig(file)

    def __plot_pert_g2(self, ax, mus, suffix="pert"):
        g2s = []
        for mu in mus:
            g2s.append(lpd.get_g2_pert(mu*T_in_GeV, Nf=0, Nloop=3))
        ax.errorbar(mus, g2s, fmt='-.', lw=self.linewidth,  markersize=0, label=r'pert. ($\mu = \mu_\mathrm{F} \equiv 1/\sqrt{8\tau_F} $)', zorder=-1)
        lpd.save_columns_to_file(self.args.outputpath_data + "g2_muF_" + suffix + ".txt", (mus, g2s), ["mu_F/T", "g2s"])

    def __plot_lattice_data(self, ax):
        for i, Nt in enumerate(self.args.Nts):
            ax.errorbar(self.raw_data_container.muF_by_T_arr[i], self.raw_data_container.g2_arr[i], fmt='-', lw=self.linewidth, markersize=0, label=r'$N_\tau = ' + str(Nt) + r'$', zorder=-i - 10)

    def __plot_grey_flow_extr_band(self, ax):
        ax.fill_between([min_muF_by_T_in_flow_extr, max_muF_by_T_in_flow_extr],
                         [-1, -1], [100, 100], facecolor='k', alpha=0.15, zorder=-1000)

    def __plot_cont(self, ax):
        cont_data = self.cont_data_container.g2_cont[:, 0]
        ax.errorbar(self.cont_data_container.muF_by_T_cont, cont_data, fmt='-', markersize=0, lw=self.linewidth, label=r'flow', zorder=-1)  # r'cont., linear in $a^' + str(Ntexp) + r'$'
        # tmp_idx = numpy.fabs(self.cont_data_container.muF_by_T_cont-19.179).argmin()
        # print(self.cont_data_container.muF_by_T_cont[tmp_idx], cont_data[tmp_idx])

    def __plot_pert_run_starting_from_reference(self, ax, g2_0, mu_0, target_mus, color, suffix):
        alpha_s0 = g2_0 / (4 * numpy.pi)
        import rundec
        crd = rundec.CRunDec()
        g2s = []
        nf = 0
        nloop = 5
        plot_mus = []
        for mu in target_mus:
            if mu >= mu_0:
                Alphas = crd.AlphasExact(alpha_s0, mu_0, mu, nf, nloop)
                g2 = 4. * numpy.pi * Alphas
                g2s.append(g2)
                plot_mus.append(mu)
        plot_mus = numpy.asarray(plot_mus)
        ax.errorbar(plot_mus, g2s, fmt='--', markersize=0, lw=self.linewidth, zorder=(int(mu_0))+20, color=color,
                    label=r'\begin{flushleft}$\text{flow}\rightarrow\overline{\text{MS}}$ at ' + r' $\frac{\mu_\text{F}}{T}=' + lpd.format_float(mu_0, 1) + r'$ \newline with ' + str(nloop)+r'-loop run \end{flushleft}')
        lpd.save_columns_to_file(self.args.outputpath_data + "g2_muF_" + suffix + ".txt", (plot_mus, g2s), ["mu_F/T", "g2s"])

    def __nonpert_ref_with_pert_run_wrapper(self, ax, cont, target_muF_by_T, color, ref_idx):
        g2_0 = cont[ref_idx, 0]
        mu_0 = target_muF_by_T[ref_idx]
        self.__plot_pert_run_starting_from_reference(ax, g2_0, mu_0, target_muF_by_T, color, "mu0_"+lpd.format_float(mu_0,2)+"_pertrun")

    def __plot_various_pert_runs(self, ax):
        mask = max_muF_by_T_for_running >= self.cont_data_container.muF_by_T_cont
        target_muF_by_T = self.cont_data_container.muF_by_T_cont[mask]
        cont = self.cont_data_container.g2_cont[mask]
        self.__nonpert_ref_with_pert_run_wrapper(ax, cont, target_muF_by_T, "tab:cyan", ref_idx=0)  # smallest mu

        ref_idx = numpy.fabs(self.cont_data_container.muF_by_T_cont - min_muF_by_T_in_flow_extr).argmin()
        self.__nonpert_ref_with_pert_run_wrapper(ax, cont, target_muF_by_T, "k", ref_idx)  # largest mu
        self.__plot_pert_g2(ax, numpy.geomspace(min_muF_by_T_in_flow_extr, target_muF_by_T[-1], 200))

    def __plot_pert_run_from_converted_g2MSbar(self, ax):

        for reference_muF_by_T, color in zip((4, 8), ("black", "red")):

            mask = max_muF_by_T_for_running >= self.msbar_data_container.mubar_by_T
            target_muF_by_T = self.msbar_data_container.mubar_by_T[mask]

            ref_idx = numpy.fabs(self.msbar_data_container.mubar_by_T - reference_muF_by_T).argmin()
            g2_0 = self.msbar_data_container.g2_MSBAR[ref_idx]
            mu_0 = self.msbar_data_container.mubar_by_T[ref_idx]

            self.__plot_pert_run_starting_from_reference(ax, g2_0, mu_0, target_muF_by_T, color, "mu0_"+lpd.format_float(mu_0,2)+"_pertrun")

    def __plot_MSBAR(self, ax):
        ax.errorbar(*self.msbar_data_container, fmt='-', markersize=0, lw=self.linewidth, label=r'$\text{flow}\rightarrow\overline{\text{MS}}$')

    def plot_full(self, file_suffix="_all"):
        fig, ax = self.__setup_plot()

        self.__plot_grey_flow_extr_band(ax)
        self.__plot_lattice_data(ax)
        self.__plot_cont(ax)
        self.__plot_various_pert_runs(ax)
        self.__plot_MSBAR(ax)

        self.__finalize_plot(fig, ax, file_suffix)

    def plot_selected(self, file_suffix=""):
        fig, ax = self.__setup_plot()

        self.__plot_grey_flow_extr_band(ax)
        self.__plot_cont(ax)
        self.__plot_MSBAR(ax)
        self.__plot_pert_run_from_converted_g2MSbar(ax)
        # self.__plot_pert_g2(ax, numpy.geomspace(2, muB_by_T, 200))

        self.__finalize_plot(fig, ax, file_suffix)


def extrapolation_ansatz(x, m, b):
    return m * x + + b


def chisqdof(fitparams, xdata, ydata, edata):  # remove edata?
    return numpy.sum(((ydata - extrapolation_ansatz(xdata, *fitparams))/edata)**2) / (len(ydata)-len(fitparams))


def fit_sample(ydata, xdata, edata):
    fitparams = scipy.optimize.minimize(chisqdof, x0=numpy.asarray([-1, 2]), args=(xdata, ydata, edata))
    fitparams = fitparams.x

    this_chisqdof = chisqdof(fitparams, xdata, ydata, edata)

    return [*fitparams, this_chisqdof]


@lpd.typed_frozen_data
class ContDataContainer:
    muF_by_T_cont: NDArray[Shape["*"], Float64]
    g2_latt_from_int: NDArray[Shape["* nflow, * nNts"], Float64]
    g2_latt_from_int_err: NDArray[Shape["* nflow, * nNts"], Float64]
    g2_cont: NDArray[Shape["* nflow, [g2, g2err]"], Float64]
    g2_extr_slope: NDArray[Shape["* nflow, [slope, slope_err]"], Float64]


def do_cont_extr(args, raw_data_container):
    muF_by_T_cont = numpy.sort(numpy.unique(numpy.concatenate((
        numpy.geomspace(raw_data_container.largest_min_muF_by_T, max_muF_by_T_for_running, 200),
        necessary_muF_by_T,
        [muB_by_T, ]))))
    nflow = len(muF_by_T_cont)

    nNts = len(args.Nts)
    g2_latt_from_int = numpy.empty((nflow, nNts))
    g2_latt_from_int_err = numpy.empty((nflow, nNts))
    cont = numpy.empty((nflow, 2))
    slope = numpy.empty((nflow, 2))
    chisqdof = numpy.empty((nflow, 2))

    for i, muF_by_T in enumerate(muF_by_T_cont):
        g2_latt_from_int[i] = [spline(muF_by_T) for spline in raw_data_container.g2_ints]
        g2_latt_from_int_err[i] = [spline(muF_by_T) for spline in raw_data_container.g2_err_ints]

    output_file = args.outputpath_data + "g2_muF" + "_cont_extr.txt"

    if args.calc_cont:
        for i in range(nflow):
            fitparams, fitparams_err = bootstr.bootstr_from_gauss(fit_sample, data=g2_latt_from_int[i], data_std_dev=g2_latt_from_int_err[i], numb_samples=100,
                                                                  sample_size=1, return_sample=False,
                                                                  args=[1 / numpy.asarray(args.Nts) ** Ntexp * factor, g2_latt_from_int_err[i]],
                                                                  parallelize=True, nproc=20)
            chisqdof[i, 0] = fitparams[2]
            chisqdof[i, 1] = fitparams_err[2]
            cont[i, 0] = fitparams[1]
            cont[i, 1] = fitparams_err[1]
            slope[i, 0] = fitparams[0]
            slope[i, 1] = fitparams_err[0]

        columns = (muF_by_T_cont, cont, slope, chisqdof)
        labels = ["mu_F/T", "g2", "err", "slope", "err", "chisqdof", "err"]

        lpd.save_columns_to_file(output_file, columns, labels)

    else:
        try:
            _, b, berr, m, merr, chisqdof, chisqdof_err = numpy.loadtxt(output_file, unpack=True)
        except OSError:
            print("Error: could not find ", output_file)
            exit(1)
        for i in range(nflow):
            cont[i, 0] = b[i]
            cont[i, 1] = berr[i]
            slope[i, 0] = m[i]
            slope[i, 1] = merr[i]

    print("max relative error in %:", numpy.amax(g2_latt_from_int_err / g2_latt_from_int) * 100)

    return ContDataContainer(muF_by_T_cont, g2_latt_from_int, g2_latt_from_int_err, cont, slope)


def flow_coupling(tauF2E: float, a2bytauF: float) -> float:
    # some constants
    prefactor = numpy.pi**2 / 8
    c_0 = 3/128
    # c_2 = -1/256
    # c_4 = 343/327680
    c_2 = 0
    c_4 = 0

    flow_coupling = prefactor * tauF2E / (c_0 + c_2 * a2bytauF + c_4 * a2bytauF ** 2)

    return flow_coupling


@lpd.typed_frozen_data
class RawDataContainer:
    muF_by_T_arr: NDArray[Shape["*, *"], Float64]
    g2_arr: NDArray[Shape["*, *"], Float64]
    g2_err_arr: NDArray[Shape["*, *"], Float64]
    g2_ints: List[scipy.interpolate._fitpack2.InterpolatedUnivariateSpline]
    g2_err_ints: List[scipy.interpolate._fitpack2.InterpolatedUnivariateSpline]
    largest_min_muF_by_T: float


def load_data_and_interpolate(args: argparse.Namespace) -> RawDataContainer:
    def virtual_nt(beta):
        r0Tc = 0.7457

        # below we implicitly use the scale setting of our finite temp lattices. Note that we do not account for the slight mistuning of the temperatures
        # (i.e., we set T=1.5Tc for all lattices instead of 1.47, 1.51...) since what we actually want here is to continuum extrapolate g^2 at a fixed scale.
        TbyTc = 1.5

        r0_T = r0Tc * TbyTc
        inv_T_a = rS.r0_div_a(
            beta) / r0_T  # 1/(a T). a is lattice spacing. this is a virtual Nt (= the Nt, that we should have chosen to get T=1.5Tc while keeping beta fixed.)
        # print(inv_T_a)
        return inv_T_a

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
        tauFbya2, tauF2E, tauF2E_err = numpy.loadtxt(args.input_basepath + "/" + args.input_files[i], unpack=True)

        # convert to a more intuitive scale than lattice spacing. It's wrong to simply use the Nt of the finite temp lattices here to get, e.g. tauFT^2, because
        # then you actually have a slightly different scale for each lattice since there is some slight temperature mistuning in the scale setting.
        # for the continuum extrapolation you want to have a fixed scale.
        # Choices:
        # tauF/r0**2 or tauF/t0 or taufF_T^2 (via r0Tc and T/Tc==1.500)
        tauFbyScale = tauFbya2 / virtual_nt(args.betas[i])**2
        muF_by_T = numpy.flip(convert_taufT2_to_muFByT(tauFbyScale))
        muF_by_T_arr.append(muF_by_T)

        largest_min_muF_by_T = max(largest_min_muF_by_T, numpy.min(muF_by_T))

        # coupling
        thisg2 = numpy.flip(flow_coupling(tauF2E, 1/tauFbya2))
        thisg2_err = numpy.flip(numpy.fabs(flow_coupling(tauF2E_err, 1/tauFbya2)))
        g2_arr.append(thisg2)
        g2_err_arr.append(thisg2_err)

    for i in range(len(args.input_files)):
        muF_by_T = muF_by_T_arr[i]
        thisg2 = g2_arr[i]
        thisg2_err = g2_err_arr[i]
        min_idx, max_idx = get_min_and_max_indices(muF_by_T, min_muF_by_T_for_running, max_muF_by_T_for_running)
        # print(min_idx, max_idx, min_muF_by_T_in_plot, max_muF_by_T_in_plot, muF_by_T[max_idx])
        g2_ints.append(scipy.interpolate.InterpolatedUnivariateSpline(muF_by_T[0:max_idx], thisg2[0:max_idx], k=order, ext=2))
        g2_err_ints.append(scipy.interpolate.InterpolatedUnivariateSpline(muF_by_T[0:max_idx], thisg2_err[0:max_idx], k=order, ext=2))

    return RawDataContainer(numpy.asarray(muF_by_T_arr), numpy.asarray(g2_arr), numpy.asarray(g2_err_arr), g2_ints, g2_err_ints, largest_min_muF_by_T)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('--calc_cont', help='calc continuum extrapolation and save to file instead of reading it from the file', action="store_true")
    parser.add_argument('--input_basepath', default="/work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/coupling/")
    parser.add_argument('--input_files', nargs='*', type=str, required=True)
    parser.add_argument('--outputpath_plot', default="/work/home/altenkort/work/correlators_flow/plots/quenched_1.50Tc_zeuthenFlow/coupling/")
    parser.add_argument('--outputpath_data', default="/work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/coupling/")
    parser.add_argument('--Nts', nargs='*', type=int)
    parser.add_argument('--betas', nargs='*', type=float)
    args = parser.parse_args()

    if len(args.betas) != len(args.input_files) or len(args.Nts) != len(args.input_files):
        print("ERROR: --betas and --input_files must have same number of arguments.")
        exit(1)

    return args


def convert_g2_flow_to_MSBAR(g2_flow: float, nf: int = 0) -> float:
    # Harlander, R.V., Neumann, T. The perturbative QCD gradient flow to three loops. J. High Energ. Phys. 2016, 161 (2016).
    # https://doi.org/10.1007/JHEP06(2016)161

    k2 = -0.982 - 0.070 * nf + 0.002 * nf**2
    k1 = 1.098 + 0.008 * nf

    alpha_flow = g2_flow / (4 * numpy.pi)

    coefficients = [k2, k1, 1, -alpha_flow]

    solutions = numpy.roots(coefficients)

    # Filter the real solutions
    real_solution = [sol*4*numpy.pi for sol in solutions if numpy.isreal(sol) and sol >= 0 and sol*4*numpy.pi <= 5][0]

    return real_solution


@lpd.typed_frozen_data
class MSBarDataContainer:
    mubar_by_T: NDArray[Shape["*"], Float64]
    g2_MSBAR: NDArray[Shape["*"], Float64]


def calc_g2_MSBAR(cont_data_container: ContDataContainer) -> MSBarDataContainer:

    g2_MSBAR_arr = []

    for g2flow in cont_data_container.g2_cont:
        g2_MSBAR = convert_g2_flow_to_MSBAR(g2flow[0])
        g2_MSBAR_arr.append(g2_MSBAR)

    return MSBarDataContainer(cont_data_container.muF_by_T_cont, numpy.asarray(g2_MSBAR_arr))


def main():
    # this script load tf^2 <E> from some files, then calculates the tree-level improved coupling g^2 and does a continuum extrapolation of g^2

    args = parse_args()

    raw_data_container = load_data_and_interpolate(args)
    cont_data_container = do_cont_extr(args, raw_data_container)
    MSBar_data_container = calc_g2_MSBAR(cont_data_container)

    plotter_g2_vs_mu = Plotter_g2_vs_mu(args, raw_data_container, cont_data_container, MSBar_data_container)
    plotter_g2_vs_mu.plot_full()
    plotter_g2_vs_mu.plot_selected()

    plot_g2_vs_1overNt2(args, cont_data_container)


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
