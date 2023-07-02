#!/usr/bin/env python3
import lib_process_data as lpd
import numpy
from latqcdtools.statistics import bootstr
import argparse
import scipy.optimize
import scipy.interpolate
import latqcdtools.physics.referenceScales as sq

try:
    import colored_traceback.auto
except ImportError:
    pass

Ntexp=2
factor = 96**Ntexp


def convert_sqrt8taufTByTau_to_muFByT(sqrt8taufTByTau, tauT):
    muFByT = 1/(sqrt8taufTByTau*tauT)
    return muFByT


def get_necessary_muF_by_T():
    necessary_muF_by_T = numpy.asarray([])

    for tauT in numpy.arange(0.25, 0.51, 1 / 36):
        tmp = convert_sqrt8taufTByTau_to_muFByT(lpd.get_relflow_range(), tauT)
        necessary_muF_by_T = numpy.concatenate((necessary_muF_by_T, tmp))

    return necessary_muF_by_T


min_muF_by_T_in_extr = convert_sqrt8taufTByTau_to_muFByT(sqrt8taufTByTau=0.3, tauT=0.5)
max_muF_by_T_in_extr = convert_sqrt8taufTByTau_to_muFByT(sqrt8taufTByTau=0.25, tauT=0.25)

necessary_muF_by_T = get_necessary_muF_by_T()
max_muF_by_T_in_plot = numpy.max(necessary_muF_by_T)
min_muF_by_T_in_plot = numpy.min(necessary_muF_by_T)

muB_by_T = 19.179
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


def plot2(args, data, data_err, cont, slope, muF_by_T_samples, plotlabel, label):

    # plot continuum extrapolation at different flow times
    xlabel = r'$N_\tau^{-'+str(Ntexp)+r'}$'
    fig3, ax3, _ = lpd.create_figure(xlabel=xlabel, ylabel=get_ylabel())

    xpoints = numpy.linspace(0, 1 / args.Nts[-1] ** Ntexp, 10)

    plot_muF_by_T = numpy.geomspace(min_muF_by_T_in_extr, muB_by_T, 10)

    for j, muF_by_T in enumerate(plot_muF_by_T):
        i = (numpy.fabs(muF_by_T_samples - muF_by_T)).argmin()
        color = lpd.get_color(plot_muF_by_T, j)
        ax3.errorbar(numpy.insert(1 / numpy.asarray(args.Nts)**Ntexp, 0, 0),
                     numpy.insert(data[i], 0, cont[i, 0]),
                     numpy.insert(data_err[i], 0, cont[i, 1]),
                     color=color,
                     fmt='|', label=lpd.format_float(muF_by_T_samples[i], 2), zorder=-i)
        ax3.errorbar(xpoints, extrapolation_ansatz(xpoints, slope[i, 0] * factor, cont[i, 0]),
                     **lpd.fitlinestyle, color=color)
        # counter += 1
    ax3.legend(title=r'$\displaystyle \mu_\mathrm{F}/ ' + plotlabel + r'$', loc="center left", bbox_to_anchor=(1, 0.5), **lpd.leg_err_size())

    # thisylims = ax3.get_ylim()
    # ax3.set_ylim((0, thisylims[1]*1.1))
    ax3.set_ylim(0.75, 2.85)
    thisxlims = ax3.get_xlim()
    ax3.set_xlim((thisxlims[0], thisxlims[1] * 1.025))

    ax3.figure.canvas.draw()
    offset = ax3.xaxis.get_major_formatter().get_offset()
    print(offset)
    ax3.xaxis.offsetText.set_visible(False)
    ax3.xaxis.set_label_text(xlabel + " " + offset)

    file = args.outputpath_plot + "g2_" + label + "_cont_extr_imp.pdf"
    print("saving", file)
    fig3.savefig(file)


def plot1(args, muF_by_T_arr, g2_arr, g2_err_arr, muF_by_T_samples, cont, plotlabel, label):

    def plot_pert_g2(ax, mus, suffix="pert"):
        g2s = []
        for mu in mus:
            g2s.append(lpd.get_g2_pert(mu*T_in_GeV, Nf=0))
        ax.errorbar(mus, g2s, fmt='-.', lw=0.5,  markersize=0, label=r'pert. ($\bar{\mu} = \mu_\mathrm{F} \equiv 1/\sqrt{8\tau_F} $)', zorder=-1)
        lpd.save_columns_to_file(args.outputpath_data + "g2_muF_" + suffix + ".txt", (mus, g2s), ["mu_F/T", "g2s"])

    def plot_nonpert_g2_with_pert_running(ax, g2_0, mu0, mus, color, suffix):
        alpha_s0 = g2_0 / (4 * numpy.pi)
        import rundec
        crd = rundec.CRunDec()
        g2s = []
        nf = 0
        nloop = 5
        plot_mus = []
        for mu in mus:
            if mu >= mu0:
                Alphas = crd.AlphasExact(alpha_s0, mu0, mu, nf, nloop)
                g2 = 4. * numpy.pi * Alphas
                g2s.append(g2)
                plot_mus.append(mu)
        plot_mus = numpy.asarray(plot_mus)
        ax.errorbar(plot_mus, g2s, fmt='--', markersize=0, lw=0.75,  label='cont. + pert. run', zorder=-1, color=color)
        lpd.save_columns_to_file(args.outputpath_data + "g2_muF_" + suffix + ".txt", (plot_mus, g2s), ["mu_F/T", "g2s"])

    def plot_lattice_data(ax):
        for i, Nt in enumerate(args.Nts):
            ax.errorbar(muF_by_T_arr[i], g2_arr[i], fmt='-', lw=0.75, markersize=0, label=r'$N_\tau = ' + str(Nt) + r'$', zorder=-i - 10)
        ax.fill_between([min_muF_by_T_in_extr, max_muF_by_T_in_extr],
                         [-1, -1], [100, 100], facecolor='k', alpha=0.15, zorder=-1000)

    def plot_cont(ax2, cont, muF_by_T_samples):
        cont_data = cont[:, 0]
        ax2.errorbar(muF_by_T_samples, cont_data, fmt='-', markersize=0, lw=0.75, label=r'cont., linear in $a^' + str(Ntexp) + r'$', zorder=-1)

    def nonpert_ref_with_pert_run_wrapper(ax2, cont, muF_by_T_samples, color, ref_idx):
        g2_0 = cont[ref_idx, 0]
        mu_0 = muF_by_T_samples[ref_idx]
        plot_nonpert_g2_with_pert_running(ax2, g2_0, mu_0, muF_by_T_samples, color, "mu0_"+lpd.format_float(mu_0,2)+"_pertrun")

    # plot coupling for all lattices
    fig2, ax2, _ = lpd.create_figure(xlabel=r'$\displaystyle \mu_{\mathrm{F}} / ' + plotlabel + r'$', ylabel=get_ylabel())
    ax2.set_ylim(0, 4.5)
    ax2.set_xlim(1, 100)
    ax2.set_xscale('log')

    plot_lattice_data(ax2)

    mask = muB_by_T >= muF_by_T_samples
    target_muF_by_T = muF_by_T_samples[mask]
    cont = cont[mask]

    plot_cont(ax2, cont, target_muF_by_T)
    nonpert_ref_with_pert_run_wrapper(ax2, cont, target_muF_by_T, "tab:cyan", ref_idx=0)  # smallest mu
    nonpert_ref_with_pert_run_wrapper(ax2, cont, target_muF_by_T, "k", ref_idx=numpy.fabs(muF_by_T_samples - min_muF_by_T_in_extr).argmin())  # largest mu
    plot_pert_g2(ax2, numpy.geomspace(min_muF_by_T_in_extr, target_muF_by_T[-1], 200))

    ax2.axvline(muB_by_T, alpha=1, dashes=(2, 2), zorder=-10000, lw=0.5, color='k')

    ax2.legend(**lpd.leg_err_size(), loc="lower left", bbox_to_anchor=(0, 0), handlelength=1, fontsize=8, framealpha=0)

    file = args.outputpath_plot + "g2_" + label + "_imp.pdf"
    print("saving", file)
    fig2.savefig(file)


def extrapolation_ansatz(x, m, b):
    return m * x + + b


def chisqdof(fitparams, xdata, ydata, edata):  # remove edata?
    return numpy.sum(((ydata - extrapolation_ansatz(xdata, *fitparams))/edata)**2) / (len(ydata)-len(fitparams))


def fit_sample(ydata, xdata, edata):
    fitparams = scipy.optimize.minimize(chisqdof, x0=numpy.asarray([-1, 2]), args=(xdata, ydata, edata))
    fitparams = fitparams.x

    this_chisqdof = chisqdof(fitparams, xdata, ydata, edata)

    return [*fitparams, this_chisqdof]


def do_cont_extr(args, g2_ints, g2_err_ints, output_file, largest_min_muF_by_T):

    muF_by_T_samples = numpy.geomspace(largest_min_muF_by_T, max_muF_by_T_in_plot, 200)
    muF_by_T_samples = numpy.sort(numpy.unique(numpy.concatenate((muF_by_T_samples, necessary_muF_by_T, [muB_by_T, ]))))
    nflow = len(muF_by_T_samples)

    nNts = len(args.Nts)
    data = numpy.empty((nflow, nNts))
    data_err = numpy.empty((nflow, nNts))
    cont = numpy.empty((nflow, 2))
    slope = numpy.empty((nflow, 2))
    chisqdof = numpy.empty((nflow, 2))

    for i, muF_by_T in enumerate(muF_by_T_samples):
        data[i] = [spline(muF_by_T) for spline in g2_ints]
        data_err[i] = [spline(muF_by_T) for spline in g2_err_ints]

    if args.calc_cont:
        for i in range(nflow):
            fitparams, fitparams_err = bootstr.bootstr_from_gauss(fit_sample, data=data[i], data_std_dev=data_err[i], numb_samples=100,
                                                                  sample_size=1, return_sample=False,
                                                                  args=[1 / numpy.asarray(args.Nts) ** Ntexp * factor, data_err[i]],
                                                                  parallelize=True, nproc=20)
            chisqdof[i, 0] = fitparams[2]
            chisqdof[i, 1] = fitparams_err[2]
            cont[i, 0] = fitparams[1]
            cont[i, 1] = fitparams_err[1]
            slope[i, 0] = fitparams[0]
            slope[i, 1] = fitparams_err[0]

        columns = (muF_by_T_samples, cont, slope, chisqdof)
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

    return muF_by_T_samples, data, data_err, cont, slope


def flow_coupling(tauF2E, a2bytauF):
    # some constants
    prefactor = numpy.pi**2 / 8
    c_0 = 3/128
    # c_2 = -1/256
    # c_4 = 343/327680
    c_2 = 0
    c_4 = 0

    flow_coupling = prefactor * tauF2E / (c_0 + c_2 * a2bytauF + c_4 * a2bytauF ** 2)

    return flow_coupling


def load_data_and_interpolate(args, scalefunc):
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
        tauFbyScale = tauFbya2 / scalefunc(args.betas[i])**2
        muF_by_T = numpy.flip(1/numpy.sqrt(8*tauFbyScale))
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
        min_idx, max_idx = get_min_and_max_indices(muF_by_T, min_muF_by_T_in_plot, max_muF_by_T_in_plot)
        # print(min_idx, max_idx, min_muF_by_T_in_plot, max_muF_by_T_in_plot, muF_by_T[max_idx])
        g2_ints.append(scipy.interpolate.InterpolatedUnivariateSpline(muF_by_T[0:max_idx], thisg2[0:max_idx], k=order, ext=2))
        g2_err_ints.append(scipy.interpolate.InterpolatedUnivariateSpline(muF_by_T[0:max_idx], thisg2_err[0:max_idx], k=order, ext=2))

    return muF_by_T_arr, g2_arr, g2_err_arr, g2_ints, g2_err_ints, largest_min_muF_by_T


def virtual_nt(beta):
    r0Tc = 0.7457

    # below we implicitly use the scale setting of our finite temp lattices. Note that we do not account for the slight mistuning of the temperatures
    # (i.e., we set T=1.5Tc for all lattices instead of 1.47, 1.51...) since what we actually want here is to continuum extrapolate g^2 at a fixed scale.
    TbyTc = 1.5

    r0_T = r0Tc * TbyTc
    inv_T_a = sq.r0_div_a(beta) / r0_T  # 1/(a T). a is lattice spacing. this is a virtual Nt (= the Nt, that we should have chosen to get T=1.5Tc while keeping beta fixed.)
    # print(inv_T_a)
    return inv_T_a


def get_scale_params(scale):
    if scale == "T_via_r0Tc":
        plotlabel = r'T'
        label = ""
        scalefunc = virtual_nt
    else:
        print("Scales other than T_via_r0Tc are deprecated")
        exit(1)
    # if scale == "r0":
    #     plotlabel = r'r_0'
    #     label = "r0"
    #     scalefunc = sq.r0_div_a
    # elif scale == "t0":
    #     plotlabel = r'\sqrt{t_0}'
    #     label = "t0"
    #     scalefunc = sq.sqrtt0

    return plotlabel, label, scalefunc


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--calc_cont', help='calc continuum extrapolation and save to file instead of reading it from the file', action="store_true")
    parser.add_argument('--ref_scale', choices=["r0", "t0", "T_via_r0Tc"], required=True)
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


def main():

    # this script load tf^2 <E> from some files, then calculates the tree-level improved coupling g^2 and does a continuum extrapolation of g^2

    args = parse_args()
    plotlabel, label, scalefunc = get_scale_params(args.ref_scale)
    output_file = args.outputpath_data + "g2_muF" + label + "_cont_extr.txt"
    muF_by_T_arr, g2_arr, g2_err_arr, g2_ints, g2_err_ints, largest_min_muF_by_T = load_data_and_interpolate(args, scalefunc)
    muF_by_T_samples, data, data_err, cont, slope = do_cont_extr(args, g2_ints, g2_err_ints, output_file, largest_min_muF_by_T)

    # TODO. now that plot is correct, extract the four possibilities of g^2 before plotting, save it to file, and then pass them to the plot function.
    # I need g^2 as a function of mu_F/T in the range min_muF_by_T_in_extr to muB_by_T

    plot1(args, muF_by_T_arr, g2_arr, g2_err_arr, muF_by_T_samples, cont, plotlabel, label)
    plot2(args, data, data_err, cont, slope, muF_by_T_samples, plotlabel, label)


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
