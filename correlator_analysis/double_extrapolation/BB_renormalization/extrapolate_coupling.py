#!/usr/bin/env python3
import lib_process_data as lpd
import numpy
from latqcdtools.statistics import bootstr
import argparse
import scipy.optimize
import scipy.interpolate
import latqcdtools.physics.referenceScales as sq


Ntexp=2
factor = 96**Ntexp


def extrapolation_ansatz(x, m, b):
    return m * x + + b


def chisqdof(fitparams, xdata, ydata, edata):  # remove edata?
    return numpy.sum(((ydata - extrapolation_ansatz(xdata, *fitparams))/edata)**2) / (len(ydata)-len(fitparams))


def fit_sample(ydata, xdata, edata):
    fitparams = scipy.optimize.minimize(chisqdof, x0=numpy.asarray([-1, 2]), args=(xdata, ydata, edata))
    fitparams = fitparams.x

    this_chisqdof = chisqdof(fitparams, xdata, ydata, edata)

    return [*fitparams, this_chisqdof]


def load_data(args, scalefunc):
    # declare a few arrays
    tauFbyScale_arr = []
    g2_arr = []
    g2_err_arr = []
    g2_ints = []
    g2_err_ints = []

    # some constants
    prefactor = numpy.pi**2 / 8
    c_0 = 3/128
    # c_2 = -1/256
    # c_4 = 343/327680
    c_2 = 0
    c_4 = 0
    order = 3

    # load data and save various "versions" of it
    for i in range(len(args.input_files)):
        tauFbya2, tauF2E, tauF2E_err = numpy.loadtxt(args.input_basepath + "/" + args.input_files[i], unpack=True)

        a2bytauF = 1/tauFbya2

        # flow times
        tauFbyScale_arr.append(tauFbya2 / scalefunc(args.betas[i])**2)  # tauF/r0**2 or tauF/t0 or taufF_T^2
        # tauFbyScale_arr.append(tauFbya2 / args.Nts[i]**2)  # tauF/r0**2 or tauF/t0

        # coupling
        thisg2 = prefactor * tauF2E / (c_0 + c_2 * a2bytauF + c_4 * a2bytauF**2)
        thisg2_err = numpy.fabs(prefactor * tauF2E_err / (c_0 + c_2 * a2bytauF + c_4 * a2bytauF**2))
        g2_arr.append(thisg2)
        g2_err_arr.append(thisg2_err)
        g2_ints.append(scipy.interpolate.InterpolatedUnivariateSpline(tauFbyScale_arr[-1], g2_arr[-1], k=order, ext=2))
        g2_err_ints.append(scipy.interpolate.InterpolatedUnivariateSpline(tauFbyScale_arr[-1], g2_err_arr[-1], k=order, ext=2))

    return tauFbyScale_arr, g2_arr, g2_err_arr, g2_ints, g2_err_ints


def temperature(beta):
    r0Tc = 0.7457
    TbyTc = 1.5
    r0_T = r0Tc * TbyTc
    inv_T_a = sq.r0_div_a(beta) / r0_T  # 1/(a T). a is lattice spacing.
    return inv_T_a


def get_scale_params(scale):
    if scale == "r0":
        plotlabel = r'r_0'
        label = "r0"
        scalefunc = sq.r0_div_a
    elif scale == "T_via_r0Tc":
        plotlabel = r'T'
        label = "T_via_r0Tc"
        scalefunc = temperature
    elif scale == "t0":
        plotlabel = r'\sqrt{t_0}'
        label = "t0"
        scalefunc = sq.sqrtt0

    return plotlabel, label, scalefunc


def get_ylabel():
    # return r'$ g^{2 \text{(flow)}} = \frac{\pi^2}{8} \tau_F^2 \langle E \rangle / \frac{3}{128}$'
    # return r'$ \scriptstyle g^{2 \text{(flow)}} = \frac{\pi^2\tau_F^2 \langle E\rangle}{8}  / \left( \frac{3}{128} + \scriptstyle \underset{n=1,2}{\sum} \scriptscriptstyle  c_{2n} \frac{a^{2n}}{\tau_\mathrm{F}^{n}} \right)$'
    return r'$ \displaystyle g^{2}$'


def convert_taufT2_to_sqrt8tauFByr0(taufT2, TbyTc=1.5, r0Tc=0.7457):
    return numpy.sqrt(8*taufT2)/(TbyTc*r0Tc)


def convert_sqrt8taufT_to_sqrt8tauFByr0(sqrt8taufT, TbyTc=1.5, r0Tc=0.7457):
    return sqrt8taufT / (TbyTc*r0Tc)


def plot1(args, tauFbyScale_arr, g2_arr, g2_err_arr, cont, flowstart, flowend, plotlabel, label):

    def plot_pert_g2(ax, mus):
        g2s = []
        for mu in mus:
            g2s.append(lpd.get_g2_pert(mu, Nf=0))

        ax.errorbar(mus/0.472, g2s, fmt='-.', lw=0.5,  markersize=0, label=r'pert. ($\bar{\mu} = \mu_\mathrm{F} \equiv 1/\sqrt{8\tau_F} $)', zorder=-1)

    def plot_nonpert_g2_with_pert_running(ax, g2_0, mu0, mus, color):
        alpha_s0 = g2_0 / (4 * numpy.pi)
        import rundec
        crd = rundec.CRunDec()
        g2s = []
        nf = 0
        nloop = 5
        plot_mus = []
        for mu in mus:
            if mu > mu0:
                Alphas = crd.AlphasExact(alpha_s0, mu0, mu, nf, nloop)
                g2 = 4. * numpy.pi * Alphas
                g2s.append(g2)
                plot_mus.append(mu)

        ax.errorbar(numpy.asarray(plot_mus)/0.472, g2s, fmt='--', markersize=0, lw=0.75,  label='cont. + pert. run', zorder=-1, color=color)
        # ax.axvline(mu0/0.472, alpha=1, dashes=(2,2), zorder=-10000, lw=0.5, color=color)

    def plot_lattice_data(ax):
        for i, Nt in enumerate(args.Nts):
            muF_scale = 1 / numpy.sqrt(8 * tauFbyScale_arr[i])
            ax.errorbar(muF_scale, g2_arr[i], fmt='-', lw=0.75, markersize=0, label=r'$N_\tau = ' + str(Nt) + r'$', zorder=-i - 10)
        ax.fill_between([1 / (0.5 * 0.3), 1 / (0.25 * 0.25)],
                         [-1, -1], [4, 4], facecolor='k', alpha=0.15, zorder=-1000)

    # plot coupling for all lattices
    fig2, ax2, _ = lpd.create_figure(xlabel=r'$\displaystyle \mu_{\mathrm{F}} / ' + plotlabel + r'$', ylabel=get_ylabel())
    ax2.set_ylim(0, 4.5)
    ax2.set_xlim(1, 900)

    mu_inv = numpy.sqrt(8*tauFbyScale_arr[0])

    cont_data = cont[:,0]
    muF_scale_cont = 1/numpy.sqrt(8*tauFbyScale_arr[0][:len(cont_data)])
    ax2.errorbar(muF_scale_cont, cont[:, 0], fmt='-', markersize=0, lw=0.75, label=r'cont., linear in $a^'+str(Ntexp)+r'$', zorder=-1)

    plot_lattice_data(ax2)

    # ax2.set_xticks((0, 0.05, 0.1, 0.15, 0.2, 0.25))

    plot_pert_g2(ax2, numpy.logspace(0, 4, 200))

    # smallest mu
    ref_idx = min(numpy.abs(mu_inv - 0.4).argmin(), len(cont[:, 0])-1)
    g2_0 = cont[ref_idx, 0]
    mu_0 = muF_scale_cont[ref_idx]*0.472
    plot_nonpert_g2_with_pert_running(ax2, g2_0, mu_0, muF_scale_cont*0.472, "tab:cyan")

    # largest mu
    ref_idx = numpy.abs(mu_inv - 0.5*0.3).argmin()
    g2_0 = cont[ref_idx, 0]
    mu_0 = muF_scale_cont[ref_idx]*0.472
    plot_nonpert_g2_with_pert_running(ax2, g2_0, mu_0, muF_scale_cont*0.472, "k")

    ax2.set_xscale('log')

    ax2.legend(**lpd.leg_err_size(), loc="upper right", bbox_to_anchor=(1, 1), handlelength=1, fontsize=8, framealpha=0)

    file = args.outputpath_plot + "g2_" + label + "_imp.pdf"
    print("saving", file)
    fig2.savefig(file)


def plot2(args, data, data_err, cont, slope, tauFbyScale_arr, flowstart, flowend, plotlabel, label):
    # plot continuum extrapolation at different flow times
    nflow = flowend - flowstart

    xlabel = r'$N_\tau^{-'+str(Ntexp)+r'}$'

    fig3, ax3, _ = lpd.create_figure(xlabel=xlabel, ylabel=get_ylabel())

    min_idx = numpy.abs(numpy.sqrt(8*tauFbyScale_arr[0])-convert_sqrt8taufT_to_sqrt8tauFByr0(0.25*0.25)).argmin()-1
    max_idx = numpy.abs(numpy.sqrt(8*tauFbyScale_arr[0])-convert_sqrt8taufT_to_sqrt8tauFByr0(0.5*0.3)).argmin()+1

    counter = 0
    xpoints = numpy.linspace(0, 1 / args.Nts[-1] ** Ntexp, 10)
    for i in range(flowend - 1, flowstart, -1):
        if (counter % 5 == 0 and min_idx < i < max_idx) or (i == min_idx) or (i == max_idx):
            color = lpd.get_color(range(nflow), i, min_idx, max_idx, 0.95)
            ax3.errorbar(numpy.insert(1 / numpy.asarray(args.Nts)**Ntexp, 0, 0),
                         numpy.insert(data[i], 0, cont[i, 0]),
                         numpy.insert(data_err[i], 0, cont[i, 1]),
                         color=color,
                         fmt='|', label=lpd.format_float(numpy.sqrt(8*tauFbyScale_arr[0][i])), zorder=-i)
            ax3.errorbar(xpoints, extrapolation_ansatz(xpoints, slope[i, 0] * factor, cont[i, 0]),
                         **lpd.fitlinestyle, color=color)
        counter += 1
    ax3.legend(title=r'$\displaystyle \frac{\sqrt{8\tau_\mathrm{F}}}{' + plotlabel + r'}$', loc="center left", bbox_to_anchor=(1, 0.5), **lpd.leg_err_size())

    # thisylims = ax3.get_ylim()
    # ax3.set_ylim((0, thisylims[1]*1.1))
    ax3.set_ylim(0.8, 2.99)
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


def do_cont_extr(args, tauFbyScale_arr, g2_ints, g2_err_ints, output_file):
    nNts = len(args.Nts)
    smallest_last_value = numpy.min([arr[-1] for arr in tauFbyScale_arr])
    largest_initial_value = numpy.max([arr[0] for arr in tauFbyScale_arr])
    flowend = (numpy.abs(tauFbyScale_arr[0] - smallest_last_value)).argmin() - 1
    flowstart = (numpy.abs(tauFbyScale_arr[0] - largest_initial_value)).argmin() + 1
    nflow = flowend - flowstart
    data = numpy.empty((nflow, nNts))
    data_err = numpy.empty((nflow, nNts))
    cont = numpy.empty((nflow, 2))
    slope = numpy.empty((nflow, 2))
    chisqdof = numpy.empty((nflow, 2))

    for i in range(flowstart, flowend):
        data[i - flowstart] = [spline(tauFbyScale_arr[0][i]) for spline in g2_ints]
        data_err[i - flowstart] = [spline(tauFbyScale_arr[0][i]) for spline in g2_err_ints]

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

        print("saving", output_file)
        numpy.savetxt(output_file, numpy.column_stack((tauFbyScale_arr[0][flowstart:flowend], cont, slope, chisqdof)),
                      header="mu_F/T g2_cont    err     slope     err     chisqdof    err")
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

    return flowstart, flowend, data, data_err, cont, slope


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
    output_file = args.outputpath_data + "g2_" + label + "_cont_extr.txt"
    tauFbyScale_arr, g2_arr, g2_err_arr, g2_ints, g2_err_ints = load_data(args, scalefunc)
    flowstart, flowend, data, data_err, cont, slope = do_cont_extr(args, tauFbyScale_arr, g2_ints, g2_err_ints, output_file)
    plot1(args, tauFbyScale_arr, g2_arr, g2_err_arr, cont, flowstart, flowend, plotlabel, label)
    plot2(args, data, data_err, cont, slope, tauFbyScale_arr, flowstart, flowend, plotlabel, label)


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
