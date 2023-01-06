#!/usr/bin/env python3
import lib_process_data as lpd
import numpy
from latqcdtools.statistics import bootstr
import argparse
import scipy.optimize
import scipy.interpolate
import latqcdtools.physics.referenceScales as sq


def extrapolation_ansatz(x, m, b):
    return m * x + b


def fit_sample(ydata, xdata, edata, start_params=None):
    fitparams, _ = scipy.optimize.curve_fit(extrapolation_ansatz, xdata, ydata, p0=start_params, sigma=edata)
    return fitparams


def load_data(args, scalefunc):
    # declare a few arrays
    tauFbyScale_arr = []
    g2_arr = []
    g2_err_arr = []
    g2_ints = []
    g2_err_ints = []

    # some constants
    prefactor = 128 * numpy.pi ** 2 / 24
    order = 3

    # load data and save various "versions" of it
    for i in range(len(args.input_files)):
        tauFbya2, tauF2E, tauF2E_err = numpy.loadtxt(args.input_basepath + "/" + args.input_files[i], unpack=True)

        # flow times
        tauFbyScale_arr.append(tauFbya2 / scalefunc(args.betas[i]) ** 2)  # tauF/r0**2 or tauF/t0

        # coupling
        g2_arr.append(prefactor * tauF2E)
        g2_err_arr.append(prefactor * tauF2E_err)
        g2_ints.append(scipy.interpolate.InterpolatedUnivariateSpline(tauFbyScale_arr[-1], g2_arr[-1], k=order, ext=2))
        g2_err_ints.append(scipy.interpolate.InterpolatedUnivariateSpline(tauFbyScale_arr[-1], g2_err_arr[-1], k=order, ext=2))

    return tauFbyScale_arr, g2_arr, g2_err_arr, g2_ints, g2_err_ints


def get_scale_params(scale):
    if scale == "r0":
        plotlabel = r'r_0'
        label = "r0"
        scalefunc = sq.r0_div_a
    elif scale == "t0":
        plotlabel = r'\sqrt{t_0}'
        label = "t0"
        scalefunc = sq.sqrtt0

    return plotlabel, label, scalefunc


def get_ylabel():
    return r'$ g^{2 \text{(flow)}} = \frac{128\pi^2}{24} \tau_F^2 \langle E \rangle$'


def plot1(args, tauFbyScale_arr, g2_arr, g2_err_arr, cont, flowstart, flowend, plotlabel, label):
    # plot coupling for all lattices
    fig2, ax2, _ = lpd.create_figure(xlabel=r'$\frac{\sqrt{8\tau_\mathrm{F}}}{' + plotlabel + r'}$', ylabel=get_ylabel())
    # ax2.errorbar(numpy.sqrt(8*tauFbyScale_arr[0][flowstart:flowend]), cont[:, 0], cont[:, 1], fmt='-', markersize=0, label='cont', zorder=-1)
    ax2.errorbar(numpy.sqrt(8 * tauFbyScale_arr[0][flowstart:flowend]), cont[:, 0], fmt='-', markersize=0, label='cont', zorder=-1)
    # ax2.fill_between(numpy.sqrt(8 * tauFbyScale_arr[0][flowstart:flowend]), cont[:, 0]-cont[:, 1], cont[:, 0]+cont[:, 1], label='cont', zorder=-1)
    for i, Nt in enumerate(args.Nts):
        # ax2.fill_between(numpy.sqrt(8 * tauFbyScale_arr[i]), g2_arr[i]-g2_err_arr[i], g2_arr[i]+g2_err_arr[i], label=str(Nt), zorder=-i - 10)
        # ax2.errorbar(numpy.sqrt(8*tauFbyScale_arr[i]), g2_arr[i], g2_err_arr[i], fmt='-', markersize=0, label=str(Nt), zorder=-i - 10)
        ax2.errorbar(numpy.sqrt(8 * tauFbyScale_arr[i]), g2_arr[i], fmt='-', markersize=0, label=str(Nt), zorder=-i - 10)
    ax2.legend(title=r'$N_\tau$', **lpd.leg_err_size(), loc="center right", bbox_to_anchor=(1, 0.5), handlelength=1)
    file = args.outputpath_plot + "g2_" + label + ".pdf"
    print("saving", file)
    fig2.savefig(file)


def plot2(args, data, data_err, cont, slope, tauFbyScale_arr, flowstart, flowend, plotlabel, label):
    # plot continuum extrapolation at different flow times
    nflow = flowend - flowstart
    fig3, ax3, _ = lpd.create_figure(xlabel=r'$N_\tau^{-2}$', ylabel=get_ylabel())
    min_idx = 20
    xpoints = numpy.linspace(0, 1 / args.Nts[-1] ** 2, 10)
    for i in range(flowend - 1, flowstart, -1):
        if i % 10 == 0 and i >= min_idx:
            ax3.errorbar(numpy.insert(1 / numpy.asarray(args.Nts) ** 2, 0, 0),
                         numpy.insert(data[i], 0, cont[i, 0]),
                         numpy.insert(data_err[i], 0, cont[i, 1]),
                         color=lpd.get_color(range(nflow), i, min_idx, nflow - 1, 1),
                         fmt='|', label=lpd.format_float(numpy.sqrt(8*tauFbyScale_arr[0][i])), zorder=-i)
            ax3.errorbar(xpoints, extrapolation_ansatz(xpoints, slope[i, 0], cont[i, 0]),
                         fmt='--', alpha=0.5, color=lpd.get_color(range(nflow), i, min_idx, nflow-1, 1))
    ax3.legend(title=r'$\frac{\sqrt{8\tau_\mathrm{F}}}{' + plotlabel + r'}$', loc="center left", bbox_to_anchor=(1, 0.5), **lpd.leg_err_size())

    thisylims = ax3.get_ylim()
    ax3.set_ylim((0, thisylims[1]*1.1))
    thisxlims = ax3.get_xlim()
    ax3.set_xlim((thisxlims[0], thisxlims[1] * 1.025))

    file = args.outputpath_plot + "g2_" + label + "_cont_extr.pdf"
    print("saving", file)
    fig3.savefig(file)


def do_cont_extr(args, tauFbyScale_arr, g2_ints, g2_err_ints, output_file):
    nNts = len(args.Nts)
    smallest_last_value = numpy.min([arr[-1] for arr in tauFbyScale_arr])
    largest_initial_value = numpy.max([arr[0] for arr in tauFbyScale_arr])
    flowend = (numpy.abs(tauFbyScale_arr[0] - smallest_last_value)).argmin() + 1
    flowstart = (numpy.abs(tauFbyScale_arr[0] - largest_initial_value)).argmin() + 1
    nflow = flowend - flowstart
    data = numpy.empty((nflow, nNts))
    data_err = numpy.empty((nflow, nNts))
    cont = numpy.empty((nflow, 2))
    slope = numpy.empty((nflow, 2))

    for i in range(flowstart, flowend):
        data[i - flowstart] = [spline(tauFbyScale_arr[0][i]) for spline in g2_ints]
        data_err[i - flowstart] = [spline(tauFbyScale_arr[0][i]) for spline in g2_err_ints]

    if args.calc_cont:
        for i in range(nflow):
            fitparams, fitparams_err = bootstr.bootstr_from_gauss(fit_sample, data=data[i], data_std_dev=data_err[i], numb_samples=1000,
                                                                  sample_size=1, return_sample=False,
                                                                  args=[1 / numpy.asarray(args.Nts) ** 2, data_err[i], [1, 0]],
                                                                  parallelize=False)
            cont[i, 0] = fitparams[1]
            cont[i, 1] = fitparams_err[1]
            slope[i, 0] = fitparams[0]
            slope[i, 1] = fitparams_err[0]
        print("saving", output_file)
        numpy.savetxt(output_file, numpy.column_stack((tauFbyScale_arr[0][flowstart:flowend], cont, slope)),
                      header="tau_F/r_0^2 g2_cont    err     slope     err")
    else:
        try:
            _, b, berr, m, merr = numpy.loadtxt(output_file, unpack=True)
        except OSError:
            print("Error: could not find ", output_file)
            exit(1)
        for i in range(nflow):
            cont[i, 0] = b[i]
            cont[i, 1] = berr[i]
            slope[i, 0] = m[i]
            slope[i, 1] = merr[i]

    print(cont[-1, 1]/cont[-1, 0])

    return flowstart, flowend, data, data_err, cont, slope


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--calc_cont', help='calc continuum extrapolation and save to file instead of reading it from the file', action="store_true")
    parser.add_argument('--ref_scale', choices=["r0", "t0"], required=True)
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

    # # plot both cont extr
    # if args.ref_scale == "r0":
    #     otherlabel = "t0"  # yes this is the other one!
    # elif args.ref_scale == "t0":
    #     otherlabel = "r0"
    # otherdatafile = outputpathdata + "g2_" + otherlabel + "_cont_extr.txt"
    # tauFbytau0, b, berr, m, merr = numpy.loadtxt(otherdatafile, unpack=True)
    #
    # fig4, ax4, _ = lpd.create_figure(xlabel=r'$\frac{\tau_F}{'+plotlabel+r'}$', xlabelpos=(1.1, 0.05), ylabelpos=(0.1, 0.9), figsize=((16/9)*(3+3/8), 3+3/8 - 1/2.54))  # r'$\tau_F T^2\big|_{T=1.5T_c}$',
    # ax4.set_title(r'$ g^2_\mathrm{gf}\equiv \frac{128\pi^2 \tau_F^2 \langle E \rangle}{24}$', y=0.6)
    # ax4.errorbar(tauFbytau0_arr[0], cont[:, 0], cont[:, 1], fmt='x-', lw=0.5, markersize=2, mew=0.3, capsize=1, label=r'$'+plotlabel+r'$', zorder=-1)
    # ax4.errorbar(tauFbytau0_arr[0], cont[:, 0], cont[:, 1], fmt='x-', lw=0.5, markersize=2, mew=0.3, capsize=1, label='cont', zorder=-1)
    # ax4.legend(**lpd.legendstyle)
    # fig4.savefig(outputpath+"g2_"+label+".pdf")
