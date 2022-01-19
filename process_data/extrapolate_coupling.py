#!/usr/local/bin/python3.7 -u
import lib_process_data as lpd
import numpy
from latqcdtools import bootstr
import argparse
import scipy.optimize


def extrapolation_ansatz(x, m, b):
    return m * x + b


def fit_sample(ydata, xdata, edata, start_params=None):
    fitparams, _ = scipy.optimize.curve_fit(extrapolation_ansatz, xdata, ydata, p0=start_params, sigma=edata)
    return fitparams


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--load_cont', help='load continuum extrapolation instead of calculating it again', action="store_true")
    args = parser.parse_args()

    inputpath = "/work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/coupling/"
    outputpath = "/work/home/altenkort/work/correlators_flow/plots/quenched_1.50Tc_zeuthenFlow/coupling/"
    outputpathdata = "/work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/coupling/"
    lpd.create_folder(outputpath)
    confdict = {"64": "s064t64_b0687361", "80": "s080t80_b0703500", "96": "s096t96_b0719200", "120": "s096t120_b0739400", "144": "s096t144_b0754400"}
    Nts = numpy.asarray((144, 120, 96, 80))  #, 64))
    Nts_finite_temp = (36, 30, 24, 20)

    tauFT2_arr = []
    g2_arr = []
    g2_err_arr = []
    prefactor = 128 * numpy.pi ** 2 / 24
    for i, Nt in enumerate(Nts):
        tauFbya2, tauF2E, tauF2E_err = numpy.loadtxt(inputpath + "flow_t2E_" + confdict[str(Nt)] + ".dat", unpack=True)
        tauFT2_arr.append(tauFbya2 / Nts_finite_temp[i] ** 2)  # use the correct scale
        g2_arr.append(prefactor * tauF2E)
        g2_err_arr.append(prefactor * tauF2E_err)

    nNts = len(Nts)
    nflow = len(tauFT2_arr[0])

    data = numpy.empty((nflow, nNts))
    data_err = numpy.empty((nflow, nNts))
    cont = numpy.empty((nflow, 2))
    slope = numpy.empty((nflow, 2))

    datafile = outputpathdata + "g2_cont_extr.txt"

    for i in range(nflow):
        data[i] = [entry[i] for entry in g2_arr]
        data_err[i] = [entry[i] for entry in g2_err_arr]

    if not args.load_cont:
        for i in range(nflow):
            fitparams, fitparams_err = bootstr.bootstr_from_gauss(fit_sample, data=data[i], data_std_dev=data_err[i], numb_samples=100,
                                                                  sample_size=1, return_sample=False, args=[1/Nts**2, data_err[i]], parallelize=False)
            cont[i, 0] = fitparams[1]
            cont[i, 1] = fitparams_err[1]
            slope[i, 0] = fitparams[0]
            slope[i, 1] = fitparams_err[0]
        numpy.savetxt(datafile, numpy.column_stack((tauFT2_arr[0], cont, slope)), header="tfT2 |_T=1.5Tc  g2_cont    err     slope     err")
    elif args.load_cont:
        try:
            _, b, berr, m, merr = numpy.loadtxt(datafile, unpack=True)
        except OSError:
            print("Error: could not find ", datafile)
            exit(1)
        for i in range(nflow):
            cont[i, 0] = b[i]
            cont[i, 1] = berr[i]
            slope[i, 0] = m[i]
            slope[i, 1] = merr[i]

    # plot coupling for all lattices
    fig2, ax2, _ = lpd.create_figure(xlabel=r'$\tau_F T^2\big|_{T=1.5T_c}$', xlabelpos=(1.1, 0.05), ylabelpos=(0.1, 0.9), figsize=((16/9)*(3+3/8), 3+3/8 - 1/2.54))
    ax2.set_title(r'$ g^2_\mathrm{gf}$', y=0.6)
    ax2.errorbar(tauFT2_arr[0], cont[:, 0], cont[:, 1], fmt='x-', lw=0.5, markersize=2, mew=0.3, capsize=1, label='cont', zorder=-1)
    for i, Nt in enumerate(Nts):
        ax2.errorbar(tauFT2_arr[i], g2_arr[i], g2_err_arr[i], fmt='x-', lw=0.5, markersize=2, mew=0.3, capsize=1, label=str(Nt), zorder=-Nt)
    ax2.legend(**lpd.legendstyle, title=r'$N_\tau$')
    fig2.savefig(outputpath+"g2.pdf")

    # plot continuum extrapolation at different flow times
    fig3, ax3, _ = lpd.create_figure(xlabel=r'$N_\tau^{-2}$', xlabelpos=(1.1, 0.05), ylabelpos=(0.1, 0.9), figsize=(1.3*(3+3/8), 3+3/8 - 1/2.54), UseTex=False)
    ax3.set_title(r'$ g^2_\mathrm{gf}$ ', y=0.6)
    min_idx = 5
    xpoints = numpy.linspace(0, 1.05*1/Nts[-1]**2, 10)
    for i in range(nflow, -1, -1):
        if i % 10 == 0 and i > min_idx:
            print(data[i])
            ax3.errorbar(numpy.insert(1/Nts**2, 0, 0), numpy.insert(data[i], 0, cont[i, 0]), numpy.insert(data_err[i], 0, cont[i, 1]),
                         color=lpd.get_color(range(nflow), i, min_idx, nflow-1), fmt='x', lw=0.5, markersize=2, mew=0.3, capsize=1,
                         label='{0:.6f}'.format(tauFT2_arr[0][i]), zorder=-i)
            ax3.errorbar(xpoints, extrapolation_ansatz(xpoints, slope[i, 0], cont[i, 0]), fmt='--', lw=0.5, alpha=0.5, color=lpd.get_color(range(nflow), i, min_idx, nflow-1))
    ax3.legend(**lpd.legendstyle, title=r'$\tau_F T^2\big|_{T=1.5T_c} $')
    fig3.savefig(outputpath+"g2_cont_extr.pdf")


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
