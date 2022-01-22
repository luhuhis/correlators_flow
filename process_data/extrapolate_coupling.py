#!/usr/local/bin/python3.7 -u
import lib_process_data as lpd
import numpy
from latqcdtools import bootstr
import argparse
import scipy.optimize
import scipy.interpolate


def extrapolation_ansatz(x, m, b):
    return m * x + b


def fit_sample(ydata, xdata, edata, start_params=None):
    fitparams, _ = scipy.optimize.curve_fit(extrapolation_ansatz, xdata, ydata, p0=start_params, sigma=edata)
    return fitparams


# TODO add source for these parameters
def sqrtt0bya(beta):
    b0=11./(4*numpy.pi)**2
    b1=102/(4*numpy.pi)**4
    c1=-9.945
    c2=24.191
    c3=-5.334
    c4=1.452
    return numpy.exp((beta/(12*b0)+b1/(2.*b0**2)*numpy.log(6*b0/beta))*(1+c1/beta+c2/beta**2)/(1+c3/beta+c4/beta**2))


# TODO add source for these parameters
def r0bya(beta):
    b0=11./(4*numpy.pi)**2
    b1=102/(4*numpy.pi)**4
    c1=-8.9664
    c2=19.21
    c3=-5.25217
    c4=0.606828
    return numpy.exp( (beta/(12*b0)+b1/(2.*b0**2)*numpy.log(6*b0/beta))*(1+c1/beta+c2/beta**2)/(1+c3/beta+c4/beta**2))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--calc_cont', help='calc continuum extrapolation and save to file instead of reading it from the file', action="store_true")
    parser.add_argument('--ref_scale', choices=["r0", "t0"], required=True)
    args = parser.parse_args()

    if args.ref_scale == "r0":
        plotlabel = r'r_0^2'
        label = "r0"
        scalefunc = r0bya
    elif args.ref_scale == "t0":
        plotlabel = r't_0'
        label = "t0"
        scalefunc = sqrtt0bya

    # input/output paths
    inputpath = "/work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/coupling/"
    outputpath = "/work/home/altenkort/work/correlators_flow/plots/quenched_1.50Tc_zeuthenFlow/coupling/"
    outputpathdata = "/work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/coupling/"
    lpd.create_folder(outputpath)
    datafile = outputpathdata + "g2_"+label+"_cont_extr.txt"

    # info about the lattices
    confdict = {"64": "s064t64_b0687361", "80": "s080t80_b0703500", "96": "s096t96_b0719200", "120": "s096t120_b0739400", "144": "s096t144_b0754400"}
    Nts = numpy.asarray((144, 120, 96, 80))  # , 64))
    # Nts_finite_temp = (36, 30, 24, 20)
    betas = (7.544, 7.394, 7.192, 7.035, 6.87361)

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
    for i, Nt in enumerate(Nts):
        tauFbya2, tauF2E, tauF2E_err = numpy.loadtxt(inputpath + "flow_t2E_" + confdict[str(Nt)] + ".dat", unpack=True)

        # flow times
        tauFbyScale_arr.append(tauFbya2/scalefunc(betas[i])**2)
        # tauFT2_arr.append(tauFbya2 / Nts_finite_temp[i] ** 2)

        # coupling
        g2_arr.append(prefactor * tauF2E)
        g2_err_arr.append(prefactor * tauF2E_err)
        g2_ints.append(scipy.interpolate.InterpolatedUnivariateSpline(tauFbyScale_arr[-1], g2_arr[-1], k=order))
        g2_err_ints.append(scipy.interpolate.InterpolatedUnivariateSpline(tauFbyScale_arr[-1], g2_err_arr[-1], k=order))

    nNts = len(Nts)
    nflow = len(tauFbyScale_arr[0])

    data = numpy.empty((nflow, nNts))
    data_err = numpy.empty((nflow, nNts))
    cont = numpy.empty((nflow, 2))
    slope = numpy.empty((nflow, 2))

    for i in range(nflow):
        data[i] = [spline(tauFbyScale_arr[0][i]) for spline in g2_ints]
        data_err[i] = [spline(tauFbyScale_arr[0][i]) for spline in g2_err_ints]

    if args.calc_cont:
        for i in range(nflow):
            fitparams, fitparams_err = bootstr.bootstr_from_gauss(fit_sample, data=data[i], data_std_dev=data_err[i], numb_samples=1000,
                                                                  sample_size=1, return_sample=False, args=[1/Nts**2, data_err[i], [1, 0]], parallelize=False)
            cont[i, 0] = fitparams[1]
            cont[i, 1] = fitparams_err[1]
            slope[i, 0] = fitparams[0]
            slope[i, 1] = fitparams_err[0]
        numpy.savetxt(datafile, numpy.column_stack((tauFbyScale_arr[0], cont, slope)), header="tfT2 |_T=1.5Tc  g2_cont    err     slope     err")
    else:
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
    fig2, ax2, _ = lpd.create_figure(xlabel=r'$\frac{\tau_F}{'+plotlabel+r'}$', xlabelpos=(1.03, 0.05), ylabelpos=(0.1, 0.9), figsize=((16/9)*(3+3/8), 3+3/8 - 1/2.54))  # r'$\tau_F T^2\big|_{T=1.5T_c}$',
    ax2.set_title(r'$ g^2_\mathrm{gf}\equiv \frac{128\pi^2 \tau_F^2 \langle E \rangle}{24}$', y=0.6)
    ax2.errorbar(tauFbyScale_arr[0], cont[:, 0], cont[:, 1], fmt='x-', lw=0.5, markersize=2, mew=0.3, capsize=1, label='cont', zorder=-1)
    for i, Nt in enumerate(Nts):
        ax2.errorbar(tauFbyScale_arr[i], g2_arr[i], g2_err_arr[i], fmt='x-', lw=0.5, markersize=2, mew=0.3, capsize=1, label=str(Nt), zorder=-i-10)
    ax2.legend(**lpd.legendstyle, title=r'$N_\tau$')
    fig2.savefig(outputpath+"g2_"+label+".pdf")

    # plot continuum extrapolation at different flow times
    fig3, ax3, _ = lpd.create_figure(xlabel=r'$N_\tau^{-2}$', xlabelpos=(1.03, 0.05), ylabelpos=(0.1, 0.9), figsize=(1.5*(3+3/8), 3+3/8 - 1/2.54), UseTex=False)
    ax3.set_title(r'$ g^2_\mathrm{gf}\equiv \frac{128\pi^2 \tau_F^2 \langle E \rangle}{24}$ ', y=0.6)
    min_idx = 5
    xpoints = numpy.linspace(0, 1.05*1/Nts[-1]**2, 10)
    for i in range(nflow, -1, -1):
        if i % 10 == 0 and i > min_idx:
            ax3.errorbar(numpy.insert(1/Nts**2, 0, 0), numpy.insert(data[i], 0, cont[i, 0]), numpy.insert(data_err[i], 0, cont[i, 1]),
                         color=lpd.get_color(range(nflow), i, min_idx, nflow-1), fmt='x', lw=0.5, markersize=2, mew=0.3, capsize=1,
                         label='{0:.6f}'.format(tauFbyScale_arr[0][i]), zorder=-i)
            ax3.errorbar(xpoints, extrapolation_ansatz(xpoints, slope[i, 0], cont[i, 0]), fmt='--', lw=0.5, alpha=0.5, color=lpd.get_color(range(nflow), i, min_idx, nflow-1))
    ax3.legend(**lpd.legendstyle, title=r'$\frac{\tau_F}{'+plotlabel+r'}$')
    fig3.savefig(outputpath+"g2_"+label+"_cont_extr.pdf")

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


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
