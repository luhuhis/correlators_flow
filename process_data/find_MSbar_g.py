#!/usr/local/bin/python3.7 -u
import lib_process_data as lpd
import numpy
import scipy.interpolate
import scipy.optimize
import extrapolate_coupling


def Eq27(alpha_MS, alpha_gf, L, nf=0):
    beta0=11/(4*numpy.pi)-1/(6*numpy.pi)*nf
    beta1=51/(8*numpy.pi**2)-19/(24*numpy.pi**2)*nf
    k1 = 1.098 + 0.008 * nf
    k2 = -0.982 - 0.070 * nf + 0.002 * nf**2
    return (alpha_MS + alpha_MS**2 * (k1 + beta0*L) + alpha_MS**3 *(k2+(beta1+2*k1*beta0)*L+ beta0**2 * L**2) ) - alpha_gf


def running_coupling(mu, nf, r0):
    C3 = 1.202057

    beta0 = (1 / 4) * (11 - (2 / 3) * nf)
    beta1 = (1 / 16) * (102 - (38 / 3) * nf)
    beta2 = (1 / 64) * ((2857 / 2) - (5033 / 18) * nf + (325 / 54) * nf ** 2)
    beta3 = (1 / 256) * (
                (149753 / 6) + 3564 * C3 + (-(1078361 / 162) - (6508 / 27) * C3) * nf + ((50065 / 162) + (6472 / 81) * C3) * nf ** 2 + (1093 / 729) * nf ** 3)

    b1 = beta1 / beta0
    b2 = beta2 / beta0
    b3 = beta3 / beta0

    # Lambda MS taken from S. Capitani et al. [ALPHA Collaboration], Nucl. Phys. B544, 669 (1999) [hep-lat/9810063]
    lambd = 0.602 / r0

    L = numpy.log(mu ** 2 / lambd ** 2)
    lnL = numpy.log(L)

    b0L = 1 / (beta0 * L)

    a_s = (b0L - (b1 * lnL) / ((beta0 * L) ** 2) + b0L ** 3 * (b1 ** 2 * (lnL ** 2 - lnL - 1) + b2) + b0L ** 4 * (
                b1 ** 3 * (-lnL ** 3 + (5 / 2) * lnL ** 2 + 2 * lnL - (1 / 2)) - 3 * b1 * b2 * lnL + (b3 / 2))) * numpy.pi

    return a_s


def TbyTc(Nt, beta):
    r0Tc = 0.7457
    return 1/(r0Tc * Nt) * extrapolate_coupling.r0bya(beta)


def main():

    # load data
    inputpathdata = "/work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/coupling/"
    outputpathplot = "/work/home/altenkort/work/correlators_flow/plots/quenched_1.50Tc_zeuthenFlow/coupling/"
    datafile = inputpathdata + "g2_cont_extr.txt"
    try:
        tauFByScale, g2_flow_data, g2_flow_data_err, _, _ = numpy.loadtxt(datafile, unpack=True)
    except OSError:
        print("Error: could not find ", datafile)
        exit(1)
    # nflow = len(tauFT2)

    # interpolate continuum extrapolated g2 flow to any tauFT2
    order = 3
    g2_flow = scipy.interpolate.InterpolatedUnivariateSpline(tauFByScale, g2_flow_data, k=order)
    g2_flow_err = scipy.interpolate.UnivariateSpline(tauFByScale, g2_flow_data_err, k=order, s=0.00001)  # smooth the errors a tiny bit

    xpoints = numpy.linspace(0, 0.0175, 10000)

    # plot interpolated g2 flow
    fig2, ax2, _ = lpd.create_figure(xlabel=r'$\tau_F T^2\big|_{T=1.5T_c}$', xlabelpos=(1.1, 0.05), ylabelpos=(0.1, 0.9), figsize=((16/9)*(3+3/8), 3+3/8 - 1/2.54))
    ax2.set_title(r'$ g^2_\mathrm{gf}$', y=0.6)
    ax2.errorbar(xpoints, g2_flow(xpoints), fmt='-', lw=0.5, markersize=0, mew=0.3, capsize=1, label='cont', zorder=-1)
    ax2.legend(**lpd.legendstyle, title=r'$N_\tau$')
    fig2.savefig(outputpathplot+"g2_cont_int.pdf")

    # get interesting OmegaByT rang
    PhiUV = numpy.loadtxt("/work/data/htshu/ee_spf/PHIUV_a.dat")
    PhiUV = PhiUV[:, 0:2]
    OmegaByT_ = PhiUV[0:, 0]
    # PhiuvByT3 = scipy.interpolate.InterpolatedUnivariateSpline(PhiUV[:, 0], PhiUV[:, 1], k=order)
    # MaxOmegaByT = PhiUV[-1][0]

    # plot MSbar coupling
    fig3, ax3, _ = lpd.create_figure(xlabel=r'$\omega/T_c$', xlabelpos=(1.1, 0.05), ylabelpos=(0.1, 0.9), figsize=((16/9)*(3+3/8), 3+3/8 - 1/2.54))
    ax3.set_title(r'$ g^2_{\overline{MS}}$', y=0.6)
    ax3.set_xscale('log')
    ax3.axvline(x=numpy.pi, **lpd.verticallinestyle)

    Nts_finite_temp = (36, 30, 24, 20)
    betas = (7.544, 7.394, 7.192, 7.035, 6.87361)

    # Tcsqrtt0 = 0.2490
    # sqrtt0Byr0 = 0.333
    r0Tc = 0.7457
    nf = 0
    min_scale_byT = numpy.pi
    for i, Nt in enumerate(Nts_finite_temp):
        beta = betas[i]
        # compute MSbar coupling using four-loop running
        OmegaByTc = []
        muByTc = []
        g2_MSbar = []
        thisTbyTc = TbyTc(Nt, beta)
        for OmegaByT in OmegaByT_:
            OmegaByTc.append(OmegaByT*thisTbyTc)
            muByTc.append(numpy.fmax(min_scale_byT*thisTbyTc, OmegaByT*thisTbyTc))
            g2_MSbar.append(4*numpy.pi*running_coupling(muByTc[-1], nf, r0Tc))
        ax3.errorbar(OmegaByTc, g2_MSbar, fmt='x', lw=0.5, markersize=1, mew=0.3, capsize=1, label=str(Nt)+", "+str(beta), zorder=-1)

    # converted MSbar coupling from gradient flow
    L_ = (0.347,) # -numpy.log(2), 0, numpy.log(2))
    for L in L_:
        g2_MS = []
        for OmegaByT in OmegaByT_:
            muByT = numpy.fmax(numpy.pi, OmegaByT)

            # tauFByScale = numpy.exp(L) / (8 * muByScale**2)
            #TODO: what to do here?
            tauFT2 = numpy.exp(L) / (8 * muByT**2)
            this_g2_flow = g2_flow(tauFT2)
            g2_MS.append(scipy.optimize.fsolve(Eq27, 1, args=(this_g2_flow/(4*numpy.pi), L))[0] * 4 * numpy.pi)
        ax3.errorbar(tauFT2, g2_MS, fmt='x', lw=0.5, markersize=1, mew=0.3, capsize=1, label=str(Nt) + ", " + str(beta), zorder=-1)

    ax3.legend(**lpd.legendstyle, title=r'$N_\tau, \beta$')
    fig3.savefig(outputpathplot+"g2_MSbar.pdf")


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
