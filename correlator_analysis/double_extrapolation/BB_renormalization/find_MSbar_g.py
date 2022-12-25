#!/usr/local/bin/python3.7 -u
import lib_process_data as lpd
import numpy
import scipy.interpolate
import scipy.optimize
import argparse


def Eq27(alpha_MS, alpha_gf, L, nf=0):
    beta0 = 11/(4*numpy.pi)-1/(6*numpy.pi)*nf
    beta1 = 51/(8*numpy.pi**2)-19/(24*numpy.pi**2)*nf
    k1 = 1.098 + 0.008 * nf
    k2 = -0.982 - 0.070 * nf + 0.002 * nf**2
    return (alpha_MS + alpha_MS**2 * (k1 + beta0*L) + alpha_MS**3 * (k2+(beta1+2*k1*beta0)*L + beta0**2 * L**2)) - alpha_gf


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--ref_scale', choices=["r0", "t0"], required=True, help="scale in which g2_cont_extr is provided")
    parser.add_argument('--TbyTc', type=float, help='needed to convert to the right scale. for continuum data just use the value of the finest lattice.')
    parser.add_argument('--r0Tc', default=0.7457, help="default r_0 T_c from https://arxiv.org/pdf/1503.05652.pdf", type=float)
    parser.add_argument('--g2_pert', default="Nf0_g2_1.5Tc.dat", help="path to perturbative g^2 file. column format: omega/T, g^2")
    args = parser.parse_args()

    label = None
    if args.ref_scale == "r0":
        label = "r0"
    elif args.ref_scale == "t0":
        label = "t0"

    inputpathdata = "/work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/coupling/"
    outputpathplot = "/work/home/altenkort/work/correlators_flow/plots/quenched_1.50Tc_zeuthenFlow/coupling/"
    datafile = inputpathdata + "g2_"+label+"_cont_extr.txt"

    # load data
    tauFByScale_ = None
    g2_flow_data = None
    try:
        tauFByScale_, g2_flow_data, g2_flow_data_err, _, _ = numpy.loadtxt(datafile, unpack=True)
    except OSError:
        print("Error: could not find ", datafile)
        exit(1)

    # interpolate continuum extrapolated g2 flow to any tauFT2
    order = 3
    g2_flow = scipy.interpolate.InterpolatedUnivariateSpline(tauFByScale_, g2_flow_data, k=order, ext=2)
    # g2_flow_err = scipy.interpolate.UnivariateSpline(tauFByScale_, g2_flow_data_err, k=order, s=0.00001, ext=2)  # smooth the errors a tiny bit

    xpoints = numpy.linspace(tauFByScale_[0], tauFByScale_[-1], 10000)

    # plot interpolated g2 flow
    fig2, ax2, _ = lpd.create_figure(xlabel=r'$\tau_F T^2\big|_{T=1.5T_c}$', xlabelpos=(0.9, 0.1), ylabelpos=(0.1, 0.9), figsize=((16/9)*(3+3/8), 3+3/8 - 1/2.54))
    ax2.set_title(r'$ g^2_\mathrm{gf}$', y=0.6)
    ax2.errorbar(xpoints, g2_flow(xpoints), fmt='-', lw=0.5, markersize=0, mew=0.3, capsize=1, label='cont', zorder=-1)
    ax2.legend(**lpd.legendstyle, title=r'$N_\tau$')
    fig2.savefig(outputpathplot+"g2_cont_int.pdf")

    # plot MSbar coupling
    fig3, ax3, _ = lpd.create_figure(xlabel=r'$\frac{\omega}{T}\big|_{T=1.5T_c}$', xlims=(1, 1100), ylims=(-0.1, 4.5), xlabelpos=(0.9, 0.1), ylabelpos=(0.1, 0.9),
                                     figsize=((16/9)*(3+3/8), 3+3/8 - 1/2.54))
    ax3.set_title(r'', y=0.6)
    ax3.set_xscale('log')
    ax3.axvline(x=numpy.pi, **lpd.verticallinestyle)

    # load pert coupling and plot it
    OmegaByT_MS, g2_MS = numpy.loadtxt(args.g2_pert, unpack=True)
    ax3.errorbar(OmegaByT_MS, g2_MS, fmt='-', lw=0.5, markersize=1, mew=0.3, capsize=1, zorder=-1,
                 label=r'$ g^2_{\overline{MS}}({\mu}), \mathrm{(4-loop)}$')

    # converted MSbar coupling from gradient flow
    L_ = (-numpy.log(2), 0, 0.347, numpy.log(2), 1)  #
    for L in L_:
        validOmegaByT_ = []
        g2_converted = []
        for OmegaByT in OmegaByT_MS:
            tauFT2 = numpy.exp(L) / (8 * OmegaByT**2)  # match the scale to the gradient flow
            tauFByr02 = tauFT2/args.TbyTc**2/args.r0Tc**2  # TODO abstract this to other scales than r_0 !
            try:
                this_g2_flow = g2_flow(tauFByr02)
                validOmegaByT_.append(OmegaByT)
                sol = (scipy.optimize.fsolve(Eq27, numpy.asarray(0), args=(this_g2_flow / (4 * numpy.pi), L)) * 4 * numpy.pi)
                g2_converted.append(sol[0])
            except ValueError:
                pass
        ax3.errorbar(validOmegaByT_, g2_converted, fmt='-', lw=0.5, markersize=1, mew=0.3, capsize=1, zorder=-1,
                     label="$L="+'{0:.3f}'.format(L)+"$")
    ax3.set_title("$\\mu\\equiv  \\omega, \\quad g^2_{\\overline{MS}}(g^2_\\mathrm{gf}(\\tau_F = \\frac{e^{L}}{8\\mu^2}))$")
    ax3.legend(**lpd.chmap(lpd.legendstyle, labelspacing=1))
    fig3.savefig(outputpathplot+"g2_MSbar.pdf")


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
