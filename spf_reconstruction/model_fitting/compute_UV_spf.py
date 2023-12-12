#!/usr/bin/env python3
import rundec
import numpy as np
import argparse
import lib_process_data as lpd
import BB_NLO_SPF_term as BB
from nptyping import NDArray, Float64
from typing import Literal as Shape
from scipy import interpolate, integrate


def smooth_max(a, b):
    return np.sqrt(a ** 2 + b ** 2)


def get_minscale(min_scale, T_in_GeV, Nc, Nf):
    effective_theory_thermal_scale = 4 * np.pi * T_in_GeV * np.exp(-np.euler_gamma) * np.exp(
        (-Nc + 4 * Nf * np.log(4)) / (22 * Nc - 4 * Nf))

    # check min_scale
    if min_scale == "piT":
        min_scale = np.pi * T_in_GeV
    elif min_scale == "2piT":
        min_scale = 2 * np.pi * T_in_GeV
    elif min_scale == "eff":
        min_scale = effective_theory_thermal_scale
    elif min_scale == "mu_IR_NLO":
        min_scale = 4 * np.pi * np.exp(1 - np.euler_gamma) * T_in_GeV
    else:
        min_scale = float(min_scale)

    print("min_scale = ", lpd.format_float(min_scale / T_in_GeV), "T")

    return min_scale


def get_omega_prefactor(omega_prefactor, Nc, Nf, T_in_GeV, min_scale):
    mu_opt_omega_term = 2 * np.exp(
        ((24 * np.pi ** 2 - 149) * Nc + 20 * Nf) / (6 * (11 * Nc - 2 * Nf)))  # see eq. 4.17 of arXiv:1006.0867
    omega_exponent = 1.0
    if omega_prefactor == "opt":
        omega_prefactor = mu_opt_omega_term  # ~14.743 for Nf=3, ~7.6 for Nf=0
    elif omega_prefactor == "optBB":
        gamma_0 = 3 / (8 * np.pi ** 2)  # color-magnetic anomalous dimension
        b_0 = 11 / (16 * np.pi ** 2)
        omega_exponent = (1 - gamma_0 / b_0)
        omega_prefactor = min_scale ** (gamma_0 / b_0)
    elif omega_prefactor == "optBBpiT":
        gamma_0 = 3 / (8 * np.pi ** 2)  # color-magnetic anomalous dimension
        b_0 = 11 / (16 * np.pi ** 2)
        omega_exponent = (1 - gamma_0 / b_0)
        omega_prefactor = (np.pi * T_in_GeV) ** (gamma_0 / b_0)
    else:
        omega_prefactor = float(omega_prefactor)

    print("omega_prefactor = ", lpd.format_float(omega_prefactor))

    return omega_prefactor, omega_exponent


def print_min_scale(scale, min_scale, first, OmegaByT):
    if scale > min_scale and first:
        print("scale > min_scale at OmegaByT=", lpd.format_float(OmegaByT))
        first = False
    return first


def get_maxfunc(max_type):
    if max_type == "smooth":
        maxfunc = smooth_max
    elif max_type == "hard":
        maxfunc = np.fmax
    return maxfunc


def get_lambdamsbar(Nf):
    if Nf == 0:
        Lambda_MSbar = 0.253  # (JHEP11(2017)206)
    elif Nf == 3:
        Lambda_MSbar = 0.339  # 2111.09849 (p.11)
    else:
        print("ERROR: I don't know any LambdaMSbar for this Nf")
        exit(1)
    return Lambda_MSbar


def get_g2_and_alphas_and_lospf(crd, Lambda_MSbar, mu, Nf, Nloop, C_F, OmegaByT):
    Alphas = crd.AlphasLam(Lambda_MSbar, mu, Nf, Nloop)
    g2 = 4. * np.pi * Alphas
    lo_spf = g2 * C_F * OmegaByT ** 3 / 6. / np.pi
    return g2, Alphas, lo_spf


@lpd.typed_frozen_data
class SpfParams:
    omega_prefactor_input: str
    OmegaByT_values: NDArray[Shape["*"], Float64]
    max_type: str
    T_in_GeV: float
    min_scale: float
    Lambda_MSbar: float
    Nf: int
    Nloop: int
    omega_prefactor: float
    omega_exponent: float
    Nc: int
    r20: float
    r21: float
    C_F: float
    b_0: float
    gamma_0: float
    Npoints: int
    order: str
    corr: str
    mu_IR_by_T: float


def get_cBsq(params: SpfParams):
    # Given input parameters This function returns
    # c_B^2(mu,mubarir) =  exp(int^{\mu^2}_{\mubarir^2}\gamma_0 \gsqms(\mubar)\frac{d\mubar^2}{\mubar^2})
    # Note that in the output function for cBsq, the input variable mu is expected to have units GeV.

    mu_IR_in_GeV = params.mu_IR_by_T * params.T_in_GeV

    # Here we just want array that covers the whole range of possible scales in GeV up to the deep UV (we use at most 2*omega and integrate up to 1000T, so 2500 should be enough),
    # but not unnecessary small values as crd then complains.
    minimum_mu = mu_IR_in_GeV if mu_IR_in_GeV <= params.min_scale else params.min_scale
    mu_in_GeV_arr = np.asarray([*np.linspace(0.9*minimum_mu, 2500*params.T_in_GeV, params.Npoints)])

    crd = rundec.CRunDec()

    def g2(mu_in_GeV):
        # run g^2 from the MSBAR coupling we found by measuring the flow coupling and converting it at mu_flow/T=mu_MSBAR/T=4
        mu0 = 4.000425807761551766e+00 * params.T_in_GeV  # TODO remove hard coding of these values and instead read them from the corresponding file
        alphas_0 = 2.742374650650232670e+00 / (4 * np.pi)
        return 4. * np.pi * crd.AlphasExact(alphas_0, mu0, mu_in_GeV, params.Nf, params.Nloop)
    g2_vec = np.vectorize(g2)

    integrand_arr = 2 * params.gamma_0 * g2_vec(mu_in_GeV_arr) / mu_in_GeV_arr
    integrand_spline = interpolate.InterpolatedUnivariateSpline(mu_in_GeV_arr, integrand_arr, k=3, ext=2)

    cBsq_arr = np.asarray(  # TODO make nproc an input
        lpd.parallel_function_eval(calc_cBsq, mu_in_GeV_arr, 100, integrand_spline, mu_IR_in_GeV))

    cBsq = interpolate.InterpolatedUnivariateSpline(mu_in_GeV_arr, cBsq_arr, k=3, ext=2)

    return cBsq


def calc_cBsq(target_scale, integrand_spline, ir_scale):
    # Both start_scale and ir_scale are expected to be in GeV
    try:
        integral = integrate.quad(integrand_spline, ir_scale, target_scale)[0]
    except ValueError:
        print("error at", target_scale, ir_scale)
        exit(1)
    return np.exp(integral)


def inner_loop(index, params: SpfParams, cBsq):
    crd = rundec.CRunDec()
    OmegaByT = params.OmegaByT_values[index]
    maxfunc = get_maxfunc(params.max_type)

    scale = params.omega_prefactor * (OmegaByT * params.T_in_GeV) ** params.omega_exponent
    mu = maxfunc(params.min_scale, scale)

    if params.corr == "BB":
        # Note that "omega_prefactor" is actually not used here.
        g2, Alphas, PhiUVByT3 = get_g2_and_alphas_and_lospf(crd, params.Lambda_MSbar, mu, params.Nf, params.Nloop,
                                                            params.C_F, OmegaByT)

        if params.order == "NLO":
            # === NLO ===
            # The following equations follow C.23 of https://arxiv.org/pdf/2204.14075.pdf
            first_paren = 5./3.*np.log(mu**2/(2*OmegaByT*params.T_in_GeV)**2) + 134./9. - 2*np.pi**2/3.  # Guy corrected the last term: -8*np.pi**2/3 to -2*np.pi**2/3.
            second_paren = 2./3*np.log(mu**2/(2*OmegaByT*params.T_in_GeV)**2) + 26./9.
            mu_dep_part = PhiUVByT3*(1 + g2/(4*np.pi)**2*(params.Nc*first_paren - params.Nf*second_paren))
            mu_free_part = (g2**2*params.C_F)/(12*np.pi**3) * BB.get_mu_free(OmegaByT, params.Nc, params.Nf)  # from line 3 to the end of eq.(C.23)
            # mu_free_part = 0
            PhiUVByT3 = cBsq(mu) * (mu_dep_part + mu_free_part)

    elif params.corr == "EE":
        g2, Alphas, PhiUVByT3 = get_g2_and_alphas_and_lospf(crd, params.Lambda_MSbar, mu, params.Nf, params.Nloop,
                                                            params.C_F, OmegaByT)
        if params.order == "NLO":
            l_ = np.log((scale / params.T_in_GeV) ** 2 / OmegaByT ** 2)
            PhiUVByT3 = PhiUVByT3 * (1 + (params.r20 + params.r21 * l_) * Alphas / np.pi)

    else:
        print("Error: unknown corr", params.corr)
        exit(1)

    return OmegaByT, PhiUVByT3, g2


def get_parameters(Nf, max_type, min_scale_str, T_in_GeV, omega_prefactor_input, Npoints, Nloop, order, corr, mu_IR_by_T_str):
    if order not in ["LO", "NLO"]:
        print("Error: unknown UV spf order", order)
        exit(1)

    # set some parameters
    Lambda_MSbar = get_lambdamsbar(Nf)
    Nc = 3
    min_scale = get_minscale(min_scale_str, T_in_GeV, Nc, Nf)
    omega_prefactor, omega_exponent = get_omega_prefactor(omega_prefactor_input, Nc, Nf, T_in_GeV, min_scale)
    OmegaByT_values = np.logspace(-5, 3, Npoints, base=10)

    mu_IR_by_T = 4 * np.pi * np.exp(1-np.euler_gamma) if mu_IR_by_T_str == "NLO" else (2 * np.pi if mu_IR_by_T_str == "LO" else np.nan)
    print("mu_IR_by_T=", mu_IR_by_T)

    r20 = Nc * (149. / 36. - 11. * np.log(2.) / 6. - 2 * np.pi ** 2 / 3.) - Nf * (5. / 9. - np.log(2.) / 3.)
    r21 = (11. * Nc - 2. * Nf) / 12.
    C_F = (Nc ** 2 - 1) / 2 / Nc
    b_0 = 11 / (16 * np.pi ** 2)
    gamma_0 = 3 / (8 * np.pi ** 2)

    additional_params = SpfParams(omega_prefactor_input, OmegaByT_values, max_type, T_in_GeV, min_scale, Lambda_MSbar,
                                  Nf, Nloop, omega_prefactor, omega_exponent, Nc, r20, r21, C_F, b_0, gamma_0, Npoints,
                                  order, corr, mu_IR_by_T)

    return additional_params


def get_spf(Nf: int, max_type: str, min_scale_str: str, T_in_GeV: float, omega_prefactor_input, Npoints: int,
            Nloop: int, order: str, corr: str, mu_IR_by_T):
    params = get_parameters(Nf, max_type, min_scale_str, T_in_GeV, omega_prefactor_input, Npoints, Nloop, order, corr,
                            mu_IR_by_T)

    cBsq = get_cBsq(params)

    OmegaByT_arr, PhiUVByT3, g2_arr = lpd.parallel_function_eval(inner_loop,
                                                                 range(len(params.OmegaByT_values)),
                                                                 128,
                                                                 params,
                                                                 cBsq)

    return np.asarray(OmegaByT_arr), np.asarray(g2_arr), np.asarray(PhiUVByT3)


def save_UV_spf(args, OmegaByT_arr, g2_arr, PhiUVByT3):
    # prepare file names
    try:
        float(args.min_scale)
        min_scale_label = '{0:.2f}'.format(args.min_scale)
    except ValueError:
        min_scale_label = args.min_scale
    try:
        float(args.omega_prefactor)
        omega_prefactor_label = '{0:.1f}'.format(args.omega_prefactor)
    except ValueError:
        omega_prefactor_label = args.omega_prefactor
    if args.max_type == "smooth":
        max_label = "smax"
    elif args.max_type == "hard":
        max_label = "hmax"
    if args.suffix != "":
        args.suffix = "_" + args.suffix
    if args.prefix != "":
        args.prefix = args.prefix + "_"
    prefix = args.outputpath + "/" + args.prefix
    middle_part = "Nf" + str(args.Nf) + "_" + '{0:.3f}'.format(
        args.T_in_GeV) + "_" + min_scale_label + "_" + omega_prefactor_label + "_" + max_label

    # save files

    file = prefix + "g2_" + middle_part + args.suffix + ".npy"
    print("saving ", file)
    np.save(file, np.column_stack((OmegaByT_arr, g2_arr)))

    file = prefix + "SPF_" + args.corr + "_" + args.order + "_" + middle_part + args.suffix + ".npy"  # format: 'omega/T g^2'
    print("saving ", file)
    np.save(file, np.column_stack((OmegaByT_arr, PhiUVByT3)))  # format: 'omega/T rho/T^3'


def add_args(parser):
    """ arguments that are relevant for get_spf() """
    parser.add_argument("--T_in_GeV", type=float, required=True)
    parser.add_argument("--Nf", help="number of flavors", type=int, required=True)

    parser.add_argument("--min_scale", default=1.3,
                        help="choices: piT, 2piT, or any number in GeV. limits how far the coupling will be run down. "
                             "From lattice calculations of the moments of quarkonium "
                             "correlators we know that going to lower than 1.3 GeV is probably unreasonable for Nf=3.")
    parser.add_argument("--omega_prefactor", help="mu = prefactor * omega. choices: opt, or any number.", default=1)
    parser.add_argument("--max_type", choices=["smooth", "hard"],
                        help="whether to use max(a,b) (hard) or sqrt(a^2+b^2) (smooth) maximum to determine the scale",
                        default="hard")

    parser.add_argument("--Npoints", default=10000, type=int, help="number of points between min and max OmegaByT")
    parser.add_argument("--Nloop", help="number of loops", type=int, default=5)
    parser.add_argument("--order", help="perturbative order. decide between LO or NLO.", choices=["LO", "NLO"],
                        type=str, required=True)
    parser.add_argument("--corr",
                        help="decide between EE or BB correlator function forms. LO is identical but the scale we use may change.",
                        type=str, choices=["EE", "BB"], required=True)
    parser.add_argument("--mu_IR_by_T", help="In the BB case, this is the scale to which the spectral function should be run", default=np.nan)

    return


def main():
    parser = argparse.ArgumentParser()
    add_args(parser)

    parser.add_argument("--suffix", type=str, help="string to append to end of output file name", default="")
    parser.add_argument("--prefix", type=str, help="string to prepend to end of output file name", default="")
    parser.add_argument("--outputpath", default="/work/home/altenkort/work/correlators_flow/data/merged/spf_coupling/")
    # parser.add_argument("--nproc", default=1, type=int, help="number of processes")

    args = parser.parse_args()

    OmegaByT_arr, g2_arr, PhiUVByT3 = get_spf(args.Nf, args.max_type, args.min_scale, args.T_in_GeV,
                                              args.omega_prefactor, args.Npoints, args.Nloop, args.order, args.corr, args.mu_IR_by_T)

    save_UV_spf(args, OmegaByT_arr, g2_arr, PhiUVByT3)


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
