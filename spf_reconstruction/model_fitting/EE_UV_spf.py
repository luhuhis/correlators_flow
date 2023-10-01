#!/usr/bin/env python3
import rundec
import numpy as np
import argparse
import lib_process_data as lpd
import BB_UV_functions as BB


def smooth_max(a, b):
    return np.sqrt(a**2+b**2)


def add_args(parser):
    """ arguments that are relevant for get_spf() """
    parser.add_argument("--T_in_GeV", type=float, required=True)
    parser.add_argument("--Nf", help="number of flavors", type=int, required=True)

    parser.add_argument("--min_scale", default=1.3,
                        help="choices: piT, 2piT, or any number in GeV. limits how far the coupling will be run down. "
                             "From lattice calculations of the moments of quarkonium "
                             "correlators we know that going to lower than 1.3 GeV is probably unreasonable for Nf=3.")
    parser.add_argument("--omega_prefactor", help="mu = prefactor * omega. choices: opt, or any number.", default=1)
    parser.add_argument("--max_type", choices=["smooth", "hard"], help="whether to use max(a,b) (hard) or sqrt(a^2+b^2) (smooth) maximum to determine the scale",
                        default="hard")

    parser.add_argument("--Npoints", default=100000, type=int, help="number of points between min and max OmegaByT")
    parser.add_argument("--Nloop", help="number of loops", type=int, default=5)
    return


def get_minscale(min_scale, T_in_GeV, Nc, Nf):

    effective_theory_thermal_scale = 4 * np.pi * T_in_GeV * np.exp(-np.euler_gamma) * np.exp((-Nc + 4 * Nf * np.log(4)) / (22 * Nc - 4 * Nf))

    # check min_scale
    if min_scale == "piT":
        min_scale = np.pi * T_in_GeV
    elif min_scale == "2piT":
        min_scale = 2 * np.pi * T_in_GeV
    elif min_scale == "eff":
        min_scale = effective_theory_thermal_scale
    else:
        min_scale = float(min_scale)

    print("min_scale = ", lpd.format_float(min_scale))

    return min_scale


def get_omega_prefactor(omega_prefactor, Nc, Nf, T_in_GeV, min_scale):

    mu_opt_omega_term = 2 * np.exp(((24 * np.pi ** 2 - 149) * Nc + 20 * Nf) / (6 * (11 * Nc - 2 * Nf)))  # see eq. 4.17 of arXiv:1006.0867
    omega_exponent = 1
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


def get_constants(Nf):
    Nc = 3
    r20 = Nc * (149. / 36. - 11. * np.log(2.) / 6. - 2 * np.pi ** 2 / 3.) - Nf * (5. / 9. - np.log(2.) / 3.)
    r21 = (11. * Nc - 2. * Nf) / 12.
    C_F = (Nc ** 2 - 1) / 2 / Nc
    b_0 = 11 / (16 * np.pi**2)
    gamma_0 = 3 / (8 * np.pi**2)

    return Nc, r20, r21, C_F, b_0, gamma_0


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


def inner_loop(index, omega_prefactor, omega_exponent, omega_prefactor_input, OmegaByT, maxfunc, r20, r21, T_in_GeV, min_scale, crd, Lambda_MSbar, Nf, Nc, Nloop, C_F, cBsq):
    if omega_prefactor_input == "optBB":
        mu = np.sqrt(4 * (OmegaByT * T_in_GeV) ** 2 + min_scale ** 2)
        g2, Alphas, lo_spf = get_g2_and_alphas_and_lospf(crd, Lambda_MSbar, mu, Nf, Nloop, C_F, OmegaByT)

        # === NLO ===
        # The following equations follow C.23 of https://arxiv.org/pdf/2204.14075.pdf

        first_paren = 5./3. * np.log(mu**2/(2*OmegaByT*T_in_GeV)**2) + 134./9. - 2*np.pi**2/3.  # Guy changed the last term: -8*np.pi**2/3 to -2*np.pi**2/3.
        second_paren = 2./3 * np.log(mu**2/(2*OmegaByT*T_in_GeV)**2) + 26./9.
        mu_dep_part = lo_spf * (1 + g2 / (4 * np.pi)**2 * (Nc * first_paren - Nf * second_paren))
        mu_free_part = (g2**2*C_F)/(12*np.pi**3) * BB.get_mu_free(OmegaByT, Nc, Nf)  # from line 3 to the end of eq.(C.23)
        nlo_spf = cBsq[index] * (mu_dep_part + mu_free_part)

    else:
        scale = omega_prefactor * (OmegaByT * T_in_GeV) ** omega_exponent
        mu = maxfunc(min_scale, scale)
        g2, Alphas, lo_spf = get_g2_and_alphas_and_lospf(crd, Lambda_MSbar, mu, Nf, Nloop, C_F, OmegaByT)

        # === NLO ===
        # Technically there are two scales here: "mu" and "scale".
        # "mu" is a scale which is bounded from below (lower bound = thermal scale) and used to evaluate g2.
        # Above the thermal scale it is given by "scale".
        # "scale" also occurs in the NLO terms of the spf and is (if omega_prefactor_input=="opt") defined to exactly
        # cancel the NLO terms (including the log), even below the thermal scale. Therefore, it is not bounded from below.
        # Effectively, we recover the LO spf with an "NLO scale choice" this way since the NLO terms are always canceled.
        # However, all of this is a technicality since this distinction only matters at thermal scales where we do not
        # use the UV spf anyway, but plots are less confusing this way :)
        l_ = np.log((scale / T_in_GeV) ** 2 / OmegaByT ** 2)
        nlo_spf = lo_spf * (1 + (r20 + r21 * l_) * Alphas / np.pi)

    return OmegaByT, lo_spf, nlo_spf, g2


def get_spf(Nf: int, max_type: str, min_scale, T_in_GeV, omega_prefactor_input, Npoints, Nloop):
    Lambda_MSbar = get_lambdamsbar(Nf)
    maxfunc = get_maxfunc(max_type)
    Nc, r20, r21, C_F, b_0, gamma_0 = get_constants(Nf)
    min_scale = get_minscale(min_scale, T_in_GeV, Nc, Nf)
    omega_prefactor, omega_exponent = get_omega_prefactor(omega_prefactor_input, Nc, Nf, T_in_GeV, min_scale)

    crd = rundec.CRunDec()

    cBsq = BB.get_cBsq(min_scale, T_in_GeV, Npoints, gamma_0, Lambda_MSbar, Nf, Nloop)

    OmegaByT_arr, LO_SPF, NLO_SPF, g2_arr = zip(*[
        inner_loop(
            index, omega_prefactor, omega_exponent, omega_prefactor_input, OmegaByT, maxfunc,
            r20, r21, T_in_GeV, min_scale, crd, Lambda_MSbar, Nf, Nc, Nloop, C_F, cBsq
        )
        for index, OmegaByT in enumerate(np.logspace(-6, 3, Npoints, base=10))
    ])

    return np.asarray(OmegaByT_arr), np.asarray(g2_arr), np.asarray(LO_SPF), np.asarray(NLO_SPF)


def save_UV_spf(args, OmegaByT_arr, g2_arr, LO_SPF, NLO_SPF):
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
    middle_part = "Nf" + str(args.Nf) + "_" + '{0:.3f}'.format(args.T_in_GeV) + "_" + min_scale_label + "_" + omega_prefactor_label + "_" + max_label

    # save files

    file = prefix + "g2_" + middle_part + args.suffix + ".npy"
    print("saving ", file)
    np.save(file, np.column_stack((OmegaByT_arr, g2_arr)))

    file = prefix + "SPF_LO_" + middle_part + args.suffix + ".npy"  # format: 'omega/T g^2'
    print("saving ", file)
    np.save(file, np.column_stack((OmegaByT_arr, LO_SPF)))  # format: 'omega/T rho/T^3'

    file = prefix + "SPF_NLO_" + middle_part + args.suffix + ".npy"
    print("saving ", file)
    np.save(file, np.column_stack((OmegaByT_arr, NLO_SPF)))  # format: 'omega/T rho/T^3'


def main():
    parser = argparse.ArgumentParser()
    add_args(parser)

    parser.add_argument("--suffix", type=str, help="string to append to end of output file name", default="")
    parser.add_argument("--prefix", type=str, help="string to prepend to end of output file name", default="")
    parser.add_argument("--outputpath", default="/work/home/altenkort/work/correlators_flow/data/merged/spf_coupling/")

    args = parser.parse_args()

    OmegaByT_arr, g2_arr, LO_SPF, NLO_SPF = get_spf(args.Nf, args.max_type, args.min_scale, args.T_in_GeV, args.omega_prefactor, args.Npoints, args.Nloop)

    save_UV_spf(args, OmegaByT_arr, g2_arr, LO_SPF, NLO_SPF)


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
