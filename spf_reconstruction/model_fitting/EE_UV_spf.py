#!/usr/bin/env python3
import rundec
import numpy as np
import argparse
import lib_process_data as lpd


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


def set_minscale(min_scale, T_in_GeV, Nc, Nf):

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


def set_omega_prefactor(omega_prefactor, Nc, Nf, T_in_GeV, min_scale):

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


def get_spf(Nf: int, max_type: str, min_scale, T_in_GeV, omega_prefactor_input, Npoints, Nloop):
    # check Nf
    if Nf == 0:
        Lambda_MSbar = 0.253  # (JHEP11(2017)206)
    elif Nf == 3:
        Lambda_MSbar = 0.339  # 2111.09849 (p.11)
    else:
        print("ERROR: I don't know any LambdaMSbar for this Nf")
        exit(1)

    if max_type == "smooth":
        maxfunc = smooth_max
    elif max_type == "hard":
        maxfunc = np.fmax

    # some constants
    Nc = 3
    r20 = Nc * (149. / 36. - 11. * np.log(2.) / 6. - 2 * np.pi ** 2 / 3.) - Nf * (5. / 9. - np.log(2.) / 3.)
    r21 = (11. * Nc - 2. * Nf) / 12.
    C_F = (Nc ** 2 - 1) / 2 / Nc
    b_0 = 11 / (16 * np.pi**2)
    gamma_0 = 3 / (8 * np.pi**2)

    min_scale = set_minscale(min_scale, T_in_GeV, Nc, Nf)
    omega_prefactor, omega_exponent = set_omega_prefactor(omega_prefactor_input, Nc, Nf, T_in_GeV, min_scale)

    crd = rundec.CRunDec()

    # calculation
    OmegaByT_arr = []
    LO_SPF = []
    NLO_SPF = []
    g2_arr = []
    first = True
    for OmegaByT in np.logspace(-6, 3, Npoints, base=10):
        scale = omega_prefactor * (OmegaByT*T_in_GeV)**omega_exponent
        mu = maxfunc(min_scale, scale)
        if scale > min_scale and first:
            print("scale > min_scale at OmegaByT=", lpd.format_float(OmegaByT))
            first = False
        Alphas = crd.AlphasLam(Lambda_MSbar, mu, Nf, Nloop)
        g2 = 4. * np.pi * Alphas
        g2_arr.append(g2)
        lo_spf = g2 * C_F * OmegaByT ** 3 / 6. / np.pi
        l = np.log((scale / T_in_GeV) ** 2 / OmegaByT ** 2)  # TODO should this be scale or mu?

        if omega_prefactor_input == "optBB":
            BB_inner_bracket = (b_0-gamma_0) * np.log(scale**2 / (OmegaByT*T_in_GeV)**2) + gamma_0 * np.log(scale**2/min_scale**2)
            BB_bracket = (1 + g2 * BB_inner_bracket)
            NLO_SPF.append(lo_spf * BB_bracket)
        else:
            NLO_SPF.append(lo_spf * (1 + (r20 + r21 * l) * Alphas / np.pi))
        OmegaByT_arr.append(OmegaByT)
        LO_SPF.append(lo_spf)

    OmegaByT_arr = np.asarray(OmegaByT_arr)
    g2_arr = np.asarray(g2_arr)
    LO_SPF = np.array(LO_SPF)
    NLO_SPF = np.array(NLO_SPF)
    return OmegaByT_arr, g2_arr, LO_SPF, NLO_SPF


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
