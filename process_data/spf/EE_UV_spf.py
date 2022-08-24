#!/usr/local/bin/python3.7m -u
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
                        default="smooth")

    parser.add_argument("--Npoints", default=100000, type=int, help="number of points between min and max OmegaByT")
    parser.add_argument("--Nloop", help="number of loops", type=int, default=5)
    return


def get_spf(Nf: int, max_type: str, min_scale, T_in_GeV, omega_prefactor, Npoints, Nloop):
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

    Nc = 3

    mu_opt_T_term = np.exp(np.log(4*np.pi*T_in_GeV) - np.euler_gamma - (Nc - 8 * np.log(2) * Nf) / (2*(11 * Nc - 2 * Nf)))
    # check min_scale
    if min_scale == "piT":
        min_scale = np.pi * T_in_GeV
    elif min_scale == "2piT":
        min_scale = 2 * np.pi * T_in_GeV
    elif min_scale == "3piT":
        min_scale = 3 * np.pi * T_in_GeV
    elif min_scale == "4piT":
        min_scale = 4 * np.pi * T_in_GeV
    elif min_scale == "opt":
        min_scale = mu_opt_T_term
    elif min_scale == "0.25opt":
        min_scale = mu_opt_T_term / 4
    elif min_scale == "0.5opt":
        min_scale = mu_opt_T_term / 2
    elif min_scale == "2opt":
        min_scale = mu_opt_T_term * 2
    else:
        min_scale = float(min_scale)

    print("min_scale = ", min_scale)

    # check omega_prefactor
    mu_opt_omega_term = 2 * np.exp(((24 * np.pi ** 2 - 149) * Nc + 20 * Nf) / (6 * (11 * Nc - 2 * Nf)))  # see eq. 4.17 of arXiv:1006.0867
    if omega_prefactor == "opt":
        omega_prefactor = mu_opt_omega_term  # ~14 for Nf=3, ~7.6 for Nf=0
    elif omega_prefactor == "0.25opt":
        omega_prefactor = mu_opt_omega_term / 4
    elif omega_prefactor == "0.5opt":
        omega_prefactor = mu_opt_omega_term / 2
    elif omega_prefactor == "2opt":
        omega_prefactor = mu_opt_omega_term * 2
    else:
        omega_prefactor = float(omega_prefactor)

    print("omega_prefactor = ", omega_prefactor)

    crd = rundec.CRunDec()
    r20 = Nc * (149. / 36. - 11. * np.log(2.) / 6. - 2 * np.pi ** 2 / 3.) - Nf * (5. / 9. - np.log(2.) / 3.)
    r21 = (11. * Nc - 2. * Nf) / 12.
    C_F = (Nc ** 2 - 1) / 2 / Nc

    # mu_opt_T_term = - np.euler_gamma - (Nc - 8 * np.log(2) * args.Nf) / (2*(11 * Nc - 2 * args.Nf))

    # calculation
    LO_SPF = []
    NLO_SPF = []
    g2_arr = []
    first = True
    for OmegaByT in np.logspace(-6, 3, Npoints, base=10):
        scale = omega_prefactor * OmegaByT * T_in_GeV
        mu = maxfunc(min_scale, scale)
        if scale > min_scale and first:
            print("scale > min_scale at OmegaByT=", OmegaByT)
            first = False
        Alphas = crd.AlphasLam(Lambda_MSbar, mu, Nf, Nloop)
        g2 = 4. * np.pi * Alphas
        g2_arr.append([OmegaByT, g2])
        lo_spf = g2 * C_F * OmegaByT ** 3 / 6. / np.pi
        l = np.log((mu / T_in_GeV) ** 2 / OmegaByT ** 2)
        LO_SPF.append([OmegaByT, lo_spf])
        NLO_SPF.append([OmegaByT, lo_spf * (1 + (r20 + r21 * l) * Alphas / np.pi)])
    g2_arr = np.asarray(g2_arr)
    LO_SPF = np.array(LO_SPF)
    NLO_SPF = np.array(NLO_SPF)
    return g2_arr, LO_SPF, NLO_SPF


def main():
    parser = argparse.ArgumentParser()
    add_args(parser)

    parser.add_argument("--suffix", type=str, help="string to append to end of output file name", default="")
    parser.add_argument("--prefix", type=str, help="string to prepend to end of output file name", default="")
    parser.add_argument("--outputpath", default="/work/home/altenkort/work/correlators_flow/data/merged/spf_coupling/")

    args = parser.parse_args()

    g2_arr, LO_SPF, NLO_SPF = get_spf(args.Nf, args.max_type, args.min_scale, args.T_in_GeV, args.omega_prefactor, args.Npoints, args.Nloop)

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
        args.suffix = "_"+args.suffix
    if args.prefix != "":
        args.prefix = args.prefix+"_"
    prefix = args.outputpath+"/" + args.prefix
    middle_part = "Nf" + str(args.Nf) + "_" + '{0:.3f}'.format(args.T_in_GeV) + "_" + min_scale_label + "_" + omega_prefactor_label + "_" + max_label

    # save files
    file = prefix + "g2_" + middle_part + args.suffix + ".dat"
    print("saving ", file)
    np.savetxt(file, np.column_stack((g2_arr[:, 0], g2_arr[:, 1])), fmt='%10.9e', header='omega/T g^2')

    file = prefix + "SPF_LO_"+middle_part + args.suffix+".dat"
    print("saving ", file)
    np.savetxt(file, np.column_stack((LO_SPF[:, 0], LO_SPF[:, 1])), fmt='%10.9e', header='omega/T rho/T^3')

    file = prefix + "SPF_NLO_"+middle_part + args.suffix+".dat"
    print("saving ", file)
    np.savetxt(file, np.column_stack((NLO_SPF[:, 0], NLO_SPF[:, 1])), fmt='%10.9e', header='omega/T rho/T^3')


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
