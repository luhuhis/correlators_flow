#!/usr/local/bin/python3.7m -u
import rundec
import numpy as np
import argparse
import lib_process_data as lpd


def smooth_max(a, b):
    return np.sqrt(a**2+b**2)


def main():

    #TODO rename file to "EE_pert_spf.py"

    parser = argparse.ArgumentParser()

    parser.add_argument("--T_in_GeV", type=float, required=True)
    parser.add_argument("--Nf", help="number of flavors", type=int, required=True)

    parser.add_argument("--min_scale", default=1.3,
                        help="choices: piT, 2piT, or any number in GeV. limits how far the coupling will be run down. "
                             "From lattice calculations of the moments of quarkonium "
                             "correlators we know that going to lower than 1.3 GeV is probably unreasonable for Nf=3.")
    parser.add_argument("--omega_prefactor", help="mu = prefactor * omega. choices: opt, or any number.", default=1)
    parser.add_argument("--max_type", choices=["smooth", "hard"], help="whether to use max(a,b) (hard) or sqrt(a^2+b^2) (smooth) maximum to determine the scale",
                        default="smooth")

    parser.add_argument("--outputpath", default="/work/home/altenkort/work/correlators_flow/data/merged/spf_coupling/")
    parser.add_argument("--suffix", type=str, help="string to append to end of output file name", default="")
    parser.add_argument("--prefix", type=str, help="string to prepend to end of output file name", default="")

    parser.add_argument("--min_OmegaByT_exp", help="exponent to base 10 for minimum OmegaByT", type=float, default=-6)
    parser.add_argument("--max_OmegaByT_exp", help="exponent to base 10 for maximum OmegaByT", type=float, default=3)
    parser.add_argument("--Npoints", default=100000, type=int, help="number of points between min and max OmegaByT")
    parser.add_argument("--Nloop", help="number of loops", type=int, default=5)

    args = parser.parse_args()

    # check Nf
    if args.Nf == 0:
        Lambda_MSbar = 0.253  # (JHEP11(2017)206)
    elif args.Nf == 3:
        Lambda_MSbar = 0.339  # 2111.09849 (p.11)
    else:
        print("ERROR: I don't know any LambdaMSbar for this Nf")
        exit(1)

    if args.max_type == "smooth":
        maxfunc = smooth_max
        max_label = "smax"
    elif args.max_type == "hard":
        maxfunc = np.fmax
        max_label = "hmax"

    # check min_scale
    if args.min_scale == "piT":
        args.min_scale = np.pi * args.T_in_GeV
        min_scale_label = "piT"
    elif args.min_scale == "2piT":
        args.min_scale = 2 * np.pi * args.T_in_GeV
        min_scale_label = "2piT"
    elif args.min_scale == "3piT":
        args.min_scale = 3 * np.pi * args.T_in_GeV
        min_scale_label = "3piT"
    else:
        args.min_scale = float(args.min_scale)
        min_scale_label = '{0:.2f}'.format(args.min_scale)

    Nc = 3

    # check omega_prefactor
    if args.omega_prefactor == "opt":
        mu_opt_omega_term = ((24 * np.pi ** 2 - 149) * Nc + 20 * args.Nf) / (6 * (11 * Nc - 2 * args.Nf))  # see eq. 4.17 of arXiv:1006.0867
        args.omega_prefactor = 2*np.exp(mu_opt_omega_term)  # ~14
        omega_prefactor_label = "opt"
    else:
        args.omega_prefactor = float(args.omega_prefactor)
        omega_prefactor_label = '{0:.1f}'.format(args.omega_prefactor)

    crd = rundec.CRunDec()
    r20 = Nc * (149. / 36. - 11. * np.log(2.) / 6. - 2 * np.pi ** 2 / 3.) - args.Nf * (5. / 9. - np.log(2.) / 3.)
    r21 = (11. * Nc - 2. * args.Nf) / 12.
    C_F = (Nc ** 2 - 1) / 2 / Nc

    # mu_opt_T_term = - np.euler_gamma - (Nc - 8 * np.log(2) * args.Nf) / (2*(11 * Nc - 2 * args.Nf))

    # calculation
    LO_SPF = []
    NLO_SPF = []
    g2_arr = []
    for OmegaByT in np.logspace(-6, 3, args.Npoints, base=10):
        mu = maxfunc(args.min_scale, args.omega_prefactor*OmegaByT*args.T_in_GeV)
        Alphas = crd.AlphasLam(Lambda_MSbar, mu, args.Nf, args.Nloop)
        g2 = 4. * np.pi * Alphas
        g2_arr.append([OmegaByT, g2])
        lo_spf = g2 * C_F * OmegaByT**3 / 6. / np.pi
        l = np.log((mu / args.T_in_GeV) ** 2 / OmegaByT ** 2)
        LO_SPF.append([OmegaByT, lo_spf])
        NLO_SPF.append([OmegaByT, lo_spf * (1 + (r20 + r21 * l) * Alphas / np.pi)])
    g2_arr = np.asarray(g2_arr)
    LO_SPF = np.array(LO_SPF)
    NLO_SPF = np.array(NLO_SPF)

    # save files
    if args.suffix != "":
        args.suffix = "_"+args.suffix
    if args.prefix != "":
        args.prefix = args.prefix+"_"

    prefix = args.outputpath+"/" + args.prefix
    middle_part = "Nf" + str(args.Nf) + "_" + '{0:.3f}'.format(args.T_in_GeV) + "_" + min_scale_label + "_" + omega_prefactor_label + "_" + max_label
    file = prefix + "g2_" + middle_part + args.suffix + ".dat"
    print("saving ", file)
    np.savetxt(file, np.column_stack((g2_arr[:, 0], g2_arr[:, 1])), fmt='%10.9e', header='omega/T g^2')

    file = prefix + "SPF_LO_"+middle_part +args.suffix+".dat"
    print("saving ", file)
    np.savetxt(file, np.column_stack((LO_SPF[:, 0], LO_SPF[:, 1])), fmt='%10.9e', header='omega/T rho/T^3')

    file = prefix + "SPF_NLO_"+middle_part +args.suffix+".dat"
    print("saving ", file)
    np.savetxt(file, np.column_stack((NLO_SPF[:, 0], NLO_SPF[:, 1])), fmt='%10.9e', header='omega/T rho/T^3')


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
