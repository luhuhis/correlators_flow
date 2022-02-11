#!/usr/local/bin/python3.7m -u
import rundec
import numpy as np
import argparse


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("--TinGeV", type=float, required=True)
    parser.add_argument("--suffix", type=str, help="file will be called ..._suffix.dat")
    parser.add_argument("--qcdtype", choices=["quenched", "hisq"])

    args = parser.parse_args()

    if args.qcdtype == "quenched":
        Nf = 0
        Nloop = 4
        Lambda_MSbar = 0.253  # from JHEP11(2017)206
    elif args.qcdtype == "hisq":
        Nf = 3  # quark mass effects neglected. from ALPHA collaboration 1701.03075
        Nloop = 5
        Lambda_MSbar = 0.332

    crd = rundec.CRunDec()
    Nc = 3
    r20 = Nc * (149. / 36. - 11. * np.log(2.) / 6. - 2 * np.pi ** 2 / 3.) - Nf * (5. / 9. - np.log(2.) / 3.)
    r21 = (11. * Nc - 2. * Nf) / 12.
    C_F = (Nc ** 2 - 1) / 2 / Nc

    # T=1.5*1.24*Lambda_MSbar ######quenched at T=1.5T_c. T_c=1.24*Lambda_MSbar

    delta_omega = 0.1 * args.TinGeV  # binsize=0.1 from Mikko

    LO_SPF = []
    NLO_SPF = []
    g2_arr = []
    for i in range(1, 10001):
        mu = np.pi * args.TinGeV if i * delta_omega < np.pi * args.TinGeV else i * delta_omega
        a = crd.AlphasLam(Lambda_MSbar, mu, Nf, Nloop)  # 3 flavors, 5 loop
        g2 = 4. * np.pi * a
        OmegaByT = 0.1 * i
        g2_arr.append([OmegaByT, g2])
        lo_spf = g2 * C_F * OmegaByT ** 3. / 6. / np.pi
        l = np.log((mu / args.TinGeV) ** 2 / OmegaByT ** 2)
        LO_SPF.append([OmegaByT, lo_spf])
        NLO_SPF.append([OmegaByT, lo_spf * (1 + (r20 + r21 * l) * a / np.pi)])
    g2_arr = np.asarray(g2_arr)
    LO_SPF = np.array(LO_SPF)
    NLO_SPF = np.array(NLO_SPF)
    np.savetxt("Nf" + str(Nf) + "_g2_" + args.suffix + ".dat", np.column_stack((g2_arr[:, 0], g2_arr[:, 1])), fmt='%10.9e', header='omega/T g^2')
    np.savetxt("Nf"+str(Nf)+"_LO_SPF_"+args.suffix+".dat", np.column_stack((LO_SPF[:, 0], LO_SPF[:, 1])), fmt='%10.9e', header='omega/T rho/T^3')
    np.savetxt("Nf"+str(Nf)+"_NLO_SPF_"+args.suffix+".dat", np.column_stack((NLO_SPF[:, 0], NLO_SPF[:, 1])), fmt='%10.9e', header='omega/T rho/T^3')
