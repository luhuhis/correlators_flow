#!/usr/bin/env python3
import lib_process_data as lpd
import numpy
import argparse



def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--kappaByT3', type=float)
    parser.add_argument('--r0_in_fm', type=float, default=0.47, help="default from https://arxiv.org/abs/1401.3270")
    parser.add_argument('--r0Tc', type=float, default=0.7457, help="default from https://arxiv.org/abs/1503.05652")
    parser.add_argument('--T_in_GeV', type=float)
    parser.add_argument('--TbyTc', type=float)
    parser.add_argument('--M_in_GeV', type=float)
    parser.add_argument('--output_digits', type=int, default=3)
    args = parser.parse_args()
    return args


def taukin_in_fm_quenched(r0_in_fm, r0Tc, fm_by_GeV, kappaByT3, M_in_GeV, TbyTc):
    return 2 * r0_in_fm**2 / r0Tc**2 / fm_by_GeV * M_in_GeV / TbyTc**2 / kappaByT3


def taukin_in_fm_hisq(fm_by_GeV, T_in_GeV, kappa_byT3, M_in_GeV):
    return 2 / kappa_byT3 * M_in_GeV / T_in_GeV**2 * fm_by_GeV



def print_hisq(args, kappas, T_in_GeV):
    fm_by_GeV = 0.19733
    n = args.output_digits

    taukin_charm_low = taukin_in_fm_hisq(fm_by_GeV, T_in_GeV, kappas[1], 1.275)
    taukin_charm_high = taukin_in_fm_hisq(fm_by_GeV, T_in_GeV, kappas[0], 1.275)
    taukin_bottom_low = taukin_in_fm_hisq(fm_by_GeV, T_in_GeV, kappas[1], 4.18)
    taukin_bottom_high = taukin_in_fm_hisq(fm_by_GeV, T_in_GeV, kappas[0], 4.18)
    twopiTD_low = 2 * numpy.pi * 2 / kappas[1]
    twopiTD_high = 2 * numpy.pi * 2 / kappas[0]

    print("$", lpd.format_float(kappas[0], n), "\dots", lpd.format_float(kappas[1], n), "$ &",
          "$", lpd.format_float(twopiTD_low, n), "\dots", lpd.format_float(twopiTD_high, n), "$ &",
          "$", lpd.format_float(taukin_charm_low, n), "\dots", lpd.format_float(taukin_charm_high, n), "$ &",
          "$", lpd.format_float(taukin_bottom_low, n), "\dots", lpd.format_float(taukin_bottom_high, n), r'$ \\')


def main():

    args = parse_args()
    if not args.T_in_GeV and args.TbyTc:
        print("using r0Tc, r0 and TbyTc to convert")
        taukin = taukin_in_fm_quenched(args.r0_in_fm, args.r0Tc, fm_by_GeV, args.kappaByT3, args.M_in_GeV, args.TbyTc)
    else:

        Ts=[0.195, 0.220, 0.251, 0.293]
        kappas=[
            [8.544, 13.448], #[8.458, 13.419],
            [6.032, 10.759], #[5.834, 10.702],
            [4.722, 9.039], #[4.966, 8.913],
            [3.85, 7.70] #[3.991, 7.610]
        ]
        print("T [MeV] & kappa/T^3 & 2piTD & tau_kin [fm]")
        for i, T in enumerate(Ts):
            print(T, end=" & ")
            print_hisq(args, kappas[i], T)

    return


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
