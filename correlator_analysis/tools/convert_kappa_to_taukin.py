#!/usr/bin/env python3
import lib_process_data as lpd
import numpy
import argparse



def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--kappaByT3', type=float)
    parser.add_argument('--r0_in_fm', type=float, default=0.47, help="default from https://arxiv.org/abs/1401.3270")
    parser.add_argument('--r0Tc', type=float, default=0.7457, help="default from https://arxiv.org/abs/1503.05652")
    parser.add_argument('--TbyTc', type=float)
    parser.add_argument('--M_in_GeV', type=float)
    parser.add_argument('--output_digits', type=int, default=3)
    args = parser.parse_args()
    return args


def taukin_in_fm(r0_in_fm, r0Tc, fm_by_GeV, kappaByT3, M_in_GeV, TbyTc):

    return 2 * r0_in_fm**2 / r0Tc**2 / fm_by_GeV * M_in_GeV / TbyTc**2 / kappaByT3


def main():

    args = parse_args()

    fm_by_GeV = 0.19733

    print("kappa/T^3= ", lpd.format_float(args.kappaByT3, args.output_digits))

    taukin = taukin_in_fm(args.r0_in_fm, args.r0Tc, fm_by_GeV, args.kappaByT3, args.M_in_GeV, args.TbyTc)
    print("tau_kin  = ", lpd.format_float(taukin, args.output_digits), "fm")

    twopiTD = 2*numpy.pi * 2/ args.kappaByT3
    print("2piTD    = ", lpd.format_float(twopiTD, args.output_digits))


    return


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
