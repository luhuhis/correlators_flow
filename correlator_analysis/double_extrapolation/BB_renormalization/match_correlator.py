#!/usr/bin/env python3
import lib_process_data as lpd
import numpy
import argparse



def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--test')
    args = parser.parse_args()
    return args


def load_data(args, basepath, flowsteps):  # TODO flowtimes
    samples = numpy.load(basepath + "/cont_extr/" + args.corr + "_cont_relflow_samples.npy")
    nt_finest_half = int(args.finest_Nt/2)
    cont_samples = samples[:, :, :nt_finest_half]
    print(cont_samples.shape)
    n_samples = len(cont_samples)
    cont_samples = numpy.swapaxes(cont_samples, 1, 2)
    data_std = lpd.dev_by_dist(cont_samples, axis=0)
    # shape: (1000, 18, 221)

    return cont_samples, data_std, n_samples

def convert_Gflow_to_GB_MSBAR_UV(Gflow, g2MSbar):
    pass


def save_data(GB_PHYS):
    pass


def run_GBMSBAR(GB_MSBAR_UV, g2MSbar, muUV, muIR):
    pass


def convert_GBMSBAR_IR_to_GB_PHYS(GB_MSBAR_IR, g2MSbar, muIR):
    pass






def compute_Zf():


def main():

    args = parse_args()

    # input: Gflow, g^2 MSbar, muUV



    Gflow, g2MSbar = load_data()

    GB_MSBAR_UV = convert_Gflow_to_GB_MSBAR_UV(Gflow, g2MSbar, muUV)
    GB_MSBAR_IR = run_GBMSBAR(GB_MSBAR_UV, g2MSbar, muUV, muIR)
    GB_PHYS = convert_GBMSBAR_IR_to_GB_PHYS(GB_MSBAR_IR, g2MSbar, muIR)
    save_data(GB_PHYS)


    return


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
