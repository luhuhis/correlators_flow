#!/usr/bin/env python3
import lib_process_data as lpd
import numpy
import argparse



def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--test')
    args = parser.parse_args()
    return args


def main():
    data = [[64, 20, 5899],
    [64, 24, 3435],
    [96, 36, 2256],
    [64, 20, 7923],
    [64, 24, 2715],
    [96, 32, 912],
    [64, 20, 6786],
    [64, 24, 5325],
    [96, 28, 1680],
    [64, 20, 6534],
    [64, 22, 9101],
    [96, 24, 688]]

    sum = 0

    BYTES = 8
    NSU3 = 18
    NDIM = 4
    COMPRESSION = 2/3

    for dat in data:
        latfactor = COMPRESSION * NDIM * NSU3 * BYTES * dat[0]**3 * dat[1]
        print(f"lattice {dat[0]} {dat[1]}: {latfactor/10**9} GB")
        sum += latfactor * dat[2]
        print(sum/10**12)

    args = parse_args()

    return


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
