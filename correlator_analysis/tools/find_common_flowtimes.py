#!/usr/bin/env python3
import lib_process_data as lpd
import numpy
import argparse


def float_intersection(a, b):
    return b[numpy.isclose(a[:, None], b).any(0)]


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--files', nargs='*', type=str)
    parser.add_argument('--basepath', type=str, default="")
    parser.add_argument('--output', type=str, required=True, help="output flowradii file")  # TODO CHANGE THIS TO FLOWTIMES, NOT RADII
    args = parser.parse_args()


    # TODO add option for different nt, otherwise this is useless since you've removed the flowradii files
    flowtimes = []
    for file in args.files:
        tmp = numpy.loadtxt(args.basepath+file)
        flowtimes.append(tmp)

    intersection = flowtimes[0]
    print(intersection.shape)
    for i in range(1, len(flowtimes)):
        intersection = float_intersection(flowtimes[i], intersection)
        print(intersection.shape)

    numpy.savetxt(args.output, intersection)
    return


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
