#!/usr/local/bin/python3.7 -u
import lib_process_data as lpd
import numpy
import argparse


def main():

    parser = argparse.ArgumentParser("script to convert flow times at fixed temperature from one lattice spacing to another")
    parser.add_argument('--input', type=str, help="path to input file that contains dimensionless lattice flow times t_F = tau_F / a^2", required=True)
    parser.add_argument('--output', type=str, help="path to output file", required=True)
    parser.add_argument('--input_Nt', type=float, help="Nt of input lattice at temperature T", required=True)
    parser.add_argument('--output_Nt', type=float, help="Nt of output lattice at temperature T", required=True)
    parser.add_argument('--threshold', default=0.3, type=float, help="maximum flow RADIUS in units of temperature")
    args = parser.parse_args()

    flowtimes_input = numpy.loadtxt(args.input)
    flowtimes_output = []
    one_more = True
    for tf in flowtimes_input:
        tf_prime = tf/args.input_Nt**2 * args.output_Nt**2
        if numpy.sqrt(8*tf_prime)/args.output_Nt <= args.threshold:
            flowtimes_output.append(tf_prime)
        elif one_more:
            flowtimes_output.append(tf_prime)
            one_more = False

    print(numpy.asarray(flowtimes_output))
    numpy.savetxt(args.output+"_"+str(args.output_Nt)+".txt", flowtimes_output, fmt='%.9f', newline=' ')

    return


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
