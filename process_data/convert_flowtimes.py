#!/usr/local/bin/python3.7 -u
import lib_process_data as lpd
import numpy
import argparse


def main():

    parser = argparse.ArgumentParser("script to convert flow times at fixed temperature from one lattice spacing to another")
    parser.add_argument('--input', type=str, help="path to input file that contains dimensionless lattice flow times t_F = tau_F / a^2", required=True)
    parser.add_argument('--output', type=str, help="path to output file", required=True)
    parser.add_argument('--input_Nt', type=float, help="Nt of input lattice at fixed temperature T")
    parser.add_argument('--output_Nt', type=float, help="Nt of output lattice at fixed temperature T", required=True)
    parser.add_argument('--input_a', type=float, help="lattice spacing of input lattice")
    parser.add_argument('--output_a', type=float, help="lattice spacing of output lattice")
    parser.add_argument('--threshold', default=0.3, type=float, help="maximum flow RADIUS in units of temperature")
    parser.add_argument('--type', choices=["fixed_temperature", "different_temperature"])
    args = parser.parse_args()

    if args.type == "fixed_temperature":
        if not args.input_Nt or not args.output_Nt:
            print("ERROR: For type=fixed_temperature need input_Nt and output_Nt")
            exit(1)
    if args.type == "different_temperature":
        if not args.input_a or not args.output_a:
            print("ERROR: For type=different_temperature need input_a and output_a")
            exit(1)
        if args.input_a > args.output_a:
            print("WARN: input_a should be smaller than output_a, to ensure that stepsizes do not increase.")

    flowtimes_input = numpy.loadtxt(args.input)
    flowtimes_output = []
    one_more = True
    for tf in flowtimes_input:
        if args.type == "fixed_temperature":
            tf_prime = tf/args.input_Nt**2 * args.output_Nt**2
        elif args.type == "different_temperature":
            tf_prime = tf * args.input_a ** 2 / args.output_a ** 2
        if numpy.sqrt(8*tf_prime)/args.output_Nt <= args.threshold:
            flowtimes_output.append(tf_prime)
        elif one_more:
            flowtimes_output.append(tf_prime)
            one_more = False

    print(numpy.asarray(flowtimes_output))
    numpy.savetxt("flowtimes_"+args.output+"_Nt"+str(int(args.output_Nt))+".txt", flowtimes_output, fmt='%.9f', newline=' ')

    return


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
