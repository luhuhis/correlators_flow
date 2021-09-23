import lib_process_data as lpd
import sys
import numpy


def get_args():
    """ parse cmd line arguments """
    parser, requiredNamed = lpd.get_parser()

    #parser.add_argument('--only_plot_no', type=int, nargs='*', default=[1, 2, 3])

    requiredNamed.add_argument('--conftype1', help="format: s064t64_b0824900_m002022_m01011", required=True)
    requiredNamed.add_argument('--conftype2', help="format: s064t64_b0824900_m002022_m01011", required=True)

    args = parser.parse_args()

    #if 'conftype_2' in vars(args) and 'reconstruct' not in vars(args):
        #parser.error('The --conftype_2 argument requires the --reconstruct argument!')

    return args


def load_data(args, inputfolder: str, prefix_load: str, nt: int):

    flow_radius = numpy.loadtxt(inputfolder + "flowradii_" + args.conftype + ".dat")
    flow_times = numpy.loadtxt(inputfolder + "flowtimes_" + args.conftype + ".dat")
    XX = numpy.loadtxt(inputfolder + args.corr + prefix_load + "_" + args.conftype + ".dat")
    XX_err = numpy.loadtxt(inputfolder + args.corr + prefix_load + "_err_" + args.conftype + ".dat")

    tauT = lpd.get_tauTs(nt)

    valid_flowtimes = None
    delete_indices = []
    # filter out flowtimes based on the file input
    if args.flowselectionfile:
        valid_flowtimes = numpy.loadtxt(args.flowselectionfile)
        for i, val in enumerate(flow_times):
            if val not in valid_flowtimes:
                delete_indices.append(i)
        flow_times = numpy.delete(flow_times, delete_indices, 0)
        flow_radius = numpy.delete(flow_radius, delete_indices, 0)
        XX = numpy.delete(XX, delete_indices, 0)
        XX_err = numpy.delete(XX_err, delete_indices, 0)

    return flow_radius, flow_times, XX, XX_err, tauT, valid_flowtimes


def main():

    # this script plots a custom ratio of a correlator at Nt1 with RECONSTRUCTED correlator from Nt2

    # parse arguments
    args = get_args()
    beta1, ns1, nt1, nt_half1 = lpd.parse_conftype(args.conftype1)
    beta2, ns2, nt2, nt_half2 = lpd.parse_conftype(args.conftype2)
    
    if nt1 != nt_half2:
        sys.exit("Nt of conftype1 has to be half of Nt of conftype2")
    
    #fermions, temp, flowtype = lpd.parse_qcdtype(args.qcdtype)
    #if fermions not in ("quenched", "hisq"):
        #sys.exit("couldn't parse fermion type (either fermions=hisq or fermions=quenched) from qcdtype. use this syntax for qcdtype: <fermions>_<other_stuff>")

    outputfolder = lpd.get_plot_path(args.qcdtype, args.corr, args.conftype)
    lpd.create_folder(outputfolder)
    inputfolder = lpd.get_merged_data_path(args.qcdtype, args.corr, args.conftype)

    # load data
    flow_radius, flow_times, XX, XX_err, tauT, valid_flowtimes = load_data(args, inputfolder, prefix_load, nt)
    XX_bak = numpy.copy(XX)
    XX_err_bak = numpy.copy(XX_err)


if __name__ == '__main__':
    main()
