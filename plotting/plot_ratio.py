import lib_process_data as lpd
import sys
import numpy



def get_args():
    """ parse cmd line arguments """
    parser, requiredNamed = lpd.get_parser()

    # parser.add_argument('--only_plot_no', type=int, nargs='*', default=[1, 2, 3])

    requiredNamed.add_argument('--conftype1', help="format: s064t64_b0824900_m002022_m01011", required=True)
    requiredNamed.add_argument('--conftype2', help="format: s064t64_b0824900_m002022_m01011", required=True)

    args = parser.parse_args()

    # if 'conftype_2' in vars(args) and 'reconstruct' not in vars(args):
        # parser.error('The --conftype_2 argument requires the --reconstruct argument!')

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
    
    outputfolder = lpd.get_plot_path(args.qcdtype, args.corr, args.conftype1)
    lpd.create_folder(outputfolder)
    inputfolder1 = lpd.get_merged_data_path(args.qcdtype, args.corr, args.conftype1)
    inputfolder2 = lpd.get_merged_data_path(args.qcdtype, args.corr, args.conftype2)

    # load data
    flow_radius1, flow_times1, num1, num_err1, tauT1, valid_flowtimes1 = load_data(args, inputfolder1, "_numerator", nt1)
    flow_radius2, flow_times2, num2, num_err2, tauT2, valid_flowtimes2 = load_data(args, inputfolder2, "_numerator", nt2)
    flow_radius1, flow_times1, ploop1, ploop_err1, tauT1, valid_flowtimes1 = load_data(args, inputfolder1, "_polyakovloop", nt1)
    flow_radius2, flow_times2, ploop2, ploop_err2, tauT2, valid_flowtimes2 = load_data(args, inputfolder2, "_polyakovloop", nt2)

    for i, flowtime1 in enumerate(flow_times1):
        if flowtime1 != flow_times2[i]:
            print("ERROR: flow times don't seem to match")
            exit(1)

    num2_bak = numpy.copy(num2)
    num_err2_bak = numpy.copy(num_err2)

    # create reconstructed correlator
    for i in range(len(flow_radius2)):
        for j in range(nt_half1):  # from 0 to 16 for a reconstructed Nt=32 correlator based on Nt=64
            new_index = ((nt2-1) - (j + nt_half1))
            num2[i, j] += num2_bak[i, new_index]
            num_err2[i, j] = numpy.sqrt(num_err2_bak[i, j] ** 2 + num_err2_bak[i, new_index] ** 2)

    # calculate things to plot
    for i in range(len(flow_radius1)):
        for j in range(nt_half1):
            num_err1[i, j] = numpy.sqrt((num1[i, j]/num_err1[i, j])**2 + (num2[i, j]/num_err2[i, j])**2)
            num1[i, j] = num1[i, j] / num2[i, j]
            ploop1[i, j] = ploop1[i, j] / ploop2[i, j]

    # options for this plot
    figsize = (1.5 * (3 + 3 / 8), 1.5 * (3 + 3 / 8) / 16 * 9)
    xlims = [0, 0.505]
    xlabel = r'$\tau T$'
    ylabelpos = (0.01, 0.97)
    xlabelpos = (0.95, 0.06)
    ylims = ([1, 100000])
    ylabel = r'$\displaystyle\frac{G^\mathrm{latt }_{\tau_\mathrm{F}}(\tau)}{T_\mathrm{norm}^4 }$'

    # create the figure and axes objects
    fig, ax, plots = lpd.create_figure(xlims=xlims, ylims=ylims, xlabel=xlabel, xlabelpos=xlabelpos, ylabel=ylabel,
                                           ylabelpos=ylabelpos, UseTex=True, figsize=figsize)

    # save pdf
    filename = outputfolder + args.conftype1 + "_" + args.corr + "_" + args.conftype2 + ".pdf"
    fig.savefig(filename)
    print("saved correlator plot", filename)




if __name__ == '__main__':
    main()
