#!/usr/local/bin/python3.7m
import lib_process_data as lpd
import numpy


def main():
    parser, requiredNamed = lpd.get_parser()

    parser.add_argument('--PathOutputFolder', help='the path of output folder like /a/b', type=str)

    # which quantitiy to plot
    requiredNamed.add_argument('--obs', help='which observable to plot', choices=["kappa", "chisqdof"], type=str, required=True)

    # identify the spf model
    requiredNamed.add_argument('--model', help='input file identifier. which model was used', choices=[2], type=int, required=True)
    requiredNamed.add_argument('--PhiUVtype', help='input file identifier. which PhiUV was used [ LO (a) or NLO (b) ]', type=str, choices=["a", "b"], required=True)
    parser.add_argument('--constrain', help='input file identifier. whether the spf was forced to reach the UV limit at large omega', action="store_true")
    requiredNamed.add_argument('--mu', help='input file identifier. which "en" function was used', choices=["alpha", "beta"], type=str, required=True)
    requiredNamed.add_argument('--nmax', help='input file identifier. which nmax was used. valid only for model 1,2.', choices=[1, 2, 3, 4, 5, 6], type=int, required=True)
    requiredNamed.add_argument('--nsamples', help='input file identifier. number of bootstrap samples that were used.', type=int)
    requiredNamed.add_argument('--tol', help='input file identifier. tolerance that was used', type=float)
    parser.add_argument('--start_from_mean_fit', help='input file identifier. whether the inital fit params where determined by a fit of the mean.', action="store_true")
    # parser.add_argument('--omega_idx', help='index which omega to plot', default=0, type=int)

    # plot args
    parser.add_argument('--nbins', help='the number of bins to plot. its also possible to specify a more sophisticated method like fd (Freedman-Diaconis).', default="fd")

    args = parser.parse_args()

    constrainstr = "s1" if not args.constrain else "s2"  # s1 = dont constrain, s2 = constrain
    startstr = "" if not args.start_from_mean_fit else "_m"
    if args.model == 2:
        modelidentifier = str(args.model) + "_" + args.PhiUVtype + "_" + constrainstr + "_" + str(args.mu) + "_" + str(args.nmax)
    elif args.model == 3:
        modelidentifier = str(args.model) + "_" + args.PhiUVtype
    fileidentifier = modelidentifier + "_" + str(args.nsamples) + "_" + '{:.0e}'.format(args.tol) + startstr

    # read in the normalized correlator
    inputfolder = "../"+lpd.get_merged_data_path(args.qcdtype, args.corr, "")+"/spf/"+fileidentifier+"/"
    outputfolder = inputfolder if not args.PathOutputFolder else args.PathOutputFolder


    index = None
    xlabel = None
    if args.obs == "kappa":
        index = 0
        xlabel = r'$\kappa / T^3$'
    # if args.obs == "spf":
    #     index = 1
    # if args.obs == "corr":
    #     index = 2
    if args.obs == "chisqdof":
        index = 3
        xlabel = r'$\chi^2/\mathrm{d.o.f.}$'

    samples = numpy.load(inputfolder+"samples_"+fileidentifier+".npy")
    structure = numpy.loadtxt(inputfolder+"samples_structure_"+fileidentifier+".dat", dtype=int)
    ab = structure[index]

    samples = samples[:, ab[0]:ab[1]]
    if args.obs == "kappa":
        samples = samples[:, 0]
    # if args.obs == "spf":
    #     samples = samples[:, args.omega_idx]

    fig, ax, plots = lpd.create_figure(xlims=(0, 15), ylabel=r'$n$', xlabel=xlabel, xlabelpos=(0.95, 0.2), ylabelpos=(0.02, 0.97), UseTex=False)

    _, bins, _ = ax.hist(samples, args.nbins)

    ax.set_title(str(len(bins))+' bins, '+fileidentifier)

    fig.savefig(outputfolder+"/"+args.corr+"_"+args.obs+"_"+fileidentifier+".pdf")


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
