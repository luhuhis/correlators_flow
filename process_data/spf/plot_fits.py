#!/usr/local/bin/python3.7m
import lib_process_data as lpd
import numpy


def main():
    parser, requiredNamed = lpd.get_parser()

    parser.add_argument('--PathOutputFolder', help='the path of output folder like /a/b', type=str)
    # parser.add_argument('--oldformat', help='read in data files of the old format', action="store_true")

    # identify the run
    requiredNamed.add_argument('--model', help='identify the spf model', type=int, choices=[2, 3], required=True)
    requiredNamed.add_argument('--nsamples', help='identify file input. number of bootstrap samples', type=int, required=True)
    requiredNamed.add_argument('--tol', help='identify file input. which tolerance was used', type=float)
    requiredNamed.add_argument('--start_from_mean_fit', help="identify file input.", action="store_true")
    parser.add_argument('--obs', help='whether to plot kappa or chisqdof', choices=["kappa", "chisqdof", "corr", "spf"], default="kappa")

    # TODO add support for choosing which models to plot

    # plot settings
    parser.add_argument('--wideaspect', action="store_true",
                        help="use a wide aspect ratio (basically make the plots larger)")
    parser.add_argument('--xlims', help='custom xlims for plot', nargs=2, type=float)
    parser.add_argument('--ylims', help='custom ylims for plot', nargs=2, type=float)
    parser.add_argument('--usetex', help='use latex to plot labels', action="store_true")

    args = parser.parse_args()

    # read in the normalized correlator
    inputfolder = "../"+lpd.get_merged_data_path(args.qcdtype, args.corr, "")+"/spf/"
    outputfolder = inputfolder if not args.PathOutputFolder else args.PathOutputFolder

    # file identifier
    startstr = "" if not args.start_from_mean_fit else "_m"
    suffix = "_" + str(args.nsamples) + "_" + '{:.0e}'.format(args.tol) + startstr

    # plot settings
    ylabel = None
    if args.obs == "kappa" or args.obs == "chisqdof":
        if args.obs == "kappa":
            xlims = (1.5, 12.5)
            idx = 0
            obslabel = r'\kappa/ T^3'
            filelabel = "params_"
            ylims = (0, 5)
        if args.obs == "chisqdof":
            xlims = (0, 9)
            idx = -1
            obslabel = r'\chi^2/\mathrm{d.o.f.}'
            filelabel = "params_"
            ylims = (0, 5)
        xlabel = r'$\mathrm{median}(' + obslabel + r') \pm 34^{\mathrm{th}}\, \%$'
        ylabelpos = None
    if args.obs == "corr" or args.obs == "spf":
        ylabelpos = (0.2, 0.9)
        if args.obs == "corr":
            obslabel = r'G^\mathrm{model}'
            filelabel = "corrfit_"
            ylims = (0.99, 1.03)
            ylabel = r'$\frac{\mathrm{median}('+obslabel+r') \pm 34^{\mathrm{th}}\, \%}{G^\mathrm{cont}}$'
            xlabel = r'$\tau T$'
            xlims = (0.2, 0.52)
        if args.obs == "spf":
            obslabel = r'\frac{\rho}{\omega^2 T}'
            xlabel = r'$\omega/T$'
            ylims = (0.1, 1000)
            filelabel = "spffit_"
            ylabel = r'$\mathrm{median}(' + obslabel + r') \pm 34^{\mathrm{th}}\, \%$'

    if args.model == 2:
        labels = ("2_a_s2_alpha_4",    "2_a_s2_alpha_5",    "2_a_s1_alpha_4",   "2_a_s1_alpha_5",     "2_a_s2_beta_4",    "2_a_s2_beta_5",    "2_a_s1_beta_4",
                  "2_a_s1_beta_5")
        labels_plot = (r'$s_2 \alpha a 4$', r'$s_2 \alpha a 5$', r'$s_1 \alpha a 4$', r'$s_1 \alpha a 5$', r'$s_2 \beta a 4$', r'$s_2 \beta a 5$',
                       r'$s_1 \beta a 4$', r'$s_1 \beta a 5$')
    elif args.model == 3:
        labels = ("3_a",)
        labels_plot = (r'$3a$',)
    labels = [label+suffix for label in labels]

    xdata = []
    ydata = []
    errorsleft = []
    errorsright = []
    for label in labels:
        try:
            file = inputfolder + "/"+label+"/" + filelabel + label + ".dat"
            print("read in", file)
            data = numpy.loadtxt(file)
        except OSError:
            print(label, "not found, skipping")
            xdata.append(numpy.nan)
            ydata.append(numpy.nan)
            errorsleft.append(numpy.nan)
            errorsright.append(numpy.nan)
            continue
        if args.obs == "kappa" or args.obs == "chisqdof":
            xdata.append(data[idx, 0])
            errorsleft.append(data[idx, 1])
            errorsright.append(data[idx, 2])
        elif args.obs == "corr":
            ydata.append(data[:, 3]/data[:, 1])
            xdata.append(data[:, 0])
            errorsleft.append(data[:, 4]/data[:, 1])
            errorsright.append(data[:, 5]/data[:, 1])
        elif args.obs == "spf":
            ydata.append(data[:, 1])
            xdata.append(data[:, 0])
            errorsleft.append(data[:, 2])
            errorsright.append(data[:, 3])

    figsize = (3 + 3 / 8, 3 + 3 / 8 - 1 / 2.54)
    if args.wideaspect:
        figsize = (1.5 * (3 + 3 / 8), 1.5 * (3 + 3 / 8) / 16 * 9)
    if args.xlims:
        xlims = args.xlims
    if args.ylims:
        ylims = args.ylims
    fig, ax, plots = lpd.create_figure(figsize=figsize, xlabel=xlabel, ylabel=ylabel, ylabelpos=ylabelpos, xlims=xlims, xlabelpos=(0.5, -0.1), ylims=ylims, UseTex=args.usetex)

    if args.obs == "kappa" or args.obs == "chisqdof":
        pos = 0.3+numpy.asarray((0.45, 0.85, 1.25, 1.65, 2.65, 3.05, 3.45, 3.85))
        colors = ('red', 'red', 'blue', 'blue', 'red', 'red', 'blue', 'blue')
        for i, label in enumerate(labels):
            ax.errorbar(xdata[i], pos[i], xerr=[[errorsleft[i]], [errorsright[i]]], color=colors[i], fmt='x-', fillstyle='none', markersize=5, mew=0.25, lw=0.8, elinewidth=0.5, capsize=1.2, zorder=-10)
        # custom yaxis
        ax.set_yticks(pos)
        ax.set_yticklabels(labels_plot)
        # vertical dashed lines
        currentxlims = ax.get_xlim()
        for i in range(0, int(currentxlims[1])):
            ax.axvline(i, **(lpd.chmap(lpd.verticallinestyle, alpha=0.4)))

    fmts = [val + "-" for val in lpd.markers]

    if args.obs == "corr" or args.obs == "spf":
        if args.obs == "corr":
            for i, label in enumerate(labels):
                plots.append(ax.errorbar(xdata[i]+i*0.003, ydata[i], yerr=[errorsleft[i], errorsright[i]], fmt=fmts[i], fillstyle='none', markersize=2, mew=0.25, lw=0.2, elinewidth=0.2, capsize=1.2, zorder=-10))
        if args.obs == "spf":
            PhiUV = numpy.loadtxt("/work/data/htshu/ee_spf/PHIUV_a.dat")
            PhiUV = PhiUV[:, 0:2]
            ax.set_yscale('log')
            ax.set_xscale('log')
            plots.append(ax.errorbar(PhiUV[:, 0], PhiUV[:, 1]/PhiUV[:, 0], fmt='-', label=r'$\Phi_a^\mathrm{UV}$', lw=0.4))
            for i, label in enumerate(labels):
                plots.append(ax.errorbar(xdata[i]+i*0.003, ydata[i]/xdata[i],  fmt=fmts[i], fillstyle='none', markersize=0, mew=0.25, lw=0.2, elinewidth=0.2, capsize=1.2, zorder=-10))  # yerr=[errorsleft[i], errorsright[i]],
        ax.legend(handles=plots, labels=(r'$\Phi^\mathrm{UV}$', *labels_plot), **lpd.legendstyle)
        ax.axhline(y=1, **lpd.horizontallinestyle)

    # title
    startstrplot = "naive" if not args.start_from_mean_fit else "mean"
    ax.set_title(r'$\mathrm{'+args.corr+r'},n_\mathrm{bs}='+str(args.nsamples)+r', \mathrm{tol}=\mathrm{'+'{:.0e}'.format(args.tol)+r'}, \mathrm{start}=\mathrm{'+startstrplot+r'}$')

    fig.savefig(outputfolder+"/"+args.obs+"_"+args.corr+suffix+".pdf")


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
