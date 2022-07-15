#!/usr/local/bin/python3.7m
import lib_process_data as lpd
import numpy
import argparse


def load_data(files, obs: str, verbose):
    # === load data
    xdata = []
    ydata = []
    errorsleft = []
    errorsright = []
    for file in files:
        try:
            data = numpy.loadtxt(file)
            if verbose:
                print("success: read in", file, "\n")
        except OSError:
            print("fail: could not find ", file, "\n")
            xdata.append(numpy.nan)
            ydata.append(numpy.nan)
            errorsleft.append(numpy.nan)
            errorsright.append(numpy.nan)
            continue
        if obs == "kappa" or obs == "chisqdof":
            if obs == "kappa":
                idx = 0
            elif obs == "chisqdof":
                idx = -1
            xdata.append(data[idx, 0])
            errorsleft.append(data[idx, 1])
            errorsright.append(data[idx, 2])
        elif obs == "corr":
            ydata.append((data[:, 3]/data[:, 1]))
            xdata.append(data[:, 0])
            errorsleft.append(numpy.abs(ydata[-1]) * numpy.sqrt((data[:, 4]/data[:, 3])**2 + (data[:, 2]/data[:, 1])**2))
            errorsright.append(numpy.abs(ydata[-1]) * numpy.sqrt((data[:, 5]/data[:, 3])**2 + (data[:, 2]/data[:, 1])**2))
        elif obs == "spf" or "spf_phys":
            ydata.append(data[:, 1])
            xdata.append(data[:, 0])
            errorsleft.append(data[:, 2])
            errorsright.append(data[:, 3])
        else:
            print("ERROR: UNKNOWN OBS")
            exit(1)
    return xdata, ydata, errorsleft, errorsright


def get_plot_settings(obs: str):
    ylabel = None
    ylims = None
    xlims = None
    xlabel = None
    ylabelpos = None
    filelabel = ""

    if obs == "kappa" or obs == "chisqdof":
        xlims = (150, 400)
        filelabel = "params_"
        ylabelpos = (0.05,0.9)
        xlabel = r'$T \mathrm{[MeV]}$'
        if obs == "kappa":
            obslabel = r'\kappa \mathrm{[GeV}^3\mathrm{]}'
            ylims = (0, 0.5)
            ylabel = r'$' + obslabel + '$'
        if obs == "chisqdof":
            obslabel = r'\chi^2/\mathrm{d.o.f.}'
            ylims = (0, 5)
            ylabel = r'$'+obslabel+'$'
    elif obs == "corr":
        filelabel = "corrfit_"
        ylabel = r'$\frac{G^\mathrm{fit}}{G^\mathrm{data}}$'
        xlabel = r'$\tau T$'
        xlims = (0.2, 0.52)
        ylims = (0.95, 1.05)
        ylabelpos = (0.05, 0.95)
    elif obs == "spf":
        obslabel = r'\frac{2\rho}{ \omega T^2}'
        xlabel = r'$\omega/T$'
        # xlims = ()
        # ylims = (1, 100)
        # --ylims  --xlims 0.001 1
        xlims = (0.1, 100)
        ylims = (1, 1000)
        filelabel = "spffit_"
        ylabel = r'$' + obslabel + r'$'
        ylabelpos = (0.05, 0.95)
    elif obs == "spf_phys":
        obslabel = r'\rho [\mathrm{GeV}^3]'
        xlabel = r'$\omega [\mathrm{GeV}]$'
        xlims = (0, 2)
        ylims = (0, 2)
        filelabel = "spffit_"
        ylabel = r'$' + obslabel + r'$'
        ylabelpos = (0.05, 0.95)
    else:
        print("ERROR: UNKNOWN OBS")
        exit(1)

    return filelabel, xlabel, ylabel, xlims, ylims, ylabelpos


def plot(obs, nfiles, xdata, ydata, errorsleft, errorsright, ax, plots, labels, PhiUV_file, plot_spf_err, pos, colors, Ntaus, kappa_swap_axes, omegaUV, legtitle):
    if obs == "kappa":

        for i in range(nfiles):
            if not numpy.isnan(xdata[i]).any():
                if kappa_swap_axes:
                    factor = 1
                    x = xdata[i] * factor
                    y = pos[i]
                    ax.errorbar(x, y, xerr=[[errorsleft[i] * factor], [errorsright[i] * factor]], color=colors[i], fmt='x-', fillstyle='none', markersize=5,
                                mew=0.25, lw=0.8, elinewidth=0.5, capsize=1.2, zorder=-10)
                    # custom yaxis
                    ax.set_yticks(pos[0:nfiles])
                    ax.set_yticklabels(labels)
                    # currentxlims = ax.get_xlim()
                    # for i in range(0, int(currentxlims[1])):
                    #     ax.axvline(i, **(lpd.chmap(lpd.verticallinestyle, alpha=0.4)))
                else:
                    factor = (pos[i]/1000)**3
                    x = pos[i]
                    y = xdata[i]*factor
                    ax.errorbar(x, y, yerr=[[errorsleft[i] * factor], [errorsright[i] * factor]], color=colors[i], fmt='x-', fillstyle='none', markersize=5,
                                mew=0.25, lw=0.8, elinewidth=0.5, capsize=1.2, zorder=-10)

        if obs == "kappa" and kappa_swap_axes:
            ax.axvline(14, **lpd.verticallinestyle)
            ax.tick_params(axis='y', which='major', labelsize=6)

        if obs == "kappa" and not kappa_swap_axes:
            xpoints = numpy.linspace(0, 500, 100)
            ax.errorbar(xpoints, 14 * (xpoints / 1000) ** 3, fmt='--', color='grey', alpha=0.8)
            # linear fit
            if not numpy.isnan(numpy.asarray(xdata)).any():
                import scipy.optimize
                xdatafit = pos
                ydatafit = numpy.asarray(xdata)*(numpy.asarray(pos)/1000)**3
                edatafit = numpy.asarray([numpy.fmax(a,b)*(c/1000)**3 for a,b,c in zip(errorsleft,errorsright,pos)])
                result = scipy.optimize.curve_fit(lambda x, m, b: m * x + b, xdata=xdatafit, ydata=ydatafit)[0]  # p0=(0.001, -0.1)
                m = result[0]
                b = result[1]
                xline = numpy.linspace(0, 500, 5)
                ax.errorbar(xline, m * xline + b, fmt=':', lw=0.5, color='gray', alpha=0.5)
    elif obs == "chisqdof":
        for i in range(nfiles):
            if not numpy.isnan(xdata[i]).any():
                ax.errorbar(pos[i], xdata[i], yerr=[[errorsleft[i]], [errorsright[i]]], color=colors[i], fmt='x-', fillstyle='none', markersize=5,
                            mew=0.25, lw=0.8, elinewidth=0.5, capsize=1.2, zorder=-10)
    elif obs == "corr" or obs == "spf" or obs == "spf_phys":
        fmts = [val + "-" for val in lpd.markers]
        if obs == "corr":
            for i in range(nfiles):
                if not numpy.isnan(xdata[i]).any():
                    plots.append(ax.errorbar(xdata[i]+i*0.003, ydata[i], yerr=[errorsleft[i], errorsright[i]], fmt=fmts[i], fillstyle='none', markersize=2, mew=0.25, lw=0.2, elinewidth=0.2, capsize=1.2, zorder=-10))
            ax.axhline(y=1, **lpd.horizontallinestyle)
        if obs == "spf":
            if PhiUV_file is not None:
                PhiUV = numpy.loadtxt(PhiUV_file)
                PhiUV = PhiUV[:, 0:2]
                labels = (r'$\Phi^\mathrm{UV}$', *labels)
                plots.append(ax.errorbar(PhiUV[:, 0], 2*PhiUV[:, 1]/PhiUV[:, 0], fmt='-', lw=0.4)) # /PhiUV[:, 0]
            ax.set_yscale('log')
            ax.set_xscale('log')
            for i in range(nfiles):
                if not numpy.isnan(xdata[i]).any():
                    if plot_spf_err:
                        yerr = [errorsleft[i], errorsright[i]]
                    else:
                        yerr = None

                    # plots.append(ax.errorbar(xdata[i] / Ntaus[i], 2*ydata[i]/xdata[i] / Ntaus[i]**2,  yerr=yerr, fmt='--', fillstyle='none', markersize=0, mew=0.25, lw=0.4, elinewidth=0.2, capsize=1.2, zorder=-10))  # yerr=[errorsleft[i], errorsright[i]],
                    plots.append(
                        ax.errorbar(xdata[i], 2 * ydata[i] / xdata[i], yerr=yerr, fmt='--', fillstyle='none', markersize=0, mew=0.25,
                                    lw=0.4, elinewidth=0.2, capsize=1.2, zorder=-10))  # yerr=[errorsleft[i], errorsright[i]],
                    # /xdata[i]
            ax.axvline(x=1, **lpd.verticallinestyle)
            ax.axvline(x=omegaUV, **lpd.verticallinestyle)
            # ax.axvline(x=2.817735677972127, **lpd.verticallinestyle)
            # ax.axhline(y=4*numpy.pi / 0.9, **lpd.horizontallinestyle)
        # ax.axhline(y=1, **lpd.horizontallinestyle)
        if obs == "spf_phys":
            # ax.axvline(x=1, **lpd.verticallinestyle)
            # ax.axvline(x=omegaUV, **lpd.verticallinestyle)
            for i in range(nfiles):
                if not numpy.isnan(xdata[i]).any():
                    if plot_spf_err:
                        yerr = [errorsleft[i], errorsright[i]]
                    else:
                        yerr = None
                    plots.append(ax.errorbar(xdata[i]*(pos[i]/1000), ydata[i]*(pos[i]/1000)**3,  yerr=yerr, fmt='--', fillstyle='none', markersize=0, mew=0.25, lw=0.4, elinewidth=0.2, capsize=1.2, zorder=-10))  # yerr=[errorsleft[i], errorsright[i]],
        ax.legend(handles=plots, labels=labels, title=legtitle, title_fontsize=6, **lpd.chmap(lpd.legendstyle, fontsize=6))
    else:
        print("ERROR: UNKNOWN OBS")
        exit(1)

def main():
    parser = argparse.ArgumentParser()
    requiredNamed = parser.add_argument_group('required named arguments')

    parser.add_argument('--outputpath', help='the path of output folder like /a/b', type=str, required=True)
    parser.add_argument('--suffix', help='append this to the output file name', type=str, default="")

    parser.add_argument('--file_basepath', help="prepend this to each argument passed via --files and to --PhiUV_file", type=str, default="")
    requiredNamed.add_argument('--files', help="files to plot", nargs='*', type=str, required=True)
    parser.add_argument('--PhiUV_file', help='path to PhiUV file if you want to include it in the spf plot', type=str, default=None)

    requiredNamed.add_argument('--labels', help="labels for the plots", nargs='*', type=str, required=True)
    requiredNamed.add_argument('--obs', help='whether to plot kappa, chisqdof, corr or spf', choices=["kappa", "chisqdof", "corr", "spf", "spf_phys"], required=True)

    parser.add_argument('--title', help='title of plot', type=str)

    parser.add_argument('--colors', help="colors in order for each file", default=('k', 'C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9'), nargs='*', type=str)
    parser.add_argument('--pos', help="y positions of kappas. only used when obs=kappa", nargs='*', type=float)

    parser.add_argument('--wideaspect', action="store_true", help="use a wide aspect ratio (basically make the plots larger)")
    parser.add_argument('--xlims', help='custom xlims for plot', nargs=2, type=float, default=None)
    parser.add_argument('--ylims', help='custom ylims for plot', nargs=2, type=float, default=None)
    parser.add_argument('--xlabel', type=str, default=None)
    parser.add_argument('--ylabel', type=str, default=None)
    parser.add_argument('--plot_spf_err', help='usually one cant see anything when plotting these errors, but sometimes it may be helpful.', action="store_true")
    parser.add_argument('--Nts', type=int, help="List of Nt of the lattices. necessary for spf plot.", nargs='*')
    parser.add_argument('--kappa_swap_axes', action="store_true", help="swap axes of kappa plot")
    parser.add_argument('--OmegaByT_UV', type=float, default=-1, help="a dashed vertical line is drawn at this point for the spf plots.")
    parser.add_argument('--legtitle', default="", help="title of legend in spf plots.")
    parser.add_argument('--verbose', action="store_true", help="log debug output")

    args = parser.parse_args()

    for i,  file in enumerate(args.files):
        args.files[i] = args.file_basepath + "/" + file
    nfiles = len(list(args.files))
    if args.file_basepath and args.PhiUV_file:
        args.PhiUV_file = args.file_basepath + "/" + args.PhiUV_file

    filelabel, xlabel, ylabel, xlims, ylims, ylabelpos = get_plot_settings(args.obs)
    xdata, ydata, errorsleft, errorsright = load_data(args.files, args.obs, args.verbose)

    # setup figure
    figsize = (3 + 3 / 8, 3 + 3 / 8 - 1 / 2.54)
    if args.wideaspect:
        figsize = (1.5 * (3 + 3 / 8), 1.5 * (3 + 3 / 8) / 16 * 9)
    if args.xlims:
        xlims = args.xlims
    if args.ylims:
        ylims = args.ylims
    if args.xlabel:
        xlabel = args.xlabel
    if args.ylabel:
        ylabel = args.ylabel
    fig, ax, plots = lpd.create_figure(figsize=figsize, xlabel=xlabel, ylabel=ylabel, ylabelpos=ylabelpos, xlims=xlims, xlabelpos=(0.5, -0.1), ylims=ylims)

    plot(args.obs, nfiles, xdata, ydata, errorsleft, errorsright, ax, plots, args.labels, args.PhiUV_file, args.plot_spf_err, args.pos, args.colors, args.Nts, args.kappa_swap_axes, args.OmegaByT_UV, args.legtitle)

    if args.title is not None:
        ax.set_title(args.title, fontsize=8)

    # save figure
    lpd.create_folder(args.outputpath)
    outfile = args.outputpath+"/"+args.obs+"_"+args.suffix+".pdf"
    fig.savefig(outfile)
    print("saved ", outfile)


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
