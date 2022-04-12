#!/usr/local/bin/python3.7m
import lib_process_data as lpd
import numpy
import argparse


def load_data(files, obs: str):
    # === load data
    xdata = []
    ydata = []
    errorsleft = []
    errorsright = []
    for file in files:
        try:
            data = numpy.loadtxt(file)
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
        elif obs == "spf":
            ydata.append(data[:, 1])
            xdata.append(data[:, 0])
            errorsleft.append(data[:, 2])
            errorsright.append(data[:, 3])
    return xdata, ydata, errorsleft, errorsright


def get_plot_settings(obs: str):
    ylabel = None
    ylims = None
    xlims = None
    xlabel = None
    ylabelpos = None
    filelabel = ""

    if obs == "kappa" or obs == "chisqdof":
        if obs == "kappa":
            xlims = (1.5, 12.5)
            obslabel = r'\kappa/ T^3'
            filelabel = "params_"
            ylims = (0, 5)
        if obs == "chisqdof":
            xlims = (0, 9)
            obslabel = r'\chi^2/\mathrm{d.o.f.}'
            filelabel = "params_"
            ylims = (0, 5)
        xlabel = r'$\mathrm{median}(' + obslabel + r') \pm 34^{\mathrm{th}}\, \%$'
        ylabelpos = None

    if obs == "corr" or obs == "spf":
        ylabelpos = (0.05, 0.95)
        if obs == "corr":
            filelabel = "corrfit_"
            # ylabel = r'$\frac{\mathrm{median}(G^\mathrm{model}) \pm 34^{\mathrm{th}}\, \%}{G^\mathrm{cont} \pm \delta G}$'
            ylabel = r'$\frac{G^\mathrm{model}}{G^\mathrm{cont}}$'
            xlabel = r'$\tau T$'
            xlims = (0.2, 0.52)
            ylims = (0.95, 1.05)
        if obs == "spf":
            obslabel = r'\frac{\rho}{\omega T^2}'
            xlabel = r'$\omega/T$'
            xlims = (0.1, 100)
            ylims = (1, 100)
            filelabel = "spffit_"
            # ylabel = r'$\mathrm{median}(' + obslabel + r') \pm 34^{\mathrm{th}}\, \%$'
            ylabel = r'$' + obslabel + r'$'

    return filelabel, xlabel, ylabel, xlims, ylims, ylabelpos


def plot(obs, nfiles, xdata, ydata, errorsleft, errorsright, ax, plots, labels, PhiUV_file, plot_spf_err, pos, colors):
    if obs == "kappa" or obs == "chisqdof":

        for i in range(nfiles):
            if not numpy.isnan(xdata[i]).any():
                ax.errorbar(xdata[i], pos[i], xerr=[[errorsleft[i]], [errorsright[i]]], color=colors[i], fmt='x-', fillstyle='none', markersize=5, mew=0.25, lw=0.8, elinewidth=0.5, capsize=1.2, zorder=-10)
        # custom yaxis
        ax.set_yticks(pos)
        ax.set_yticklabels(labels)
        # vertical dashed lines
        currentxlims = ax.get_xlim()
        for i in range(0, int(currentxlims[1])):
            ax.axvline(i, **(lpd.chmap(lpd.verticallinestyle, alpha=0.4)))

    fmts = [val + "-" for val in lpd.markers]

    if obs == "corr" or obs == "spf":
        if obs == "corr":
            for i in range(nfiles):
                if not numpy.isnan(xdata[i]).any():
                    plots.append(ax.errorbar(xdata[i]+i*0.003, ydata[i], yerr=[errorsleft[i], errorsright[i]], fmt=fmts[i], fillstyle='none', markersize=2, mew=0.25, lw=0.2, elinewidth=0.2, capsize=1.2, zorder=-10))
        if obs == "spf":
            if PhiUV_file is not None:
                PhiUV = numpy.loadtxt(PhiUV_file)
                PhiUV = PhiUV[:, 0:2]
                labels = (r'$\Phi^\mathrm{UV}$', *labels)
                plots.append(ax.errorbar(PhiUV[:, 0], PhiUV[:, 1]/PhiUV[:, 0], fmt='-', lw=0.4))
            ax.set_yscale('log')
            ax.set_xscale('log')
            for i in range(nfiles):
                if not numpy.isnan(xdata[i]).any():
                    if plot_spf_err:
                        yerr = [errorsleft[i], errorsright[i]]
                    else:
                        yerr = None
                    plots.append(ax.errorbar(xdata[i]+i*0.003, ydata[i]/xdata[i],  yerr=yerr, fmt=fmts[i], fillstyle='none', markersize=0, mew=0.25, lw=0.2, elinewidth=0.2, capsize=1.2, zorder=-10))  # yerr=[errorsleft[i], errorsright[i]],
        ax.legend(handles=plots, labels=labels, **lpd.legendstyle)
        ax.axhline(y=1, **lpd.horizontallinestyle)


def main():
    parser = argparse.ArgumentParser()
    requiredNamed = parser.add_argument_group('required named arguments')

    parser.add_argument('--outputpath', help='the path of output folder like /a/b', type=str, required=True)

    parser.add_argument('--file_basepath', help="prepend this to each argument passed via --files", type=str, default="")
    requiredNamed.add_argument('--files', help="files to plot", nargs='*', type=str, required=True)
    parser.add_argument('--PhiUV_file', help='path to PhiUV file if you want to include it in the spf plot', type=str, default=None)

    requiredNamed.add_argument('--labels', help="labels for the plots", nargs='*', type=str, required=True)
    requiredNamed.add_argument('--obs', help='whether to plot kappa, chisqdof, corr or spf', choices=["kappa", "chisqdof", "corr", "spf"], required=True)
    parser.add_argument('--suffix', help='append this to the output file name', type=str, default="")
    parser.add_argument('--title', help='title of plot', type=str)

    parser.add_argument('--colors', help="colors in order for each file", default=('k', 'C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9'), nargs='*', type=str)
    parser.add_argument('--pos', help="y positions of kappas. only used when obs=kappa", default=0.3 + numpy.asarray((0.45, 0.85, 1.25, 1.65, 2.65, 3.05, 3.45, 3.85)), nargs='*')

    parser.add_argument('--wideaspect', action="store_true", help="use a wide aspect ratio (basically make the plots larger)")
    parser.add_argument('--xlims', help='custom xlims for plot', nargs=2, type=float, default=None)
    parser.add_argument('--ylims', help='custom ylims for plot', nargs=2, type=float, default=None)
    parser.add_argument('--plot_spf_err', help='usually one cant see anything when plotting these errors, but sometimes it may be helpful.', action="store_true")

    args = parser.parse_args()

    for i,  file in enumerate(args.files):
        args.files[i] = args.file_basepath + "/" + file
    nfiles = len(list(args.files))

    filelabel, xlabel, ylabel, xlims, ylims, ylabelpos = get_plot_settings(args.obs)
    xdata, ydata, errorsleft, errorsright = load_data(args.files, args.obs)

    # setup figure
    figsize = (3 + 3 / 8, 3 + 3 / 8 - 1 / 2.54)
    if args.wideaspect:
        figsize = (1.5 * (3 + 3 / 8), 1.5 * (3 + 3 / 8) / 16 * 9)
    if args.xlims:
        xlims = args.xlims
    if args.ylims:
        ylims = args.ylims
    fig, ax, plots = lpd.create_figure(figsize=figsize, xlabel=xlabel, ylabel=ylabel, ylabelpos=ylabelpos, xlims=xlims, xlabelpos=(0.5, -0.1), ylims=ylims)

    plot(args.obs, nfiles, xdata, ydata, errorsleft, errorsright, ax, plots, args.labels, args.PhiUV_file, args.plot_spf_err, args.pos, args.colors)

    if args.title is not None:
        ax.set_title(args.title)

    # save figure
    lpd.create_folder(args.outputpath)
    outfile = args.outputpath+"/"+args.obs+"_"+args.suffix+".pdf"
    fig.savefig(outfile)
    print("saved ", outfile)


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
