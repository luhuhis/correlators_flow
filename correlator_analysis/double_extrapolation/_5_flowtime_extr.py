#!/usr/bin/env python3
import numpy
import matplotlib
import lib_process_data as lpd
import scipy.optimize
import scipy.interpolate


import warnings
warnings.filterwarnings('ignore', r'.*All-NaN slice encountered.*')


def extrapolation_ansatz(x, m, b):
    return m * x + b


def chisq(params, args):
    xdata = args[0]
    ydata = args[1]
    if args[2] is not None:
        edata = args[2]
    else:
        edata = 1
    if args[3] is not None:
        corr_inv = args[3]
    else:
        corr_inv = numpy.identity(len(ydata))

    r = (ydata - extrapolation_ansatz(xdata, params[0], params[1]))/edata
    return r.T @ corr_inv @ r


def plot_extrapolation(args, xdata, ydata, edata, ydata_extr, edata_extr, indices, plotbasepath):

    finest_Nt_half = int(args.finest_Nt/2)
    tauTs = lpd.get_tauTs(args.finest_Nt)

    mintauTindex = None

    # setup figure
    displaystyle = '' if not args.use_tex else r'\displaystyle'
    ylims = (2.5, 3.8) if not args.custom_ylims else args.custom_ylims
    ylabel_prefix = ""
    if args.Zf2_file is not None:
        ylabel_prefix = r'Z_f^2'
    fig, ax, _ = lpd.create_figure(ylims=ylims, ylabel=r'$' + displaystyle + r'\frac{'+ylabel_prefix+r'G'+lpd.get_corr_subscript(args.corr)+r'}{G^\text{norm}}$',
                                       xlabel=r'$'+displaystyle+r'{8\tau_\mathrm{F}}/{\tau^2}$')
                                       # r'\frac{Z_f^2 G_B}{G^\text{norm}}$', UseTex=args.use_tex)
    ax.set_xlim([-0.005, args.max_FlowradiusBytauT**2 * 1.15])
    plots = []

    # fill figure with data
    for i in range(len(ydata)):
        tauT = tauTs[i]
        if tauT >= args.min_tauT_plot:

            if mintauTindex is None:
                mintauTindex = i

            if args.relflow_file:
                xdataplot = xdata**2
            else:
                xdataplot = xdata*8/tauT**2
            mycolor = lpd.get_color(tauTs, i, mintauTindex, finest_Nt_half - 1)

            if not args.no_extr and not numpy.isnan(ydata_extr[i][1]):
                # plot fit result at tf=0
                plots.append(ax.errorbar(0, ydata_extr[i][1], edata_extr[i][1], fmt='|', color=mycolor, zorder=1, label='{0:.3f}'.format(tauT)))

                # plot linear fit line
                x = numpy.asarray([0, xdataplot[-1]])
                ax.errorbar(x, extrapolation_ansatz(x, *ydata_extr[i]), color=mycolor, **lpd.fitlinestyle, zorder=-10000)

                # plot data points used in extrapolation
                ax.errorbar(xdataplot[indices[i]], ydata[i][indices[i]], edata[i][indices[i]], fmt='|', color=mycolor)

            # plot data error band
            if args.no_extr:
                label = '{0:.3f}'.format(tauT)
                alpha = 1
            else:
                label = None
                alpha = 0.5
            ax.fill_between(xdataplot, ydata[i]-edata[i], ydata[i]+edata[i], facecolor=mycolor, alpha=alpha, zorder=-300+i, label=label)

    # second x-axis for flow radius
    ax2 = ax.twiny()
    new_tick_locations = numpy.array([args.min_FlowradiusBytauT**2, args.max_FlowradiusBytauT**2])
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(new_tick_locations)
    ax2.set_xticklabels(["%.1f" % z if z == 0 else "%.2f" % z for z in numpy.sqrt(new_tick_locations)])
    ax2.set_xlabel(r'$'+displaystyle+r'{\sqrt{8\tau_\mathrm{F}}}/{\tau}$', horizontalalignment='right', verticalalignment='top', zorder=999999, bbox=lpd.labelboxstyle)
    ax2.xaxis.set_label_coords(0.97, 0.97)
    ax2.tick_params(direction='in', pad=0, width=0.5)

    # legend
    ax.legend(handles=plots)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1], title=r'$\tau T$', loc="center left", bbox_to_anchor=(1.01, 0.5), **lpd.leg_err_size(1, 0.2))

    # save plot
    file = plotbasepath+"/"+args.corr+"_flow_extr_quality"+args.output_suffix+".pdf"
    print("saving", file)
    fig.savefig(file)
    matplotlib.pyplot.close(fig)


def plot_corr(args, xdata, ydata, edata, plotbasepath):
    # plot final double-extrapolated correlator in its own plot
    ylims = [1.5, 4] if not args.custom_ylims else args.custom_ylims
    fig, ax, plots = lpd.create_figure(xlims=[0, 0.51], ylims=ylims, xlabel=r'$\tau T$',
                                       ylabel=r'$\frac{G'+lpd.get_corr_subscript(args.corr)+r'}'
                                              r'{ G^\text{norm}}$',
                                       UseTex=args.use_tex)
    ax.errorbar(xdata, ydata, edata, color='black', markersize=0)
    file = plotbasepath+"/"+args.corr+"_flow_extr"+args.output_suffix+".pdf"
    print("saving", file)
    fig.savefig(file)
    matplotlib.pyplot.close(fig)


def do_flow_extr(index, tauTs, cont_samples, data_std, n_samples, args, flowsteps):
    """ do the extrapolation for all samples of one tauT"""

    # make a flag for rel flows and change flowtimes / flowradii in this function accordingly. maybe make a general flow variable instead of calling flow"times" and flow"radii".

    tauT = tauTs[index]

    results = numpy.empty((n_samples, 2))
    results[:] = numpy.nan

    if args.relflow_file:
        relflows = flowsteps
        flowend = numpy.fabs(relflows - args.max_FlowradiusBytauT).argmin()
        flowstart = numpy.fabs(relflows - args.min_FlowradiusBytauT).argmin()
        midpoint = numpy.abs(relflows - (args.max_FlowradiusBytauT+args.min_FlowradiusBytauT)/2).argmin()
        indices = numpy.unique([flowstart, midpoint, flowend])
        xdatatmp = relflows[indices]**2
    else:
        flowtimes = flowsteps
        flowradii = numpy.sqrt(8*flowtimes)
        flowend = numpy.abs(flowradii - lpd.upper_flowradius_limit_(tauT, args.max_FlowradiusBytauT, 0)).argmin()
        flowstart = numpy.abs(flowradii - lpd.upper_flowradius_limit_(tauT, args.min_FlowradiusBytauT, 0)).argmin()
        geom_mean = (numpy.abs(flowtimes - numpy.sqrt(flowtimes[flowstart] * flowtimes[flowend]))).argmin()
        indices = numpy.unique([flowstart, geom_mean, flowend])
        xdatatmp = flowtimes[indices] * 8 / tauT ** 2

    if len(indices) >= 3:

        for n in range(n_samples):

            ydatatmp = cont_samples[n][index][indices]
            edatatmp = data_std[index][indices]
            mask = numpy.isnan(ydatatmp)
            xdata = xdatatmp[~mask]
            ydata = ydatatmp[~mask]
            edata = edatatmp[~mask]

            # TODO undo this again when working on corr matrix again
            # this_cov = cov[i][index[:, numpy.newaxis], index]

            if len(xdata) >= 3:  # minimum amount for linear fit

                if not args.no_extr:
                    fitparams = scipy.optimize.minimize(chisq, x0=numpy.asarray([0, 3]), bounds=(args.slope_bounds, (None, None)), args=[xdata, ydata, edata, None])
                    fitparams = fitparams.x

                    results[n] = fitparams
    print("done ", '{0:.3f}'.format(tauT))  # TODO , ", taufT2=[", lpd.format_float(numpy.sqrt(8*flowtimes[flowstart])), ", ", lpd.format_float(numpy.sqrt(8*flowtimes[flowend])), "]", sep="")
    return results, indices


def load_data(args, basepath):  # TODO flowtimes
    if not args.relflow_file:
        samples = numpy.load(basepath + "/cont_extr/" + args.corr + "_cont_samples.npy")
        cont_samples = samples[:, :, :, 0]
    else:
        samples = numpy.load(basepath + "/cont_extr/" + args.corr + "_cont_relflow_samples.npy")
        nt_finest_half = int(args.finest_Nt/2)
        cont_samples = samples[:, :, :nt_finest_half]
    n_samples = len(cont_samples)
    # TODO for BB, also take care of Zf2 in relflow units
    # if args.Zf2_file is not None:
    #     tfT2, Zf2 = numpy.loadtxt(args.Zf2_file, unpack=True)
    #     Zf2_int = scipy.interpolate.InterpolatedUnivariateSpline(tfT2, Zf2, k=3, ext=1)
    #     for j in range(len(flowtimes)):
    #         if Zf2_int(flowtimes[j]) == 0:
    #             print("warn", flowtimes[j])
    #         cont_samples[:, j, :] *= Zf2_int(flowtimes[j])
    cont_samples = numpy.swapaxes(cont_samples, 1, 2)
    data_std = lpd.dev_by_dist(cont_samples, axis=0)
    # shape: (1000, 18, 221)

    return cont_samples, data_std, n_samples


def parse_args():
    # parse cmd line arguments
    parser, requiredNamed = lpd.get_parser()

    parser.add_argument('--flowtimes_finest', help='file with list of all dimensionless flowtimes of the finest lattice', type=str)
    parser.add_argument('--finest_Nt', help='finest Nt used in previous cont extr.', type=int, required=True)
    parser.add_argument('--use_tex', action="store_true", default=False)
    parser.add_argument('--max_FlowradiusBytauT', type=float, default=numpy.sqrt(8 * 0.014),
                        help='modify the flow filter based on tauT to be more/less strict. default value of 0.33 means that for each flow radius the tauT '
                             'cannot be greater than 3*flowradius, or that for each tauT the flow radius must be less than 0.33*tauT.')
    parser.add_argument('--min_FlowradiusBytauT', type=float, default=0.2)
    parser.add_argument('--min_tauT_plot', default=0, type=float)
    parser.add_argument('--no_extr', help='do NOT perform a flow-time-to-zero extrapolation, just plot the data', action="store_true")
    parser.add_argument('--custom_ylims', help="custom y-axis limits for both plots", type=float, nargs=2)
    parser.add_argument('--temp_subfolder', help="suffix for input path, to incorporate subfolders for different temperatures", default="", type=str)
    parser.add_argument('--basepath', help='where to look for files', type=str)
    parser.add_argument('--basepath_plot', help='where to save plots', type=str)
    parser.add_argument('--Zf2_file', help="Z_f^2 for BB correlator renormalization", default=None, type=str)
    parser.add_argument('--output_suffix', help="add this string to the end of plot file names", default="", type=str)
    parser.add_argument('--slope_bounds', help="bound the slope fit parameter between these values", default=(None, None), nargs=2, type=float)
    parser.add_argument('--nproc', type=int, default=20)
    parser.add_argument('--relflow_file', help="if provided, this indicates that input samples are already in relative flow units, and this is then the path to the relative flow times file.")
    parser.add_argument('--combined_fit', action="store_true")
    parser.add_argument('--n_samples', type=int)
    args = parser.parse_args()

    # check if given arguments are valid
    if args.relflow_file and args.flowtimes_finest:
        parser.error("use either --relflow_file or --flowtimes_finest")

    return args


def independent_extrapolation(args, ntauT, finest_tauTs, cont_samples, data_std, n_samples, flowsteps):
    fitparams, indices = lpd.parallel_function_eval(do_flow_extr, range(0, ntauT), args.nproc, finest_tauTs, cont_samples, data_std, n_samples, args, flowsteps)

    results = numpy.swapaxes(numpy.asarray(fitparams), 0, 1)

    return results, indices


def plot_combined_extrapolation(args, xdata, ydata, edata, ydata_extr, edata_extr, indices, plotbasepath):

    finest_Nt_half = int(args.finest_Nt/2)
    tauTs = lpd.get_tauTs(args.finest_Nt)

    mintauTindex = None

    # setup figure
    displaystyle = '' if not args.use_tex else r'\displaystyle'
    ylims = (2.5, 3.8) if not args.custom_ylims else args.custom_ylims
    ylabel_prefix = ""
    if args.Zf2_file is not None:
        ylabel_prefix = r'Z_f^2'
    fig, ax, _ = lpd.create_figure(ylims=ylims, ylabel=r'$' + displaystyle + r'\frac{'+ylabel_prefix+r'G'+lpd.get_corr_subscript(args.corr)+r'}{G^\text{norm}}$',
                                       xlabel=r'$'+displaystyle+r'{8\tau_\mathrm{F}}/{\tau^2}$')
                                       # r'\frac{Z_f^2 G_B}{G^\text{norm}}$', UseTex=args.use_tex)
    ax.set_xlim([-0.005, args.max_FlowradiusBytauT**2 * 1.15])
    plots = []

    # fill figure with data
    for i in range(len(ydata)):
        tauT = tauTs[i]
        if tauT >= args.min_tauT_plot:

            if mintauTindex is None:
                mintauTindex = i

            if args.relflow_file:
                xdataplot = xdata**2
            else:
                xdataplot = xdata*8/tauT**2
            mycolor = lpd.get_color(tauTs, i, mintauTindex, finest_Nt_half - 1)

            print(ydata_extr[i], edata_extr[i])

            if not args.no_extr and not numpy.isnan(ydata_extr[i]):
                # plot fit result at tf=0
                plots.append(ax.errorbar(0, ydata_extr[i], edata_extr[i], fmt='|', color=mycolor, zorder=1, label='{0:.3f}'.format(tauT)))

                # plot linear fit line
                x = numpy.asarray([0, xdataplot[-1]])
                ax.errorbar(x, combined_fit_ansatz(x, tauT, ydata_extr[i], *ydata_extr[len(ydata):-1]), color=mycolor, **lpd.fitlinestyle, zorder=-10000)

                # plot data points used in extrapolation
                ax.errorbar(xdataplot[indices], ydata[i][indices], edata[i][indices], fmt='|', color=mycolor)

            # plot data error band
            if args.no_extr:
                label = '{0:.3f}'.format(tauT)
                alpha = 1
            else:
                label = None
                alpha = 0.5
            ax.fill_between(xdataplot[indices[0]:indices[-1]+1], ydata[i][indices[0]:indices[-1]+1]-edata[i][indices[0]:indices[-1]+1], ydata[i][indices[0]:indices[-1]+1]+edata[i][indices[0]:indices[-1]+1], facecolor=mycolor, alpha=alpha, zorder=-300+i, label=label)

    # second x-axis for flow radius
    ax2 = ax.twiny()
    new_tick_locations = numpy.array([args.min_FlowradiusBytauT**2, args.max_FlowradiusBytauT**2])
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(new_tick_locations)
    ax2.set_xticklabels(["%.1f" % z if z == 0 else "%.2f" % z for z in numpy.sqrt(new_tick_locations)])
    ax2.set_xlabel(r'$'+displaystyle+r'{\sqrt{8\tau_\mathrm{F}}}/{\tau}$', horizontalalignment='right', verticalalignment='top', zorder=999999, bbox=lpd.labelboxstyle)
    ax2.xaxis.set_label_coords(0.97, 0.97)
    ax2.tick_params(direction='in', pad=0, width=0.5)

    # legend
    ax.legend(handles=plots)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1], title=r'$\tau T$', loc="center left", bbox_to_anchor=(1.01, 0.5), **lpd.leg_err_size(1, 0.2))

    # save plot
    file = plotbasepath+"/"+args.corr+"_flow_extr_quality"+args.output_suffix+".pdf"
    print("saving", file)
    fig.savefig(file)
    matplotlib.pyplot.close(fig)


n_additional_fitparams = 1
def combined_fit_ansatz(x, tauT, cont, base_slope):
    return cont + (base_slope)*x


def combined_chisqdof(fitparams, ydata, xdata, edata, tauTs):
    ndata = ydata.size
    ntauT = ydata.shape[0]
    chisq = 0

    combined_fitparams = fitparams[ntauT:]

    for i in range(ntauT):
        chisq += numpy.nansum(((combined_fit_ansatz(xdata, tauTs[i], fitparams[i], *combined_fitparams) - ydata[i])/edata[i]) ** 2)

    nfitparams = len(fitparams)
    chisqdof = chisq / (ndata-nfitparams)
    return chisqdof


def perform_combined_fit(ydata, xdata, tauTs, edata, nparams):
    # fitparams = scipy.optimize.minimize(chisq, x0=numpy.asarray([0, 3]), bounds=(args.slope_bounds, (None, None)), args=[xdata, ydata, edata, None])
    # fitparams = fitparams.x
    #
    # results[n] = fitparams

    ntauT = len(tauTs)
    # print(xdata, ydata, edata)
    fitparams = scipy.optimize.minimize(combined_chisqdof, x0=numpy.asarray([*[1 for _ in range(ntauT)], *[1 for _ in range(nparams)]]),
                                        bounds=([*[(0, 20) for _ in range(ntauT)], (None, 0), *[(None, None) for _ in range(nparams-2)]]), args=(ydata, xdata, edata, tauTs))
    fitparams = fitparams.x
    chisqdof = combined_chisqdof(fitparams, ydata, xdata, edata, tauTs)
    return [*fitparams, chisqdof]


def count_falses_from_start(arr):
    count = 0
    for value in arr:
        if not value:
            count += 1
        else:
            break
    return count


def combined_extrapolation(args, ntauT, finest_tauTs, cont_samples, data_std, n_samples, flowsteps):

    """ do the extrapolation for all samples of one tauT"""

    # make a flag for rel flows and change flowtimes / flowradii in this function accordingly. maybe make a general flow variable instead of calling flow"times" and flow"radii".
    results = numpy.empty((n_samples, ntauT+n_additional_fitparams+1))
    results[:] = numpy.nan

    relflows = flowsteps
    flowend = numpy.fabs(relflows - args.max_FlowradiusBytauT).argmin()
    flowstart = numpy.fabs(relflows - args.min_FlowradiusBytauT).argmin()
    midpoint = numpy.abs(relflows - (args.max_FlowradiusBytauT+args.min_FlowradiusBytauT)/2).argmin()
    indices = numpy.unique([flowstart, midpoint, flowend])
    xdata = relflows[indices]**2

    if len(indices) >= 3:

        edata = numpy.asarray([data_std[i][indices] for i in range(ntauT)])
        mask = numpy.isnan(edata).any(axis=1)  # this should have shape tauT
        edata = edata[~mask]
        offset = count_falses_from_start(~mask)
        for n in range(args.n_samples):
            ydatatmp = []
            for t in range(ntauT):
                ydatatmp.append(cont_samples[n][t][indices])
            ydata = numpy.asarray(ydatatmp)[~mask]  # ydata has now length of valid tauTs, and contains three flow points. same as edata should have

            if len(xdata) >= 3:  # minimum amount for linear fit
                fitresults = perform_combined_fit(ydata, xdata, finest_tauTs[~mask], edata, n_additional_fitparams)
                for i, val in enumerate(fitresults):
                    results[n][offset + i] = val
            if n % 100 == 0:
                print(n)

    return results, indices


def plot_final_corr_and_savetxt(args, results_mean, results_std, ntauT, finest_tauTs, plotbasepath, basepath, suffix):
    # plot and save final correlator

    if args.combined_fit:
        ydata_extr = results_mean[:ntauT]
        edata_extr = results_std[:ntauT]
    else:
        ydata_extr = results_mean[:, 1]
        edata_extr = results_std[:, 1]
    plot_corr(args, finest_tauTs, ydata_extr, edata_extr, plotbasepath)
    correlator = numpy.column_stack((finest_tauTs, ydata_extr, edata_extr))
    print(correlator)
    file = basepath + "/" + args.corr + "_flow_extr"+suffix+".txt"
    print("saving mean and error in ", file)
    numpy.savetxt(file, correlator, header="                tauT                G/Gnorm                    err", fmt='%22.15e')


def save_extr_samples(args, basepath, suffix, results):
    file = basepath + "/" + args.corr + "_flow_extr"+suffix+".npy"
    print("saving all samples in", file)
    numpy.save(file, results)
    print("results shape:", results.shape)


def plot_flow_extr(args, flowsteps, cont_samples, data_std, results_mean, results_std, indices, plotbasepath):

    data_mean = numpy.nanmedian(cont_samples, axis=0)  # original data points
    if args.combined_fit:
        plot_combined_extrapolation(args, flowsteps, data_mean, data_std, results_mean, results_std, indices, plotbasepath)
    else:
        plot_extrapolation(args, flowsteps, data_mean, data_std, results_mean, results_std, indices, plotbasepath)


def main():

    # TODO how are the flowtimes actually used in this script?

    args = parse_args()

    basepath = lpd.get_merged_data_path(args.qcdtype, args.corr, "", args.basepath) + args.temp_subfolder

    # parse some more arguments
    if not args.relflow_file:
        flowsteps = numpy.loadtxt(args.flowtimes_finest)/args.finest_Nt**2  # flowsteps = flowtimes T^2
    else:
        flowsteps = numpy.loadtxt(args.relflow_file)  # flowsteps = relative flow radii

    finest_tauTs = lpd.get_tauTs(args.finest_Nt)
    ntauT = len(finest_tauTs)

    cont_samples, data_std, n_samples = load_data(args, basepath)  # TODO flowtimes for BB

    if args.n_samples is None:
        args.nsamples = cont_samples.shape[0]

    if args.combined_fit:
        results, indices = combined_extrapolation(args, ntauT, finest_tauTs, cont_samples, data_std, n_samples, flowsteps)
    else:
        results, indices = independent_extrapolation(args, ntauT, finest_tauTs, cont_samples, data_std, n_samples, flowsteps)

    suffix = ""
    if args.relflow_file:
        suffix = "_relflow"

    save_extr_samples(args, basepath, suffix, results)

    # plots
    results_mean = numpy.nanmedian(results[:args.n_samples], axis=0)
    results_std = lpd.dev_by_dist(results[:args.n_samples], axis=0)
    print(results_mean[-1])
    if not args.combined_fit:
        print(numpy.column_stack((results_mean[:,0], results_std[:,0])))
    plotbasepath = lpd.get_plot_path(args.qcdtype, args.corr, "", args.basepath_plot) + args.temp_subfolder
    plot_flow_extr(args, flowsteps, cont_samples, data_std, results_mean, results_std, indices, plotbasepath)
    plot_final_corr_and_savetxt(args, results_mean, results_std, ntauT, finest_tauTs, plotbasepath, basepath, suffix)


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()