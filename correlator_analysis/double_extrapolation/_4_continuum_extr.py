#!/usr/bin/env python3

import lib_process_data as lpd
import numpy
import matplotlib
import scipy.optimize
from matplotlib.backends.backend_pdf import PdfPages
import warnings

warnings.filterwarnings('ignore', r'All-NaN (slice|axis) encountered.*')
warnings.filterwarnings('ignore', r'Calling figure.constrained_layout, but figure not setup to do constrained layout.*')
warnings.filterwarnings('ignore', r'constrained_layout not applied.*')
warnings.filterwarnings('ignore', r'.*converting a masked element to nan.*')
warnings.filterwarnings('ignore', r'.*Covariance of the parameters could not be estimated.*')
matplotlib.rcParams['figure.max_open_warning'] = 0  # suppress warning due to possibly large number of figures...
numpy.set_printoptions(linewidth=numpy.inf, precision=3, threshold=numpy.inf, floatmode="fixed")


def main():
    args = parse_args()
    Nts, nt_finest, nt_coarsest = parse_nts(args)

    samples = load_data(args)
    edatas = get_weights(samples)

    if not args.relflow:
        traditional_extr(args, samples, edatas, nt_finest, Nts)
    else:
        new_extr(args, samples, edatas, nt_finest, Nts)


def parse_args():
    # parse cmd line arguments
    parser, requiredNamed = lpd.get_parser()

    requiredNamed.add_argument('--conftypes', nargs='*', required=True,
                               help="ORDERED list of conftypes (from coarse to fine), e.g. s080t20_b0703500 s096t24_b0719200 s120t30_b0739400 s144t36_b0754400")

    parser.add_argument('--use_tex', action="store_true", help="use LaTeX when plotting")
    parser.add_argument('--custom_ylims', help="custom y-axis limits for both plots", type=float, nargs=2)
    parser.add_argument('--max_FlowradiusBytauT', type=float, default=0.3001)
    parser.add_argument('--min_FlowradiusBytauT', type=float, default=0.2499)
    parser.add_argument('--output_suffix', default="", help="append this to the output folder name")
    parser.add_argument('--verbose', help='print out progress information', action="store_true")
    parser.add_argument('--basepath', type=str, help="where to look for the data")
    parser.add_argument('--basepath_plot', type=str, help="where to save the plots")
    parser.add_argument('--min_flowradius', type=float, help="minimum flowradius for the extrapolation. default=1/min(Nts)", default=None)
    parser.add_argument('--nproc', type=int, default=20, help="number of processes for parallelization")
    parser.add_argument('--flow_index_range', default=(0, -1), type=int, nargs=2,
                        help="which flow indices to consider (default considers all). useful if youre only interested in some speficic ones.")
    parser.add_argument('--ansatz', type=str, choices=["linear", "constant", "custom"], default="linear")
    parser.add_argument('--relflow', action="store_true")
    parser.add_argument('--combined_fit', action="store_true")
    parser.add_argument('--nsamples', type=int, default=None)
    parser.add_argument('--nterms', help='how many terms to add in the slope of the combined fit', default=1, choices=[1, 2], type=int)

    args = parser.parse_args()

    global n_additional_fitparams
    n_additional_fitparams = args.nterms
    global combined_fit_ansatz
    if args.nterms == 1:
        combined_fit_ansatz = combined_fit_ansatz_1
    elif args.nterms == 2:
        combined_fit_ansatz = combined_fit_ansatz_2

    return args


def parse_nts(args):
    # parse Nts
    Nts = []
    for conftype in args.conftypes:
        _, _, tmp_nt, _ = lpd.parse_conftype(conftype)
        Nts.append(tmp_nt)
    Nts = numpy.asarray(Nts)
    nt_finest = numpy.max(Nts)
    nt_coarsest = numpy.min(Nts)

    return Nts, nt_finest, nt_coarsest


def load_data(args):
    print("load data...")
    samples = []

    suffix = ""
    if args.relflow:
        suffix = "relflow_"

    for conftype in args.conftypes:
        path = lpd.get_merged_data_path(args.qcdtype, args.corr, conftype,
                                        args.basepath) + "/" + args.corr + "_" + conftype + "_interpolation_" + suffix + "samples.npy"
        tmp = numpy.load(path)
        samples.append(tmp)

    samples = numpy.stack([conftype_samples for conftype_samples in samples], 0)

    if args.nsamples is None:
        args.nsamples = samples.shape[2]

    print("Done. Data layout: (n_conftypes, n_flowtimes, n_samples, Nt/2): ", samples.shape)

    return samples


def get_weights(samples):
    # we use the sample deviation as weights for each sample fit
    edatas = []
    for sample_set in samples:
        tmp = lpd.dev_by_dist(sample_set, axis=1)
        edatas.append(tmp)
    edatas = numpy.asarray(edatas)
    return edatas


def linear_ansatz(x, b, m, chisq_dof=None):
    return m * x + b


def fit_sample(ydata, xdata, edata):
    def chisq_dof(fitparams, ydata, xdata, edata, extrapolation_ansatz):
        ndata = len(ydata)
        chisq = 0
        for i in range(ndata):
            chisq += (extrapolation_ansatz(xdata[i], *fitparams) - ydata[i]) ** 2 / edata[i] ** 2

        nparams = 2
        chisqdof = chisq / (ndata - nparams)
        return chisqdof

    start_params = [3, 0]

    if not numpy.isnan(ydata).any():
        fitparams, _ = scipy.optimize.curve_fit(linear_ansatz, xdata, ydata, p0=start_params, sigma=edata)
        chisqdof = chisq_dof(fitparams, ydata, xdata, edata, linear_ansatz)
        return [*fitparams, chisqdof]
    else:
        return [numpy.nan] * (len(start_params) + 1)


def traditional_extr(args, samples, edatas, nt_finest, Nts):
    def load_flowtimes(args):
        # load flow times from finest lattice
        flowtimes = numpy.loadtxt(
            lpd.get_merged_data_path(args.qcdtype, args.corr, args.conftypes[-1], args.basepath) + "/flowtimes_" + args.conftypes[-1] + ".dat")
        if args.flow_index_range == (0, -1):
            indices = range(0, len(flowtimes))
        else:
            indices = range(*args.flow_index_range)

        return flowtimes, indices

    print("calculate extrapolation and create figures...")
    print("      ", lpd.get_tauTs(nt_finest))
    flowtimes, indices = load_flowtimes(args)
    fitparams, figs_corr, figs_extr = lpd.parallel_function_eval(wrapper, indices, args.nproc, flowtimes, samples, args, Nts, edatas)

    save_figs(args, figs_corr, "_cont")
    save_figs(args, figs_extr, "_cont_quality")
    save_data(args, fitparams, flowtimes, nt_finest)

    return


def wrapper(index, flowtimes, samples, args, Nts, edatas):
    nt_finest = numpy.max(Nts)
    nt_coarsest = numpy.min(Nts)

    flowtime = flowtimes[index]
    flowradius = numpy.sqrt(flowtime * 8) / nt_finest
    flowradius_str = r'{0:.4f}'.format(flowradius)

    # define some parameters
    lower_tauT_lim = flowradius / args.max_FlowradiusBytauT
    if args.min_flowradius is None:
        args.min_flowradius = 1 / nt_coarsest
    upper_tauT_lim = flowradius / args.min_FlowradiusBytauT
    tauTs_fine = lpd.get_tauTs(nt_finest)
    valid_tauTs = [tauT for tauT in tauTs_fine if lower_tauT_lim <= tauT <= upper_tauT_lim]
    n_valid_tauTs = len(valid_tauTs)
    nt_half_fine = int(nt_finest / 2)
    offset = nt_half_fine - n_valid_tauTs

    # declarations
    nsamples = samples[0].shape[1]

    results = numpy.empty((nsamples, nt_half_fine, 3))
    results[:] = numpy.nan

    # at fixed flowtime, perform cont extr for each tauT for each sample
    if flowradius >= args.min_flowradius:
        for j, tauT in enumerate(tauTs_fine):
            xdata = numpy.asarray([1 / Ntau ** 2 for k, Ntau in enumerate(Nts)])

            if tauT in valid_tauTs and len(xdata) >= 2:

                for m in range(nsamples):
                    ydata = [conftype_samples[index][m, j] for conftype_samples in samples]
                    edata = [conftype_edata[index, j] for conftype_edata in edatas]

                    results[m][j] = fit_sample(ydata, xdata, edata)

    # TODO else: load data

    # assemble data and calculate means in order to produce some intermediate plots

    # mean and std of the original data samples
    data_mean = numpy.nanmedian(samples, axis=2)[:, index, :]
    data_std = lpd.dev_by_dist(samples, axis=2)[:, index, :]

    # mean and std of the fit parameters
    results_mean = numpy.nanmedian(results, axis=0)
    results_std = lpd.dev_by_dist(results, axis=0)

    continuum_mean = results_mean[:, 0]
    continuum_std = results_std[:, 0]
    chisqdof_mean = results_mean[:, -1]

    if not numpy.isnan(chisqdof_mean).all():
        print(flowradius_str, chisqdof_mean)

    # merge continuum and data together into one object to pass it to the plot function
    plot_data = numpy.transpose(numpy.stack([continuum_mean, *[tmp for tmp in data_mean]]))
    plot_std = numpy.transpose(numpy.stack([continuum_std, *[tmp for tmp in data_std]]))

    fig_corr = plot_corr(args, lpd.get_tauTs(nt_finest), continuum_mean, continuum_std, flowradius_str, lpd.lower_tauT_limit_(flowradius, 0.3))
    fig_extr = plot_extrapolation(args, Nts, plot_data, plot_std, results_mean, flowradius_str, tauTs_fine,
                                  [tauT for tauT in tauTs_fine if flowradius / 0.301 <= tauT <= flowradius / 0.249],
                                  offset)  # TODO remove hard code for the plot

    return results, fig_corr, fig_extr


def plot_extrapolation(args, xdata, ydata, edata, fitparams, flowradius_str, tauTs_fine, valid_tauTs, offset):
    nt_coarse = numpy.min(xdata)

    xdata = numpy.asarray([0, *[1 / tmp ** 2 for tmp in xdata]])

    ylims = (2.55, 3.75) if not args.custom_ylims else args.custom_ylims

    ylabel = r'$\displaystyle \frac{G' + lpd.get_corr_subscript(args.corr) + r'}{G^\mathrm{norm}}$'
    fig, ax, _ = lpd.create_figure(xlims=[-0.0001, 1 / nt_coarse ** 2 * 1.05], ylims=ylims, xlabel=r'$N_\tau^{-2}$', ylabel=ylabel)
    plots = []
    ax.text(0.99, 0.99, r'$ \sqrt{8\tau_\mathrm{F}}T = ' + flowradius_str + "$", ha='right', va='top', transform=ax.transAxes, zorder=-1000,
            bbox=lpd.labelboxstyle)

    maxtauTindex_plot = 0
    mintauTindex = None

    for j in range(len(ydata)):
        if tauTs_fine[j] in valid_tauTs:
            # print(xdata.shape, ydata[j].shape, edata[j].shape)

            # plot extrapolation
            if mintauTindex is None:
                mintauTindex = j
            if j > maxtauTindex_plot:
                maxtauTindex_plot = j
            mycolor = lpd.get_color(valid_tauTs, len(valid_tauTs) - 1 - (j - mintauTindex))
            plots.append(ax.errorbar(xdata, ydata[j], edata[j], fmt='|', color=mycolor, zorder=-100 + j,
                                     label='{0:.3f}'.format(tauTs_fine[j])))  # lpd.markers[j - nt_half_fine + len(lpd.markers)]
            x = numpy.linspace(0, 0.1, 100)
            ax.errorbar(x, linear_ansatz(x, *fitparams[j]), color=mycolor, **lpd.fitlinestyle, zorder=-100)

    ax.legend(handles=plots)
    handles, labels = ax.get_legend_handles_labels()
    # reverse ordering of legend
    ax.legend(handles[::-1], labels[::-1], title=r'$\tau T$', loc="center left", bbox_to_anchor=(1, 0.5), **lpd.leg_err_size(1, 0.25), ncol=1,
              columnspacing=0.5)

    return fig


def plot_corr(args, xdata, ydata, edata, flowradius_str, lower_tauT_lim):
    # save plot of single continuum extrapolation for this flow time
    if args.use_tex:
        ylabel = r'$\displaystyle \frac{G' + lpd.get_corr_subscript(args.corr) + r'}{G^\mathrm{norm}}$'
    else:
        ylabel = r'$\frac{G}{G^\mathrm{norm}}$'
    ylims = (1.4, 4) if not args.custom_ylims else args.custom_ylims
    fig, ax, plots = lpd.create_figure(xlims=[0, 0.55], ylims=ylims, xlabel=r'$\tau T$',
                                       ylabel=ylabel, UseTex=args.use_tex, constrained_layout=True)
    ax.text(0.99, 0.99, r'$ \sqrt{8\tau_\mathrm{F}}T = ' + flowradius_str + "$", ha='right', va='top', transform=ax.transAxes, zorder=-1000)
    ax.axvline(x=lower_tauT_lim, **lpd.verticallinestyle)
    ax.set_xticks((0.0, 0.1, 0.2, 0.3, 0.4, 0.5))
    # filter out nan's to suppress warnings in ax.errorbar
    mask = ~numpy.isnan(ydata)
    xdata = xdata[mask]
    ydata = ydata[mask]
    edata = edata[mask]
    ax.errorbar(xdata, ydata, edata, color="black", markersize=0, fmt='x')
    return fig


def save_figs(args, figs, suffix):
    # save figures to pdf
    print("save figures...")
    lpd.set_rc_params()  # for some reason we need to repeat this here...
    plotpath = lpd.get_plot_path(args.qcdtype, args.corr, args.output_suffix, args.basepath_plot)
    lpd.create_folder(plotpath)
    with PdfPages(plotpath + "/" + args.corr + suffix + ".pdf") as pdf:
        for fig in figs:
            if fig is not None:
                pdf.savefig(fig)
                matplotlib.pyplot.close(fig)


def save_data(args, fitparams, flowtimes, nt_finest):
    # save fitparams to npy file
    print("save extrapolation data...")
    fitparams = numpy.asarray(fitparams)
    fitparams = numpy.swapaxes(fitparams, 0, 1)

    print("fitparams shape:", fitparams.shape)
    folder = lpd.get_merged_data_path(args.qcdtype, args.corr, args.output_suffix, args.basepath) + "/cont_extr/"
    lpd.create_folder(folder)
    numpy.save(folder + args.corr + "_cont_samples.npy", fitparams)

    with open(folder + args.corr + "_cont.dat", 'w') as outfile:
        outfile.write('# bootstrap mean of continuum ' + args.corr + ' correlator for ' + args.qcdtype + '\n')
        outfile.write('# rows correspond to flow times, columns to dt = {1, ... , Ntau/2} of the finest lattice that entered the extrapolation\n')
        cont_mean = numpy.mean(fitparams, axis=0)[:, :, 0]
        numpy.savetxt(outfile, cont_mean)
    with open(folder + args.corr + "_cont_err.dat", 'w') as outfile:
        outfile.write('# bootstrap err of mean of continuum ' + args.corr + ' correlator for ' + args.qcdtype + '\n')
        outfile.write('# rows correspond to flow times, columns to dt = {1, ... , Ntau/2} of the finest lattice that entered the extrapolation\n')
        cont_std = numpy.std(fitparams, axis=0)[:, :, 0]
        numpy.savetxt(outfile, cont_std)
    numpy.savetxt(folder + args.corr + "_cont_flowtimesT2.dat", flowtimes / nt_finest ** 2, header='# tau_F T^2')


# ==============================================================================================================================================================


combined_fit_ansatz = None
n_additional_fitparams = None


def combined_fit_ansatz_1(x, tauT, cont, a):
    return cont + (- (a / tauT) ** 2) * x


def combined_fit_ansatz_2(x, tauT, cont, a, b):
    return cont + (- (a / tauT) ** 2 - (b / tauT) ** 4) * x


def perform_combined_fit(ydata, xdata, tauTs, edata, nparams):
    def combined_chisqdof(fitparams, ydata, xdata, edata, tauTs):
        ndata = ydata.size
        ntauT = ydata.shape[1]
        chisq = 0

        combined_fitparams = fitparams[ntauT:]

        for i in range(ntauT):
            chisq += numpy.nansum(((combined_fit_ansatz(xdata, tauTs[i], fitparams[i], *combined_fitparams) - ydata[:, i]) / edata[:, i]) ** 2)

        nfitparams = len(fitparams)
        chisqdof = chisq / (ndata - nfitparams)
        return chisqdof

    ntauT = len(tauTs)
    fitparams = scipy.optimize.minimize(combined_chisqdof, x0=numpy.asarray([*[1 for _ in range(ntauT)], *[1 for _ in range(nparams)]]),
                                        bounds=([*[(0, 20) for _ in range(ntauT)], (None, None), *[(None, None) for _ in range(nparams - 1)]]),
                                        args=(ydata, xdata, edata, tauTs))
    fitparams = fitparams.x
    chisqdof = combined_chisqdof(fitparams, ydata, xdata, edata, tauTs)
    return [*fitparams, chisqdof]


def new_extr(args, samples, edatas, nt_finest, Nts):
    merged_data_path = lpd.get_merged_data_path(args.qcdtype, args.corr, args.conftypes[-1], args.basepath)
    relflows = numpy.loadtxt(merged_data_path + args.corr + "_" + args.conftypes[-1] + "_relflows.txt")
    nflow = len(relflows)

    nfitparams = int(nt_finest / 2) + n_additional_fitparams

    if args.combined_fit:
        results = lpd.parallel_function_eval(combined_extr_at_relflow, range(nflow), args.nproc, args, samples, edatas, nfitparams, nt_finest, Nts)
        results = numpy.asarray(results)
        results = results.swapaxes(0, 1)
    else:
        results = lpd.parallel_function_eval(individual_extr_at_relflow, range(nflow), args.nproc, args, samples, edatas, nt_finest, Nts)
        results = numpy.asarray(results)
        results = results.swapaxes(0, 1)

    print("results shape", results.shape)
    # samples[:, :, :nt_finest_half]

    # save data
    save_relflow_data(args, results, relflows, nt_finest)

    # plot extr
    xdata = 1 / numpy.asarray(Nts) ** 2
    figs = lpd.parallel_function_eval(plot_relflow_extr, range(nflow), args.nproc, args, relflows, results, xdata, samples, edatas, lpd.get_tauTs(nt_finest))
    save_figs(args, figs, "_cont_quality_relflow")


def combined_extr_at_relflow(index, args, samples, edatas, nfitparams, nt_finest, Nts):
    tauTs_finest = lpd.get_tauTs(nt_finest)
    nt_finest_half = int(nt_finest / 2)

    results = numpy.empty((args.nsamples, nfitparams + 1))
    results[:] = numpy.nan

    xdata = numpy.asarray([1 / Ntau ** 2 for Ntau in Nts])

    valid_tauT_indices = ~numpy.isnan(edatas[0, index])
    valid_tauT_indices = valid_tauT_indices & [*[False for _ in range(int(nt_finest_half / 2) - 1)], *[True for _ in range(int(nt_finest_half / 2) + 1)]]
    offset = count_falses_from_start(valid_tauT_indices)

    edata = numpy.asarray([conftype_edata[index][valid_tauT_indices] for conftype_edata in edatas])

    for m in range(args.nsamples):
        ydata = numpy.asarray([conftype_samples[index][m][valid_tauT_indices] for conftype_samples in samples])
        fitresults = perform_combined_fit(ydata, xdata, tauTs_finest[valid_tauT_indices], edata, n_additional_fitparams)
        for i, val in enumerate(fitresults):
            results[m][offset + i] = val

    print("done", index)

    return results


def individual_extr_at_relflow(index, args, samples, edatas, nt_finest, Nts):
    nt_finest_half = int(nt_finest / 2)

    # we sort the fit results like this: first the nt/2 y-intercepts, then the nt/2 slopes, then the chisqdofs.
    # this is such that the flow extrapolation finds the intercepts in the first part, as it expects it from the combined fit as well.
    results = numpy.empty((args.nsamples, int(nt_finest_half * 3)))
    results[:] = numpy.nan

    xdata = numpy.asarray([1 / Ntau ** 2 for Ntau in Nts])

    valid_tauT_indices = ~numpy.isnan(edatas[0, index])
    valid_tauT_indices = valid_tauT_indices & [*[False for _ in range(int(nt_finest_half / 2) - 1)], *[True for _ in range(int(nt_finest_half / 2) + 1)]]
    offset = count_falses_from_start(valid_tauT_indices)

    for j in range(offset, nt_finest_half):

        for m in range(args.nsamples):
            ydata = [conftype_samples[index][m, j] for conftype_samples in samples]
            edata = [conftype_edata[index, j] for conftype_edata in edatas]

            fitparams = fit_sample(ydata, xdata, edata)

            results[m][j] = fitparams[0]
            results[m][j + nt_finest_half] = fitparams[1]
            results[m][j + int(2 * nt_finest_half)] = fitparams[2]

    print("done", index)

    return results


def count_falses_from_start(arr):
    count = 0
    for value in arr:
        if not value:
            count += 1
        else:
            break
    return count


def plot_relflow_extr(index, args, relflows, results, xdata, samples, edatas, tauTs):
    relflow = relflows[index]

    if 0.2 <= relflow <= 0.33:

        orig_ydata = numpy.asarray([numpy.nanmedian(conftype_samples[index], axis=0) for conftype_samples in samples]).swapaxes(0, 1)
        orig_edata = numpy.asarray([conftype_edata[index] for conftype_edata in edatas]).swapaxes(0, 1)

        xlims = (-0.0002, numpy.amax(xdata) * 1.05)
        fig, ax, _ = lpd.create_figure(xlims=xlims, ylims=args.custom_ylims, xlabel=r'$1/N_\tau^2$',
                                       ylabel=r'$\displaystyle \frac{G' + lpd.get_corr_subscript(args.corr) + r'}{G^\mathrm{norm}}$')

        fitparams = numpy.nanmedian(results, axis=0)[index]
        if 0.25 <= relflows[index] <= 0.3:
            print(lpd.format_float(relflows[index]), fitparams[:len(tauTs)], fitparams[len(tauTs):])
        fitparams_err = lpd.dev_by_dist(results, axis=0)[index]
        xpoints = numpy.linspace(0, numpy.amax(xdata), 5)
        plots = []

        mintauTindex = None

        for i, tauT in enumerate(tauTs):
            if numpy.count_nonzero(~numpy.isnan(orig_ydata[i])) == len(xdata) and tauT >= 0.249:
                if mintauTindex is None:
                    mintauTindex = i
                mycolor = lpd.get_color(tauTs, len(tauTs) - 1 - (i - mintauTindex), mintauTindex)
                plots.append(
                    ax.errorbar(xdata, orig_ydata[i], orig_edata[i], fmt='|', label=lpd.format_float(tauT), color=mycolor))
                # print(lpd.format_float(tauT), "-> chisqdof=", fitparams[-1])

                if args.combined_fit:
                    ax.errorbar(xpoints, combined_fit_ansatz(xpoints, tauT, fitparams[i], *fitparams[len(tauTs):-1]), **lpd.fitlinestyle, color=mycolor)
                else:
                    ax.errorbar(xpoints, linear_ansatz(xpoints, fitparams[i], fitparams[i + len(tauTs)]), **lpd.fitlinestyle, color=mycolor)
                ax.errorbar(0, fitparams[i], fitparams_err[i], fmt='|', color=mycolor)

        ax.text(0.99, 0.99, r'$\sqrt{8\tau_\mathrm{F}}/\tau=' + lpd.format_float(relflows[index]) + r'$', ha='right', va='top', transform=ax.transAxes)
        ax.legend(handles=plots)
        handles, labels = ax.get_legend_handles_labels()
        # reverse ordering of legend
        ax.legend(handles[::-1], labels[::-1], title=r'$\tau T$', loc="center left", bbox_to_anchor=(1, 0.5), **lpd.leg_err_size(1, 0.25), ncol=1,
                  columnspacing=0.5)
        return fig
    else:
        return None


def save_relflow_data(args, fitparams, relflows, nt_finest):
    # we need samples flowtime tauT as the shape

    folder = lpd.get_merged_data_path(args.qcdtype, args.corr, args.output_suffix, args.basepath) + "/cont_extr/"
    lpd.create_folder(folder)
    numpy.save(folder + args.corr + "_cont_relflow_samples.npy", fitparams)

    nt_finest_half = int(nt_finest / 2)
    with open(folder + args.corr + "_cont_relflow.dat", 'w') as outfile:
        outfile.write('# bootstrap mean of continuum ' + args.corr + ' correlator for ' + args.qcdtype + '\n')
        outfile.write('# rows correspond to flow times, columns to dt = {1, ... , Ntau/2} of the finest lattice that entered the extrapolation\n')
        cont_mean = numpy.nanmedian(fitparams, axis=0)[:, :nt_finest_half]
        numpy.savetxt(outfile, cont_mean)
    with open(folder + args.corr + "_cont_relflow_err.dat", 'w') as outfile:
        outfile.write('# bootstrap err of mean of continuum ' + args.corr + ' correlator for ' + args.qcdtype + '\n')
        outfile.write('# rows correspond to flow times, columns to dt = {1, ... , Ntau/2} of the finest lattice that entered the extrapolation\n')
        cont_std = lpd.dev_by_dist(fitparams, axis=0)[:, :nt_finest_half]
        numpy.savetxt(outfile, cont_std)
    numpy.savetxt(folder + args.corr + "_cont_relflows.dat", relflows, header='# sqrt(8tau_F)/tau')


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
