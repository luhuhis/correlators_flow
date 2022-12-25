#!/usr/bin/env python3

import lib_process_data as lpd
import numpy
import matplotlib
import scipy.optimize
from matplotlib.backends.backend_pdf import PdfPages


# hide warning for fit when only having 2 data points.
import warnings
# warnings.simplefilter("ignore", OptimizeWarning)
# warnings.simplefilter('error', UserWarning)
warnings.filterwarnings('ignore', r'All-NaN (slice|axis) encountered')
warnings.filterwarnings('ignore', r'Calling figure.constrained_layout, but figure not setup to do constrained layout.*')
warnings.filterwarnings('ignore', r'constrained_layout not applied*')


def extrapolation_ansatz(x, m, b):
    return m * x + b


def fit_sample(ydata, xdata, edata, start_params=None):
    if start_params is None:
        start_params = [0, 3]
    fitparams, _ = scipy.optimize.curve_fit(extrapolation_ansatz, xdata, ydata, p0=start_params, sigma=edata)
    return fitparams


def wrapper(index, flowtimes, samples, args, Nts, edatas):

    nt_finest = numpy.max(Nts)
    nt_coarsest = numpy.min(Nts)

    flowtime = flowtimes[index]
    flowradius = numpy.sqrt(flowtime*8)/nt_finest
    flowradius_str = r'{0:.4f}'.format(flowradius)

    # define some parameters
    lower_tauT_lim = lpd.lower_tauT_limit_(flowradius, args.max_FlowradiusBytauT, args.max_FlowradiusBytauT_offset)
    tauTs_fine = lpd.get_tauTs(nt_finest)
    valid_tauTs = [tauT for tauT in tauTs_fine if tauT > lower_tauT_lim]
    n_valid_tauTs = len(valid_tauTs)
    nt_half_fine = int(nt_finest/2)
    offset = nt_half_fine-n_valid_tauTs

    # declarations
    nsamples = samples[0].shape[1]
    # TODO if compute:

    results = numpy.empty((nsamples, nt_half_fine, 2))
    results[:] = numpy.nan

    # at fixed flowtime, perform cont extr for each tauT for each sample
    if flowradius >= args.min_flowradius:
        for j, tauT in enumerate(tauTs_fine):
            xdata = numpy.asarray([1/Ntau**2 for k, Ntau in enumerate(Nts)])

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

    # split the fit parameters into slope and intercept
    slope_mean, continuum_mean = numpy.squeeze(numpy.split(results_mean, 2, axis=1))
    slope_std, continuum_std = numpy.squeeze(numpy.split(results_std, 2, axis=1))

    # merge continuum and data together into one object to pass it to the plot function
    plot_data = numpy.transpose(numpy.stack([continuum_mean, *[tmp for tmp in data_mean]]))
    plot_std = numpy.transpose(numpy.stack([continuum_std, *[tmp for tmp in data_std]]))

    fig_corr = plot_corr(args, lpd.get_tauTs(nt_finest), continuum_mean, continuum_std, flowradius_str, lower_tauT_lim)
    fig_extr = plot_extrapolation(args, Nts, plot_data, plot_std, results_mean, flowradius_str, tauTs_fine, valid_tauTs, offset)

    return results, fig_corr, fig_extr


def plot_extrapolation(args, xdata, ydata, edata, fitparams, flowradius_str, tauTs_fine, valid_tauTs, offset):

    nt_half_fine = int(numpy.max(xdata)/2)
    nt_coarse = numpy.min(xdata)

    xdata = numpy.asarray([0, *[1/tmp**2 for tmp in xdata]])

    ylims = (2.55, 3.75) if not args.custom_ylims else args.custom_ylims
    fig, ax, plots = lpd.create_figure(xlims=[-0.0001, 1/nt_coarse**2*1.05], ylims=ylims, xlabel=r'$N_\tau^{-2}$', ylabel=r'$\displaystyle\frac{G}{G^\mathrm{norm}}$')
    ax.text(0.99, 0.99, r'$ \sqrt{8\tau_\mathrm{F}}T = ' + flowradius_str + "$", ha='right', va='top', transform=ax.transAxes, zorder=-1000, bbox=lpd.labelboxstyle)

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
            mycolor = lpd.get_color(tauTs_fine, nt_half_fine - 1 - j + offset, mintauTindex, nt_half_fine - 1)
            plots.append(ax.errorbar(xdata, ydata[j], edata[j], fmt='|', color=mycolor, zorder=-100 + j, label='{0:.3f}'.format(tauTs_fine[j])))  # lpd.markers[j - nt_half_fine + len(lpd.markers)]
            x = numpy.linspace(0, 0.1, 100)
            ax.errorbar(x, extrapolation_ansatz(x, *fitparams[j]), color=mycolor, fmt=':', zorder=-100)

    # numpy.savetxt(lpd.get_merged_data_path(args.qcdtype, args.corr, args.output_suffix)+"/cont_extr/"+args.corr+"_"+flowradius_str+"_cont.txt", results,
    #               header="tauT    G/G_norm    err")

    # save continuum extrapolation quality plot for this flow time

    ax.legend(handles=plots)
    handles, labels = ax.get_legend_handles_labels()
    # reverse ordering of legend
    ax.legend(handles[::-1], labels[::-1], title=r'$\tau T$', loc="lower left", bbox_to_anchor=(-0.01, -0.01), **lpd.leg_err_size(1, 0.25), ncol=4,
              columnspacing=0.5)

    return fig


def plot_corr(args, xdata, ydata, edata, flowradius_str, lower_tauT_lim):
    # save plot of single continuum extrapolation for this flow time
    if args.use_tex:
        displaystyle = r'\displaystyle'
        ylabel = r'$'+displaystyle+r'\frac{G}{G^\mathrm{norm}}$'
    else:
        ylabel = r'$\frac{G}{G^\mathrm{norm}}$'
    ylims = (1.4, 4) if not args.custom_ylims else args.custom_ylims
    fig, ax, plots = lpd.create_figure(xlims=[0, 0.55], ylims=ylims, xlabel=r'$\tau T$',
                                       ylabel=ylabel, UseTex=args.use_tex, constrained_layout=True)
    ax.text(0.99, 0.99, r'$ \sqrt{8\tau_\mathrm{F}}T = '+flowradius_str+"$", ha='right', va='top', transform=ax.transAxes, zorder=-1000)
    ax.axvline(x=lower_tauT_lim, **lpd.verticallinestyle)
    ax.set_xticks((0.0, 0.1, 0.2, 0.3, 0.4, 0.5))
    # filter out nan's to suppress warnings in ax.errorbar
    mask = ~numpy.isnan(ydata)
    xdata = xdata[mask]
    ydata = ydata[mask]
    edata = edata[mask]
    ax.errorbar(xdata, ydata, edata, color="black", markersize=0, fmt='x')
    return fig


def main():

    # parse cmd line arguments
    parser, requiredNamed = lpd.get_parser()

    requiredNamed.add_argument('--conftypes', nargs='*', required=True,
                               help="ORDERED list of conftypes (from coarse to fine), e.g. s080t20_b0703500 s096t24_b0719200 s120t30_b0739400 s144t36_b0754400")

    parser.add_argument('--use_tex', action="store_true", help="use LaTeX when plotting")
    parser.add_argument('--custom_ylims', help="custom y-axis limits for both plots", type=float, nargs=2)
    parser.add_argument('--max_FlowradiusBytauT', type=float, default=numpy.sqrt(8*0.014),
                        help='modify the tauT filter based on flow time to be more/less strict. default value of 0.33 means that for each tauT the flow radius '
                             'cannot be greater than 0.33*tauT, or that for each flow radius the tauT must be atleast 3*flowradius.')
    parser.add_argument('--max_FlowradiusBytauT_offset', type=float, default=1/20,
                        help='fixed offset to make lower_tauT_limit stricter (by e.g. one lattice spacing 1/Nt), as the 0.33 criterion is only valid in the '
                             'continuum. on the lattice one has to be stricter. 1/Nt_coarsest is a good value.')
    parser.add_argument('--output_suffix', default="", help="append this to the output folder name")
    parser.add_argument('--verbose', help='print out progress information', action="store_true")
    parser.add_argument('--basepath', type=str, help="where to look for the data")
    parser.add_argument('--basepath_plot', type=str, help="where to save the plots")
    parser.add_argument('--min_flowradius', type=float, help="minimum flowradius for the extrapolation. default=1/min(Nts)", default=None)
    parser.add_argument('--nproc', type=int, default=20, help="number of processes for parallelization")
    parser.add_argument('--flow_index_range', default=(0,-1), type=int, nargs=2,
                        help="which flow indices to consider (default considers all). useful if youre only interested in some speficic ones.")

    args = parser.parse_args()

    # load flow times from finest lattice
    flowtimes = numpy.loadtxt(lpd.get_merged_data_path(args.qcdtype, args.corr, args.conftypes[-1], args.basepath)+"/flowtimes_"+args.conftypes[-1]+".dat")
    if args.flow_index_range == (0,-1):
        indices = range(0, len(flowtimes))
    else:
        indices = range(*args.flow_index_range)

    _, _, _, gaugeaction, flowaction = lpd.parse_qcdtype(args.qcdtype)

    # parse Nts
    Nts = []
    for conftype in args.conftypes:
        _, _, tmp_nt, _ = lpd.parse_conftype(conftype)
        Nts.append(tmp_nt)
    Nts = numpy.asarray(Nts)
    nt_finest = numpy.max(Nts)
    nt_coarsest = numpy.min(Nts)

    if args.min_flowradius is None:
        args.min_flowradius = 1/nt_coarsest

    # load data
    print("load data...")
    samples = []
    for conftype in args.conftypes:
        if conftype != args.conftypes[-1]:
            path = lpd.get_merged_data_path(args.qcdtype, args.corr, conftype, args.basepath) + "/" + args.corr + "_" + conftype + "_interpolation_samples.npy"
            tmp = numpy.load(path)
            samples.append(tmp)
        else:
            path = lpd.get_merged_data_path(args.qcdtype, args.corr, conftype, args.basepath) + "/" + args.corr + "_" + conftype + "_samples.npy"
            tmp = numpy.load(path)
            # add tree-level imp for finest lattice
            for j in range(tmp.shape[1]):
                for k in range(tmp.shape[2]):
                    tmp[:, j, k] *= nt_finest ** 4 / lpd.G_latt_LO_flow(k, flowtimes[j], args.corr, nt_finest, flowaction, gaugeaction)
            tmp = numpy.swapaxes(tmp, 0, 1)
            samples.append(tmp)

    samples = numpy.stack([conftype_samples for conftype_samples in samples], 0)
    print("Done. Data layout: (n_conftypes, n_flowtimes, n_samples, Nt/2): ", samples.shape)

    # we use the sample deviation as weights for each sample fit
    edatas = []
    for sample_set in samples:
        tmp = lpd.dev_by_dist(sample_set, axis=1)
        edatas.append(tmp)
    edatas = numpy.asarray(edatas)

    matplotlib.rcParams['figure.max_open_warning'] = 0  # suppress warning due to possibly large number of figures...

    print("calculate extrapolation and create figures...")

    # fitparams, figs_corr, figs_extr = lpd.serial_function_eval(wrapper, indices, flowtimes, samples, args, Nts, edatas)
    fitparams, figs_corr, figs_extr = lpd.parallel_function_eval(wrapper, indices, args.nproc, flowtimes, samples, args, Nts, edatas)

    # save figures to pdf
    print("save figures...")
    lpd.set_rc_params()  # for some reason we need to repeat this here...
    plotpath = lpd.get_plot_path(args.qcdtype, args.corr, args.output_suffix, args.basepath_plot)
    lpd.create_folder(plotpath)
    with PdfPages(plotpath + "/" + args.corr + "_cont.pdf") as pdf:
        for fig in figs_corr:
            pdf.savefig(fig)
            matplotlib.pyplot.close(fig)
    with PdfPages(plotpath + "/" + args.corr + "_cont_quality.pdf") as pdf:
        for fig in figs_extr:
            pdf.savefig(fig)
            matplotlib.pyplot.close(fig)

    # save fitparams to npy file
    print("save extrapolation data...")
    fitparams = numpy.asarray(fitparams)
    fitparams = numpy.swapaxes(fitparams, 0, 1)

    print("fitparams shape:", fitparams.shape)
    folder = lpd.get_merged_data_path(args.qcdtype, args.corr, args.output_suffix, args.basepath) + "/cont_extr/"
    lpd.create_folder(folder)
    numpy.save(folder + args.corr + "_cont_samples.npy", fitparams)

    with open(folder + args.corr + "_cont.dat", 'w') as outfile:
        outfile.write('# bootstrap mean of continuum '+args.corr+' correlator for '+args.qcdtype+'\n')
        outfile.write('# rows correspond to flow times, columns to dt = {1, ... , Ntau/2} of the finest lattice that entered the extrapolation\n')
        cont_mean = numpy.mean(fitparams, axis=0)[:, :, 1]
        numpy.savetxt(outfile, cont_mean)
    with open(folder + args.corr + "_cont_err.dat", 'w') as outfile:
        outfile.write('# bootstrap err of mean of continuum '+args.corr+' correlator for '+args.qcdtype+'\n')
        outfile.write('# rows correspond to flow times, columns to dt = {1, ... , Ntau/2} of the finest lattice that entered the extrapolation\n')
        cont_std = numpy.std(fitparams, axis=0)[:, :, 1]
        numpy.savetxt(outfile, cont_std)
    numpy.savetxt(folder+args.corr+"_cont_flowtimesT2.dat", flowtimes/nt_finest**2, header='# tau_F T^2')


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
