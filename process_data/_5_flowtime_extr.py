#!/usr/bin/env python3
import numpy
import matplotlib
from matplotlib import container, legend_handler
import lib_process_data as lpd
import scipy.optimize


def extrapolation_ansatz(x, m, b):
    return m * x + b


def extr_jac(x, m, b):
    return x, 1


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

    density_weight = []
    for i in range(len(xdata)):
        if i == len(xdata)-1:
            density_weight.append((xdata[i]-xdata[i-1]))
        elif i == 0:
            density_weight.append((xdata[i + 1]-xdata[i]))
        else:
            density_weight.append((xdata[i + 1] - xdata[i - 1]) / 2)

    density_weight = numpy.asarray(density_weight)
    r = (ydata - extrapolation_ansatz(xdata, params[0], params[1]))/edata
    return r.T @ corr_inv @ (r * density_weight)


def do_fit(ydata, xdata, edata, this_cov_inv):
    fitparams = scipy.optimize.minimize(chisq, x0=numpy.asarray([-65, 3]), args=[xdata, ydata, edata, this_cov_inv])
    return fitparams.x


def fit_sample(ydata, xdata, edata):
    fitparams, _ = scipy.optimize.curve_fit(extrapolation_ansatz, xdata, ydata, p0=[-64, 3], sigma=edata, maxfev=2000, method='lm', ftol=1e-12, xtol=1e-12, gtol=1e-12)
    return fitparams


def plot_extrapolation(xdata, ydata, edata, ydata_extr, edata_extr, args, mintauTindex, plotbasepath):

    finest_Nt_half = int(args.finest_Nt/2)
    tauTs = lpd.get_tauTs(args.finest_Nt)

    # setup figure
    displaystyle = '' if not args.use_tex else r'\displaystyle'
    ylims = (2.5, 3.8) if not args.custom_ylims else args.custom_ylims
    lpd.labelboxstyle.update(dict(alpha=0.8))
    fig, ax, plots = lpd.create_figure(ylims=ylims,
                                       ylabel=r'$' + displaystyle +
                                              r'\frac{G}{G^\text{norm}}$',
                                       # r'\frac{Z_f^2 G_B}{G^\text{norm}}$',
                                       ylabelpos=(0.35, 0.96), xlabelpos=(0.5, 0.01), UseTex=args.use_tex, figsize=(3, 2.5))
    ax.set_xlabel(r'$\frac{8\tau_\mathrm{F}}{\tau^2}$')  # horizontalalignment='center', verticalalignment='bottom'
    ax.set_xlim([-0.005, args.max_FlowradiusBytauT**2 * 1.15])

    # fill figure with data
    for i in range(len(ydata)):
        if not numpy.isnan(ydata_extr[i][1]):
            tauT = tauTs[i]
            xdataplot = xdata*8/tauT**2
            mycolor = lpd.get_color(tauTs, i, mintauTindex, finest_Nt_half - 1)

            # plot fit result at tf=0
            plots.append(ax.errorbar(0, ydata_extr[i][1], edata_extr[i][1], fmt='o', markersize=0, lw=0.75, mew=0.75, capsize=1.5, alpha=1, color=mycolor, zorder=1, label='{0:.3f}'.format(tauT)))

            # plot linear fit line
            x = numpy.asarray([0, xdataplot[-1]])
            ax.errorbar(x, extrapolation_ansatz(x, *ydata_extr[i]), color='black', alpha=0.25, fmt='-', lw=0.5, zorder=-10000)

            # plot data error band
            ax.fill_between(xdataplot, ydata[i]-edata[i], ydata[i]+edata[i], facecolor=mycolor, alpha=0.5, zorder=-300+i)

            # plot data points used in extrapolation
            ax.errorbar([x for x in xdataplot if args.max_FlowradiusBytauT ** 2 > x > args.min_FlowradiusBytauT ** 2],
                        [y for i, y in enumerate(ydata[i]) if args.max_FlowradiusBytauT ** 2 > xdataplot[i] > args.min_FlowradiusBytauT ** 2],
                        fmt='o', markersize=0.75, color=mycolor, alpha=1, lw=1)

    # second x-axis for flow radius
    ax2 = ax.twiny()
    new_tick_locations = numpy.array([args.min_FlowradiusBytauT**2, args.max_FlowradiusBytauT**2])
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(new_tick_locations)
    ax2.set_xticklabels(["%.1f" % z if z == 0 else "%.2f" % z for z in numpy.sqrt(new_tick_locations)])
    ax2.set_xlabel(r'$\frac{\sqrt{8\tau_\mathrm{F}}}{\tau}$', horizontalalignment='right', verticalalignment='top', bbox=lpd.labelboxstyle, zorder=999999)
    ax2.xaxis.set_label_coords(0.97, 0.97)
    ax2.tick_params(direction='in', pad=0, width=0.5)

    # legend
    lpd.legendstyle.update(dict(loc="center left", bbox_to_anchor=(1.01, 0.5), columnspacing=0.5, handlelength=0.75, labelspacing=0.1, handletextpad=0.2,
                                borderpad=0, handler_map={matplotlib.container.ErrorbarContainer: matplotlib.legend_handler.HandlerErrorbar(yerr_size=0.2)}))
    ax.legend(handles=plots, **lpd.legendstyle)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1], title=r'$\tau T$', **lpd.legendstyle).set_zorder(-1)
    ax.xaxis.set_label_coords(0.98, 0.02)

    # TODO uncomment for BB
    # ax.text(0.5, 0.9, 'preliminary', transform=ax.transAxes,
    #         fontsize=9, color='k', alpha=0.4,
    #         ha='center', va='center', rotation='0', zorder=-1000000)

    # save plot
    file = plotbasepath+"/"+args.corr+"_flow_extr_quality.pdf"
    print("saving", file)
    fig.savefig(file)
    matplotlib.pyplot.close(fig)


def plot_corr(xdata, ydata, edata, args, plotbasepath):
    # plot final double-extrapolated correlator in its own plot
    ylims = [1.5, 4] if not args.custom_ylims else args.custom_ylims
    fig, ax, plots = lpd.create_figure(xlims=[0, 0.51], ylims=ylims, xlabel=r'$\tau T$',
                                       ylabel=r'$\frac{G}'
                                              r'{ G^\text{norm}}$',
                                       UseTex=args.use_tex)
    ax.errorbar(xdata, ydata, edata, color='black', markersize=0, lw=0.75, mew=0.75)
    file = plotbasepath+"/"+args.corr+"_flow_extr.pdf"
    print("saving", file)
    fig.savefig(file)
    matplotlib.pyplot.close(fig)


def main():

    # parse cmd line arguments
    parser, requiredNamed = lpd.get_parser()

    parser.add_argument('--flowtimes_finest', help='file with list of all dimensionless flowtimes of the finest lattice', type=str)
    parser.add_argument('--finest_Nt', help='finest Nt used in previous cont extr.', type=int)
    parser.add_argument('--coarsest_Nt', help='coarsest Nt used in previous cont extr. needed for lower flow time limit.', type=int)
    parser.add_argument('--use_tex', type=bool, default=True)
    parser.add_argument('--max_FlowradiusBytauT', type=float, default=numpy.sqrt(8*0.014),
                        help='modify the flow filter based on tauT to be more/less strict. default value of 0.33 means that for each flow radius the tauT '
                             'cannot be greater than 3*flowradius, or that for each tauT the flow radius must be less than 0.33*tauT.')
    parser.add_argument('--min_FlowradiusBytauT', type=float, default=0.2)
    parser.add_argument('--max_FlowradiusBytauT_offset', type=float, default=0,
                        help='fixed offset to make upper_flowradius_limit_ stricter (by e.g. one lattice spacing 1/Nt), as the 0.33 criterion is only valid in the '
                             'continuum. on the lattice one has to be stricter. 1/Nt_coarsest is a good value.')
    parser.add_argument('--no_extr', help='do NOT perform a flow-time-to-zero extrapolation, just plot the data', action="store_true")
    parser.add_argument('--custom_ylims', help="custom y-axis limits for both plots", type=float, nargs=2)
    parser.add_argument('--plot_all_flowtimes', help="plot all flowtimes instead of only the ones that are used in the extrapolation", action="store_true")
    parser.add_argument('--input_type', help="cont: one input file per flow time, as provided by continuum_extr.py. latt: read from single file that "
                                             "contains all flow times, as provided by reduce_data.py", choices=["cont", "latt"], default="cont")
    parser.add_argument('--conftype', help="format: s064t64_b0824900*", type=str)
    parser.add_argument('--flow_cov_file', help="load covariance between flow times from this file and use it in the fits", type=str, default=None)
    parser.add_argument('--temp_subfolder', help="suffix for input path, to incorporate subfolders for different temperatures", default="", type=str)
    parser.add_argument('--basepath', help='where to look for files', type=str)
    parser.add_argument('--basepath_plot', help='where to save plots', type=str)

    args = parser.parse_args()

    # check if given arguments are valid
    if (args.input_type == "cont" and not args.coarsest_Nt) or (args.input_type == "cont" and not args.finest_Nt):
        parser.error("--input_type cont requires --coarsest_Nt and --finest_Nt")
    if args.input_type == "latt" and not args.conftype:
        parser.error("--input_type latt requires --conftype")
    if args.input_type == "latt" and args.flowtimes_finest:
        parser.error("--input_type latt prohibits --flowtimes_finest")
    if args.input_type == "cont" and not args.flowtimes_finest:
        parser.error("--input_type cont requires --flowtimes_finest")
    if (args.input_type == "latt" and args.finest_Nt) or (args.input_type == "latt" and args.coarsest_Nt):
        parser.error("--input_type latt prohibits --finest_Nt or --coarsest_Nt")
    if args.input_type == "cont":
        basepath = lpd.get_merged_data_path(args.qcdtype, args.corr, "", args.basepath)+args.temp_subfolder
        plotbasepath = lpd.get_plot_path(args.qcdtype, args.corr, "", args.basepath_plot)+args.temp_subfolder
    elif args.input_type == "latt":
        basepath = lpd.get_merged_data_path(args.qcdtype, args.corr, args.conftype, args.basepath)+args.temp_subfolder
        plotbasepath = lpd.get_plot_path(args.qcdtype, args.corr, args.conftype, args.basepath_plot)+args.temp_subfolder
        args.flowtimes_finest = lpd.get_merged_data_path(args.qcdtype, args.corr, args.conftype, args.basepath)+"flowtimes_"+args.conftype+".dat"
        _, _, nt, _ = lpd.parse_conftype(args.conftype)
        args.finest_Nt = nt
        args.coarsest_Nt = nt

    # parse some more arguments
    flowtimes = numpy.loadtxt(args.flowtimes_finest)/args.finest_Nt**2
    flowradii = numpy.sqrt(8*flowtimes)
    finest_tauTs = lpd.get_tauTs(args.finest_Nt)
    ntauT = len(finest_tauTs)

    if args.input_type == "cont":
        cont_samples = numpy.load(basepath+"/cont_extr/" + args.corr + "_cont_samples.npy")[:, :, :, 1]
        n_samples = len(cont_samples)
        cont_samples = numpy.swapaxes(cont_samples, 1, 2)
        data_std = lpd.dev_by_dist(cont_samples, axis=0)
        # shape: (1000, 18, 221)
        # TODO figure out weights for the fit. covariance matrix?
        # TODO: since we now have samples in the continuum, we could calculate the covariance across flowtime in the continuum, right? or does this not work somehow?

    elif args.input_type == "latt":
        print("ERROR currently not implemented")
        exit(0)
    # TODO reimplement this

    #     # load single lattice data of all flow times
    #     try:
    #         filename = lpd.get_merged_data_path(args.qcdtype, args.corr, args.conftype) + args.corr + "_" + args.conftype + ".dat"
    #         filename_err = lpd.get_merged_data_path(args.qcdtype, args.corr, args.conftype) + args.corr + "_err_" + args.conftype + ".dat"
    #         _, _, nt, _ = lpd.parse_conftype(args.conftype)
    #         tmp = numpy.loadtxt(filename)
    #         tmp_err = numpy.loadtxt(filename_err)
    #     except OSError:
    #         print("ERROR: could not read in \n", filename, "\n", filename_err)
    #         exit(1)
    #     print("SUCCESS: read in \n", filename, "\n", filename_err)
    #
    #     fermions, _, flowtype, gaugeaction, flowaction = lpd.parse_qcdtype(args.qcdtype)
    #
    #     for i, flowradius in enumerate(flowradii):
    #         for j in range(ntauT):
    #             XX[j][i][0] = tmp[i, j] * nt**4 / lpd.G_latt_LO_flow(j, (flowradius*nt)**2/8, args.corr, nt, flowaction, gaugeaction) #
    #             XX[j][i][1] = tmp_err[i, j] * nt**4 / lpd.G_latt_LO_flow(j, (flowradius*nt)**2/8, args.corr, nt, flowaction, gaugeaction) # (flowradius*nt)**2/8

    fermions, _, flowtype, gaugeaction, flowaction = lpd.parse_qcdtype(args.qcdtype)

    # if args.flow_cov_file:
    #     _, _, nt, _ = lpd.parse_conftype(args.conftype)
    # cov = None
    # if args.flow_cov_file is not None:
    #     cov = numpy.load(args.flow_cov_file)

    # for t in range(ntauT):
    #     for i in range(nflow):
    #         for j in range(nflow):
    #             cov[t][i,j] = cov[t][i,j] * nt**8 / lpd.G_latt_LO_flow(t, (flowradii[i]*nt)**2/8, args.corr, nt, flowaction, gaugeaction) / lpd.G_latt_LO_flow(t, (flowradii[j]*nt)**2/8, args.corr, nt, flowaction, gaugeaction)


    # declarations
    results = numpy.empty((n_samples, ntauT, 2))
    results[:] = numpy.nan
    mintauTindex = None

    # TODO parallelize over tauT?
    for n in range(n_samples):
        for i, tauT in enumerate(finest_tauTs):
            maxflowradius = lpd.upper_flowradius_limit_(tauT, args.max_FlowradiusBytauT, args.max_FlowradiusBytauT_offset)
            flowend = (numpy.abs(flowradii - maxflowradius)).argmin()
            flowstart = numpy.abs(flowradii - lpd.upper_flowradius_limit_(tauT, args.min_FlowradiusBytauT, 0)).argmin()

            # index = numpy.unique(numpy.linspace(flowstart, flowend, 10, dtype=int))
            index = numpy.arange(flowstart, flowend, 1, dtype=int)
            # flowmiddle_l = int((flowend-flowstart)/3+flowstart)
            # flowmiddle_r = int(2*(flowend-flowstart) / 3+flowstart)
            # flowmiddle = int((flowend+flowstart)/2)
            # index = numpy.unique(numpy.array([flowstart,flowmiddle_l, flowmiddle_r,flowend]))
            # index = numpy.unique(numpy.array([flowstart, flowend]))
            # index = numpy.array(range(flowstart, flowend+1))
            # print(index)
            xdata_all = flowtimes * 8/tauT**2
            xdatatmp = flowtimes[index] * 8/tauT**2
            ydatatmp = cont_samples[n][i][index]
            edatatmp = data_std[i][index]
            mask = numpy.isnan(ydatatmp)
            index = index[~mask]
            xdata = xdatatmp[~mask]
            ydata = ydatatmp[~mask]
            edata = edatatmp[~mask]

            # TODO undo this again when working on corr matrix again
            # this_cov = cov[i][index[:, numpy.newaxis], index]

            # todo get correlation matrix in here

            if len(xdata) >= 2:  # minimum amount for linear fit

                if mintauTindex is None:
                    mintauTindex = i

                if not args.no_extr:
                    # perform extrapolation
                    # print(edata)
                    # print(numpy.linalg.eig(this_cov)[0])
                    # diag = []
                    # for t in range(flowend-flowstart):
                    #     for u in range(flowend - flowstart):
                    #         if t != u:
                    #             this_cov[t, u] = 0
                    # print(this_cov)
                    # print(this_cov)
                    # this_cov_inv = numpy.linalg.inv(this_cov)
                    # print(this_cov_inv @ this_cov)
                    # fitparams = scipy.optimize.minimize(chisq, x0=numpy.asarray([-65, 3]), args=[xdata, ydata, edata, this_cov_inv])
                    # print(xdata, ydata)
                    fitparams = scipy.optimize.minimize(chisq, x0=numpy.asarray([0, 3]), args=[xdata, ydata, edata, None])
                    fitparams = fitparams.x
                    # fitparams_err = [0, 0]
                    # fitparams, fitparams_err = scipy.optimize.curve_fit(extrapolation_ansatz, xdata, ydata, p0=[-65, 3.5], sigma=this_cov)  # , maxfev=2000, method='lm', ftol=1e-12, xtol=1e-12, gtol=1e-12
                    # for i in range(2):
                    #     for j in range(2):
                    #         fitparams_err[i,j] = fitparams_err[i,j] / numpy.sqrt(fitparams_err[i,i] * fitparams_err[j,j])
                    # print(fitparams_err)
                    # fitparams_err = numpy.sqrt(numpy.diag(fitparams_err))

                    # fitparams, fitparams_err = bootstr.bootstr_from_gauss(do_fit, data=ydata, data_std_dev=edata, numb_samples=args.nsamples,
                    #                                                       sample_size=1, return_sample=False, args=[xdata, edata, this_cov_inv], nproc=20, seed=0)
                    #
                    # fitparams, fitparams_err = scipy.optimize.curve_fit(extrapolation_ansatz, xdata, ydata, p0=[-65, 3.5], sigma=edata)  # , maxfev=2000, method='lm', ftol=1e-12, xtol=1e-12, gtol=1e-12
                    # fitparams_err = numpy.sqrt(numpy.diag(fitparams_err))
                    # print(fitparams, fitparams_err)
                    # fitparams, fitparams_err = bootstr.bootstr_from_gauss(fit_sample, data=ydata, data_std_dev=edata, numb_samples=args.nsamples,
                    #                                                       sample_size=1, return_sample=False, args=[xdata, edata], nproc=20, seed=0)

                    results[n][i] = fitparams
            # print("done ", '{0:.3f}'.format(tauT))

    file = basepath + "/" + args.corr + "_flow_extr.npy"
    print("saving", file)
    numpy.save(file, results)

    # original data points
    data_mean = numpy.nanmedian(cont_samples, axis=0)

    results_mean = numpy.median(results, axis=0)
    results_std = lpd.dev_by_dist(results, axis=0)

    xdata = finest_tauTs

    plot_extrapolation(flowtimes, data_mean, data_std, results_mean, results_std, args, mintauTindex, plotbasepath)

    # plot and save final correlator
    ydata_extr = results_mean[:, 1]
    edata_extr = results_std[:, 1]
    plot_corr(finest_tauTs, ydata_extr, edata_extr, args, plotbasepath)
    correlator = numpy.column_stack((xdata, ydata_extr, edata_extr))
    print(correlator)
    file = basepath + "/" + args.corr + "_flow_extr.txt"
    print("saving", file)
    numpy.savetxt(file , correlator, header="                tauT                G/Gnorm                    err", fmt='%22.15e')


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
