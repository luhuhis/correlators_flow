#!/usr/bin/env python3
import numpy
import matplotlib
import lib_process_data as lpd
import scipy.optimize
import scipy.interpolate


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

    density_weight = []
    for i in range(len(xdata)):
        if i == len(xdata)-1:
            density_weight.append((xdata[i]-xdata[i-1])/xdata[i])
        elif i == 0:
            density_weight.append((xdata[i + 1]-xdata[i])/xdata[i])
        else:
            density_weight.append((xdata[i + 1] - xdata[i - 1]) / 2/ xdata[i])

    density_weight = numpy.asarray(density_weight)
    density_weight = 1  # TODO undo
    r = (ydata - extrapolation_ansatz(xdata, params[0], params[1]))/edata
    return r.T @ corr_inv @ (r * density_weight)


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
    fig, ax, plots = lpd.create_figure(ylims=ylims, ylabel=r'$' + displaystyle + r'\frac{'+ylabel_prefix+r'G}{G^\text{norm}}$', xlabel=r'$'+displaystyle+r'{8\tau_\mathrm{F}}/{\tau^2}$')
                                       # r'\frac{Z_f^2 G_B}{G^\text{norm}}$', UseTex=args.use_tex)
    ax.set_xlim([-0.005, args.max_FlowradiusBytauT**2 * 1.15])

    # fill figure with data
    for i in range(len(ydata)):
        tauT = tauTs[i]
        if tauT >= args.min_tauT_plot:

            if mintauTindex is None:
                mintauTindex = i

            xdataplot = xdata*8/tauT**2
            mycolor = lpd.get_color(tauTs, i, mintauTindex, finest_Nt_half - 1)

            if not args.no_extr and not numpy.isnan(ydata_extr[i][1]):
                # plot fit result at tf=0
                plots.append(ax.errorbar(0, ydata_extr[i][1], edata_extr[i][1], fmt='|', color=mycolor, zorder=1, label='{0:.3f}'.format(tauT)))

                # plot linear fit line
                x = numpy.asarray([0, xdataplot[-1]])
                ax.errorbar(x, extrapolation_ansatz(x, *ydata_extr[i]), color='black', alpha=0.25, fmt='-',  zorder=-10000)  #lw=0.5,

                # plot data points used in extrapolation
                ax.errorbar(xdataplot[indices[i]], ydata[i][indices[i]], edata[i][indices[i]], fmt='|', color=mycolor)

            # plot data error band
            # TODO add option to only plot data where flowradius > 1.5a
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
    # ax.xaxis.set_label_coords(0.98, 0.02)

    # TODO uncomment for BB
    # ax.text(0.5, 0.9, 'preliminary', transform=ax.transAxes,
    #         fontsize=9, color='k', alpha=0.4,
    #         ha='center', va='center', rotation='0', zorder=-1000000)

    # save plot
    file = plotbasepath+"/"+args.corr+"_flow_extr_quality"+args.output_suffix+".pdf"
    print("saving", file)
    fig.savefig(file)
    matplotlib.pyplot.close(fig)


def plot_corr(args, xdata, ydata, edata, plotbasepath):
    # plot final double-extrapolated correlator in its own plot
    ylims = [1.5, 4] if not args.custom_ylims else args.custom_ylims
    fig, ax, plots = lpd.create_figure(xlims=[0, 0.51], ylims=ylims, xlabel=r'$\tau T$',
                                       ylabel=r'$\frac{G}'
                                              r'{ G^\text{norm}}$',
                                       UseTex=args.use_tex)
    ax.errorbar(xdata, ydata, edata, color='black', markersize=0)
    file = plotbasepath+"/"+args.corr+"_flow_extr"+args.output_suffix+".pdf"
    print("saving", file)
    fig.savefig(file)
    matplotlib.pyplot.close(fig)


def do_flow_extr(index, tauTs, cont_samples, data_std, n_samples, args, flowtimes, flowradii):
    """ do the extrapolation for all samples of one tauT"""

    tauT = tauTs[index]

    results = numpy.empty((n_samples, 2))
    results[:] = numpy.nan

    maxflowradius = lpd.upper_flowradius_limit_(tauT, args.max_FlowradiusBytauT, args.max_FlowradiusBytauT_offset)
    flowend = (numpy.abs(flowradii - maxflowradius)).argmin()
    flowstart = numpy.abs(flowradii - lpd.upper_flowradius_limit_(tauT, args.min_FlowradiusBytauT, 0)).argmin()
    # indices = numpy.arange(flowstart, flowend, 1, dtype=int)

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
    print("done ", '{0:.3f}'.format(tauT))
    return results, indices


def main():

    # parse cmd line arguments
    parser, requiredNamed = lpd.get_parser()

    parser.add_argument('--flowtimes_finest', help='file with list of all dimensionless flowtimes of the finest lattice', type=str)
    parser.add_argument('--finest_Nt', help='finest Nt used in previous cont extr.', type=int)
    parser.add_argument('--coarsest_Nt', help='coarsest Nt used in previous cont extr. needed for lower flow time limit.', type=int)
    parser.add_argument('--use_tex', action="store_true", default=False)
    parser.add_argument('--max_FlowradiusBytauT', type=float, default=numpy.sqrt(8*0.014),
                        help='modify the flow filter based on tauT to be more/less strict. default value of 0.33 means that for each flow radius the tauT '
                             'cannot be greater than 3*flowradius, or that for each tauT the flow radius must be less than 0.33*tauT.')
    parser.add_argument('--min_FlowradiusBytauT', type=float, default=0.2)
    parser.add_argument('--max_FlowradiusBytauT_offset', type=float, default=0,
                        help='fixed offset to make upper_flowradius_limit_ stricter (by e.g. one lattice spacing 1/Nt), as the 0.33 criterion is only valid in the '
                             'continuum. on the lattice one has to be stricter. 1/Nt_coarsest is a good value.')
    parser.add_argument('--min_tauT_plot', default=0, type=float)
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
    parser.add_argument('--Zf2_file', help="Z_f^2 for BB correlator renormalization", default=None, type=str)
    parser.add_argument('--output_suffix', help="add this string to the end of plot file names", default="", type=str)
    parser.add_argument('--slope_bounds', help="bound the slope fit parameter between these values", default=(None, None), nargs=2, type=float)

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
        cont_samples = numpy.load(basepath+"/cont_extr/" + args.corr + "_cont_samples.npy")[:, :, :, 0]
        print(cont_samples.shape)
        print(len(flowtimes))
        n_samples = len(cont_samples)
        if args.Zf2_file is not None:
            tfT2, Zf2 = numpy.loadtxt(args.Zf2_file, unpack=True)
            Zf2_int = scipy.interpolate.InterpolatedUnivariateSpline(tfT2, Zf2, k=3, ext=1)
            for j in range(len(flowtimes)):
                if Zf2_int(flowtimes[j]) == 0:
                    print("warn", flowtimes[j])
                cont_samples[:, j, :] *= Zf2_int(flowtimes[j])
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

    fitparams, indices = lpd.parallel_function_eval(do_flow_extr, range(0, ntauT), 20, finest_tauTs, cont_samples, data_std, n_samples, args, flowtimes, flowradii)

    results = numpy.swapaxes(numpy.asarray(fitparams), 0, 1)

    file = basepath + "/" + args.corr + "_flow_extr.npy"
    print("saving all samples in", file)
    numpy.save(file, results)
    print("results shape:", results.shape)

    # TODO save correlation matrices
    # XX_samples = numpy.swapaxes(XX_samples, 0, 2)  # change data structure such that numpy.cov understands it
    # pcov = []
    # for i in range(int(nt/2)):
    #     pcov.append(numpy.corrcoef(XX_samples[i]))
    # pcov = numpy.asarray(pcov)
    # numpy.save(lpd.print_var("write", file_prefix + "_flow_cov_" + conftype + ".npy"), pcov)
    # print(pcov.shape)

    # original data points
    data_mean = numpy.nanmedian(cont_samples, axis=0)

    results_mean = numpy.median(results, axis=0)
    results_std = lpd.dev_by_dist(results, axis=0)
    print(results_mean[-1], results_std[-1])
    xdata = finest_tauTs

    plot_extrapolation(args, flowtimes, data_mean, data_std, results_mean, results_std, indices, plotbasepath)

    # plot and save final correlator
    ydata_extr = results_mean[:, 1]
    edata_extr = results_std[:, 1]
    plot_corr(args, finest_tauTs, ydata_extr, edata_extr, plotbasepath)
    correlator = numpy.column_stack((xdata, ydata_extr, edata_extr))
    print(correlator)
    file = basepath + "/" + args.corr + "_flow_extr.txt"
    print("saving mean and error in ", file)
    numpy.savetxt(file, correlator, header="                tauT                G/Gnorm                    err", fmt='%22.15e')


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()

# index = numpy.unique(numpy.linspace(flowstart, flowend, 10, dtype=int))
# flowmiddle_l = int((flowend-flowstart)/3+flowstart)
# flowmiddle_r = int(2*(flowend-flowstart) / 3+flowstart)
# flowmiddle = int((flowend+flowstart)/2)
# index = numpy.unique(numpy.array([flowstart,flowmiddle_l, flowmiddle_r,flowend]))
# index = numpy.unique(numpy.array([flowstart, flowend]))
# index = numpy.array(range(flowstart, flowend+1))
# print(index)
# xdata_all = flowtimes * 8 / tauT ** 2


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