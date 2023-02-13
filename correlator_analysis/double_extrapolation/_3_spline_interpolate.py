#!/usr/bin/env python3

# this script interpolates the XX correlator using interpolating (not smoothing!) third order splines.

import numpy
import lib_process_data as lpd
import scipy.interpolate
import matplotlib.pyplot
from matplotlib.backends.backend_pdf import PdfPages


import warnings
warnings.filterwarnings('ignore', r'All-NaN slice encountered.*')


def interpolate_XX_flow(xdata, ydata, ydata_norm, output_xdata1, output_xdata2):
    spline = scipy.interpolate.CubicSpline(xdata, ydata, bc_type=((2, 0.0), (1, 0.0)))  # true interpolating spline
    norm = scipy.interpolate.CubicSpline(xdata, ydata_norm, bc_type=((2, 0.0), (1, 0.0)))  # true interpolating spline
    return spline(output_xdata1)/norm(output_xdata1), spline(output_xdata2)/norm(output_xdata2)


def plot(args, flowradius, x, y, yerr, merged_data_path, index, flowtime, nt, flowaction, gaugeaction):

    # plot interpolations and underlying data points
    ylabel = r'$\displaystyle \frac{G' + lpd.get_corr_subscript(args.corr) + r'}{G^\mathrm{norm}}$'
    fig, ax, plots = lpd.create_figure(xlims=[0, 0.52], ylims=args.ylims, xlabel=r'$\tau T$', ylabel=ylabel, constrained_layout=True)
    ax.text(0.99, 0.99, r'$ \sqrt{8\tau_F}T = ' + '{0:.3f}'.format(flowradius)+'$', ha='right', va='top', transform=ax.transAxes, zorder=-1000, bbox=lpd.labelboxstyle)
    lower_limit = lpd.lower_tauT_limit_(flowradius, 0.33, 0)
    min_index = numpy.fabs(x-lower_limit).argmin()
    ax.axvline(x=lower_limit, **lpd.verticallinestyle)
    ax.text(lower_limit, args.lower_limit_text_pos, r'$ 3\sqrt{8\tau_F}T$', ha='center', va='center', zorder=-1000, bbox=lpd.labelboxstyle)
    ax.set_xticks((0.0, 0.1, 0.2, 0.3, 0.4, 0.5))

    # interpolations
    ax.fill_between(x[min_index:], y[min_index:] - yerr[min_index:], y[min_index:] + yerr[min_index:], label="int.")
    ax.errorbar(x[min_index:], y[min_index:], fmt='-', zorder=-10)

    #data
    XX = numpy.loadtxt(merged_data_path + "/" + args.corr + "_" + args.conftype + ".dat")
    XX_err = numpy.loadtxt(merged_data_path + "/" + args.corr + "_err_" + args.conftype + ".dat")
    for a in range(len(XX[index])):
        XX[index, a] = XX[index, a] * nt ** 4 / lpd.G_latt_LO_flow(a, flowtime, args.corr, nt, flowaction, gaugeaction)
        XX_err[index, a] = XX_err[index, a] * nt ** 4 / numpy.fabs(lpd.G_latt_LO_flow(a, flowtime, args.corr, nt, flowaction, gaugeaction))
    x = lpd.get_tauTs(nt)
    min_index = numpy.ceil(numpy.fabs(x*nt - lower_limit*nt)).argmin()+1
    ax.errorbar(x[min_index:], XX[index][min_index:], XX_err[index][min_index:], fmt='|', label="data")

    ax.legend(loc="center right", bbox_to_anchor=(0.99, 0.25), **lpd.leg_err_size())
    return fig


def tau_interpolation(index, flowtimes, XX_samples, args, merged_data_path):

    beta, ns, nt, nt_half = lpd.parse_conftype(args.conftype)
    _, _, _, gaugeaction, flowaction = lpd.parse_qcdtype(args.qcdtype)

    flowtime = flowtimes[index]
    flowradius = numpy.sqrt(flowtimes[index]*8)/nt

    # get interpolation of the tree-level improvement factor
    ydata_norm = []
    for a in range(nt_half):
        ydata_norm.append(
            lpd.G_latt_LO_flow(a, flowtime, args.corr, nt, flowaction, gaugeaction) / lpd.G_latt_LO_flow(a, 0, args.corr, nt, flowaction, gaugeaction))

    theoutputdata = []
    theplotdata = []
    xpoints = lpd.get_tauTs(args.int_Nt)
    xpointsplot = numpy.linspace(0, xpoints[-1], 1000)

    # perform spline fits for each sample
    for m in range(args.nsamples):
        ydata = numpy.asarray(XX_samples[m, index])
        for a in range(nt_half):
            # We use tf=0 here so that the interpolation works better. Afterwards we divide it out again.
            ydata[a] = ydata[a] * nt ** 4 / lpd.G_latt_LO_flow(a, 0, args.corr, nt, flowaction, gaugeaction)
        output, plot_output = interpolate_XX_flow(lpd.get_tauTs(nt), ydata, ydata_norm, xpoints, xpointsplot)
        theoutputdata.append(output)
        theplotdata.append(plot_output)

    yplot = numpy.mean(theplotdata, axis=0)
    eplot = numpy.std(theplotdata, axis=0)

    fig = plot(args, flowradius, xpointsplot, yplot, eplot, merged_data_path, index, flowtime, nt, flowaction, gaugeaction)

    return theoutputdata, fig, yplot


def interpolate_data(xdata, ydata, output_xdata, bc=((2, 0.0), (2, 0.0))):
    spline = scipy.interpolate.CubicSpline(xdata, ydata, bc_type=bc, extrapolate=False)
    return spline(output_xdata)


def parse_args():
    # parse cmd line arguments
    parser, requiredNamed = lpd.get_parser()

    requiredNamed.add_argument('--conftype', help="format: s096t20_b0824900 for quenched or s096t20_b0824900_m002022_m01011 for hisq", required=True)
    parser.add_argument('--nsamples', help="number of gaussian bootstrap samples that are contained in the input files", type=int, default=1000)
    parser.add_argument('--int_Nt', help='use tauT of this Nt as xdata for the interpolation output', type=int, default=36)
    parser.add_argument('--max_FlowradiusBytauT', type=float, default=numpy.sqrt(8*0.014),
                        help='modify the tauT filter based on flow time to be more/less strict. default value of 0.33 means that for each tauT the flow radius '
                             'cannot be greater than 0.33*tauT, or that for each flow radius the tauT must be atleast 3*flowradius.')
    parser.add_argument('--max_FlowradiusBytauT_offset', type=float, default=1/20,
                        help='fixed offset to make lower_tauT_limit stricter (by e.g. one lattice spacing 1/Nt), as the 0.33 criterion is only valid in the '
                             'continuum. on the lattice one has to be stricter. 1/Nt_coarsest is a good value.')
    parser.add_argument('--ylims', default=[-0.2, 4], nargs=2, type=float, help="custom ylims for plot")
    parser.add_argument('--nproc', type=int, default=20, help="number of parallel processes (different flow times can be interpolated independently)")
    parser.add_argument('--basepath', type=str, help="where to look for the data")
    parser.add_argument('--basepath_plot', type=str, help="where to save the output plots")
    parser.add_argument('--lower_limit_text_pos', type=float)

    args = parser.parse_args()

    return args


def save_data(args, merged_data_path, interpolations, interpolation_mean_for_plotting, suffix="interpolation"):
    interpolations = numpy.asarray(interpolations)
    interpolation_mean_for_plotting = numpy.asarray(interpolation_mean_for_plotting)
    file = merged_data_path + args.corr + "_" + args.conftype + "_" + suffix + "_samples.npy"
    print("save ", file)
    numpy.save(file, interpolations)
    file = merged_data_path + args.corr + "_" + args.conftype + "_" + suffix + "_mean.npy"
    print("save ", file)
    numpy.save(file, interpolation_mean_for_plotting)


def save_figs(args, figs, suffix="interpolation"):
    outputfolder_plot = lpd.get_plot_path(args.qcdtype, args.corr, args.conftype, args.basepath_plot)
    lpd.create_folder(outputfolder_plot)
    filepath = outputfolder_plot + "/" + args.corr + "_" + suffix + ".pdf"
    print("saving ", filepath)

    lpd.set_rc_params()  # for some reason we need to call this here...
    with PdfPages(filepath) as pdf:
        for fig in figs:
            pdf.savefig(fig)
            matplotlib.pyplot.close(fig)


def convert_sqrt8tauFT_to_tauFBya2(sqrt8tauFT, nt):
    return (sqrt8tauFT*nt)**2/8


def traditional_interpolation(args, merged_data_path, flowtimes, XX_samples):
    matplotlib.rcParams['figure.max_open_warning'] = 0  # suppress warning due to possibly large number of figures...
    indices = range(len(flowtimes))
    interpolations, figs, interpolation_mean_for_plotting = lpd.parallel_function_eval(tau_interpolation, indices, args.nproc, flowtimes, XX_samples, args, merged_data_path)

    save_figs(args, figs)
    save_data(args, merged_data_path, interpolations, interpolation_mean_for_plotting)


def interpolate_to_relative_flowtimes(tauT_index, tauTs, XX_samples, args, absolute_flowtimes, relative_flowradii):

    _, _, nt, _ = lpd.parse_conftype(args.conftype)
    _, _, _, gaugeaction, flowaction = lpd.parse_qcdtype(args.qcdtype)

    new_XX_samples = numpy.empty((args.nsamples, len(relative_flowradii)))

    # only use well-behaved part for interpolation
    flow_indices = (absolute_flowtimes >= 1 / 8) & (absolute_flowtimes <= convert_sqrt8tauFT_to_tauFBya2(0.35 * 0.5, nt))
    new_absolute_flowtimes = convert_sqrt8tauFT_to_tauFBya2(relative_flowradii * tauTs[tauT_index], nt)
    factor = nt**4 / lpd.G_latt_LO_flow(tauT_index, absolute_flowtimes[flow_indices], args.corr, nt, flowaction, gaugeaction)
    for m in range(args.nsamples):
        ydata_old = XX_samples[m, flow_indices, tauT_index] * factor
        new_XX_samples[m] = interpolate_data(absolute_flowtimes[flow_indices], ydata_old, new_absolute_flowtimes)

    return new_XX_samples


def interpolate_tauTs_at_relative_flowtime(flow_index, args, XX_samples, int_tauTs):
    def first_non_nan_index(arr):
        for i, value in enumerate(arr):
            if not numpy.isnan(value):
                return i
        return None

    def first_non_nan_index_from_end(arr):
        for i in range(len(arr) - 1, -1, -1):
            if not numpy.isnan(arr[i]):
                return i
        return None

    _, _, nt, _ = lpd.parse_conftype(args.conftype)
    orig_tauTs = lpd.get_tauTs(nt)

    new_XX_samples = numpy.empty((args.nsamples, len(int_tauTs)))

    first_not_nan_index = first_non_nan_index(XX_samples[0, flow_index])
    last_not_nan_index = first_non_nan_index_from_end(XX_samples[0, flow_index])

    for m in range(args.nsamples):
        ydata_old = XX_samples[m, flow_index][first_not_nan_index:last_not_nan_index+1]
        new_XX_samples[m] = interpolate_data(orig_tauTs[first_not_nan_index:last_not_nan_index+1], ydata_old, int_tauTs, bc=((2, 0.0), (1, 0.0)))

    return new_XX_samples


def plot_relative_flow_ints(flowindex, args, relflow_range, orig_XX_samples, int_XX_samples, int_xdata):
    relflow = relflow_range[flowindex]
    _, _, nt, _ = lpd.parse_conftype(args.conftype)

    orig_xdata = lpd.get_tauTs(nt)
    orig_ydata = numpy.nanmedian(orig_XX_samples, axis=0)[flowindex]
    orig_edata = lpd.dev_by_dist(orig_XX_samples, axis=0)[flowindex]

    int_ydata = numpy.nanmedian(int_XX_samples, axis=0)[flowindex]
    int_edata = lpd.dev_by_dist(int_XX_samples, axis=0)[flowindex]

    # plot interpolations and underlying data points
    ylabel = r'$\displaystyle \frac{G' + lpd.get_corr_subscript(args.corr) + r'}{G^\mathrm{norm}}$'
    fig, ax, plots = lpd.create_figure(xlims=[0, 0.52], ylims=args.ylims, xlabel=r'$\tau T$', ylabel=ylabel, constrained_layout=True)
    ax.text(0.99, 0.99, r'$ \sqrt{8\tau_F}/\tau = ' + '{0:.3f}'.format(relflow) + '$', ha='right', va='top', transform=ax.transAxes, zorder=-1000,
            bbox=lpd.labelboxstyle)
    ax.set_xticks((0.0, 0.1, 0.2, 0.3, 0.4, 0.5))

    # interpolations
    ax.errorbar(orig_xdata, orig_ydata, orig_edata, label="data", fmt='|', color="C0", zorder=0)
    ax.errorbar(int_xdata, int_ydata, lw=0.5, color="C1", zorder=1)
    ax.fill_between(int_xdata, int_ydata-int_edata, int_ydata+int_edata, label="int.", facecolor="C1", zorder=2)

    ax.legend(loc="center right", bbox_to_anchor=(0.99, 0.25), **lpd.leg_err_size())
    return fig


def new_interpolation(args, merged_data_path, flowtimes, XX_samples):
    _, _, nt, nt_half = lpd.parse_conftype(args.conftype)

    min_relflow = 0.2
    max_relflow = 0.35
    stepsize = 0.005
    relflow_range = numpy.arange(min_relflow, max_relflow + stepsize, stepsize)
    nflow = len(relflow_range)

    orig_xdata = lpd.get_tauTs(nt)

    # interpolate_to_relative_flowtimes
    flow_int_XX_samples = lpd.parallel_function_eval(interpolate_to_relative_flowtimes, range(nt_half), args.nproc, orig_xdata, XX_samples, args, flowtimes, relflow_range)
    flow_int_XX_samples = numpy.asarray(flow_int_XX_samples).swapaxes(0, 2)
    flow_int_XX_samples = flow_int_XX_samples.swapaxes(0, 1)

    # interpolate to tauTs of finest lattice
    tauT_int_flow_int_XX_samples = lpd.parallel_function_eval(interpolate_tauTs_at_relative_flowtime, range(len(relflow_range)), args.nproc,
                                                              args, flow_int_XX_samples, lpd.get_tauTs(args.int_Nt))
    tauT_int_flow_int_XX_samples = numpy.asarray(tauT_int_flow_int_XX_samples).swapaxes(0, 1)

    # interpolate to all the tauTs for the plot
    int_xdata = numpy.linspace(0, 0.5, 200)
    plot_tauT_int_flow_int_XX_samples = lpd.parallel_function_eval(interpolate_tauTs_at_relative_flowtime, range(len(relflow_range)), args.nproc,
                                                                   args, flow_int_XX_samples, int_xdata)
    plot_tauT_int_flow_int_XX_samples = numpy.asarray(plot_tauT_int_flow_int_XX_samples).swapaxes(0, 1)

    numpy.savetxt(merged_data_path + args.corr + "_" + args.conftype + "_relflows.txt", relflow_range, fmt='%.4f')
    save_data(args, merged_data_path, tauT_int_flow_int_XX_samples, numpy.median(plot_tauT_int_flow_int_XX_samples, axis=0), suffix="interpolation_relflow")

    # create figs
    figs = lpd.parallel_function_eval(plot_relative_flow_ints, range(nflow), args.nproc, args, relflow_range,
                                      flow_int_XX_samples, plot_tauT_int_flow_int_XX_samples, int_xdata)

    save_figs(args, figs, suffix="interpolation_relflow")


def main():

    args = parse_args()

    # load data
    merged_data_path = lpd.get_merged_data_path(args.qcdtype, args.corr, args.conftype, args.basepath)
    flowtimes = numpy.loadtxt(merged_data_path + "/flowtimes_" + args.conftype + ".dat")
    XX_samples = numpy.load(merged_data_path + "/" + args.corr + "_" + args.conftype + "_samples.npy")

    traditional_interpolation(args, merged_data_path, flowtimes, XX_samples)
    new_interpolation(args, merged_data_path, flowtimes, XX_samples)


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
