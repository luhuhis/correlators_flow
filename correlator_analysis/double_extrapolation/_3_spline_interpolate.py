#!/usr/bin/env python3

# this script interpolates the XX correlator using interpolating (not smoothing!) third order splines.

import numpy
import lib_process_data as lpd
import scipy.interpolate
import matplotlib.pyplot
from matplotlib.backends.backend_pdf import PdfPages
# import warnings
# warnings.simplefilter("ignore", OptimizeWarning)
# warnings.simplefilter('error', UserWarning)
# warnings.filterwarnings('ignore', r'Calling figure.constrained_layout, but figure not setup to do constrained layout.*')


def interpolate_XX_flow(xdata, ydata, ydata_norm, output_xdata1, output_xdata2):
    spline = scipy.interpolate.CubicSpline(xdata, ydata, bc_type=((2, 0.0), (1, 0.0)))  # true interpolating spline
    norm = scipy.interpolate.CubicSpline(xdata, ydata_norm, bc_type=((2, 0.0), (1, 0.0)))  # true interpolating spline
    return spline(output_xdata1)/norm(output_xdata1), spline(output_xdata2)/norm(output_xdata2)


def plot(args, flowradius, x, y, yerr, merged_data_path, index, flowtime, nt, flowaction, gaugeaction):

    # plot interpolations and underlying data points
    ylabel = r'$\displaystyle \frac{G}{G^\mathrm{norm}}$'
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

    # data
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


def wrapper(index, flowtimes, XX_samples, args, merged_data_path):

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


def main():

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

    # file paths
    merged_data_path = lpd.get_merged_data_path(args.qcdtype, args.corr, args.conftype, args.basepath)
    outputfolder_plot = lpd.get_plot_path(args.qcdtype, args.corr, args.conftype, args.basepath_plot)
    lpd.create_folder(outputfolder_plot)

    # load flow times
    flowtimes = numpy.loadtxt(merged_data_path + "/flowtimes_" + args.conftype + ".dat")
    indices = range(0, len(flowtimes))

    # load data
    XX_samples = numpy.load(merged_data_path+"/"+args.corr + "_" + args.conftype + "_samples.npy")

    matplotlib.rcParams['figure.max_open_warning'] = 0  # suppress warning due to possibly large number of figures...
    interpolations, figs, interpolation_mean_for_plotting = lpd.parallel_function_eval(wrapper, indices, args.nproc, flowtimes, XX_samples, args, merged_data_path)

    filepath = outputfolder_plot + "/" + args.corr + "_interpolation.pdf"
    print("saving ", filepath)

    lpd.set_rc_params()  # for some reason we need to call this here...
    with PdfPages(filepath) as pdf:
        for fig in figs:
            pdf.savefig(fig)
            matplotlib.pyplot.close(fig)

    # save data
    interpolations = numpy.asarray(interpolations)
    interpolation_mean_for_plotting = numpy.asarray(interpolation_mean_for_plotting)
    print(interpolation_mean_for_plotting.shape)
    numpy.save(merged_data_path + args.corr + "_" + args.conftype + "_interpolation_samples.npy", interpolations)
    numpy.save(merged_data_path + args.corr + "_" + args.conftype + "_interpolation_mean.npy", interpolation_mean_for_plotting)


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()