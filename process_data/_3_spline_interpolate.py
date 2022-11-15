#!/usr/bin/env python3

# this script interpolates the XX correlator using interpolating (not smoothing!) third order splines.

import numpy
import lib_process_data as lpd
import scipy.interpolate
import matplotlib.pyplot
from matplotlib.backends.backend_pdf import PdfPages


def interpolate_XX_flow(xdata, ydata, ydata_norm, output_xdata1, output_xdata2):
    spline = scipy.interpolate.CubicSpline(xdata, ydata, bc_type=((2, 0.0), (1, 0.0)))  # true interpolating spline
    norm = scipy.interpolate.CubicSpline(xdata, ydata_norm, bc_type=((2, 0.0), (1, 0.0)))  # true interpolating spline
    return spline(output_xdata1)/norm(output_xdata1), spline(output_xdata2)/norm(output_xdata2)


def plot(args, flowradius, xpointsplot, theplotdata, merged_data_path, index, flowtime, nt, flowaction, gaugeaction):
    ypointsplot = numpy.median(theplotdata, axis=0)
    epointsplot = lpd.dev_by_dist(theplotdata, axis=0)

    # plot interpolations and underlying data points
    ylabel = 'G'
    fig, ax, plots = lpd.create_figure(xlims=[0, 0.51], ylims=args.ylims, xlabel=r'$\tau T$', ylabel=ylabel, xlabelpos=(0.95, 0.05), constrained_layout=False)
    ax.set_title(r'$ \sqrt{8\tau_F}T = $' + '{0:.3f}'.format(flowradius))
    ax.axvline(x=lpd.lower_tauT_limit_(flowradius, args.max_FlowradiusBytauT, args.max_FlowradiusBytauT_offset), **lpd.verticallinestyle)
    ax.fill_between(xpointsplot, ypointsplot - epointsplot, ypointsplot + epointsplot, alpha=0.5)
    ax.errorbar(xpointsplot, ypointsplot, fmt='-', lw=0.5, mew=0)
    XX = numpy.loadtxt(merged_data_path + "/" + args.corr + "_" + args.conftype + ".dat")
    XX_err = numpy.loadtxt(merged_data_path + "/" + args.corr + "_err_" + args.conftype + ".dat")
    for a in range(len(XX[index])):
        XX[index, a] = XX[index, a] * nt ** 4 / lpd.G_latt_LO_flow(a, flowtime, args.corr, nt, flowaction, gaugeaction)
        XX_err[index, a] = XX_err[index, a] * nt ** 4 / lpd.G_latt_LO_flow(a, flowtime, args.corr, nt, flowaction, gaugeaction)
    ax.errorbar(lpd.get_tauTs(nt), XX[index], XX_err[index], **lpd.plotstyle_add_point_single)

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

    fig = plot(args, flowradius, xpointsplot, theplotdata, merged_data_path, index, flowtime, nt, flowaction, gaugeaction)

    return theoutputdata, fig


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

    args = parser.parse_args()

    # file paths
    merged_data_path = lpd.get_merged_data_path(args.qcdtype, args.corr, args.conftype, args.basepath)
    outputfolder_data = merged_data_path+"/interpolations/"
    outputfolder_plot = lpd.get_plot_path(args.qcdtype, args.corr, args.conftype, args.basepath_plot)
    lpd.create_folder(outputfolder_plot, outputfolder_data)

    # load flow times
    flowtimes = numpy.loadtxt(merged_data_path + "/flowtimes_" + args.conftype + ".dat")
    indices = range(0, len(flowtimes))

    # load data
    XX_samples = numpy.load(merged_data_path+"/"+args.corr + "_" + args.conftype + "_samples.npy")

    matplotlib.rcParams['figure.max_open_warning'] = 0  # suppress warning due to possibly large number of figures...
    interpolations, figs = lpd.parallel_function_eval(wrapper, indices, args.nproc, flowtimes, XX_samples, args, merged_data_path)

    filepath = outputfolder_plot + "/" + args.corr + "_interpolation.pdf"
    print("saving ", filepath)

    lpd.set_rc_params()  # for some reason we need to call this here...
    with PdfPages(filepath) as pdf:
        for fig in figs:
            pdf.savefig(fig, bbox_inches='tight', pad_inches=0.05)
            matplotlib.pyplot.close(fig)

    # save data
    numpy.save(merged_data_path + args.corr + "_" + args.conftype + "_interpolation_samples.npy", numpy.asarray(interpolations))


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
