#!/usr/local/bin/python3.7m -u

# this script interpolates the XX correlator using interpolating (not smoothing!) third order splines. it does that for one single given flow time.

import numpy
import lib_process_data as lpd
import scipy.interpolate
import sys


def filter_tauTs(flowradius, xdata, max_FlowradiusBytauT, max_FlowradiusBytauT_offset, int_left_tauT_offset):
    """ remove datapoints that the flow has destroyed """
    indices = numpy.asarray([k for k, j in enumerate(xdata) if j >= lpd.lower_tauT_limit_(flowradius, max_FlowradiusBytauT, max_FlowradiusBytauT_offset) - int_left_tauT_offset])
    if len(indices) == 0:
        sys.exit("Error at "+"flow_radius="+'{0:.4f}'.format(flowradius)+". There are no datapoints that the flow hasn't destroyed according to pert. flow limits.")
    # min_ind = numpy.min(indices)
    # add one helper point outside of flow limit to stabilize lower part of interpolation
    # if min_ind >= 1:
    #     indices = numpy.asarray([min_ind-1, *indices])
    # if min_ind >= 2:
    #     indices = numpy.asarray([min_ind-2, *indices])
    xdata = numpy.asarray(xdata)
    xdata = xdata[indices]
    return xdata, indices


def filter_corr_data(flowradius, xdata, ydata, max_FlowradiusBytauT, max_FlowradiusBytauT_offset, int_left_tauT_offset):
    """ remove datapoints that the flow has destroyed """
    xdata, indices = filter_tauTs(flowradius, xdata, max_FlowradiusBytauT, max_FlowradiusBytauT_offset, int_left_tauT_offset)
    ydata = numpy.asarray(ydata)
    ydata = ydata[indices]
    return xdata, ydata


# TODO find a way to do smoothing AND constrain the first derivative at 0.5 to zero
def interpolate_XX_flow(xdata, ydata, weights, output_xdata1, output_xdata2):
    # spline = scipy.interpolate.CubicSpline(xdata, ydata, bc_type=((2, 0.0), (1, 0.0)))  # true interpolating spline
    spline = scipy.interpolate.UnivariateSpline(xdata, ydata, w=weights, k=3, s=0.001)  # smoothing spline with minimal smoothing and fit weights
    return spline(output_xdata1), spline(output_xdata2)


def main():

    # parse cmd line arguments
    parser, requiredNamed = lpd.get_parser()

    requiredNamed.add_argument('--conftype', help="format: s096t20_b0824900 for quenched or s096t20_b0824900_m002022_m01011 for hisq", required=True)
    requiredNamed.add_argument('--flow_index', help='which flow time to interpolate', type=int, required=True)

    parser.add_argument('--use_imp', help='whether to use tree-level improvement', type=bool, default=True)
    parser.add_argument('--nsamples', help="number of gaussian bootstrap samples that are contained in the input files", type=int, default=1000)
    parser.add_argument('--int_Nt', help='use tauT of this Nt as xdata for the interpolation output', type=int, default=36)
    parser.add_argument('--int_left_tauT_offset', help='include points <int_left_tauT_offset> below the lower tauT limit in the spline calculation. this helps '
                                                       'to make the interpolation smooth/consistent across varying flow time.  a good value is ~2/Nt_coarsest.',
                        default=2/20, type=float)
    parser.add_argument('--max_FlowradiusBytauT', type=float, default=numpy.sqrt(8*0.014),
                        help='modify the tauT filter based on flow time to be more/less strict. default value of 0.33 means that for each tauT the flow radius '
                             'cannot be greater than 0.33*tauT, or that for each flow radius the tauT must be atleast 3*flowradius.')
    parser.add_argument('--max_FlowradiusBytauT_offset', type=float, default=1/20,
                        help='fixed offset to make lower_tauT_limit stricter (by e.g. one lattice spacing 1/Nt), as the 0.33 criterion is only valid in the '
                             'continuum. on the lattice one has to be stricter. 1/Nt_coarsest is a good value.')

    args = parser.parse_args()

    beta, ns, nt, nt_half = lpd.parse_conftype(args.conftype)

    # file paths
    merged_data_path = lpd.get_merged_data_path(args.qcdtype, args.corr, args.conftype)
    outputfolder_data = merged_data_path+"/interpolations/"
    outputfolder_plot = lpd.get_plot_path(args.qcdtype, args.corr, args.conftype)+"/interpolations/"
    lpd.create_folder(outputfolder_plot, outputfolder_data)

    # load and set flowradius
    flowradius = numpy.loadtxt(merged_data_path+"/flowradii_"+args.conftype+".dat")[args.flow_index]
    flowradiusstr = '{0:.4f}'.format(flowradius)

    # load data
    XX_samples = numpy.loadtxt(merged_data_path+"/btstrp_samples/"+args.corr+"_"+flowradiusstr+"_Nt"+str(nt)+"_btstrp_samples.dat")

    # result variables
    theoutputdata = []
    interpolation_points = numpy.arange(0, 0.5, 0.001)
    xpoints = numpy.asarray([x for x in numpy.sort(numpy.unique([*interpolation_points, *lpd.get_tauTs(args.int_Nt), 0.5,
                                                                 lpd.lower_tauT_limit_(flowradius, args.max_FlowradiusBytauT, args.max_FlowradiusBytauT_offset)]))
                             if x >= lpd.lower_tauT_limit_(flowradius, args.max_FlowradiusBytauT, args.max_FlowradiusBytauT_offset)])  # xdata[0]

    # this is only relevant for the plot
    tauTs_used_in_int, _ = filter_tauTs(flowradius, lpd.get_tauTs(nt), args.max_FlowradiusBytauT, args.max_FlowradiusBytauT_offset, args.int_left_tauT_offset)
    xpointsplot = numpy.linspace(tauTs_used_in_int[0], xpoints[-1], 1000)
    theplotdata = []

    weights = ((tauTs_used_in_int / flowradius) * args.max_FlowradiusBytauT)**3  # TODO justify the exponent here?

    # perform spline fits for each sample
    for m in range(args.nsamples):
        ydata = numpy.asarray(XX_samples[m*nt_half:(m+1)*nt_half, 1])
        for a in range(nt_half):
            ydata[a] = ydata[a] * lpd.improve_corr_factor(a, nt, args.flow_index, args.use_imp)
        xdata, ydata = filter_corr_data(flowradius, lpd.get_tauTs(nt), ydata, args.max_FlowradiusBytauT, args.max_FlowradiusBytauT_offset, args.int_left_tauT_offset)
        output, plot_output = interpolate_XX_flow(xdata, ydata, weights, xpoints, xpointsplot)
        theoutputdata.append(output)
        theplotdata.append(plot_output)

    # compute spline average over all samples at all desired tauT (xpoints) and save
    ypoints = numpy.mean(theoutputdata, axis=0)
    epoints = numpy.std(theoutputdata, axis=0)
    numpy.savetxt(outputfolder_data+args.corr+"_"+'{0:.4f}'.format(flowradius)+"_interpolation.txt",
                  numpy.stack((xpoints, ypoints, epoints), axis=-1), header="tauT     G_tauF(tau)/Gnorm      err")

    ypointsplot = numpy.mean(theplotdata, axis=0)
    epointsplot = numpy.std(theplotdata, axis=0)

    # plot interpolations and underlying data points
    UseTex = False
    if UseTex:
        displaystyle = r'\displaystyle'
        ylabel = r'$'+displaystyle+r'\frac{ G^\mathrm{latt }_{\tau_F} (\tau)}{G^{\substack{ \text{\tiny  norm} \\[-0.4ex] \text{\tiny latt } } }_{\tau_F = 0} (\tau)}$'
    else:
        ylabel = 'G'
    fig, ax, plots = lpd.create_figure(xlims=[0, 0.51], ylims=[1.4, 4], xlabel=r'$\tau T$', ylabel=ylabel, xlabelpos=(0.95, 0.05), UseTex=UseTex)
    ax.set_title(r'$ \sqrt{8\tau_F}T = $'+'{0:.3f}'.format(flowradius))  # +", nknots = "+nknot_str)
    ax.axvline(x=lpd.lower_tauT_limit_(flowradius, args.max_FlowradiusBytauT, args.max_FlowradiusBytauT_offset), **lpd.verticallinestyle)
    ax.fill_between(xpointsplot, ypointsplot-epointsplot, ypointsplot+epointsplot, alpha=0.5)
    ax.errorbar(xpointsplot, ypointsplot, fmt='-', lw=0.5, mew=0)
    XX = numpy.loadtxt(merged_data_path+"/"+args.corr+"_"+args.conftype+".dat")
    XX_err = numpy.loadtxt(merged_data_path+"/"+args.corr+"_err_"+args.conftype+".dat")
    for a in range(len(XX[args.flow_index])):
        XX[args.flow_index, a] = XX[args.flow_index, a] * lpd.improve_corr_factor(a, nt, args.flow_index, args.use_imp)
        XX_err[args.flow_index, a] = XX_err[args.flow_index, a] * lpd.improve_corr_factor(a, nt, args.flow_index, args.use_imp)
    ax.errorbar(lpd.get_tauTs(nt), XX[args.flow_index], XX_err[args.flow_index], **lpd.plotstyle_add_point_single)
    filepath = outputfolder_plot+"/"+args.corr+"_"+'{0:.4f}'.format(flowradius)+"_interpolation.pdf"
    print("saving ", filepath)
    fig.savefig(filepath)


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
