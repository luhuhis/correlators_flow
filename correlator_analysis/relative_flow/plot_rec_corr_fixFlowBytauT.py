#!/usr/bin/env python3

import lib_process_data as lpd
import numpy
import scipy.integrate
import scipy.interpolate
import scipy.optimize
import argparse
import matplotlib
from spf_reconstruction.model_fitting.spf_reconstruct import Gnorm, Kernel


import warnings
warnings.filterwarnings('ignore', r'.*The maximum number of subdivisions.*')
warnings.filterwarnings('ignore', r'.*The occurrence of roundoff error is detected.*')


def set_title(args, ax):

    title = args.title
    if title != "":
        ax.text(0.5, 1, title, transform=ax.transAxes, ha='center', va='bottom')


def set_text(args, ax):
    if args.custom_text is not None:
        x = float(args.custom_text[0])
        y = float(args.custom_text[1])
        text = args.custom_text[2]
        ha = args.custom_text[3]
        va = args.custom_text[4]
        transform = None
        if args.custom_text[5] == "rel":
            transform = ax.transAxes

        ax.text(x, y, text, ha=ha, va=va, transform=transform)


def set_legend(args, ax):

    # Legend
    for i in range(args.leg_n_dummies):
        ax.errorbar(0, 0, label=' ', markersize=0, alpha=0, lw=0)

    ax.legend()
    handles, labels = ax.get_legend_handles_labels()
    n = len(handles)
    n2 = int(n/2)

    if args.leg_interleave:
        index = numpy.asarray([(i, i+n2) for i in range(0, n2)], dtype=int).flatten()
    else:
        index = range(n)

    legend = ax.legend([handles[i] for i in index], [labels[i] for i in index],
                       loc=args.leg_loc, bbox_to_anchor=args.leg_pos, ncol=args.leg_n_col, **lpd.leg_err_size(), columnspacing=0.5, framealpha=args.leg_framealpha)

    legend.set_title(args.leg_title)


def plot_error_interpolation(xdata, ydata, edata, y_int, e_int, tauT):

    x_int = numpy.linspace(xdata[0], xdata[-1], 1000)

    fig, ax, _ = lpd.create_figure()
    ax.set_xlim(0, tauT/2)
    ax.set_ylim(0, 12)
    ax.errorbar(xdata, ydata, edata, color='k', zorder=10)
    ax.fill_between(x_int, y_int(x_int)-e_int(x_int), y_int(x_int)+e_int(x_int), zorder=9)

    ax.text(0.99, 0.99, lpd.format_float(tauT), ha='right', va='top', transform=ax.transAxes)

    return fig


def plot_relative_flow(args, ax, plot_offset):
    # find correlator at various fixed flowradiusBytauT
    if args.qcdtype is not None and args.corr is not None and args.conftype is not None:
        n = max(len(args.flowradiusBytauT), len(args.qcdtype), len(args.corr), len(args.conftype))
        color_data = args.color_data[plot_offset:plot_offset + n]
        markers = args.markers[plot_offset:plot_offset + n]
        fillstyles = args.fillstyle[plot_offset:plot_offset + n]
        leg_labels = args.leg_labels[plot_offset:plot_offset + n]
        x_scales = args.x_scales[plot_offset:plot_offset + n] if args.x_scales is not None else [1 for _ in range(n)]
        for i in range(n):
            flowradiusBytauT = args.flowradiusBytauT[i] if len(args.flowradiusBytauT) > 1 else args.flowradiusBytauT[0]
            qcdtype = args.qcdtype[i] if len(args.qcdtype) > 1 else args.qcdtype[0]
            corr = args.corr[i] if len(args.corr) > 1 else args.corr[0]
            conftype = args.conftype[i] if len(args.conftype) > 1 else args.conftype[0]
            color = color_data[i]
            marker = markers[i]
            fillstyle = fillstyles[i]
            label_suffix = leg_labels[i]
            x_scale = x_scales[i]
            inputfolder = lpd.get_merged_data_path(qcdtype, corr, conftype, args.basepath)
            these_flowtimes = numpy.loadtxt(inputfolder + "flowtimes_" + conftype + ".dat")
            XX = numpy.loadtxt(inputfolder + corr + "_" + conftype + ".dat")
            XX_err = numpy.loadtxt(inputfolder + corr + "_err_" + conftype + ".dat")
            beta, ns, nt, nt_half = lpd.parse_conftype(conftype)
            these_flowradii = numpy.sqrt(8 * these_flowtimes) / nt
            these_tauT = numpy.arange(1 / nt, 0.501, 1 / nt)

            fermions, temp, flowtype, gaugeaction, flowaction = lpd.parse_qcdtype(qcdtype)

            # interpolate between flow times
            y_int = []
            e_int = []
            for i in range(len(these_tauT)):
                G_latt_LO_flow = lpd.G_latt_LO_flow(i, these_flowtimes, corr, nt, flowaction, gaugeaction)
                ydata = numpy.asarray(XX[:, i]) * nt ** 4 / G_latt_LO_flow
                edata = numpy.asarray(XX_err[:, i]) * nt ** 4 / numpy.fabs(G_latt_LO_flow)
                xdata = these_flowradii
                min_index = numpy.fabs(these_flowradii-args.min_flowradius).argmin()-1
                y_int.append(scipy.interpolate.InterpolatedUnivariateSpline(xdata[min_index:], ydata[min_index:], k=3, ext=2, check_finite=True))
                e_int.append(scipy.interpolate.InterpolatedUnivariateSpline(xdata[min_index:], edata[min_index:], k=3, ext=2, check_finite=True))

            find_and_plot_and_save_relflow(args, ax, flowradiusBytauT, these_tauT, y_int, e_int,
                                           inputfolder, color, marker, fillstyle,
                                           int(-1000 * flowradiusBytauT),
                                           x_scale, label_suffix)


def find_and_plot_and_save_relflow(args, ax, flowradiusBytauT, possible_tauTs, y_int_list, e_int_list, out_file_path, color, marker, fillstyle,
                                   zorder=-10, xscale=1, label=""):

    relflow_edata = []
    relflow_ydata = []
    relflow_xdata = []
    for i, tauT in enumerate(possible_tauTs):
        wanted_flowradius = flowradiusBytauT * tauT
        if y_int_list[i] is not None:
            try:
                relflow_edata.append(e_int_list[i](wanted_flowradius))
                relflow_ydata.append(y_int_list[i](wanted_flowradius))
                relflow_xdata.append(tauT)
            except ValueError:
                # print("Warn: ValueError for flow radius", wanted_flowradius)
                pass
    relflow_xdata = numpy.asarray(relflow_xdata)
    relflow_ydata = numpy.asarray(relflow_ydata)
    relflow_edata = numpy.asarray(relflow_edata)
    flowstr = r'{0:.2f}'.format(flowradiusBytauT)
    lpd.create_folder(out_file_path + "/rel_flow/")
    numpy.savetxt(out_file_path + "/rel_flow/"+args.corr[0]+"_relflow_" + flowstr + ".dat",
                  numpy.column_stack((relflow_xdata, relflow_ydata, relflow_edata)),
                  header="tauT      G/Gnorm_sqrt(8tauF)/tau=" + flowstr + "       err")

    xdata = numpy.asarray(relflow_xdata) * xscale
    ax.errorbar(xdata, relflow_ydata, relflow_edata, zorder=zorder, fmt=marker, color=color, fillstyle=fillstyle,
                label=label, markeredgewidth=0.75, elinewidth=0.75, markersize=4)

    if not args.no_connection:
        ax.errorbar(xdata, relflow_ydata, zorder=-100 * zorder,  markersize=0, alpha=0.5, fmt='-', color=color)


def plot_extrapolated_data(args, ax, plot_offset):
    if args.plot_flow_extr is not None:
        for i, path in enumerate(args.plot_flow_extr):
            try:
                XX_flow_extr = numpy.loadtxt(path, unpack=True)
                ax.errorbar(XX_flow_extr[0]*args.flow_extr_custom_units[i], XX_flow_extr[1], XX_flow_extr[2],
                            zorder=-10000, color=args.color_data[i], fillstyle=args.fillstyle[i], fmt=args.markers[i], label=args.leg_labels[i], markersize=4, markeredgewidth=0.75, elinewidth=0.75)
                if not args.no_connection:
                    ax.errorbar(XX_flow_extr[0], XX_flow_extr[1], zorder=-100000, markersize=0, color=args.color_data[i], fmt='-', alpha=0.5)
                plot_offset += 1
            except OSError:
                print("WARN: skipping tf->0 extr., could not find ", path)
    return plot_offset

#
# # TODO fix this
# def plot_relflow_continuum_data(args, ax, plot_offset):
#
#     for datapath, flowpath in zip(args.plot_cont_extr, args.cont_flowtimes):
#         plot_offset += 1
#         conftype = args.conftype[0]  # TODO fix conftype
#
#         # TODO update this for new sample-based format
#         # load continuum quenched data and interpolate it
#         EE_cont_arr = []
#         flowtimes = numpy.loadtxt(flowpath)
#         flowradii = numpy.sqrt(8*flowtimes)
#         min_idx = 0
#         max_idx = 0
#         stop = False
#
#         # TODO replace this loop over flowradius with one line that loads the npy file containing all flow and all samples.
#         # TODO then compute the mean+std and use that further down the line
#         for flowradius in flowradii:
#             flowradius_str = r'{0:.4f}'.format(flowradius)
#             fitparam_file = path+"/cont_extr/EE_" + flowradius_str + "_cont.txt"
#             try:
#                 EE_cont_arr.append(numpy.loadtxt(fitparam_file, unpack=True))
#                 stop = True
#                 max_idx += 1
#             except OSError:
#                 # print("could not find ", fitparam_file)
#                 if not stop:
#                     min_idx += 1
#         max_idx += min_idx
#         flowradii = flowradii[min_idx:max_idx]  # filter which flowradii got actually read in
#         possible_tauTs = EE_cont_arr[0][0]  # doesn't matter which one we choose here, they are all the same
#
#         # interpolate between flowtimes for each tauT. for that we first need to remove the nan's
#         cont_int = []
#         cont_err_int = []
#         for i in range(len(possible_tauTs)):
#             ydata = numpy.asarray([EE_cont[1][i] for EE_cont in EE_cont_arr])
#             edata = numpy.asarray([EE_cont[2][i] for EE_cont in EE_cont_arr])
#             mask = ~numpy.isnan(ydata)
#             xdata = flowradii[mask]
#             ydata = ydata[mask]
#             edata = edata[mask]
#             if len(xdata) >= 4:
#                 cont_int.append(scipy.interpolate.InterpolatedUnivariateSpline(xdata, ydata, k=3, ext=2, check_finite=True))
#                 cont_err_int.append(scipy.interpolate.InterpolatedUnivariateSpline(xdata, edata, k=3, ext=2, check_finite=True))
#             else:
#                 cont_int.append(None)
#                 cont_err_int.append(None)
#         outfolder = path+"/cont_rel_flow/"
#         lpd.create_folder(outfolder)
#
#         find_and_plot_and_save_relflow(args, ax, args.flowradiusBytauT[0], possible_tauTs, cont_int, cont_err_int, outfolder, args.color_data[1], r'a\rightarrow 0',
#                                        None, None, zorder=-9999)
#
#     return plot_offset


def plot_fit_corrs(args, ax, xpoints, model_means, just_UVs):
    zorders = range(len(model_means))

    # plot reconstructed correlators
    for i in zorders:
        x_scale = args.x_scales[i] if args.x_scales is not None else 1
        these_xpoints = xpoints * x_scale
        model_mean = model_means[i]
        just_UV = just_UVs[i]
        color = args.color_fit[i]
        zorder = i

        ax.errorbar(these_xpoints, model_mean, color=color, lw=0.75, fmt='--', zorder=-10000+zorder)

        if args.show_UV_corrs:
            ax.errorbar(these_xpoints, just_UV, color=color, alpha=0.25, fmt='-.', zorder=-10000+zorder)


def Integrand(OmegaByT, tauT, SpfByT3):
    return 1. / numpy.pi * Kernel(OmegaByT, tauT) * SpfByT3(OmegaByT)


def calc_Gmodel(tauT, spf, MinOmegaByT, MaxOmegaByT):
    normalization = Gnorm(tauT)
    Gmodel = scipy.integrate.quad(Integrand, MinOmegaByT, MaxOmegaByT, args=(tauT, spf), limit=100)[0]
    return Gmodel / normalization


def plot_fits(args, ax):

    if args.fit_folders is not None:

        xpoints = numpy.linspace(0.25, 0.5, args.npoints)
        Gmodel_by_norm_arr = []
        GmodelUV_by_norm_arr = []

        file_prefixes = [args.fit_basepath+folder for folder in args.fit_folders]

        # calculate fit corrs for all given models
        for file_prefix in file_prefixes:
            # load PhiUV and calculate model correlators
            PhiUVByT3 = numpy.load(file_prefix + "/phiUV.npy")
            spfByT3 = numpy.load(file_prefix + "/spffit.npy")

            PhiUVByT3_spline = scipy.interpolate.InterpolatedUnivariateSpline(PhiUVByT3[:, 0], PhiUVByT3[:, 1], k=3, ext=2, check_finite=True)
            spfByT3_spline = scipy.interpolate.InterpolatedUnivariateSpline(spfByT3[:, 0], spfByT3[:, 1], k=3, ext=2, check_finite=True)
            MinOmegaByT = spfByT3[0, 0]
            MaxOmegaByT = spfByT3[-1, 0]

            Gmodel_by_norm = lpd.parallel_function_eval(calc_Gmodel, xpoints, 20, spfByT3_spline, MinOmegaByT, MaxOmegaByT)
            GmodelUV_by_norm = lpd.parallel_function_eval(calc_Gmodel, xpoints, 20, PhiUVByT3_spline, MinOmegaByT, MaxOmegaByT)

            Gmodel_by_norm_arr.append(Gmodel_by_norm)
            GmodelUV_by_norm_arr.append(GmodelUV_by_norm)

        plot_fit_corrs(args, ax, xpoints, Gmodel_by_norm_arr, GmodelUV_by_norm_arr)

    return


def prepare_plot_canvas(args):

    xlabel = args.xlabel
    subscript = r''
    if args.corr == ['EE']:
        subscript = r'_E'
    elif args.corr == ['BB']:
        subscript = r'_B'

    ylabel = r'$\dfrac{G'+subscript+r'}{G^\mathrm{norm}}$'

    fig, ax, _ = lpd.create_figure(figsize=args.figsize, UseTex=args.usetex, xlims=args.xlims, xlabel=xlabel, ylims=args.ylims, ylabel=ylabel)

    # plot dashed vlines
    if not args.no_vlines and args.tauT_vlines is not None and args.color_vlines is not None:
        for tauT, color in zip(args.tauT_vlines, args.color_vlines):
            ax.axvline(tauT, **lpd.chmap(lpd.verticallinestyle, color=color, alpha=1))

    if args.xticks != "auto":
        ax.set_xticks([float(i) for i in args.xticks])

    return fig, ax


def parse_args():

    def set_output_path(outputpath, qcdtype, corr):
        if not outputpath:
            if qcdtype:
                outputpath = lpd.get_plot_path(qcdtype[0], "", "") + "/rel_flow/"
            if qcdtype and corr:
                outputpath = lpd.get_plot_path(qcdtype[0], corr[0], "") + "/rel_flow/"
        lpd.create_folder(outputpath)
        return outputpath

    parser = argparse.ArgumentParser()

    parser.add_argument('--basepath', default="../../../../data/merged/")
    parser.add_argument('--plot_flow_extr', help='paths to files containing flow-extr data.', default=None, type=str, nargs='*')
    parser.add_argument('--flow_extr_custom_units', nargs='*', type=float)
    parser.add_argument('--plot_cont_extr', help='paths to files containing continuum-extr data.', default=None, type=str, nargs='*')
    parser.add_argument('--cont_flowtimes', help='paths to files containing continuum flowtimes in temperature units', default=None, type=str, nargs='*')

    parser.add_argument('--plot_quenched_systematics', action="store_true")  # TODO soften hard code

    parser.add_argument('--flowradiusBytauT', nargs='*', type=float, help='which sqrt(8tauF)/tauT to plot')
    parser.add_argument('--min_flowradius', default=0.075, type=float, help="set to at least one lattice spacing of the coarsest lattice.")

    # fitparam options
    parser.add_argument('--fit_folders', nargs='*', help='plot model 3 spf for these fit params')
    parser.add_argument('--fit_basepath', default="", type=str, help='shared basepath of all fitparam files, prepended to the folders.')

    # general plot options
    parser.add_argument('--usetex', action="store_true")
    parser.add_argument('--ylims', help='custom ylims', nargs=2, type=float, default=(-0.03, 4))
    parser.add_argument('--xlims', help='custom xlims', nargs=2, type=float, default=(-0.01, 0.55))
    parser.add_argument('--xticks', default="auto", nargs="*")
    parser.add_argument('--title', help='title prefix of plot', default="", type=str)
    parser.add_argument('--custom_text', help="show some custom text. modify this script to customize it. (i dont feel like adding 10 more parameters "
                                              "just for some text)", nargs='*')
    parser.add_argument('--xlabel', default=r'$\tau T$')

    parser.add_argument('--x_scales', type=float, nargs='*')
    parser.add_argument('--min_tauT', type=float, help='lower limits from which on the corr is fitted. appear in the kappa plot.', nargs='*')
    # parser.add_argument('--min_tauT_plot', type=float, help='lower limit from which tauT on the reconstructed corr is valid/should be plotted. '
    #                                                         'default is to parse it from fitparam file.', nargs='*', default=[0.05,])
    parser.add_argument('--tauT_vlines', type=float, help='where to plot vertical dashed lines', nargs='*')
    parser.add_argument('--npoints', help='how many tauT to plot between 0 and 0.5 for the integrated spfs (model correlators)', default=100, type=int)

    # legend options
    parser.add_argument('--leg_pos', nargs=2, default=(0.02, 1), type=float, help="where to put the upper left corner of the legend")
    parser.add_argument('--leg_loc', type=str, default='lower right', help="anchor point of the legend")
    parser.add_argument('--leg_n_dummies', help="how many dummy entries to put in the legend in order to align things nicely", default=0, type=int)
    parser.add_argument('--leg_n_col', help="how many columns the legend should have", default=1, type=int)
    parser.add_argument('--leg_title', help="Put this at end of legend title", type=str, default="")
    parser.add_argument('--leg_labels', help="additional strings to put after the Ntau in the labels in legend", nargs='*', type=str)
    parser.add_argument('--leg_interleave', help="whether to interleave the given data in the legend. this is useful if you want to pair extrapolated data with nonextrapolated data.",
                        default=False, action="store_true")
    parser.add_argument('--leg_framealpha', default=0.8, type=int)

    # plot styling
    parser.add_argument('--figsize', help="size of the figure", default=None, nargs='*')
    parser.add_argument('--markers', nargs='*', default='|', type=str, help="fmt keywords for the plots")
    parser.add_argument('--fillstyle', nargs='*', default='none', type=str, help="fillstyle keywords for the plots")
    parser.add_argument('--color_data', help="specify colors in order that the data should have", type=str, nargs='*',
                        default=('C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9'))
    parser.add_argument('--color_fit', help="specify colors in order that the fits should have", type=str, nargs='*',
                        default=('C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9'))
    parser.add_argument('--color_vlines', help="specify colors in order that the vlines should have", type=str, nargs='*',
                        default=('C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9'))

    # hide certain plot elements
    parser.add_argument('--no_vlines', help="hide vertical lines that indicate lower fit limit", action="store_true")
    parser.add_argument('--no_connection', help="hide connections between data points", action="store_true")
    parser.add_argument('--no_label', help="hide labels of fit corrs", action="store_true")
    parser.add_argument('--show_UV_corrs', action="store_true", help="plot UV part of model spf translated to correlator domain")
    parser.add_argument('--hide_fits', help="hide fits even if they were provided", action="store_true")

    # when you want to load finite lattice spacing
    parser.add_argument('--qcdtype', help="list of qcdtypes. format doesnt matter, only used for finding data. example: quenched_1.50Tc_zeuthenFlow", nargs='*',
                        type=str)
    parser.add_argument('--corr', help="list of corrs. choose from EE, EE_clover, BB, BB_clover", type=str, nargs='*')
    parser.add_argument('--conftype', nargs='*', type=str, help="list of conftypes, e.g. s080t20_b0703500 s096t24_b0719200 s120t30_b0739400 s144t36_b0754400")

    # file output options
    parser.add_argument('--output_path', help="where to save the plot.")
    parser.add_argument('--output_suffix', type=str, help="suffix to add to plot file output name")

    args = parser.parse_args()

    # check if arguments make sense
    if args.plot_cont_extr is not None:
        if not all(element == args.flowradiusBytauT[0] for element in args.flowradiusBytauT):
            # TODO remove this limitation
            print("ERROR: when plotting the quenched cont/flow extr, all flowradiusBytauT have to be identical")
            exit(1)

    if args.conftype is not None and args.flowradiusBytauT is None:
        print("ERROR: need flowradiusBytauT")
        exit(1)
    args.output_path = set_output_path(args.output_path, args.qcdtype, args.corr)

    return args


def main():

    # TODO load flowtime instead of flowradius file

    args = parse_args()

    fig, ax = prepare_plot_canvas(args)

    if not args.hide_fits:
        plot_fits(args, ax)

    plot_offset = plot_extrapolated_data(args, ax, plot_offset=0)
    # we return an offset, so that plot_relative_flow knows that something else has already been plotted, so that it can skip already used colors.

    # plot_offset = plot_relflow_continuum_data(args, ax, plot_offset)

    plot_relative_flow(args, ax, plot_offset)

    set_legend(args, ax)
    set_text(args, ax)
    set_title(args, ax)

    # save figure
    file = args.output_path + "/"+args.corr[0]+"_relflow" + args.output_suffix + ".pdf"
    print("saving ", file)
    fig.savefig(file)
    matplotlib.pyplot.close(fig)


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
