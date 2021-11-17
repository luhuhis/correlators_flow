#!/usr/local/bin/python

import lib_process_data as lpd
import numpy
import sys
import matplotlib
from matplotlib import pyplot, legend_handler, container
from matplotlib.ticker import MaxNLocator

# Global variables for this script
figsize = None
corr_label = ""


# input
#   inputfolder+"flowradii_"+conftype+".dat"
#   inputfolder+"XX_"+conftype+".dat"
#   inputfolder+"XX_err_"+conftype+".dat"
# output
#   outputfolder+conftype+"_XX_fill.pdf"
#   outputfolder+"/"+conftype+"_XX_flowlim_fill.pdf"


def skip_large_errors(XX, XX_err, boundary):
    for i in range(0, XX.shape[0]):
        for j in range(0, XX.shape[1]):
            if XX_err[i, j] > boundary:
                XX[i, j] = None


def plot1(XX, XX_err, args, prefix, flow_selection, flow_var, xdata_plot, nt_half, flow_str, outputfolder):
    # default options for all plots
    global corr_label
    global figsize

    # default options for this plot
    xlims = [0, 0.505]
    xlabel = r'$\tau T$'
    ylabelpos = (0.01, 0.97)
    xlabelpos = (0.95, 0.06)
    ylims = ([1, 100000])
    ylabel = r'$\displaystyle\frac{G^\mathrm{latt }}{T_\mathrm{norm}^4 }$'

    # special options from user input
    if args.show_TauByA:
        xlims = [0, nt_half * 1.01]
        xlabel = r'$\tau/a$'
    if prefix == "_numerator":
        ylims = [-50, 17]
        ylabel = r'$\displaystyle\frac{ \langle ' + corr_label + corr_label + r' \rangle}{T_\mathrm{norm}^4 }$'
    if prefix == "_polyakovloop":
        ylims = [1, 1000000]
        ylabel = r'$\displaystyle\frac{ \langle -- \rangle}{T_\mathrm{norm}^4 }$'
    if args.wideaspect:
        ylabelpos = (-0.2, 0.95)
        xlabelpos = (0.5, -0.06)
    if args.reconstruct:
        ylabel = r'$\displaystyle\frac{G^\mathrm{rec}}{T_\mathrm{norm}^4 }$'
    if args.custom_ylims:
        ylims = args.custom_ylims

    # create the figure and axes objects
    fig, ax, plots = lpd.create_figure(xlims=xlims, ylims=ylims, xlabel=xlabel, xlabelpos=xlabelpos, ylabel=ylabel,
                                       ylabelpos=ylabelpos, UseTex=True, figsize=figsize)

    # modify axes if necessary
    ax.set_yscale('log', nonposy='mask')
    if prefix == "_numerator":
        ax.set_yscale('linear')

    # draw the data
    for i in flow_selection:
        mycolor = lpd.get_color(flow_var, i, flow_selection[0], flow_selection[-1])
        ax.errorbar(list(xdata_plot[:nt_half]), XX[i, :nt_half], color=mycolor, **lpd.plotstyle_lines, zorder=-flow_selection[-1] + i + 1)
        plots.append(ax.errorbar(list(xdata_plot[:nt_half]), XX[i, :nt_half], XX_err[i, :nt_half], color=mycolor, fmt='D', linewidth=0.5,
                                 mew=0.25, mec='grey', markersize=1.5, capsize=3, zorder=-flow_selection[-1] + i + 2))

    # draw the legend
    leg = ax.legend(handles=plots, labels=['{0:.3f}'.format(j) for i, j in enumerate(flow_var) if i in flow_selection], title=flow_str, **lpd.legendstyle)

    ax.set_title(r'$N_\tau=' + str(int(nt_half * 2)) + '$', x=0.5, y=0.9)

    # change first xtick label from '0.0' to '0'
    fig.canvas.draw()
    ticks = ax.get_xticks().tolist()
    ticks = ['{0:.1f}'.format(x) for x in ticks]
    ticks[0] = '0'
    ax.set_xticklabels(ticks)

    # save pdf
    filename = outputfolder + args.conftype + "_" + args.corr + prefix + "_T.pdf"
    fig.savefig(filename)
    print("SAVED correlator plot \n", filename)

    # clear canvas
    ax.lines.clear()
    ax.collections.clear()
    plots.clear()


def plot2(XX, XX_err, args, prefix, flow_selection, flow_var, xdata_plot, nt_half, flow_str, outputfolder,
          fermions):
    # default options for all plots
    global corr_label
    global figsize

    # default options for this plot
    xlabel = r'$\tau T$'
    xlabelpos = (0.95, 0.06)
    ylabelpos = (0.01, 0.97)
    xlims = [0, 0.505]
    ylims = ([-1, 20])
    ylabel = r'$\displaystyle\frac{G^\mathrm{latt }}{G_{\tau_\mathrm{F}=0}^{\begin{subarray}{l}\mathrlap{\textrm{\tiny  norm}}\\[-0.3ex] \textrm{\tiny latt}\end{subarray}}}$'
    num_format = '{0:.3f}'
    tau_format = '{0:.2f}'

    # special options from user input
    if args.show_TauByA:
        xlims = [0, nt_half * 1.01]
        xlabel = r'$\tau/a$'
        tau_format = '{0:.0f}'
    if args.show_flowtime_instead_of_flowradii:
        num_format = '{0:.1f}'
    if fermions == "hisq":
        ylims = ([-0.5, 10])
        # if corr == "BB_clover" or corr == "BB":
        # ylims = ([-0.5,6])
    if fermions == "quenched":
        ylims = ([0, 4])
    if args.reconstruct:
        ylabel = r'$\displaystyle\frac{G^\mathrm{rec}_{T\rightarrow 2T}}{G_{\tau_\mathrm{F}=0}^{\begin{subarray}{l}\mathrlap{\textrm{\tiny  norm}}\\[-0.3ex] \textrm{\tiny latt}\end{subarray}}}$'
    if args.part_obs == "numerator":
        ylims = [-50, 17]
        ylabel = r'$\displaystyle  \langle ' + corr_label + corr_label + r' \rangle$'
    if args.part_obs == "polyakovloop":
        ylims = [1, 1000000]
        ylabel = r'$\displaystyle\frac{ \langle - \rangle}{G_{\tau_\mathrm{F}=0}^{\begin{subarray}{l}\mathrlap{\textrm{\tiny  norm}}\\[-0.3ex] \textrm{\tiny latt}\end{subarray}}}$'

    # if args.wideaspect:
        # ylabelpos = (-0.2, 0.95)
        # xlabelpos = (0.5, -0.06)
    if args.custom_ylims:
        ylims = args.custom_ylims
    if args.custom_xlims2:
        xlims = args.custom_xlims2

    # create the figure and axes objects
    fig, ax, plots = lpd.create_figure(xlims=xlims, ylims=ylims, xlabel=xlabel, xlabelpos=xlabelpos,
                                       ylabelpos=ylabelpos,
                                       ylabel=ylabel,
                                       UseTex=True, figsize=figsize)

    if args.conftype_2:
        ax.set_yscale('log', nonposy='mask')

    # draw the data
    for i in flow_selection:
        mycolor = lpd.get_color(flow_var, flow_selection[-1]-i+flow_selection[0], flow_selection[0], flow_selection[-1])
        plots.append(ax.fill_between(list(xdata_plot[:nt_half]), XX[i, :nt_half] - XX_err[i, :nt_half],
                                     XX[i, :nt_half] + XX_err[i, :nt_half], facecolor=mycolor, lw=lpd.mylinewidth, zorder=-flow_selection[-1] + i))
        ax.errorbar(list(xdata_plot[:nt_half]), XX[i, :nt_half], color=mycolor, **lpd.plotstyle_lines, zorder=-flow_selection[-1] + i + 1)
        ax.errorbar(list(xdata_plot[:nt_half]), XX[i, :nt_half], color=mycolor, fmt='D', linewidth=0.5, mew=0.25, mec='grey', markersize=1.5,
                    capsize=3, zorder=-flow_selection[-1] + i + 2)

    # draw the legend
    leg = ax.legend(handles=plots, labels=[num_format.format(j) for i, j in enumerate(flow_var) if i in flow_selection], title=flow_str, **lpd.legendstyle)

    if args.reconstruct:
        ax.set_title(r'$N_\tau=' + str(int(nt_half * 4)) + r'\xrightarrow{\mathrm{rec}}' + str(int(nt_half * 2)) + '$', x=0.5, y=0.85)
        if args.conftype_2:
            ax.set_title(r'\scriptsize [errors of $\langle--\rangle$ ignored!]', x=0.5, y=0.89)
    else:
        ax.set_title(r'$N_\tau=' + str(int(nt_half * 2)) + '$', x=0.5, y=0.89)

    # change first xtick label from '0.0' to '0'
    fig.canvas.draw()
    ticks = ax.get_xticks().tolist()
    ticks = [tau_format.format(x) for x in ticks]
    ticks[0] = '0'
    ax.set_xticklabels(ticks)

    # save pdf
    filename = outputfolder + args.conftype + "_" + args.corr + prefix + ".pdf"
    fig.savefig(filename)
    print("SAVED correlator plot\n", filename)

    # clear canvas
    ax.lines.clear()
    ax.collections.clear()
    plots.clear()


def plot3(XX, XX_err, args, prefix, flow_var, xdata_plot, nt_half: int, outputfolder: str, fermions: str, valid_flowtimes, tauT, flow_radius, flowend: int):
    # default options for all plots
    global corr_label
    global figsize

    # default options for this plot
    ylabelpos = (0.01, 0.97)
    xlabelpos = (0.94, 0.05)
    xlabel = r'$ \tau_\mathrm{F} T^2 $'
    xlabel2 = r'$\sqrt{8\tau_\mathrm{F}}T$'
    legend_title = r"$\tau T=$"
    xlims = [-0.0001, flow_var[-1] ** 2 / 8 * 1.01]
    xlabel2_format = "%.2f"
    tau_format = '{0:.3f}'
    ylabel = r'$\displaystyle\frac{G^\mathrm{latt }}{G_{{\tau_\mathrm{F}}=0}^{\substack{ \text{\tiny  norm} \\[-0.4ex] \text{\tiny latt } } } }$'

    # special options from user input
    # if args.wideaspect:
        # ylabelpos = (-0.2, 0.95)
        # xlabelpos = (0.5, -0.06)
    if fermions == "hisq":
        ylims = (-0.5, 10)
        ylabelpos = (0.05, 0.97)
        xlabelpos = (0.94, 0.1)
    if fermions == "quenched":
        ylims = (0, 4)
        xlims = (-0.0001, 0.003)
    if args.show_flowtime_instead_of_flowradii:
        xlabel = r'$ \tau_\mathrm{F} /a^2 $'
        xlabel2 = r'$\sqrt{8\tau_\mathrm{F}}/a$'
        if valid_flowtimes is not None:
            xlims = [-0.05, valid_flowtimes[-1] * 1.01]
        xlabel2_format = "%.1f"
    if args.show_TauByA:
        legend_title = r"$\tau/a =$"
        tau_format = '{0:.0f}'
    if args.reconstruct:
        ylabel = r'$\displaystyle\frac{G^\mathrm{rec}_{T\rightarrow 2T} }{G_{\tau_\mathrm{F}=0}^{\begin{subarray}{l}\mathrlap{\textrm{\tiny  norm}}\\[-0.3ex] \textrm{\tiny latt}\end{subarray}}   }$'
    if args.conftype_2:
        ylabel = r'$\displaystyle \frac { \langle E-E \rangle _{' + str(nt_half * 2) + r'} }' \
                                                                                       r'{\langle E-E \rangle_{' + str(nt_half * 4) + r' }^{\mathrm{rec} }}' \
                                                                                                                                      r'\times \scriptstyle \frac{ \langle -- \rangle_{' + str(
            nt_half * 4) + r'}}' \
                           r'{\langle -- \rangle_{' + str(nt_half * 2) + r'} }$'
    if args.part_obs == "numerator":
        ylims = [-2, 17]
        ylabelpos = (0.1, 0.97)
        ylabel = r'$\displaystyle{ \langle ' + corr_label + corr_label + r' \rangle} $'
    if args.part_obs == "polyakovloop":
        ylims = [1e-6, 1]
        ylabel = r'$\displaystyle{ \langle --- \rangle}$'
    # if args.wideaspect:
        # ylabelpos = (-0.2, 0.95)
        # xlabelpos = (0.5, -0.06)
    if args.custom_ylims:
        ylims = args.custom_ylims
    if args.custom_xlims3:
        xlims = args.custom_xlims3

    # create the figure and axes objects
    fig, ax, plots = lpd.create_figure(xlims=xlims, ylims=ylims, xlabelpos=xlabelpos, ylabelpos=ylabelpos,
                                       xlabel=xlabel, ylabel=ylabel, UseTex=True,
                                       figsize=figsize)

    # modify axes if necessary
    if prefix == "_polyakovloop":
        ax.set_yscale('log', nonposy='mask')

    # calculate positions on the curves for perturbative flow limits
    # if not args.ignore_pert_lims:
    flow_limit = numpy.zeros(len(tauT))
    interpolation_value = numpy.zeros(len(tauT))
    for i in range(nt_half):
        flow_limit[i] = lpd.upper_flow_limit(tauT[i])
        index = (numpy.abs(flow_radius[:] - flow_limit[i])).argmin()
        if index == 0:
            interpolation_value[i] = XX[index, i]
            continue
        offset = 1 if flow_limit[i] - flow_radius[index] > 0 else -1
        if (index + offset) >= len(flow_radius):
            interpolation_value[i] = numpy.nan
            continue
        offset2 = -1 if offset == -1 else 0
        interpolation_value[i] = XX[index + offset2, i] - abs(XX[index, i] - XX[index + offset, i]) / abs(
            flow_radius[index] - flow_radius[index + offset]) * abs(flow_limit[i] - flow_radius[index + offset2])

    tau_selection = numpy.arange(nt_half)
    if args.tau_selection:
        tau_selection = tau_selection[args.tau_selection]

    for idx, i in enumerate(tau_selection):
        if args.part_obs == "polyakovloop" and idx != 0:
            continue
        flow_extr_filter_high = 0  # determines up to which index the discrete error bars should be plotted for this tauT
        flow_extr_filter_low = -1  # determines from which index on the discrete error bars should be plotted
        flow_extr_filter_low_grey = -1  # determines from which index on the grey band should be plotted
        flow_extr_filter_low_greyline = -1  # determines from which index on the grey line should be plotted. this is serves as a "minimum line width" or as an indicator if the errors are too large.

        rel_err = numpy.fabs(XX_err[:, i] / XX[:, i])
        # loop over flow times
        for j, val in enumerate(rel_err):

            # only show discrete error bars up to the perturbative flow limits
            if not args.ignore_pert_lims:
                if flow_extr_filter_high == 0:
                    if flow_limit[i] < flow_radius[j]:
                        flow_extr_filter_high = j

            # if args.conftype_2:
            #     # print(j, len(rel_err) - 1 - j, rel_err[len(rel_err) - 1 - j], args.rel_err_limit)
            #     if flow_extr_filter_high == 0:
            #         if rel_err[len(rel_err)-1-j] <= args.rel_err_limit:
            #             flow_extr_filter_high = len(rel_err)-1-j

            # only show discrete error bars from here on
            if flow_extr_filter_low == -1:
                if fermions == "hisq" and val <= args.rel_err_limit or fermions == "quenched" and val <= args.rel_err_limit:  # * tauT[i]
                    flow_extr_filter_low = max(j, args.min_trusted_flow_idx)

            # only show grey line with fixed width from here on
            if flow_extr_filter_low_greyline == -1:
                if fermions == "hisq" and val <= 1 or fermions == "quenched" and val <= args.rel_err_limit * tauT[i]:
                    flow_extr_filter_low_greyline = max(j, args.min_trusted_flow_idx)

            # only show grey error band from here on
            if flow_extr_filter_low_grey == -1:
                if val <= 0.10:
                    flow_extr_filter_low_grey = j

        if args.ignore_pert_lims:
            flow_extr_filter_high = flowend
        else:
            flow_extr_filter_high = min(flow_extr_filter_high, flowend)

        if args.flowend:
            flow_extr_filter_high = args.flowend
        if args.flowstart:
            flow_extr_filter_low = args.flowstart

        shift_factor = 0.9975 if not (idx % 2 == 0) else 1.0025

        xdata = flow_var[flow_extr_filter_low:flow_extr_filter_high] ** 2 / 8 * shift_factor
        ydata = XX[flow_extr_filter_low:flow_extr_filter_high, i]
        edata = XX_err[flow_extr_filter_low:flow_extr_filter_high, i]

        # --- explicit colorful datapoints ---
        zorder = -100 * (len(tau_selection)) + i
        color = lpd.get_color(tau_selection, len(tau_selection) - 1 - idx, 0, -1) if not args.part_obs == "polyakovloop" else 'black'

        # --- light grey lines and error bands ---
        if fermions == "quenched":
            plots.append(ax.errorbar(xdata, ydata, edata, color=color, zorder=zorder,
                                     **lpd.chmap(lpd.plotstyle_add_point, fmt='x-', markersize=1.25, capsize=0.6)))

            LOW = flow_extr_filter_low_grey
            xdata = flow_var[LOW:] ** 2 / 8
            ydata = XX[LOW:, i]
            edata = XX_err[LOW:, i]

            # light grey backgrounds
            # ax.fill_between(xdata, ydata - edata, ydata + edata, linewidth=0.5, color=(0.9, 0.9, 0.9, 1),
            #                 zorder=zorder)  #
            # light grey lines for minimum width of background!
            # ax.errorbar(xdata, ydata, color='lightgrey', zorder=zorder, fmt='-', lw=0.5, markersize=0.5)

            # dark grey lines (0.7,0.7,0.7,1)
            # ax.errorbar(flow_var ** 2 / 8, XX[:, i], linewidth=0.5, color=color,
            #             zorder=zorder / 100 + 1, fmt='-')  # , label=tau_format.format(xdata_plot[i]))
            # perturbative flow limit markers
            plots.append(ax.errorbar(flow_limit[i] ** 2 / 8, interpolation_value[i],
                                     marker=lpd.markers[i % len(lpd.markers)], fillstyle='none', markersize=4,
                                     mew=0.25, color=color, label=tau_format.format(tauT[i]),
                                     capsize=0, lw=0))

        # dark grey lines that are cut off, may break perturbative flow limit markes (?)
        if fermions == "hisq":
            plots.append(ax.errorbar(xdata, ydata, edata, color=color, zorder=zorder, label=tau_format.format(xdata_plot[i]),
                                     **lpd.chmap(lpd.plotstyle_add_point, fmt='x', markersize=1.25, capsize=0.6)))

            ax.errorbar(flow_limit[i] ** 2 / 8, interpolation_value[i],
                                     marker=lpd.markers[i % len(lpd.markers)], fillstyle='none', markersize=4,
                                     mew=0.25, color=color, capsize=0, lw=0)

            xdata = flow_var[flow_extr_filter_low:] ** 2 / 8  # FIXME reset to flow_extr_filter_low_greyline
            ydata = XX[flow_extr_filter_low:, i]  # FIXME reset to flow_extr_filter_low_greyline

            # colored lines
            ax.errorbar(xdata, ydata, linewidth=0.5, color=color, alpha=0.5, zorder=zorder)

        # vertical dashed line
        if args.min_trusted_flow_idx > 0:
            ax.axvline(x=flow_var[args.min_trusted_flow_idx] ** 2 / 8, **lpd.verticallinestyle)

    ax.axhline(y=0, **lpd.horizontallinestyle)
    if args.conftype_2:
        ax.axhline(y=1, **lpd.horizontallinestyle)

    if args.part_obs == "numerator":
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))

    # draw second axes
    ax2 = ax.twiny()

    def tick_function(X):
        V = numpy.sqrt(X * 8)
        return ["%.1f" % z if z == 0 else xlabel2_format % z for z in V]

    ax2.tick_params(direction='in', pad=0, width=0.5)
    ax2.set_xlim(ax.get_xlim())
    nticks = len(ax.xaxis.get_ticklabels())
    new_tick_locations = numpy.linspace(0, xlims[1], nticks)
    ax2.set_xticks(new_tick_locations)  # TODO make this a cmd line input
    ax2.set_xticklabels(tick_function(new_tick_locations))
    ax2.set_xlabel(xlabel2)  # , horizontalalignment='right', verticalalignment='top', bbox=lpd.labelboxstyle, zorder=99)
    ax2.xaxis.set_label_coords(0.92, 0.94)

    # draw legend
    if not args.part_obs == "polyakovloop":
        xerr_size = 0.3 if fermions == "hisq" else 0
        ax.legend(handles=plots)
        handles, labels = ax.get_legend_handles_labels()
        ncol = 2 if not args.wideaspect else 1
        bbox_to_anchor = (1, 0.1) if not args.wideaspect else (1, 1)
        loc = "lower right" if not args.wideaspect else "upper left"
        fontsize = lpd.fontsize if not args.wideaspect else 8
        leg = ax.legend(handles[::-1], labels[::-1], ncol=ncol, title=legend_title, fontsize=fontsize,
                        **lpd.chmap(lpd.legendstyle, loc=loc, bbox_to_anchor=bbox_to_anchor, columnspacing=0.3, handletextpad=0.2, handlelength=0.5,
                                    handler_map={matplotlib.container.ErrorbarContainer: matplotlib.legend_handler.HandlerErrorbar(xerr_size=xerr_size)}))
    # for line in leg.get_lines():
    # line.set_linewidth(4.0)

    if args.reconstruct:
        ax.set_title(r'$N_\tau=' + str(int(nt_half * 4)) + r'\xrightarrow{\mathrm{rec}}' + str(int(nt_half * 2)) + '$', x=0.5, y=0.85)
        if args.conftype_2:
            ax.set_title(r'\scriptsize [errors of $\langle--\rangle$ ignored!]', x=0.5, y=0.89)
    else:
        ax.set_title(r'$N_\tau=' + str(int(nt_half * 2)) + '$', x=0.5, y=0.89)

    # save pdf
    filename = outputfolder + "/" + args.conftype + "_" + args.corr + prefix + "_flow.pdf"
    fig.savefig(filename)
    print("SAVED flow effect plot\n", filename)

    # clear canvas
    ax.lines.clear()
    ax.collections.clear()
    plots.clear()


def get_args():
    """ parse cmd line arguments """
    parser, requiredNamed = lpd.get_parser()

    parser.add_argument('--tau_selection', type=int, nargs='*',
                        help='list of indices which tauT to plot in third flow plot')

    parser.add_argument('--flowselectionfile', help="only consider the flowtimes given in this file")
    parser.add_argument('--flowend', type=int, help="index of the maximum flow time")
    parser.add_argument('--flowstart', type=int, help="index of the minimum flow time", default='0')
    parser.add_argument('--flowstep', type=int, help="index step size for flow time", default='10')

    parser.add_argument('--part_obs', help='choose to only show numerator or denominator (polyakovloop) of corr',
                        choices=['numerator', 'polyakovloop'])

    parser.add_argument('--reconstruct', action="store_true",
                        help="instead of plotting the actual correlator data of <conftype> at <nt>, use that data to plot the reconstructed correlator corresponding to 2x<nt>")
    parser.add_argument('--conftype_2',
                        help="<conftype_2> should have half the <nt> of <conftype>. if --reconstruct is given: load data from this <conftype_2> and plot the ratio (<conftype_2>/<conftype>) of the numerators WITH errors and normalize it with the ratio of the polyakovloops WITHOUT errors.")

    parser.add_argument('--show_flowtime_instead_of_flowradii', action="store_true")
    parser.add_argument('--show_TauByA', help="show tau/a instead of tau*T", action="store_true")

    parser.add_argument('--wideaspect', action="store_true",
                        help="use a wide aspect ratio (basically make the plots larger)")
    parser.add_argument('--custom_ylims', nargs=2, type=float, help="overwrite any default ylims")

    parser.add_argument('--custom_xlims2', nargs=2, type=float, help="overwrite any default xlims for plot2")
    parser.add_argument('--custom_xlims3', nargs=2, type=float, help="overwrite any default xlims for plot3")

    parser.add_argument('--error_reduce_factor', type=float, default='1.0', help="reduce errors of all data by this factor")
    parser.add_argument('--rel_err_limit', type=float, default=0.05,
                        help="float in [0,inf]. for plot no3: only show error bars if relative error is less than this value. default value is 0.05 (5% relative error)")

    parser.add_argument('--hide_err_above', type=float, default=2,
                        help="hide all data with relative errors larger than this. default is just above rel_err_limit.")

    parser.add_argument('--ignore_pert_lims', action="store_true", help="for third plot ignore the perturbative flow limits for deciding the plot range")

    parser.add_argument('--min_trusted_flow_idx', type=int, default=0,
                        help="minimum trusted flowtime index. Set this to the index where the flow radius is larger than one lattice spacing of the coarsest lattice used in cont extr. only data for flow times after this are plotted.")

    parser.add_argument('--only_plot_no', type=int, nargs='*', default=[2, 3])

    requiredNamed.add_argument('--conftype', help="format: s064t64_b0824900_m002022_m01011", required=True)

    args = parser.parse_args()

    if 3 in args.only_plot_no and args.rel_err_limit and args.flowstart:
        print("WARNING: plot 3: flowstart will override rel_err_limit!")

    if not args.hide_err_above:
        args.hide_err_above = args.rel_err_limit * 1.001

    if 1 not in args.only_plot_no and 2 not in args.only_plot_no and 3 not in args.only_plot_no:
        parser.error('usage: --only_plot_no 1 2 3')

    if args.conftype_2 is not None and not args.reconstruct:
        parser.error('The --conftype_2 argument requires the --reconstruct argument!')

    return args


def load_data(args, conftype: str, inputfolder: str, prefix_load: str, nt: int):
    """load the reduced data (means and error)"""
    print("READING from", args.corr + prefix_load + "_" + conftype + ".dat")
    flow_radius = numpy.loadtxt(inputfolder + "flowradii_" + conftype + ".dat")
    flow_times = numpy.loadtxt(inputfolder + "flowtimes_" + conftype + ".dat")
    XX = numpy.loadtxt(inputfolder + args.corr + prefix_load + "_" + conftype + ".dat")
    XX_err = numpy.loadtxt(inputfolder + args.corr + prefix_load + "_err_" + conftype + ".dat")

    tauT = lpd.get_tauTs(nt)

    valid_flowtimes = None
    delete_indices = []
    # filter out flowtimes based on the file input
    if args.flowselectionfile:
        valid_flowtimes = numpy.loadtxt(args.flowselectionfile)
        for i, val in enumerate(flow_times):
            if val not in valid_flowtimes:
                delete_indices.append(i)
        # print("deleting flow times with the following indices:", delete_indices)
        flow_times = numpy.delete(flow_times, delete_indices, 0)
        flow_radius = numpy.delete(flow_radius, delete_indices, 0)
        XX = numpy.delete(XX, delete_indices, 0)
        XX_err = numpy.delete(XX_err, delete_indices, 0)

    return flow_radius, flow_times, XX, XX_err, tauT, valid_flowtimes


def symmetry_index(nt, input_idx):
    input_idx += 1  # the zeroeth entry in the array corresponds to tau=1, so we have to add 1 here
    if input_idx > nt / 2:
        return nt - input_idx - 1  # remove the +1 again to obtain the correct array index
    elif input_idx <= nt / 2:
        return input_idx - 1


def main():
    # parse arguments
    args = get_args()
    beta, ns, nt, nt_half = lpd.parse_conftype(args.conftype)
    fermions, temp, flowtype = lpd.parse_qcdtype(args.qcdtype)
    if fermions not in ("quenched", "hisq"):
        sys.exit("couldn't parse fermion type (either fermions=hisq or fermions=quenched) from qcdtype. use this syntax for qcdtype: <fermions>_<other_stuff>")

    prefix = "" if not args.part_obs else "_" + args.part_obs  # file path prefix for loading different observables
    prefix_load = prefix
    if args.conftype_2:
        prefix_load = "_numerator"  # why? see --help
    outputfolder = lpd.get_plot_path(args.qcdtype, args.corr, args.conftype)
    lpd.create_folder(outputfolder)
    inputfolder = lpd.get_merged_data_path(args.qcdtype, args.corr, args.conftype)

    # load data
    flow_radius, flow_times, XX, XX_err, tauT, valid_flowtimes = load_data(args, args.conftype, inputfolder, prefix_load, nt)
    XX_bak = numpy.copy(XX)
    XX_err_bak = numpy.copy(XX_err)

    # decision between sqrt(8tauF)/a and sqrt(8tauF)T for plot labels
    flow_var = numpy.copy(flow_radius) if not args.show_flowtime_instead_of_flowradii else flow_radius * nt
    flow_str = r"$ \sqrt{8\tau_\mathrm{F}}T$" if not args.show_flowtime_instead_of_flowradii else r"$ \sqrt{8\tau_\mathrm{F}}/a$"

    # decide whether to plot all (valid) flow times or to skip some
    flowend = len(flow_var) if not args.flowend else args.flowend
    flow_selection = range(args.flowstart, flowend, args.flowstep)

    # ONLY UNCOMMENT THESE WHEN USING OLD VERSIONS OF PARALLELGPUCODE WHERE MULTIPLICITIES OF DISCRETIZATIONS ARE NOT ACCOUNTED FOR
    # if corr == "EE_clover":
    # multiplicity_correction = 2
    # if corr == "BB_clover":
    # multiplicity_correction = -1.5
    # if corr == "BB":
    # multiplicity_correction = 1
    multiplicity_correction = 1

    # global plot settings
    global corr_label
    corr_label = r'X'
    if args.corr == "EE":
        corr_label = 'E'
    if args.corr == "BB":
        corr_label = 'B'

    global figsize
    figsize = (3 + 3 / 8, 3 + 3 / 8 - 1 / 2.54)
    if args.wideaspect:
        figsize = (1.5 * (3 + 3 / 8), 1.5 * (3 + 3 / 8) / 16 * 9)

    # create reconstructed correlator
    if args.reconstruct:
        nt = int(nt / 2)
        nt_half = int(nt / 2)
        tauT = tauT[1:nt:2]
        for i in range(XX.shape[0]):
            for j in range(nt_half):  # from 0 to 15 (or tau/a=[1,..,16]) for a reconstructed Nt=32 correlator based on Nt=64
                new_index = symmetry_index(int(nt * 2), (j + nt))
                # print(j, j+nt, new_index)
                XX_err[i, j] = numpy.sqrt(XX_err_bak[i, j] ** 2 + XX_err_bak[i, new_index] ** 2)
                XX[i, j] += XX_bak[i, new_index]
        # load finite temp corr.
        if args.conftype_2:
            inputfolder_2 = lpd.get_merged_data_path(args.qcdtype, args.corr, args.conftype_2)
            flow_times_2, XX_2, XX_2_err = load_data(args, args.conftype_2, inputfolder_2, prefix_load, nt)[1:4]

            for i in range(len(flow_times)):
                if flow_times[i] != flow_times_2[i]:
                    print("bla")

            # load polyakovloop data as well
            ploop1, ploop_err1 = load_data(args, args.conftype, inputfolder, "_polyakovloop", int(nt * 2))[2:4]  # nt*2 because this is the original one
            ploop2, ploop_err2 = load_data(args, args.conftype_2, inputfolder_2, "_polyakovloop", nt)[2:4]

            # TODO reconstruct polakovloop ????
            # for i in range(len(flow_var)):
            #     for j in range(nt_half):  # from 0 to 15 (or tau/a=[1,..,16]) for a reconstructed Nt=32 correlator based on Nt=64
            #         new_index = symmetry_index(int(nt * 2), (j + nt))
            #         # print(j, j+nt, new_index)
            #         ploop_err1[i, j] = numpy.sqrt(ploop_err1[i, j] ** 2 + ploop_err1[i, new_index] ** 2)
            #         ploop1[i, j] += ploop1[i, new_index]

            for i in range(XX.shape[0]):
                for j in range(nt_half):
                    ploop1[i, j] = ploop2[i, j] / ploop1[i, j]
                    # f = a A/B , df = a A/B sqrt( (dA/A)^2 (dB/B)^2 )
                    XX_err[i, j] = numpy.fabs(XX_2[i, j] / XX[i, j]) * numpy.sqrt((XX_err[i, j] / XX[i, j]) ** 2 + (XX_2_err[i, j] / XX_2[i, j]) ** 2) / ploop1[
                        i, j]
                    XX[i, j] = XX_2[i, j] / XX[i, j] / ploop1[i, j]

    skip_large_errors(XX, XX_err, args.hide_err_above)

    # decision between tauT and tau/a for plot labels
    xdata_plot = numpy.copy(tauT)
    if args.show_TauByA:
        xdata_plot = tauT * nt

    # DEBUG OPTION: reduce errors for plotting. default is to not do this. (args.error_reduce_factor=1)
    for i in range(len(flow_var)):
        for j in range(nt_half):
            XX_err[i, j] /= args.error_reduce_factor

    # append "_rec" to output folder name if reconstruct
    if args.reconstruct:
        args.conftype += "_rec"
        if args.conftype_2:
            args.conftype += "_ratio"

    # ==================================================================================================================
    # ===================== PLOT: corr normalized to temperature, x-axis tauT, flowtimes fixed =========================
    # ==================================================================================================================

    # normalize to temperature (as long as were not computing this special ratio with two conftypes)
    if not args.conftype_2 and not args.part_obs == "polyakovloop":
        for i in range(len(flow_var)):
            for j in range(nt_half):
                XX[i, j] *= nt ** 4 * multiplicity_correction
                XX_err[i, j] *= nt ** 4 * multiplicity_correction

    # plot the first figure
    if 1 in args.only_plot_no and not args.part_obs == "polyakovloop":
        plot1(XX, XX_err, args, prefix, flow_selection, flow_var, xdata_plot, nt_half, flow_str, outputfolder)

    # ==================================================================================================================
    # ==================== PLOT: corr normalized to LO pert corr, x-axis tauT, flowtimes fixed =========================
    # ==================================================================================================================

    # normalize data to perturbative lattice correlators
    if not args.part_obs and not args.conftype_2:
        for i in range(len(flow_var)):
            for j in range(nt_half):
                XX[i, j] *= lpd.improve_corr_factor(j, nt, i) / nt ** 4
                XX_err[i, j] *= lpd.improve_corr_factor(j, nt, i) / nt ** 4

    # plot the second figure
    if 2 in args.only_plot_no and not args.part_obs == "polyakovloop":
        plot2(XX, XX_err, args, prefix, flow_selection, flow_var, xdata_plot, nt_half, flow_str, outputfolder,
              fermions)

    # ==================================================================================================================
    # ========================= PLOT: x-axis flowtimes, tauT fixed =====================================================
    # ==================================================================================================================

    # plot the third figure
    if 3 in args.only_plot_no:
        plot3(XX, XX_err, args, prefix, flow_var, xdata_plot, nt_half, outputfolder, fermions, valid_flowtimes, tauT, flow_radius, flowend)

    # ==================================================================================================================

    lpd.save_script_call(outputfolder)


if __name__ == '__main__':
    main()