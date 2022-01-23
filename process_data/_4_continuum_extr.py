#!/usr/local/bin/python3.7m -u

import lib_process_data as lpd
import numpy
import matplotlib
from matplotlib import container, legend_handler
import re
import scipy.optimize
from latqcdtools import bootstr
import sys

# input
#   XX_mean_finest = numpy.loadtxt(lpd.get_merged_data_path(qcdtype, conftypes[-1], corr)+"/"+corr+"_"+conftypes[-1]+".dat")
#   lpd.tauT(nt_fine)
#   lpd.get_merged_data_path+"/"+conftype+"/interpolations/EE_"+flowradius_str+"_interpolation.txt", for each conftype and flowradius

# output
#   lpd.inputfolder+"/cont_extr_quality/EE_"+flowradius_str+"_cont_quality.txt"
#   lpd.outputfolder+"/cont_extr_quality/EE_"+flowradius_str+"_cont_quality.pdf"
#   lpd.inputfolder+"/cont_extr/EE_"+flowradius_str+"_cont.txt"
#   lpd.outputfolder+"/cont_extr/EE_"+flowradius_str+"_cont.pdf"


def extrapolation_ansatz(x, m, b):
    return m * x + b


def fit_sample(ydata, xdata, edata, start_params=None):
    if start_params is None:
        start_params = [80, 3]
    fitparams, _ = scipy.optimize.curve_fit(extrapolation_ansatz, xdata, ydata, p0=start_params, sigma=edata)
    return fitparams


def main():

    # parse cmd line arguments
    parser, requiredNamed = lpd.get_parser()

    requiredNamed.add_argument('--flow_index', help='which flow time to interpolate', type=int, required=True)
    requiredNamed.add_argument('--conftypes', nargs='*', required=True,
                               help="ORDERED list of conftypes (from coarse to fine), e.g. s080t20_b0703500 s096t24_b0719200 s120t30_b0739400 s144t36_b0754400")

    parser.add_argument('--use_imp', help='whether to use tree-level improvement for the finest lattice. use the same setting you used for the interpolation of the coarser lattices!', type=bool, default=True)
    parser.add_argument('--nsamples', help="number of artifical gaussian bootstrap samples to generate", type=int, default=400)
    parser.add_argument('--use_tex', action="store_true", help="use LaTeX when plotting")
    parser.add_argument('--start_at_zero', action="store_true", help="replace the flow indices from which on the continuum extr shall be performed")
    parser.add_argument('--custom_ylims', help="custom y-axis limits for both plots", type=float, nargs=2)
    parser.add_argument('--int_only', help="load an interpolation of the finest lattice instead of the original non-interpolated data. useful if you want to "
                                           "look at more tauT-points than what the finest lattice can offer.", action="store_true")
    parser.add_argument('--int_Nt', help="if --int_only, then specify the Nt, for whose tauT values the interpolations were saved at", type=int)
    parser.add_argument('--grace_factor', help='make the lower_tauT_limit less strict by using grace_factor < 1 in order to allow higher flow times. '
                                               'use the same value here that was also used for the interpolation.', type=float, default=1)

    args = parser.parse_args()

    if (args.int_only and not args.int_Nt) or (args.int_Nt and not args.int_only):
        parser.error("Error while parsing arguments. --int_only requires --int_Nt and vice versa.")

    # load flow radius
    flowradius = numpy.loadtxt(lpd.get_merged_data_path(args.qcdtype, args.corr, args.conftypes[-1])+"/flowradii_"+args.conftypes[-1]+".dat")[args.flow_index]
    flowradius_str = r'{0:.4f}'.format(flowradius)

    # set some params from args
    Nts = []
    for conftype in args.conftypes:
        dummy = re.sub(r'(^.*?)t', '', conftype)
        Nts.append(int(re.sub(r'(^.*?)_b(.*)', r'\1', dummy)))
    nt_coarse = Nts[0]
    if not args.int_Nt:
        nt_fine = Nts[-1]
    else:
        nt_fine = args.int_Nt
    print("INFO: perform extrapolation using Nts:", Nts)

    if flowradius < 1/nt_coarse:
        exit("ERROR: need at least a flow radius sqrt(8t)/a of 1/"+str(nt_coarse)+" to have enough flow for cont extr")

    lower_tauT_lim = args.grace_factor*lpd.lower_tauT_limit(flowradius, nt_coarse)
    if not args.int_only:
        # load finest lattice
        XX_mean_finest = numpy.loadtxt(lpd.get_merged_data_path(args.qcdtype, args.corr, args.conftypes[-1])+"/"+args.corr+"_"+args.conftypes[-1]+".dat")
        XX_err_finest = numpy.loadtxt(lpd.get_merged_data_path(args.qcdtype, args.corr, args.conftypes[-1])+"/"+args.corr+"_err_"+args.conftypes[-1]+".dat")

        # --- only use points above the lower limit
        tauTs_fine = lpd.get_tauTs(nt_fine)
        print("INFO: lower tauT limit for this flow time and coarsest lattice: ", lower_tauT_lim)
        XX_finest = numpy.asarray(([tauT for tauT in tauTs_fine if tauT > lower_tauT_lim],
                                   [val*lpd.improve_corr_factor(j, nt_fine, args.flow_index, args.use_imp) for j, val in enumerate(XX_mean_finest[args.flow_index])
                                    if tauTs_fine[j] > lower_tauT_lim],
                                   [val*lpd.improve_corr_factor(j, nt_fine, args.flow_index, args.use_imp) for j, val in enumerate(XX_err_finest[args.flow_index])
                                    if tauTs_fine[j] > lower_tauT_lim]))
    else:
        tauTs_fine = lpd.get_tauTs(args.int_Nt)

    # --- load interpolations
    XX_ints = []
    for conftype in args.conftypes:
        if conftype == args.conftypes[-1] and not args.int_only:
            XX_ints.append(XX_finest)
        else:
            try:
                tmp = numpy.loadtxt(lpd.get_merged_data_path(args.qcdtype, args.corr, conftype)+"/interpolations/"+args.corr+"_"+flowradius_str+"_interpolation.txt", unpack=True)
            except OSError:
                sys.exit("ERROR while reading "+lpd.get_merged_data_path(args.qcdtype, args.corr, conftype)+"/interpolations/"+args.corr+"_"+flowradius_str+"_interpolation.txt")
            tmp2 = numpy.asarray(([tauT for tauT in tmp[0] if tauT in tauTs_fine],
                                  [val for j, val in enumerate(tmp[1]) if tmp[0][j] in tauTs_fine],
                                  [val for j, val in enumerate(tmp[2]) if tmp[0][j] in tauTs_fine]))
            XX_ints.append(tmp2)

    # define some parameters
    if not args.int_only:
        valid_tauTs = XX_finest[0]
    else:
        valid_tauTs = [tauT for tauT in XX_ints[-1][0] if tauT > lower_tauT_lim]
        print(valid_tauTs)
    n_valid_tauTs = len(valid_tauTs)
    nt_half_fine = int(nt_fine/2)
    offset = nt_half_fine-n_valid_tauTs

    # skip flow times (indices) that are smaller than these, as it turned out in retrospect that they are not used in the flow time extrapolation
    # flowstarts contain for each tau the minimum index of flow times. if args.flow_index is less than
    # if args.start_at_zero:
    #     flowstarts = [0 for dummy in range(0, nt_half_fine)]
    # else:
    #     flowstarts = [-1, -1, -1, -1, -1, -1, -1, 50, 53, 57, 61, 65, 68, 71, 74, 77, 82, 83]  # "hardcoded" as used in PhysRevD.103.014511
    #
    # maxtauTindex = 0
    # for j in range(nt_half_fine):
    #     if args.flow_index >= flowstarts[j]:
    #         maxtauTindex = j-offset

    if args.use_tex:
        displaystyle = r'\displaystyle'
        ylabel = r'$'+displaystyle+r'\frac{G^\mathrm{latt }} { G^{\substack{ \text{\tiny  norm} \\[-0.4ex] \text{\tiny latt } } }_{\tau_\mathrm{F} = 0} } $'
    else:
        displaystyle = ''
        ylabel = 'G'

    # plot settings
    ylims = (2.55, 3.75) if not args.custom_ylims else args.custom_ylims
    fig, ax, plots = lpd.create_figure(xlims=[-0.0001, 1/nt_coarse**2*1.05], ylims=ylims, xlabel=r'$N_\tau^{-2}$', ylabel=ylabel,
                                       xlabelpos=(0.95, 0.07), ylabelpos=(0.08, 0.98), UseTex=args.use_tex)
    lpd.titlestyle.update(dict(y=0.95))
    ax.set_title(r'$ \sqrt{8\tau_\mathrm{F}}T = '+flowradius_str+r'$', **lpd.titlestyle)

    plotpath = lpd.get_plot_path(args.qcdtype, args.corr, "")
    ceq_prefix = lpd.get_merged_data_path(args.qcdtype, args.corr, "")+"/cont_extr_quality/"+args.corr+"_"+flowradius_str
    lpd.create_folder(plotpath+"/cont_extr_quality/",
                      plotpath+"/cont_extr/",
                      lpd.get_merged_data_path(args.qcdtype, args.corr, "")+"/cont_extr_quality/",
                      lpd.get_merged_data_path(args.qcdtype, args.corr, "")+"/cont_extr/")

    # declarations
    results = numpy.empty((nt_half_fine, 3))
    results[:] = numpy.nan
    maxtauTindex_plot = 0
    mintauTindex = None

    with open(ceq_prefix+"_cont_quality.txt", 'w') as outfile:
        outfile.write('# continuum extrapolation quality at fixed flowtime for all tauT of Ntau='+str(nt_fine)+r' \n')
        outfile.write('# tauT    1/Ntau^2    G/G_norm     err \n')
        with open(ceq_prefix+"_cont_quality_slope.txt", 'w') as outfile_slope:
            outfile_slope.write('# continuum extrapolation slope at fixed flowtime for all tauT of Ntau='+str(nt_fine)+r' \n')
            outfile_slope.write('# tauT     slope     slope_err \n')

            # perform cont extr for each tauT at fixed flowtime
            for j, tauT in enumerate(tauTs_fine):
                xdata = numpy.asarray([1/Ntau**2 for k, Ntau in enumerate(Nts) if XX_ints[k] is not None])

                # skip if flow time is too small
                if tauT in valid_tauTs and len(xdata) >= 3:  # and args.flow_index >= flowstarts[j]:
                    print("working on flowtime ", flowradius_str, " tauT=", '{0:.4f}'.format(tauT))

                    # actually perform extr
                    # for m, XX_int in enumerate(XX_ints):
                    # print("XX_int no ", m, "of length", len(XX_int[1]))

                    ydata = [XX_int[1][j-offset] for XX_int in XX_ints if XX_int is not None]
                    edata = [XX_int[2][j-offset] for XX_int in XX_ints if XX_int is not None]
                    fitparams, fitparams_err = bootstr.bootstr_from_gauss(fit_sample, data=ydata, data_std_dev=edata, numb_samples=args.nsamples,
                                                                          sample_size=1, return_sample=False, args=[xdata, edata])
                    # fitparams, fitparams_err = fit_sample(ydata, xdata, edata)

                    results[j][0] = tauT
                    results[j][1] = fitparams[1]
                    results[j][2] = fitparams_err[1]
                    # noinspection PyTypeChecker
                    numpy.savetxt(outfile_slope, numpy.stack(([tauT], [fitparams[0]], [fitparams_err[0]]), axis=-1))

                    # plot extrapolation
                    if mintauTindex is None:
                        mintauTindex = j
                    if j > maxtauTindex_plot:
                        maxtauTindex_plot = j
                    mycolor = lpd.get_color(tauTs_fine, nt_half_fine-1-j+offset, mintauTindex, nt_half_fine-1)
                    lpd.plotstyle_add_point_single.update(dict(fmt=lpd.markers[j-nt_half_fine+len(lpd.markers)]))
                    plots.append(ax.errorbar(xdata, ydata, edata, **lpd.plotstyle_add_point_single, color=mycolor, zorder=-100+j, label='{0:.3f}'.format(tauT)))
                    ax.errorbar(0, fitparams[1], fitparams_err[1], **lpd.plotstyle_add_point_single, color=mycolor, zorder=1)
                    x = numpy.linspace(0, 0.1, 100)
                    ax.errorbar(x, extrapolation_ansatz(x, *fitparams), color=mycolor, alpha=1, fmt=':', lw=0.5, zorder=-100)
                else:
                    results[j][0] = tauT
                    results[j][1] = numpy.nan
                    results[j][2] = numpy.nan
                    # noinspection PyTypeChecker
                    numpy.savetxt(outfile_slope, numpy.stack(([tauT], [numpy.nan], [numpy.nan]), axis=-1))
                    continue
                # save extrapolation to file
                # noinspection PyTypeChecker
                numpy.savetxt(outfile, numpy.stack(([tauT for _ in range(len(xdata)+1)], [0, *xdata], [results[j][1], *ydata], [results[j][2], *edata]), axis=-1))
                outfile.write('# \n')
    # lpd.verticallinestyle.update(dict())#, dashes=(4,1)) )
    # noinspection PyTypeChecker
    numpy.savetxt(lpd.get_merged_data_path(args.qcdtype, args.corr, "")+"/cont_extr/"+args.corr+"_"+flowradius_str+"_cont.txt", results,
                  header="tauT    G/G_norm    err")

    # save continuum extrapolation quality plot for this flow time
    lpd.plotstyle_add_point_single.update(dict(fmt='D-'))
    ax.axvline(x=0, ymin=((results[mintauTindex, 1]-ylims[0])/(ylims[1]-ylims[0])), ymax=((results[maxtauTindex_plot, 1]-ylims[0])/(ylims[1]-ylims[0])), alpha=1, color='grey',
               zorder=-1000, lw=0.5, dashes=(5, 2))
    lpd.legendstyle.update(dict(loc="lower left", bbox_to_anchor=(-0.01, -0.01), columnspacing=0.1, labelspacing=0.25, handletextpad=0, borderpad=0,
                                framealpha=0, handler_map={matplotlib.container.ErrorbarContainer: matplotlib.legend_handler.HandlerErrorbar(xerr_size=0.4)}))
    ax.legend(handles=plots)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1], title=r'$\tau T=$', **lpd.legendstyle, ncol=3)  # reverse ordering of legend
    # matplotlib.pyplot.tight_layout(0)
    fig.savefig(plotpath+"/cont_extr_quality/"+args.corr+"_"+flowradius_str+"_cont_quality.pdf")
    ax.lines.clear()
    ax.collections.clear()
    plots.clear()

    # save plot of single continuum extrapolation for this flow time
    if args.use_tex:
        displaystyle = r'\displaystyle'
        ylabel = r'$'+displaystyle+r'\frac{G^\mathrm{cont}}{ G^{\substack{ \text{\tiny  norm} \\[-0.4ex] \text{\tiny latt } } }_{\tau_\mathrm{F} = 0} }$'
    else:
        ylabel = 'G'
    fig, ax, plots = lpd.create_figure(xlims=[0.15, 0.51], ylims=[1.4, 4], xlabel=r'$\tau T$',
                                       ylabel=ylabel, UseTex=args.use_tex)
    ax.set_title(r'$ \sqrt{8\tau_\mathrm{F}}T = $'+flowradius_str, x=0.5, y=0.89)
    ax.axvline(x=lower_tauT_lim, **lpd.verticallinestyle)
    results = numpy.swapaxes(results, 0, 1)
    ax.errorbar(results[0], results[1], results[2], color="black", **lpd.plotstyle_add_point)
    # matplotlib.pyplot.tight_layout(0.1)
    fig.savefig(plotpath+"/cont_extr/"+args.corr+"_"+flowradius_str+"_cont.pdf")


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
