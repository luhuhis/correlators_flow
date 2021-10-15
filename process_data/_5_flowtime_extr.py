#!/usr/local/bin/python
import numpy
import matplotlib
from matplotlib import pyplot, container, legend_handler
import lib_process_data as lpd
from latqcdtools import fitting
from latqcdtools import bootstr


def extrapolation_ansatz(x, m, b):
    return m * x + b 


def fit_sample(ydata, xdata, edata):
    fitter = fitting.Fitter(func=extrapolation_ansatz, xdata=xdata, ydata=ydata, edata=edata, func_sup_numpy=True, always_return=False)
    fitparams, fitparams_err, chi_dof = fitter.try_fit(algorithms=["curve_fit"], start_params=[-100, 3])
    return fitparams  # , fitparams_err


def main():
    
    # parse cmd line arguments
    parser, requiredNamed = lpd.get_parser()

    requiredNamed.add_argument('--flowradii_file', help='file with list of all flow radii', required=True)
    requiredNamed.add_argument('--finest_Nt', help='finest Nt used in previous cont extr.', required=True, type=int)
    requiredNamed.add_argument('--coarsest_Nt', help='coarsest Nt used in previous cont extr. needed for lower flow time limit.', required=True, type=int)

    parser.add_argument('--nsamples', help="number of artifical gaussian bootstrap samples to generate", type=int, default=200)
    parser.add_argument('--use_tex', type=bool, default=True)
    parser.add_argument('--error_treshhold', help='factor of tauT that limits the relative error', default='0.01', type=float)
    parser.add_argument('--flowend', help='maximum flow index for all tauT that should be used in the extrapolation', type=int, default=134)
    parser.add_argument('--no_extr', help='do NOT perform an extrapolation, just plot the continuum data', action="store_true")
    parser.add_argument('--custom_ylims', help="custom y-axis limits for both plots", type=float, nargs=2)

    args = parser.parse_args()

    basepath = lpd.get_merged_data_path(args.qcdtype, args.corr, "")
    plotbasepath = lpd.get_plot_path(args.qcdtype, args.corr, "")

    # parse some more arguments
    flowradii = numpy.loadtxt(args.flowradii_file)
    nflow = len(flowradii)
    finest_tauTs = lpd.get_tauTs(args.finest_Nt)
    ntauT = len(finest_tauTs)
    finest_Nt_half = int(args.finest_Nt/2)

    # variable to hold the continuum data
    XX = numpy.empty((ntauT, nflow, 2))
    XX[:] = numpy.nan

    # load continuum extr for each flowtime, only extract the desired tauT
    for i, flowradius in enumerate(flowradii):
        flowradius_str = '{0:.4f}'.format(flowradius)
        try:
            tmp = numpy.loadtxt(basepath+"/cont_extr/"+args.corr+"_"+flowradius_str+"_cont.txt")
            for j, row in enumerate(tmp):
                XX[j][i][0] = row[1]
                XX[j][i][1] = row[2]
        except OSError:
            pass

    # filter low distance high error data that is not really linear
    flow_extr_filter = []  # this is the firts flowindex that shall be used
    for i in range(finest_Nt_half):  # loop over tauT
        rel_err = numpy.fabs(XX[i][:, 1]/XX[i][:, 0])
        found = False
        for j, val in enumerate(rel_err):  # loop over flowtimes
            if val <= args.error_treshhold*finest_tauTs[i]:
                flow_extr_filter.append(j)
                found = True
                break
        if not found:
            flow_extr_filter.append(-1)
    flowstarts = flow_extr_filter  

    # declarations
    results = numpy.zeros((ntauT, 3))
    mintauTindex = None

    # plot settings
    displaystyle = '' if not args.use_tex else r'\displaystyle'
    ylims = (2.5, 3.8) if not args.custom_ylims else args.custom_ylims
    lpd.labelboxstyle.update(dict(alpha=0.8))

    fig, ax, plots = lpd.create_figure(xlims=[-0.0001, 0.0033], ylims=ylims, xlabel=r'$\tau_\mathrm{F} T^2$',
                                       ylabel=r'$' + displaystyle + r'\frac{G^\mathrm{cont }}{G_{\tau_\mathrm{F}=0}^{\substack{\text{\tiny norm} \\[-0.5ex] \text{\tiny cont }}}}$',
                                       xlabelpos=(0.94, 0.05), ylabelpos=(0.08, 0.96), UseTex=args.use_tex)

    with open(basepath+"/"+args.corr+"_flow_extr_quality.txt", 'w') as outfile:
        outfile.write('# flowtime zero extrapolation quality at fixed tauT for all tauT of Ntau=36 \n')
        outfile.write('# tauT    tau_F    G/G_norm     err \n')

        with open(basepath+"/"+args.corr+"_flow_extr_quality_slope.txt", 'w') as outfile_slope:
            outfile_slope.write('# flowtime zero extrapolation slope for all tauT of Ntau=36 \n')
            outfile_slope.write('# tauT     slope     slope_err \n')

            for i, tauT in enumerate(finest_tauTs):
                flowstart = flowstarts[i]

                if not numpy.isnan(flowstart):

                    # organize data
                    maxflowtime = (tauT-1/args.coarsest_Nt)*numpy.sqrt(8*0.014)
                    flowend = min((numpy.abs(flowradii - maxflowtime)).argmin(), args.flowend)
                    xdatatmp = flowradii[flowstart:flowend]**2/8
                    ydatatmp = XX[i][flowstart:flowend, 0]
                    edatatmp = XX[i][flowstart:flowend, 1]
                    mask = numpy.isnan(ydatatmp)
                    xdata = xdatatmp[~mask]
                    ydata = ydatatmp[~mask]
                    edata = edatatmp[~mask]

                    if len(xdata) > 2:  # minimum amount for linear fit

                        if mintauTindex is None:
                            mintauTindex = i

                        mycolor = lpd.get_color(finest_tauTs, finest_Nt_half-1-i+mintauTindex, mintauTindex, finest_Nt_half-1)

                        lpd.plotstyle_add_point_single.update(dict(fmt='-', markersize=2, mew=0.25))
                        # ax.errorbar([x for i,x in enumerate(xdata) if i%2==0], [y for i,y in enumerate(ydata) if i%2==0], **lpd.plotstyle_add_point_single, color='grey', zorder=-100+i)#, label='{0:.3f}'.format(tauT))
                        # ax.fill_between(xdata, ydata-edata, ydata+edata, facecolor=mycolor, alpha=1, zorder=-300+i)#, label='{0:.3f}'.format(tauT))
                        if args.no_extr:
                            plots.append(ax.errorbar(xdata, ydata, edata, **lpd.plotstyle_add_point_single, color=mycolor, zorder=-100+i, label='{0:.3f}'.format(tauT)))
                        else:
                            ax.errorbar(xdata, ydata, edata, **lpd.plotstyle_add_point_single, color=mycolor, zorder=-100 + i)

                        if not args.no_extr:
                            # perform extrapolation
                            fitparams, fitparams_err = bootstr.bootstr_from_gauss(fit_sample, data=ydata, data_std_dev=edata, numb_samples=args.nsamples, sample_size=1, return_sample=False, args=[xdata, edata])
                            results[i][0] = tauT
                            results[i][1] = fitparams[1]
                            results[i][2] = fitparams_err[1]
                            # noinspection PyTypeChecker
                            numpy.savetxt(outfile_slope, numpy.stack(([tauT], [fitparams[0]], [fitparams_err[0]]), axis=-1))

                            lpd.plotstyle_add_point_single.update(dict(fmt=lpd.markers[i - finest_Nt_half], mew=0.25, markersize=5))
                            plots.append(ax.errorbar(0, fitparams[1], fitparams_err[1], **lpd.plotstyle_add_point_single, alpha=1, color=mycolor, zorder=1, label='{0:.3f}'.format(tauT)))
                            x = numpy.linspace(0, 0.1, 100)
                            ax.errorbar(x, extrapolation_ansatz(x, *fitparams[0:2]), color=mycolor, alpha=1, fmt=':', lw=0.5, zorder=-100)

                            # noinspection PyTypeChecker
                            numpy.savetxt(outfile, numpy.stack(([tauT for _ in range(len(xdata) + 1)], [0, *xdata], [fitparams[1], *ydata], [fitparams_err[1], *edata]), axis=-1))
                            outfile.write('# \n')
                        else:
                            results[i] = None
                    else:
                        results[i] = None
                else:
                    results[i] = None
                print("done ", '{0:.3f}'.format(tauT))

    # noinspection PyTypeChecker
    numpy.savetxt(basepath+"/"+args.corr+"_final.txt", results, header="tauT    G/Gnorm    err")
    
    # second x-axis for flow radius
    ax2 = ax.twiny()
    new_tick_locations = numpy.array([0, 0.05**2/8, 0.08**2/8, 0.1**2/8, 0.13**2/8, 0.15**2/8])

    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(new_tick_locations)
    ax2.set_xticklabels(["%.1f" % z if z == 0 else "%.2f" % z for z in numpy.sqrt(new_tick_locations*8)])
    ax2.set_xlabel(r'$\sqrt{8\tau_\mathrm{F}}T$', horizontalalignment='right', verticalalignment='top', bbox=lpd.labelboxstyle, zorder=999999)
    ax2.xaxis.set_label_coords(0.99, 0.97)
    ax2.tick_params(direction='in', pad=0, width=0.5)

    ax.axvline(x=0, ymin=((results[mintauTindex, 1]-ylims[0])/(ylims[1]-ylims[0])), ymax=((results[finest_Nt_half-1, 1]-ylims[0])/(ylims[1]-ylims[0])), alpha=1, color='grey', zorder=-1000, lw=0.5, dashes=(5, 2))
    
    # lpd.legendstyle.update(dict(handlelength=0.5))
    lpd.legendstyle.update(dict(loc="lower right", bbox_to_anchor=(1.01, 0.01), columnspacing=0.5, handlelength=0.75, labelspacing=0.1, handletextpad=0.2,
                                borderpad=0, handler_map={matplotlib.container.ErrorbarContainer: matplotlib.legend_handler.HandlerErrorbar(yerr_size=0.4)}))
    ax.legend(handles=plots)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1], title=r'$\tau T=$', ncol=2, **lpd.legendstyle).set_zorder(-1)
    matplotlib.pyplot.tight_layout(0)
    fig.savefig(plotbasepath+"/"+args.corr+"_flow_extr_quality.pdf")
    ax.lines.clear()
    ax.collections.clear()
    plots.clear()

    # plot final double-extrapolated correlator in its one plot
    fig, ax, plots = lpd.create_figure(xlims=[0.15, 0.51], ylims=[1.5, 4], xlabel=r'$\tau T$',
                                       ylabel=r'$'+displaystyle+r'\frac{G^\mathrm{cont}_{\tau_\mathrm{F}\rightarrow 0}}{ G^{\substack{ \text{\tiny  norm} \\[-0.4ex] \text{\tiny cont } } }_{\tau_\mathrm{F} = 0} }$', UseTex=args.use_tex)
    lpd.plotstyle_add_point.update(dict(fmt='D-'))
    results = numpy.swapaxes(results, 0, 1)
    ax.errorbar(results[0], results[1], results[2], color='black', **lpd.plotstyle_add_point)
    matplotlib.pyplot.tight_layout(0.2)
    fig.savefig(plotbasepath+"/"+args.corr+"_final.pdf")


if __name__ == '__main__':
    main()