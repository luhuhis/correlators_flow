#!/usr/local/bin/python3.7m
import numpy
import lib_process_data as lpd


def load_merged_data(qcdtype, corr, conftype):
    inputfolder = lpd.get_merged_data_path(qcdtype, corr, conftype)

    print("read  " + inputfolder + "flowtimes_" + conftype + ".dat")
    flow_times = numpy.loadtxt(inputfolder + "flowtimes_" + conftype + ".dat")
    n_flow = len(flow_times)

    print("read  " + inputfolder + "n_datafiles_" + conftype + ".dat")
    metadata = numpy.loadtxt(inputfolder + "n_datafiles_" + conftype + ".dat")
    n_datafiles, n_streams = [int(i) for i in metadata[0:2]]
    n_discarded = [int(i) for i in metadata[2:2+n_streams]]
    n_conf_per_stream = [int(i) for i in metadata[2+n_streams:]]

    print("read  " + inputfolder + "polyakov_real_" + conftype + "_merged.dat")
    polyakov_real = numpy.loadtxt(inputfolder + "polyakov_real_" + conftype + "_merged.dat")

    return flow_times, n_flow, n_datafiles, n_streams, n_conf_per_stream, n_discarded, polyakov_real


def main():

    parser, requiredNamed = lpd.get_parser()
    requiredNamed.add_argument('--conftype', help='format example: s096t32_b0824900_m002022_m01011', type=str, required=True)
    parser.add_argument('--outputpath', help='where to store the plot', type=str)
    parser.add_argument('--streamids', help='how to call the streams in the plot', nargs='*', type=str)
    parser.add_argument('--show_id', help="list of indices that shall be plotted. e.g. 0 1 2 3 for the first 4 streams", type=int, nargs='*')
    parser.add_argument('--vlines', help="where to plot vlines", nargs='*', type=float)
    parser.add_argument('--offsets', help="trajctory offset for each stream.", nargs='*', type=int)
    parser.add_argument('--suffix', help="suffix for output file", default="")
    args = parser.parse_args()
    beta, ns, nt, nt_half = lpd.parse_conftype(args.conftype)

    outputfolder = lpd.get_plot_path(args.qcdtype, args.corr, args.conftype) if not args.outputpath else args.outputpath
    lpd.create_folder(outputfolder)

    flow_times, n_flow, n_datafiles, n_streams, n_conf_per_stream, n_discarded, polyakov_real = load_merged_data(args.qcdtype, args.corr, args.conftype)

    xlabel = r'n traj.'
    xlabel2 = r'n traj. discarded from left'
    xlabel3 = r'n traj. discarded from left'
    ylabel = r'$\mathrm{Re}($P.loop$)$,$\tau_F/a^2='+r'{0:.3f}'.format(flow_times[-1])+r', \sqrt{8\tau_F}T \approx '+r'{0:.2f}'.format(numpy.sqrt(8*flow_times[-1])/nt)+r'$'
    ylabel2 = r'Re(P.loop) stddev without first n traj.,\\ $\sqrt{8\tau_F}T \approx'+r'{0:.2f}'.format(numpy.sqrt(8*flow_times[-1])/nt)+r'$'
    ylabel3 = r'Re(P.loop) mean, $\sqrt{8\tau_F}T \approx'+r'{0:.2f}'.format(numpy.sqrt(8*flow_times[-1])/nt)+r'$'
    ylabelpos = (-0.1, 0.5)
    stepsize = 5  # number of trajectories between written out confs

    streamids = ("a", "b", "c", "d", "e", "f", "g", "h", "i", "j") if not args.streamids else args.streamids
    flowidx = n_flow-1
    usetex = False
    fig, ax, _ = lpd.create_figure(xlabel=xlabel, ylabel=ylabel, xlabelpos=(0.5, -0.1), ylabelpos=ylabelpos, UseTex=usetex, figsize=(16/9*(3+3/8), 3+3/8))  # ylims=(0, 5), xlims=(0,5),
    fig2, ax2, _ = lpd.create_figure(xlabel=xlabel2, ylabel=ylabel2, xlabelpos=(0.5, -0.1), ylabelpos=ylabelpos, UseTex=usetex, figsize=(16/9*(3+3/8), 3+3/8))
    fig3, ax3, plots = lpd.create_figure(xlabel=xlabel3, ylabel=ylabel3, xlabelpos=(0.5, -0.1), ylabelpos=ylabelpos, UseTex=usetex, figsize=(16/9*(3+3/8), 3+3/8))

    for thisax in ax, ax2, ax3:
        thisax.set_ylabel(thisax.get_ylabel(), fontsize=8, rotation=90)

    for i, nconf in enumerate(n_conf_per_stream):
        if i in args.show_id:
            print(nconf, streamids[i])
            offset = numpy.sum(n_conf_per_stream[:i]) if i > 0 else 0
            x = range(args.offsets[i]+n_discarded[i]*stepsize, n_discarded[i]*stepsize+args.offsets[i]+nconf*stepsize, stepsize)
            y = [polyakov_real[int(flowidx+n_flow*(k+offset))] for k in range(nconf)]
            ax.errorbar(x, y, fmt='.-', lw=0.3, fillstyle='full', markersize=2.5, mew=0, label=streamids[i], alpha=1)
            errs = []
            means = []
            for j in range(0,nconf):
                try:
                    errs.append(numpy.std(y[j:]))
                    means.append(numpy.mean(y[j:]))
                except:
                    continue
            errs = numpy.asarray(errs)
            means = numpy.asarray(means)
            ax2.errorbar(x, errs, fmt='.', fillstyle='full', markersize=1, mew=0.25, lw=0.5, elinewidth=0.25, capsize=1.2, alpha=0.5, label=streamids[i])
            ax3.errorbar(x, means, fmt='.', fillstyle='full', markersize=1, mew=0.25, lw=0.5, elinewidth=0.25, capsize=1.2, alpha=0.5, label=streamids[i])

    legendoptions = dict(edgecolor='none', fancybox=False, facecolor="w", columnspacing=0.1,
                   labelspacing=0.1, borderpad=0.1, handletextpad=0.4, handlelength=1)
    ax.legend(**legendoptions)
    ax2.legend(**legendoptions)
    ax3.legend(**legendoptions)
    if args.vlines:
        for x in args.vlines:
         for thisax in ax, ax2, ax3:
            thisax.axvline(x=x, **lpd.verticallinestyle)
    ax.set_title(r'$'+str(ns)+r'^3 \times '+str(nt)+r'$, ms/5, HISQ', x=0.5)
    ax2.set_title(r'$'+str(ns)+r'^3 \times '+str(nt)+r'$')
    ax3.set_title(r'$' + str(ns) + r'^3 \times ' + str(nt) + r'$')

    for thisax in ax, ax2, ax3:
        thisax.set_xlim(xmin = 0)

    fig.savefig(outputfolder+"/"+"polyakovloop_"+args.suffix+".pdf")
    fig2.savefig(outputfolder+"/"+"polyakovloop_err_"+args.suffix+".pdf")
    fig3.savefig(outputfolder+"/"+"polyakovloop_mean_"+args.suffix+".pdf")



if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
