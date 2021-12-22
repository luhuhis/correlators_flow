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
    n_datafiles, n_streams, n_discarded = [int(i) for i in metadata[0:3]]
    n_conf_per_stream = [int(i) for i in metadata[3:]]

    print("read  " + inputfolder + "polyakov_real_" + conftype + "_merged.dat")
    polyakov_real = numpy.loadtxt(inputfolder + "polyakov_real_" + conftype + "_merged.dat")

    return flow_times, n_flow, n_datafiles, n_streams, n_conf_per_stream, n_discarded, polyakov_real


def main():

    parser, requiredNamed = lpd.get_parser()
    requiredNamed.add_argument('--conftype', help='format example: s096t32_b0824900_m002022_m01011', type=str, required=True)
    parser.add_argument('--outputpath', help='where to store the plot', type=str)
    parser.add_argument('--streamids', help='how to call the streams in the plot', nargs='*', type=str)
    args = parser.parse_args()
    beta, ns, nt, nt_half = lpd.parse_conftype(args.conftype)

    outputfolder = lpd.get_plot_path(args.qcdtype, args.corr, args.conftype) if not args.outputpath else args.outputpath
    lpd.create_folder(outputfolder)

    flow_times, n_flow, n_datafiles, n_streams, n_conf_per_stream, n_discarded, polyakov_real = load_merged_data(args.qcdtype, args.corr, args.conftype)

    xlabel = r'no. traj.'
    ylabel = r'$\mathrm{Re}($P.loop$)$,$\tau_F/a^2='+r'{0:.3f}'.format(flow_times[-1])+r', \sqrt{8\tau_F}T \approx '+r'{0:.2f}'.format(numpy.sqrt(8*flow_times[-1])/nt)+r'$'
    ylabel2 = r'Re(P.loop) stddev without first n traj.,\\ $\sqrt{8\tau_F}T \approx'+r'{0:.2f}'.format(numpy.sqrt(8*flow_times[-1])/nt)+r'$'
    ylabel3 = r'Re(P.loop) mean without first n traj., \\$\sqrt{8\tau_F}T \approx'+r'{0:.2f}'.format(numpy.sqrt(8*flow_times[-1])/nt)+r'$'
    ylabelpos = (0.1, 0.5)
    stepsize = 5  # number of trajectories between written out confs

    streamids = ("a", "b", "c", "d", "e", "f", "g", "h", "i", "j") if not args.streamids else args.streamids
    flowidx = n_flow-1
    usetex = False
    fig, ax, _ = lpd.create_figure(xlabel=xlabel, ylabel=ylabel, xlabelpos=(0.5, -0.1), ylabelpos=ylabelpos, UseTex=usetex)  # ylims=(0, 5), xlims=(0,5),
    fig2, ax2, _ = lpd.create_figure(xlabel=xlabel, ylabel=ylabel2, xlabelpos=(0.5, -0.1), ylabelpos=ylabelpos, UseTex=usetex)
    fig3, ax3, plots = lpd.create_figure(xlabel=xlabel, ylabel=ylabel3, xlabelpos=(0.5, -0.1), ylabelpos=ylabelpos, UseTex=usetex)

    for i, nconf in enumerate(n_conf_per_stream):
        if i != 0:
            continue
        print(nconf, streamids[i])
        offset = numpy.sum(n_conf_per_stream[:i]) if i > 0 else 0
        x = range(0, nconf*stepsize, stepsize)
        y = [polyakov_real[int(flowidx+n_flow*(k+offset))] for k in range(nconf)]
        ax.errorbar(x, y, fmt='.', fillstyle='full', markersize=3, mew=0, label=streamids[i], alpha=0.5)
        errs = []
        means = []
        for j in range(0,nconf):
            errs.append(numpy.std(y[j:]))
            means.append(numpy.mean(y[j:]))
        errs = numpy.asarray(errs)
        means = numpy.asarray(means)
        ax2.errorbar(x, errs, fmt='.', fillstyle='full', markersize=1, mew=0.25, lw=0.5, elinewidth=0.25, capsize=1.2, alpha=0.5, label=streamids[i])
        ax3.errorbar(x, means, fmt='.', fillstyle='full', markersize=1, mew=0.25, lw=0.5, elinewidth=0.25, capsize=1.2, alpha=0.5, label=streamids[i])

    legendoptions = dict(edgecolor='none', fancybox=False, facecolor="w", columnspacing=0.1,
                   labelspacing=0.1, borderpad=0.1, handletextpad=0.4, handlelength=1)
    ax.legend(**legendoptions)
    ax2.legend(**legendoptions)
    ax3.legend(**legendoptions)
    ax.axvline(x=2000, **lpd.verticallinestyle)
    ax2.axvline(x=2000, **lpd.verticallinestyle)
    ax3.axvline(x=2000, **lpd.verticallinestyle)
    ax.set_title(r'$'+str(ns)+r'^3 \times '+str(nt)+r'$, ms/5, HISQ', x=0.5)
    ax2.set_title(r'$'+str(ns)+r'^3 \times '+str(nt)+r'$')
    ax3.set_title(r'$' + str(ns) + r'^3 \times ' + str(nt) + r'$')
    fig.savefig(outputfolder+"/"+str(n_discarded)+"_polyakovloop.pdf")
    fig2.savefig(outputfolder+"/"+str(n_discarded)+"_polyakovloop_err.pdf")
    fig3.savefig(outputfolder+"/"+str(n_discarded)+"_polyakovloop_mean.pdf")



if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
