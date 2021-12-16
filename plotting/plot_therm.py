#!/usr/local/bin/python3.7m
import numpy
import lib_process_data as lpd
from latqcdtools import statistics as stat


def load_merged_data(qcdtype, corr, conftype):
    inputfolder = lpd.get_merged_data_path(qcdtype, corr, conftype)

    print("read  " + inputfolder + "flowtimes_" + conftype + ".dat")
    flow_times = numpy.loadtxt(inputfolder + "flowtimes_" + conftype + ".dat")
    n_flow = len(flow_times)

    print("read  " + inputfolder + "n_datafiles_" + conftype + ".dat")
    metadata = numpy.loadtxt(inputfolder + "n_datafiles_" + conftype + ".dat")
    n_datafiles, n_streams = [int(i) for i in metadata[0:2]]
    n_conf_per_stream = [int(i) for i in metadata[2:]]

    print("read  " + inputfolder + "polyakov_real_" + conftype + "_merged.dat")
    polyakov_real = numpy.loadtxt(inputfolder + "polyakov_real_" + conftype + "_merged.dat")

    EE_real = numpy.loadtxt(inputfolder + "EE_real_" + conftype + "_merged.dat")
    print("readin done")

    return flow_times, n_flow, n_datafiles, n_streams, n_conf_per_stream, polyakov_real, EE_real


def main():

    parser, requiredNamed = lpd.get_parser()
    requiredNamed.add_argument('--conftype', help='format example: s096t32_b0824900_m002022_m01011', type=str, required=True)
    parser.add_argument('--outputpath', help='where to store the plot', type=str)
    args = parser.parse_args()
    beta, ns, nt, nt_half = lpd.parse_conftype(args.conftype)

    outputfolder = lpd.get_plot_path(args.qcdtype, args.corr, args.conftype) if not args.outputpath else args.outputpath

    flow_times, n_flow, n_datafiles, n_streams, n_conf_per_stream, polyakov_real, EE_real = load_merged_data(args.qcdtype, args.corr, args.conftype)

    xlabel = r'no. traj.'
    ylabel = r'$\mathrm{Re}($P.loop$)$,$\tau_F/a^2='+r'{0:.3f}'.format(flow_times[-1])+r', \sqrt{8\tau_F}T \approx 0.3$'
    stepsize = 5  # number of trajectories between written out confs

    streamids = ("a", "b")
    flowidx = 60
    fig, ax, plots = lpd.create_figure(xlabel=xlabel, ylabel=ylabel, xlabelpos=(0.5, -0.1), ylabelpos=(0.1, 0.9), UseTex=True)  # ylims=(0, 5), xlims=(0,5),

    # print([polyakov_real[int(flowidx+k*n_flow)] for k in range(int(1200))])

    # ax.errorbar(range(0, n_datafiles), [polyakov_real[int(flowidx + k * n_flow)] for k in range(n_datafiles)]

    for i, nconf in enumerate(n_conf_per_stream):
        offset = n_conf_per_stream[i-1] if i > 0 else 0
        x = range(0, nconf*stepsize, stepsize)
        # y = [polyakov_real[int(flowidx+n_flow*(k+offset))] for k in range(nconf)]
        y = [EE_real[int(flowidx+n_flow*(k+offset))][-1] for k in range(nconf)]
        ax.errorbar(x, y, fmt='.', fillstyle='full', markersize=1, label=streamids[i])
    # ax.axvline(2500, **lpd.verticallinestyle)
    ax.set_title(r'$'+str(ns)+r'^3 \times '+str(nt)+r'$, ms/5, HISQ', x=0.5)
    ax.legend(loc="center right")
    fig.savefig(outputfolder+"/polyakovloop.pdf")

    fig2, ax2, plots = lpd.create_figure(xlabel=xlabel, ylabel=ylabel, xlabelpos=(0.5, -0.1), UseTex=True)

    for i, nconf in enumerate(n_conf_per_stream):
        offset = n_conf_per_stream[i - 1] if i > 0 else 0
        x = range(0, nconf, 1)

        y = [polyakov_real[int(flowidx + n_flow * (k + offset))] for k in range(nconf)]
        means = []
        errs = []
        for j in x:
            means.append(numpy.mean(y[j:]))
            errs.append(numpy.std(y[j:]))
        means = numpy.asarray(means)
        errs = numpy.asarray(errs)
        x = range(0, nconf*5, 5)
        ax2.errorbar(x, errs, fmt='.', fillstyle='full', markersize=1, mew=0.25, lw=0.5, elinewidth=0.25, capsize=1.2, label=streamids[i])
    # ax2.axvline(2500, **lpd.verticallinestyle)
    ax2.set_title(r'$'+str(ns)+r'^3 \times '+str(nt)+r'$')

    fig2.savefig(outputfolder+"/polyakovloop2.pdf")


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
