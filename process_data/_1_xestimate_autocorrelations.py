#!/usr/bin/env python3

import lib_process_data as lpd
import numpy
import _2_reduce_data as rd
from latqcdtools.statistics import statistics as stat
import matplotlib


def myapp(vec, val):
    # print(val)
    vec.append(val)


def get_equally_spaced_timeseries(x, y, MC_stepsize):
    # clean data such that trajectory spacing is constant=args.MC_stepsize.
    # 1. filter out confnums that are not divisible by 5 TODO fix this hard code
    # 2. average pairs of confs spaced out by only by 5 traj to get one conf with spacing args.MC_stepsize.
    # 3. whenever there is one conf missing, simply insert the average between the previous and next conf instead.
    # (this only happens very rarely, for example when some data file was corrupted).
    mask = numpy.asarray(numpy.mod(x, 5), dtype=bool)
    x = x[~mask]
    y = y[~mask]

    y_bin = []
    x_bin = []
    i = 0
    while i < len(x) - 1:
        diff = x[i + 1] - x[i]
        xi_is_multiple_of_ten = numpy.mod(x[i], MC_stepsize) == 0
        diff_is_multiple_of_ten = numpy.mod(diff, MC_stepsize) == 0
        if xi_is_multiple_of_ten and diff == 5:

            val = (y[i])  # TODO use mean here (y[i+1]+y[i])/2  ?
            y_bin.append(val)

            myapp(x_bin, x[i])
            if i < len(x) - 2 and x[i + 2] - x[i + 1] == MC_stepsize:
                # print("add missing", end="->")
                myapp(x_bin, x[i + 1] + 5)
                y_bin.append(y[i + 2])
                i += 1
            i += 2
        elif xi_is_multiple_of_ten and diff == MC_stepsize:
            myapp(x_bin, x[i])
            y_bin.append(y[i])
            i += 1
        elif diff > MC_stepsize and xi_is_multiple_of_ten and diff_is_multiple_of_ten:
            # print(x[i], "WARNa: diff greater args.MC_stepsize:", diff, end = " ---> ")
            nconfs_missing = int(diff / MC_stepsize) - 1
            for j in range(1, nconfs_missing + 1):
                myapp(x_bin, x[i] + j * MC_stepsize)
                y_bin.append((y[i + 1] + y[i]) / 2)
            i += 1
        elif diff > MC_stepsize and not xi_is_multiple_of_ten and diff_is_multiple_of_ten:
            # print("WARNb: diff greater args.MC_stepsize:", diff, end = " ---> ")
            nconfs_missing = int(diff / MC_stepsize)
            for j in range(1, nconfs_missing + 1):
                myapp(x_bin, x[i] + 5 + (j - 1) * MC_stepsize)
                y_bin.append((y[i + 1] + y[i]) / 2)
            i += 1
        elif diff > MC_stepsize and not xi_is_multiple_of_ten and not diff_is_multiple_of_ten:
            # print("WARNc: diff greater args.MC_stepsize:", diff, end = " ---> ")
            nconfs_missing = diff // MC_stepsize
            for j in range(1, nconfs_missing + 1):
                myapp(x_bin, x[i] + 5 + (j - 1) * MC_stepsize)
                y_bin.append((y[i + 1] + y[i]) / 2)
            i += 1
        elif diff > MC_stepsize and xi_is_multiple_of_ten and not diff_is_multiple_of_ten:
            # print("WARNd: diff greater args.MC_stepsize:", diff, end=" ---> ")
            nconfs_missing = diff // MC_stepsize
            myapp(x_bin, x[i])
            y_bin.append(y[i])
            for j in range(1, nconfs_missing + 1):
                myapp(x_bin, x[i] + j * MC_stepsize)
                y_bin.append((y[i + 1] + y[i]) / 2)
            i += 1
        elif i > 0 and x[i] - x_bin[-1] > MC_stepsize:
            print("WARN: There is a large gap in the time series! x[i]=", x[i], "x_bin[i-1]", x_bin[-1])
            i += 1
        else:
            # print("skipping", x[i], "diff=", diff)
            i += 1

    if x_bin[-1] - x_bin[0] != (len(x_bin) - 1) * MC_stepsize:
        print("ERROR: something went wrong when filtering the data, the MC time stepsize is not constant=args.MC_stepsize:")
        print("x_bin[-1]-x_bin[0]=", x_bin[-1] - x_bin[0], "but (len(x_bin)-1)*args.MC_stepsize=", (len(x_bin) - 1) * MC_stepsize)

    return numpy.asarray(x_bin), numpy.asarray(y_bin)


def get_start_and_tauint(x, y, min_conf, MC_stepsize, streamid, flowidx, tauidx, blocksize, tpickmax_increment):

    y = y[:, flowidx, tauidx]
    streamoffset = numpy.argmin(numpy.abs(x-min_conf*MC_stepsize))
    nconftotal = len(x)
    best_start = streamoffset
    best_nconfeff = 0
    bestvals = [0, 0, 0, 0]
    stepsize = blocksize
    for start in range(streamoffset, nconftotal-stepsize, stepsize):
        nconf = nconftotal-start
        diff = 0
        tpickmax = 10
        too_small = False
        while diff < 1:
            nblocks = int(nconf/tpickmax)
            if nblocks < 3:
                too_small = True
                break

            tau_int, tau_inte, tau_intbias, itpick = stat.getTauInt(y[start:], nblocks, tpickmax, acoutfileName='acor.d', showPlot=False)

            diff = tpickmax-itpick
            tpickmax += tpickmax_increment

        nconf_eff = nconf/(tau_int+tau_intbias)  # maximize this quantity
        if nconf_eff > best_nconfeff:
            best_nconfeff = nconf_eff
            best_start = start
            bestvals = [tau_int, tau_inte, tau_intbias, itpick]
            best_too_small = too_small

    print("stream=", streamid, ", nconf=[", int(x[best_start]/MC_stepsize), ",", int(x[-1]/MC_stepsize), "]", " ---> tau_int=", r'{0:.0f}'.format(bestvals[0]), "+-", r'{0:.0f}'.format(bestvals[1]),
          " (+", r'{0:.0f}'.format(bestvals[2]), ") at n=", bestvals[3], ", nconfeff=", int(best_nconfeff), sep="", end="")

    if best_too_small:
        print(", WARN: data set too small to find a decrease of tau_int inside the largest possible jackknife block. very small number of blocks (3).")
    else:
        print("")

    return best_start, tau_int, tau_inte, tau_intbias, itpick


def plot_MC_time(x, y, flowidx, dataindex, results, args, flowtime, nt, ns, beta):

    results = numpy.asarray(results)
    n_streams = results.shape[0]

    xlabel = r'ntraj/'+str(args.MC_stepsize)
    ylabel = r'$\mathrm{Re}[U_{(\beta,0)}]$'
    ylabelpos = (0.01, 0.99)
    fig, ax, _ = lpd.create_figure(xlabel=xlabel, ylabel=ylabel, xlabelpos=(0.99, 0.01), ylabelpos=ylabelpos, figsize=(5, 2.5))
    ax.set_xlabel(ax.get_xlabel(), fontsize=8, rotation=0, horizontalalignment='right', verticalalignment='bottom')
    ax.set_ylabel(ax.get_ylabel(), fontsize=8, rotation=0, horizontalalignment='left', verticalalignment='top')

    for k in range(n_streams):

        thisx = x[k]
        thisy = y[k][:, flowidx, dataindex]

        best_start = int(results[k][0])
        bestvals = results[k][1:]

        if k in args.show_id:
            ref = ax.errorbar(thisx[:best_start + 1] / args.MC_stepsize, thisy[:best_start + 1], fmt='-', lw=0.4, fillstyle='full', markersize=2.5, mew=0, alpha=0.25)
            label = str(k) + ", " + r'{0:.0f}'.format(bestvals[0]) + ", " + str(int(thisx[best_start] / args.MC_stepsize))
            ax.errorbar(thisx[best_start:] / args.MC_stepsize, thisy[best_start:], fmt='-', lw=0.4, fillstyle='full', markersize=2.5, mew=0, alpha=1,
                        color=ref.lines[0].get_color(), label=label)

    legendoptions = dict(edgecolor='none', fancybox=False, facecolor="w", columnspacing=0.1,
                         labelspacing=0.1, borderpad=0.1, handletextpad=0.4, handlelength=1, loc="center left", bbox_to_anchor=[1, 0.5])
    ax.legend(**legendoptions, title=r'id, $\tau_{int}$, start')
    if args.vlines:
        for x in args.vlines:
            ax.axvline(x=x, **lpd.verticallinestyle)
    ax.set_title(r'$' + str(ns) + r'^3 \times ' + str(nt) + r'$, $\beta=' + str(beta) + '$,' + r'$\sqrt{8\tau_F}/a=' + r'{0:.3f}'.format(
        numpy.sqrt(8 * flowtime)) + r', \sqrt{8\tau_F}T \approx ' + r'{0:.2f}'.format(
        numpy.sqrt(8 * flowtime) / nt) + r'$' + ', ms/5, HISQ', x=0.5, fontsize=8)

    ax.set_xlim(xmin=-50)

    thisylims = ax.get_ylim()
    diff = thisylims[1] - thisylims[0]
    ax.set_ylim(numpy.fmax(0, thisylims[0] - diff / 2), numpy.fmin(0.5,thisylims[1] + diff / 2))

    outputfolder = lpd.get_plot_path(args.qcdtype, args.corr, args.conftype, args.basepath_plot) if not args.outputpath else args.outputpath
    lpd.create_folder(outputfolder)
    fig.savefig(outputfolder + "/" + "polyakovloop_MCtime" + args.suffix + ".pdf")
    matplotlib.pyplot.close(fig)


def compute_XX_corr(data):
    """
    function for bootstrap routine that computes an XX correlator (with XX=EE or BB) normalized by the polyakov loop. numerator data (i.e. --X--X--) is first index, polyakov loop is second index of data.
    """
    nt = data.shape[-1]-1
    numerator, denominator = numpy.split(data, [nt,], axis=2)
    numerator_mean = numpy.mean(numerator, axis=0)
    denominator_mean = numpy.mean(denominator, axis=0)
    XX = numerator_mean/denominator_mean

    return XX


def main():
    parser, requiredNamed = lpd.get_parser()
    requiredNamed.add_argument('--conftype', help='format example: s096t32_b0824900_m002022_m01011', type=str, required=True)
    parser.add_argument('--outputpath', help='where to store the plot', type=str)
    parser.add_argument('--show_id', help="list of indices that shall be plotted. e.g. 0 1 2 3 for the first 4 streams. default=show all streams.", type=int, nargs='*', default=None)
    parser.add_argument('--vlines', help="where to plot vlines", nargs='*', type=float)
    parser.add_argument('--suffix', help="suffix for output file", default="")
    parser.add_argument('--basepath', type=str, default="")
    parser.add_argument('--basepath_plot', type=str, default="")
    parser.add_argument('--min_conf', help="ignore data below this conf number. nargs needs to be equal to nstreams.", nargs='*', default=[0,], type=int)
    parser.add_argument('--MC_stepsize', default=10, type=int)
    parser.add_argument('--blocksize', default=50, type=int, help="iterate over removing i*<blocksize> configurations from the beginning, "
                                                                   "then estimate tau_int and see if it improves the number of effective configurations.")
    parser.add_argument('--tpickmax_increment', default=2, type=int, help="if the jackknife bin size (to estimate the tau_int error and bias) is too small, "
                                                                           "then we never observe a decrease in the tau_int estimate and may underestimate tau_int."
                                                                           "so we keep increasing the binsize by <tpickmax_increment> and try again until we do.")
    parser.add_argument('--flowradius', default=0.075, type=float, help="at which flowtime tau_int should be calculated")
    args = parser.parse_args()

    beta, ns, nt, nt_half = lpd.parse_conftype(args.conftype)

    # load and reorganize data
    flow_times, n_flow, n_datafiles, n_streams, n_files_per_stream, XX_data, confnums = rd.load_merged_data(args.qcdtype, args.corr, args.conftype, args.basepath, None)
    polyakov_real = XX_data[1]
    numerator_real = XX_data[0]

    # TODO add option to just resample without estimating autocorrelations
    data = numpy.concatenate((numerator_real, polyakov_real), axis=2)
    index_selection = dict(polyakovloop=-1, numerator_at_largest_tau=-2)
    dataindex = index_selection["polyakovloop"]

    if args.show_id is None:
        args.show_id = range(n_streams)
    if args.min_conf == (0,):
        args.min_conf = [0 for _ in range(n_streams)]

    flowradii = numpy.sqrt(8*flow_times)/nt
    flowidx = numpy.argmin(numpy.abs(flowradii-args.flowradius))

    # separate data into streams
    x_eq_spaced_array = []
    y_eq_spaced_array = []
    results = []
    offset = 0

    y_binned = []

    for k in range(n_streams):
        x = numpy.asarray(confnums[k])
        y = data[offset:offset + n_files_per_stream[k]]
        offset += n_files_per_stream[k]

        x_eq_spaced, y_eq_spaced = get_equally_spaced_timeseries(x, y, args.MC_stepsize)
        x_eq_spaced_array.append(x_eq_spaced)
        y_eq_spaced_array.append(y_eq_spaced)
        best_start, tau_int, tau_inte, tau_intbias, itpick = get_start_and_tauint(x_eq_spaced, y_eq_spaced, args.min_conf[k], args.MC_stepsize, k, flowidx, dataindex, args.blocksize, args.tpickmax_increment)
        results.append([best_start, tau_int, tau_inte, tau_intbias, itpick])

        # reverse order
        y_eq_spaced = numpy.flip(y_eq_spaced, axis=0)

        # bin data
        ndata = len(y_eq_spaced)-best_start
        binlength = int(2*numpy.ceil(tau_int+tau_intbias))
        nbins = int(ndata/binlength)
        for b in range(nbins):
            y_binned.append(numpy.mean(y_eq_spaced[b*binlength:(b+1)*binlength], axis=0))

    y_binned = numpy.asarray(y_binned)
    print(args.conftype, "reduced ", n_datafiles, "to ", len(y_binned))

    plot_MC_time(x_eq_spaced_array, y_eq_spaced_array, flowidx, dataindex, results, args, flow_times[flowidx], nt, ns, beta)

    outputfolder = lpd.get_merged_data_path(args.qcdtype, args.corr, args.conftype, args.basepath)
    rd.resample_and_save_data(compute_XX_corr, y_binned, 1000, len(y_binned), outputfolder+"/"+args.corr, flow_times, args.qcdtype, args.conftype, args.corr, nt, False, 20, 0)

    print("done", args.conftype)


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
