#!/usr/bin/env python3

import lib_process_data as lpd
import numpy
from latqcdtools.statistics import statistics as stat
from latqcdtools.statistics import bootstr
import matplotlib


def myapp(vec, val):
    # print(val)
    vec.append(val)


def resample_and_save_data(reduce_function, XX_data, n_samples, n_datafiles, file_prefix, flow_times, qcdtype, conftype, corr, nt, BB_renorm, nproc, conf_axis):
    XX_samples, XX, XX_err = stat.bootstr(reduce_function, XX_data, numb_samples=n_samples, sample_size=n_datafiles, conf_axis=conf_axis, return_sample=True, same_rand_for_obs=True, parallelize=True, nproc=nproc, seed=0, err_by_dist=True)
    XX_samples = numpy.asarray(XX_samples)
    print(XX_samples.shape)
    # XX:     this is the bootstrap estimate for the mean of the correlator
    # XX_err: this is the bootstrap estimate for the error of the mean of the correlator

    # save samples
    numpy.save(lpd.print_var("write", file_prefix + "_" + conftype + "_samples.npy"), XX_samples)

    # add renormalization factor for BB correlator. errors for it are extremely small. FIXME abstract this
    if BB_renorm and corr == "BB":
        path = "/work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/coupling/"
        tfT2, Z2 = numpy.loadtxt(path + "Z2_" + str(nt) + ".dat", unpack=True)
        for i in range(XX.shape[0]):
            for j in range(XX.shape[1]):
                XX[i, j] *= Z2[i]
                XX_err[i, j] *= Z2[i]

    # write XX and XX_err to file
    with open(lpd.print_var("write", file_prefix+"_"+conftype+".dat"), 'w') as outfile:
        outfile.write('# bootstrap mean of '+str(n_datafiles)+' measurements of '+corr+' correlator for '+qcdtype+' '+conftype+'\n')
        outfile.write('# rows correspond to flow times, columns to dt = {1, ... , Ntau/2}\n')
        lpd.write_flow_times(outfile, flow_times)
        numpy.savetxt(outfile, XX)

    with open(lpd.print_var("write", file_prefix+"_err_"+conftype+".dat"), 'w') as outfile:
        outfile.write('# bootstrap mean err of '+str(n_datafiles)+' measurements of '+corr+' correlator '+qcdtype+' '+conftype+'\n')
        outfile.write('# rows correspond to flow times, columns to dt = {1, ... , Ntau/2}\n')
        lpd.write_flow_times(outfile, flow_times)
        numpy.savetxt(outfile, XX_err)

    # for each tau, compute flow correlation matrix.
    XX_samples = numpy.swapaxes(XX_samples, 0, 2)  # change data structure such that numpy.cov understands it
    print(XX_samples.shape)
    pcov = []
    for i in range(int(nt/2)):
        pcov.append(numpy.corrcoef(XX_samples[i]))
    pcov = numpy.asarray(pcov)
    numpy.save(lpd.print_var("write", file_prefix + "_flow_cov_" + conftype + ".npy"), pcov)
    print(pcov.shape)


# TODO incorporate here how the data should be binned?
def load_merged_data(qcdtype, corr, conftype, basepath, n_discard_per_stream=None, only_metadata=False):
    """
    load and reorganize data

    organize the data in the format for analysistoolbox general jackknife routine
    XX_numerator_real and polyakov_real: first index = different measurements, second index = flowtime, third index = tauT
    XX_data: first index [0]=XX_numerator_real, [1]=polyakov_real, second index = different measurements (-> pairs of the form (XX_numerator_real[x], XX_polyakov_real[x]))
    the bootstrap leaves out pairs of (XX_numerator, Polyakov), this is done by the last argument (1), which specifies the axis on which the data pairs lie.
    """

    inputfolder = lpd.get_merged_data_path(qcdtype, corr, conftype, basepath)

    flow_times = numpy.loadtxt(lpd.print_var("read", inputfolder+"flowtimes_"+conftype+".dat"))
    n_flow = len(flow_times)

    metadata = numpy.loadtxt(lpd.print_var("read", inputfolder + "n_datafiles_" + conftype + ".dat"))
    n_datafiles, n_streams = [int(i) for i in metadata[0:2]]
    if n_discard_per_stream is None:
        print("INFO: setting n_discard_per_stream to 0 for all streams")
        n_discard_per_stream = [0 for _ in range(n_streams)]

    n_files_per_stream = [int(i) for i in metadata[2:2+n_streams]]
    tmp = [int(i) for i in metadata[2+n_streams:]]
    confnums = []
    offset = 0
    for i in range(n_streams):
        confnums.append(tmp[offset+n_discard_per_stream[i]:offset+n_files_per_stream[i]])
        offset += n_files_per_stream[i]

    if n_streams != len(n_discard_per_stream):
        print("ERROR: n_streams=", n_streams, "but len(n_discard_per_stream)=", len(n_discard_per_stream), "(should be equal)")
        exit(1)

    if not only_metadata:
        XX_numerator_real_tmp = numpy.load(lpd.print_var("read", inputfolder+corr+"_real_"+conftype+"_merged.npy"))
        polyakov_real_tmp = numpy.load(lpd.print_var("read", inputfolder + "polyakov_real_" + conftype + "_merged.npy"))

        # discard data that is not thermalized
        slices = []
        counter_up = 0
        counter_dn = n_discard_per_stream[0]
        for i in range(len(n_files_per_stream)):
            if i > 0:
                counter_dn = counter_up + n_discard_per_stream[i]
            counter_up += n_files_per_stream[i]
            slices.append(slice(counter_dn, counter_up))

        XX_numerator_real = numpy.concatenate([XX_numerator_real_tmp[thisslice] for thisslice in slices])
        polyakov_real = numpy.expand_dims(numpy.concatenate([polyakov_real_tmp[thisslice] for thisslice in slices]), axis=2)
        XX_data = (XX_numerator_real, polyakov_real)
    else:
        XX_data = None

    print("readin done")


    return flow_times, n_flow, n_datafiles, n_streams, n_files_per_stream, XX_data, confnums




def get_equally_spaced_timeseries(x, y, MC_stepsize):
    # clean data such that trajectory spacing is constant=args.MC_stepsize.
    # 1. filter out confnums that are not divisible by *half* the MC stepsize
    # 2. average pairs of confs spaced out by only by *half* the stepsize to get one conf with spacing one full MC_stepsize.
    # 3. whenever there is one conf missing, simply insert the average between the previous and next conf instead.
    # (this should only happen very rarely, for example when some data file was corrupted).
    half_MC_stepsize = int(MC_stepsize / 2)

    mask = numpy.asarray(numpy.mod(x, int(half_MC_stepsize)), dtype=bool)
    x = x[~mask]
    y = y[~mask]

    y_bin = []
    x_bin = []
    i = 0
    while i < len(x) - 1:
        diff = x[i + 1] - x[i]
        xi_is_multiple_of_stepsize = numpy.mod(x[i], MC_stepsize) == 0
        diff_is_multiple_of_stepsize = numpy.mod(diff, MC_stepsize) == 0
        if xi_is_multiple_of_stepsize and diff == half_MC_stepsize:

            val = (y[i])  # TODO use mean here (y[i+1]+y[i])/2  ?
            y_bin.append(val)

            myapp(x_bin, x[i])
            if i < len(x) - 2 and x[i + 2] - x[i + 1] == MC_stepsize:
                print("add missing", end=", ")
                myapp(x_bin, x[i + 1] + half_MC_stepsize)
                y_bin.append(y[i + 2])
                i += 1
            i += 2
        elif xi_is_multiple_of_stepsize and diff == MC_stepsize:
            myapp(x_bin, x[i])
            y_bin.append(y[i])
            i += 1
        elif diff > MC_stepsize and xi_is_multiple_of_stepsize and diff_is_multiple_of_stepsize:
            print(x[i], "WARNa: diff greater MC_stepsize:", diff, end=", ")
            nconfs_missing = int(diff / MC_stepsize) - 1
            for j in range(1, nconfs_missing + 1):
                myapp(x_bin, x[i] + j * MC_stepsize)
                y_bin.append((y[i + 1] + y[i]) / 2)
            i += 1
        elif diff > MC_stepsize and not xi_is_multiple_of_stepsize and diff_is_multiple_of_stepsize:
            print("WARNb: diff greater MC_stepsize:", diff, end=", ")
            nconfs_missing = int(diff / MC_stepsize)
            for j in range(1, nconfs_missing + 1):
                myapp(x_bin, x[i] + half_MC_stepsize + (j - 1) * MC_stepsize)
                y_bin.append((y[i + 1] + y[i]) / 2)
            i += 1
        elif diff > MC_stepsize and not xi_is_multiple_of_stepsize and not diff_is_multiple_of_stepsize:
            print("WARNc: diff greater MC_stepsize:", diff, end=", ")
            nconfs_missing = diff // MC_stepsize
            for j in range(1, nconfs_missing + 1):
                myapp(x_bin, x[i] + half_MC_stepsize + (j - 1) * MC_stepsize)
                y_bin.append((y[i + 1] + y[i]) / 2)
            i += 1
        elif diff > MC_stepsize and xi_is_multiple_of_stepsize and not diff_is_multiple_of_stepsize:
            print("WARNd: diff greater MC_stepsize:", diff, end=", ")
            nconfs_missing = diff // MC_stepsize
            myapp(x_bin, x[i])
            y_bin.append(y[i])
            for j in range(1, nconfs_missing + 1):
                myapp(x_bin, x[i] + j * MC_stepsize)
                y_bin.append((y[i + 1] + y[i]) / 2)
            i += 1
        elif i > 0 and x[i] - x_bin[-1] > MC_stepsize:
            print("WARN: There is a large gap in the time series! x[i]=", x[i], ", x_bin[i-1]=", x_bin[-1], sep="")
            i += 1
        else:
            print("skipping ", x[i], ", diff=", diff, sep="", end=", ")
            i += 1

    # don't forget the last entry
    if x[-1] - x_bin[-1] == MC_stepsize:
        y_bin.append(y[-1])
        myapp(x_bin, x[-1])
    print("")

    if x_bin[-1] - x_bin[0] != (len(x_bin) - 1) * MC_stepsize:
        print("ERROR: something went wrong when filtering the data, the MC time stepsize is not constant=MC_stepsize:")
        print("x_bin[-1]-x_bin[0]=", x_bin[-1] - x_bin[0], "but (len(x_bin)-1)*MC_stepsize=", (len(x_bin) - 1) * MC_stepsize)

    return numpy.asarray(x_bin), numpy.asarray(y_bin)


def find_reliable_tauint(index_max, index_start, y, tpickmax_increment):
    nconf = index_max - index_start
    diff = 0
    tpickmax = 10
    too_small = False
    while diff < 1:
        nblocks = int(nconf / tpickmax)
        if nblocks < 3:
            too_small = True
            break

        tau_int, tau_inte, tau_intbias, itpick = stat.getTauInt(y[index_start:], nblocks, tpickmax, acoutfileName='acor.d', showPlot=False)

        diff = tpickmax - itpick
        tpickmax += tpickmax_increment

    return tau_int, tau_inte, tau_intbias, itpick, too_small, nconf


def get_start_and_tauint(x, y, index_offset, MC_stepsize, streamid, blocksize, tpickmax_increment, include_bias):
    index_max = len(x)
    best_start_index = index_offset
    best_nconfeff = 0
    bestvals = [0, 0, 0, 0]
    for index_start in range(index_offset, index_max - blocksize, blocksize):
        tau_int, tau_inte, tau_intbias, itpick, too_small, nconf = find_reliable_tauint(index_max, index_start, y, tpickmax_increment)

        # maximize this quantity
        if include_bias:
            # this can somewhat reliably find the point where the time series is stationary (ie, fluctuates around a fixed value).
            # don't use if you know that your time series is already stationary, as the tau_intbias estimate is a bit unreliable.
            nconf_eff = nconf / (tau_int + tau_intbias)
        else:
            nconf_eff = nconf / tau_int

        if nconf_eff > best_nconfeff:
            best_nconfeff = nconf_eff
            best_start_index = index_start
            bestvals = [tau_int, tau_inte, tau_intbias, itpick]
            best_too_small = too_small

    print_result(streamid, x, best_start_index, MC_stepsize, tau_int, tau_inte, tau_intbias, itpick, index_max)

    if best_too_small:
        print_too_small_warning()
    else:
        print("")

    return best_start_index, tau_int, tau_inte, tau_intbias, itpick


def print_result(streamid, x, best_start_index, MC_stepsize, tau_int, tau_inte, tau_intbias, itpick, index_max):
    print("stream=", streamid, ", nconf=[", int(x[best_start_index] / MC_stepsize), ",", int(x[-1] / MC_stepsize), "]", " ---> tau_int=",
          r'{0:.0f}'.format(tau_int), "+-", r'{0:.0f}'.format(tau_inte),
          " (+", r'{0:.0f}'.format(tau_intbias), ") at n=", itpick, ", nconfeff=", int((index_max - best_start_index) / tau_int), sep="", end="")




def print_too_small_warning():
    print(", WARN: data set too small to find a decrease of tau_int inside the largest possible jackknife block. very small number of blocks (3).")


def plot_MC_time(x, y, flowidx, dataindex, results, args, flowtime, nt, ns, beta):
    results = numpy.asarray(results)
    n_streams = results.shape[0]

    xlabel = r'$n_\mathrm{traj}/ ' + str(args.MC_stepsize) + r'$'
    ylabel = r'$\mathrm{Re}\mathrm{Tr}[U_{(\beta,0)}]$'
    fig, ax, _ = lpd.create_figure(xlabel=xlabel, ylabel=ylabel, figsize="fullwidth_slim")

    for k in range(n_streams):

        thisx = x[k]
        thisy = y[k][:, flowidx, dataindex]

        best_start = int(results[k][0])
        bestvals = results[k][1:]

        if k in args.show_id:
            ref = ax.errorbar(thisx[:best_start + 1] / args.MC_stepsize, thisy[:best_start + 1], fmt='-', lw=0.4, fillstyle='full', markersize=2.5, mew=0,
                              alpha=0.25)
            label = str(k)  # + ", " + str(int(bestvals[0])) + ", " + str(int(thisx[best_start] / args.MC_stepsize)
            ax.errorbar(thisx[best_start:] / args.MC_stepsize, thisy[best_start:], fmt='-', lw=0.4, fillstyle='full', markersize=2.5, mew=0, alpha=1,
                        color=ref.lines[0].get_color(), label=label)

    legendoptions = dict(edgecolor='none', fancybox=False, facecolor="w", columnspacing=0.1,
                         labelspacing=0.1, borderpad=0.1, handletextpad=0.4, handlelength=1, loc="center left", bbox_to_anchor=[1, 0.5])
    ax.legend(**legendoptions, title=r'stream')  # , $\tau_{int}$, start
    if args.vlines:
        for x in args.vlines:
            ax.axvline(x=x, **lpd.verticallinestyle)
    ax.text(0.99, 0.99, r'$' + str(ns) + r'^3 \times ' + str(nt) + r'$, $\beta=' + str(beta) + '$,' + r'$\sqrt{8\tau_F}/a=' + r'{0:.1f}'.format(
        numpy.sqrt(8 * flowtime)) + r'$', ha='right', va='top', transform=ax.transAxes)  # + ', ms/5, HISQ'
    # r', \sqrt{8\tau_F}T \approx ' + r'{0:.2f}'.format(numpy.sqrt(8 * flowtime) / nt) +
    ax.set_xlim(xmin=-50)

    thisylims = ax.get_ylim()
    diff = thisylims[1] - thisylims[0]
    ax.set_ylim(numpy.fmax(0, thisylims[0] - diff / 4), numpy.fmin(0.5, thisylims[1] + diff / 4))

    outputfolder = lpd.get_plot_path(args.qcdtype, args.corr, args.conftype, args.basepath_plot) if not args.outputpath else args.outputpath
    lpd.create_folder(outputfolder)
    fig.savefig(outputfolder + "/" + "polyakovloop_MCtime" + args.suffix + ".pdf")
    matplotlib.pyplot.close(fig)


def compute_XX_corr(data):
    """
    function for bootstrap routine that computes an XX correlator (with XX=EE or BB) normalized by the polyakov loop. numerator data (i.e. --X--X--) is first index, polyakov loop is second index of data.
    """
    # TODO the way this is written FORCES same_rand_for_obs=true. maybe we dont want this.
    nt = data.shape[-1] - 1
    numerator, denominator = numpy.split(data, [nt, ], axis=2)
    numerator_mean = numpy.mean(numerator, axis=0)
    denominator_mean = numpy.mean(denominator, axis=0)
    XX = numerator_mean / denominator_mean

    return XX


def main():
    # TODO sort arguments according to their purpose

    # TODO add custom ylims here
    parser, requiredNamed = lpd.get_parser()
    requiredNamed.add_argument('--conftype', help='format example: s096t32_b0824900_m002022_m01011', type=str, required=True)
    parser.add_argument('--outputpath', help='where to store the plot', type=str)
    parser.add_argument('--show_id', help="list of indices that shall be plotted. e.g. 0 1 2 3 for the first 4 streams. default=show all streams.", type=int,
                        nargs='*', default=None)
    parser.add_argument('--vlines', help="where to plot vlines", nargs='*', type=float)
    parser.add_argument('--suffix', help="suffix for output file", default="")
    parser.add_argument('--basepath', type=str, default="")
    parser.add_argument('--basepath_plot', type=str, default="")
    parser.add_argument('--min_conf', help="ignore data below this conf number. nargs needs to be equal to nstreams.", nargs='*', default=[0, ], type=int)
    parser.add_argument('--already_equally_spaced', action="store_true", default=False, help="if data is already perfectly equally spaced across MC time.")
    parser.add_argument('--skip_binning', action="store_true", default=False,
                        help="dont bin the data according to tau_int, just resample. useful if you know the data is independent.")
    parser.add_argument('--include_bias', help="maximize nconf/(tau_int+bias) instead of nconf/tau_int. this can help you find a point where the "
                                               "chain is relatively stationary, but the bias estimate is generally a bit unreliable, so use with care.",
                        default=False, action="store_true")
    parser.add_argument('--MC_stepsize', default=10, type=int)
    parser.add_argument('--blocksize', default=0, type=int, help="iterate over removing i*<blocksize> configurations from the beginning, "
                                                                 "then estimate tau_int and see if it improves the number of effective configurations (nconf/tau_int).")
    parser.add_argument('--tpickmax_increment', default=2, type=int, help="if the jackknife bin size (to estimate the tau_int error and bias) is too small, "
                                                                          "then we never observe a decrease in the tau_int estimate and may underestimate tau_int."
                                                                          "so we keep increasing the binsize by <tpickmax_increment> and try again until we do.")
    parser.add_argument('--flowradius', default=0.075, type=float, help="at which flowtime tau_int should be calculated")
    parser.add_argument('--n_samples', default=10000, type=int, help="number of bootstrap samples to draw")
    parser.add_argument('--n_proc', default=20, type=int, help="number of processes for parallelization")

    args = parser.parse_args()

    beta, ns, nt, nt_half = lpd.parse_conftype(args.conftype)

    # TODO add option to choose between full correlator or only numerator or denominator

    # load and reorganize data
    flow_times, n_flow, n_datafiles, n_streams, n_files_per_stream, XX_data, confnums = load_merged_data(args.qcdtype, args.corr, args.conftype,
                                                                                                            args.basepath, None)
    polyakov_real = XX_data[1]
    numerator_real = XX_data[0]

    data = numpy.concatenate((numerator_real, polyakov_real), axis=2)
    index_selection = dict(polyakovloop=-1, numerator_at_largest_tau=-2)
    dataindex = index_selection["polyakovloop"]

    if args.show_id is None:
        args.show_id = range(n_streams)
    if args.min_conf == (0,):
        args.min_conf = [0 for _ in range(n_streams)]

    flowradii = numpy.sqrt(8 * flow_times) / nt
    flowidx = numpy.argmin(numpy.abs(flowradii - args.flowradius))

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

        if args.already_equally_spaced:
            x_eq_spaced = x
            y_eq_spaced = y
        else:
            x_eq_spaced, y_eq_spaced = get_equally_spaced_timeseries(x, y, args.MC_stepsize)
        x_eq_spaced_array.append(x_eq_spaced)
        y_eq_spaced_array.append(y_eq_spaced)

        index_offset = numpy.argmin(numpy.abs(x - args.min_conf[k] * args.MC_stepsize))

        if args.blocksize:
            print("trying to optimize effective number of independent measurements, that is nconf/tau_int.")
            best_start, tau_int, tau_inte, tau_intbias, itpick = get_start_and_tauint(
                x_eq_spaced, y_eq_spaced[:, flowidx, dataindex], index_offset, args.MC_stepsize, k, args.blocksize, args.tpickmax_increment, args.include_bias)
        else:
            tau_int, tau_inte, tau_intbias, itpick, too_small, nconf = find_reliable_tauint(len(y_eq_spaced), index_offset, y[:, flowidx, dataindex], args.tpickmax_increment)
            print_result(k, x, index_offset, args.MC_stepsize, tau_int, tau_inte, tau_intbias, itpick, len(y_eq_spaced))
            if too_small:
                print_too_small_warning()
            else:
                print("")
            best_start = index_offset

        results.append([best_start, tau_int, tau_inte, tau_intbias, itpick])

        # reverse order
        y_eq_spaced = numpy.flip(y_eq_spaced, axis=0)

        if args.skip_binning:  # don't bin
            for val in y_eq_spaced:
                y_binned.append(val)
        else:  # bin data
            ndata = len(y_eq_spaced) - best_start
            binlength = int(tau_int)
            nbins = int(ndata / binlength)
            for b in range(nbins):
                y_binned.append(numpy.mean(y_eq_spaced[b * binlength:(b + 1) * binlength], axis=0))

    y_binned = numpy.asarray(y_binned)
    if len(y_binned) < n_datafiles:
        print(args.conftype, "binned ", n_datafiles, "files into ", len(y_binned), "independent bins")

    plot_MC_time(x_eq_spaced_array, y_eq_spaced_array, flowidx, dataindex, results, args, flow_times[flowidx], nt, ns, beta)

    outputfolder = lpd.get_merged_data_path(args.qcdtype, args.corr, args.conftype, args.basepath)
    # TODO add BB renorm
    resample_and_save_data(compute_XX_corr, y_binned, args.n_samples, len(y_binned),
                              outputfolder + "/" + args.corr, flow_times, args.qcdtype, args.conftype, args.corr, nt, False, args.n_proc, 0)

    print("done", args.conftype)


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
