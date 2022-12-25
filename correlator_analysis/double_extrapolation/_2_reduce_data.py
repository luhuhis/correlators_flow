#!/usr/bin/env python3
import numpy
import lib_process_data as lpd
from latqcdtools.statistics import bootstr



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


def compute_XX_corr(data):
    """
    function for bootstrap routine that computes an XX correlator (with XX=EE or BB) normalized by the polyakov loop. numerator data (i.e. --X--X--) is first index, polyakov loop is second index of data.
    """
    numerator_mean = numpy.mean(data[0], axis=0)
    denominator_mean = numpy.mean(data[1], axis=0)  # median?
    XX = numerator_mean/denominator_mean
    return XX


def compute_only_numerator(data):
    return numpy.mean(data[0], axis=0), numpy.std(data[0], axis=0)


def compute_only_denominator(data):
    return numpy.mean(data[1], axis=0), numpy.std(data[1], axis=0)


def resample_and_save_data(reduce_function, XX_data, n_samples, n_datafiles, file_prefix, flow_times, qcdtype, conftype, corr, nt, BB_renorm, nproc, conf_axis):
    XX_samples, XX, XX_err = bootstr.bootstr(reduce_function, XX_data, numb_samples=n_samples, sample_size=n_datafiles, conf_axis=conf_axis, return_sample=True, same_rand_for_obs=True, parallelize=True, nproc=nproc, seed=0, err_by_dist=True)
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


def main():

    # parse cmd line arguments
    parser, requiredNamed = lpd.get_parser()
    requiredNamed.add_argument('--conftype', help="format: s096t20_b0824900 for quenched or s096t20_b0824900_m002022_m01011 for hisq", required=True)
    parser.add_argument('--part_obs', help="only compute part of the observable", choices=['numerator', 'polyakovloop'])
    parser.add_argument('--n_samples', help='number of bootstrap resamplings', type=int, default='1000')
    parser.add_argument('--BB_renorm', action="store_true", help="multiply corr with Z factor")
    parser.add_argument('--basepath', default="", type=str, help="where to look for the data")
    parser.add_argument('--n_discard', help="number of configurations, counted from the lowest conf_num, in each stream that should be ignored", type=int,
                        nargs='*', required=True)
    parser.add_argument('--nproc', default=10, type=int, help="how many processes to use for parallelization")
    args = parser.parse_args()

    beta, ns, nt, nt_half = lpd.parse_conftype(args.conftype)

    inputfolder = lpd.get_merged_data_path(args.qcdtype, args.corr, args.conftype, args.basepath)
    outputfolder = inputfolder

    flow_times, n_flow, n_datafiles, n_streams, n_files_per_stream, XX_data, confnums = load_merged_data(args.qcdtype, args.corr, args.conftype, args.basepath, args.n_discard)
    XX_data = numpy.asarray(XX_data)
    print(XX_data.shape)

    reduce_function = compute_XX_corr
    file_prefix = args.corr
    if args.part_obs:
        if args.part_obs == "numerator":
            reduce_function = compute_only_numerator
            file_prefix += "_numerator"
        if args.part_obs == "polyakovloop":
            reduce_function = compute_only_denominator
            file_prefix += "_polyakovloop"

    file_prefix = outputfolder+file_prefix

    # TODO add BB renorm Z factor in this step, not afterwards!
    # Z = 1
    # if args.BB_renorm and args.corr == "BB":
    #     Z = Z2[i]

    resample_and_save_data(reduce_function, XX_data, args.n_samples, n_datafiles, file_prefix, flow_times, args.qcdtype, args.conftype, args.corr, nt, args.BB_renorm, args.nproc, 1)

    print("done reducing data for "+args.conftype+" "+args.corr)


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
