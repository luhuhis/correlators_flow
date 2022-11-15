#!/usr/bin/env python3
import numpy
import lib_process_data as lpd
from latqcdtools.statistics import bootstr


def load_merged_data(qcdtype, corr, conftype, basepath):
    """
    load and reorganize data

    organize the data in the format for analysistoolbox general jackknife routine
    XX_numerator_real and polyakov_real: first index = different measurements, second index = flowtime, third index = tauT
    XX_data: first index [0]=XX_numerator_real, [1]=polyakov_real, second index = different measurements (-> pairs of the form (XX_numerator_real[x], XX_polyakov_real[x]))
    the bootstrap leaves out pairs of (XX_numerator, Polyakov), this is done by the last argument (1), which specifies the axis on which the data pairs lie.
    """
    inputfolder = lpd.get_merged_data_path(qcdtype, corr, conftype, basepath)
    nt, nt_half = lpd.parse_conftype(conftype)[2:]

    print("read  "+inputfolder+"flowtimes_"+conftype+".dat")
    flow_times = numpy.loadtxt(inputfolder+"flowtimes_"+conftype+".dat")
    n_flow = len(flow_times)

    print("read  " + inputfolder + "n_datafiles_" + conftype + ".dat")
    metadata = numpy.loadtxt(inputfolder + "n_datafiles_" + conftype + ".dat")
    n_datafiles, n_streams = [int(i) for i in metadata[0:2]]

    print("read  "+inputfolder+corr+"_real_"+conftype+"_merged.dat")
    XX_numerator_real_tmp = numpy.loadtxt(inputfolder+corr+"_real_"+conftype+"_merged.dat")
    XX_numerator_real = []
    for i in range(n_datafiles):
        XX_numerator_real.append(XX_numerator_real_tmp[i*n_flow:(i+1)*n_flow, :])

    print("read  "+inputfolder+"polyakov_real_"+conftype+"_merged.dat")
    polyakov_real_tmp = numpy.loadtxt(inputfolder+"polyakov_real_"+conftype+"_merged.dat")

    print("readin done")

    polyakov_real = []
    for i in range(n_datafiles):
        polyakov_real.append(polyakov_real_tmp[i*n_flow:(i+1)*n_flow])
    polyakov_real = [[nt_half*[j] for j in i] for i in polyakov_real]  # copy the polyakov into every cell in order to match shapes for analysistoolbox bootstrap

    XX_data = [XX_numerator_real, polyakov_real]

    return flow_times, n_flow, n_datafiles, n_streams, XX_data


def compute_XX_corr(data):
    """
    function for bootstrap routine that computes an XX correlator (with XX=EE or BB) normalized by the polyakov loop. numerator data (i.e. --X--X--) is first index, polyakov loop is second index of data.
    """
    numerator_mean = numpy.mean(data[0], axis=0)
    denominator_mean = numpy.mean(data[1], axis=0)
    XX = numerator_mean/denominator_mean
    return XX


def compute_only_numerator(data):
    return numpy.mean(data[0], axis=0), numpy.std(data[0], axis=0)


def compute_only_denominator(data):
    return numpy.mean(data[1], axis=0), numpy.std(data[1], axis=0)


def main():

    # parse cmd line arguments
    parser, requiredNamed = lpd.get_parser()
    requiredNamed.add_argument('--conftype', help="format: s096t20_b0824900 for quenched or s096t20_b0824900_m002022_m01011 for hisq", required=True)
    parser.add_argument('--part_obs', help="only compute part of the observable", choices=['numerator', 'polyakovloop'])
    parser.add_argument('--n_samples', help='number of bootstrap samples', type=int, default='1000')
    parser.add_argument('--BB_renorm', action="store_true", help="multiply corr with Z factor")
    parser.add_argument('--basepath', default="", type=str, help="where to look for the data")
    parser.add_argument('--nproc', default=10, type=int, help="how many processes to use for parallelization")
    args = parser.parse_args()

    beta, ns, nt, nt_half = lpd.parse_conftype(args.conftype)

    inputfolder = lpd.get_merged_data_path(args.qcdtype, args.corr, args.conftype, args.basepath)
    outputfolder = inputfolder

    flow_times, n_flow, n_datafiles, n_streams, XX_data = load_merged_data(args.qcdtype, args.corr, args.conftype, args.basepath)
    # TODO: add ndiscard here, and not in merge data
    XX_data = numpy.asarray(XX_data)
    print(XX_data.shape)
    tauTs = numpy.arange(1/nt, 0.501, 1/nt)
    n_samples = args.n_samples

    reduce_function = compute_XX_corr
    file_prefix = args.corr
    if args.part_obs:
        if args.part_obs == "numerator":
            reduce_function = compute_only_numerator
            file_prefix += "_numerator"
        if args.part_obs == "polyakovloop":
            reduce_function = compute_only_denominator
            file_prefix += "_polyakovloop"

    # TODO add BB renorm Z factor in this step, not afterwards!
    # Z = 1
    # if args.BB_renorm and args.corr == "BB":
    #     Z = Z2[i]

    # Use same_rand_for_obs=False in order to break up correlations between the observables
    XX_samples, XX, XX_err = bootstr.bootstr(reduce_function, XX_data, numb_samples=n_samples, sample_size=n_datafiles, conf_axis=1, return_sample=True, same_rand_for_obs=True, parallelize=True, nproc=args.nproc, seed=0, err_by_dist=True)
    XX_samples = numpy.asarray(XX_samples)
    print(XX_samples.shape)
    # XX:     this is the bootstrap estimate for the mean of the correlator
    # XX_err: this is the bootstrap estimate for the error of the mean of the correlator

    # save samples
    filename = outputfolder + file_prefix + "_" + args.conftype + "_samples.npy"
    print("write " + filename)
    numpy.save(filename, XX_samples)

    # add renormalization factor for BB correlator. errors for it are extremely small. FIXME abstract this
    if args.BB_renorm and args.corr == "BB":
        path = "/work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/coupling/"
        tfT2, Z2 = numpy.loadtxt(path + "Z2_" + str(nt) + ".dat", unpack=True)
        for i in range(XX.shape[0]):
            for j in range(XX.shape[1]):
                XX[i, j] *= Z2[i]
                XX_err[i, j] *= Z2[i]

    # write XX and XX_err to file
    filename = outputfolder+file_prefix+"_"+args.conftype+".dat"
    print("write "+filename)
    with open(filename, 'w') as outfile:
        outfile.write('# bootstrap mean of '+str(n_datafiles)+' measurements of '+args.corr+' correlator for '+args.qcdtype+' '+args.conftype+'\n')
        outfile.write('# rows correspond to flow times, columns to dt = {1, ... , Ntau/2}\n')
        lpd.write_flow_times(outfile, flow_times)
        numpy.savetxt(outfile, XX)

    filename = outputfolder+file_prefix+"_err_"+args.conftype+".dat"
    print("write "+filename)
    with open(filename, 'w') as outfile:
        outfile.write('# bootstrap mean err of '+str(n_datafiles)+' measurements of '+args.corr+' correlator '+args.qcdtype+' '+args.conftype+'\n')
        outfile.write('# rows correspond to flow times, columns to dt = {1, ... , Ntau/2}\n')
        lpd.write_flow_times(outfile, flow_times)
        numpy.savetxt(outfile, XX_err)

    # for each tau, compute flow correlation matrix.
    XX_samples = numpy.swapaxes(XX_samples, 0, 2)  # change data structure such that numpy.cov understands it
    pcov = []
    for i in range(nt_half):
        pcov.append(numpy.corrcoef(XX_samples[i]))
    pcov = numpy.asarray(pcov)
    filename = outputfolder + file_prefix + "_flow_cov_" + args.conftype + ".npy"
    print("write " + filename)
    numpy.save(filename, pcov)
    print(pcov.shape)

    print("done reducing data for "+args.conftype+" "+args.corr)


if __name__ == '__main__':
    main()
    lpd.save_script_call()
