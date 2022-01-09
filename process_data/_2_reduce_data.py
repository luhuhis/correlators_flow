#!/usr/local/bin/python3.7m -u
import numpy
import lib_process_data as lpd
from latqcdtools import bootstr


def load_merged_data(qcdtype, corr, conftype):
    """
    load and reorganize data

    organize the data in the format for analysistoolbox general jackknife routine
    XX_numerator_real and polyakov_real: first index = different measurements, second index = flowtime, third index = tauT
    XX_data: first index [0]=XX_numerator_real, [1]=polyakov_real, second index = different measurements (-> pairs of the form (XX_numerator_real[x], XX_polyakov_real[x]))
    the bootstrap leaves out pairs of (XX_numerator, Polyakov), this is done by the last argument (1), which specifies the axis on which the data pairs lie.
    """
    inputfolder = lpd.get_merged_data_path(qcdtype, corr, conftype)
    nt, nt_half = lpd.parse_conftype(conftype)[2:]

    print("read  "+inputfolder+"flowtimes_"+conftype+".dat")
    flow_times = numpy.loadtxt(inputfolder+"flowtimes_"+conftype+".dat")
    n_flow = len(flow_times)

    print("read  " + inputfolder + "n_datafiles_" + conftype + ".dat")
    metadata = numpy.loadtxt(inputfolder + "n_datafiles_" + conftype + ".dat")
    n_datafiles, n_streams, n_discarded = [int(i) for i in metadata[0:3]]
    n_conf_per_stream = [int(i) for i in metadata[3:]]

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
    normalization_factor = 1  # the multiplicities are now correctly taken care of inside of the ParallelGPUCode
    numerator_mean = numpy.mean(data[0], axis=0)
    denominator_mean = numpy.mean(data[1], axis=0)
    XX = numerator_mean/normalization_factor / denominator_mean
    numerator_std_dev = numpy.std(data[0], axis=0)
    denominator_std_dev = numpy.std(data[1], axis=0)
    XX_std_dev = numpy.sqrt((numerator_std_dev/(normalization_factor*denominator_mean))**2+(numerator_mean*denominator_std_dev/(normalization_factor*denominator_mean**2))**2)
    return XX, XX_std_dev


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
    args = parser.parse_args()

    beta, ns, nt, nt_half = lpd.parse_conftype(args.conftype)

    inputfolder = lpd.get_merged_data_path(args.qcdtype, args.corr, args.conftype)
    outputfolder = inputfolder

    lpd.create_folder(outputfolder+"/btstrp_samples/")
    flow_times, n_flow, n_datafiles, n_streams, XX_data = load_merged_data(args.qcdtype, args.corr, args.conftype)
    XX_data = numpy.asarray(XX_data)
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
    # Use same_rand_for_obs=False in order to break up correlations between the observables
    XX_samples, XX, XX_err = bootstr.bootstr(reduce_function, XX_data, numb_samples=n_samples, sample_size=n_datafiles, conf_axis=1, return_sample=True, same_rand_for_obs=False, parallelize=True, nproc=40, seed=0)

    # extract the right data
    XX_samples_mean = numpy.asarray([i[0] for i in XX_samples])
    XX_samples_err = numpy.asarray([i[1] for i in XX_samples])
    XX_mean = XX[0]  # this is the bootstrap estimate for the mean of the correlator
    XX_err = XX_err[0]  # this is the bootstrap estimate for the error of the mean of the correlator

    # add renormalization factor for BB correlator. errors for it are extremely small. FIXME abstract this
    path = "/work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/coupling/"
    if args.corr == "BB":
        tfT2, Z2 = numpy.loadtxt(path + "Z2_" + str(nt) + ".dat", unpack=True)
        for i in range(XX_mean.shape[0]):
            for j in range(XX_mean.shape[1]):
                XX_mean[i, j] *= Z2[i]
                XX_err[i, j] *= Z2[i]

    # write XX_mean and XX_err to file
    filename = outputfolder+file_prefix+"_"+args.conftype+".dat"
    print("write "+filename)
    with open(filename, 'w') as outfile:
        outfile.write('# bootstrap mean of '+str(n_datafiles)+' measurements of '+args.corr+' correlator for '+args.qcdtype+' '+args.conftype+'\n')
        outfile.write('# rows correspond to flow times, columns to dt = {1, ... , Ntau/2}\n')
        lpd.write_flow_times(outfile, flow_times)
        numpy.savetxt(outfile, XX_mean)

    filename = outputfolder+file_prefix+"_err_"+args.conftype+".dat"
    print("write "+filename)
    with open(filename, 'w') as outfile:
        outfile.write('# bootstrap mean err of '+str(n_datafiles)+' measurements of '+args.corr+' correlator '+args.qcdtype+' '+args.conftype+'\n')
        outfile.write('# rows correspond to flow times, columns to dt = {1, ... , Ntau/2}\n')
        lpd.write_flow_times(outfile, flow_times)
        numpy.savetxt(outfile, XX_err)

    # write bootstrap samples in files for each flowtime separately
    for i in range(n_flow):
        filename = outputfolder+"/btstrp_samples/"+file_prefix+"_"+'{0:.4f}'.format(numpy.sqrt(flow_times[i]*8)/nt)+"_Nt"+str(nt)+"_btstrp_samples.dat"
        print("write "+filename)

        # add BB renorm
        Z = 1
        if args.corr == "BB":
            Z = Z2[i]

        with open(filename, 'w') as outfile:
            outfile.write("# "+str(n_samples)+" bootstrap samples for "+args.corr+" correlator "+args.qcdtype+" "+args.conftype+"\n# one sample consists of "+str(n_datafiles)+" draws(measurements), which are then reduced\n")
            outfile.write("# flowtime="+str(flow_times[i])+" or flowradius="+str(numpy.sqrt(flow_times[i]*8)/nt)+'\n')
            outfile.write("# first column: \\tau_T, second column: G(tauT), third column: std_dev\n")
            for k in range(n_samples):
                XX_sample_single_flowtime = numpy.empty((nt_half, 3))
                for j in range(nt_half):
                    XX_sample_single_flowtime[j, 0] = tauTs[j]
                    XX_sample_single_flowtime[j, 1] = XX_samples_mean[k][i, j] * Z
                    XX_sample_single_flowtime[j, 2] = XX_samples_err[k][i, j] * Z
                numpy.savetxt(outfile, XX_sample_single_flowtime)
                outfile.write('#\n')

    print("done reducing data for "+args.conftype+" "+args.corr)


if __name__ == '__main__':
    main()
    lpd.save_script_call()
