#!/usr/bin/env python3
import numpy
from os import listdir
import lib_process_data as lpd


# TODO write basepath and conf numbers to file for future reference where the data came from


def main():

    # parse cmd line arguments
    parser, requiredNamed = lpd.get_parser()
    parser.add_argument('--acc_sts', help="accuracy and stepsize. format: acc0.000010_sts0.000010")
    requiredNamed.add_argument('--conftype', help="format: s096t20_b0824900 for quenched or s096t20_b0824900_m002022_m01011 for hisq", required=True)
    parser.add_argument('--basepath', help="override default base input path with this one", type=str)
    parser.add_argument('--legacy', help="use legacy file names and legacy multiplicity factor of -3", action="store_true")
    parser.add_argument('--excess_workaround', help="ignore additional flow times at the end of later files", action="store_true")
    parser.add_argument('--reference_flowradii', default=None, type=str,
                        help="only consider flowradii contained in this file. "
                             "the flowradii in this file have to be a subset of the ones you're trying to read in."
                             "useful if you want to combine different runs which different max flow time, or when you made noncritical copy-paste mistake for a few"
                             "flow times.")
    parser.add_argument('--min_conf_nr', help="ignore datafiles with conf number less than this.", default=0, type=int)
    parser.add_argument('--output_basepath', type=str, default="")

    args = parser.parse_args()

    beta, ns, nt, nt_half = lpd.parse_conftype(args.conftype)
    fermions, temp, flowtype, _, _ = lpd.parse_qcdtype(args.qcdtype)

    multiplicity = 1
    if args.corr == "EE":
        if args.legacy:
            XX_label = "ColElecCorrTimeSlices_s"
            multiplicity = -6
        else:
            XX_label = "ColElecCorrTimeSlices_naive_s"
    elif args.corr == "BB":
        XX_label = "ColMagnCorrTimeSlices_naive_s"
    elif args.corr == "EE_clover":
        XX_label = "ColElecCorrTimeSlices_clover_s"
    elif args.corr == "BB_clover":
        XX_label = "ColMagnCorrTimeSlices_clover_s"

    if not args.acc_sts:
        flow_prefix = flowtype+"_"
    else:
        flow_prefix = flowtype+"_"+args.acc_sts + "_"

    if args.basepath:
        inputfolder = args.basepath + "/" + args.qcdtype + "/" + args.conftype + "/"
    else:
        inputfolder = lpd.get_raw_data_path(args.qcdtype, args.conftype)

    outputfolder = lpd.get_merged_data_path(args.qcdtype, args.corr, args.conftype, args.output_basepath)
    lpd.create_folder(outputfolder)

    XX_numerator_real, XX_numerator_imag, polyakov_real, polyakov_imag, flow_times = ([] for _ in range(5))
    n_datafiles, n_streams = (0, 0)
    n_files_per_stream = []
    streamids = []

    full_prefix = flow_prefix+XX_label
    print("searching for files named "+inputfolder+full_prefix+"*")

    # read in data from many files
    shape = (0, 0)
    folders = listdir(inputfolder)
    folders.sort()
    corrupt_files = []
    conf_nums = []
    for stream_folder in folders:
        n_files_this_stream = 0
        if stream_folder.startswith(args.conftype+"_"):
            print(stream_folder, end=', ')
            datafiles = listdir(inputfolder+"/"+stream_folder)
            datafiles.sort()
            for datafile in datafiles:
                if datafile.startswith(full_prefix):
                    path = inputfolder+"/"+stream_folder+"/"+datafile
                    # print("reading "+datafile)
                    confnum = int(lpd.remove_left_of_last('_U', path))
                    if confnum >= args.min_conf_nr:
                        tmp = numpy.loadtxt(path)
                        if shape == (0, 0):
                            shape = tmp.shape
                        if shape != tmp.shape:
                            if args.excess_workaround and tmp.shape != (0,) and tmp.shape[1] == shape[1] and tmp.shape[0] > shape[0]:
                                print("INFO: discarding", tmp.shape[0]-shape[0], " flow times from the end")
                                tmp = tmp[:shape[0]]
                            else:
                                print("error! shapes of input files don't match. ignoring this file.")
                                print(shape, " (previous) vs ", tmp.shape, " (current file)")
                                corrupt_files.append(datafile)
                                continue
                        if n_datafiles == 0:
                            flow_times = tmp[:, 0]
                            if args.reference_flowradii is not None:
                                flowradii = numpy.sqrt(8*flow_times)/nt
                                flowradii_ref = numpy.loadtxt(args.reference_flowradii)
                                indices = []
                                for f in range(len(flow_times)):
                                    for g in range(len(flowradii_ref)):
                                        if numpy.isclose(flowradii[f], flowradii_ref[g]):  # TODO check whether tolerance of isclose is small enough for our stepsizes.
                                            indices.append(f)
                                            break
                                indices = numpy.asarray(indices)
                            else:
                                indices = numpy.asarray(range(0, len(flow_times)))

                            flow_times = flow_times[indices]
                            # print("nflow=", len(flow_times))
                        polyakov_real.append(tmp[indices, 1])
                        polyakov_imag.append(tmp[indices, 2])
                        XX_numerator_real.append(tmp[indices, 3:int((3 + nt_half))] / multiplicity)
                        XX_numerator_imag.append(tmp[indices, int((3 + nt_half)):] / multiplicity)
                        conf_nums.append(confnum)
                        n_datafiles += 1
                        n_files_this_stream += 1
            print("n_files_this_stream: ", n_files_this_stream)
            if n_files_this_stream > 0:
                n_files_per_stream.append(n_files_this_stream)
                n_streams += 1
                streamid = lpd.remove_left_of_last(r'_', stream_folder)
                streamids.append(streamid)
    if len(corrupt_files) > 0:
        print("\n============================= \nWARNING!!!!!!!!!!")
        print("These files are corrupt and were skipped:")
        for file in corrupt_files:
            print(file)
        print("================================\n")
    if n_datafiles == 0:
        print("Didn't find any files! Are the input parameters correct?", args.conftype, beta, ns, nt, nt_half, args.qcdtype, fermions, temp, flowtype, args.corr, args.acc_sts)
        exit()

    # TODO this could be simplified by pickling a dict to disk
    filename = outputfolder+'n_datafiles_'+args.conftype+'.dat'
    print("write "+filename)
    with open(filename, 'w') as outfile:
        outfile.write('# number of datafiles (i.e. confs) for '+args.qcdtype+' '+args.conftype+'\n')
        outfile.write(str(n_datafiles)+'\n')
        outfile.write('# number of streams for '+args.qcdtype+' '+args.conftype+'\n')
        outfile.write(str(n_streams)+'\n')
        # outfile.write('# number of discarded confs from the beginning of each stream:\n')
        # for ndisc, strid in zip(args.n_discard, streamids):
        #     outfile.write(str(ndisc) + '  # ' + strid + '\n')
        numpy.savetxt(outfile, numpy.asarray(n_files_per_stream), header='number of confs contained in each stream respectively', fmt='%i')
        outfile.write('# confnums, ordered by stream\n')
        numpy.savetxt(outfile, numpy.asarray(conf_nums), fmt='%i')
    flow_times = [i for i in flow_times]
    numpy.savetxt(outputfolder+'flowtimes_'+args.conftype+'.dat', flow_times, header=r'flow times \tau_F for '+args.qcdtype+' '+args.conftype)

    # these have the following shape (nconf, nflow, Ntau/2)
    numpy.save(lpd.print_var("write", outputfolder + args.corr + '_real_' + args.conftype + '_merged.npy'), XX_numerator_real)
    numpy.save(lpd.print_var("write", outputfolder + args.corr + '_imag_' + args.conftype + '_merged.npy'), XX_numerator_imag)
    numpy.save(lpd.print_var("write", outputfolder+'polyakov_real_'+args.conftype+'_merged.npy'), polyakov_real)
    numpy.save(lpd.print_var("write", outputfolder+'polyakov_imag_'+args.conftype+'_merged.npy'), polyakov_imag)

    # with open(filename, 'w') as outfile:
    #     outfile.write('# real part of numerator of '+args.corr+' correlator '+args.qcdtype+' '+args.conftype+'\n')
    #     outfile.write('# '+str(n_datafiles)+' confs in one file\n')
    #     outfile.write('# rows correspond to flow times, columns to dt = {1, ... , Ntau/2}\n')
    #     lpd.write_flow_times(outfile, flow_times)
    #     for conf in XX_numerator_real:
    #         numpy.savetxt(outfile, conf)
    #         outfile.write('# \n')


    # with open(filename, 'w') as outfile:
    #     outfile.write('# imag part of numerator of '+args.corr+' correlator '+args.qcdtype+' '+args.conftype+'\n')
    #     outfile.write('# '+str(n_datafiles)+' confs in one file\n')
    #     outfile.write('# rows correspond to flow times, columns to dt = {1, ... , Ntau/2}\n')
    #     lpd.write_flow_times(outfile, flow_times)
    #     for conf in XX_numerator_imag:
    #         numpy.savetxt(outfile, conf)
    #         outfile.write('# \n')

    # filename =
    # print("write "+filename)
    # with open(filename, 'w') as outfile:
    #     outfile.write('# real part of polyakov loop '+args.qcdtype+' '+args.conftype+'\n')
    #     outfile.write('# '+str(n_datafiles)+' confs in one file\n')
    #     outfile.write('# rows correspond to flow times\n')
    #     lpd.write_flow_times(outfile, flow_times)
    #     for conf in polyakov_real:
    #         numpy.savetxt(outfile, conf)
    #         outfile.write('# \n')
    #
    # filename = outputfolder+'polyakov_imag_'+args.conftype+'_merged.dat'
    # print("write "+filename)
    # with open(filename, 'w') as outfile:
    #     outfile.write('# imag part of polyakov loop '+args.qcdtype+' '+args.conftype+'\n')
    #     outfile.write('# '+str(n_datafiles)+' confs in one file\n')
    #     outfile.write('# rows correspond to flow times\n')
    #     lpd.write_flow_times(outfile, flow_times)
    #     for conf in polyakov_imag:
    #         numpy.savetxt(outfile, conf)
    #         outfile.write('# \n')

    print("done with "+args.qcdtype+" "+args.conftype)


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
