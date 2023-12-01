#!/usr/bin/env python3
import numpy
from os import listdir
import lib_process_data as lpd
from natsort import natsorted

# TODO write basepath and conf numbers to file for future reference where the data came from


def get_flow_indices(flowradii, flowradii_ref):
    indices = []
    for f in range(len(flowradii)):
        for g in range(len(flowradii_ref)):
            if numpy.isclose(flowradii[f], flowradii_ref[g]):  # TODO check whether tolerance of isclose is small enough for our stepsizes.
                indices.append(f)
                break
    return numpy.asarray(indices)


def are_all_floats_of_a_in_b(a, b):
    return numpy.array([numpy.isclose(val, b).any() for val in a]).all()


def main():

    # it is assumed that the first data file that is read in is not corrupted and has the correct shape.

    # parse cmd line arguments
    parser, requiredNamed = lpd.get_parser()
    parser.add_argument('--acc_sts', help="accuracy and stepsize. format: acc0.000010_sts0.000010")
    requiredNamed.add_argument('--conftype', help="format: s096t20_b0824900 for quenched or s096t20_b0824900_m002022_m01011 for hisq", required=True)
    parser.add_argument('--basepath', help="override default base input path with this one", type=str)
    parser.add_argument('--legacy', help="use legacy file names and legacy multiplicity factor", action="store_true")
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
    print("searching in subfolders of "+inputfolder + "   for  ", full_prefix+"*")

    flowradii_ref = None
    if args.reference_flowradii is not None:
        flowradii_ref = numpy.loadtxt(args.reference_flowradii)

    # read in data from many files
    shape = (0, 0)
    folders = listdir(inputfolder)
    folders = natsorted(folders)
    corrupt_files = []
    conf_nums = []
    for stream_folder in folders:
        n_files_this_stream = 0
        if stream_folder.startswith(args.conftype+"_"):
            print(stream_folder, end=', ')
            datafiles = listdir(inputfolder+"/"+stream_folder)
            datafiles = natsorted(datafiles)
            for datafile in datafiles:
                if datafile.startswith(full_prefix):
                    path = inputfolder+"/"+stream_folder+"/"+datafile
                    # print("reading "+datafile)
                    confnum = int(lpd.remove_left_of_last('_U', path))
                    if confnum >= args.min_conf_nr:
                        tmp = numpy.loadtxt(path)
                        if tmp.shape != (0,):
                            # print(tmp.shape)
                            if shape == (0, 0):
                                shape = tmp.shape
                            shape_wrong_but_all_flowtimes_exist = False
                            if shape != tmp.shape:
                                flow_times = tmp[:, 0]
                                these_flowradii = numpy.sqrt(8 * flow_times) / nt
                                if args.excess_workaround and tmp.shape != (0,) and tmp.shape[1] == shape[1] and tmp.shape[0] > shape[0]:
                                    print("INFO: discarding", tmp.shape[0]-shape[0], " flow times from the end")
                                    tmp = tmp[:shape[0]]
                                elif args.reference_flowradii and are_all_floats_of_a_in_b(these_flowradii, flowradii_ref):
                                    shape_wrong_but_all_flowtimes_exist = True
                                    # TODO check that all flowradii in flowradii are also in the ref file.
                                else:
                                    print("error! shapes of input files don't match. ignoring this file.")
                                    print(shape, " (previous) vs ", tmp.shape, " (current file)")
                                    corrupt_files.append(datafile)
                                    continue

                            flow_times = tmp[:, 0]
                            if args.reference_flowradii is not None or shape_wrong_but_all_flowtimes_exist:
                                flowradii = numpy.sqrt(8*flow_times)/nt
                                indices = get_flow_indices(flowradii, flowradii_ref)
                                if len(indices) != len(flowradii_ref):
                                    print("error! could not find all ref flowradii. skipping", datafile)
                                    continue
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

    filename = outputfolder+'n_datafiles_'+args.conftype+'.dat'
    print("write "+filename)
    with open(filename, 'w') as outfile:
        outfile.write('# number of datafiles (i.e. confs) for '+args.qcdtype+' '+args.conftype+'\n')
        outfile.write(str(n_datafiles)+'\n')
        outfile.write('# number of streams for '+args.qcdtype+' '+args.conftype+'\n')
        outfile.write(str(n_streams)+'\n')
        numpy.savetxt(outfile, numpy.asarray(n_files_per_stream), header='number of confs contained in each stream respectively', fmt='%i')
        outfile.write('# confnums, ordered by stream\n')
        numpy.savetxt(outfile, numpy.asarray(conf_nums), fmt='%i')
    flowtime_filename = outputfolder+'flowtimes_'+args.conftype+'.dat'
    print("write ", flowtime_filename)
    numpy.savetxt(flowtime_filename, flow_times, header=r'flow times \tau_F for '+args.qcdtype+' '+args.conftype)

    # these have the following shape (nconf, nflow, Ntau/2)
    numpy.save(lpd.print_var("write", outputfolder + args.corr + '_real_' + args.conftype + '_merged.npy'), XX_numerator_real)
    numpy.save(lpd.print_var("write", outputfolder + args.corr + '_imag_' + args.conftype + '_merged.npy'), XX_numerator_imag)
    numpy.save(lpd.print_var("write", outputfolder+'polyakov_real_'+args.conftype+'_merged.npy'), polyakov_real)
    numpy.save(lpd.print_var("write", outputfolder+'polyakov_imag_'+args.conftype+'_merged.npy'), polyakov_imag)

    print("done with "+args.qcdtype+" "+args.conftype)


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
