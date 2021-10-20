#!/usr/local/bin/python
import numpy
from os import listdir 
import lib_process_data as lpd


def main():

  
    #parse cmd line arguments
    parser, requiredNamed = lpd.get_parser()
    requiredNamed.add_argument('--acc_sts', help="accuracy and stepsize. format: acc0.000010_sts0.000010", default="acc0.000010_sts0.000010", required=True)
    requiredNamed.add_argument('--conftype', help="format: s096t20_b0824900 for quenched or s096t20_b0824900_m002022_m01011 for hisq", required=True)
    
    
    
    args = parser.parse_args()
   
    beta, ns, nt, nt_half = lpd.parse_conftype(args.conftype)
    fermions, temp, flowtype = lpd.parse_qcdtype(args.qcdtype)
   
   
    if args.corr == "EE": 
        XX_label = "ColElecCorrTimeSlices_naive_s"
    elif args.corr == "BB":
        XX_label = "ColMagnCorrTimeSlices_naive_s"
    elif args.corr == "EE_clover":
        XX_label = "ColElecCorrTimeSlices_clover_s"
    elif args.corr == "BB_clover":
        XX_label = "ColMagnCorrTimeSlices_clover_s"

    flow_prefix = flowtype+"_"+args.acc_sts+"_"
    
    inputfolder = lpd.get_raw_data_path(args.qcdtype, args.conftype)
    outputfolder = lpd.get_merged_data_path(args.qcdtype, args.corr, args.conftype)
    lpd.create_folder(outputfolder)

    XX_numerator_real, XX_numerator_imag, polyakov_real, polyakov_imag, flow_times = ([] for i in range(5))
    n_datafiles, n_streams = (0,0)
    
    full_prefix = flow_prefix+XX_label
    print("searching for files named "+inputfolder+full_prefix+"*")

    """read in data from many files"""
    shape = (0,0)
    for stream_folder in listdir(inputfolder):
        print(stream_folder)
        if stream_folder.startswith(args.conftype+"_"):
            n_streams += 1
            datafiles = listdir(inputfolder+"/"+stream_folder)
            datafiles.sort()
            for datafile in datafiles:
                if datafile.startswith(full_prefix): 
                    path = inputfolder+"/"+stream_folder+"/"+datafile
                    print("reading "+datafile)
                    tmp = numpy.loadtxt(path) 
                    if shape == (0,0):
                        shape = tmp.shape
                    if shape != tmp.shape:
                        print("error! shapes of input files don't match")
                        print(shape, " (previous) vs ", tmp.shape, " (current file)")
                        exit()
                    if n_datafiles == 0:
                        flow_times = tmp[:,0]
                    polyakov_real.append(tmp[:,1])
                    polyakov_imag.append(tmp[:,2])
                    XX_numerator_real.append(tmp[:,3:int((3 + nt_half))])
                    XX_numerator_imag.append(tmp[:,int((3 + nt_half)):])
                    n_datafiles += 1                
    if n_datafiles == 0:
        print("Didn't find any files! Are the input parameters correct?", args.conftype, beta, ns, nt, nt_half, args.qcdtype, fermions, temp, flowtype, args.corr, args.acc_sts)
        exit()

    """write data to one file"""
    filename = outputfolder+'n_datafiles_'+args.conftype+'.dat'
    print("write "+filename)
    with open(filename, 'w') as outfile:
        outfile.write('# number of datafiles (i.e. confs) for '+args.qcdtype+' '+args.conftype+'\n')
        outfile.write(str(n_datafiles)+'\n')
        outfile.write('# number of streams for '+args.qcdtype+' '+args.conftype+'\n')
        outfile.write(str(n_streams)+'\n')
    flow_times=[i for i in flow_times]
    flow_radii = [numpy.sqrt(i*8)/nt for i in flow_times]
    numpy.savetxt(outputfolder+'flowtimes_'+args.conftype+'.dat', flow_times, header='flow times \\tau_F for '+args.qcdtype+' '+args.conftype)
    numpy.savetxt(outputfolder+'flowradii_'+args.conftype+'.dat', flow_radii, header='flow radii (\sqrt{8\\tau_F}T) for '+args.qcdtype+' '+args.conftype)
    filename = outputfolder+args.corr+'_real_'+args.conftype+'_merged.dat'
    print("write "+filename)
    with open(filename, 'w') as outfile:
        outfile.write('# real part of numerator of '+args.corr+' correlator '+args.qcdtype+' '+args.conftype+'\n')
        outfile.write('# '+str(n_datafiles)+' confs in one file\n')
        outfile.write('# rows correspond to flow times, columns to dt = {1, ... , Ntau/2}\n')
        lpd.write_flow_times(outfile, flow_times)
        for conf in XX_numerator_real: 
            numpy.savetxt(outfile, conf)
            outfile.write('# \n')
    filename=outputfolder+args.corr+'_imag_'+args.conftype+'_merged.dat'
    print("write "+filename)
    with open(filename, 'w') as outfile:
        outfile.write('# imag part of numerator of '+args.corr+' correlator '+args.qcdtype+' '+args.conftype+'\n')
        outfile.write('# '+str(n_datafiles)+' confs in one file\n')
        outfile.write('# rows correspond to flow times, columns to dt = {1, ... , Ntau/2}\n')
        lpd.write_flow_times(outfile, flow_times)
        for conf in XX_numerator_imag: 
            numpy.savetxt(outfile, conf)
            outfile.write('# \n')
    filename=outputfolder+'polyakov_real_'+args.conftype+'_merged.dat'
    print("write "+filename)
    with open(filename, 'w') as outfile:
        outfile.write('# real part of polyakov loop '+args.qcdtype+' '+args.conftype+'\n')
        outfile.write('# '+str(n_datafiles)+' confs in one file\n')
        outfile.write('# rows correspond to flow times\n')
        lpd.write_flow_times(outfile, flow_times)
        for conf in polyakov_real: 
            numpy.savetxt(outfile, conf)
            outfile.write('# \n')
    filename=outputfolder+'polyakov_imag_'+args.conftype+'_merged.dat'
    print("write "+filename)
    with open(filename, 'w') as outfile:
        outfile.write('# imag part of polyakov loop '+args.qcdtype+' '+args.conftype+'\n')
        outfile.write('# '+str(n_datafiles)+' confs in one file\n')
        outfile.write('# rows correspond to flow times\n')
        lpd.write_flow_times(outfile, flow_times)
        for conf in polyakov_imag: 
            numpy.savetxt(outfile, conf)
            outfile.write('# \n')
            
    print("done with "+args.qcdtype+" "+args.conftype)
    
if __name__ == '__main__':
    main()
    lpd.save_script_call()
