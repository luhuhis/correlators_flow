#!/usr/local/bin/python
import numpy
from os import listdir 
import lib_process_data as lpd


def main():

    conftype, beta, ns, nt, nt_half, qcdtype, fermions, temp, flowtype, corr, add_params = lpd.read_args()

    accuracy_and_stepsize = "acc0.000010_sts0.000010_"
    if add_params is not None:
        accuracy_and_stepsize = add_params[0]

    if corr == "EE": 
        XX_label = "ColElecCorrTimeSlices_s"
    elif corr == "BB":
        XX_label = "ColMagCorrTimeSlices_s"
    elif corr == "EE_clover":
        XX_label = "ColElecCorrTimeSlices_clover_s"
    elif corr == "BB_clover":
        XX_label = "ColMagCorrTimeSlices_clover_s"

    flow_prefix = flowtype+"_"+accuracy_and_stepsize
    
    inputfolder = lpd.get_raw_data_path(qcdtype, conftype)
    outputfolder = lpd.get_merged_data_path(qcdtype, corr, conftype)
    lpd.create_folder(outputfolder)

    XX_numerator_real, XX_numerator_imag, polyakov_real, polyakov_imag, flow_times = ([] for i in range(5))
    n_datafiles, n_streams = (0,0)
    
    full_prefix = flow_prefix+XX_label
    print("searching for files named "+inputfolder+full_prefix+"*")

    """read in data from many files"""
    for stream_folder in listdir(inputfolder):
        print(stream_folder)
        if stream_folder.startswith(conftype+"_"):
            n_streams += 1
            datafiles = listdir(inputfolder+"/"+stream_folder)
            datafiles.sort()
            for datafile in datafiles:
                if datafile.startswith(full_prefix): 
                    path = inputfolder+"/"+stream_folder+"/"+datafile
                    tmp = numpy.loadtxt(path) 
                    if n_datafiles == 0:
                        flow_times = tmp[:,0]
                    polyakov_real.append(tmp[:,1])
                    polyakov_imag.append(tmp[:,2])
                    XX_numerator_real.append(tmp[:,3:int((3 + nt_half))])
                    XX_numerator_imag.append(tmp[:,int((3 + nt_half)):])
                    print("read "+datafile)
                    n_datafiles += 1                
    if n_datafiles == 0:
        print("Didn't find any files! Are the input parameters correct?", conftype, beta, ns, nt, nt_half, qcdtype, fermions, temp, flowtype, corr, accuracy_and_stepsize)
        exit()

    """write data to one file"""
    print("write "+outputfolder+'n_datafiles_'+conftype+'.dat')
    with open(outputfolder+'n_datafiles_'+conftype+'.dat', 'w') as outfile:
        outfile.write('# number of datafiles (i.e. confs) for '+qcdtype+' '+conftype+'\n')
        outfile.write(str(n_datafiles)+'\n')
        outfile.write('# number of streams for '+qcdtype+' '+conftype+'\n')
        outfile.write(str(n_streams)+'\n')
    flow_times=[i for i in flow_times]
    flow_radii = [numpy.sqrt(i*8)/nt for i in flow_times]
    numpy.savetxt(outputfolder+'flowtimes_'+conftype+'.dat', flow_times, header='flow times \\tau_F for '+qcdtype+' '+conftype)
    numpy.savetxt(outputfolder+'flowradii_'+conftype+'.dat', flow_radii, header='flow radii (\sqrt{8\\tau_F}T) for '+qcdtype+' '+conftype)
    print("write "+outputfolder+corr+'_real_'+conftype+'_merged.dat')
    with open(outputfolder+corr+"_real_"+conftype+'_merged.dat', 'w') as outfile:
        outfile.write('# real part of numerator of '+corr+' correlator '+qcdtype+' '+conftype+'\n')
        outfile.write('# '+str(n_datafiles)+' confs in one file\n')
        outfile.write('# rows correspond to flow times, columns to dt = {1, ... , Ntau/2}\n')
        lpd.write_flow_times(outfile, flow_times)
        for conf in XX_numerator_real: 
            numpy.savetxt(outfile, conf)
            outfile.write('# \n')
    print("write "+outputfolder+corr+'_imag_'+conftype+'_merged.dat')
    with open(outputfolder+corr+'_imag_'+conftype+'_merged.dat', 'w') as outfile:
        outfile.write('# imag part of numerator of '+corr+' correlator '+qcdtype+' '+conftype+'\n')
        outfile.write('# '+str(n_datafiles)+' confs in one file\n')
        outfile.write('# rows correspond to flow times, columns to dt = {1, ... , Ntau/2}\n')
        lpd.write_flow_times(outfile, flow_times)
        for conf in XX_numerator_imag: 
            numpy.savetxt(outfile, conf)
            outfile.write('# \n')
    print("write "+outputfolder+'polyakov_real_'+conftype+'_merged.dat')
    with open(outputfolder+'polyakov_real_'+conftype+'_merged.dat', 'w') as outfile:
        outfile.write('# real part of polyakov loop '+qcdtype+' '+conftype+'\n')
        outfile.write('# '+str(n_datafiles)+' confs in one file\n')
        outfile.write('# rows correspond to flow times\n')
        lpd.write_flow_times(outfile, flow_times)
        for conf in polyakov_real: 
            numpy.savetxt(outfile, conf)
            outfile.write('# \n')
    print("write "+outputfolder+'polyakov_imag_'+conftype+'_merged.dat')
    with open(outputfolder+'polyakov_imag_'+conftype+'_merged.dat', 'w') as outfile:
        outfile.write('# imag part of polyakov loop '+qcdtype+' '+conftype+'\n')
        outfile.write('# '+str(n_datafiles)+' confs in one file\n')
        outfile.write('# rows correspond to flow times\n')
        lpd.write_flow_times(outfile, flow_times)
        for conf in polyakov_imag: 
            numpy.savetxt(outfile, conf)
            outfile.write('# \n')
            
    print("done with "+qcdtype+" "+conftype)
    
if __name__ == '__main__':
    main()
