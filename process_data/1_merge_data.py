#!/usr/local/bin/python
import numpy
from os import listdir 
import lib_process_data as lpd

qcdtype, conftype, beta, ns, nt, nt_half = lpd.read_args()

outputfolder="../../data_merged/"+qcdtype+"/"+conftype+"/"
inputfolder="../../data_raw/"+qcdtype+"/"+conftype+"/"


EE_numerator_real, EE_numerator_imag, polyakov_real, polyakov_imag, flow_times = ([] for i in range(5))
n_datafiles, n_streams = (0,0)

"""read in"""
for stream_folder in listdir(inputfolder):
    if stream_folder.startswith(conftype+"_"):
        n_streams += 1
        for datafile in listdir(inputfolder+"/"+stream_folder):
            if datafile.startswith("zeuthenFlow"): 
                path = inputfolder+"/"+stream_folder+"/"+datafile
                tmp = numpy.loadtxt(path) 
                if n_datafiles == 0:
                    flow_times = tmp[:,0]
                polyakov_real.append(tmp[:,1])
                polyakov_imag.append(tmp[:,2])
                EE_numerator_real.append(tmp[:,3:int((3 + nt_half))])
                EE_numerator_imag.append(tmp[:,int((3 + nt_half)):])
                n_datafiles += 1

"""write to files"""
with open(outputfolder+'n_datafiles_'+conftype+'.dat', 'w') as outfile:
    outfile.write('# number of datafiles (i.e. confs) for '+qcdtype+' '+conftype+'\n')
    outfile.write(str(n_datafiles)+'\n')
    outfile.write('# number of streams for '+qcdtype+' '+conftype+'\n')
    outfile.write(str(n_streams)+'\n')
flow_times=[i for i in flow_times]
flow_radii = [numpy.sqrt(i*8)/nt for i in flow_times]
numpy.savetxt(outputfolder+'flowtimes_'+conftype+'.dat', flow_times, header='flow times \\tau_F for '+qcdtype+' '+conftype)
numpy.savetxt(outputfolder+'flowradii_'+conftype+'.dat', flow_radii, header='flow radii (\sqrt{8\\tau_F}T) for '+qcdtype+' '+conftype)
with open(outputfolder+'EE_real_'+conftype+'_merged.dat', 'w') as outfile:
    outfile.write('# real part of numerator of EE correlator '+qcdtype+' '+conftype+'\n')
    outfile.write('# '+str(n_datafiles)+' confs in one file\n')
    outfile.write('# rows correspond to flow times, columns to dt = {1, ... , Ntau/2}\n')
    lpd.write_flow_times(outfile, flow_times)
    for conf in EE_numerator_real: 
        numpy.savetxt(outfile, conf)
        outfile.write('# \n')
with open(outputfolder+'EE_imag_'+conftype+'_merged.dat', 'w') as outfile:
    outfile.write('# imag part of numerator of EE correlator '+qcdtype+' '+conftype+'\n')
    outfile.write('# '+str(n_datafiles)+' confs in one file\n')
    outfile.write('# rows correspond to flow times, columns to dt = {1, ... , Ntau/2}\n')
    lpd.write_flow_times(outfile, flow_times)
    for conf in EE_numerator_imag: 
        numpy.savetxt(outfile, conf)
        outfile.write('# \n')
with open(outputfolder+'polyakov_real_'+conftype+'_merged.dat', 'w') as outfile:
    outfile.write('# real part of polyakov loop '+qcdtype+' '+conftype+'\n')
    outfile.write('# '+str(n_datafiles)+' confs in one file\n')
    outfile.write('# rows correspond to flow times\n')
    lpd.write_flow_times(outfile, flow_times)
    for conf in polyakov_real: 
        numpy.savetxt(outfile, conf)
        outfile.write('# \n')
with open(outputfolder+'polyakov_imag_'+conftype+'_merged.dat', 'w') as outfile:
    outfile.write('# imag part of polyakov loop '+qcdtype+' '+conftype+'\n')
    outfile.write('# '+str(n_datafiles)+' confs in one file\n')
    outfile.write('# rows correspond to flow times\n')
    lpd.write_flow_times(outfile, flow_times)
    for conf in polyakov_imag: 
        numpy.savetxt(outfile, conf)
        outfile.write('# \n')
        
print("done with "+conftype)
