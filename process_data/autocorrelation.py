#!/usr/local/bin/python
import numpy
from os import listdir  
from matplotlib import pyplot as plt
from latqcdtools import autocorrelation


qcdtype = "quenched"
conftype = "s144t36_b0754400"
nt_half=18
flowtype = "zeuthenFlow_acc0.000010_sts0.000010_ColElecCorrTimeSlices" 
outputfolder="../../data_merged/"+qcdtype+"/"+conftype+"/"
inputfolder="../../data_raw/"+qcdtype+"/"+conftype+"/"+conftype+"_a/"
#lpd.create_folder(outputfolder)

single_polyakov, single_EE, EE_numerator_real, EE_numerator_imag, polyakov_real, polyakov_imag, flow_times = ([] for i in range(7))
n_datafiles = 0

files = listdir(inputfolder)

files.sort()

flowindex = 133
tauTindex = 17

"""read in"""
for datafile in files:
    if n_datafiles < 100:
        if datafile.startswith(flowtype): 
            path = inputfolder+"/"+datafile
            tmp = numpy.loadtxt(path)
            if n_datafiles == 0:
                flow_times = tmp[:,0]
            polyakov_real.append(tmp[:,1])
            polyakov_imag.append(tmp[:,2])
            EE_numerator_real.append(tmp[:,3:int((3 + nt_half))])
            EE_numerator_imag.append(tmp[:,int((3 + nt_half)):])
            single_EE.append(EE_numerator_real[n_datafiles][flowindex][tauTindex])
            single_polyakov.append(polyakov_real[n_datafiles][flowindex])
            n_datafiles += 1

t_arr, corr_arr, corr_err_arr, int_corr_arr, int_corr_err_arr = autocorrelation.corr(single_EE)

fig = plt.figure()
ax = fig.add_subplot(1,1,1) 
ax.set_xlabel("autocorrelation time in 500 sweeps")
ax.set_ylabel("autocorrelation of -E-E- at largest separation and used flow time")
plt.errorbar(t_arr, corr_arr)

plt.savefig("autocorrelation_EE.pdf")

t_arr, corr_arr, corr_err_arr, int_corr_arr, int_corr_err_arr = autocorrelation.corr(single_polyakov)

fig = plt.figure()
ax = fig.add_subplot(1,1,1) 
ax.set_xlabel("autocorrelation time in 500 sweeps")
ax.set_ylabel("autocorrelation of polyakovloop at largest used flow time")
plt.errorbar(t_arr, corr_arr)

plt.savefig("autocorrelation_polyakovloop.pdf")
