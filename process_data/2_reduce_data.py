#!/usr/local/bin/python
import numpy
import math 
from os import listdir 
import lib_process_data as lpd

from latqcdtools import jackknife

def normalize_EE( data, data_err, tauT_imp):
    datashape = data.shape
    Ntau = datashape[1]*2
    for i in range(0, datashape[0]):
        for j in range(0, datashape[1]):
            data[i,j] = data[i,j] / lpd.norm_corr(tauT_imp[j]) * Ntau**4
            data_err[i,j] = data_err[i,j] / lpd.norm_corr(tauT_imp[j]) * Ntau**4

qcdtype, conftype, beta, ns, nt, nt_half = lpd.read_args()
 
inputfolder=lpd.inputfolder(qcdtype,conftype)
outputfolder=inputfolder

flow_times, n_flow, n_datafiles, n_streams, EE_data = lpd.load_merged_data(qcdtype, conftype)

"""perform a jackknife analysis of the data where one stream is a block or a bootstrap analysis with 1000 samples. This gives the mean and error of the data"""
if qcdtype == "quenched":
    EE, EE_err = jackknife.jackknife(lpd.compute_EE_mean, EE_data, n_streams, conf_axis=1) 
if qcdtype == "hisq":
    EE, EE_err = bootstr.bootstr(func=lpd.compute_EE_mean, data=EE_data, numb_samples=1000, sample_size=n_datafiles, same_rand_for_obs = True, conf_axis=1, return_sample = False)

normalize_EE(EE, EE_err, lpd.tauT_imp[nt]) 

with open(outputfolder+"EE_"+conftype+".dat", 'w') as outfile:
    outfile.write('# reduced data of '+str(n_datafiles)+' measurements of EE correlator, normalized to G_norm=G_free/(g^2*C_F) with improved distances for '+qcdtype+' '+conftype+'\n')
    outfile.write('# rows correspond to flow times, columns to dt = {1, ... , Ntau/2}\n')
    lpd.write_flow_times(outfile, flow_times)
    numpy.savetxt(outfile, EE)
with open(outputfolder+"EE_err_"+conftype+".dat", 'w') as outfile:
    outfile.write('# jackknife err data of '+str(n_datafiles)+' measurements of EE correlator, normalized to G_norm=G_free/(g^2*C_F) with improved distances for '+qcdtype+' '+conftype+'\n')
    outfile.write('# rows correspond to flow times, columns to dt = {1, ... , Ntau/2}\n')
    lpd.write_flow_times(outfile, flow_times)
    numpy.savetxt(outfile, EE_err)
