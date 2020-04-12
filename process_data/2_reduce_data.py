#!/usr/local/bin/python
import numpy
import math 
from os import listdir 
import lib_process_data as lpd

from latqcdtools import jackknife
from latqcdtools import bootstr

def normalize_EE( data, data_err, tauT_imp):
    datashape = data.shape
    Ntau = datashape[1]*2
    for i in range(0, datashape[0]):
        for j in range(0, datashape[1]):
            data[i,j] = data[i,j] / lpd.norm_corr(tauT_imp[j]) * Ntau**4
            data_err[i,j] = data_err[i,j] / lpd.norm_corr(tauT_imp[j]) * Ntau**4


def compute_means_separate(data): #mean function for Haukes jackknife routine
    mean = []
    for observable in data:
        mean.append(numpy.mean(observable, axis=0))
    return mean

#numerator data is first index, polyakov loop is second index of data
def compute_EE_correlator(data):
    EE_correlator = numpy.mean(data[0], axis=0) / (-6* numpy.mean(data[1], axis=0))
    Ntau = EE_correlator.shape[1]*2
    tauT_imp = lpd.tauT_imp[Ntau]
    for i in range(len(EE_correlator)):
        for j in range(len(EE_correlator[0])):
            EE_correlator[i,j] = EE_correlator[i,j] / lpd.norm_corr(tauT_imp[j]) * Ntau**4
    #print(EE_correlator)
    return EE_correlator

def compute_EE(data):
    EE_correlator = numpy.mean(data[0], axis=0) / (-6* numpy.mean(data[1], axis=0))
    Ntau = EE_correlator.shape[1]*2
    tauT_imp = lpd.tauT_imp[Ntau]
    for i in range(len(EE_correlator)):
        for j in range(len(EE_correlator[0])):
            EE_correlator[i,j] = EE_correlator[i,j] / lpd.norm_corr(tauT_imp[j]) * Ntau**4
    #print(EE_correlator)
    return EE_correlator

qcdtype, conftype, beta, ns, nt, nt_half = lpd.read_args()
 
inputfolder=lpd.inputfolder(qcdtype,conftype)
outputfolder=inputfolder

flow_times, n_flow, n_datafiles, n_streams, EE_data = lpd.load_merged_data(qcdtype, conftype)

nsamples = 10000

"""perform a jackknife analysis of the data where one stream is a block or a bootstrap analysis with 1000 samples. This gives the mean and error of the data"""
if qcdtype == "quenched":
    #generate bootstrap samples with same random numbers for the numerator and polyakov loop
    #EE_samples_separate, dummy, dummy = bootstr.bootstr(compute_means_separate, EE_data, numb_samples=nsamples, sample_size=n_datafiles, conf_axis=1, return_sample=True, same_rand_for_obs=True)
    #rearrange the samples for a second bootstrap
    #EE_numerator_samples = [k[0] for k in EE_samples_separate]
    #polyakov_samples = [k[1] for k in EE_samples_separate]
    #EE_samples_separate = [EE_numerator_samples, polyakov_samples]
    #do the second bootstrap based on the samples on the first one to get samples for the correlator and mean+error estimates
    #EE_samples_combined, EE, EE_err = bootstr.bootstr(compute_EE_correlator, EE_samples_separate, numb_samples=nsamples, sample_size=nsamples, conf_axis=1, return_sample=True, same_rand_for_obs=True)
    EE_samples, EE, EE_err = bootstr.bootstr(compute_EE, EE_data, numb_samples=nsamples, sample_size=n_datafiles, conf_axis=1, return_sample=True, same_rand_for_obs=True)
    
    #EE, EE_err = jackknife.jackknife(lpd.compute_EE_mean, EE_data, n_streams, conf_axis=1) 
    #hier bootstrap für die beiden größen einzeln performen (dafür jeweils eine funktion schreiben), samples zurückgeben lassen. aus den samples dann ein neuen datenblock machen und daraus wieder bootstrappen und samples zurückgeben (auch hierfür wird eine funktion benötigt)
    
    
#if qcdtype == "hisq":
    #EE, EE_err = bootstr.bootstr(func=lpd.compute_EE_mean, data=EE_data, numb_samples=1000, sample_size=n_datafiles, same_rand_for_obs = True, conf_axis=1, return_sample = False)

#normalize_EE(EE, EE_err, lpd.tauT_imp[nt]) 

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


for i in range(n_flow):
    with open(outputfolder+"/btstrp_samples/EE_"+'{0:.4f}'.format(numpy.sqrt(flow_times[i]*8)/nt)+"_Nt"+str(nt)+"_btstrp_samples.dat", 'w') as outfile:
        outfile.write("# "+str(nsamples)+" bootstrap samples for EE correlator "+qcdtype+" "+conftype+"\n# one sample consists of "+str(n_datafiles)+" draws(measurements), which are then reduced\n")
        outfile.write("# flowtime="+str(flow_times[i])+" or flowradius="+str(numpy.sqrt(flow_times[i]*8)/nt)+'\n')
        outfile.write("# first column: \\tau_T, second column: G/G_norm\n")
        for k in range(nsamples):
            EE_sample_single_flowtime = numpy.empty((nt_half, 2))
            for j in range(nt_half):
                EE_sample_single_flowtime[j,0] = lpd.tauT_imp[nt][j]
                EE_sample_single_flowtime[j,1] = EE_samples[k][i,j]
            numpy.savetxt(outfile, EE_sample_single_flowtime)
            outfile.write('#\n')

print("done reducing data for "+conftype)
