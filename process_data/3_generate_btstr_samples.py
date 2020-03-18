#!/usr/local/bin/python
import numpy
import lib_process_data as lpd
from latqcdtools import bootstr
from pathlib import Path


def normalize_EE_samples( data, tauT_imp):
    for k in range(len(data)):
        datashape = data[k].shape
        Ntau = datashape[1]*2
        for i in range(datashape[0]):
            for j in range(datashape[1]):
                data[k][i,j] = data[k][i,j] / lpd.norm_corr(tauT_imp[j]) * Ntau**4

qcdtype, conftype, beta, ns, nt, nt_half = lpd.read_args()
 
inputfolder=lpd.inputfolder(qcdtype,conftype)
outputfolder=inputfolder
Path(outputfolder+"/btstrp_samples/").mkdir(parents=True, exist_ok=True)

flow_times, n_flow, n_datafiles, n_streams, EE_data = lpd.load_merged_data(qcdtype, conftype)

"""generate bootstrap samples, which are used in the cont extr"""
nsamples=1000
EE_samples, EE_bootstrap, EE_err_bootstrap = bootstr.bootstr(func=lpd.compute_EE_mean, data=EE_data, numb_samples=nsamples, sample_size=10000, same_rand_for_obs = True, conf_axis=1, return_sample = True)
normalize_EE_samples(EE_samples, lpd.tauT_imp[nt])

if conftype in ("s064t16_b0687361", "s080t20_b0703500", "s096t24_b0719200", "s120t30_b0739400"):
    for i in range(n_flow):
        EE_sample_single_flowtime = numpy.empty((nt_half*nsamples, 2))
        for k in range(nsamples):
            for j in range(nt_half):
                EE_sample_single_flowtime[k*nt_half+j,0] = lpd.tauT_imp[nt][j]
                EE_sample_single_flowtime[k*nt_half+j,1] = EE_samples[k][i,j]
        numpy.savetxt(outputfolder+"/btstrp_samples/EE_"+'{0:.4f}'.format(numpy.sqrt(flow_times[i]*8)/nt)+"_Nt"+str(nt)+"_btstrp_samples.dat", EE_sample_single_flowtime, header=str(nsamples)+"bootstrap samples for EE correlator "+qcdtype+" "+conftype+"\n # one sample consists of 10000 draws(measurements), which are then reduced")

print("done with", conftype)
