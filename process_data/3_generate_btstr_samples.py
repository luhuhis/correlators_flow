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

"""own custom bootstrap routine"""
#def bootstrap(func, data, n_samples = 10, sample_size = 10000, conf_axis = 1)
    
    #has to return the individual samples with their respective covariance matrices



"""use analysistoolbox"""
"""generate bootstrap samples, which are used in the cont extr"""
nsamples=10
EE_samples, EE_bootstrap, EE_err_bootstrap = bootstr.bootstr(func=lpd.compute_EE_mean, data=EE_data, numb_samples=nsamples, sample_size=10000, same_rand_for_obs = True, conf_axis=1, return_sample = True)

"""calculate covariance matrix"""
#ich habe alle einzelnen samples und ich habe den mittelwert und die standardabweichung
#cov(X,Y)=1/n sum_(i=1)^n (x_i - E(X))(y_i - E(Y))


normalize_EE_samples(EE_samples, lpd.tauT_imp[nt])

if conftype in ("s064t16_b0687361", "s080t20_b0703500", "s096t24_b0719200", "s120t30_b0739400"):
    for i in range(n_flow):
        with open(outputfolder+"/btstrp_samples/EE_"+'{0:.4f}'.format(numpy.sqrt(flow_times[i]*8)/nt)+"_Nt"+str(nt)+"_btstrp_samples.dat", 'w') as outfile:
            outfile.write("# "+str(nsamples)+" bootstrap samples for EE correlator "+qcdtype+" "+conftype+"\n# one sample consists of 10000 draws(measurements), which are then reduced\n")
            outfile.write("# flowtime="+str(flow_times[i])+" or flowradius="+str(numpy.sqrt(flow_times[i]*8)/nt)+'\n')
            outfile.write("# first column: \\tau_T, second column: G/G_norm\n")
            for k in range(nsamples):
                EE_sample_single_flowtime = numpy.empty((nt_half, 2))
                for j in range(nt_half):
                    EE_sample_single_flowtime[j,0] = lpd.tauT_imp[nt][j]
                    EE_sample_single_flowtime[j,1] = EE_samples[k][i,j]
                numpy.savetxt(outfile, EE_sample_single_flowtime)
                outfile.write('#\n')

print("done with", conftype)
