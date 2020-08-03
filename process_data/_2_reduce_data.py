#!/usr/local/bin/python
import numpy
import math 
from os import listdir 
import lib_process_data as lpd

#from latqcdtools import jackknife
from latqcdtools import bootstr
from latqcdtools import logger

#function for bootstrap routine that computes the EE correlator. numerator data is first index, polyakov loop is second index of data.
def compute_EE(data):
    numerator_mean = numpy.mean(data[0], axis=0)
    denominator_mean = numpy.mean(data[1], axis=0)
    EE = numerator_mean / (-6* denominator_mean)
    numerator_std_dev = numpy.std(data[0], axis=0)
    denominator_std_dev = numpy.std(data[1], axis=0)
    EE_err = numpy.sqrt((numerator_std_dev/(6*denominator_mean))**2+(numerator_mean*denominator_std_dev/(6*denominator_mean**2))**2)
    Ntau_half = EE.shape[1]
    Ntau = Ntau_half*2
    nflow = EE.shape[0]
    tauT_imp = lpd.tauT_imp[Ntau]
    #compute time slice covariance matrix for every flowtime
    #nt_cov = numpy.zeros((nflow,Ntau_half,Ntau_half))
    #for i in range(nflow):
        #for j in range(Ntau_half):
            #for k in range(Ntau_half):
                #for l in range(len(data)):
                    #EEj = data[0,l,i,j] / (-6 * data[1,l,i,j]) 
                    #EEk = data[0,l,i,k] / (-6 * data[1,l,i,k]) 
                    #nt_cov[i,j,k] += 1/len(data) * (EEj - EE[i,j]) * (EEk - EE[i,k])
    return EE, EE_err


#def calc_mean(data):
    #mean = numpy.mean(data, axis=0)
    #return mean

qcdtype, conftype, beta, ns, nt, nt_half = lpd.read_args() 
inputfolder=lpd.inputfolder(qcdtype,conftype)
outputfolder=inputfolder
lpd.create_folder(outputfolder+"/btstrp_samples/")
flow_times, n_flow, n_datafiles, n_streams, EE_data = lpd.load_merged_data(qcdtype, conftype)

#do separate bootstrap because only then we can get the covariance matrices for a single sample

n_samples = 10000

EE_data = numpy.asarray(EE_data)


#numerator_samples, numerator, numerator_err = bootstr.bootstr(calc_mean, EE_data[0], numb_samples=n_samples_1, sample_size=n_datafiles, conf_axis=0, return_sample=True, same_rand_for_obs=True)
#polyakov_samples, polyakov, polyakov_err = bootstr.bootstr(calc_mean, EE_data[1], numb_samples=n_samples_1, sample_size=n_datafiles, conf_axis=0, return_sample=True, same_rand_for_obs=True)

#combined_samples = numpy.asarray([numerator_samples, polyakov_samples])

EE_samples, EE, EE_err = bootstr.bootstr(compute_EE, EE_data, numb_samples=n_samples, sample_size=n_datafiles, conf_axis=1, return_sample=True, same_rand_for_obs=False)

#from the samples and the mean now compute the covariance matrix

#extract the right data
EE_samples_mean = numpy.asarray([i[0] for i in EE_samples])
EE_samples_err = numpy.asarray([i[1] for i in EE_samples])
#EE_samples_cov = numpy.asarray([i[1] for i in EE_samples])
EE_mean = EE[0] #this is the bootstrap estimate for the average of the correlator 
EE_err = EE_err[0] #this is the bootstrap estimate for the error of the average of the correlator

#EE_cov = numpy.zeros((n_flow,nt_half,nt_half)) #calculate the covariance matrix of the average of the correlator from the samples
#for i in range(n_flow):
        #for j in range(nt_half):
            #for k in range(nt_half):
                #for l in range(n_samples):
                    #EE_cov[i,j,k] += 1/(n_samples-1) * (EE_samples_mean[l,i,j] - EE_mean[i,j]) * (EE_samples_mean[l,i,k] - EE_mean[i,k])



#write EE_mean and EE_err to file
with open(outputfolder+"EE_"+conftype+".dat", 'w') as outfile:
    outfile.write('# reduced data of '+str(n_datafiles)+' measurements of EE correlator, normalized to G_norm=G_free/(g^2*C_F) with improved distances for '+qcdtype+' '+conftype+'\n')
    outfile.write('# rows correspond to flow times, columns to dt = {1, ... , Ntau/2}\n')
    lpd.write_flow_times(outfile, flow_times)
    numpy.savetxt(outfile, EE_mean)
#with open(outputfolder+"EE_cov_"+conftype+".dat", 'w') as outfile:
    #outfile.write('# reduced covariance matrix of '+str(n_datafiles)+' measurements of EE correlator, normalized to G_norm=G_free/(g^2*C_F) with improved distances for '+qcdtype+' '+conftype+'\n')
    #outfile.write('# covariance matrix of the time slices for each flow time, rows and columns correspond to dt = {1, ... , Ntau/2}\n')
    #lpd.write_flow_times(outfile, flow_times)
    #for i in range(n_flow):
        #numpy.savetxt(outfile, EE_cov[i])
        #outfile.write('#\n')
with open(outputfolder+"EE_err_"+conftype+".dat", 'w') as outfile:
    outfile.write('# jackknife err data of '+str(n_datafiles)+' measurements of EE correlator, normalized to G_norm=G_free/(g^2*C_F) with improved distances for '+qcdtype+' '+conftype+'\n')
    outfile.write('# rows correspond to flow times, columns to dt = {1, ... , Ntau/2}\n')
    lpd.write_flow_times(outfile, flow_times)
    numpy.savetxt(outfile, EE_err)

#write bootstrap samples in files for each flowtime separately 
for i in range(n_flow):
    with open(outputfolder+"/btstrp_samples/EE_"+'{0:.4f}'.format(numpy.sqrt(flow_times[i]*8)/nt)+"_Nt"+str(nt)+"_btstrp_samples.dat", 'w') as outfile:
        outfile.write("# "+str(n_samples)+" bootstrap samples for EE correlator "+qcdtype+" "+conftype+"\n# one sample consists of "+str(n_datafiles)+" draws(measurements), which are then reduced\n")
        outfile.write("# flowtime="+str(flow_times[i])+" or flowradius="+str(numpy.sqrt(flow_times[i]*8)/nt)+'\n')
        outfile.write("# first column: \\tau_T, second column: G(tauT), third column: std_dev\n")
        for k in range(n_samples):
            EE_sample_single_flowtime = numpy.empty((nt_half, 3))
            for j in range(nt_half):
                EE_sample_single_flowtime[j,0] = lpd.tauT_imp[nt][j]
                EE_sample_single_flowtime[j,1] = EE_samples_mean[k][i,j]
                EE_sample_single_flowtime[j,2] = EE_samples_err[k][i,j]
            numpy.savetxt(outfile, EE_sample_single_flowtime)
            outfile.write('#\n')
#for i in range(n_flow):
    #with open(outputfolder+"/btstrp_samples/EE_"+'{0:.4f}'.format(numpy.sqrt(flow_times[i]*8)/nt)+"_Nt"+str(nt)+"_btstrp_samples_cov.dat", 'w') as outfile:
        #outfile.write("# covariance matrices for the "+str(n_samples)+" bootstrap samples for EE correlator "+qcdtype+" "+conftype+"\n# one sample consists of "+str(n_datafiles)+" draws(measurements), which are then reduced\n")
        #outfile.write("# flowtime="+str(flow_times[i])+" or flowradius="+str(numpy.sqrt(flow_times[i]*8)/nt)+'\n')
        #outfile.write("# rows/columns: \\tau_T \n")
        #for k in range(n_samples):
            #EE_sample_single_flowtime = numpy.empty((nt_half, nt_half))
            #for j in range(nt_half):
                #for l in range(nt_half):
                    #EE_sample_single_flowtime[j,l] = EE_samples_cov[k][i][j,l]
            #numpy.savetxt(outfile, EE_sample_single_flowtime)
            #outfile.write('#\n')

print("done reducing data for "+conftype)
