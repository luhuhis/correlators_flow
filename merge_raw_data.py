#!/usr/bin/python3
import math
import numpy
from os import listdir
#import own stuff
import sys
from latqcdtools import jackknife
from latqcdtools import bootstr


def compute_EE_mean(EE_data): #mean function for Haukes jackknife routine
    mean_EE_numerator_real = numpy.mean(EE_data[0], axis=0)
    mean_polyakov_real      = numpy.mean(EE_data[1], axis=0)
    for k in range(0,mean_polyakov_real.shape[0]):
            mean_EE_numerator_real[k,:] /= -6*mean_polyakov_real[k][0]
    return mean_EE_numerator_real

def free_corr( tauT ):
    norm = math.pi**2 * ( math.cos(math.pi*tauT)**2 / math.sin(math.pi*tauT)**4 + 1/(3*math.sin(math.pi*tauT)**2))
    return norm


def normalize_EE( data, data_err, tauT_imp):
    datashape = data.shape
    Ntau = datashape[1]*2
    for i in range(0, datashape[0]):
        for j in range(0, datashape[1]):
            data[i,j] = data[i,j] / free_corr(tauT_imp[j]) * Ntau**4
            data_err[i,j] = data_err[i,j] / free_corr(tauT_imp[j]) * Ntau**4


def normalize_EE_samples( data, tauT_imp):
    for k in range(0,len(data)):
        datashape = data[k].shape
        Ntau = datashape[1]*2
        for i in range(0, datashape[0]):
            for j in range(0, datashape[1]):
                data[k][i,j] = data[k][i,j] / free_corr(tauT_imp[j]) * Ntau**4
                
#----
#main
#----
try:
    qcdtype = sys.argv[1]
    conftype = sys.argv[2]
    beta = float(sys.argv[3])
    ns = int(sys.argv[4])
    nt = int(sys.argv[5])
    nt_half = int(nt/2)
except IndexError:
    exit("Invalid Arguments: merge_raw_data.py <qcdtype> <conftype> <beta> <ns> <nt>")

print("merging", conftype)
outputfolder="../data_merged/"+qcdtype+"/"+conftype+"/"
inputfolder="../data_raw/"+qcdtype+"/"+conftype+"/"

#---------
#load data
#---------
EE_numerator_real, EE_numerator_imag, polyakov_real, polyakov_imag, flow_times = ([] for i in range(5))
n_datafiles, n_streams = (0,0)

for stream_folder in listdir(inputfolder):
    if stream_folder.startswith(conftype+"_"):
        n_streams += 1
        for datafile in listdir(inputfolder+"/"+stream_folder):
            if datafile.startswith("zeuthenFlow"): #_acc0.000010_sts0.001000_ColElecCorrTimeSlices_"):#+conftype+"_U"):
                path = inputfolder+"/"+stream_folder+"/"+datafile
                tmp = numpy.loadtxt(path) # read in values as a numpy.ndarray
                if n_datafiles == 0:
                    flow_times = tmp[:,0]
                polyakov_real.append(tmp[:,1])
                polyakov_imag.append(tmp[:,2])
                EE_numerator_real.append(tmp[:,3:int((3 + nt_half))])
                EE_numerator_imag.append(tmp[:,int((3 + nt_half)):])
                n_datafiles += 1
                
                

flow_radius = numpy.sqrt(flow_times*8)/nt
for i in range(0,len(flow_radius)):
    flow_radius[i] = round(flow_radius[i],4)    
    
tauT_imp = {16:(6.883194e-02, 1.085070e-01, 1.616768e-01, 2.247303e-01, 2.914720e-01, 3.571978e-01, 4.178695e-01, 4.527297e-01), 
            20:(5.506555e-02, 8.680590e-02, 1.293465e-01, 1.798279e-01, 2.334154e-01, 2.867992e-01, 3.389318e-01, 3.894339e-01, 4.362603e-01, 4.632426e-01), 
            24:(4.588796e-02, 7.233829e-02, 1.077897e-01, 1.498636e-01, 1.945462e-01, 2.391203e-01, 2.828289e-01, 3.257494e-01, 3.679753e-01, 4.092103e-01, 4.476142e-01, 4.697792e-01),
            28:(0.039332543108858704, 0.06200429713148807, 0.09239146163313365, 0.12845583966778362, 0.16676106589862263, 0.20498393407330573, 0.2424922712244153, 0.2793926631132959, 0.31587145378525494, 0.3520157792617203, 0.38777053007225404, 0.4227840462554289, 0.45543361809657334, 0.4742844242530861),
            30:(3.671037e-02, 5.787068e-02, 8.623207e-02, 1.198924e-01, 1.556450e-01, 1.913227e-01, 2.263377e-01, 2.607949e-01, 2.948799e-01, 3.287019e-01, 3.622833e-01, 3.955382e-01, 4.281210e-01, 4.585112e-01, 4.760585e-01),
            32:(0.034415975319896056,0.054253764481208354,0.08084259567163814,0.11239933766750176,0.14591798888855687,0.17936738975712976,0.21219781707809077,0.24450999410099183,0.2764833450186914,0.3082316375131209,0.3398021996562756,0.37118383234030866,0.4022790685779857,0.43275510318199356,0.4611842090525582,0.4775996080475127),
            64:(0.015625, 0.03125, 0.046875, 0.0625, 0.078125, 0.09375, 0.109375, 0.125, 0.140625, 0.15625, 0.171875, 0.1875, 0.203125, 0.21875, 0.234375, 0.25, 0.265625, 0.28125, 0.296875, 0.3125, 0.328125, 0.34375, 0.359375, 0.375, 0.390625, 0.40625, 0.421875, 0.4375, 0.453125, 0.46875, 0.484375, 0.5)}

#nonimproved tauTs for testing purposes
#tauT_imp = {16:list(numpy.arange(1/16, 0.51, 1/16)), 
            #20:list(numpy.arange(1/20, 0.51, 1/20)), 
            #24:list(numpy.arange(1/24, 0.51, 1/24)), 
            #30:list(numpy.arange(1/30, 0.51, 1/30))}

#-------------------------------------------
#do mean&jackknife and save results to disk
#-------------------------------------------
#EE_numerator_real/polyakov_real: first index = different measurements, second index = flowtime, third index = tauT
#EE_data: first index = [0]=EE_numerator_real, [1]=EE_polyakov_real, second index = different measurements (-> pairs of the form (EE_numerator_real[x], EE_polyakov_real[x]))
#the jackknife leaves out pairs of (EE_numerator, Polyakov), this is done by the last argument (1), which specifies the axis on which the data pairs lie.

if qcdtype == "quenched":
    EE_data=[EE_numerator_real, [ [ nt_half*[j] for j in i ] for i in polyakov_real]] #get data in the format for Haukes general jackknife routine
    EE, EE_err = jackknife.jackknife(compute_EE_mean, EE_data, n_streams, 1) 
if qcdtype == "hisq":
    EE_data=[EE_numerator_real, [ [ nt_half*[j] for j in i ] for i in polyakov_real]] #get data in the format for Haukes general jackknife routine
    EE, EE_err = bootstr.bootstr(func=compute_EE_mean, data=EE_data, numb_samples=1000, sample_size=n_datafiles, same_rand_for_obs = True, conf_axis=1, return_sample = False)


#EE, EE_err = jackknife_EE(EE_numerator_real, polyakov_real, n_streams) #old way via my own jackknife_EE function

normalize_EE(EE, EE_err, tauT_imp[nt]) 

numpy.savetxt(outputfolder+"EE_"+conftype+".dat", EE) 
numpy.savetxt(outputfolder+"EE_err_"+conftype+".dat", EE_err)
numpy.savetxt(outputfolder+"tauT_imp"+conftype+".dat", numpy.array(tauT_imp[nt]))
numpy.savetxt(outputfolder+"n_datafiles_"+conftype+".dat", numpy.array(n_datafiles).reshape(1,))
numpy.savetxt(outputfolder+"flowradius_"+conftype+".dat", flow_radius)

#----------------
#generate bootstrap samples, which are used in the cont extr
nsamples=1000
EE_samples, EE_bootstrap, EE_err_bootstrap = bootstr.bootstr(func=compute_EE_mean, data=EE_data, numb_samples=nsamples, sample_size=10000, same_rand_for_obs = True, conf_axis=1, return_sample = True)
normalize_EE_samples(EE_samples, tauT_imp[nt])

if conftype in ("s064t16_b0687361", "s080t20_b0703500", "s096t24_b0719200", "s120t30_b0739400"):
    for i in range(0,len(flow_radius)):
        EE_sample_single_flowtime = numpy.empty((nt_half*nsamples, 2))
        for k in range(0,nsamples):
            for j in range(0, nt_half):
                EE_sample_single_flowtime[k*nt_half+j,0] = tauT_imp[nt][j]
                EE_sample_single_flowtime[k*nt_half+j,1] = EE_samples[k][i,j]
        numpy.savetxt(outputfolder+"/single_flow/EE_"+'{0:.4f}'.format(flow_radius[i])+"_Nt"+str(nt)+"_btstrp_samples.dat", EE_sample_single_flowtime)
#----------------
print("done with", conftype)
