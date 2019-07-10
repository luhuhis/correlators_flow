#!/usr/bin/python3
import math
import numpy
from os import listdir
#import own stuff
import sys
from latqcdtools import jackknife

#---------
#functions
#---------
def compute_EE_mean(EE_data): #mean function for Haukes jackknife routine
    mean_EE_numerator_real = numpy.mean(EE_data[0], axis=0)
    mean_polyakov_real      = numpy.mean(EE_data[1], axis=0)
    for k in range(0,mean_polyakov_real.shape[0]):
            mean_EE_numerator_real[k,:] /= -6*mean_polyakov_real[k][0]
    return mean_EE_numerator_real

def free_corr( tauT ):
    norm = math.pi**2 * ( math.cos(math.pi*tauT)**2 / math.sin(math.pi*tauT)**4 + 1/(3*math.sin(math.pi*tauT)**2))
    return norm

def free_corr_flow( tauT, flowtime, Nt ):
    result = 0
    for n in range(-50,50):
        x = tauT*Nt + n*Nt
        xi_sq = x**2 / (8 * flowtime)
        result += 1/x**4  * ((xi_sq**2 + xi_sq + 1)*math.exp(-xi_sq) - 1)
    result *= - 1 / math.pi**2 * Nt**4
    return result

def normalize_EE( data, data_err, tauT_imp, beta ):
    datashape = data.shape
    Ntau = datashape[1]*2
    for i in range(0, datashape[0]):
        for j in range(0, datashape[1]):
            data[i,j] = data[i,j] / free_corr(tauT_imp[j]) * Ntau**4 #* (1+0.079*6/beta) 
            data_err[i,j] = data_err[i,j] / free_corr(tauT_imp[j]) * Ntau**4 #* (1+0.079*6/beta) 

def normalize_EE_flow( data, data_err, tauT_imp, beta, flowtimes ):
    datashape = data.shape
    Ntau = datashape[1]*2    
    for i in range(0, datashape[0]):
        for j in range(0, datashape[1]):
            if i == 0 :
                data[i,j] = data[i,j] / free_corr(tauT_imp[j]) * Ntau**4 #* (1+0.079*6/beta) 
                data_err[i,j] = data_err[i,j] / free_corr(tauT_imp[j]) * Ntau**4 #* (1+0.079*6/beta) 
            else:
                data[i,j] = data[i,j] / free_corr_flow(tauT_imp[j], flowtimes[i], Ntau ) * Ntau**4 #* (1+0.079*6/beta) 
                data_err[i,j] = data_err[i,j] / free_corr_flow(tauT_imp[j], flowtimes[i], Ntau ) * Ntau**4 #* (1+0.079*6/beta) 

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

print(qcdtype, conftype, beta, ns, nt)
outputfolder="../data_merged/"+qcdtype+"/"+conftype+"/"
inputfolder="../../data_raw/"+qcdtype+"/"+conftype+"/"

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
                EE_numerator_imag.append(tmp[:,int((4 + nt_half)):])
                n_datafiles += 1
                
                

flow_radius = numpy.sqrt(flow_times*8)/int(nt)
for i in range(0,len(flow_radius)):
    flow_radius[i] = round(flow_radius[i],4)    
    
tauT_imp = {16:(6.883194e-02, 1.085070e-01, 1.616768e-01, 2.247303e-01, 2.914720e-01, 3.571978e-01, 4.178695e-01, 4.527297e-01), 
            20:(5.506555e-02, 8.680590e-02, 1.293465e-01, 1.798279e-01, 2.334154e-01, 2.867992e-01, 3.389318e-01, 3.894339e-01, 4.362603e-01, 4.632426e-01), 
            24:(4.588796e-02, 7.233829e-02, 1.077897e-01, 1.498636e-01, 1.945462e-01, 2.391203e-01, 2.828289e-01, 3.257494e-01, 3.679753e-01, 4.092103e-01, 4.476142e-01, 4.697792e-01),
            64:(0.015625, 0.03125, 0.046875, 0.0625, 0.078125, 0.09375, 0.109375, 0.125, 0.140625, 0.15625, 0.171875, 0.1875, 0.203125, 0.21875, 0.234375, 0.25, 0.265625, 0.28125, 0.296875, 0.3125, 0.328125, 0.34375, 0.359375, 0.375, 0.390625, 0.40625, 0.421875, 0.4375, 0.453125, 0.46875, 0.484375, 0.5)}

#-------------------------------------------
#do mean&jackknife and save results to disk
#-------------------------------------------
EE_data=[EE_numerator_real, [ [ nt_half*[j] for j in i ] for i in polyakov_real]] #get data in the format for Haukes general jackknife routine
EE, EE_err = jackknife.jackknife(compute_EE_mean, EE_data, n_streams, 1)
#EE_flow = numpy.copy(EE); EE_flow_err = numpy.copy(EE_err)
normalize_EE(EE, EE_err, tauT_imp[nt], beta) 

numpy.savetxt(outputfolder+"EE_"+conftype+".dat", EE) 
numpy.savetxt(outputfolder+"EE_err_"+conftype+".dat", EE_err)
numpy.savetxt(outputfolder+"tauT_imp"+conftype+".dat", numpy.array(tauT_imp[nt]))
numpy.savetxt(outputfolder+"n_datafiles_"+conftype+".dat", numpy.array(n_datafiles).reshape(1,))
numpy.savetxt(outputfolder+"flowradius_"+conftype+".dat", flow_radius)

#normalize_EE_flow(EE_flow, EE_flow_err, tauT_imp[nt], beta, flow_times) 
#numpy.savetxt(outputfolder+"EE_flow_"+conftype+".dat", EE_flow) 
#numpy.savetxt(outputfolder+"EE_flow_err_"+conftype+".dat", EE_flow_err)


#save data in correct format for continuum limit extrapolation via analysistoolbox's extrapolator.py
for i in range(0,len(flow_radius)):
    EE_singleflowradius = numpy.empty((nt_half,3))
    for j in range(0, nt_half):
        #TODO: only put in tauT's for which the flowtime is smaller than tauT*numpy.sqrt(8*0.014)=flowlimit 
        if (tauT_imp[nt][j]*numpy.sqrt(8*0.014)*0.5 >= flow_radius[i]) and (EE_err[i,j] < 0.1):
            EE_singleflowradius[j,0] = tauT_imp[nt][j]
            EE_singleflowradius[j,1] = EE[i,j]
            EE_singleflowradius[j,2] = EE_err[i,j]
    numpy.savetxt(outputfolder+"/continuum_limit/EE_"+str(flow_radius[i])+"_Nt"+str(nt_half)+".dat", EE_singleflowradius)
