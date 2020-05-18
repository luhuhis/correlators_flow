import numpy
import sys
import re
import math
from latqcdtools import statistics

def read_args():
    try:
        qcdtype = sys.argv[1]
        conftype = sys.argv[2]
        
        if len(sys.argv) > 3:
            add_params = sys.argv[3:]
        else:
            add_params = None

        beta=re.sub(r'(^.*?)_b', '', conftype) ; beta=float(beta)/100000
        ns=re.sub(r'(^.*?)s', '', conftype) ; ns=int(re.sub(r'(^.*?)t(.*)', r'\1', ns))
        nt=re.sub(r'(^.*?)t', '', conftype) ; nt=int(re.sub(r'(^.*?)_b(.*)', r'\1', nt))
        nt_half=int(nt/2)
        if add_params is not None:
            return qcdtype, conftype, beta, ns, nt, nt_half, add_params
        else:
            return qcdtype, conftype, beta, ns, nt, nt_half
    except IndexError:
        exit("Invalid Arguments. Usage: script.py <qcdtype> <conftype>, e.g. script.py quenched s064t16_b0687361")

def norm_corr( tauT ):
    norm = math.pi**2 * ( math.cos(math.pi*tauT)**2 / math.sin(math.pi*tauT)**4 + 1/(3*math.sin(math.pi*tauT)**2))
    return norm

def compute_EE_mean(EE_data): #mean function for Haukes jackknife routine
    mean_EE_numerator_real = numpy.mean(EE_data[0], axis=0)
    mean_polyakov_real = numpy.mean(EE_data[1], axis=0)
    correlator_sample = numpy.zeros(mean_EE_numerator_real.shape)
    nflowtimes = mean_polyakov_real.shape[0]
    nconf = len(EE_data[0])
    #single_flow_corr = numpy.zeros(correlator_sample[0].shape)
    for k in range(nflowtimes):
        correlator_sample[k,:] = mean_EE_numerator_real[k,:] / (-6*mean_polyakov_real[k][0])
        
        
        #statistics.calc_cov(data)
    #here also compute covariance matrix and return it!
    #maybe us latqcdtools covariance matrix implementation? from latqcdtools import statistics
    return correlator_sample

#tauT_imp = {16:(6.883194e-02, 1.085070e-01, 1.616768e-01, 2.247303e-01, 2.914720e-01, 3.571978e-01, 4.178695e-01, 4.527297e-01), 
            #20:(5.506555e-02, 8.680590e-02, 1.293465e-01, 1.798279e-01, 2.334154e-01, 2.867992e-01, 3.389318e-01, 3.894339e-01, 4.362603e-01, 4.632426e-01), 
            #24:(4.588796e-02, 7.233829e-02, 1.077897e-01, 1.498636e-01, 1.945462e-01, 2.391203e-01, 2.828289e-01, 3.257494e-01, 3.679753e-01, 4.092103e-01, 4.476142e-01, 4.697792e-01),
            #28:(0.039332543108858704, 0.06200429713148807, 0.09239146163313365, 0.12845583966778362, 0.16676106589862263, 0.20498393407330573, 0.2424922712244153, 0.2793926631132959, 0.31587145378525494, 0.3520157792617203, 0.38777053007225404, 0.4227840462554289, 0.45543361809657334, 0.4742844242530861),
            #30:(3.671037e-02, 5.787068e-02, 8.623207e-02, 1.198924e-01, 1.556450e-01, 1.913227e-01, 2.263377e-01, 2.607949e-01, 2.948799e-01, 3.287019e-01, 3.622833e-01, 3.955382e-01, 4.281210e-01, 4.585112e-01, 4.760585e-01),
            #36:(0.030591978346158538, 0.048225570112454395, 0.07186010884057005, 0.0999106894769821, 0.12970558033931176, 0.15943980536845198, 0.1886256825235924, 0.2173546715249435, 0.24579006824642974, 0.2740405351256582, 0.302163599109069, 0.3301833888944985, 0.3580983972194321, 0.385875234263269, 0.41341447770969214, 0.44041363885624496, 0.46560294274945346, 0.480147821801901)}

#nonimproved tauTs for testing purposes
tauT_imp = {16:list(numpy.arange(1/16, 0.501, 1/16)), 
            20:list(numpy.arange(1/20, 0.501, 1/20)), 
            24:list(numpy.arange(1/24, 0.501, 1/24)), 
            30:list(numpy.arange(1/30, 0.501, 1/30)),
            36:list(numpy.arange(1/36, 0.501, 1/36))} 

def inputfolder(qcdtype,conftype):
    return "../../data_merged/"+str(qcdtype)+"/"+str(conftype)+"/"

def create_folder(*paths):
    from pathlib import Path
    for path in paths:
        Path(path).mkdir(parents=True, exist_ok=True)
    return

"""load and reorganize data
organize the data in the format for analysistoolbox general jackknife routine
EE_numerator_real and polyakov_real: first index = different measurements, second index = flowtime, third index = tauT
EE_data: first index [0]=EE_numerator_real, [1]=polyakov_real, second index = different measurements (-> pairs of the form (EE_numerator_real[x], EE_polyakov_real[x]))
the jackknife leaves out pairs of (EE_numerator, Polyakov), this is done by the last argument (1), which specifies the axis on which the data pairs lie."""
def load_merged_data(qcdtype, conftype):
    inputfolder="../../data_merged/"+qcdtype+"/"+conftype+"/"
    nt=re.sub(r'(^.*?)t', '', conftype) ; nt=int(re.sub(r'(^.*?)_b(.*)', r'\1', nt))
    nt_half=int(nt/2)
        
    flow_times=numpy.loadtxt(inputfolder+"flowtimes_"+conftype+".dat")
    n_flow=len(flow_times)
    n_datafiles, n_streams=[int(i) for i in numpy.loadtxt(inputfolder+"n_datafiles_"+conftype+".dat")]
    EE_numerator_real_tmp=numpy.loadtxt(inputfolder+"EE_real_"+conftype+"_merged.dat")
    EE_numerator_real = []
    for i in range(n_datafiles):
        EE_numerator_real.append(EE_numerator_real_tmp[i*n_flow:(i+1)*n_flow,:])
    polyakov_real_tmp=numpy.loadtxt(inputfolder+"polyakov_real_"+conftype+"_merged.dat")
    polyakov_real = []
    for i in range(n_datafiles):
        polyakov_real.append(polyakov_real_tmp[i*n_flow:(i+1)*n_flow])
    polyakov_real = [ [ nt_half*[j] for j in i ] for i in polyakov_real] #copy the polyakov into every cell in order to match shapes for analysistoolbox bootstrap 
    EE_data = [EE_numerator_real, polyakov_real] 
    return flow_times, n_flow, n_datafiles, n_streams, EE_data 

def write_flow_times(outfile, flow_times):
    outfile.write('# flow times: ')
    for i in flow_times:
        outfile.write(str(i)+" ")
    outfile.write('\n')
    
ToverTc = {16:1.510424, 20:1.473469, 24:1.484769, 30:1.511761, 36: 1.504171}   
