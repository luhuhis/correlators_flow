import numpy
import sys
import re

def read_args():
    try:
        qcdtype = sys.argv[1]
        conftype = sys.argv[2]

        beta=re.sub(r'(^.*?)_b', '', conftype) ; beta=float(beta)/100000
        ns=re.sub(r'(^.*?)s', '', conftype) ; ns=int(re.sub(r'(^.*?)t(.*)', r'\1', ns))
        nt=re.sub(r'(^.*?)t', '', conftype) ; nt=int(re.sub(r'(^.*?)_b(.*)', r'\1', nt))
        nt_half=int(nt/2)
        return qcdtype, conftype, beta, ns, nt, nt_half
    except IndexError:
        exit("Invalid Arguments: script.py <qcdtype> <conftype>, e.g. script.py quenched s064t16_b0687361")

def compute_EE_mean(EE_data): #mean function for Haukes jackknife routine
    mean_EE_numerator_real = numpy.mean(EE_data[0], axis=0)
    mean_polyakov_real      = numpy.mean(EE_data[1], axis=0)
    for k in range(0,mean_polyakov_real.shape[0]):
            mean_EE_numerator_real[k,:] /= -6*mean_polyakov_real[k][0]
    return mean_EE_numerator_real

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

def inputfolder(qcdtype,conftype):
    return "../../data_merged/"+str(qcdtype)+"/"+str(conftype)+"/"

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
    polyakov_real = [ [ nt_half*[j] for j in i ] for i in polyakov_real]
    EE_data=[EE_numerator_real, polyakov_real] 
    return flow_times, n_flow, n_datafiles, n_streams, EE_data 

def write_flow_times(outfile, flow_times):
    outfile.write('# flow times: ')
    for i in flow_times:
        outfile.write(str(i)+" ")
    outfile.write('\n')
