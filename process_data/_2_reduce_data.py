#!/usr/local/bin/python
import numpy
import math 
from os import listdir 
import lib_process_data as lpd

#from latqcdtools import jackknife
from latqcdtools import bootstr
from latqcdtools import logger

"""
load and reorganize data

organize the data in the format for analysistoolbox general jackknife routine
XX_numerator_real and polyakov_real: first index = different measurements, second index = flowtime, third index = tauT
XX_data: first index [0]=XX_numerator_real, [1]=polyakov_real, second index = different measurements (-> pairs of the form (XX_numerator_real[x], XX_polyakov_real[x]))
the bootstrap leaves out pairs of (XX_numerator, Polyakov), this is done by the last argument (1), which specifies the axis on which the data pairs lie.
"""
def load_merged_data(qcdtype, corr, conftype):
    inputfolder=lpd.get_merged_data_path(qcdtype, corr, conftype)
    nt, nt_half = lpd.parse_conftype(conftype)[2:]
    
    print("read  "+inputfolder+"flowtimes_"+conftype+".dat")
    flow_times=numpy.loadtxt(inputfolder+"flowtimes_"+conftype+".dat")
    n_flow=len(flow_times)
    
    print("read  "+inputfolder+"n_datafiles_"+conftype+".dat")
    n_datafiles, n_streams=[int(i) for i in numpy.loadtxt(inputfolder+"n_datafiles_"+conftype+".dat")]
    print("read  "+inputfolder+corr+"_real_"+conftype+"_merged.dat")
    XX_numerator_real_tmp=numpy.loadtxt(inputfolder+corr+"_real_"+conftype+"_merged.dat")
    XX_numerator_real = []
    for i in range(n_datafiles):
        XX_numerator_real.append(XX_numerator_real_tmp[i*n_flow:(i+1)*n_flow,:])
        
    print("read  "+inputfolder+"polyakov_real_"+conftype+"_merged.dat")
    polyakov_real_tmp=numpy.loadtxt(inputfolder+"polyakov_real_"+conftype+"_merged.dat")
    polyakov_real = []
    for i in range(n_datafiles):
        polyakov_real.append(polyakov_real_tmp[i*n_flow:(i+1)*n_flow])
    polyakov_real = [ [ nt_half*[j] for j in i ] for i in polyakov_real] #copy the polyakov into every cell in order to match shapes for analysistoolbox bootstrap 
    XX_data = [XX_numerator_real, polyakov_real] 
    return flow_times, n_flow, n_datafiles, n_streams, XX_data 


"""
function for bootstrap routine that computes an XX correlator (with XX=EE or BB) normalized by the polyakov loop. numerator data (i.e. --X--X--) is first index, polyakov loop is second index of data.
"""
def compute_XX_corr(data):
    normalization_factor = 1 # the multiplicities are now correctly taken care of inside of the ParallelGPUCode
    numerator_mean = numpy.mean(data[0], axis=0)
    denominator_mean = numpy.mean(data[1], axis=0)
    XX = numerator_mean/normalization_factor / denominator_mean
    numerator_std_dev = numpy.std(data[0], axis=0)
    denominator_std_dev = numpy.std(data[1], axis=0)
    XX_std_dev = numpy.sqrt((numerator_std_dev/(normalization_factor*denominator_mean))**2+(numerator_mean*denominator_std_dev/(normalization_factor*denominator_mean**2))**2)
    #Ntau_half = XX.shape[1]
    #Ntau = Ntau_half*2
    #nflow = XX.shape[0]
    #compute time slice covariance matrix for every flowtime
    #nt_cov = numpy.zeros((nflow,Ntau_half,Ntau_half))
    #for i in range(nflow):
        #for j in range(Ntau_half):
            #for k in range(Ntau_half):
                #for l in range(len(data)):
                    #XXj = data[0,l,i,j] / (-6 * data[1,l,i,j]) 
                    #XXk = data[0,l,i,k] / (-6 * data[1,l,i,k]) 
                    #nt_cov[i,j,k] += 1/len(data) * (XXj - XX[i,j]) * (XXk - XX[i,k])
    return XX, XX_std_dev

def compute_only_numerator(data):
    return numpy.mean(data[0], axis=0), numpy.std(data[0], axis=0)

def compute_only_denominator(data):
    return numpy.mean(data[1], axis=0), numpy.std(data[1], axis=0)

def main():
    
    #parse cmd line arguments
    parser, requiredNamed = lpd.get_parser()
    parser.add_argument('--part_obs', help="only compute part of the observable", choices=['numerator','polyakovloop'])
    parser.add_argument('--n_samples', help='number of bootstrap samples', type=int, default='200')
    args = parser.parse_args()
   
    beta, ns, nt, nt_half = lpd.parse_conftype(args.conftype)
    fermions, temp, flowtype = lpd.parse_qcdtype(args.qcdtype)
    
     
    inputfolder = lpd.get_merged_data_path(args.qcdtype,args.corr,args.conftype)
    outputfolder = inputfolder
    
    lpd.create_folder(outputfolder+"/btstrp_samples/")
    flow_times, n_flow, n_datafiles, n_streams, XX_data = load_merged_data(args.qcdtype, args.corr, args.conftype)
    XX_data = numpy.asarray(XX_data)
    tauTs = numpy.arange(1/nt, 0.501, 1/nt)
    n_samples = args.n_samples

    ### do separate bootstrap because only then we can get the covariance matrices for a single sample
    #numerator_samples, numerator, numerator_err = bootstr.bootstr(calc_mean, XX_data[0], numb_samples=n_samples_1, sample_size=n_datafiles, conf_axis=0, return_sample=True, same_rand_for_obs=True)
    #polyakov_samples, polyakov, polyakov_err = bootstr.bootstr(calc_mean, XX_data[1], numb_samples=n_samples_1, sample_size=n_datafiles, conf_axis=0, return_sample=True, same_rand_for_obs=True)
    #combined_samples = numpy.asarray([numerator_samples, polyakov_samples])


    reduce_function = compute_XX_corr
    file_prefix = args.corr
    if args.part_obs:
        if args.part_obs == "numerator":
            reduce_function = compute_only_numerator
            file_prefix += "_numerator"
        if args.part_obs == "polyakovloop":
            reduce_function = compute_only_denominator
            file_prefix += "_polyakovloop"
    #Use same_rand_for_obs=False in order to break up correlations between the observables
    XX_samples, XX, XX_err = bootstr.bootstr(reduce_function, XX_data, numb_samples=n_samples, sample_size=n_datafiles, conf_axis=1, return_sample=True, same_rand_for_obs=False)

        
    #from the samples and the mean now compute the covariance matrix

    ### extract the right data
    XX_samples_mean = numpy.asarray([i[0] for i in XX_samples])
    XX_samples_err = numpy.asarray([i[1] for i in XX_samples])
    #XX_samples_cov = numpy.asarray([i[1] for i in XX_samples])
    XX_mean = XX[0] #this is the bootstrap estimate for the mean of the correlator 
    XX_err = XX_err[0] #this is the bootstrap estimate for the error of the mean of the correlator

    #XX_cov = numpy.zeros((n_flow,nt_half,nt_half)) #calculate the covariance matrix of the average of the correlator from the samples
    #for i in range(n_flow):
            #for j in range(nt_half):
                #for k in range(nt_half):
                    #for l in range(n_samples):
                        #XX_cov[i,j,k] += 1/(n_samples-1) * (XX_samples_mean[l,i,j] - XX_mean[i,j]) * (XX_samples_mean[l,i,k] - XX_mean[i,k])



    #write XX_mean and XX_err to file
    filename = outputfolder+file_prefix+"_"+args.conftype+".dat"
    print("write "+filename)
    with open(filename, 'w') as outfile:
        outfile.write('# bootstrap mean of '+str(n_datafiles)+' measurements of '+args.corr+' correlator for '+args.qcdtype+' '+args.conftype+'\n')
        outfile.write('# rows correspond to flow times, columns to dt = {1, ... , Ntau/2}\n')
        lpd.write_flow_times(outfile, flow_times)
        numpy.savetxt(outfile, XX_mean)
    #with open(outputfolder+"XX_cov_"+args.conftype+".dat", 'w') as outfile:
        #outfile.write('# reduced covariance matrix of '+str(n_datafiles)+' measurements of '+args.corr+' correlator for '+args.qcdtype+' '+args.conftype+'\n')
        #outfile.write('# covariance matrix of the time slices for each flow time, rows and columns correspond to dt = {1, ... , Ntau/2}\n')
        #lpd.write_flow_times(outfile, flow_times)
        #for i in range(n_flow):
            #numpy.savetxt(outfile, XX_cov[i])
            #outfile.write('#\n')
    filename = outputfolder+file_prefix+"_err_"+args.conftype+".dat"
    print("write "+filename)
    with open(filename, 'w') as outfile:
        outfile.write('# bootstrap mean err of '+str(n_datafiles)+' measurements of '+args.corr+' correlator '+args.qcdtype+' '+args.conftype+'\n')
        outfile.write('# rows correspond to flow times, columns to dt = {1, ... , Ntau/2}\n')
        lpd.write_flow_times(outfile, flow_times)
        numpy.savetxt(outfile, XX_err)

    #write bootstrap samples in files for each flowtime separately 
    for i in range(n_flow):
        filename = outputfolder+"/btstrp_samples/"+file_prefix+"_"+'{0:.4f}'.format(numpy.sqrt(flow_times[i]*8)/nt)+"_Nt"+str(nt)+"_btstrp_samples.dat"
        print("write "+filename)
        with open(filename, 'w') as outfile:
            outfile.write("# "+str(n_samples)+" bootstrap samples for "+args.corr+" correlator "+args.qcdtype+" "+args.conftype+"\n# one sample consists of "+str(n_datafiles)+" draws(measurements), which are then reduced\n")
            outfile.write("# flowtime="+str(flow_times[i])+" or flowradius="+str(numpy.sqrt(flow_times[i]*8)/nt)+'\n')
            outfile.write("# first column: \\tau_T, second column: G(tauT), third column: std_dev\n")
            for k in range(n_samples):
                XX_sample_single_flowtime = numpy.empty((nt_half, 3))
                for j in range(nt_half):
                    XX_sample_single_flowtime[j,0] = tauTs[j]
                    XX_sample_single_flowtime[j,1] = XX_samples_mean[k][i,j]
                    XX_sample_single_flowtime[j,2] = XX_samples_err[k][i,j]
                numpy.savetxt(outfile, XX_sample_single_flowtime)
                outfile.write('#\n')
    #for i in range(n_flow):
        #with open(outputfolder+"/btstrp_samples/"+args.corr+"_"+'{0:.4f}'.format(numpy.sqrt(flow_times[i]*8)/nt)+"_Nt"+str(nt)+"_btstrp_samples_cov.dat", 'w') as outfile:
            #outfile.write("# covariance matrices for the "+str(n_samples)+" bootstrap samples for "+args.corr+" correlator "+args.qcdtype+" "+args.conftype+"\n# one sample consists of "+str(n_datafiles)+" draws(measurements), which are then reduced\n")
            #outfile.write("# flowtime="+str(flow_times[i])+" or flowradius="+str(numpy.sqrt(flow_times[i]*8)/nt)+'\n')
            #outfile.write("# rows/columns: \\tau_T \n")
            #for k in range(n_samples):
                #XX_sample_single_flowtime = numpy.empty((nt_half, nt_half))
                #for j in range(nt_half):
                    #for l in range(nt_half):
                        #XX_sample_single_flowtime[j,l] = XX_samples_cov[k][i][j,l]
                #numpy.savetxt(outfile, XX_sample_single_flowtime)
                #outfile.write('#\n')

    print("done reducing data for "+args.conftype+" "+args.corr)

if __name__ == '__main__':
    main()
