#!/usr/local/bin/python
import numpy
import sys
import lib_plot as lp
import lib_process_data as lpd

Ntaus = [36,30,24,20,16]
max_nt_half = int(Ntaus[0]/2)
conftypes = {36:'s144t36_b0754400', 30:'s120t30_b0739400', 24:'s096t24_b0719200', 20:'s080t20_b0703500', 16:'s064t16_b0687361'}
n_lat = len(Ntaus)
n_samples = 10000
n_flow = 221

all_data = numpy.zeros((n_samples,n_flow,n_lat,2,max_nt_half))
all_data[:] = numpy.nan
outputfolder = "../../data_merged/quenched/samples_for_full_analysis/"


for flowindex,flowradius in enumerate(lp.flow_radius):
    flowradius_str = '{0:.4f}'.format(flowradius)
    
    for i,Ntau in enumerate(Ntaus):
        tmp = numpy.loadtxt("../../data_merged/quenched/"+conftypes[Ntau]+"/btstrp_samples/EE_"+flowradius_str+"_Nt"+str(Ntau)+"_btstrp_samples.dat")
        nt_half = int(Ntau/2)
        for n in range(n_samples):
            if n % 1000 == 0:
                print(n)
            all_data[n,flowindex,i,0,:nt_half] = tmp[n*nt_half:(n+1)*nt_half,0]
            all_data[n,flowindex,i,1,:nt_half] = tmp[n*nt_half:(n+1)*nt_half,1]
        print(Ntau) 
    print("loaded flowindex ", flowindex)
    
for n in range(n_samples):
    numpy.save(outputfolder+"EE_sample_"+str(n)+".npy", all_data[n])
    print("saved ", n)
