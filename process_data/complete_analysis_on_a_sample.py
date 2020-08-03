#!/usr/local/bin/python
# perform complete analysis at a fixed flowtime for <nsamples> samples
import numpy
#import itertools
#import math
import sys
#from matplotlib import pyplot
#from latqcdtools import fitting
#from latqcdtools import spline_interpolate
import lib_plot as lp
import lib_process_data as lpd 

import _3_spline_interpolate
import _4_continuum_extr
import _5_flowtime_extr

try:
    current_sample_idx = int(sys.argv[1])
except IndexError:
    exit("Invalid Arguments. Usage: script.py <flowindex>, e.g. script.py 100")

use_imp = True
if use_imp:
    impstr = ""
else:
    impstr = "_unimp"

data = numpy.load("../../data_merged/quenched/samples_for_full_analysis/EE_sample_"+str(current_sample_idx)+".npy")

Ntaus = [36,30,24,20] #16
conftypes = {36:'s144t36_b0754400', 30:'s120t30_b0739400', 24:'s096t24_b0719200', 20:'s080t20_b0703500', 16:'s064t16_b0687361'} #16:'s064t16_b0687361',
max_nt_half = int(Ntaus[0]/2)

nflow = 221
flow_start_index = 50
flow_end_index = 118


EE        = numpy.loadtxt("../../data_merged/quenched/"+conftypes[Ntaus[0]]+"/EE_"+conftypes[Ntaus[0]]+".dat")
EE_err = numpy.loadtxt("../../data_merged/quenched/"+conftypes[Ntaus[0]]+"/EE_err_"+conftypes[Ntaus[0]]+".dat")
flow_extr_filter = []
for i in range(max_nt_half):
    rel_err = numpy.fabs(EE_err[:,i]/EE[:,i])
    for j,val in enumerate(rel_err):
        if val <= 0.0075 *lpd.tauT_imp[Ntaus[0]][i]:
            flow_extr_filter.append(j)
            break

cont_corr = numpy.empty((max_nt_half, 2, nflow)) #tauTs x (flowtime, corrvalue)
cont_corr[:] = numpy.nan

for flowindex in range(flow_start_index,flow_end_index):
    flowradius = lp.flow_radius[flowindex]
    flowtime = flowradius**2 / 8
    flowradius_str = '{0:.4f}'.format(flowradius)
    
    ### interpolation settings
    order = 3
    constraints = numpy.asarray([[_3_spline_interpolate.lower_tauT_limit(flowradius)-1/20,2,0],[0.5,1,0]])
    output_xdata, indices = _3_spline_interpolate.filter_tauTs(flowradius, lpd.tauT_imp[Ntaus[0]])
    output_xdata = output_xdata[1:] #remove the additional point that is only valid for helping the interpolation fit
    
    ### interpolate all coarser lattices to tauT of the finest lattice
    interpolations = {}
    tauTs = lpd.tauT_imp[Ntaus[0]]
    for i,Ntau in enumerate(Ntaus):
        nt_half = int(Ntau/2)
        xdata = data[flowindex][i][0]
        ydata = data[flowindex][i][1]
        for a in range(nt_half):
            ydata[a] = ydata[a] * lpd.improve_corr_factor(a,Ntau,flowindex, use_imp) 
        xdata, ydata, edata = _3_spline_interpolate.filter_corr_data(flowradius, nt_half, xdata, ydata)     
        if Ntau == Ntaus[0]: ### for the finest lattice we do not interpolate. exclude first datapoint to match interpolations!!
            interpolations[Ntau] = ydata[1:]
        else:
            try:
                knots = _3_spline_interpolate.get_knots(xdata)
                interpolations[Ntau] = _3_spline_interpolate.interpolate_EE_flow(xdata, ydata, edata, knots, order, constraints, output_xdata, return_params = False)[0]
            except:
                interpolations[Ntau] = None
    
    ### continuum extrapolation
    xdata = numpy.asarray([1/Ntau**2 for k,Ntau in enumerate(Ntaus) if interpolations[Ntau] is not None])
    offset = max_nt_half-len(output_xdata)
    for j in range(len(output_xdata)):
        ydata = [interpolations[Ntau][j] for Ntau in interpolations if interpolations[Ntau] is not None]
        cont_corr[j+offset,0,flowindex] = flowtime
        try:
            cont_corr[j+offset,1,flowindex] = _4_continuum_extr.fit_sample(ydata, xdata)[1]
        except:
            cont_corr[j+offset,1,flowindex] = numpy.nan

final_corr = numpy.empty((2,max_nt_half))
final_corr[:] = numpy.nan



### flow time zero extrapolation
for j,tauT in enumerate(lpd.tauT_imp[Ntaus[0]]):
    xdatatmp = cont_corr[j,0][flow_extr_filter[j]:]
    ydatatmp = cont_corr[j,1][flow_extr_filter[j]:]
    mask = numpy.isnan(ydatatmp)
    xdata = xdatatmp[~mask]
    ydata = ydatatmp[~mask]
    try:
        final_corr[0,j] = tauT
        final_corr[1,j] = _5_flowtime_extr.fit_sample(ydata, xdata)[1]
    except:
        final_corr[0,j] = numpy.nan
        final_corr[1,j] = numpy.nan
    
#save corr sample estimate in file.
#create script that evaluates all results and plots final corr estimate.
lpd.create_folder("../../data_merged/quenched/final_corr_samples"+impstr+"/")
numpy.savetxt("../../data_merged/quenched/final_corr_samples"+impstr+"/EE_final_"+str(current_sample_idx)+".txt", numpy.swapaxes(final_corr,0,1), header="#tauT    corr_sample_val")
fig, ax, plots = lp.create_figure(xlims=[0.15,0.51], ylims=[1.5, 4], xlabel=r'$\tau T$', ylabel=r'$ \frac{G (\tau)}{G^\mathrm{ norm } (\tau)}$', UseTex = False)
lp.plotstyle_add_point.update(dict(fmt='D-'))
ax.errorbar(final_corr[0], final_corr[1], color='black', **lp.plotstyle_add_point)
lpd.create_folder("../../plots/quenched/final_corr_samples"+impstr)
fig.savefig("../../plots/quenched/final_corr_samples"+impstr+"/EE_final_"+str(current_sample_idx)+".pdf", **lp.figurestyle)
        
    
