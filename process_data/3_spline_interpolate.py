#!/usr/local/bin/python
### this script interpolates the EE correlator using the average of third order splines. it does that for one single given flow time.

import numpy
import itertools
import math
from matplotlib import pyplot
from latqcdtools import fitting
from latqcdtools import spline_interpolate
import lib_plot as lp
import lib_process_data as lpd

### manual knot function
def symmetric_knots(xdata, points_per_knot):
    knots = []
    ndata = len(xdata)
    nknots = int(numpy.round(ndata/points_per_knot, 0))
    interval_len_mean = ndata/(nknots+1)
    
    knot_slots = []
    for i in range(ndata-1):
        knot_slots.append((xdata[i]+xdata[i+1])/2)
    
    smallest_std_dev = numpy.inf
    knots_with_smallest_avg_dist = []
    for knot_combination in itertools.product(knot_slots, repeat=nknots):
        if len(numpy.unique(knot_combination)) != nknots:
            continue 
        if len(knot_combination) > 1:
            if numpy.any(numpy.diff(knot_combination) <= 0):
                continue
        
        interval_len_std_dev = 0
        for j in range(nknots):
            add = 0
            if j == 0:
                for x in xdata:
                    if x < knot_combination[j]:
                       add += 1
                interval_len_std_dev += ( add - interval_len_mean )**2
            elif j == nknots - 1:
                for x in xdata:
                    if x > knot_combination[j]:
                        add += 1
                interval_len_std_dev += ( add - interval_len_mean )**2
            else:
                for x in xdata:
                    if ( x > knot_combination[j] ) and ( x < knot_combination[j+1] ):
                        add += 1
                interval_len_std_dev += ( add - interval_len_mean )**2
        if numpy.isclose(interval_len_std_dev, smallest_std_dev):
            knots_with_smallest_avg_dist.append(knot_combination)
        else:
            if interval_len_std_dev < smallest_std_dev:
                smallest_std_dev = interval_len_std_dev
                knots_with_smallest_avg_dist = [knot_combination]
    return knots_with_smallest_avg_dist


### fit ansatz
def mysplinefunc(x, params, knots, order, constraints):
    return spline_interpolate.constraint_spline(x, knots=knots, coeffs=params, base_point=0.001, order=order, constraints=constraints)

### called for every set of spline parameters
def perform_spline_fit(knots, xdata, ydata, edata, constraints):
    fitter = fitting.Fitter(func=mysplinefunc, args = [knots, order, constraints], xdata=xdata, ydata=ydata, edata=edata, expand=False, always_return = True)
    mystart_params=[1 for i in range(len(knots)+order-1)]
    fitparams, fitparams_err, chi_dof, cov = fitter.do_fit(start_params = mystart_params, ret_pcov = True) 
    return fitparams


### compute average of all splines with different nknots in one sample. then compute average and std_dev of all samples.
def average_spline(x, all_fitparams, all_knots, order, constraints):
    nsamples = len(all_fitparams)
    values = []
    for i,fitparams,knots in zip(range(nsamples),all_fitparams,all_knots):
        counter = 0
        mean_single_sample = 0
        for params,knots in zip(fitparams, knots):
            if params is not None and knots is not None:
                    mean_single_sample += mysplinefunc(x, params, knots, order, constraints)
                    counter += 1
        try:
            mean_single_sample /= counter
        except:
            exit("Error: in at least one sample there are not enough points to fit any given spline")
        values.append(mean_single_sample)
    mean = numpy.mean(values, axis=0)
    std_dev = numpy.std(values, axis=0)
    return mean, std_dev


### start of script
qcdtype, conftype, beta, ns, nt, nt_half, add_params = lpd.read_args() 
inputfolder=lpd.inputfolder(qcdtype,conftype)
outputfolder="../../plots/"+qcdtype+"/"+conftype+"/interpolations/"
lpd.create_folder(outputfolder, inputfolder+"/interpolations/")

### load and set flowradius
flow_radii = numpy.loadtxt(lpd.inputfolder(qcdtype,conftype)+"/flowradii_"+conftype+".dat")
nflow = len(flow_radii)
flowindex = int(add_params[0])
flowradius = flow_radii[flowindex]
flowradiusstr = '{0:.4f}'.format(flow_radii[flowindex])
if flowradius < 1/nt:
    exit("skipped "+flowradiusstr+" for nt = "+str(nt))

### load data
EE_samples=numpy.loadtxt(inputfolder+"/btstrp_samples/EE_"+flowradiusstr+"_Nt"+str(nt)+"_btstrp_samples.dat")
nsamples = 250
### nsamples = int(len(EE_samples)/nt_half)

### spline parameters
order = 3
constraints = numpy.asarray([[flowradius/numpy.sqrt(8*0.014),2,0],[0.5,1,0]])
#knotconstant = {16:1, 20:2, 24:3, 30:4, 36:5}

### result variables
fitparams = [[None, None, None] for dummy in range(nsamples)]
theknots = [[None, None, None] for dummy in range(nsamples)]

### select trustworthy tauTs for this flow time
indices = numpy.asarray([k for k,j in enumerate(lpd.tauT_imp[nt]) if j > flowradius/numpy.sqrt(8*0.014)]) #r_F T < tauT/3
min_ind = numpy.min(indices)
if min_ind >= 1:
    ### add one helper point outside of flow limit to stabilize lower part of interpolation
    indices = numpy.asarray([min_ind-1, *indices])
xdata = numpy.asarray(lpd.tauT_imp[nt])
xdata = xdata[indices]
#if len(xdata) < 4:
    #exit("skip "+str(nt)+" "+str(flowradius))

### perform spline fits for each sample
for m in range(nsamples):
    
    ### select corr values
    ydata = numpy.asarray(EE_samples[m*nt_half:(m+1)*nt_half,1])
    for a in range(nt_half):
        ydata[a] = ydata[a] * lpd.improve_corr_factor(a,nt,flowindex) 
    ydata = ydata[indices]
    
    ### drastically but smoothly reduce weight of the one helper point outside of the fit interval based on its distance to the flow limit
    edata = numpy.fmax([1 for dummy in range(len(xdata))], [(flowradius/numpy.sqrt(8*0.014)/tauT)**(nt_half) for tauT in xdata])
      
    ### calc knot positions and actually perform fits
        #for nknots in range(2,5):
        #if (order+nknots+1 -2 +1) >= (len(xdata) ): 
            #break
    #nknots = len(xdata)-(order+1-2+knotconstant[nt]) #spline has order+nknots+1 params. minus 2 for the two constraints. there should always be atleast one more datapoint than parameters for the smallest lattice. 
                                                     #for larger lattices we require more free datapoints
    
    
    #add a new not every 4 datapoints. for 5 datapoints, one averages over splines with 0,1,2 knots. for 7 data points it is 1,2,3. for 11 it is 2,3,4. all equally distributed.
    nknots = max(1,int(round(len(xdata[:])/4,0)))
    ### if number of parameters is gr/eq the number of datapoints, then skip, because we want atleast one degree of freedom in the fit
    if (nknots+1)+1+1 >= len(xdata):
        exit("skip "+str(nt)+" "+str(flowradius))
    
    #nknots = max(1,int(round(len(xdata[:])/3,0)))

    #theknots[m][0] = spline_interpolate.random_knots(xdata, nknots, constraints, randomization_factor = 0.5)
    if nknots >= 0:
        theknots[m][2] = spline_interpolate.auto_knots(xdata[:], nknots+1)
        fitparams[m][2] = perform_spline_fit(theknots[m][2], xdata, ydata, edata, constraints)
        theknots[m][1] = spline_interpolate.auto_knots(xdata[:], nknots)
        fitparams[m][1] = perform_spline_fit(theknots[m][1], xdata, ydata, edata, constraints)
    if nknots >= 1:
        theknots[m][0] = spline_interpolate.auto_knots(xdata[:], nknots-1)
        fitparams[m][0] = perform_spline_fit(theknots[m][0], xdata, ydata, edata, constraints)
        


### compute final spline average at all desired tauT (xpoints) and save
interpolation_points = numpy.arange(0,0.5,0.01)
xpoints = numpy.asarray([x for x in numpy.sort(numpy.unique([*interpolation_points, *lpd.tauT_imp[36], 0.5, flowradius/numpy.sqrt(8*0.014)])) if x >= flowradius/numpy.sqrt(8*0.014)]) #
ypoints, epoints = average_spline(xpoints, fitparams, theknots, order, constraints)
numpy.savetxt(inputfolder+"/interpolations/"+"EE_"+'{0:.3f}'.format(flowradius)+"_interpolation.txt", numpy.stack((xpoints, ypoints, epoints), axis=-1), header="tauT     G_tauF(tau)/Gnorm      err")


### plot interpolations and underlying data points
fig, ax, plots = lp.create_figure(xlims=[0,0.51], ylims=[1.4, 4], xlabel=r'$\tau T$', ylabel=r'$ \frac{G_{\tau_F} (\tau)}{G^\mathrm{ norm }_{\tau_F = 0} (\tau)}$', UseTex = False)
ax.set_title(r'$ \sqrt{8\tau_F}T = $'+'{0:.3f}'.format(flowradius)+", nknots = "+str(nknots))
ax.axvline(x=(flowradius/numpy.sqrt(8*0.014)), **lp.verticallinestyle)
ax.fill_between(xpoints, ypoints-epoints, ypoints+epoints, alpha=0.5)
ax.errorbar(xpoints, ypoints, fmt = '-', lw = 0.5, mew=0)
EE=numpy.loadtxt(inputfolder+"/EE_"+conftype+".dat")
EE_err=numpy.loadtxt(inputfolder+"/EE_err_"+conftype+".dat")
for a in range(len(EE[flowindex])):
        EE[flowindex,a] = EE[flowindex,a] * lpd.improve_corr_factor(a,nt,flowindex) 
        EE_err[flowindex,a] = EE_err[flowindex,a] * lpd.improve_corr_factor(a,nt,flowindex)  
ax.errorbar(lpd.tauT_imp[nt], EE[flowindex], EE_err[flowindex], **lp.plotstyle_add_point_single)
fig.savefig(outputfolder+"/EE_"+'{0:.3f}'.format(flowradius)+"_interpolation.pdf" , **lp.figurestyle)#"_"+str(nknots)+ 

print("done ", nt, flowradius)
