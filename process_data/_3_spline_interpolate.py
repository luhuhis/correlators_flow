#!/usr/local/bin/python
### this script interpolates the EE correlator using the average of third order splines. it does that for one single given flow time.

import numpy
import itertools
import math
import matplotlib
from matplotlib import pyplot
from latqcdtools import fitting
from latqcdtools import spline_interpolate
import lib_plot as lp
import lib_process_data as lpd
import scipy.interpolate 

def lower_tauT_limit(flowradius):
    return flowradius/numpy.sqrt(8*0.014)+1/20  #r_F T < tauT/3 - 1/20, the last term is to shift the flow limit one point of the coarsest lattice that we have.

### fit ansatz
def mysplinefunc(x, params, knots, order, constraints):
    return spline_interpolate.constraint_spline(x, knots=knots, coeffs=params, base_point=0.001, order=order, constraints=constraints)

### called for every set of spline parameters
def perform_spline_fit(knots, xdata, ydata, edata, order, constraints):
    fitter = fitting.Fitter(func=mysplinefunc, args = [knots, order, constraints], xdata=xdata, ydata=ydata, edata=edata, expand=False, func_sup_numpy=True)
    #some_params = [1.22070268e+00,  2.05907247e+00,  3.45373663e+01,  3.30913109e+02,
#-6.39897695e+02,  2.00868085e+02,  1.38420525e+03, -1.68906472e+03]
    #some_params = [1.77350694,  -5.62974886,  42.12872541, -64.39893562]
    #mystart_params=[some_params[i] for i in range(len(knots)+order-1)]
    mystart_params=[1 for i in range(len(knots)+order-1)]
    
    fitparams, fitparams_err, chi_dof= fitter.try_fit(algorithms=["curve_fit"], start_params = mystart_params) 
    return fitparams

### remove datapoints that the flow has destroyed
def filter_tauTs(flowradius, xdata):
    indices = numpy.asarray([k for k,j in enumerate(xdata) if j > lower_tauT_limit(flowradius) ])
    min_ind = numpy.min(indices)
    ### add one helper point outside of flow limit to stabilize lower part of interpolation
    if min_ind >= 1:
        indices = numpy.asarray([min_ind-1, *indices])
    #if min_ind >= 2:
        #indices = numpy.asarray([min_ind-2, *indices])
    xdata = numpy.asarray(xdata)
    xdata = xdata[indices]
    return xdata, indices

### remove datapoints that the flow has destroyed
def filter_corr_data(flowradius, nt_half, xdata, ydata):
    xdata, indices = filter_tauTs(flowradius, xdata) 
    ydata = numpy.asarray(ydata)
    ydata = ydata[indices]
    ### drastically but smoothly reduce weight of the one helper point outside of the fit interval based on its distance to the flow limit
    #edata = numpy.fmax([1 for dummy in range(len(xdata))], [ 1+5*(lower_tauT_limit(flowradius)-tauT)/(1/20) for tauT in xdata])
    edata = [1 for dummy in range(len(xdata))]
    return xdata, ydata, edata


### compute average of all splines with different nknots in one sample. then compute average and std_dev of all samples if there are more than one.
def average_spline(x, all_fitparams, theknots, order, constraints):
    nsamples = len(all_fitparams)
    values = []
    for i,fitparams in enumerate(all_fitparams):
        counter = 0
        mean_single_sample = 0
        for params,knots in zip(fitparams, theknots):
            if params is not None and knots is not None:
                mean_single_sample += mysplinefunc(x, params, knots, order, constraints)
                counter += 1
        try:
            mean_single_sample /= counter
        except:
            exit("Error: in at least one sample there are not enough points to fit any given spline")
        values.append(mean_single_sample)
    if len(values) == 1:
        return values[0], numpy.nan
    mean = numpy.mean(values, axis=0)
    std_dev = numpy.std(values, axis=0)
    return mean, std_dev

def get_knots(xdata):
    knots = [[None]]
    nknots = len(xdata)-2
    if nknots >= 0:
        knots[0] = spline_interpolate.auto_knots(xdata, nknots)
    else:
        raise Exception("skip, too few datapoints")
    return knots

def interpolate_EE_flow(xdata, ydata, edata, knots, order, constraints, output_xdata=None, return_params=False, fit=True):
    if not fit and output_xdata is not None:
        spline = scipy.interpolate.CubicSpline(xdata, ydata, bc_type=((2, 0.0), (1, 0.0)))
        output_ydata = spline(output_xdata)
        return output_ydata 
    else:        
        ### result variables
        fitparams = len(knots)*[None]
        ### if number of parameters is gr/eq the number of datapoints, then skip, because we want atleast one degree of freedom in the fit
        ### perform spline fits for different nknots
        for i,current_knots in enumerate(knots):
            if current_knots != [None]:
                fitparams[i] = perform_spline_fit(current_knots, xdata, ydata, edata, order, constraints)
        if output_xdata is not None:
            output_ydata, output_edata = average_spline(output_xdata, [fitparams], knots, order, constraints)
            if not return_params:
                return output_ydata, output_edata
        if return_params:
            if output_xdata is None:
                return fitparams
            else:
                return output_ydata, output_edata, fitparams
    

def main():
    use_imp = True
    
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
    #if flowradius < 1/nt:
        #exit("skipped "+flowradiusstr+" for nt = "+str(nt))

    ### load data
    EE_samples=numpy.loadtxt(inputfolder+"/btstrp_samples/EE_"+flowradiusstr+"_Nt"+str(nt)+"_btstrp_samples.dat")
    nsamples = 10000
    ### nsamples = int(len(EE_samples)/nt_half)

    order = 3

    ### result variables
    thefitparams = [[None, None, None] for dummy in range(nsamples)]
    theoutputdata = []
    xdata, dummy = filter_tauTs(flowradius, lpd.tauT_imp[nt])
    theknots = get_knots(xdata)
    interpolation_points = numpy.arange(0,0.5,0.001)
    xpoints = numpy.asarray([x for x in numpy.sort(numpy.unique([*interpolation_points, *lpd.tauT_imp[36], 0.5, lower_tauT_limit(flowradius)])) if x >= lower_tauT_limit(flowradius) ]) #xdata[0]
    
    ### perform spline fits for each sample
    for m in range(nsamples):
        ydata = numpy.asarray(EE_samples[m*nt_half:(m+1)*nt_half,1])
        for a in range(nt_half):
            ydata[a] = ydata[a] * lpd.improve_corr_factor(a,nt,flowindex,use_imp)
        xdata, ydata, edata = filter_corr_data(flowradius, nt_half, lpd.tauT_imp[nt], ydata)
        constraints = numpy.asarray([[xdata[0],2,0],[0.5,1,0]])        
        #thefitparams[m] = interpolate_EE_flow(xdata, ydata, edata, theknots, order, constraints, output_xdata = None, return_params = True)
        theoutputdata.append(interpolate_EE_flow(xdata, ydata, edata, theknots, order, constraints, output_xdata=xpoints, return_params =False, fit=False))
        #print("done sample ", m)
    
    
    ### compute spline average over all samples at all desired tauT (xpoints) and save
    ypoints = numpy.mean(theoutputdata, axis=0)
    epoints = numpy.std(theoutputdata, axis=0)
    #ypoints, epoints = average_spline(xpoints, thefitparams, theknots, order, constraints)
    numpy.savetxt(inputfolder+"/interpolations/"+"EE_"+'{0:.3f}'.format(flowradius)+"_interpolation.txt", numpy.stack((xpoints, ypoints, epoints), axis=-1), header="tauT     G_tauF(tau)/Gnorm      err")


    ### plot interpolations and underlying data points
    UseTex = False
    if UseTex:
        displaystyle = r'\displaystyle'
        ylabel = r'$'+displaystyle+r'\frac{ G^\mathrm{latt }_{\tau_F} (\tau)}{G^{\substack{ \text{\tiny  norm} \\[-0.4ex] \text{\tiny latt } } }_{\tau_F = 0} (\tau)}$'
    else:
        ylabel = 'G'
        displaystyle = ''
    fig, ax, plots = lp.create_figure(xlims=[0,0.51], ylims=[1.4, 4], xlabel=r'$\tau T$', ylabel=ylabel, xlabelpos=(0.95,0.05), UseTex=UseTex) 
    nknot_str = ""
    for knots in theknots:
        if knots != [None]:
            nknot_str = nknot_str+str(len(knots))+","    
    ax.set_title(r'$ \sqrt{8\tau_F}T = $'+'{0:.3f}'.format(flowradius)+", nknots = "+nknot_str)
    ax.axvline(x=lower_tauT_limit(flowradius), **lp.verticallinestyle)
    ax.fill_between(xpoints, ypoints-epoints, ypoints+epoints, alpha=0.5)
    ax.errorbar(xpoints, ypoints, fmt = '-', lw = 0.5, mew=0)
    EE=numpy.loadtxt(inputfolder+"/EE_"+conftype+".dat")
    EE_err=numpy.loadtxt(inputfolder+"/EE_err_"+conftype+".dat")
    for a in range(len(EE[flowindex])):
            EE[flowindex,a] = EE[flowindex,a] * lpd.improve_corr_factor(a,nt,flowindex,use_imp) 
            EE_err[flowindex,a] = EE_err[flowindex,a] * lpd.improve_corr_factor(a,nt,flowindex, use_imp)  
    ax.errorbar(lpd.tauT_imp[nt], EE[flowindex], EE_err[flowindex], **lp.plotstyle_add_point_single)
    matplotlib.pyplot.tight_layout(0.2)
    fig.savefig(outputfolder+"/EE_"+'{0:.3f}'.format(flowradius)+"_interpolation.pdf")#"_"+str(nknots)+ 

    print("done ", nt, flowradius)

if __name__ == '__main__':
    main()
