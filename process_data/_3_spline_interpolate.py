#!/usr/local/bin/python
"""
this script interpolates the XX correlator using the average of third order splines. it does that for one single given flow time.
"""

import numpy
import itertools
import math
import matplotlib
from matplotlib import pyplot
import lib_process_data as lpd
import scipy.interpolate 


"""
remove datapoints that the flow has destroyed
"""
def filter_tauTs(flowradius, xdata):
    indices = numpy.asarray([k for k,j in enumerate(xdata) if j > lpd.lower_tauT_limit(flowradius) ])
    min_ind = numpy.min(indices)
    ### add one helper point outside of flow limit to stabilize lower part of interpolation
    if min_ind >= 1:
        indices = numpy.asarray([min_ind-1, *indices])
    xdata = numpy.asarray(xdata)
    xdata = xdata[indices]
    return xdata, indices

""" 
remove datapoints that the flow has destroyed
"""
def filter_corr_data(flowradius, nt_half, xdata, ydata):
    xdata, indices = filter_tauTs(flowradius, xdata) 
    ydata = numpy.asarray(ydata)
    ydata = ydata[indices]
    ### drastically but smoothly reduce weight of the one helper point outside of the fit interval based on its distance to the flow limit
    #edata = numpy.fmax([1 for dummy in range(len(xdata))], [ 1+5*(lpd.lower_tauT_limit(flowradius)-tauT)/(1/20) for tauT in xdata])
    edata = [1 for dummy in range(len(xdata))]
    return xdata, ydata, edata


def interpolate_XX_flow(xdata, ydata, output_xdata):
    spline = scipy.interpolate.CubicSpline(xdata, ydata, bc_type=((2, 0.0), (1, 0.0)))
    return spline(output_xdata)


def main():
    
    #parse cmd line arguments
    parser, requiredNamed = lpd.get_parser()
    #parser.add_argument('--part_obs', help="only compute part of the observable", choices=['numerator','polyakovloop'])
    parser.add_argument('--use_imp', help='whether to use tree-level improvement', type=bool, default=True)
    parser.add_argument('--nsamples', help="number of artifical gaussian bootstrap samples to generate", type=int, default=1000)
    parser.add_argument('--int_Nt', help='use tauT of this Nt as xdata for the interpolation output', type=int, default=36)
    parser.add_argument('--min_Nt', help="only perform intepolation if sqrt(8tF)T>1/min_Nt", default=16)
    requiredNamed.add_argument('--flow_index', help='which flow time to interpolate', type=int, required=True)
    
    
    args = parser.parse_args()
   
    beta, ns, nt, nt_half = lpd.parse_conftype(args.conftype)
    fermions, temp, flowtype = lpd.parse_qcdtype(args.qcdtype)
    
    merged_data_path = lpd.get_merged_data_path(args.qcdtype,args.corr,args.conftype)
    outputfolder_data = merged_data_path+"/interpolations/"
    outputfolder_plot = lpd.get_plot_path(args.qcdtype,args.corr,args.conftype)+"/interpolations/"
    lpd.create_folder(outputfolder_plot, outputfolder_data)

    ### load and set flowradius
    flow_radii = numpy.loadtxt(merged_data_path+"/flowradii_"+args.conftype+".dat")
    nflow = len(flow_radii)
    
    flowindex = int(args.flow_index)
    flowradius = flow_radii[flowindex]
    flowradiusstr = '{0:.4f}'.format(flow_radii[flowindex])
    
    ### load data
    XX_samples=numpy.loadtxt(merged_data_path+"/btstrp_samples/"+args.corr+"_"+flowradiusstr+"_Nt"+str(nt)+"_btstrp_samples.dat")
    ### args.nsamples = int(len(XX_samples)/nt_half)

    #order = 3 #order of the spline polynomials when using fit

    ### result variables
    #thefitparams = [[None, None, None] for dummy in range(args.nsamples)]
    theoutputdata = []
    xdata, dummy = filter_tauTs(flowradius, lpd.get_tauTs(nt))
    #theknots = get_knots(xdata)
    interpolation_points = numpy.arange(0,0.5,0.001)
    xpoints = numpy.asarray([x for x in numpy.sort(numpy.unique([*interpolation_points, *lpd.get_tauTs(args.int_Nt), 0.5, lpd.lower_tauT_limit(flowradius)])) if x >= lpd.lower_tauT_limit(flowradius, args.min_Nt) ]) #xdata[0]
    
    ### perform spline fits for each sample
    for m in range(args.nsamples):
        ydata = numpy.asarray(XX_samples[m*nt_half:(m+1)*nt_half,1])
        for a in range(nt_half):
            ydata[a] = ydata[a] * lpd.improve_corr_factor(a,nt,flowindex,args.use_imp)
        xdata, ydata, edata = filter_corr_data(flowradius, nt_half, lpd.get_tauTs(nt), ydata)
        #constraints = numpy.asarray([[xdata[0],2,0],[0.5,1,0]])        
        #thefitparams[m] = interpolate_XX_flow(xdata, ydata, edata, theknots, order, constraints, output_xdata = None, return_params = True)
        #theoutputdata.append(fit_interpolate_XX_flow(xdata, ydata, edata, theknots, order, constraints, output_xdata=xpoints, return_params =False))
        theoutputdata.append(interpolate_XX_flow(xdata, ydata, output_xdata=xpoints))
        #print("done sample ", m)
    
    
    ### compute spline average over all samples at all desired tauT (xpoints) and save
    ypoints = numpy.mean(theoutputdata, axis=0)
    epoints = numpy.std(theoutputdata, axis=0)
    #ypoints, epoints = average_spline(xpoints, thefitparams, theknots, order, constraints)
    numpy.savetxt(outputfolder_data+args.corr+"_"+'{0:.4f}'.format(flowradius)+"_interpolation.txt", numpy.stack((xpoints, ypoints, epoints), axis=-1), header="tauT     G_tauF(tau)/Gnorm      err")


    ### plot interpolations and underlying data points
    UseTex = False
    if UseTex:
        displaystyle = r'\displaystyle'
        ylabel = r'$'+displaystyle+r'\frac{ G^\mathrm{latt }_{\tau_F} (\tau)}{G^{\substack{ \text{\tiny  norm} \\[-0.4ex] \text{\tiny latt } } }_{\tau_F = 0} (\tau)}$'
    else:
        ylabel = 'G'
        displaystyle = ''
    fig, ax, plots = lpd.create_figure(xlims=[0,0.51], ylims=[1.4, 4], xlabel=r'$\tau T$', ylabel=ylabel, xlabelpos=(0.95,0.05), UseTex=UseTex) 
    nknot_str = ""
    #for knots in theknots:
        #if knots != [None]:
            #nknot_str = nknot_str+str(len(knots))+","    
    ax.set_title(r'$ \sqrt{8\tau_F}T = $'+'{0:.3f}'.format(flowradius))#+", nknots = "+nknot_str)
    ax.axvline(x=lpd.lower_tauT_limit(flowradius), **lpd.verticallinestyle)
    ax.fill_between(xpoints, ypoints-epoints, ypoints+epoints, alpha=0.5)
    ax.errorbar(xpoints, ypoints, fmt = '-', lw = 0.5, mew=0)
    XX=numpy.loadtxt(merged_data_path+"/"+args.corr+"_"+args.conftype+".dat")
    XX_err=numpy.loadtxt(merged_data_path+"/"+args.corr+"_err_"+args.conftype+".dat")
    for a in range(len(XX[flowindex])):
            XX[flowindex,a] = XX[flowindex,a] * lpd.improve_corr_factor(a,nt,flowindex,args.use_imp) 
            XX_err[flowindex,a] = XX_err[flowindex,a] * lpd.improve_corr_factor(a,nt,flowindex, args.use_imp)  
    ax.errorbar(lpd.get_tauTs(nt), XX[flowindex], XX_err[flowindex], **lpd.plotstyle_add_point_single)
    matplotlib.pyplot.tight_layout(0.2)
    filepath=outputfolder_plot+"/"+args.corr+"_"+'{0:.4f}'.format(flowradius)+"_interpolation.pdf"
    print("saving ", filepath)
    fig.savefig(filepath)#"_"+str(nknots)+ 

    #print("done ", nt, '{0:.3f}'.format(flowradius))

if __name__ == '__main__':
    main()
    
    
"""unused and probably deprecated code (used for averaged fitting smoothing splines using analysistoolbox constrained spline):"""
    
#from latqcdtools import spline_interpolate
#from latqcdtools import fitting
    
### fit ansatz
def mysplinefunc(x, params, knots, order, constraints):
    return spline_interpolate.constraint_spline(x, knots=knots, coeffs=params, base_point=0.001, order=order, constraints=constraints)

### called for every set of spline parameters
def perform_spline_fit(knots, xdata, ydata, edata, order, constraints):
    fitter = fitting.Fitter(func=mysplinefunc, args = [knots, order, constraints], xdata=xdata, ydata=ydata, edata=edata, expand=False, func_sup_numpy=True)
    mystart_params=[1 for i in range(len(knots)+order-1)]
    
    fitparams, fitparams_err, chi_dof= fitter.try_fit(algorithms=["curve_fit"], start_params = mystart_params) 
    return fitparams
    

def get_knots(xdata):
    knots = [[None]]
    nknots = len(xdata)-2
    if nknots >= 0:
        knots[0] = spline_interpolate.auto_knots(xdata, nknots)
    else:
        raise Exception("skip, too few datapoints")
    return knots

def fit_interpolate_XX_flow(xdata, ydata, edata, knots, order, constraints, output_xdata=None, return_params=False):
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

    
""" 
compute average of all splines with different nknots in one sample. then compute average and std_dev of all samples if there are more than one.
"""
def average_spline(x, all_fitparams, theknots, order, constraints):
    args.nsamples = len(all_fitparams)
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
