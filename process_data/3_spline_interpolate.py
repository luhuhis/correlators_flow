#!/usr/local/bin/python
#this script interpolates the correlator using the average of third order splines

from latqcdtools import spline_interpolate
import lib_plot as lp
import lib_process_data as lpd
import numpy
from latqcdtools import fitting
from latqcdtools import plotting
from matplotlib import pyplot
import lib_plot as pl

def mysplinefunc(x, params, knots, order, leftboundary):
    return spline_interpolate.constraint_spline(x, knots=knots, coeffs=params, base_point=0, order=order, constraints=[[leftboundary,2,0],[0.5,1,0]])

def perform_spline_fit(i, flowradius, nknots, xdata, ydata, edata, leftboundary):
    myknots=spline_interpolate.auto_knots(xdata, nknots)
    fitter = fitting.Fitter(func=mysplinefunc, args = [myknots, order, leftboundary], xdata=xdata, ydata=ydata, edata=edata, expand=False, always_return = True)#numpy.log(edata+1.1)
    mystart_params=[1 for i in range(len(myknots)+order-1)]
    fitparams, fitparams_err, chi_dof, cov = fitter.do_fit(start_params = mystart_params, ret_pcov = True) 
    return fitparams

def average_spline(x, xdata, fitparams, order, leftboundary):
    nsamples = len(fitparams)
    values = []
    for i,sampleparams in enumerate(fitparams):
        counter = 0
        mean_single_sample = 0
        for k,params in enumerate(sampleparams):
            if params is not None:
                myknots = spline_interpolate.auto_knots(xdata, k)
                mean_single_sample += mysplinefunc(x, params, myknots, order, leftboundary)
                counter += 1
        try:
            mean_single_sample /= counter
        except:
            exit("Error: in at least one sample there are not enough points to fit any given spline")
        values.append(mean_single_sample)
    mean = numpy.mean(values, axis=0)
    std_dev = numpy.std(values, axis=0)
    return mean, std_dev


qcdtype, conftype, beta, ns, nt, nt_half, add_params = lpd.read_args() 
flowindex = int(add_params[0])
inputfolder=lpd.inputfolder(qcdtype,conftype)
outputfolder="../../plots/"+qcdtype+"/"+conftype+"/interpolations/"

lpd.create_folder(outputfolder, inputfolder+"/interpolations/")

flow_radii = numpy.loadtxt(lpd.inputfolder(qcdtype,conftype)+"/flowradii_"+conftype+".dat")
nflow = len(flow_radii)
flowradius = flow_radii[flowindex]
flowradiusstr = '{0:.4f}'.format(flow_radii[flowindex])

EE=numpy.loadtxt(inputfolder+"/EE_"+conftype+".dat")
EE_err=numpy.loadtxt(inputfolder+"/EE_err_"+conftype+".dat")
for a in range(len(EE[flowindex])):
        EE[flowindex,a] = EE[flowindex,a] / lpd.norm_corr(lpd.tauT_imp[nt][a]) * nt**4
        EE_err[flowindex,a] = EE_err[flowindex,a] / lpd.norm_corr(lpd.tauT_imp[nt][a]) * nt**4

EE_samples=numpy.loadtxt(inputfolder+"/btstrp_samples/EE_"+flowradiusstr+"_Nt"+str(nt)+"_btstrp_samples.dat")
#nsamples = int(len(EE_samples)/nt_half)
nsamples = 1000

#cheatfactor = lpd.ToverTc[nt]/lpd.ToverTc[36]
cheatfactor = 1

order = 3
fitparams = [[None, None, None, None, None] for dummy in range(nsamples)]
indices = numpy.asarray([k for k,j in enumerate(lpd.tauT_imp[nt]) if j > flowradius/numpy.sqrt(8*0.014)]) #r_F T < tauT/3
#indices = numpy.asarray(range(nt_half)) #r_F T < tauT/3
min_ind = numpy.min(indices)
#if min_ind >= 3:
    #indices = numpy.asarray([min_ind-3, min_ind-2, min_ind-1, *indices])
#if min_ind >= 2:
    #indices = numpy.asarray([min_ind-2, min_ind-1, *indices])
if min_ind >= 1:
    indices = numpy.asarray([min_ind-1, *indices])
xdata = numpy.asarray(lpd.tauT_imp[nt])
xdata = xdata[indices]
print(xdata)
for m in range(nsamples):
    ydata = numpy.asarray(cheatfactor * EE_samples[m*nt_half:(m+1)*nt_half,1])
    for v in range(nt_half):
        ydata[v] = ydata[v] / lpd.norm_corr(lpd.tauT_imp[nt][v]) * nt**4
    ydata = ydata[indices]
    edata = numpy.fmax([1.0 for dummy in range(len(xdata))], [(flowradius/numpy.sqrt(8*0.014)/tauT)**11 for tauT in xdata])
    for nknots in range(1,4):
        if (order + nknots +1 -2 +1) >= (len(xdata) ): #spline has order+nknots+1 params. minus 2 for the two constraints. there should always be atleast one more datapoint than parameters.
            break
        fitparams[m][nknots] = perform_spline_fit(flowindex, flowradius, nknots, xdata, ydata, edata, leftboundary = flowradius/numpy.sqrt(8*0.014))
        #print("done sample no ", m, flowindex, nknots)

interpolation_points = numpy.arange(0,0.5,0.0001)
xpoints = numpy.asarray([x for x in numpy.sort(numpy.unique([*interpolation_points, *lpd.tauT_imp[36], 0.5])) if x >= flowradius/numpy.sqrt(8*0.014)]) #
ypoints, epoints = average_spline(xpoints, xdata, fitparams, order, leftboundary = flowradius/numpy.sqrt(8*0.014))
numpy.savetxt(inputfolder+"/interpolations/"+"EE_"+'{0:.3f}'.format(flowradius)+"_interpolation.txt", numpy.stack((xpoints, ypoints, epoints), axis=-1), header="tauT     G_tauF(tau)/Gnorm      err")


#plot interpolations
fig, ax, plots = pl.create_figure(xlims=[0,0.51], ylims=[1.4, 4], xlabel=r'$\tau T$', ylabel=r'$ \frac{G_{\tau_F} (\tau)}{G^\mathrm{ norm }_{\tau_F = 0} (\tau)}$', UseTex = False)
ax.set_title(r'$ \sqrt{8\tau_F}T = $'+'{0:.3f}'.format(flowradius))
ax.axvline(x=(flowradius/numpy.sqrt(8*0.014)), **pl.verticallinestyle)
ax.fill_between(xpoints, ypoints-epoints, ypoints+epoints, alpha=0.5)
ax.errorbar(xpoints, ypoints, fmt = '-', lw = 0.5, mew=0)
ax.errorbar(lpd.tauT_imp[nt], cheatfactor*EE[flowindex], cheatfactor*EE_err[flowindex], **pl.plotstyle_add_point_single)
fig.savefig(outputfolder+"/EE_"+'{0:.3f}'.format(flowradius)+"_interpolation.pdf" , **pl.figurestyle)#"_"+str(nknots)+ 


#EE_cov_tmp=numpy.loadtxt(inputfolder+"/EE_cov_"+conftype+".dat")
#EE_cov=numpy.zeros((nflow,nt_half,nt_half))

#for i in range(nflow):
    #EE_cov[i] = EE_cov_tmp[i*nt_half:(i+1)*nt_half,:]

#samples=numpy.loadtxt("../../data_merged/quenched/s144t36_b0754400/btstrp_samples/EE_0.0500_Nt36_btstrp_samples.dat")
#samples_cov=numpy.loadtxt("../../data_merged/quenched/s144t36_b0754400/btstrp_samples/EE_0.0500_Nt36_btstrp_samples_cov.dat")

#EE = samples[0,0:nt_half

#delete nth row: numpy.delete(mymatrix, n-1, 0)
#delete nth column: numpy.delete(mymatrix, n-1, 1)

#print("flowradius, lowestLOOCVscore, bestorder, bestnknots")
#for i,flowradius in enumerate(flow_radii):
    #bestorder=0
    #bestnknots=0
    #lowestLOOCVscore=numpy.inf
    #if i in (50,60,70,80,90,100):
        #for nknots in range(1,3):
            #for myorder in range(2,5):
                #xdata = list(lpd.tauT_imp[nt])
                #myknots=spline_interpolate.auto_knots(xdata, nknots)
                #loocv=0
                #for k in range(nt_half):
                    #ydata = list(EE[i,:])
                    #edata = list(EE_err[i,:])
                    #cov   = EE_cov[i]
                    #cov   = numpy.delete(cov, k, 0)
                    #cov   = numpy.delete(cov, k, 1)
                    #xdata = list(lpd.tauT_imp[nt])
                    
                    #xtest = xdata[k]
                    #ytest = ydata[k]
                    #etest = edata[k]
                    #del xdata[k]
                    #del ydata[k]
                    #del edata[k]
                    
                    #def mysplinefunc(x, params):
                        #return spline_interpolate.constraint_spline(x, knots=myknots, coeffs=params, base_point=0, order=myorder, constraints=[[0.5,1,0]])

                    #mystart_params=[1 for i in range(len(myknots)+myorder)]

                    #fitter = fitting.Fitter(func=mysplinefunc, xdata=xdata, ydata=ydata, edata=cov, expand=False, always_return = True)
                    #fitparams, fitparams_err, chi_dof, cov = fitter.do_fit(start_params = mystart_params, ret_pcov = True) 
                    #loocv += numpy.fabs(mysplinefunc(xtest, fitparams) - ytest)

                #loocv /= nt_half 
                #if lowestLOOCVscore > loocv: 
                    #lowestLOOCVscore = loocv
                    #bestorder = myorder
                    #bestnknots = nknots
                #print(i, nknots, myorder, '{0:.3f}'.format(loocv))
        #print(flowradius, lowestLOOCVscore, bestorder, bestnknots)
        
        ##do fit for whole datset with best nknots and order and show plot
        #xdata = list(lpd.tauT_imp[nt])
        #myknots=spline_interpolate.auto_knots(xdata, bestnknots)
        #def mysplinefunc(x, params):
            #return spline_interpolate.constraint_spline(x, knots=myknots, coeffs=params, base_point=0, order=bestorder, constraints=[[0.5,1,0]])
        #ydata = list(EE[i,:])
        #edata = list(EE_err[i,:])
        #cov   = EE_cov[i]
        #fitter = fitting.Fitter(func=mysplinefunc, xdata=xdata, ydata=ydata, edata=cov, expand=False, always_return = True)
        #mystart_params=[1 for i in range(len(myknots)+bestorder)]
        #fitparams, fitparams_err, chi_dof, cov = fitter.do_fit(start_params = mystart_params, ret_pcov = True) 
        #plotting.plot_dots(xdata, ydata, edata)
        #plotting.plot_func(mysplinefunc, args = (fitparams,), xmax=0.5)
        #pyplot.savefig(outputfolder+"/"+'{0:.3f}'.format(flowradius)+".pdf") 
        #pyplot.cla()
