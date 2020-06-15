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
 
