#!/usr/local/bin/python
import numpy
import matplotlib
import lib_plot as lp
import lib_process_data as lpd
from latqcdtools import fitting
from latqcdtools import bootstr

def extrapolation_ansatz(x, m, b):
    return m * x + b 

def fit_sample(ydata, xdata, edata):
    fitter = fitting.Fitter(func=extrapolation_ansatz, xdata=xdata, ydata=ydata, edata=edata, func_sup_numpy=True, always_return = False)
    fitparams, fitparams_err, chi_dof= fitter.try_fit(algorithms = ["curve_fit"], start_params = [-100,3]) 
    return fitparams#, fitparams_err

def main():


    max_nt_half = 8
    max_nt = 16
    nflow = len(lp.flow_radius)
    ntauT = len(lpd.tauT_imp[max_nt])
    nsamples = 10
    
    results = {}
    conftype = "s064t16_b0687361"
    qcdtypes = ["quenched", "quenched_wilsonflow"]

    fmts = ['o-', 's-'] 
    
    ymax=3.8
    ymin=2.5
    UseTex = True
    
    ylabels = [r'$\displaystyle\frac{G^{\begin{subarray}{l}\mathrlap{\textrm{\tiny  Zeuthen}}\\[-0.3ex] \textrm{\tiny latt}\end{subarray}}_\tau {\scriptstyle(\tau_F)}}{G_{\tau,\tau_F=0}^{\begin{subarray}{l}\mathrlap{\textrm{\tiny  norm}}\\[-0.4ex] \textrm{\tiny latt}\end{subarray}}  } $',
               r'$\displaystyle\frac{G^{\begin{subarray}{l}\mathrlap{\textrm{\tiny  Wilson}}\\[-0.3ex] \textrm{\tiny latt}\end{subarray}}_\tau {\scriptstyle(\tau_F)}}{G_{\tau}^{\begin{subarray}{l}\mathrlap{\textrm{\tiny  norm}}\\[-0.4ex] \textrm{\tiny latt}\end{subarray}} {\scriptstyle(\tau_F)}} $']
    
    for qcdtype,ylabel in zip(qcdtypes, ylabels):
        fig, ax, plots = lp.create_figure(xlims=[-0.0001,0.0023], ylims=[ymin,ymax], xlabel=r'$\tau_F$', ylabel=ylabel, 
                                            xlabelpos=(0.96,0.05), ylabelpos=(0.08,0.96), UseTex = UseTex)
        
        
        inputfolder = lpd.inputfolder(qcdtype, conftype)   
        outputfolder = lpd.outputfolder(qcdtype, conftype)
        
        paths = [inputfolder+"/EE_"+conftype+".dat", 
                 inputfolder+"/EE_err_"+conftype+".dat", 
                 inputfolder+"/EE_"+conftype+"_flow_extr_quality.txt", 
                 inputfolder+"/EE_"+conftype+"_flow_extr.txt", 
                 outputfolder+"/EE_"+conftype+"_flow_extr_quality.pdf", 
                 outputfolder+"/EE_"+conftype+"_flow_extr.pdf"]
        
        EE = numpy.empty((ntauT,nflow,2))
        EE[:] = numpy.nan
        
        ### load data
        tmp1 = numpy.loadtxt(paths[0])
        tmp2 = numpy.loadtxt(paths[1])
        
        if qcdtype == "quenched":
            improve_with_flow = False
        elif qcdtype == "quenched_wilsonflow":
            improve_with_flow = True
            
        #apply improvement
        for i in range(150):
            for j in range(tmp1.shape[1]):
                tmp1[i,j] = tmp1[i,j] * lpd.improve_corr_factor(j, max_nt, i, improve=True, improve_with_flow=improve_with_flow)
                tmp2[i,j] = tmp2[i,j] * lpd.improve_corr_factor(j, max_nt, i, improve=True, improve_with_flow=improve_with_flow)
            
        #reorganize data
        for j,row in enumerate(tmp1):
            for i,val in enumerate(row):
                EE[i][j][0] = val
        for j,row in enumerate(tmp2):
            for i,val in enumerate(row):
                EE[i][j][1] = val
            
        if qcdtype == "quenched":
            ### filter low distance high error data that is not really linear
            flow_extr_filter = [] #this is the firts flowindex that shall be used
            for i in range(ntauT): #loop over tauT
                rel_err = numpy.fabs(EE[i,:,1]/EE[i,:,0])
                found = False
                for j,val in enumerate(rel_err): #loop over flowtimes
                    if val <= 0.006 *lpd.tauT_imp[max_nt][i]:
                        flow_extr_filter.append(j)
                        found = True
                        break
                if not found:
                    flow_extr_filter.append(-1)
            flowstarts = flow_extr_filter  
        
        results[qcdtype] = numpy.zeros((ntauT,3))
        mintauTindex = None
        with open(paths[2], 'w') as outfile:
            outfile.write('# flowtime zero extrapolation quality at fixed tauT for all tauT of Ntau=36 \n')
            outfile.write('# tauT    tau_F    G/G_norm     err \n')
            for i,tauT in enumerate(lpd.tauT_imp[max_nt]):
                flowstart = flowstarts[i]
                if not numpy.isnan(flowstart):
                    ### organize data
                    maxflowtime = (tauT-1/16)*numpy.sqrt(8*0.014)
                    flowend = min((numpy.abs(lp.flow_radius - maxflowtime)).argmin(), 134)
                    xdatatmp = lp.flow_radius[flowstart:flowend]**2 /8
                    ydatatmp = EE[i][flowstart:flowend,0]    
                    edatatmp = EE[i][flowstart:flowend,1]    
                    mask = numpy.isnan(ydatatmp)
                    xdata = xdatatmp[~mask]
                    ydata = ydatatmp[~mask]
                    edata = edatatmp[~mask]
                    if len(xdata) > 2:
                        ### perform extrapolation
                        fitparams, fitparams_err = bootstr.bootstr_from_gauss(fit_sample, data=ydata, data_std_dev=edata, numb_samples = nsamples, sample_size = 1, return_sample=False, args=[xdata, edata] )
                        results[qcdtype][i][0] = tauT
                        results[qcdtype][i][1] = fitparams[1]
                        results[qcdtype][i][2] = fitparams_err[1]
                        
                        ##plot extr quality
                        if mintauTindex is None:
                            mintauTindex = i
                        mycolor = lp.get_color(lpd.tauT_imp[max_nt], i,mintauTindex,max_nt_half)
                        lp.plotstyle_add_point_single.update(dict(fmt='-', markersize=2, mew=0.25))
                        ax.errorbar(xdata, ydata, edata, **lp.plotstyle_add_point_single, color=mycolor, zorder=-100+i)#, label='{0:.3f}'.format(tauT))
                        lp.plotstyle_add_point_single.update(dict(fmt=lp.markers[i-18], mew=0.25,markersize=5))
                        x = numpy.linspace(0,0.1,100)
                        plots.append(ax.errorbar(0, fitparams[1], fitparams_err[1], **lp.plotstyle_add_point_single, alpha=1, color=mycolor, zorder=1, label='{0:.3f}'.format(tauT)))
                        ax.errorbar(x, extrapolation_ansatz(x, *fitparams[0:2]), color=mycolor, alpha=1, fmt=':', lw=0.5, zorder=-100)
                        
                        numpy.savetxt(outfile, numpy.stack(([tauT for dummy in range(len(xdata)+1)], [0, *xdata], [fitparams[1], *ydata], [fitparams_err[1], *edata]), axis=-1))
                        outfile.write('# \n')
                    else:
                        results[qcdtype][i] = None
                else:
                    results[qcdtype][i] = None
                print("done ", tauT)
        numpy.savetxt(paths[3], results[qcdtype], header="tauT    G/Gnorm    err")
        
        ## second x-axis for flow radius
        ax2=ax.twiny()
        new_tick_locations = numpy.array([0,0.05**2/8,0.07**2/8,0.09**2/8,0.1**2/8, 0.11**2/8, 0.12**2/8, 0.13**2/8]) #0.06**2/8,
        def tick_function(X):
            V = numpy.sqrt(X*8)
            return ["%.1f" % z if z ==0 else "%.2f" % z for z in V]

        ax2.set_xlim(ax.get_xlim())
        ax2.set_xticks(new_tick_locations)
        ax2.set_xticklabels(tick_function(new_tick_locations))
        ax2.set_xlabel(r'$\sqrt{8\tau_F}T$', horizontalalignment='right', verticalalignment='top', bbox=lp.labelboxstyle, zorder=999999)
        ax2.xaxis.set_label_coords(0.99,0.97)        
        ax2.tick_params(direction='in', pad=0, width=0.5)

        ax.axvline(x=0, ymin=((results[qcdtype][mintauTindex,1]-ymin)/(ymax-ymin)), ymax=((results[qcdtype][max_nt_half-1,1]-ymin)/(ymax-ymin)), alpha=1, color='grey', zorder=-1000, lw=0.5, dashes=(5,2))
        
        lp.legendstyle.update(dict(handlelength=0.5))
        lp.legendstyle.update(dict(loc="lower right", bbox_to_anchor=(1.01,0.05), columnspacing=0.5, handlelength=0.75, labelspacing=0.1, handletextpad=0.2, borderpad=0, handler_map={matplotlib.container.ErrorbarContainer: matplotlib.legend_handler.HandlerErrorbar(yerr_size=0.4)} ))
        ax.legend(handles=plots)
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[::-1], labels[::-1],title=r'$\tau T=$', ncol=2, **lp.legendstyle).set_zorder(-1)
        matplotlib.pyplot.tight_layout(0)
        fig.savefig(paths[4])
        ax.lines.clear() ; ax.collections.clear() ; plots.clear();
        
        results[qcdtype] = numpy.swapaxes(results[qcdtype],0,1)
    
    
    lp.labelboxstyle.update(dict(alpha=0.8))
    ylabel = r'$\displaystyle \frac{\left.G^{\begin{subarray}{l}\mathrlap{\textrm{\tiny  Wilson}}\\[-0.3ex] \textrm{\tiny latt}\end{subarray}}{\scriptstyle(\tau,\tau_F)}/G^{\begin{subarray}{l}\mathrlap{\textrm{\tiny  norm}}\\[-0.4ex] \textrm{\tiny latt}\end{subarray}}{\scriptstyle(\tau, \tau_F)} \right|_{\tau_F \rightarrow 0} }{\left.G^{\begin{subarray}{l}\mathrlap{\textrm{\tiny  Zeuthen}}\\[-0.3ex] \textrm{\tiny latt}\end{subarray}}{\scriptstyle(\tau,\tau_F)}/G_{ {\tau}_{ F }=0}^{\begin{subarray}{l}\mathrlap{\textrm{\tiny  norm}}\\[-0.4ex] \textrm{\tiny latt}\end{subarray}} {\scriptstyle(\tau)}\right|_{\tau_F \rightarrow 0} } $'
    fig_final, ax_final, plots_final = lp.create_figure(xlims=[0.25,0.51], ylims=[0.985, 1.04], xlabelpos=(0.95,0.07), ylabelpos=(0.04,0.97), xlabel=r'$\tau T$', ylabel=ylabel, UseTex=UseTex)
    lp.legendstyle.update(dict(loc="lower left") , bbox_to_anchor=(0.1,0.1))
    lp.plotstyle_add_point.update(dict(fmt='o-'))


    #for qcdtype,label,fmt in zip(qcdtypes, labels, fmts):
        #lp.plotstyle_add_point.update(dict(fmt=fmt))
        #ax_final.errorbar(results[qcdtype][0], results[qcdtype][1], results[qcdtype][2], label=label, **lp.plotstyle_add_point)
    
    ax_final.axhline(y=1, xmin=0, xmax=1, alpha=1, color='grey', zorder=-1000, lw=0.5, dashes=(5,2))
    ydata = results["quenched_wilsonflow"][1]/results["quenched"][1]
    edata = ydata * numpy.sqrt( (results["quenched"][2]/results["quenched"][1])**2 + (results["quenched_wilsonflow"][2]/results["quenched_wilsonflow"][1])**2 )
    ax_final.errorbar(results[qcdtype][0], ydata, edata, **lp.plotstyle_add_point)
    
    #ax_final.legend(**lp.legendstyle)
    lp.titlestyle.update(dict(x=0.91, y=0.95))
    ax_final.set_title(r'$N_\tau = 16$', **lp.titlestyle)
    matplotlib.pyplot.tight_layout(0)
    fig_final.savefig(paths[5])


if __name__ == '__main__':
    main()
