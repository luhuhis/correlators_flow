#!/usr/local/bin/python
import numpy
import matplotlib
import lib_plot as lp
import lib_process_data as lpd
from latqcdtools import fitting
from latqcdtools import bootstr

#def find_nearest(array, value):
    #array = np.asarray(array)
    #return idx = (np.abs(array - value)).argmin()


def extrapolation_ansatz(x, m, b):
    return m * x + b 

def linear(x, m, b):
    return m * x + b 

def fit_sample(ydata, xdata):
    fitter = fitting.Fitter(func=extrapolation_ansatz, xdata=xdata, ydata=ydata, always_return = True)
    fitparams, fitparams_err, chi_dof= fitter.do_fit(start_params = [-100,3]) 
    return fitparams

fig, ax, plots = lp.create_figure(xlims=[-0.00002,0.0015], ylims=[1.5, 4], xlabel=r'$\tau_F$', ylabel=r'$ \frac{G_\tau (\tau_F)}{G^\mathrm{ norm } }$', UseTex = False)
lpd.create_folder(lp.outputfolder+"/flow_extr_quality/")

nflow = len(lp.flow_radius)
ntauT = len(lpd.tauT_imp[36])
nsamples = 1000

### load continuum extr for each flowtime, only extract the desired tauT
EE_conts = [] 
for flowradius in lp.flow_radius:
    flowradius_str = '{0:.3f}'.format(flowradius)
    try:
        tmp = numpy.loadtxt(lp.inputfolder+"/cont_extr/EE_"+flowradius_str+"_cont.txt")
        indices = []
        for j,tauT in enumerate(tmp[:,0]):
            if tauT in lpd.tauT_imp[36]:
                indices.append(j)
        indices = numpy.asarray(indices)
        ### if any of the values in a row (tauT corr corr_err) is not a number, append None instead
        if numpy.isnan(tmp[indices]).any():
            EE_conts.append(None)
        else:
            EE_conts.append(tmp[indices])       
    except:
        EE_conts.append(None)


### rearrange data and put in nans where data is missing
EE = numpy.empty((ntauT,nflow,2))
EE[:] = numpy.nan

for j in range(nflow):
    if EE_conts[j] is not None:
        this_ntauT = len(EE_conts[j])
        offset = ntauT-this_ntauT
        for k in range(this_ntauT):
            EE[offset+k][j] = EE_conts[j][k][1:]
    else:
        for k in range(ntauT):
            EE[k][j] = [numpy.nan, numpy.nan]
        
        
results = numpy.zeros((ntauT,3))
flowstart=0
lpd.create_folder(lp.inputfolder+"/flow_extr_quality/")
mintauTindex = None
with open(lp.inputfolder+"/flow_extr_quality/EE_flow_extr_quality.txt", 'w') as outfile:
    outfile.write('# flowtime zero extrapolation quality at fixed tauT for all tauT of Ntau=36 \n')
    outfile.write('# tauT    tau_F    G/G_norm     err \n')
    for i,tauT in enumerate(lpd.tauT_imp[36]):
        maxflowtime = tauT*numpy.sqrt(8*0.014)
        flowend = min((numpy.abs(lp.flow_radius - maxflowtime)).argmin(), 150)
        xdata = lp.flow_radius[flowstart:flowend]**2 /8
        ydata = EE[i][flowstart:flowend,0]    
        edata = EE[i][flowstart:flowend,1]    
        mask = numpy.isnan(ydata)
        xdatatmp = xdata[~mask]
        ydatatmp = ydata[~mask]
        edatatmp = edata[~mask]
        xdata = [x for j,x in enumerate(xdatatmp) if edatatmp[j]/ydatatmp[j] < 0.025]
        ydata = [x for j,x in enumerate(ydatatmp) if edatatmp[j]/ydatatmp[j] < 0.025]
        edata = [x for j,x in enumerate(edatatmp) if edatatmp[j]/ydatatmp[j] < 0.025]
        if len(xdata) > 2:
            try:
                fitparams, fitparams_err = bootstr.bootstr_from_gauss(fit_sample, data=ydata, data_std_dev=edata*3, numb_samples = nsamples, sample_size = 1, return_sample=False, args=[xdata] )
                results[i][0] = tauT
                results[i][1] = fitparams[1]
                results[i][2] = fitparams_err[1]
                #plot extr quality
                if mintauTindex is None:
                    mintauTindex = i
                mycolor = lp.get_color(lpd.tauT_imp[36], i,mintauTindex,18)
                plots.append(ax.errorbar(xdata, ydata, edata, **lp.plotstyle_add_point_single, color=mycolor, zorder=-100+i, label='{0:.3f}'.format(tauT)))
                x = numpy.linspace(0,0.1,100)
                ax.errorbar(0, fitparams[1], fitparams_err[1], **lp.plotstyle_add_point_single, color=mycolor, zorder=1)
                ax.errorbar(x, linear(x, *fitparams[0:2]), color=mycolor, alpha=0.8, fmt='--', lw=0.5, zorder=-100)
                
                numpy.savetxt(outfile, numpy.stack(([tauT for dummy in range(len(xdata)+1)], [0, *xdata], [fitparams[1], *ydata], [fitparams_err[1], *edata]), axis=-1))
                outfile.write('# \n')
            except:
                results[i] = None
        else:
            results[i] = None
numpy.savetxt(lp.inputfolder+"/EE_final.txt", results, header="tauT    G/Gnorm    err")

lp.legendstyle.update(dict(handlelength=1))
lp.legendstyle.update(handlelength=0.5)
ax.legend(handles=plots, title=r'$\tau T$', **lp.legendstyle)
fig.savefig(lp.outputfolder+"/EE_flow_extr_quality.pdf" , **lp.figurestyle)
ax.lines.clear() ; ax.collections.clear() ; plots.clear();

fig, ax, plots = lp.create_figure(xlims=[0.15,0.51], ylims=[1.5, 4], xlabel=r'$\tau T$', ylabel=r'$ \frac{G (\tau)}{G^\mathrm{ norm } (\tau)}$', UseTex = False)
lp.plotstyle_add_point.update(dict(fmt='D-'))
results = numpy.swapaxes(results,0,1)
ax.errorbar(results[0], results[1], results[2], color='black', **lp.plotstyle_add_point)
fig.savefig(lp.outputfolder+"/EE_final.pdf" , **lp.figurestyle)
