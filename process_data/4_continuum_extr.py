#!/usr/local/bin/python
import numpy
import matplotlib
import sys
from latqcdtools import fitting
from latqcdtools import bootstr
import lib_plot as lp
import lib_process_data as lpd

def extrapolation_ansatz(x, m, b):
    return m * x + b

def fit_sample(ydata, xdata):
    fitter = fitting.Fitter(func=extrapolation_ansatz, xdata=xdata, ydata=ydata, always_return = True)
    fitparams, fitparams_err, chi_dof= fitter.do_fit(start_params = [80,3]) 
    return fitparams

try:
    flowindex = sys.argv[1]
    i = int(flowindex) #create alias
    flowradius = lp.flow_radius[i]
    flowradius_str = '{0:.3f}'.format(flowradius)
    if flowradius < 1/20:
        exit("need at least a flow radius of 1/20 to have enough lattices for cont extr")
except IndexError:
    exit("Invalid Arguments. Usage: script.py <flowindex>, e.g. script.py 100")
   

### load largest lattice
EE36 = numpy.asarray(([tauT for tauT in lpd.tauT_imp[36] if tauT > flowradius/numpy.sqrt(8*0.014)], 
                     [val*lpd.improve_corr_factor(j, 36, i) for j,val in enumerate(lp.EE_36[i]) if lpd.tauT_imp[36][j] > flowradius/numpy.sqrt(8*0.014)], [val*lpd.improve_corr_factor(j, 36, i) for j,val in enumerate(lp.EE_err_36[i]) if lpd.tauT_imp[36][j] > flowradius/numpy.sqrt(8*0.014)]))



### load interpolations
EE_ints = []
for conftype in ("s064t16_b0687361", "s080t20_b0703500", "s096t24_b0719200", "s120t30_b0739400", "s144t36_b0754400"):
    if conftype == "s064t16_b0687361":
        EE_ints.append(None)
        continue
    if conftype == "s144t36_b0754400":
        EE_ints.append(EE36)
    else:
        try:
            tmp = numpy.loadtxt(lp.inputfolder+"/"+conftype+"/interpolations/EE_"+flowradius_str+"_interpolation.txt", unpack=True)
            tmp2 = numpy.asarray(( [tauT for tauT in tmp[0] if tauT in lpd.tauT_imp[36]], 
                                    [val for j,val in enumerate(tmp[1]) if tmp[0][j] in lpd.tauT_imp[36]],
                                    [val for j,val in enumerate(tmp[2]) if tmp[0][j] in lpd.tauT_imp[36]] ))
            EE_ints.append(tmp2)
        except:
            EE_ints.append(None)
if EE_ints[1] is None:
    exit("skip")


### define some variables
#interpolation_tauTs = EE_ints[1][0]
nsamples = 1000
Ntaus = numpy.asarray([16,20,24,30,36])
n_tauTs = len(EE36[0])
tauTs = EE36[0]

results = numpy.empty((n_tauTs,3))
results[:] = numpy.nan
fitparameters = []
mintauTindex =  lpd.tauT_imp[36].index(EE36[0][0])
extrapolation_results = []
   
### plot settings
fig, ax, plots = lp.create_figure(xlims=[-0.0001,0.004], ylims=[1.4, 4], xlabel=r'$N_\tau^{-2}$', ylabel=r'$ \frac{G_{\tau_F} (\tau)}{G^\mathrm{ norm }_{\tau_F = 0} (\tau)}$', UseTex = False)
ax.set_title(r'$ \sqrt{8\tau_F}T = $'+flowradius_str)
lp.legendstyle.update(dict(labelspacing=0.25, handletextpad=0, borderpad=0, handler_map={matplotlib.container.ErrorbarContainer: matplotlib.legend_handler.HandlerErrorbar(xerr_size=0.4)} ))
lpd.create_folder(lp.outputfolder+"/cont_extr_quality/")


with open(lp.inputfolder+"/cont_extr_quality/EE_"+flowradius_str+"_cont_quality.txt", 'w') as outfile:
    outfile.write('# continuum extrapolation quality at fixed flowtime for all tauT of Ntau=36 \n')
    outfile.write('# tauT    1/Ntau^2    G/G_norm     err \n')
    for j,tauT in enumerate(tauTs):
        xdata = numpy.asarray([1/Ntau**2 for k,Ntau in enumerate(Ntaus) if EE_ints[k] is not None])
        if len(xdata) < 2:
            results[j][0] = tauT
            results[j][1] = numpy.nan
            results[j][2] = numpy.nan
            continue
        #tauTindex = interpolation_tauTs.index(tauT)
        ydata = [EE_int[1][j] for EE_int in EE_ints if EE_int is not None]
        edata = [EE_int[2][j]*3 for EE_int in EE_ints if EE_int is not None]
        fitparams, fitparams_err = bootstr.bootstr_from_gauss(fit_sample, data=ydata, data_std_dev=edata, numb_samples = nsamples, sample_size = 1, return_sample=False, args=[xdata] )
        results[j][0] = tauT
        results[j][1] = fitparams[1]
        results[j][2] = fitparams_err[1]
        mycolor = lp.get_color(lpd.tauT_imp[36], j,mintauTindex,18)
        plots.append(ax.errorbar(xdata, ydata, edata, **lp.plotstyle_add_point_single, color=mycolor, zorder=-100+j, label='{0:.3f}'.format(tauT)))
        x = numpy.linspace(0,0.1,100)
        ax.errorbar(0, fitparams[1], fitparams_err[1], **lp.plotstyle_add_point_single, color=mycolor, zorder=1)
        ax.errorbar(x, extrapolation_ansatz(x, *fitparams), color=mycolor, alpha=0.8, fmt='--', lw=0.5, zorder=-100)
        
        numpy.savetxt(outfile, numpy.stack(([tauT for dummy in range(len(xdata)+1)], [0, *xdata], [fitparams[1], *ydata], [fitparams_err[1], *edata]), axis=-1))
        outfile.write('# \n')
numpy.savetxt(lp.inputfolder+"/cont_extr/EE_"+flowradius_str+"_cont.txt", results, header="tauT    G/G_norm    err")

### save continuum extrapolation quality plot for this flow time
ax.legend(handles=plots, title=r'$\tau T$', **lp.legendstyle)
lpd.create_folder(lp.outputfolder+"/cont_extr_quality/")
fig.savefig(lp.outputfolder+"/cont_extr_quality/EE_"+flowradius_str+"_cont_quality.pdf" , **lp.figurestyle)
ax.lines.clear() ; ax.collections.clear() ; plots.clear();

### save plot of continuum extrapolation for this flow time
fig, ax, plots = lp.create_figure(xlims=[0.15,0.51], ylims=[1.4, 4], xlabel=r'$\tau T$', ylabel=r'$ \frac{G_{\tau_F} (\tau)}{G^\mathrm{ norm }_{\tau_F = 0} (\tau)}$', UseTex = False)
ax.set_title(r'$ \sqrt{8\tau_F}T = $'+flowradius_str)
ax.axvline(x=(flowradius/numpy.sqrt(8*0.014)), **lp.verticallinestyle)
#ax.fill_between(results[0], results[1]-results[2], results[1]+results[2], lw = 0.5)
results = numpy.swapaxes(results, 0, 1)
ax.errorbar(results[0], results[1], results[2], color="black", **lp.plotstyle_add_point_single)
lpd.create_folder(lp.outputfolder+"/cont_extr/")
fig.savefig(lp.outputfolder+"/cont_extr/EE_"+flowradius_str+"_cont.pdf" , **lp.figurestyle)
