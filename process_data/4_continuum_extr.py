#!/usr/local/bin/python
import numpy
import matplotlib
import sys
import itertools
#import lib_plot as lp
#import lib_process_data as lpd
from latqcdtools import fitting
from latqcdtools import bootstr
import lib_plot as lp
import lib_process_data as lpd

try:
    flowindex = sys.argv[1]
    i = int(flowindex) #create alias
    flowradius = lp.flow_radius[i]
    flowradius_str = '{0:.3f}'.format(flowradius)
except IndexError:
    exit("Invalid Arguments. Usage: script.py <flowindex>, e.g. script.py 100")
    
    
def extrapolation_ansatz(x, m, b):
    return m * x + b

def fit_sample(ydata, xdata):
    fitter = fitting.Fitter(func=extrapolation_ansatz, xdata=xdata, ydata=ydata, always_return = True)
    fitparams, fitparams_err, chi_dof= fitter.do_fit(start_params = [80,3]) 
    return fitparams

#load interpolations
EE_ints = []
for conftype in ("s064t16_b0687361", "s080t20_b0703500", "s096t24_b0719200", "s120t30_b0739400", "s144t36_b0754400"):
    try:
        EE_ints.append(numpy.loadtxt(lp.inputfolder+"/"+conftype+"/interpolations/EE_"+flowradius_str+"_interpolation.txt", unpack=True))
    except:
        EE_ints.append(None)

nsamples = 2
#EE_ints[0] = None
Ntaus = numpy.asarray([16,20,24,30,36])

tauTs = EE_ints[4][0]
n_tauTs = len(tauTs)
results = numpy.zeros((3,n_tauTs))
fitparameters = []

fig, ax, plots = lp.create_figure(xlims=[-0.0001,0.004], ylims=[1.4, 4], xlabel=r'$N_\tau^{-2}$', ylabel=r'$ \frac{G_{\tau_F} (\tau)}{G^\mathrm{ norm }_{\tau_F = 0} (\tau)}$', UseTex = False)
ax.set_title(r'$ \sqrt{8\tau_F}T = $'+flowradius_str)
lp.legendstyle.update(dict(labelspacing=0.25, handletextpad=0, borderpad=0, handler_map={matplotlib.container.ErrorbarContainer: matplotlib.legend_handler.HandlerErrorbar(xerr_size=0.4)} ))
lpd.create_folder(lp.outputfolder+"/cont_extr_quality/")

mintauTindex = lpd.tauT_imp[36].index(numpy.min(numpy.intersect1d(tauTs, lpd.tauT_imp[36])))
extrapolation_results = []

for j,tauT in enumerate(tauTs):
    #if tauT in lpd.tauT_imp[36]:
        xdata = numpy.asarray([Ntau for k,Ntau in enumerate(Ntaus) if EE_ints[k] is not None])
        xdata = 1/xdata**2
        if len(xdata) < 2:
            results[0][j] = tauT
            results[1][j] = numpy.nan
            results[2][j] = numpy.nan
            continue
        ydata = [EE_int[1][j] for EE_int in EE_ints if EE_int is not None]
        edata = [EE_int[2][j] for EE_int in EE_ints if EE_int is not None]
        
        fitparams, fitparams_err = bootstr.bootstr_from_gauss(fit_sample, data=ydata, data_std_dev=edata, numb_samples = nsamples, sample_size = 1, return_sample=False, args=[xdata] )
        results[0][j] = tauT
        results[1][j] = fitparams[1]
        results[2][j] = fitparams_err[1]
        if tauT in lpd.tauT_imp[36]:
            tauTindex = lpd.tauT_imp[36].index(tauT)
            mycolor = lp.get_color(lpd.tauT_imp[36], tauTindex,mintauTindex,18)
            plots.append(ax.errorbar(xdata, ydata, edata, **lp.plotstyle_add_point_single, color=mycolor, zorder=-100+tauTindex, label='{0:.3f}'.format(tauT)))
            x = numpy.linspace(0,0.1,100)
            ax.errorbar(0, fitparams[1], fitparams_err[1], **lp.plotstyle_add_point_single, color=mycolor, zorder=1)
            ax.errorbar(x, extrapolation_ansatz(x, *fitparams), color=mycolor, alpha=0.8, fmt='--', lw=0.75, zorder=-100)
            #numpy.savetxt(lp.inputfolder+"/cont_extr_quality/EE_"+flowradius_str+"_cont_quality.txt", numpy.stack(([tauT for dummy in range(len(xdata)+1)], [0, *xdata], [fitparams[1], *ydata], [fitparams_err[1], *edata]), axis=-1))
        
numpy.savetxt(lp.inputfolder+"/cont_extr/EE_"+flowradius_str+"_cont.txt", results)

ax.legend(handles=plots, title=r'$\tau T$', **lp.legendstyle)
lpd.create_folder(lp.outputfolder+"/cont_extr_quality/")
fig.savefig(lp.outputfolder+"/cont_extr_quality/EE_"+flowradius_str+"_cont_quality.pdf" , **lp.figurestyle)
ax.lines.clear() ; ax.collections.clear() ; plots.clear();

fig, ax, plots = lp.create_figure(xlims=[0,0.5], ylims=[1.4, 4], xlabel=r'$\tau T$', ylabel=r'$ \frac{G_{\tau_F} (\tau)}{G^\mathrm{ norm }_{\tau_F = 0} (\tau)}$', UseTex = False)
ax.set_title(r'$ \sqrt{8\tau_F}T = $'+flowradius_str)
ax.axvline(x=(flowradius/numpy.sqrt(8*0.014)), **lp.verticallinestyle)
ax.fill_between(results[0], results[1]-results[2], results[1]+results[2], lw = 0.5)
#ax.errorbar(results[0], results[1], results[2], **lp.plotstyle_add_point_single)
lpd.create_folder(lp.outputfolder+"/cont_extr/")
fig.savefig(lp.outputfolder+"/cont_extr/EE_"+flowradius_str+"_cont.pdf" , **lp.figurestyle)
