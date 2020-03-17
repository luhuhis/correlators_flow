#!/usr/local/bin/python
import numpy
import matplotlib
import plotlib as pl
from latqcdtools import fitting as fit
#from latqcdtools import bootstr

fig, ax, plots = pl.create_figure(xlims=[-0.00002,0.0016], ylims=[2.7,3.9], xlabel=r'$\tau_F T^2$', ylabel=r'$\displaystyle \frac{G_{\tau}^\mathrm{ cont } (\tau_F)}{G_{\tau,\tau_F=0}^\mathrm{norm } } $')
pl.legendstyle.update(dict(loc='center right', bbox_to_anchor=(1,0.5), handler_map={matplotlib.container.ErrorbarContainer: matplotlib.legend_handler.HandlerErrorbar(xerr_size=0.4)}, handlelength=1.5 ))

flowstart=10
flowend=21
start=2
end=14

#second x-axis for flow radius
ax2=ax.twiny()
new_tick_locations = numpy.array([0,0.05**2/8,0.06**2/8,0.07**2/8,0.08**2/8,0.09**2/8,0.1**2/8])
def tick_function(X):
    V = numpy.sqrt(X*8)
    return ["%.2f" % z for z in V]

ax2.set_xlim(ax.get_xlim())
ax2.set_xticks(new_tick_locations)
ax2.set_xticklabels(tick_function(new_tick_locations))
ax2.set_xlabel(r'$\sqrt{8\tau_F}T$', horizontalalignment='right', verticalalignment='top', bbox=pl.labelboxstyle, zorder=999999)
ax2.xaxis.set_label_coords(0.99,0.98)

EE_cont_arr= []
EE_cont_arr_backup= []

for i in range(flowstart,flowend):
    EE_cont_arr.append(numpy.loadtxt(pl.inputfolder+"/continuum_limit/EE_"+'{0:.4f}'.format(pl.flow_radius[i])+"_cont.txt", max_rows=pl.n_skiprows))
    EE_cont_arr_backup.append(numpy.loadtxt(pl.inputfolder+"/continuum_limit/EE_"+'{0:.4f}'.format(pl.flow_radius[i])+"_cont.txt", max_rows=pl.n_skiprows))
    for j in range(0,len(EE_cont_arr[i-flowstart])):
        if EE_cont_arr[i-flowstart][j,0] < pl.flow_radius[i]/numpy.sqrt(8*0.014)+pl.offset:
            EE_cont_arr[i-flowstart][j,:] = None
#tauT_30_ext=numpy.loadtxt(inputfolder+"/continuum_limit/EE_0.0500_cont.txt", skiprows=17, max_rows=1000)
#tauT_30_ext=tauT_30_ext[:,0]
zero_flow_extr = []

#t->0 extrapolations

#bootstrap method
#def flow_fit(_ydata, _xdata, _xmax):
    #func = lambda x,a,b:a*x+b
    #fitter = fit.Fitter(func, xdata = _xdata, ydata = _ydata)
    #fit_params, fit_params_err, chi_dof = fitter.do_fit(start_params = [-0.1, 3], xmin=0.049, xmax=_xmax)
    #return [fit_params, fit_params_err]

#for i in range(start,end):
    #plots.append(ax.errorbar(flow_radius[flowstart:flowend], [j[i,1] for j in EE_cont_arr], [j[i,2] for j in EE_cont_arr], **plotstyle_points, zorder=-end+start+i, color=get_color(tauT_30_ext, i, start, end)))
    #res, res_err = bootstr.bootstr_from_gauss(func=flow_fit, data=[j[i,1] for j in EE_cont_arr_backup], data_std_dev=[j[i,2] for j in EE_cont_arr_backup],
                                                       #args={'_xdata':flow_radius[flowstart:flowend], '_xmax':(tauT_30_ext[i]-0.02)*numpy.sqrt(8*0.014)}, numb_samples=100)
    #zero_flow_extr.append([tauT_30_ext[i], res[0][1], res_err[0][1]])
    #x = numpy.linspace(0,0.1,1000)
    #func = lambda x,a,b:a*x+b
    #func_err = lambda x,da,db:numpy.sqrt((da*x)**2+db**2)
    ##ax.fill_between(x, func(x, *res)-func_err(x, *res_err), func(x, *res)+func_err(x, *res_err), facecolor=get_color(tauT_30_ext, i, start, end), alpha=0.6, zorder=(-2*end+start+i))
    #ax.errorbar(x, func(x, *res[0]), color=get_color(tauT_30_ext, i, start, end), alpha=0.6, fmt='-', zorder=(-2*end+start+i))
    #ax.errorbar(0, res[0][1], res_err[0][1], color=get_color(tauT_30_ext, i, start, end), alpha=0.6, zorder=(-2*end+start+i), **plotstyle_points)
    #print("done with ", i-start, "of", end-start) 

#simple method
for i in range(start,end):
    plots.append(ax.errorbar(pl.flow_radius[flowstart:flowend]**2/8, [j[i,1] for j in EE_cont_arr], [j[i,2] for j in EE_cont_arr], **pl.plotstyle_points, zorder=-i, color=pl.get_color(pl.tauT_30_ext, i, start, end)))
    func = lambda x,a,b:a*x+b
    #func_err = lambda x,da,db:numpy.sqrt((da*x)**2+db**2)
    fitter = fit.Fitter(func, xdata = pl.flow_radius[flowstart:flowend]**2/8, ydata = [j[i,1] for j in EE_cont_arr_backup], edata = [j[i,2] for j in EE_cont_arr_backup])
    res, res_err, chi_dof = fitter.do_fit(start_params = [-3, 3], xmin=0.049**2/8, xmax=((pl.tauT_30_ext[i]-0.02)*numpy.sqrt(8*0.014))**2/8)
    zero_flow_extr.append([pl.tauT_30_ext[i], res[1], res_err[1]])
    x = numpy.linspace(0,0.1,1000)
    #ax.errorbar((tauT_30_ext[i]-0.02)*numpy.sqrt(8*0.014), 2.3+tauT_30_ext[i], color=get_color(tauT_30_ext, i, start, end), **plotstyle_points)
    #ax.fill_between(x, func(x, *res)-func_err(x, *res_err), func(x, *res)+func_err(x, *res_err), facecolor=get_color(tauT_30_ext, i, start, end), alpha=0.6, zorder=-2*i)
    ax.errorbar(x, func(x, *res), color=pl.get_color(pl.tauT_30_ext, i, start, end), alpha=0.6, fmt='--', lw=0.75, zorder=(-2*end+start+i))
    ax.errorbar(0, res[1], res_err[1], color=pl.get_color(pl.tauT_30_ext, i, start, end), alpha=0.6, zorder=(-2*end+start+i), **pl.plotstyle_points)
    
ax.legend(handles=plots, labels=['{0:.3f}'.format(i) for i in pl.tauT_30_ext[start:end]], title=r'$\tau T$', **pl.legendstyle)
fig.savefig(pl.outputfolder+"/EE_aeq0_t_extr.pdf", **pl.figurestyle) 

#save extrapolations to file
finalresult=[[k[0],k[1],k[2]] for k in zero_flow_extr]
numpy.savetxt(pl.inputfolder+"EE_final_cont_gradient_flow.txt", finalresult) 
print("saved flowtime extr plot", pl.inputfolder+"EE_final_cont_gradient_flow.txt")

