#!/usr/local/bin/python
import numpy
import matplotlib
import lib_plot as pl
from latqcdtools import fitting as fit
from pathlib import Path

Path(pl.outputfolder+"cont_lim_quality/").mkdir(parents=True, exist_ok=True)

pl.ylabelstyle.update(dict(ha='left', va='top'))
pl.xlabelstyle.update(dict(ha='right', va='bottom'))
fig, ax, plots = pl.create_figure(xlims=[-0.0001,1/15**2+0.001], ylims=[2.9,3.8], xlabel=r'$N_\tau^{-2}$', ylabel=r'$\displaystyle \frac{G_{\tau,\tau_F} }{G^{\mathrm{free}}_{\tau,\tau_F=0} }$')

ax.xaxis.set_label_coords(0.99,0.01)
ax.yaxis.set_label_coords(0.01,0.97)
pl.legendstyle.update(dict(loc='center right', labelspacing=0.25, handletextpad=0, borderpad=0, handler_map={matplotlib.container.ErrorbarContainer: matplotlib.legend_handler.HandlerErrorbar(xerr_size=0.4)} ))
pl.titlestyle.update(dict(y=0.95))
pl.plotstyle_points.update(dict(fmt='D', markersize=2))

start=6
end=14
#be careful with the 

flowstart=20
flowend=21

for i in range(flowstart,flowend):
    for tauT_index in range(start,end):
        EE_16 = numpy.loadtxt(pl.inputfolder+"continuum_limit/EE_"+'{0:.4f}'.format(pl.flow_radius[i])+"_interpolation_16.txt", max_rows=14)
        EE_20 = numpy.loadtxt(pl.inputfolder+"continuum_limit/EE_"+'{0:.4f}'.format(pl.flow_radius[i])+"_interpolation_20.txt", max_rows=14)
        EE_24 = numpy.loadtxt(pl.inputfolder+"continuum_limit/EE_"+'{0:.4f}'.format(pl.flow_radius[i])+"_interpolation_24.txt", max_rows=14)
        EE_30 = numpy.loadtxt(pl.inputfolder+"continuum_limit/EE_"+'{0:.4f}'.format(pl.flow_radius[i])+"_interpolation_30.txt", max_rows=14)
        EE_cont = numpy.loadtxt(pl.inputfolder+"continuum_limit/EE_"+'{0:.4f}'.format(pl.flow_radius[i])+"_cont.txt")
        ax.set_title(r'$\sqrt{8\tau_F}T=$'+'{0:.2f}'.format(pl.flow_radius[i]),**pl.titlestyle) #r'$ \tau T = $' + '{0:.2f}'.format(EE_cont[tauT_index,1]) +
        
        data = []
        Nts = [16, 20, 24, 30]
        data.append([0, EE_cont[tauT_index,1], EE_cont[tauT_index,2]])
        data.append([1/16**2, EE_16[tauT_index,1], EE_16[tauT_index,2]])
        data.append([1/20**2, EE_20[tauT_index,1], EE_20[tauT_index,2]])
        data.append([1/24**2, EE_24[tauT_index,1], EE_24[tauT_index,2]])
        data.append([1/30**2, EE_30[tauT_index,1], EE_30[tauT_index,2]])
        
        _xdata = [j[0] for j in data[:]]
        _ydata = [j[1] for j in data[:]]
        _edata = [j[2] for j in data[:]]
        
        plots.append(ax.errorbar(_xdata, _ydata, _edata, **pl.plotstyle_points, color=pl.get_color(EE_cont[:,0], tauT_index,start,end), zorder=-tauT_index, label='{0:.3f}'.format(EE_cont[tauT_index,0])))
        
        #fit
        func = lambda x,a,b:a*x+b
        fitter = fit.Fitter(func, xdata = _xdata, ydata = _ydata, edata = _edata)
        res, res_err, chi_dof = fitter.do_fit(start_params = [0, 3])
        x = numpy.linspace(0,0.1,1000)
        ax.errorbar(x, func(x, *res), color='grey', alpha=0.8, fmt='--', lw=0.75, zorder=0)
    ax.legend(handles=plots, title=r'$\tau T$', **pl.legendstyle)
    fig.savefig(pl.outputfolder+"cont_lim_quality/EE_rF"+'{0:.4f}'.format(pl.flow_radius[i])+"_interpolation.pdf", **pl.figurestyle) 
    print("saved extrapolation quality plot", pl.outputfolder+"continuum_limit/EE_rF"+'{0:.4f}'.format(pl.flow_radius[i])+"_interpolation.pdf")
    ax.lines.clear() ; ax.collections.clear() ; plots.clear();
