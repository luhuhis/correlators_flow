#!/usr/local/bin/python
import numpy
import matplotlib
import lib_plot as pl

pl.ylabelstyle.update(dict(ha='left', va='top'))
pl.xlabelstyle.update(dict(ha='right', va='bottom'))
fig, ax, plots = pl.create_figure(xlims=[0,0.5], ylims=[1,4], xlabel=r'$\tau T$', ylabel=r'$\displaystyle \frac{G^\mathrm{ cont }(\tau)}{G^{\mathrm{norm } }(\tau)}$')
pl.legendstyle.update(dict(loc='upper right', bbox_to_anchor=(1,0.45), labelspacing=0.4, handler_map={matplotlib.container.ErrorbarContainer: matplotlib.legend_handler.HandlerErrorbar(xerr_size=0.5)} ))
ax.xaxis.set_label_coords(0.99,0.01)
ax.yaxis.set_label_coords(0.01,0.97)

EE_final = numpy.loadtxt(pl.inputfolder+"EE_final_cont_gradient_flow.txt")
EE_cont_2015 = numpy.loadtxt("/home/altenkort/master/work/data_merged/quenched/multi-level_2015/cont_thomas.dat")
EE_cont_2015_new = numpy.loadtxt(pl.inputfolder+"/multi-level_2015/EE_2015_cont.txt")

plots.append(ax.errorbar(pl.tauT_30_ext[pl.start:pl.end], [k[1] for k in EE_final], [k[2] for k in EE_final], **pl.plotstyle_points, label="Gradient flow method", color=matplotlib.cm.gnuplot(0.8)))
plots.append(ax.errorbar(EE_cont_2015_new[:,0], EE_cont_2015_new[:,1], EE_cont_2015_new[:,2], color=matplotlib.cm.gnuplot(0.45), alpha=1, **pl.plotstyle_points, label=r'Multi-level method$^{1,2}$'+"\n"+r'w/ $2016$ renorm. update'))#\beta_g (\mathcal{Z}_\mathrm{pert} -1)/6 \approx 0.138
plots.append(ax.errorbar(EE_cont_2015[:,0], EE_cont_2015[:,1], EE_cont_2015[:,2], **pl.plotstyle_points, label=r'Multi-level method$^1$'+"\n"+r'2015', color="black"))#\beta_g (\mathcal{Z}_\mathrm{pert} -1)/6 = 0.079 

ax.legend(handles=plots, title="", **pl.legendstyle)
fig.savefig(pl.outputfolder+"/EE_aeq0_teq0.pdf", **pl.figurestyle) 
print("saved final corr plot", pl.outputfolder+"/EE_aeq0_teq0.pdf")
 
