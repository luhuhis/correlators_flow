#!/usr/local/bin/python
import numpy
import sys
import lib_plot as lp
import lib_process_data as lpd 
import matplotlib

#import _3_spline_interpolate
#import _4_continuum_extr
#import _5_flowtime_extr 

impstr = ""
#impstr = "_unimp"

inputfolder = "../../data_merged/quenched/final_corr_samples"+impstr+"/"
nsamples=10000
EE_final_samples = numpy.ndarray((nsamples,18))
EE_final_samples[:] = numpy.nan

valid_samples = 0
for i in range(nsamples):
    try:
        tmp = numpy.loadtxt(inputfolder+"EE_final_"+str(i)+".txt", unpack=True)
        EE_final_samples[i] = tmp[1]
        valid_samples += 1
    except:
        continue

print(valid_samples)
xdata = lpd.tauT_imp[36]

percentile = 68
EE_final_median = numpy.nanmedian(EE_final_samples, axis=0)
EE_final_err_up = numpy.nanpercentile(EE_final_samples, percentile, axis=0) - EE_final_median
EE_final_err_down = numpy.fabs(numpy.nanpercentile(EE_final_samples, 100-percentile, axis=0) - EE_final_median)
EE_final_err_up_all = numpy.nanpercentile(EE_final_samples, 100, axis=0) - EE_final_median
EE_final_err_down_all = numpy.fabs(numpy.nanpercentile(EE_final_samples, 0, axis=0) - EE_final_median)

EE_final_mean = numpy.nanmean(EE_final_samples,axis=0)
EE_final_stddev = numpy.nanstd(EE_final_samples,axis=0)


numpy.savetxt("../../data_merged/quenched/EE_final_btstrp"+impstr+".txt", numpy.stack((xdata,EE_final_mean,EE_final_stddev) , axis=-1))

EE_cont_2015 = numpy.loadtxt("../../data_merged/quenched/multi-level_2015/cont_thomas.dat")


if impstr == "_unimp":
    ylabel = r'$ \frac{G(\tau)}{G_^\mathrm{latt}{\tau_F = 0}^\mathrm{ pert,cont } (\tau)}$'
else:
    ylabel = r'$ \frac{G(\tau)}{G_{\tau_F = 0}^\mathrm{ pert,latt } (\tau)}$'
fig, ax, plots = lp.create_figure(xlims=[0.2,0.51], ylims=[2, 4], UseTex = False, xlabel=r'$\tau T$', ylabel=ylabel) 
lp.plotstyle_add_point.update(dict(fmt='D-'))
ax.errorbar(xdata, EE_final_median, [EE_final_err_up_all, EE_final_err_down_all], color='black', label="median, 100th percentile", **lp.plotstyle_add_point)
#ax.errorbar(xdata, EE_final_median, [EE_final_err_up, EE_final_err_down], color='red', label="median, "+str(percentile)+"th percentile", **lp.plotstyle_add_point_single)
ax.errorbar(xdata, EE_final_mean, EE_final_stddev, color='blue', label="mean, stddev", **lp.plotstyle_add_point_single)
plots.append(ax.errorbar(EE_cont_2015[:,0], EE_cont_2015[:,1], EE_cont_2015[:,2], **lp.plotstyle_add_point, label=r'Multi-level method$^1$'+"\n"+r'2015', color="black"))#\beta_g (\mathcal{Z}_\mathrm{pert} -1)/6 = 0.079  #*EE_final_mean[-1]/EE_cont_2015[:,1][-1]
#ax.errorbar(xdata, ydata, std_dev, color='grey', label="mean, std_dev", **lp.plotstyle_add_point)
lp.legendstyle.update(dict(handler_map={matplotlib.container.ErrorbarContainer: matplotlib.legend_handler.HandlerErrorbar(yerr_size=0.3, xerr_size=0.3)}, bbox_to_anchor=None, loc='lower right'))
ax.legend(**lp.legendstyle)
fig.savefig("../../plots/quenched/EE_final_btstrp"+impstr+".pdf", **lp.figurestyle)
