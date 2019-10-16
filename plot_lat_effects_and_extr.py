#!/usr/local/bin/python
import math
import numpy
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from matplotlib.legend_handler import HandlerErrorbar
from os import listdir
import sys
from latqcdtools import fitting as fit

def get_color(myarray, i, start, end):
    return cm.gnuplot((myarray[i]-myarray[start])/(myarray[end-1]-myarray[start]))

qcdtype="quenched"
outputfolder="../plots/"+qcdtype+"/"
inputfolder="../data_merged/"+qcdtype+"/"

EE_16 = numpy.loadtxt(inputfolder+"/s064t16_b0687361/EE_s064t16_b0687361.dat")
EE_20 = numpy.loadtxt(inputfolder+"/s080t20_b0703500/EE_s080t20_b0703500.dat")
EE_24 = numpy.loadtxt(inputfolder+"/s096t24_b0719200/EE_s096t24_b0719200.dat")
EE_30 = numpy.loadtxt(inputfolder+"/s120t30_b0739400/EE_s120t30_b0739400.dat")
EE_err_16 = numpy.loadtxt(inputfolder+"/s064t16_b0687361/EE_err_s064t16_b0687361.dat")
EE_err_20 = numpy.loadtxt(inputfolder+"/s080t20_b0703500/EE_err_s080t20_b0703500.dat")
EE_err_24 = numpy.loadtxt(inputfolder+"/s096t24_b0719200/EE_err_s096t24_b0719200.dat")
EE_err_30 = numpy.loadtxt(inputfolder+"/s120t30_b0739400/EE_err_s120t30_b0739400.dat")
EE_cont_2015 = numpy.loadtxt("/home/altenkort/master/work/data_merged/quenched/multi-level_2015/cont_thomas.dat")
EE_cont_2015_new = numpy.loadtxt(inputfolder+"/multi-level_2015/EE_2015_cont.txt")
EE_zero_flow_16 = numpy.loadtxt(inputfolder+"/s064t16_b0687361/EE_s064t16_b0687361_zero_flow_extr.dat")
EE_zero_flow_20 = numpy.loadtxt(inputfolder+"/s080t20_b0703500/EE_s080t20_b0703500_zero_flow_extr.dat")
EE_zero_flow_24 = numpy.loadtxt(inputfolder+"/s096t24_b0719200/EE_s096t24_b0719200_zero_flow_extr.dat")
EE_zero_flow_30 = numpy.loadtxt(inputfolder+"/s120t30_b0739400/EE_s120t30_b0739400_zero_flow_extr.dat")

flow_radius = numpy.loadtxt(inputfolder+"/s064t16_b0687361/flowradius_s064t16_b0687361.dat")

tauT_16 = list(numpy.around(numpy.loadtxt(inputfolder+"/s064t16_b0687361/tauT_imps064t16_b0687361.dat"), 4))
tauT_20 = list(numpy.around(numpy.loadtxt(inputfolder+"/s080t20_b0703500/tauT_imps080t20_b0703500.dat"), 4))
tauT_24 = list(numpy.around(numpy.loadtxt(inputfolder+"/s096t24_b0719200/tauT_imps096t24_b0719200.dat"), 4))
tauT_30 = list(numpy.around(numpy.loadtxt(inputfolder+"/s120t30_b0739400/tauT_imps120t30_b0739400.dat"), 4))
#tauT_30_ext = (*tauT_30,0.5)
tauT_30_ext = (0.229167, 0.250000, 0.270833, 0.291667, 0.312500, 0.333333, 0.354167, 0.375000, 0.395833, 0.416667, 0.437500, 0.458333, 0.479167, 0.500000)

#-------------------
#set plot parameters
#-------------------
figurestyle        = dict(bbox_inches="tight")
plotstyle_points   = dict(fmt='D-', linewidth=1, markersize=4, capsize=2, mew=0.5, fillstyle='none')
labelboxstyle      = dict(boxstyle="square", fc="w", ec='none', alpha=0.7, pad=0.15, zorder=999999)
legendstyle        = dict(loc="center left", frameon=True, framealpha=0.8, edgecolor='none', fancybox=False, facecolor="w", title="", labelspacing=0.1, borderpad=0.1, handletextpad=0.4, prop={'size': 9})#, handlelength=1)#, prop={'size': 9})
figurestyle        = dict(bbox_inches="tight", pad_inches=0)
plotstyle_fill     = dict(linewidth=0.5)
xlabelstyle        = dict(horizontalalignment='right', verticalalignment='bottom',  bbox=labelboxstyle, zorder=999998)
ylabelstyle        = dict(horizontalalignment='left', verticalalignment='top', rotation=0, bbox=labelboxstyle, zorder=999999)
titlestyle         = dict(x=0.5, y=0.9, bbox=labelboxstyle, verticalalignment='top', zorder=999999)


plt.rc('text', usetex=True)
plt.rc('text.latex')
plt.rc('font', family='serif', size=8)
aspect=1/3.8
fig = plt.figure() 
ax = fig.add_subplot(1,1,1) 
ax.xaxis.set_label_coords(0.99,0.01)
ax.yaxis.set_label_coords(0.01,0.97)
#ax.set_ylabel(r'$\displaystyle \frac{G_{EE}(\tau T, \tau_F)}{G_{EE}^{\,\mathrm{free}}(\tau T)}$', bbox=labelboxstyle, horizontalalignment='left', verticalalignment='top', rotation=0, zorder=999999)
#legendstyle["title"]=r"$N_\tau$"

offset=0.02
flowstart=10
flowend=21


#PLOT LATTICE EFFECTS FOR ANIMATION
#for i in range(10,21):
ax.set_xlim([0,0.51])
ax.set_ylim([0.75,4.25])
ax.set_aspect(1.0/ax.get_data_ratio()*aspect)
plots = []
#legendstyle['bbox_to_anchor']=(0,0.4)
#legendstyle['loc']='upper left'
#legendstyle['borderpad']=0.25
ax.set_ylabel(r'$\displaystyle \frac{G_{r_F}(\tau T)}{G_{\mathrm{norm}}(\tau T)}$', **ylabelstyle)
ax.set_xlabel(r'$\tau T$', **xlabelstyle)
ax.get_xaxis().set_visible(False)
for i in range(0,len(flow_radius)):
#for i in (0,10,20,22):
    #ax.set_title(r'quenched $|$ $T\approx 1.5 T_C$ $|$ imp. dist. $|$ $\sqrt{\tilde{\tau_F}} =$ '+'{0:.4f}'.format(flow_radius[i]))
    #ax.set_title(r'$r_F =$ '+'{0:.3f}'.format(flow_radius[i]), x=0.5, y=0.925, bbox=labelboxstyle, zorder=999999)
    ax.set_title(r'$r_F =$ '+'{0:.2f}'.format(flow_radius[i]), **titlestyle)
    #text = ax.text(.5,.95,r'$r_F =$ '+'{0:.4f}'.format(flow_radius[i]), horizontalalignment='center', transform=ax.transAxes, bbox=labelboxstyle, zorder=999999)
    plots.append(ax.errorbar(tauT_16, EE_16[i,:], EE_err_16[i,:], label=r'$N_\tau = 16$', **plotstyle_points, color=cm.gnuplot(0.03), zorder=-3))
    plots.append(ax.errorbar(tauT_20, EE_20[i,:], EE_err_20[i,:], label=r'$N_\tau = 20$', **plotstyle_points, color=cm.gnuplot(0.40), zorder=-2))
    plots.append(ax.errorbar(tauT_24, EE_24[i,:], EE_err_24[i,:], label=r'$N_\tau = 24$', **plotstyle_points, color=cm.gnuplot(0.65), zorder=-1))
    plots.append(ax.errorbar(tauT_30, EE_30[i,:], EE_err_30[i,:], label=r'$N_\tau = 30$', **plotstyle_points, color=cm.gnuplot(0.90), zorder=0))
    try:
        EE_cont = numpy.loadtxt(inputfolder+"/continuum_limit/EE_"+'{0:.4f}'.format(flow_radius[i])+"_cont.txt", skiprows=14)
        for j in range(0,len(EE_cont)):
            if EE_cont[j,0] < flow_radius[i]/numpy.sqrt(8*0.014)+offset:
                EE_cont[j,:] = None
        plots.append(ax.fill_between(EE_cont[:,0], EE_cont[:,1]-EE_cont[:,2], EE_cont[:,1]+EE_cont[:,2], facecolor='grey', zorder=-6))
        ax.legend(handles=plots, labels=[r'$N_\tau = 16$', r'$N_\tau = 20$', r'$N_\tau = 24$', r'$N_\tau = 30$', r'cont.'], loc="upper left", frameon=True, framealpha=0.5, edgecolor='none', fancybox=False, facecolor="w", title="", labelspacing=0.1, borderpad=0.1, handletextpad=0.4, handler_map={type(plots[0]): HandlerErrorbar(xerr_size=0.4)}, bbox_to_anchor=(0,0.65))#, prop={'size': myfontsize})
        ax.errorbar(EE_cont[:,0], EE_cont[:,1], color='grey', lw=1.5, fmt='-', zorder=-5, solid_capstyle='butt')
    except:
        #plots.append(ax.plot([],marker="", ls="", label="s")[0])
        ax.legend(handles=plots, labels=[r'$N_\tau = 16$', r'$N_\tau = 20$', r'$N_\tau = 24$', r'$N_\tau = 30$'], loc="upper left", frameon=True, framealpha=0.5, edgecolor='none', fancybox=False, facecolor="w", title="", labelspacing=0.1, borderpad=0.1, handletextpad=0.4, handler_map={type(plots[0]): HandlerErrorbar(xerr_size=0.4)}, bbox_to_anchor=(0,0.65))#, prop={'size': myfontsize})
        #pass
    if i != 0:
        ax.axvline(x=(flow_radius[i]/numpy.sqrt(8*0.014)+offset), ymin=0, ymax=1, color='grey', alpha=0.8, zorder=-100, dashes=(4,4), lw=0.5)
    if i == 22:
        ax.get_xaxis().set_visible(True)
    fig.savefig(outputfolder+"/single_flow/EE_flow_"+'{0:.4f}'.format(flow_radius[i])+".pdf", **figurestyle) 
    ax.lines.clear() ; ax.collections.clear() ; plots.clear();
    ax.set_title(r'')
legendstyle['borderpad']=0.1

#PLOT all continuum extr. in one plot
plt.rc('font', family='serif', size=10)
fig = plt.figure() 
ax = fig.add_subplot(1,1,1) 
ax.xaxis.set_label_coords(0.99,0.01)
ax.yaxis.set_label_coords(0.01,0.97)
ax.set_ylabel(r'$\displaystyle \frac{G_{r_F,\mathrm{cont}}(\tau T)}{G_{\mathrm{norm}}(\tau T)}$', **ylabelstyle)
ax.set_xlabel(r'$\tau T$', **xlabelstyle)
aspect=1/2
ax.set_ylim([2,4])
ax.set_xlim([0.15,0.5])
ax.set_aspect(1.0/ax.get_data_ratio()*aspect)
legendstyle['loc']='upper left'
legendstyle['bbox_to_anchor']=(1,1)
#ax.set_title(r'$|$ $T\approx 1.50 T_C$ $|$ imp. dist. $|$ continuum extrapolations')
for i in range(0,len(flow_radius)):
    flow_radius[i] = round(flow_radius[i],4)
for i in range(flowstart,flowend):
    EE_cont = numpy.loadtxt(inputfolder+"/continuum_limit/EE_"+'{0:.4f}'.format(flow_radius[i])+"_cont.txt")
    for j in range(0,len(EE_cont)):
        if EE_cont[j,0] < flow_radius[i]/numpy.sqrt(8*0.014)+offset:
            EE_cont[j,:] = None
    plots.append(ax.fill_between(EE_cont[:,0], EE_cont[:,1]-EE_cont[:,2], EE_cont[:,1]+EE_cont[:,2], label='{0:.3f}'.format(flow_radius[i]), facecolor=get_color(flow_radius, i, flowstart, flowend), zorder=i))
legendstyle["title"]=r"$r_F$"
ax.legend(handles=plots, labels=['{0:.3f}'.format(i) for i in flow_radius[flowstart:flowend]], **legendstyle)
fig.savefig(outputfolder+"/EE_aeq0_tneq0.pdf", **figurestyle) 
ax.lines.clear() ; ax.collections.clear() ; plots.clear()


#PLOT CONT LIMITS AXES FLIPPED WITH t->0 EXTRAPOLATIONS
start=1
end=14
legendstyle['prop']={'size': 10}
#ax.set_title(r'pure SU(3) $|$ $T\approx 1.5 T_C$ $|$ imp. dist. $|$ $a = 0$ $|$ $\tau_F\rightarrow 0$ ')
ax.set_ylabel(r'$\displaystyle \frac{G_{\tau T, \mathrm{cont}}(r_F)}{G_{\tau T, \mathrm{norm}}}$')
ax.set_xlabel(r'$r_F$')
ax.set_xlim([0,0.1025])
ax.set_ylim([2.45,4.1])
ax.set_aspect(1.0/ax.get_data_ratio()*aspect)
EE_cont_arr= []
EE_cont_arr_backup= []

for i in range(flowstart,flowend):
    EE_cont_arr.append(numpy.loadtxt(inputfolder+"/continuum_limit/EE_"+'{0:.4f}'.format(flow_radius[i])+"_cont.txt"))#, max_rows=17))
    EE_cont_arr_backup.append(numpy.loadtxt(inputfolder+"/continuum_limit/EE_"+'{0:.4f}'.format(flow_radius[i])+"_cont.txt"))#, max_rows=17))
    for j in range(0,len(EE_cont_arr[i-flowstart])):
        if EE_cont_arr[i-flowstart][j,0] < flow_radius[i]/numpy.sqrt(8*0.014)+offset:
            EE_cont_arr[i-flowstart][j,:] = None
#tauT_30_ext=numpy.loadtxt(inputfolder+"/continuum_limit/EE_0.0500_cont.txt", skiprows=17, max_rows=1000)
#tauT_30_ext=tauT_30_ext[:,0]
zero_flow_extr = []
for i in range(start,end):
    plots.append(ax.errorbar(flow_radius[flowstart:flowend], [j[i,1] for j in EE_cont_arr], [j[i,2] for j in EE_cont_arr], **plotstyle_points, zorder=-i, color=get_color(tauT_30_ext, i, start, end)))
    func = lambda x,a,b:a*x+b
    func_err = lambda x,da,db:numpy.sqrt((da*x)**2+db**2)
    fitter = fit.Fitter(func, xdata = flow_radius[flowstart:flowend], ydata = [j[i,1] for j in EE_cont_arr_backup], edata = [j[i,2] for j in EE_cont_arr_backup])
    res, res_err, chi_dof = fitter.do_fit(start_params = [-0.1, 3], xmin=0.049, xmax=(tauT_30_ext[i]-0.02)*numpy.sqrt(8*0.014))
    zero_flow_extr.append([tauT_30_ext[i], res[1], res_err[1]])
    x = numpy.linspace(0,0.1,1000)
    #ax.errorbar((tauT_30_ext[i]-0.02)*numpy.sqrt(8*0.014), 2.3+tauT_30_ext[i], color=get_color(tauT_30_ext, i, start, end), **plotstyle_points)
    ax.fill_between(x, func(x, *res)-func_err(x, *res_err), func(x, *res)+func_err(x, *res_err), facecolor=get_color(tauT_30_ext, i, start, end), alpha=0.6, zorder=-2*i)
legendstyle["title"]=r"$\tau T$"
legendstyle["loc"]="upper left"
legendstyle["bbox_to_anchor"]=(1,1)
ax.legend(handles=plots, labels=['{0:.3f}'.format(i) for i in tauT_30_ext[start:end]], **legendstyle, handler_map={type(plots[0]): HandlerErrorbar(xerr_size=0.4)}, handlelength=1.5)
fig.savefig(outputfolder+"/EE_aeq0_t_extr.pdf", **figurestyle) 
ax.lines.clear() ; ax.collections.clear() ; plots.clear()


finalresult=[[k[0],k[1],k[2]] for k in zero_flow_extr]
numpy.savetxt(inputfolder+"EE_final_cont_gradient_flow.txt", finalresult)
#PLOT FINAL RESULT
plt.rc('font', family='serif', size=9)
fig = plt.figure() 
ax = fig.add_subplot(1,1,1) 
ax.xaxis.set_label_coords(0.99,0.01)
ax.yaxis.set_label_coords(0.01,0.97)
#ax.set_title(r'pure SU(3) $|$ $T\approx 1.5 T_C$ $|$ imp. dist.')
#ax.set_ylabel(r'$\displaystyle \frac{G_{EE}(\tau T)}{G_{EE}^{\,\mathrm{free}}(\tau T)}$')
ax.set_ylabel(r'$\displaystyle \frac{G_\mathrm{cont}(\tau T)}{G_{\mathrm{norm}}(\tau T)}$', **ylabelstyle)
ax.set_xlabel(r'$\tau T$', **xlabelstyle)
ax.set_xlim([0,0.505])
ax.set_ylim([1,4])
ax.set_aspect(1.0/ax.get_data_ratio()*aspect)
plots.append(ax.errorbar(tauT_30_ext[start:end], [k[1] for k in zero_flow_extr], [k[2] for k in zero_flow_extr], **plotstyle_points, label="Gradient Flow", color=cm.gnuplot(0.8)))
#plots.append(ax.fill_between(EE_cont_2015_new[:,0], EE_cont_2015_new[:,1]-EE_cont_2015_new[:,2], EE_cont_2015_new[:,1]+EE_cont_2015_new[:,2], facecolor='grey', alpha=1, **plotstyle_fill, label="2015 ML \n new cont. extr. \n w/ updated Z"))
plots.append(ax.errorbar(EE_cont_2015_new[:,0], EE_cont_2015_new[:,1], EE_cont_2015_new[:,2], color=cm.gnuplot(0.45), alpha=1, **plotstyle_points, label=r'Multi-Level'+"\n"+r'$\beta_g (\mathcal{Z}_\mathrm{pert} -1)/6 \approx 0.138$'))
#plots.append(ax.errorbar(EE_cont_2015_new[:,0], EE_cont_2015_new[:,1], EE_cont_2015_new[:,2], color=cm.gnuplot(0.45), alpha=1, **plotstyle_points, label=r'Multi-Level'+"\n"+r'2-loop renorm.'))
plots.append(ax.errorbar(EE_cont_2015[:,0], EE_cont_2015[:,1], EE_cont_2015[:,2], **plotstyle_points, label=r'Multi-Level'+"\n"+r'$\beta_g (\mathcal{Z}_\mathrm{pert} -1)/6 = 0.079 $', color="black"))
#plots.append(ax.errorbar(EE_cont_2015[:,0], EE_cont_2015[:,1], EE_cont_2015[:,2], **plotstyle_points, label=r'Multi-Level'+"\n"+r'1-loop renorm.', color="black"))
#plots.append(ax.errorbar(tauT_16[3:], EE_zero_flow_16[:,0], EE_zero_flow_16[:,1], label=r'$N_\tau = 16$'+"\nflow extr.", **plotstyle_points))
#plots.append(ax.errorbar(tauT_20[3:], EE_zero_flow_20[:,0], EE_zero_flow_20[:,1], label=r'$N_\tau = 20$'+"\nflow extr.", **plotstyle_points))
#plots.append(ax.errorbar(tauT_24[3:], EE_zero_flow_24[:,0], EE_zero_flow_24[:,1], label=r'$N_\tau = 24$'+"\nflow extr.", **plotstyle_points))
#plots.append(ax.errorbar(tauT_30[3:], EE_zero_flow_30[:,0], EE_zero_flow_30[:,1], label=r'$N_\tau = 30$'+"\nflow extr.", **plotstyle_points))
ax.legend(handles=plots, loc="upper right", frameon=True, framealpha=0.8, edgecolor='none', fancybox=False, facecolor="w", title="", labelspacing=0.4, borderpad=0.1, handletextpad=0.4, bbox_to_anchor=(1,0.5))
fig.savefig(outputfolder+"/EE_aeq0_teq0.pdf", **figurestyle) 
