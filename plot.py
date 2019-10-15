#!/usr/bin/python3
import math
import numpy
from matplotlib import pyplot as plt
import matplotlib.cm as cm
import matplotlib.patheffects as path_effects
from os import listdir
import sys
from latqcdtools import scales_quenched
from latqcdtools import fitting as fit

try:
    qcdtype = sys.argv[1]
    conftype = sys.argv[2]
    beta = sys.argv[3]
    ns = int(sys.argv[4])
    nt = int(sys.argv[5])
except IndexError:
    exit("Invalid Arguments: plot.py <qcdtype> <conftype> <beta> <ns> <nt>")

print(qcdtype, conftype, beta, ns, nt)
outputfolder="../plots/"+qcdtype+"/"+conftype+"/"
inputfolder="../data_merged/"+qcdtype+"/"+conftype+"/"
outputfoldermisc="../plots/"+qcdtype+"/"

#---------
#load data
#---------
flow_radius = numpy.loadtxt(inputfolder+"flowradius_"+conftype+".dat")
EE        = numpy.loadtxt(inputfolder+"EE_"+conftype+".dat")
EE_backup = numpy.loadtxt(inputfolder+"EE_"+conftype+".dat")
EE_err = numpy.loadtxt(inputfolder+"EE_err_"+conftype+".dat")
EE_err_backup = numpy.loadtxt(inputfolder+"EE_err_"+conftype+".dat")
#EE     = EE_backup     = numpy.loadtxt(inputfolder+"EE_flow_"+conftype+".dat")
#EE_err = EE_err_backup = numpy.loadtxt(inputfolder+"EE_flow_err_"+conftype+".dat")
n_datafiles = str(int(numpy.loadtxt(inputfolder+"n_datafiles_"+conftype+".dat")))
tauT = list(numpy.loadtxt(inputfolder+"tauT_imp"+conftype+".dat"))
#EE_ML_2015 = numpy.loadtxt(inputfolder+"EE_ML_2015_"+conftype+".dat")[:int(nt/2),1]
#EE_ML_2015_err = numpy.loadtxt(inputfolder+"EE_ML_2015_"+conftype+".dat")[:int(nt/2),2]
#EE_ML_2019 = numpy.loadtxt(inputfolder+"EE_ML_2019_"+conftype+".dat")[:int(nt/2),3]
#EE_ML_2019_err = numpy.loadtxt(inputfolder+"EE_ML_2019_"+conftype+".dat")[:int(nt/2),4]
EE_cont = []
for flowradius in flow_radius:
    for datafile in listdir("../data_merged/quenched/continuum_limit/"):
            if str(flowradius)+"_cont" in datafile: 
                path = "../data_merged/quenched/continuum_limit/"+datafile
                tmp = numpy.loadtxt(path) 
                EE_cont.append(tmp)
                

#def free_corr( tauT ):
    #norm = math.pi**2 * ( math.cos(math.pi*tauT)**2 / math.sin(math.pi*tauT)**4 + 1/(3*math.sin(math.pi*tauT)**2))
    #return norm
    
#def normalize_EE( data, data_err, tauT_imp, beta ):
    #for i in range(0, int(nt/2)):
        #data[i] = data[i] / free_corr(tauT_imp[i]) * nt**4 * (1+0.079*6/beta) / (-6) 
        #data_err[i] = data_err[i] / free_corr(tauT_imp[i]) * nt**4 * (1+0.079*6/beta) /6 

#normalize_EE(EE_ML_2019, EE_ML_2019_err, tauT, 6.87361)

#---------------------------------------
#perturbative flow limits from Darmstadt
#---------------------------------------
flow_limit=numpy.zeros(len(tauT))
interpolation_value=numpy.zeros(len(tauT))
for i in range(0,len(tauT)):
    flow_limit[i] = tauT[i]*numpy.sqrt(8*0.014)
    index = (numpy.abs(flow_radius[:] - flow_limit[i])).argmin()
    offset = 1 if flow_limit[i]-flow_radius[index] > 0 else -1
    offset2 = -1 if offset == -1 else 0
    interpolation_value[i] = EE[index+offset2,i]-abs(EE[index,i]-EE[index+offset,i])/abs(flow_radius[index]-flow_radius[index+offset])*abs(flow_limit[i]-flow_radius[index+offset2])
flow_limit_lo=numpy.zeros(len(tauT))
interpolation_value_lo=numpy.zeros(len(tauT))
for i in range(0,len(tauT)):
    flow_limit_lo[i] = tauT[i]*numpy.sqrt(8*0.014)/3
    index = (numpy.abs(flow_radius[:] - flow_limit_lo[i])).argmin()
    offset = 1 if flow_limit_lo[i]-flow_radius[index] > 0 else -1
    offset2 = -1 if offset == -1 else 0
    interpolation_value_lo[i] = EE[index+offset2,i]-abs(EE[index,i]-EE[index+offset,i])/abs(flow_radius[index]-flow_radius[index+offset])*abs(flow_limit_lo[i]-flow_radius[index+offset2])

for i in range(0,len(tauT)):
    tauT[i] = round(tauT[i],4)

def skip_large_errors(EE, EE_err, boundary):
    for i in range(0,EE.shape[0]):
        for j in range(0,EE.shape[1]):
            if EE_err[i,j] > boundary :
                EE[i,j] = None

def get_color(myarray, i, start, end, offset=1.0):
    return cm.gnuplot((myarray[i]-myarray[start])/(myarray[end-1]-myarray[start]) * offset)

#-------------------
#set plot parameters
#-------------------
plots = []

mylinewidth = 0.5
mymarkeredgewidth = 0.5

plotstyle          = dict(fmt='D-', linewidth=mylinewidth, mew=0, markersize=2, capsize=3)
plotstyle_add_point= dict(fmt='D-', linewidth=mylinewidth, mew=0.25, mec='grey', markersize=2, capsize=3)
plotstyle_points   = dict(fmt='D-', linewidth=mylinewidth, markersize=0.5)
plotstyle_2015     = dict(fmt='D-', linewidth=1, color='#7CFC00', markersize=2, capsize=5, mew=0.5)
#plotstyle_2019     = dict(fmt='D-', linewidth=1, color='#458e00', markersize=2, capsize=5, mew=0.5)
flowlimitplotstyle = dict(fmt='|', mew=0.7, markersize=15)
legendstyle        = dict(bbox_to_anchor=(1,1), loc="upper left", frameon=True, framealpha=0, title="", labelspacing=0.1, borderpad=0.1, handletextpad=0.4, prop={'size': 9})
figurestyle        = dict(bbox_inches="tight",  pad_inches=0)#, pad_inches=0.05)
labelboxstyle      = dict(boxstyle="square", fc="w", ec='none', alpha=0.7, pad=0.15, zorder=999999)
xlabelstyle        = dict(horizontalalignment='right', verticalalignment='bottom',  bbox=labelboxstyle, zorder=999998)
ylabelstyle        = dict(horizontalalignment='left', verticalalignment='top', rotation=0, bbox=labelboxstyle, zorder=999999)
titlestyle         = dict(zorder=999999)#bbox=labelboxstyle, #verticalalignment='bottom', horizontalalignment='center', x=0.5, y=0.9375, 

plt.rc('text', usetex=True)
plt.rc('text.latex')
plt.rc('font', family='serif', size=10)
aspect=1/2

fig = plt.figure() 
ax = fig.add_subplot(1,1,1) 
#if qcdtype == "quenched":
    #ax.set_title(r'$'+str(ns)+'^3 \\times $'+str(nt), **titlestyle)
#if conftype == "s064t16_b07054_m00129_m0348":
    #ax.set_title(r'$T\approx 0.97\, T_C$, $'+srt(ns)+'^3 \\times $'+str(nt)+', $\#_\mathrm{conf.}=\,$'+n_datafiles+'', **titlestyle)#x=0.975, y=0.95, horizontalalignment='right', 
#if conftype == "s064t16_b07188_m00113_m0306":
    #ax.set_title(r'$T\approx 1.07\, T_C$, $'+str(ns)+'^3 \\times $'+str(nt)+', $\#_\mathrm{conf.}=\,$'+n_datafiles+'', **titlestyle)#x=0.975, y=0.95, horizontalalignment='right', 
    
ax.xaxis.set_label_coords(0.99,0.01)
ax.yaxis.set_label_coords(0.015,0.97)
ax.set_ylabel(r'$\displaystyle \frac{G_{r_F}(\tau T)}{G_{\mathrm{norm}}(\tau T)}$', 
              bbox=labelboxstyle, horizontalalignment='left', verticalalignment='top', rotation=0, zorder=999999)

#-----------------------------------
#PLOT: x-axis tauT, flowtimes legend 
#-----------------------------------
#skip_large_errors(EE, EE_err, EE_backup, EE_err_backup, 10000)
start=0; end=int(nt/2); 
flowstart=0; 
ax.set_xlim([0,0.5]); 
if qcdtype == "quenched":
    flowend=16
    ax.set_ylim([1,4])
if qcdtype == "hisq":
    flowstart=5
    flowend=21
    ax.set_ylim([0,12])
    
ax.set_aspect(1.0/ax.get_data_ratio()*aspect)
ax.set_xlabel(r'$\tau T$', horizontalalignment='right', verticalalignment='bottom', bbox=labelboxstyle, zorder=999999)
legendstyle["title"]=r"$\displaystyle r_F$"

#for i in range(flowstart,flowend):
    #plots.append(ax.errorbar(tauT[start:end], EE[i,start:end], EE_err[i,start:end], color=get_color(flow_radius, i, flowstart, flowend), **plotstyle))
    #for cap in plots[i-flowstart][1]: cap.set_markeredgewidth(mymarkeredgewidth)
##plots.append(ax.plot([],marker="", ls="")[0])
##plots.append(ax.errorbar(tauT, EE_ML_2015, EE_ML_2015_err, **plotstyle_2015))
##plots.append(ax.errorbar(tauT, EE_ML_2019, EE_ML_2019_err, **plotstyle_2019))
#leg = ax.legend(handles=plots, labels=[str(j) for j in flow_radius[flowstart:flowend]]+list(("multi-level","2015","2019")), **legendstyle)
##leg._legend_handle_box.get_children()[0].get_children()[flowend].get_children()[0].set_width(0)
#fig.savefig(outputfolder+"/"+conftype+"_EE.pdf", **figurestyle) 
#ax.lines.clear() ; ax.collections.clear() ; plots.clear()

for i in range(flowstart,flowend):
    plots.append(ax.fill_between(list(tauT), EE[i,:]-EE_err[i,:], EE[i,:]+EE_err[i,:], facecolor=get_color(flow_radius, i, flowstart, flowend), linewidth=mylinewidth))
    ax.errorbar(list(tauT), EE[i,:], color=get_color(flow_radius, i, flowstart, flowend), **plotstyle_add_point)
#plots.append(ax.plot([],marker="", ls="")[0])
#plots.append(ax.errorbar(tauT, EE_ML_2015, EE_ML_2015_err, **plotstyle_2015))
#plots.append(ax.errorbar(tauT, EE_ML_2019, EE_ML_2019_err, **plotstyle_2019))
leg = ax.legend(handles=plots, labels=['{0:.3f}'.format(j) for j in flow_radius[flowstart:flowend]]+list(("multi-level","2015","2019")), **legendstyle)
#leg._legend_handle_box.get_children()[0].get_children()[flowend].get_children()[0].set_width(0)
fig.savefig(outputfolder+"/"+conftype+"_EE_fill.pdf", **figurestyle) 
ax.lines.clear() ; ax.collections.clear() ; plots.clear()


#-----------------------------------
#PLOT: CONTINUUM, x-axis tauT, flowtimes legend 
#-----------------------------------
#for i in range(flowstart,flowend):
    #plots.append(ax.fill_between(EE_cont[i][:,0], EE_cont[i][:,1]-EE_cont[i][:,2], EE_cont[i][:,1]+EE_cont[i][:,2],  facecolor=get_color(flow_radius, i, flowstart, flowend), linewidth=mylinewidth))
#plots.append(ax.plot([],marker="", ls="")[0])
#plots.append(ax.errorbar(tauT, EE_ML_2015, EE_ML_2015_err, **plotstyle_2015))
#leg = ax.legend(handles=plots, labels=[str(j) for j in flow_radius[flowstart:flowend]]+list(("multi-level","2015")), **legendstyle)
#leg._legend_handle_box.get_children()[0].get_children()[flowend].get_children()[0].set_width(0)
#fig.savefig(outputfoldermisc+"/EE_cont_fill.pdf", **figurestyle) 
#ax.lines.clear() ; ax.collections.clear() ; plots.clear()


#PLOT: 


#------------------------------------
#PLOT: x-axis flowtimes, tautT legend 
#------------------------------------
#skip_large_errors(EE, EE_err, EE_backup, EE_err_backup, 10000)
flowstart=0; flowend=24
ax.set_xlim([0,0.16]); 
if qcdtype == "quenched":
    ax.set_ylim([0.7,3.8])
if qcdtype == "hisq":
    ax.set_ylim([0,12])
ax.set_aspect(1.0/ax.get_data_ratio()*aspect)
ax.set_xlabel(r'$\displaystyle r_F $')
ax.set_ylabel(r'$\displaystyle \frac{G_{\tau T}(r_F)}{G_{\tau T,\mathrm{norm}}}$', 
              bbox=labelboxstyle, horizontalalignment='left', verticalalignment='top', rotation=0, zorder=999999)
legendstyle["title"]=r"$\tau T$"

#for i in range(start,end):
    #plots.append(ax.errorbar(flow_radius[flowstart:flowend], EE[flowstart:flowend,i], EE_err[flowstart:flowend,i], color=get_color(tauT, i, start, end), zorder=(-i), **plotstyle))
    #for cap in plots[i-start][1]: cap.set_markeredgewidth(mymarkeredgewidth)
    #ax.errorbar(flow_limit[i], interpolation_value[i], color=get_color(tauT, i, start, end), **flowlimitplotstyle)
#ax.legend(handles=plots, labels=[str(j) for j in tauT[start:end]], **legendstyle)
#fig.savefig(outputfolder+"/"+conftype+"_EE_flowlim.pdf", **figurestyle) 
#ax.lines.clear() ; ax.collections.clear() ; plots.clear()

for i in range(start,end):
    plots.append(ax.fill_between(flow_radius[flowstart:flowend], EE[flowstart:flowend,i]-EE_err[flowstart:flowend,i], 
                                 EE[flowstart:flowend,i]+EE_err[flowstart:flowend,i], facecolor=get_color(tauT, i, start, end), zorder=(-(20+i)), linewidth=mylinewidth))
    ax.errorbar(flow_radius[flowstart:flowend], EE[flowstart:flowend,i], color=get_color(tauT, i, start, end), zorder=(-i), **plotstyle_add_point)
    ax.errorbar(flow_limit[i], interpolation_value[i], color=get_color(tauT, i, start, end), **flowlimitplotstyle)
ax.legend(handles=plots, labels=['{0:.3f}'.format(j) for j in tauT[start:end]], **legendstyle)
fig.savefig(outputfolder+"/"+conftype+"_EE_flowlim_fill.pdf", **figurestyle) 
ax.lines.clear() ; ax.collections.clear() ; plots.clear()


ax.set_ylim([0.7,4.2])
ax.set_aspect(1.0/ax.get_data_ratio()*aspect)
skip_large_errors(EE, EE_err, 0.1)
zero_flow_extr = []
for i in range(start,end):
    plots.append(ax.fill_between(flow_radius[flowstart:flowend], EE[flowstart:flowend,i]-EE_err[flowstart:flowend,i], 
                                 EE[flowstart:flowend,i]+EE_err[flowstart:flowend,i], facecolor=get_color(tauT, i, start, end), zorder=(-i), linewidth=mylinewidth))
    if i > 2:
        func = lambda x,a,b:a*x+b
        func_err = lambda x,da,db:numpy.sqrt((da*x)**2+db**2)
        fitter = fit.Fitter(func, xdata = flow_radius, ydata = EE_backup[:,i], edata = EE_err_backup[:,i])
        res, res_err, chi_dof = fitter.do_fit(start_params = [-0.1, 3], xmin=flow_limit[i]/3, xmax=flow_limit[i])
        zero_flow_extr.append([res[1], res_err[1]])
        x = numpy.linspace(0,flow_limit[i],1000)
        ax.fill_between(x, func(x, *res)-func_err(x, *res_err), func(x, *res)+func_err(x, *res_err), facecolor=get_color(tauT, i, start, end), alpha=0.4)
    ax.errorbar(flow_limit_lo[i], interpolation_value_lo[i], color=get_color(tauT, i, start, end), **flowlimitplotstyle)
    ax.errorbar(flow_limit[i], interpolation_value[i], color=get_color(tauT, i, start, end), **flowlimitplotstyle)
numpy.savetxt("../data_merged/"+qcdtype+"/"+conftype+"/EE_"+conftype+"_zero_flow_extr.dat", zero_flow_extr)
ax.legend(handles=plots, labels=[str(j) for j in tauT[start:end]], **legendstyle)
fig.savefig(outputfolder+"/"+conftype+"_EE_flowlim_fill_w_0flow_extr.pdf", **figurestyle) 
ax.lines.clear() ; ax.collections.clear() ; plots.clear()

#--------------------
#PLOT: relative error
#--------------------
#skip_large_errors(EE, EE_err, EE_backup, EE_err_backup, 100000)
plt.rc('font', family='serif', size=16)
fig = plt.figure() 
ax = fig.add_subplot(1,1,1) 
ax.xaxis.set_label_coords(0.975,0.025)
ax.yaxis.set_label_coords(0.025,0.975)
ax.set_aspect('auto')
mycols=1
if ns == 120:
    mycols=2
flowstart=0; flowend=21
ax.set_yscale('log')
ax.set_xlim([-0.0025,flow_radius[flowend-1]+0.0025]); #ax.autoscale(enable=True, axis='y')
ax.set_ylim([0.01,1000])
#ax.set_title(r'$'+str(ns)+'^3 \\times $'+str(nt), **titlestyle)
ax.set_xlabel(r'$\displaystyle r_F $', **xlabelstyle)
ax.set_ylabel(r'$\displaystyle \frac{\delta G_{\tau T}(r_F)}{\left| G_{\tau T}(r_F)\right|}$ in \%', **ylabelstyle)
for i in range(start,end):
    for j in range(flowstart, flowend):
        if flow_radius[j] > flow_limit[i]:
            EE_backup[j,i] = None
    plots.append(ax.errorbar(flow_radius[flowstart:flowend],  abs(EE_err_backup[flowstart:flowend,i]/EE_backup[flowstart:flowend,i])*100, color=get_color(tauT, i, start, end), zorder=(-i), fmt='D-', linewidth=mylinewidth, mew=0.25, mec='grey', markersize=3, capsize=3))
    for cap in plots[i][1]: cap.set_markeredgewidth(mymarkeredgewidth)
ax.legend(handles=plots, labels=['{0:.3f}'.format(j) for j in tauT[start:end]], bbox_to_anchor=(1,1), ncol=mycols, columnspacing=0.5, title=r'$\tau T$', loc="upper right", frameon=True, framealpha=0, labelspacing=0.1, borderpad=0.1, handletextpad=0.4, handlelength=1.5)
fig.savefig(outputfolder+"/"+conftype+"_EE_relerr.pdf", **figurestyle) 

print("done with", conftype)
