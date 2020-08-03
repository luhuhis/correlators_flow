#!/usr/local/bin/python
import numpy
import matplotlib
import sys
import itertools
import lib_plot as lp
import lib_process_data as lpd
from matplotlib import cm

try:
    flowindex = sys.argv[1]
    i = int(flowindex) #create alias
    flowradius = lp.flow_radius[i]
except IndexError:
    exit("Invalid Arguments. Usage: script.py <flowindex>, e.g. script.py 100")

use_imp = True

UseTex = True
if UseTex:
    displaystyle = r'\displaystyle'
    ylabel = r'$ '+displaystyle+r'\frac{G^\mathrm{ latt }_{\tau_F} (\tau)}{G_{\tau_F = 0}^{\substack{ \text{\tiny  norm} \\[-0.4ex] \text{\tiny latt } } } (\tau)} $'
else:
    ylabel = r'G'
    displaystyle = r''
    
lpd.create_folder(lp.outputfolder+"/lattice_effects/")
fig, ax, lpots = lp.create_figure(xlims=[0.15,0.51], ylims=[2,4], xlabel=r'$\tau T$', ylabel=ylabel, xlabelpos=(0.95,0.05), ylabelpos=(0.03,0.98), UseTex = UseTex)

lp.titlestyle.update(dict(y=0.95))
lp.legendstyle.update(dict(loc="lower right", bbox_to_anchor=(1.01,0.1), framealpha=0.5, handlelength=1, handler_map={matplotlib.container.ErrorbarContainer: matplotlib.legend_handler.HandlerErrorbar(xerr_size=0.4)}))

### load finest lattice and interpolations of coarser lattices into one array
EE_data = []
for conftype in ("s064t16_b0687361", "s080t20_b0703500", "s096t24_b0719200", "s120t30_b0739400", "s144t36_b0754400"): 
    ### skip coarsest lattice
    if conftype == "s064t16_b0687361":
        EE_data.append(None)
    #### load finest lattice separately because it isnt interpolated
    #elif conftype == "s144t36_b0754400":
        #EE36 = numpy.asarray(([tauT for tauT in lpd.tauT_imp[36] if tauT > lpd.lower_tauT_limit(flowradius)], 
                     #[val*lpd.improve_corr_factor(j, 36, i, use_imp) for j,val in enumerate(lp.EE_36[i]) if lpd.tauT_imp[36][j] > lpd.lower_tauT_limit(flowradius)], 
                     #[val*lpd.improve_corr_factor(j, 36, i, use_imp) for j,val in enumerate(lp.EE_err_36[i]) if lpd.tauT_imp[36][j] > lpd.lower_tauT_limit(flowradius)]))
        #EE_data.append(EE36)        
    else: 
        EE_data.append(numpy.loadtxt(lp.inputfolder+"/"+conftype+"/interpolations/EE_"+'{0:.3f}'.format(lp.flow_radius[i])+"_interpolation.txt", unpack=True))
### also load continuum extr if possible
try:
    EE_data.append(numpy.loadtxt(lp.inputfolder+"/cont_extr/EE_"+'{0:.3f}'.format(lp.flow_radius[i])+"_cont.txt", unpack=True))
except:
    EE_data.append(None)
    
### some plot settings
#colors = matplotlib.cm.gnuplot([0.05, 0.2, 0.4, 0.75, 0.9, 0])

colors=cm.tab10(numpy.linspace(0,1,10))
colors = ('', *colors)
#for c in iter(color):
    #print(c)

#colors =  ['','green', 'blue', 'red', 'brown', 'black']
labels = [r'$N_\tau = 16$', r'$N_\tau = 20$', r'$N_\tau = 24$', r'$N_\tau = 30$', r'$N_\tau = 36$', r'cont.']
plots = []
ax.set_title(r'$\sqrt{8\tau_F} T=$ '+'{0:.3f}'.format(lp.flow_radius[i]), **lp.titlestyle)
markers = ['', 's', 'o', 'D', 'H', 'h']
zorders = range(-30,0,4)
### loop over all the different data sets and plot settings
for EE_int, tauTs, EE, EE_err, label, color, marker, zorder in itertools.zip_longest(EE_data, lp.tauT, lp.EE, lp.EE_err, labels, colors, markers, zorders):
    lp.plotstyle_add_point.update(dict(fmt=marker))
    if tauTs is not None: #tauTs is None for the continuum limit
        nt = int(len(tauTs)*2)
    if EE_int is not None:
        ### remove nans
        columns = []
        nan_array = numpy.isnan(EE_int[1])
        for l,column in enumerate(EE_int):
            columns.append(column[~nan_array])
        ### plot original Data (EE, EE_err), interpolations (EE_int) and continuum extrapolation
        if tauTs is None: #continuum plot
            plots.append(ax.errorbar(columns[0], columns[1], columns[2], color=color, label=label, **lp.plotstyle_add_point, zorder=zorder+1))
            lp.plotstyle_add_point.update(dict(fmt='-'))
            ax.errorbar(columns[0], columns[1], columns[2], color=color, **lp.plotstyle_add_point, zorder=zorder+1)
        else: #interpolation plot
            ax.errorbar(columns[0], columns[1], fmt = '-', lw = 0.5, mew=0, color=color, zorder=zorder+1)
            #ax.fill_between(columns[0], columns[1]-columns[2], columns[1]+columns[2], lw=0, alpha=0.5, facecolor=color, zorder=zorder)        
            plots.append(ax.errorbar(tauTs, [val * lpd.improve_corr_factor(v, nt, i, use_imp) for v,val in enumerate(EE[i,:])], [val * lpd.improve_corr_factor(v, nt, i, use_imp) for v,val in enumerate(EE_err[i,:])], **lp.plotstyle_add_point, color=color, label=label, zorder=zorder+2))
if len(plots) < 4:
    exit("Error: need at least 4 lattices")
#try:
    #EE_cont = numpy.loadtxt(lp.inputfolder+"/continuum_limit/EE_"+'{0:.4f}'.format(lp.flow_radius[i])+"_cont.txt", skiprows=lp.n_skiprows)
    #for j in range(0,len(EE_cont)):
        #if EE_cont[j,0] < lp.flow_radius[i]/numpy.sqrt(8*0.014)+lp.offset:
            #EE_cont[j,:] = None
    #plots.append(ax.fill_between(EE_cont[:,0], EE_cont[:,1]-EE_cont[:,2], EE_cont[:,1]+EE_cont[:,2], facecolor='grey', zorder=-6))
    #ax.legend(handles=plots, labels=[r'$N_\tau = 16$', r'$N_\tau = 20$', r'$N_\tau = 24$', r'$N_\tau = 30$', r'cont.'], title="", **lp.legendstyle)
    ##plot a minimum linewidth of cont extr for visibility
    #ax.errorbar(EE_cont[:,0], EE_cont[:,1], color='grey', lw=1.5, fmt='-', zorder=-5, solid_capstyle='butt')
#except:
ax.legend(handles=plots)
handles, labels = ax.get_legend_handles_labels()
leg = ax.legend(handles[::-1], labels[::-1], title="", **lp.legendstyle)
leg.set_zorder(-1)
for line in leg.get_lines():
    line.set_linewidth(4.0)
    #pass
if i != 0:
    ax.axvline(x=lpd.lower_tauT_limit(flowradius), **lp.verticallinestyle)
matplotlib.pyplot.tight_layout(0.2)
fig.savefig(lp.outputfolder+"/lattice_effects/EE_flow_"+'{0:.3f}'.format(lp.flow_radius[i])+".pdf")#, **lp.figurestyle) 
print("saved lattice effects plot", lp.outputfolder+"/lattice_effects/EE_flow_"+'{0:.3f}'.format(lp.flow_radius[i])+".pdf") 
#ax.lines.clear() ; ax.collections.clear() ; plots.clear();
#ax.set_title(r'')

print("done with lattice effects plot") 
