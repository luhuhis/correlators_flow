#!/usr/local/bin/python
import numpy
import matplotlib
import sys
import itertools
import lib_plot as lp
import lib_process_data as lpd

try:
    flowindex = sys.argv[1]
    i = int(flowindex) #create alias
except IndexError:
    exit("Invalid Arguments. Usage: script.py <flowindex>, e.g. script.py 100")

lpd.create_folder(lp.outputfolder+"/lattice_effects/")
fig, ax, lpots = lp.create_figure(xlims=[0.1,0.505], ylims=[2.2,3.9], xlabel=r'$\tau T$', ylabel=r'$ \frac{G_{\tau_F}(\tau)}{G_{\tau_F = 0}^{\mathrm{ norm } }(\tau)}$', UseTex = False)

lp.titlestyle.update(dict(y=0.96))
lp.legendstyle.update(dict(loc="lower right", bbox_to_anchor=(0.99,0.01), framealpha=0.5, handlelength=0.5, handler_map={matplotlib.container.ErrorbarContainer: matplotlib.legend_handler.HandlerErrorbar(xerr_size=0.4)}))
lp.plotstyle_points.update(dict(fmt="D", markersize = 1.5, mew=0.25, lw=0, elinewidth=0.25, capsize=1))
#ax.xaxis.set_label_coords(0.99,0.01)
#ax.yaxis.set_label_coords(0.01,0.97)

#load interpolations
EE_ints = []
for conftype in ("s064t16_b0687361", "s080t20_b0703500", "s096t24_b0719200", "s120t30_b0739400", "s144t36_b0754400"):
    try:
        EE_ints.append(numpy.loadtxt(lp.inputfolder+"/"+conftype+"/interpolations/EE_"+'{0:.3f}'.format(lp.flow_radius[i])+"_interpolation.txt", unpack=True))
    except:
        EE_ints.append(None)
try:
    EE_ints.append(numpy.loadtxt(lp.inputfolder+"/cont_extr/EE_"+'{0:.3f}'.format(lp.flow_radius[i])+"_cont.txt"))
except:
    EE_ints.append(None)
    
colors = matplotlib.cm.gnuplot([0.05, 0.2, 0.4, 0.75, 0.9, 0])
labels = [r'$N_\tau = 16$', r'$N_\tau = 20$', r'$N_\tau = 24$', r'$N_\tau = 30$', r'$N_\tau = 36$', r'cont.']
plots = []
ax.set_title(r'$\sqrt{8\tau_F} T=$ '+'{0:.3f}'.format(lp.flow_radius[i]), **lp.titlestyle)
for k, EE_int, tauT, EE, EE_err, label, color in itertools.zip_longest(range(-30,0,4), EE_ints, lp.tauT, lp.EE, lp.EE_err, labels, colors):
    if tauT is not None:
        nt = int(len(tauT)*2)
    if EE_int is not None:
            ax.fill_between(EE_int[0], EE_int[1]-EE_int[2], EE_int[1]+EE_int[2], lw=0, alpha=0.5, facecolor=color, zorder=k)
            plots.append(ax.errorbar(EE_int[0], EE_int[1], label=label, fmt = '-', lw = 0.5, mew=0, color=color, zorder=k+1))
            if tauT is not None:
                ax.errorbar(tauT, [val / lpd.norm_corr(lpd.tauT_imp[nt][v]) * nt**4 for v,val in enumerate(EE[i,:])], [val / lpd.norm_corr(lpd.tauT_imp[nt][v]) * nt**4 for v,val in enumerate(EE_err[i,:])], **lp.plotstyle_points, color=color, zorder=k+2) 
if len(plots) <= 1:
    exit("Error: need at least 2 lattices")
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
leg = ax.legend(handles=plots, title="", **lp.legendstyle)
for line in leg.get_lines():
    line.set_linewidth(4.0)
    #pass
if i != 0:
    ax.axvline(x=(lp.flow_radius[i]/numpy.sqrt(8*0.014)), **lp.verticallinestyle)
fig.savefig(lp.outputfolder+"/lattice_effects/EE_flow_"+'{0:.3f}'.format(lp.flow_radius[i])+".pdf", **lp.figurestyle) 
print("saved lattice effects plot", lp.outputfolder+"/lattice_effects/EE_flow_"+'{0:.3f}'.format(lp.flow_radius[i])+".pdf") 
#ax.lines.clear() ; ax.collections.clear() ; plots.clear();
#ax.set_title(r'')

print("done with lattice effects plot") 
