#!/usr/local/bin/python
import numpy
import matplotlib
import lib_process_data as lpd
import lib_plot as lp

lattices_orig = ["192_48" , "144_36", "96_24", "80_20", "64_16", "cont"] 
lattices_orig_label = ["48", "36", "24", "20", "16", "cont"]#"48", 
lattices_new = ["48", "36", "24", "20", "16", "cont"]#"48", 

data_orig = []
data_new = []

for lat in lattices_orig[:-1]:
    data_orig.append(numpy.loadtxt("diffusion_"+lat+".dat"))
data_orig.append(numpy.loadtxt("cont_thomas.dat"))

for lat in lattices_new[:-1]:
    data_new.append(numpy.loadtxt("EE_2015_Nt"+lat+".dat"))
data_new.append(numpy.loadtxt("EE_2015_new.txt"))

ylabel=r'$ \displaystyle \left.\frac{G^\mathrm{latt}}{G^\mathrm{norm}} \right|_{\mathrm{imp. dist.} } \kern-2.5em (\tau T) $'

lp.plotstyle_add_point.update(dict(markersize=4))
fig, ax, plots = lp.create_figure(xlims=[0.0,0.51], ylims=[1,4], xlabelpos=(0.95,0.05), ylabelpos=(0.05,0.95), xlabel=r'$\tau T$', ylabel=ylabel,  figsize=None, UseTex=True)
for data,label,marker in zip(data_orig,lattices_orig_label,lp.markers):
    lp.plotstyle_add_point.update(dict(fmt=marker+'-'))
    ax.errorbar(data[:,0], data[:,1], data[:,2], label=label, **lp.plotstyle_add_point)
ax.legend(**lp.legendstyle)
matplotlib.pyplot.tight_layout(0.2)
fig.savefig("EE_2015_original.pdf")
ax.lines.clear() ; ax.collections.clear() ; plots.clear();


ylabel=r'$ \displaystyle \left.\frac{G^\mathrm{latt}}{G^\mathrm{norm}} \right|_{\mathrm{imp. corr. + updated Z} } \kern-7.5em (\tau T) $'
fig, ax, plots = lp.create_figure(xlims=[0.0,0.51], ylims=[1,4],  xlabelpos=(0.95,0.05), ylabelpos=(0.05,0.95),  xlabel=r'$\tau T$', ylabel=ylabel,  figsize=None, UseTex=True)
for data,label,marker in zip(data_new,lattices_new,lp.markers):
    #if label == "36":
        #ax.errorbar(data[:11,0], data[:11,1], data[:11,2], label=label, **lp.plotstyle_add_point)
    #else:
    lp.plotstyle_add_point.update(dict(fmt=marker+'-'))
    ax.errorbar(data[:,0], data[:,1], data[:,2], label=label, **lp.plotstyle_add_point)
ax.legend(**lp.legendstyle)
matplotlib.pyplot.tight_layout(0.2)
fig.savefig("EE_2015_new.pdf")

