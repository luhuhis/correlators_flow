#!/usr/local/bin/python

import lib_process_data as lpd

#input

#   inputfolder+"flowradii_"+conftype+".dat"
#   inputfolder+"XX_"+conftype+".dat"
#   inputfolder+"XX_err_"+conftype+".dat"

#output

#   outputfolder+conftype+"_XX_fill.pdf"
#   outputfolder+"/"+conftype+"_XX_flowlim_fill.pdf"


import numpy
import re
import sys
import matplotlib

def skip_large_errors(XX, XX_err, boundary):
    for i in range(0,XX.shape[0]):
        for j in range(0,XX.shape[1]):
            if XX_err[i,j] > boundary :
                XX[i,j] = None


def main():
    
    conftype, beta, ns, nt, nt_half, qcdtype, fermions, temp, flowtype, corr, add_params = lpd.read_args()
    inputfolder=lpd.get_merged_data_path(qcdtype,corr,conftype)
    outputfolder=lpd.get_plot_path(qcdtype,corr,conftype)
           
    """load data"""
    lpd.create_folder(outputfolder) 
    flow_radius = numpy.loadtxt(inputfolder+"flowradii_"+conftype+".dat")
    XX        = numpy.loadtxt(inputfolder+corr+"_"+conftype+".dat")
    XX_err = numpy.loadtxt(inputfolder+corr+"_err_"+conftype+".dat")
    
    normalization_factor = 1
    if corr == "EE_clover":
        normalization_factor = 2 #FIXME FIXME FIXME remove this for new measurements (there its already included in parallelgpucode!)
    if corr == "BB_clover":
        normalization_factor = -1.5
    if corr == "BB":
        normalization_factor = 1


    for i in range(len(flow_radius)):
    #for i in range(150):
        for j in range(nt_half):
            XX[i,j] = XX[i,j] * lpd.improve_corr_factor(j, nt, i) *normalization_factor
            XX_err[i,j] = XX_err[i,j] * lpd.improve_corr_factor(j, nt, i) *normalization_factor #/ 10 #FIXME FIXME FIXME FIXME FIXME FIXME REMOVE ERROR REDUCTION FOR TEMPORARY BETTER VISIBILITY
    tauT = lpd.get_tauTs(nt)
    
    if fermions == "hisq":
        ylims = ([-0.5,9.5])
        if corr == "BB_clover" or corr == "BB":
            ylims = ([-0.5,6])
    elif fermions == "quenched":
        ylims=([0,4]) 
        
    """PLOT: x-axis tauT, flowtimes fixed"""
    fig, ax, plots = lpd.create_figure(xlims=[0,0.505], 
                                    ylims=ylims,
                                    xlabel=r'$\tau T$', 
                                    xlabelpos=(0.95,0.06),
                                    ylabel=r'$\displaystyle\frac{G^\mathrm{latt }_{\tau_\mathrm{F}}(\tau)}{G_{\tau_\mathrm{F}=0}^{\begin{subarray}{l}\mathrlap{\textrm{\tiny  norm}}\\[-0.3ex] \textrm{\tiny latt}\end{subarray}}  (\tau)    }$',
                                    #ylabel=r'$G$',
                                    UseTex=True)
    #lpd.legendstyle.update(dict(loc="upper right"))
    #ax.yaxis.set_label_coords(0.01,0.99)

    #skip_large_errors(XX, XX_err, XX_backup, XX_err_backup, 10) 

        
    flow_selection = range(0,142,10)

    for i in flow_selection:
        mycolor = lpd.get_color(flow_radius, i, flow_selection[0], flow_selection[-1]+1)
        plots.append(ax.fill_between(list(tauT), XX[i,:]-XX_err[i,:], XX[i,:]+XX_err[i,:], facecolor=mycolor, linewidth=lpd.mylinewidth, zorder=-flow_selection[-1]+i))
        ax.errorbar(list(tauT), XX[i,:], color=mycolor, **lpd.plotstyle_lines, zorder=-flow_selection[-1]+i+1)
        ax.errorbar(list(tauT), XX[i,:], color=mycolor, fmt='D', linewidth=0.5, mew=0.25, mec='grey', markersize=1.5, capsize=3, zorder=-flow_selection[-1]+i+2)
        

    leg = ax.legend(handles=plots, labels=['{0:.3f}'.format(j) for i,j in enumerate(flow_radius) if i in flow_selection], title=r"$ \sqrt{8\tau_\mathrm{F}}T$", **lpd.legendstyle)
    matplotlib.pyplot.tight_layout(0)
    ### change first xtick label
    fig.canvas.draw()
    ticks=ax.get_xticks().tolist()
    ticks = ['{0:.1f}'.format(x) for x in ticks]
    ticks[0]='0'
    ax.set_xticklabels(ticks)
    fig.savefig(outputfolder+conftype+"_"+corr+".pdf") 
    print("saved correlator plot", outputfolder+conftype+"_"+corr+".pdf")
    ax.lines.clear() ; ax.collections.clear() ; plots.clear()


    """PLOT: x-axis flowtimes, tauT fixed"""
    xlims=[-0.0001,0.0028125+0.00015]
    if fermions == "hisq":
        xlims = [0,0.009] #FIXME reset left limit to -0.0001
    
    fig, ax, plots = lpd.create_figure(xlims=xlims, ylims=ylims, 
                                    xlabelpos=[0.94,0.05],
                                    ylabelpos=[0.4,0.97], #FIXME reset to [0.01,0.97]
                                    xlabel=r'$ \tau_\mathrm{F} T^2 $', 
                                    ylabel=r'$\displaystyle\frac{G^\mathrm{latt }_{\tau}(\tau_\mathrm{F})}{G_{\tau,{\tau_\mathrm{F}}=0}^{\substack{ \text{\tiny  norm} \\[-0.4ex] \text{\tiny latt } } } }$',
                                    #ylabel=r'$G$',
                                    UseTex=True)
    lpd.plotstyle_add_point.update(dict(fmt=',', capsize=0.6))


    """calculate flow limits"""
    flow_limit=numpy.zeros(len(tauT))
    interpolation_value=numpy.zeros(len(tauT))
    for i in range(len(tauT)):
        flow_limit[i] = lpd.upper_flow_limit(tauT[i])
        index = (numpy.abs(flow_radius[:] - flow_limit[i])).argmin()
        offset = 1 if flow_limit[i]-flow_radius[index] > 0 else -1
        offset2 = -1 if offset == -1 else 0
        interpolation_value[i] = XX[index+offset2,i]-abs(XX[index,i]-XX[index+offset,i])/abs(flow_radius[index]-flow_radius[index+offset])*abs(flow_limit[i]-flow_radius[index+offset2])



    #decide what flow times to plot how
    
    error_factor = 0.025 #FIXME change this back to 0.0075 
    minimum_trusted_flowtime_index = 0 #FIXME change this back to 50. Set this to the index where the flow radius is larger than one lattice spacing of the coarsest lattice used in cont extr

    max_flow_index = 142 #this is the max because of the coarsest lattice (too few data points to put a spline through at this flow time)
    for i in range(nt_half):
        if True: #i%2 == 0: #FIXME only use half because of large lattice
            flow_extr_filter_high = 0 #determines up to which index the discrete error bars should be plotted for this tauT
            flow_extr_filter_low = -1 # determines from which index on the discrete error bars should be plotted
            flow_extr_filter_low_grey = -1 # determines from which index on the grey band should be plotted
            flow_extr_filter_low_greyline = -1 # determines from which index on the grey line should be plotted
            
            rel_err = numpy.fabs(XX_err[:,i]/XX[:,i])
            for j,val in enumerate(rel_err):
                if flow_extr_filter_high == 0:
                    if flow_limit[i] < flow_radius[j]:
                        flow_extr_filter_high = j
                if flow_extr_filter_low == -1:
                    if val <= error_factor: #FIXME change this back to error_factor *tauT[i]:
                        flow_extr_filter_low = max(j,minimum_trusted_flowtime_index)
                if flow_extr_filter_low_greyline == -1:
                    if val <= 0.25: #FIXME change this back to error_factor *tauT[i]:
                        flow_extr_filter_low_greyline = max(j,minimum_trusted_flowtime_index)
                if flow_extr_filter_low_grey == -1:
                    if val <= 0.05:
                        flow_extr_filter_low_grey = j
            flow_extr_filter_high = min(flow_extr_filter_high,max_flow_index)
            
            #FIXME temporarily overwrite settings for hisq for testing purposes
            if fermions == "hisq":
                #flow_extr_filter_low = 40
                flow_extr_filter_high = -1
                flow_extr_filter_low_grey = 1
                
            xdata = flow_radius[flow_extr_filter_low:flow_extr_filter_high]**2/8
            ydata = XX[flow_extr_filter_low:flow_extr_filter_high,i]
            edata = XX_err[flow_extr_filter_low:flow_extr_filter_high,i]
            
            # explicit colorful datapoints
            ax.errorbar(xdata, ydata, edata, color=lpd.get_color(tauT, i, 0, nt_half), zorder=(-nt_half+i), **lpd.plotstyle_add_point)
            
            #FIXME replace flow_extr_filter_low_greyline with flow_extr_filter_low_grey here!
            xdata = flow_radius[flow_extr_filter_low_greyline:]**2/8
            ydata = XX[flow_extr_filter_low_greyline:,i]
            edata = XX_err[flow_extr_filter_low_greyline:,i]
            
            zorder = 100*(-nt_half+i)
            
            # FIXME uncomment this for light grey lines
            #ax.errorbar(xdata, ydata, color='lightgrey', zorder=zorder, fmt='-', lw=0.5, markersize=0.5)  
            
            #FIXME uncomment this for light grey backgrounds
            ax.fill_between(xdata, ydata-edata, ydata+edata, linewidth=0.5, color=(0.9,0.9,0.9,1), zorder=zorder)  #
            
            
            #FIXME remove these three lines to not cut off the dark grey lines at low end and if you uncomment the perturbative flow limit markes FIXME
            xdata = flow_radius[flow_extr_filter_low_greyline:]**2/8
            ydata = XX[flow_extr_filter_low_greyline:,i]
            edata = XX_err[flow_extr_filter_low_greyline:,i]
            
            #dark grey lines
            #FIXME reset color to color=(0.7,0.7,0.7,1)
            #FIXME remove plots.append(...) and label=... because we want perturbative markes for tauT identification
            plots.append(ax.errorbar(xdata, ydata, linewidth=0.5, color=lpd.get_color(tauT, i, 0, nt_half), zorder=zorder/100+1, label='{0:.3f}'.format(tauT[i])))
            #ax.axvline(x=flow_radius[minimum_trusted_flowtime_index]**2/8, **lpd.verticallinestyle)
            
            #FIXME uncomment this for perturbative flow limit markers
            #plots.append(ax.errorbar(flow_limit[i]**2/8, interpolation_value[i], marker=lpd.markers[i%len(lpd.markers)], fillstyle='none', markersize=4, mew=0.25, color=lpd.get_color(tauT, i), label='{0:.3f}'.format(tauT[i]), capsize=0, lw=0 ))
            #ax.errorbar(flow_limit[i]**2/8, interpolation_value[i], marker='|', fillstyle='none', markersize=4, mew=0.25, color=lpd.get_color(tauT, i, 0, nt_half), capsize=0 )
            
    ax2=ax.twiny()
    #new_tick_locations = numpy.array([0,0.05, 0.08, 0.1, 0.12, 0.13, 0.14, 0.15])**2/8
    new_tick_locations = numpy.array([0, 0.1, 0.15, 0.2, 0.25])**2/8  #FIXME
    #new_tick_locations = numpy.array([0,0.05**2/8, 0.08**2/8, 0.1**2/8, 0.12**2/8, 0.13**2/8, 0.14**2/8, 0.15**2/8])
    def tick_function(X):
        V = numpy.sqrt(X*8)
        return ["%.1f" % z if z == 0 else "%.2f" % z for z in V ]
    ax2.tick_params(direction='in', pad=0, width=0.5)
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(new_tick_locations)
    ax2.set_xticklabels(tick_function(new_tick_locations))
    ax2.set_xlabel(r'$\sqrt{8\tau_\mathrm{F}}T$')#, horizontalalignment='right', verticalalignment='top', bbox=lpd.labelboxstyle, zorder=99)
    ax2.xaxis.set_label_coords(0.92,0.92)


    #FIXME temporarily set legend outside of axes, reset to loc="center right"
    lpd.legendstyle.update(dict(loc="center right", bbox_to_anchor=(1,0.4), columnspacing=0.3, handletextpad=0.2, handlelength=0.5, handler_map={matplotlib.container.ErrorbarContainer: matplotlib.legend_handler.HandlerErrorbar(xerr_size=0)}))
    #FIXME enable legend again (uncomment three lines below)
    ax.legend(handles=plots)    
    handles, labels = ax.get_legend_handles_labels()
    leg = ax.legend(handles[::-1], labels[::-1], ncol=2, title=r"$\tau T=$", **lpd.legendstyle)
    #FIXME remove the two lines below
    for line in leg.get_lines():
        line.set_linewidth(4.0)
    
    matplotlib.pyplot.tight_layout(0)
    ### change first xtick label
    fig.canvas.draw()
    ticks=ax.get_xticks().tolist()
    ticks[0]='0.0'
    ax.set_xticklabels(ticks)
    fig.savefig(outputfolder+"/"+conftype+"_"+corr+"_flow.pdf") 
    print("saved flow effect plot", outputfolder+"/"+conftype+"_"+corr+"_flow.pdf")

if __name__ == '__main__':
    main()
