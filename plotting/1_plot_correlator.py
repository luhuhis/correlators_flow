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

def plot1(XX, XX_err, args,prefix, flow_selection, flow_var, xdata_plot, nt_norm_half, flow_str, outputfolder):
    global xlims
    global xlabel
    global ylims
    global ylabelpos
    global xlabelpos 
    global figsize 
    xlims = [0,0.505]
    xlabel=r'$\tau T$'
    ylabelpos = (0.01,0.97)
    xlabelpos = (0.95,0.06)
    figsize = (3+3/8,3+3/8-1/2.54)
    ylims=([1,100000])
    ylabel=r'$\displaystyle\frac{G^\mathrm{latt }_{\tau_\mathrm{F}}(\tau)}{T_\mathrm{norm}^4 }$'
    global corr_label
    corr_label = r'X'
    if args.corr == "EE":
        corr_label = 'E'
    if args.corr == "BB":
        corr_label = 'B'
    if args.show_TauByA:
        xlims = [0,nt_norm_half*1.01]
        xlabel = r'$\tau/a$'
    if prefix == "_numerator":
        ylims = [-50,17]
        ylabel = r'$\displaystyle\frac{ \langle '+corr_label+corr_label+r' \rangle_{\tau_\mathrm{F}}(\tau)}{T_\mathrm{norm}^4 }$'
    
    if prefix == "_polyakovloop":
        ylims = [1,1000000]
        ylabel = r'$\displaystyle\frac{ \langle - \rangle_{\tau_\mathrm{F}}(\tau)}{T_\mathrm{norm}^4 }$'
        
    if args.testmode:
        ylabelpos = (-0.2,0.95)
        xlabelpos = (0.5,-0.06) 
        figsize = (1.5*(3+3/8),1.5*(3+3/8)/16*9)
    
    if args.reconstruct:
        ylabel = r'$\displaystyle\frac{G^\mathrm{rec }_{\tau_\mathrm{F}}(\tau)}{T_\mathrm{norm}^4 }$'
    
    
    
    fig, ax, plots = lpd.create_figure(xlims=xlims, ylims=ylims, xlabel=xlabel, xlabelpos=xlabelpos, ylabel=ylabel, ylabelpos=ylabelpos, UseTex=True, figsize=figsize)
    
    #lpd.legendstyle.update(dict(loc="upper right"))
    
    ax.set_yscale('log', nonposy='mask')
    if prefix == "_numerator":
        ax.set_yscale('linear')


    #skip_large_errors(XX, XX_err, XX_backup, XX_err_backup, 10) 

    
        
    for i in flow_selection:
        mycolor = lpd.get_color(flow_var, i, flow_selection[0], flow_selection[-1]+1)
        ax.errorbar(list(xdata_plot[:nt_norm_half]), XX[i,:nt_norm_half], color=mycolor, **lpd.plotstyle_lines, zorder=-flow_selection[-1]+i+1)
        plots.append(ax.errorbar(list(xdata_plot[:nt_norm_half]), XX[i,:nt_norm_half], XX_err[i,:nt_norm_half], color=mycolor, fmt='D', linewidth=0.5, mew=0.25, mec='grey', markersize=1.5, capsize=3, zorder=-flow_selection[-1]+i+2))
        

    leg = ax.legend(handles=plots, labels=['{0:.3f}'.format(j) for i,j in enumerate(flow_var) if i in flow_selection], title=flow_str, **lpd.legendstyle)
    matplotlib.pyplot.tight_layout(0)
    ### change first xtick label
    fig.canvas.draw()
    ticks=ax.get_xticks().tolist()
    ticks = ['{0:.1f}'.format(x) for x in ticks]
    ticks[0]='0'
    ax.set_xticklabels(ticks)
    filename = outputfolder+args.conftype+"_"+args.corr+prefix+"_T.pdf"
    fig.savefig(filename) 
    print("saved correlator plot", filename)
    ax.lines.clear() ; ax.collections.clear() ; plots.clear()
    
def plot2(XX, XX_err, args,prefix, flow_selection,flow_var,xdata_plot, nt_norm_half, flow_str, outputfolder, fermions):
    global xlims
    global xlabel
    global ylims
    global ylabelpos
    global xlabelpos 
    global figsize 
    global corr_label
    
    
    # ======= CUSTOM OPTIONS ======= #
    if fermions == "hisq":
        ylims = ([-0.5,10])
        #if corr == "BB_clover" or corr == "BB":
            #ylims = ([-0.5,6])
    elif fermions == "quenched":
        ylims=([0,4]) 
    
    ylabel = r'$\displaystyle\frac{G^\mathrm{latt }_{\tau_\mathrm{F}}(\tau)}{G_{\tau_\mathrm{F}=0}^{\begin{subarray}{l}\mathrlap{\textrm{\tiny  norm}}\\[-0.3ex] \textrm{\tiny latt}\end{subarray}}  (\tau)    }$'
    if prefix == "_numerator":
        ylims = [-0.2,1.3]
        ylabel = r'$\displaystyle\frac{ \langle '+corr_label+corr_label+r' \rangle_{\tau_\mathrm{F}}(\tau)}{G_{\tau_\mathrm{F}=0}^{\begin{subarray}{l}\mathrlap{\textrm{\tiny  norm}}\\[-0.3ex] \textrm{\tiny latt}\end{subarray}}  (\tau)    }$'
    
    if prefix == "_polyakovloop":
        ylims = [1,1000000]
        ylabel = r'$\displaystyle\frac{ \langle - \rangle_{\tau_\mathrm{F}}(\tau)}{G_{\tau_\mathrm{F}=0}^{\begin{subarray}{l}\mathrlap{\textrm{\tiny  norm}}\\[-0.3ex] \textrm{\tiny latt}\end{subarray}}  (\tau)    }$'
    
    if args.reconstruct:
        ylabel = r'$\displaystyle\frac{G^\mathrm{rec }_{\tau_\mathrm{F}}(\tau)}{G_{\tau_\mathrm{F}=0}^{\begin{subarray}{l}\mathrlap{\textrm{\tiny  norm}}\\[-0.3ex] \textrm{\tiny latt}\end{subarray}}  (\tau)    }$'
    
    if args.plot_for_peter:
        ylabelpos = [0.025,0.95]
    
    
    fig, ax, plots = lpd.create_figure(xlims=xlims, ylims=ylims, xlabel=xlabel, xlabelpos=xlabelpos, ylabelpos=ylabelpos,
                                    ylabel=ylabel,
                                    UseTex=True,figsize=figsize)
    
    
    #lpd.legendstyle.update(dict(loc="upper right"))
    #skip_large_errors(XX, XX_err, XX_backup, XX_err_backup, 10) 

    if args.plot_for_peter:
        global BB
        global BB_err
        fmts=['s','o','D','v']
        for i in flow_selection:
            mycolor = lpd.get_color(flow_var, i, flow_selection[0], flow_selection[-1]+1)
            ax.errorbar(list(xdata_plot[:nt_norm_half]), XX[i,:nt_norm_half], XX_err[i,:nt_norm_half],  color=mycolor, fmt=fmts[i], linewidth=0.5, mew=0.4, mec=mycolor, markersize=4, capsize=3, zorder=-flow_selection[-1]+i+2, label='{0:.2f}'.format(flow_var[i])+r', \text{EE}', fillstyle='full')
        for i in flow_selection:
            mycolor = lpd.get_color(flow_var, i, flow_selection[0], flow_selection[-1]+1)
            ax.errorbar(list(xdata_plot[:nt_norm_half]), BB[i,:nt_norm_half], BB_err[i,:nt_norm_half], color=mycolor, fmt=fmts[i], linewidth=0.5, mew=0.4, mec=mycolor, markersize=4, capsize=3, zorder=-flow_selection[-1]+i+2, label='{0:.2f}'.format(flow_var[i])+r', \text{BB}', fillstyle='none')
        
        leg = ax.legend(title=flow_str, **lpd.legendstyle)
    
    else:       
        for i in flow_selection:
            mycolor = lpd.get_color(flow_var, i, flow_selection[0], flow_selection[-1]+1)
            plots.append(ax.fill_between(list(xdata_plot[:nt_norm_half]), XX[i,:nt_norm_half]-XX_err[i,:nt_norm_half], XX[i,:nt_norm_half]+XX_err[i,:nt_norm_half], facecolor=mycolor, linewidth=lpd.mylinewidth, zorder=-flow_selection[-1]+i))
            ax.errorbar(list(xdata_plot[:nt_norm_half]), XX[i,:nt_norm_half], color=mycolor, **lpd.plotstyle_lines, zorder=-flow_selection[-1]+i+1)
            ax.errorbar(list(xdata_plot[:nt_norm_half]), XX[i,:nt_norm_half], color=mycolor, fmt='D', linewidth=0.5, mew=0.25, mec='grey', markersize=1.5, capsize=3, zorder=-flow_selection[-1]+i+2)
        leg = ax.legend(handles=plots, labels=['{0:.3f}'.format(j) for i,j in enumerate(flow_var) if i in flow_selection], title=flow_str, **lpd.legendstyle)

    matplotlib.pyplot.tight_layout(0)
    ### change first xtick label
    fig.canvas.draw()
    ticks=ax.get_xticks().tolist()
    ticks = ['{0:.1f}'.format(x) for x in ticks]
    ticks[0]='0'
    ax.set_xticklabels(ticks)
    filename = outputfolder+args.conftype+"_"+args.corr+prefix+".pdf"
    fig.savefig(filename) 
    print("saved correlator plot", filename)
    ax.lines.clear() ; ax.collections.clear() ; plots.clear()

def plot3(XX, XX_err, args,prefix, flow_var,xdata_plot, nt_norm_half, flow_str, outputfolder, fermions, nt_half, valid_flowtimes,tauT, flow_radius, flowend, flow_selection):
    global xlims
    global xlabel
    global ylims
    global ylabelpos
    global xlabelpos 
    global figsize 
    global corr_label
    
    ylabelpos = (0.01,0.97) 
    xlabelpos = (0.94,0.05)
    
    if args.testmode:
        ylabelpos = (-0.2,0.95)
        xlabelpos = (0.5,-0.06)
    
    max_tau = nt_half
    number_format = '{0:.3f}'
    xlabel = r'$ \tau_\mathrm{F} T^2 $'
    xlabel2 = r'$\sqrt{8\tau_\mathrm{F}}T$'
    legend_title = r"$\tau T=$"
    #xlims=[-0.0001,(flow_var[flow_selection[-1]])*2/8 *1.01] #0.0028125+0.00015
    #xlims=[-0.0001,0.0028125+0.00015] 
    xlims=[-0.0001,0.02+0.00015] 
    if fermions == "quenched":
        new_tick_locations = numpy.array([0, 0.05, 0.1, 0.13])**2/8  
    if fermions == "hisq":
        new_tick_locations = numpy.array([0, 0.1, 0.15])**2/8  
        ylims=(-0.5, 12)
        ylabelpos = (0.05,0.97)
        xlabelpos = (0.94,0.1)
    
    xlabel2_format = "%.2f"
    if args.show_flowtime_instead_of_flowradii:
        xlabel = r'$ \tau_\mathrm{F} /a^2 $'
        xlabel2 = r'$\sqrt{8\tau_\mathrm{F}}/a$'
        legend_title = r"$\tau/a =$"
        if valid_flowtimes is not None:
            xlims = [-0.05,valid_flowtimes[-1]*1.01]
        number_format = '{0:.0f}'
        new_tick_locations = numpy.array([0, 3, 4, 5, 6, 7, 8, 9])**2/8  
        xlabel2_format = "%.1f"
        max_tau = 16 #FIXME hardcoded...
    
    #if fermions == "hisq":
        #xlims = [0,0.009] #FIXME reset left limit to -0.0001
    
    
    ylabel = r'$\displaystyle\frac{G^\mathrm{latt }_{\tau}(\tau_\mathrm{F})}{G_{\tau,{\tau_\mathrm{F}}=0}^{\substack{ \text{\tiny  norm} \\[-0.4ex] \text{\tiny latt } } } }$'
    if prefix == "_numerator":
        ylims = [-0.2,2]
        ylabel = r'$\displaystyle\frac{ \langle '+corr_label+corr_label+r' \rangle_{\tau_\mathrm{F}}(\tau)}{G_{\tau,{\tau_\mathrm{F}}=0}^{\substack{ \text{\tiny  norm} \\[-0.4ex] \text{\tiny latt } } } }$'
    
    if prefix == "_polyakovloop":
        ylims = [1,1000000]
        ylabel = r'$\displaystyle\frac{ \langle - \rangle_{\tau_\mathrm{F}}(\tau)}{T_\mathrm{norm}^4 }$'
    
    if args.reconstruct:
        ylabel = r'$\displaystyle\frac{G^\mathrm{rec }_{\tau_\mathrm{F}}(\tau)}{G_{\tau_\mathrm{F}=0}^{\begin{subarray}{l}\mathrlap{\textrm{\tiny  norm}}\\[-0.3ex] \textrm{\tiny latt}\end{subarray}}  (\tau)    }$'
    
    figsize = ((3+3/8),(3+3/8)-1/2.54)
    if args.testmode:
        figsize=(1.5*(3+3/8),1.5*(3+3/8)/16*9)
    fig, ax, plots = lpd.create_figure(xlims=xlims, ylims=ylims, xlabelpos=xlabelpos, ylabelpos=ylabelpos, xlabel=xlabel, ylabel=ylabel, UseTex=True, 
    figsize=figsize)

    lpd.plotstyle_add_point.update(dict(fmt=',', capsize=0.6))

    if prefix == "_polyakovloop":
        ax.set_yscale('log', nonposy='mask')



    """calculate flow limits"""
    flow_limit=numpy.zeros(len(tauT))
    interpolation_value=numpy.zeros(len(tauT))
    length=len(tauT)
    if args.conftype_2:
        length=nt_norm_half
    for i in range(length): 
        flow_limit[i] = lpd.upper_flow_limit(tauT[i])
        index = (numpy.abs(flow_radius[:] - flow_limit[i])).argmin()
        offset = 1 if flow_limit[i]-flow_radius[index] > 0 else -1
        offset2 = -1 if offset == -1 else 0
        interpolation_value[i] = XX[index+offset2,i]-abs(XX[index,i]-XX[index+offset,i])/abs(flow_radius[index]-flow_radius[index+offset])*abs(flow_limit[i]-flow_radius[index+offset2])



    #decide what flow times to plot how
    if fermions == "hisq":
        error_factor = 0.5
        minimum_trusted_flowtime_index = 0 #Set this to the index where the flow radius is larger than one lattice spacing of the coarsest lattice used in cont extr
    if fermions == "quenched":
        error_factor = 0.05 #0.0075
        minimum_trusted_flowtime_index = 20 

    max_flow_index = flowend #this is the max because of the coarsest lattice (too few data points to put a spline through at this flow time)
    
    if args.tauT_selection:
        tauT_selection = args.tauT_selection
    else:
        tauT_selection = range(max_tau)
    for i in tauT_selection:
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
                    if fermions == "hisq" and val <= error_factor or fermions == "quenched" and val <= error_factor *tauT[i]: 
                        flow_extr_filter_low = max(j,minimum_trusted_flowtime_index)
                if flow_extr_filter_low_greyline == -1:
                    if fermions == "hisq" and val <= 1 or fermions == "quenched" and val <= error_factor *tauT[i]: 
                        flow_extr_filter_low_greyline = max(j,minimum_trusted_flowtime_index)
                if flow_extr_filter_low_grey == -1:
                    if val <= 0.10:
                        flow_extr_filter_low_grey = j
            flow_extr_filter_high = min(flow_extr_filter_high,max_flow_index)
            
            #FIXME temporarily overwrite settings for hisq for testing purposes
            if fermions == "hisq":
                #flow_extr_filter_low = 0
                flow_extr_filter_high = flowend
                flow_extr_filter_low_grey = 1
                
            xdata = flow_var[flow_extr_filter_low:flow_extr_filter_high]**2/8
            ydata = XX[flow_extr_filter_low:flow_extr_filter_high,i]
            edata = XX_err[flow_extr_filter_low:flow_extr_filter_high,i]
            
            # --- explicit colorful datapoints for low error points ---
            lpd.plotstyle_add_point.update(fmt='x', markersize=1.25)
            #plots.append(
            ax.errorbar(xdata, ydata, edata, color=lpd.get_color(xdata_plot, i, 0, max_tau), zorder=(-max_tau+i), **lpd.plotstyle_add_point)
            #, label=number_format.format(xdata_plot[i])))
            
            zorder = 100*(-max_tau+i)
            
            
            # --- light grey lines and error bands ---
            if fermions == "quenched":
                LOW = flow_extr_filter_low_grey
                xdata = flow_var[LOW:]**2/8
                ydata = XX[LOW:,i]
                edata = XX_err[LOW:,i]
                
                #light grey backgrounds
                ax.fill_between(xdata, ydata-edata, ydata+edata, linewidth=0.5, color=(0.9,0.9,0.9,1), zorder=zorder)  #                
                #light grey lines for minimum width of background!
                ax.errorbar(xdata, ydata, color='lightgrey', zorder=zorder, fmt='-', lw=0.5, markersize=0.5)  
                
                #dark grey lines (0.7,0.7,0.7,1)
                ax.errorbar(flow_var**2/8, XX[:,i], linewidth=0.5, color=lpd.get_color(xdata_plot, i, 0, max_tau), zorder=zorder/100+1, fmt='-') #, label=number_format.format(xdata_plot[i]))
                #perturbative flow limit markers
                plots.append(ax.errorbar(flow_limit[i]**2/8, interpolation_value[i], marker=lpd.markers[i%len(lpd.markers)], fillstyle='none', markersize=4, mew=0.25, color=lpd.get_color(tauT, i), label='{0:.3f}'.format(tauT[i]), capsize=0, lw=0 ))
                ax.errorbar(flow_limit[i]**2/8, interpolation_value[i], marker='|', fillstyle='none', markersize=4, mew=0.25, color=lpd.get_color(tauT, i, 0, max_tau), capsize=0 )
        
            #dark grey lines that are cut off, may break perturbative flow limit markes (?)
            if fermions == "hisq": 
                xdata = flow_var[flow_extr_filter_low:]**2/8 # FIXME reset to flow_extr_filter_low_greyline
                ydata = XX[flow_extr_filter_low:,i] # FIXME reset to flow_extr_filter_low_greyline
                edata = XX_err[flow_extr_filter_low:,i] # FIXME reset to flow_extr_filter_low_greyline
            
                #colored lines
                ax.errorbar(xdata, ydata, linewidth=0.5, color=lpd.get_color(xdata_plot, i, 0, max_tau), zorder=zorder/100+1) 

            #vertical dashed line
            ax.axvline(x=flow_var[minimum_trusted_flowtime_index]**2/8, **lpd.verticallinestyle)
            
    ax2=ax.twiny()
    #new_tick_locations = numpy.array([0,0.05, 0.08, 0.1, 0.12, 0.13, 0.14, 0.15])**2/8
    
    #new_tick_locations = numpy.array([0,0.05**2/8, 0.08**2/8, 0.1**2/8, 0.12**2/8, 0.13**2/8, 0.14**2/8, 0.15**2/8])
    def tick_function(X):
        V = numpy.sqrt(X*8)
        return ["%.1f" % z if z == 0 else xlabel2_format % z for z in V ]
    ax2.tick_params(direction='in', pad=0, width=0.5)
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(new_tick_locations)
    ax2.set_xticklabels(tick_function(new_tick_locations))
    
    ax2.set_xlabel(xlabel2)#, horizontalalignment='right', verticalalignment='top', bbox=lpd.labelboxstyle, zorder=99)
    ax2.xaxis.set_label_coords(0.92,0.92)


    xerr_size=0
    if fermions == "hisq":
        xerr_size=0.3
    lpd.legendstyle.update(dict(loc="center right", bbox_to_anchor=(1,0.3), columnspacing=0.3, handletextpad=0.2, handlelength=0.5, handler_map={matplotlib.container.ErrorbarContainer: matplotlib.legend_handler.HandlerErrorbar(xerr_size=xerr_size)}))
    ax.legend(handles=plots)    
    handles, labels = ax.get_legend_handles_labels()
    leg = ax.legend(handles[::-1], labels[::-1], ncol=2, title=legend_title, **lpd.legendstyle)
    #FIXME remove the two lines below
    #for line in leg.get_lines():
        #line.set_linewidth(4.0)
    
    matplotlib.pyplot.tight_layout(0)
    ### change first xtick label
    #fig.canvas.draw()
    #ticks=ax.get_xticks().tolist()
    #ticks[0]='0.0' # FIXME remove this line
    #ax.set_xticklabels(ticks)
    filename = outputfolder+"/"+args.conftype+"_"+args.corr+prefix+"_flow.pdf"
    fig.savefig(filename) 
    print("saved flow effect plot", filename)


def main():
  
    #parse cmd line arguments
    parser, requiredNamed = lpd.get_parser()
   
    parser.add_argument('--flowselectionfile', help="only consider the flowtimes given in this file")
    parser.add_argument('--flowend', type=int, help="index of the maximum flow time")
    parser.add_argument('--flowstart', type=int, help="index of the minimum flow time", default='0')
    parser.add_argument('--flowstep', type=int, help="index step size for flow time", default='10')
    parser.add_argument('--part_obs', help='choose to only show numerator or denominator (polyakovloop) of corr', choices=['numerator','polyakovloop'])
    parser.add_argument('--conftype_2', help="if reconstruct, load data from this conftype and subtract the reconstructed corr from it and plot the result.")
    parser.add_argument('--show_flowtime_instead_of_flowradii', action="store_true")
    parser.add_argument('--testmode', action="store_true", help="turn on some temporary options (lims, labels, etc) when working on new data")
    parser.add_argument('--reconstruct', action="store_true", help="instead of plotting the actual correlator at T_1, use the data to plot the reconstructed correlator at T_2")
    parser.add_argument('--error_factor', type=float, default='1.0', help="reduce errors of all data by this factor")
    parser.add_argument('--show_TauByA', help="show tau/a instead of tau*T")
    parser.add_argument('--plot_for_peter', action="store_true", help="temporary option for incite proposal")
    parser.add_argument('--only_plot_no', type=int, nargs='*', default=[1,2,3])
    parser.add_argument('--tauT_selection', type=int, nargs='*', help='list of indices which tauT to plot in third flow plot')
   
   
    args = parser.parse_args()
   
    beta, ns, nt, nt_half = lpd.parse_conftype(args.conftype)
    fermions, temp, flowtype = lpd.parse_qcdtype(args.qcdtype)
   
    prefix = ""
    if args.part_obs:
        prefix = "_"+args.part_obs
    #if prefix == "_polyakovloop": # FIXME remove after calling reduce again
        #prefix = "polyakovloop"
        
    inputfolder=lpd.get_merged_data_path(args.qcdtype,args.corr,args.conftype)
    outputfolder=lpd.get_plot_path(args.qcdtype,args.corr,args.conftype)
           
    prefix_load = prefix
    if args.reconstruct:
        prefix_load = ""


    # ======= load data ======= 
    lpd.create_folder(outputfolder) 
    flow_radius = numpy.loadtxt(inputfolder+"flowradii_"+args.conftype+".dat")
    flow_times = numpy.loadtxt(inputfolder+"flowtimes_"+args.conftype+".dat")
    XX        = numpy.loadtxt(inputfolder+args.corr+prefix_load+"_"+args.conftype+".dat")
    XX_err = numpy.loadtxt(inputfolder+args.corr+prefix_load+"_err_"+args.conftype+".dat")
    
    tauT = lpd.get_tauTs(nt)
    xdata_plot = tauT
        

    #flow_selection = []
    #valid_flowtimes = None
    #if args.flowselectionfile:
        #valid_flowtimes = numpy.loadtxt(args.flowselectionfile)
        #for i,val in enumerate(flow_times):
            #if val in valid_flowtimes:
                #flow_selection.append(i)
        #last_entry = flow_selection[-1]
        #flow_selection = flow_selection[::10]
        #flow_selection.append(last_entry)
    #else:
        #flow_selection = range(0,flowend)


    


    valid_flowtimes = None
    #filter out flowtimes based on the file input
    if args.flowselectionfile:
        valid_flowtimes = numpy.loadtxt(args.flowselectionfile)
        delete_indices = []
        for i,val in enumerate(flow_times):
            if val not in valid_flowtimes:
                delete_indices.append(i)
        flow_times = numpy.delete(flow_times, delete_indices,0)
        flow_radius = numpy.delete(flow_radius, delete_indices,0)
        XX = numpy.delete(XX, delete_indices,0)
        XX_err = numpy.delete(XX_err, delete_indices,0)
    
    if args.plot_for_peter:
        global BB
        global BB_err 
        BB = numpy.loadtxt(lpd.get_merged_data_path(args.qcdtype,"BB",args.conftype)+"BB"+prefix_load+"_"+args.conftype+".dat")
        BB_err = numpy.loadtxt(lpd.get_merged_data_path(args.qcdtype,"BB",args.conftype)+"BB"+prefix_load+"_err_"+args.conftype+".dat")
        if args.flowselectionfile:
            BB = numpy.delete(BB, delete_indices,0)
            BB_err = numpy.delete(BB_err, delete_indices,0)
    
    flow_var = flow_radius
    
    flowend = args.flowend if args.flowend else len(flow_var)            
    flow_selection = range(args.flowstart,flowend,args.flowstep)
    
    flow_str = r"$ \sqrt{8\tau_\mathrm{F}}T$"
    if args.show_TauByA:
        xdata_plot *= nt    
    if args.show_flowtime_instead_of_flowradii:
        flow_var *= nt
        #flow_var = flow_times
        flow_str = r"$ \sqrt{8\tau_\mathrm{F}}/a$"
        
    
    multiplicity_correction = 1    
    # ONLY UNCOMMENT THESE WHEN USING OLD VERSIONS OF PARALLELGPUCODE WHERE MULTIPLICITIES OF DISCRETIZATIONS ARE NOT ACCOUNTED FOR
    #if corr == "EE_clover":
        #multiplicity_correction = 2
    #if corr == "BB_clover":
        #multiplicity_correction = -1.5
    #if corr == "BB":
        #multiplicity_correction = 1
            
        
    # ========================= PLOT: corr normalized to temperature, x-axis tauT, flowtimes fixed ========================= 
    
    nt_norm = nt
    if args.show_flowtime_instead_of_flowradii:
        nt_norm = 32
    nt_norm_half = int(nt_norm/2)
    
    XX_bak = numpy.copy(XX)
    XX_err_bak = numpy.copy(XX_err)
    
    if args.reconstruct:   
        #create reconstructed correlator 
        for i in range(len(flow_var)):
            for j in range(nt_norm_half): #from 0 to 16
                new_index = (64-(j+1+32))-1
                XX[i,j] += XX_bak[i,new_index]  
                XX_err[i,j] = numpy.sqrt( XX_err_bak[i,j]**2 + XX_err_bak[i,new_index]**2  )
        #load finite temp corr. then <finite temp corr> - <reconstructed corr>
        if args.conftype_2:    
            flow_times_2 = numpy.loadtxt(lpd.get_merged_data_path(args.qcdtype,args.corr,args.conftype_2)+"flowtimes_"+args.conftype_2+".dat")
            XX_2     = numpy.loadtxt(lpd.get_merged_data_path(args.qcdtype,args.corr,args.conftype_2)+args.corr+prefix_load+"_"+args.conftype_2+".dat")
            XX_2_err = numpy.loadtxt(lpd.get_merged_data_path(args.qcdtype,args.corr,args.conftype_2)+args.corr+prefix_load+"_err_"+args.conftype_2+".dat")
            
            #filter out flowtimes based on file
            delete_indices = []
            if args.flowselectionfile:
                valid_flowtimes = numpy.loadtxt(args.flowselectionfile)
                for i,val in enumerate(flow_times_2):
                    if val not in valid_flowtimes:
                        delete_indices.append(i)
                flow_times_2 = numpy.delete(flow_times_2, delete_indices,0)
                XX_2 = numpy.delete(XX_2, delete_indices,0)
                XX_2_err = numpy.delete(XX_2_err, delete_indices,0)     
            
            for i in range(len(flow_var)):
                for j in range(nt_norm_half): #from 0 to 16
                    XX[i,j] = XX_2[i,j] - XX[i,j] 
                    XX_err[i,j] += XX_2[i,j] 
    
    #normalize to temperature
    for i in range(len(flow_var)):
        for j in range(nt_half):
            XX[i,j] *= nt_norm**4 *multiplicity_correction
            XX_err[i,j] *= nt_norm**4 *multiplicity_correction
            if args.plot_for_peter:
                BB[i,j] *= nt_norm**4 *multiplicity_correction
                BB_err[i,j] *= nt_norm**4 *multiplicity_correction
    
            
    #TODO fix the third plot!!
            
    #reduce errors for plotting. default is to not do this. (args.error_factor=1)
    for i in range(len(flow_var)):
        for j in range(nt_half):
            XX_err[i,j] /= args.error_factor
            
    
    if 1 in args.only_plot_no:
        plot1(XX, XX_err, args,prefix, flow_selection,flow_var,xdata_plot, nt_norm_half, flow_str, outputfolder)
    
    # ========================= PLOT: corr normalized to LO pert corr, x-axis tauT, flowtimes fixed ========================= 
    
    if not prefix == "_polyakovloop":
        for i in range(len(flow_var)):
            for j in range(nt_half):
                XX[i,j] *= lpd.improve_corr_factor(j, nt_norm, i) / nt_norm**4
                XX_err[i,j] *= lpd.improve_corr_factor(j, nt_norm, i) / nt_norm**4
                if args.plot_for_peter:
                    BB[i,j] *= lpd.improve_corr_factor(j, nt_norm, i) / nt_norm**4
                    BB_err[i,j] *= lpd.improve_corr_factor(j, nt_norm, i) / nt_norm**4
    
    
    
    if 2 in args.only_plot_no:
        plot2(XX, XX_err, args,prefix, flow_selection,flow_var,xdata_plot, nt_norm_half, flow_str, outputfolder, fermions)

    # ========================= PLOT: x-axis flowtimes, tauT fixed ========================= 
    
    if 3 in args.only_plot_no:
        plot3(XX, XX_err, args,prefix, flow_var,xdata_plot, nt_norm_half, flow_str, outputfolder, fermions, nt_half, valid_flowtimes, tauT, flow_radius, flowend, flow_selection)

   
if __name__ == '__main__':
    main()
