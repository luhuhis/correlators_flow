#!/usr/local/bin/python

    #input
    #   lpd.inputfolder+conftype+"/interpolations/XX_"+'{0:.3f}'.format(lpd.flow_radius[i])+"_interpolation.txt"
    #   lpd.inputfolder+"/cont_extr/XX_"+'{0:.3f}'.format(lpd.flow_radius[i])+"_cont.txt"

    #output
    #   lpd.outputfolder+"/lattice_effects/XX_flow_"+'{0:.3f}'.format(lpd.flow_radius[i])+".pdf"

import sys
import numpy
import matplotlib
import itertools
import lib_process_data as lpd
from matplotlib import cm


def main():

    try:
        qcdtype, fermions, temp, flowtype, corr, add_params = lpd.read_args(skip_conftype=True)
        flowindex = add_params[0]
        conftypes = add_params[1:]
    except:
        exit("Invalid Arguments. Usage: script.py <qcdtype> <corr> <flowindex> <conftype1> <conftype2> <...>, e.g. script.py quenched_1.50Tc_zeuthenFlow XX s080t20_b0703500 s096t24_b0719200 s120t30_b0739400 s144t36_b0754400")

    i = int(flowindex) #create alias
    flow_radius = numpy.loadtxt(lpd.get_merged_data_path(qcdtype,corr,conftypes[0])+"flowradii_"+conftypes[0]+".dat")
    flowradius = flow_radius[i]

    outputfolder = lpd.get_plot_path(qcdtype, corr, "")+"/lattice_effects/"

    use_imp = True

    UseTex = True
    if UseTex:
        displaystyle = r'\displaystyle'
        ylabel = r'$ '+displaystyle+r'\frac{G^\mathrm{ latt }_{\tau_\mathrm{F}} (\tau)}{G_{\tau_\mathrm{F} = 0}^{\substack{ \text{\tiny  norm} \\[-0.4ex] \text{\tiny latt } } } (\tau)} $'
    else:
        ylabel = r'G'
        displaystyle = r''
        
    fig, ax, plots = lpd.create_figure(xlims=[0.15,0.51], ylims=[2,4], xlabel=r'$\tau T$', ylabel=ylabel, xlabelpos=(0.95,0.05), ylabelpos=(0.03,0.98), UseTex = UseTex)

    lpd.titlestyle.update(dict(y=0.95))
    lpd.legendstyle.update(dict(loc="lower right", bbox_to_anchor=(1.01,0.1), framealpha=0.5, handlelength=1, handler_map={matplotlib.container.ErrorbarContainer: matplotlib.legend_handler.HandlerErrorbar(xerr_size=0.4)}))


    ### load lattice data and interpolations
    XX_data = []
    XX = []
    XX_err = []
    tauT = []
    Ntaus = []
    for conftype in conftypes: 
        inputfolder = lpd.get_merged_data_path(qcdtype,corr,conftype)
        beta, ns, nt, nt_half = lpd.parse_conftype(conftype)
        Ntaus.append(nt)
        XX_data.append(numpy.loadtxt(inputfolder+"/interpolations/"+corr+"_"+'{0:.3f}'.format(flow_radius[i])+"_interpolation.txt", unpack=True))
        XX.append(numpy.loadtxt(inputfolder+"/"+corr+"_"+conftype+".dat"))
        XX_err.append(numpy.loadtxt(inputfolder+"/"+corr+"_err_"+conftype+".dat"))
        tauT.append(numpy.arange(1/nt,0.501,1/nt))
        
    labels = [r'$N_\tau = '+str(nt)+r'$' for nt in Ntaus] 
    ### also load continuum extr if possible
    try:
        XX_data.append(numpy.loadtxt(lpd.get_merged_data_path(qcdtype,corr,"")+"/cont_extr/"+corr+"_"+'{0:.3f}'.format(flow_radius[i])+"_cont.txt", unpack=True))
        labels.append(r'cont.')
    except:
        XX_data.append(None)
        
    ### some plot settings
    colors=cm.tab10(numpy.linspace(0,1,10))
    plots = []
    ax.set_title(r'$\sqrt{8\tau_\mathrm{F}} T=$ '+'{0:.3f}'.format(flow_radius[i]), **lpd.titlestyle)
    markers = ['s', 'o', 'D', 'H', 'h']
    zorders = range(-30,0,4)

    ### loop over all the different data sets and plot settings
    for XX_int, tauTs, XX, XX_err, label, color, marker, zorder in itertools.zip_longest(XX_data, tauT, XX, XX_err, labels, colors, markers, zorders):
        lpd.plotstyle_add_point.update(dict(fmt=marker))
        if tauTs is not None: #tauTs is None for the continuum limit
            nt = int(len(tauTs)*2)
        if XX_int is not None:
            ### remove nans
            columns = []
            nan_array = numpy.isnan(XX_int[1])
            for l,column in enumerate(XX_int):
                columns.append(column[~nan_array])
            ### plot original Data (XX, XX_err), interpolations (XX_int) and continuum extrapolation
            if tauTs is None: #continuum plot
                plots.append(ax.errorbar(columns[0], columns[1], columns[2], color=color, label=label, **lpd.plotstyle_add_point, zorder=zorder+1))
                lpd.plotstyle_add_point.update(dict(fmt='-'))
                ax.errorbar(columns[0], columns[1], columns[2], color=color, **lpd.plotstyle_add_point, zorder=zorder+1)
            else: #interpolation plot
                ax.errorbar(columns[0], columns[1], fmt = '-', lw = 0.5, mew=0, color=color, zorder=zorder+1)   
                plots.append(ax.errorbar(tauTs, [val * lpd.improve_corr_factor(v, nt, i, use_imp) for v,val in enumerate(XX[i,:])], [val * lpd.improve_corr_factor(v, nt, i, use_imp) for v,val in enumerate(XX_err[i,:])], **lpd.plotstyle_add_point, color=color, label=label, zorder=zorder+2))
    #if len(plots) < 4:
        #exit("Error: need at least 4 lattices")
    ax.legend(handles=plots)
    handles, labels = ax.get_legend_handles_labels()
    leg = ax.legend(handles[::-1], labels[::-1], title="", **lpd.legendstyle)
    leg.set_zorder(-1)
    for line in leg.get_lines():
        line.set_linewidth(4.0)
    if i != 0:
        ax.axvline(x=lpd.lower_tauT_limit(flowradius), **lpd.verticallinestyle)
    matplotlib.pyplot.tight_layout(0.2)
    lpd.create_folder(outputfolder)
    fig.savefig(outputfolder+"/"+corr+"_flow_"+'{0:.3f}'.format(flow_radius[i])+".pdf")#, **lpd.figurestyle) 
    print("saved lattice effects plot", outputfolder+"/"+corr+"_flow_"+'{0:.3f}'.format(flow_radius[i])+".pdf") 

if __name__ == '__main__':
    main()
