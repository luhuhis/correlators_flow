#!/usr/local/bin/python

import lib_process_data as lpd

#input
#   XX_mean_finest = numpy.loadtxt(lpd.get_merged_data_path(qcdtype, conftypes[-1], corr)+"/"+corr+"_"+conftypes[-1]+".dat")
#   lpd.tauT(Ntau_finest)
#   lpd.get_merged_data_path+"/"+conftype+"/interpolations/EE_"+flowradius_str+"_interpolation.txt", for each conftype and flowradius

#output
#   lpd.inputfolder+"/cont_extr_quality/EE_"+flowradius_str+"_cont_quality.txt"
#   lpd.outputfolder+"/cont_extr_quality/EE_"+flowradius_str+"_cont_quality.pdf"
#   lpd.inputfolder+"/cont_extr/EE_"+flowradius_str+"_cont.txt"
#   lpd.outputfolder+"/cont_extr/EE_"+flowradius_str+"_cont.pdf"

import numpy
import matplotlib
import sys
import re
from latqcdtools import fitting
from latqcdtools import bootstr

def extrapolation_ansatz(x, m, b):
    return m * x + b

def fit_sample(ydata, xdata, edata, start_params=[80,3]):
    fitter = fitting.Fitter(func=extrapolation_ansatz, xdata=xdata, ydata=ydata, edata=edata, func_sup_numpy=True, always_return=False)
    fitparams, fitparams_err, chi_dof= fitter.try_fit(algorithms=["curve_fit"], start_params=start_params ) 
    return fitparams#, fitparams_err

def main():
    
    # parse cmd line arguments
    parser, requiredNamed = lpd.get_parser()
    parser.add_argument('--use_imp', help='whether to use tree-level improvement', type=bool, default=True)
    parser.add_argument('--nsamples', help="number of artifical gaussian bootstrap samples to generate", type=int, default=200)
    requiredNamed.add_argument('--flow_index', help='which flow time to interpolate', type=int, required=True)
    requiredNamed.add_argument('--conftypes', help="list of conftypes, e.g. s080t20_b0703500 s096t24_b0719200 s120t30_b0739400 s144t36_b0754400", nargs='*')
    parser.add_argument('--use_tex', type=bool, default=False)
    
    args = parser.parse_args()
   
    fermions, temp, flowtype = lpd.parse_qcdtype(args.qcdtype)
    
    flowindex = args.flow_index  
    flowindex = int(flowindex)
    i = flowindex #create alias
    flowradii = numpy.loadtxt(lpd.get_merged_data_path(args.qcdtype, args.corr, args.conftypes[-1])+"/flowradii_"+args.conftypes[-1]+".dat")
    flowradius = flowradii[i]
    flowradius_str = r'{0:.4f}'.format(flowradius)
    
    
    # --- set some params from args
    Ntaus = []
    for conftype in args.conftypes:
        dummy = re.sub(r'(^.*?)t', '', conftype)
        Ntaus.append(int(re.sub(r'(^.*?)_b(.*)', r'\1', dummy)))
    Ntau_coarsest = Ntaus[0]
    Ntau_finest = Ntaus[-1]
    print("perform extrapolation using Ntaus:", Ntaus)

    # --- checks
    if flowradius < 1/Ntaus[0]:
        exit("need at least a flow radius of 1/"+str(Ntau_coarsest)+" to have enough lattices for cont extr")
        
    
    # --- load finest lattice
    XX_mean_finest = numpy.loadtxt(lpd.get_merged_data_path(args.qcdtype, args.corr, args.conftypes[-1])+"/"+args.corr+"_"+args.conftypes[-1]+".dat")
    XX_err_finest = numpy.loadtxt(lpd.get_merged_data_path(args.qcdtype, args.corr, args.conftypes[-1])+"/"+args.corr+"_err_"+args.conftypes[-1]+".dat")
    
    # --- only use points beyond the lower limit
    XX_finest = numpy.asarray(
        ([tauT for tauT in lpd.get_tauTs(Ntau_finest) 
          if tauT > lpd.lower_tauT_limit(flowradius, Ntau_coarsest)],
         [val*lpd.improve_corr_factor(j, Ntau_finest, i, args.use_imp) for j,val in enumerate(XX_mean_finest[i]) 
          if lpd.get_tauTs(Ntau_finest)[j] > lpd.lower_tauT_limit(flowradius, Ntau_coarsest)],
         [val*lpd.improve_corr_factor(j, Ntau_finest, i, args.use_imp) for j,val in enumerate(XX_err_finest[i]) 
          if lpd.get_tauTs(Ntau_finest)[j] > lpd.lower_tauT_limit(flowradius, Ntau_coarsest)]))

    # --- load interpolations
    XX_ints = []
    for conftype in args.conftypes:
        ### ignore coarstest lattice!
        if conftype == args.conftypes[-1]:
            XX_ints.append(XX_finest)
        else:
            tmp = numpy.loadtxt(lpd.get_merged_data_path(args.qcdtype,args.corr,conftype)+"/interpolations/"+args.corr+"_"+flowradius_str+"_interpolation.txt", unpack=True)
            tmp2 = numpy.asarray(( [tauT for tauT in tmp[0] if tauT in lpd.get_tauTs(Ntau_finest)], 
                                    [val for j,val in enumerate(tmp[1]) if tmp[0][j] in lpd.get_tauTs(Ntau_finest)],
                                    [val for j,val in enumerate(tmp[2]) if tmp[0][j] in lpd.get_tauTs(Ntau_finest)] ))
            XX_ints.append(tmp2)
    #if XX_ints[1] is None:
        #exit("skip")
        
    ### define some parameters
    #interpolation_tauTs = XX_ints[1][0]
    n_tauTs = len(XX_finest[0])
    finest_ntau_half = int(Ntau_finest/2)
    offset = finest_ntau_half-n_tauTs
    valid_tauTs = XX_finest[0]

    results = numpy.empty((finest_ntau_half,3))
    results[:] = numpy.nan
    fitparameters = []
    mintauTindex =  None
    extrapolation_results = []
    flowstarts = lpd.get_flow_start_indices()
    maxtauTindex = 0
    for j in range(n_tauTs):
        if flowindex >= flowstarts[j]:
            maxtauTindex = j+offset+1
    maxtauTindex_plot = 0
    if args.use_tex:
        displaystyle = '\displaystyle'
        ylabel = r'$'+displaystyle+r'\frac{G^\mathrm{latt }_{\tau,\tau_\mathrm{F}} } { G^{\substack{ \text{\tiny  norm} \\[-0.4ex] \text{\tiny latt } } }_{\tau,\tau_\mathrm{F} = 0} } $'
    else:
        displaystyle = ''
        ylabel = 'G'
    ### plot settings
    ymin=2.55
    ymax=3.75
    fig, ax, plots = lpd.create_figure(xlims=[-0.0001,1/16**2+0.00015], ylims=[ymin, ymax], xlabel=r'$N_\tau^{-2}$', ylabel=ylabel, 
                                      xlabelpos=(0.95,0.07), ylabelpos=(0.08,0.98), UseTex = args.use_tex)
    lpd.titlestyle.update(dict(y=0.95))
    ax.set_title(r'$ \sqrt{8\tau_\mathrm{F}}T = '+flowradius_str+r'$', **lpd.titlestyle)

    plotpath=lpd.get_plot_path(args.qcdtype,args.corr,"")
    ceq_prefix=lpd.get_merged_data_path(args.qcdtype,args.corr,"")+"/cont_extr_quality/"+args.corr+"_"+flowradius_str
    lpd.create_folder(plotpath+"/cont_extr_quality/",
                      plotpath+"/cont_extr/",
                      lpd.get_merged_data_path(args.qcdtype,args.corr,"")+"/cont_extr_quality/", 
                      lpd.get_merged_data_path(args.qcdtype,args.corr,"")+"/cont_extr/")

    with open(ceq_prefix+"_cont_quality.txt", 'w') as outfile:
        outfile.write('# continuum extrapolation quality at fixed flowtime for all tauT of Ntau=36 \n')
        outfile.write('# tauT    1/Ntau^2    G/G_norm     err \n')
        with open(ceq_prefix+"_cont_quality_slope.txt", 'w') as outfile_slope:
            outfile_slope.write('# continuum extrapolation slope at fixed flowtime for all tauT of Ntau=36 \n')
            outfile_slope.write('# tauT     slope     slope_err \n')
        
            ### perform cont extr for each tauT at fixed flowtime
            for j,tauT in enumerate(lpd.get_tauTs(Ntau_finest)):
                xdata = numpy.asarray([1/Ntau**2 for k,Ntau in enumerate(Ntaus) if XX_ints[k] is not None])
                ### skip if flow time is too small
                if tauT in valid_tauTs and len(xdata) >=3: # and flowindex >= flowstarts[j]:
                    #print("working on flowtime ", flowradius_str, " tauT=,", '{0:.4f}'.format(tauT))
                    ### actually perform extr
                    print("index: ", j-offset)
                    for m,XX_int in enumerate(XX_ints):
                        print("XX_int no ", m, "of length", len(XX_int[1]))
                    
                    ydata = [ XX_int[1][j-offset] for XX_int in XX_ints if XX_int is not None]
                    edata = [XX_int[2][j-offset] for XX_int in XX_ints if XX_int is not None]
                    fitparams, fitparams_err = bootstr.bootstr_from_gauss(fit_sample, data=ydata, data_std_dev=edata, numb_samples = args.nsamples, sample_size = 1, return_sample=False, args=[xdata,edata] )
                    #fitparams, fitparams_err = fit_sample(ydata, xdata, edata)

                    results[j][0] = tauT
                    results[j][1] = fitparams[1]
                    results[j][2] = fitparams_err[1]
                    numpy.savetxt(outfile_slope, numpy.stack(([tauT], [fitparams[0]], [fitparams_err[0]]), axis=-1))
                
                    ### plot extrapolation
                    if mintauTindex is None:
                        mintauTindex = j
                    if j > maxtauTindex_plot:
                        maxtauTindex_plot = j
                    mycolor = lpd.get_color(lpd.get_tauTs(Ntau_finest), j, mintauTindex, maxtauTindex)
                    lpd.plotstyle_add_point_single.update(dict(fmt=lpd.markers[j-18+len(lpd.markers)]))
                    plots.append(ax.errorbar(xdata, ydata, edata, **lpd.plotstyle_add_point_single, color=mycolor, zorder=-100+j, label='{0:.3f}'.format(tauT)))
                    ax.errorbar(0, fitparams[1], fitparams_err[1], **lpd.plotstyle_add_point_single, color=mycolor, zorder=1)
                    x = numpy.linspace(0,0.1,100)
                    ax.errorbar(x, extrapolation_ansatz(x, *fitparams), color=mycolor, alpha=1, fmt=':', lw=0.5, zorder=-100)
                else:
                    results[j][0] = tauT
                    results[j][1] = numpy.nan
                    results[j][2] = numpy.nan
                    numpy.savetxt(outfile_slope, numpy.stack(([tauT], [numpy.nan], [numpy.nan]), axis=-1))
                    continue
                ### save extrapolation to file
                numpy.savetxt(outfile, numpy.stack(([tauT for dummy in range(len(xdata)+1)], [0, *xdata], [results[j][1], *ydata], [results[j][2], *edata]), axis=-1))
                outfile.write('# \n')
    #lpd.verticallinestyle.update(dict())#, dashes=(4,1)) )
    numpy.savetxt(lpd.get_merged_data_path(args.qcdtype,args.corr,"")+"/cont_extr/"+args.corr+"_"+flowradius_str+"_cont.txt", results, header="tauT    G/G_norm    err")

    ### save continuum extrapolation quality plot for this flow time
    lpd.plotstyle_add_point_single.update(dict(fmt='D-'))
    ax.axvline(x=0, ymin=((results[mintauTindex,1]-ymin)/(ymax-ymin)), ymax=((results[maxtauTindex_plot,1]-ymin)/(ymax-ymin)), alpha=1, color='grey', zorder=-1000, lw=0.5, dashes=(5,2))
    lpd.legendstyle.update(dict(loc="lower left", bbox_to_anchor=(-0.01,-0.01), columnspacing=0.1, labelspacing=0.25, handletextpad=0, borderpad=0, framealpha=0, handler_map={matplotlib.container.ErrorbarContainer: matplotlib.legend_handler.HandlerErrorbar(xerr_size=0.4)} ))
    ax.legend(handles=plots)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1], title=r'$\tau T=$', **lpd.legendstyle, ncol=3) #reverse ordering of legend
    matplotlib.pyplot.tight_layout(0)
    fig.savefig(plotpath+"/cont_extr_quality/"+args.corr+"_"+flowradius_str+"_cont_quality.pdf")
    ax.lines.clear() ; ax.collections.clear() ; plots.clear();

    ### save plot of continuum extrapolation for this flow time
    fig, ax, plots = lpd.create_figure(xlims=[0.15,0.51], ylims=[1.4, 4], xlabel=r'$\tau T$', ylabel=r'$'+displaystyle+r'\frac{G_{\tau_\mathrm{F}} (\tau)}{G^\mathrm{ norm }_{\tau_\mathrm{F} = 0} (\tau)}$', UseTex = args.use_tex)
    ax.set_title(r'$ \sqrt{8\tau_\mathrm{F}}T = $'+flowradius_str)
    ax.axvline(x=lpd.lower_tauT_limit(flowradius, Ntau_coarsest), **lpd.verticallinestyle)
    results = numpy.swapaxes(results, 0, 1)
    ax.errorbar(results[0], results[1], results[2], color="black", **lpd.plotstyle_add_point)
    matplotlib.pyplot.tight_layout(0.1)
    fig.savefig(plotpath+"/cont_extr/"+args.corr+"_"+flowradius_str+"_cont.pdf")

if __name__ == '__main__':
    main()
