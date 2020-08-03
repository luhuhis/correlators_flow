#!/usr/local/bin/python
import numpy
import matplotlib
import sys
from latqcdtools import fitting
from latqcdtools import bootstr
import lib_plot as lp
import lib_process_data as lpd
import _3_spline_interpolate

def extrapolation_ansatz(x, m, b):
    return m * x + b

def fit_sample(ydata, xdata, edata, start_params=[80,3]):
    fitter = fitting.Fitter(func=extrapolation_ansatz, xdata=xdata, ydata=ydata, edata=edata, func_sup_numpy=True, always_return=False)
    fitparams, fitparams_err, chi_dof= fitter.try_fit(algorithms=["curve_fit"], start_params=start_params ) 
    return fitparams#, fitparams_err

def main():
    
    use_imp = True
    
    try:
        flowindex = sys.argv[1]
        flowindex = int(flowindex)
        i = flowindex #create alias
        flowradius = lp.flow_radius[i]
        flowradius_str = r'{0:.3f}'.format(flowradius)
        if flowradius < 1/20:
            exit("need at least a flow radius of 1/20 to have enough lattices for cont extr")
    except IndexError:
        exit("Invalid Arguments. Usage: script.py <flowindex>, e.g. script.py 100")
    

    ### load largest lattice
    EE36 = numpy.asarray(([tauT for tauT in lpd.tauT_imp[36] if tauT > _3_spline_interpolate.lower_tauT_limit(flowradius)], 
                        [val*lpd.improve_corr_factor(j, 36, i, use_imp) for j,val in enumerate(lp.EE_36[i]) if lpd.tauT_imp[36][j] > _3_spline_interpolate.lower_tauT_limit(flowradius)], [val*lpd.improve_corr_factor(j, 36, i, use_imp) for j,val in enumerate(lp.EE_err_36[i]) if lpd.tauT_imp[36][j] > _3_spline_interpolate.lower_tauT_limit(flowradius)]))



    ### load interpolations
    EE_ints = []
    for conftype in ("s064t16_b0687361", "s080t20_b0703500", "s096t24_b0719200", "s120t30_b0739400", "s144t36_b0754400"):
        ### ignore coarstest lattice!
        if conftype == "s064t16_b0687361":
            EE_ints.append(None)
            continue
        if conftype == "s144t36_b0754400":
            EE_ints.append(EE36)
        else:
            tmp = numpy.loadtxt(lp.inputfolder+"/"+conftype+"/interpolations/EE_"+flowradius_str+"_interpolation.txt", unpack=True)
            tmp2 = numpy.asarray(( [tauT for tauT in tmp[0] if tauT in lpd.tauT_imp[36]], 
                                    [val for j,val in enumerate(tmp[1]) if tmp[0][j] in lpd.tauT_imp[36]],
                                    [val for j,val in enumerate(tmp[2]) if tmp[0][j] in lpd.tauT_imp[36]] ))
            EE_ints.append(tmp2)
    #if EE_ints[1] is None:
        #exit("skip")
    
    ### define some variables
    #interpolation_tauTs = EE_ints[1][0]
    nsamples = 10000
    Ntaus = numpy.asarray([16,20,24,30,36])
    n_tauTs = len(EE36[0])
    offset = 18-n_tauTs
    valid_tauTs = EE36[0]

    results = numpy.empty((18,3))
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
    UseTex = True
    if UseTex:
        displaystyle = '\displaystyle'
        ylabel = r'$'+displaystyle+r'\frac{G^\mathrm{latt }_{\tau,\tau_F} } { G^{\substack{ \text{\tiny  norm} \\[-0.4ex] \text{\tiny latt } } }_{\tau,\tau_F = 0} } $'
    else:
        displaystyle = ''
        ylabel = 'G'
    ### plot settings
    ymin=2.55
    ymax=3.75
    fig, ax, plots = lp.create_figure(xlims=[-0.0001,1/20**2+0.00015], ylims=[ymin, ymax], xlabel=r'$N_\tau^{-2}$', ylabel=ylabel, 
                                      xlabelpos=(0.95,0.07), ylabelpos=(0.08,0.98), UseTex = UseTex)
    lp.titlestyle.update(dict(y=0.95))
    ax.set_title(r'$ \sqrt{8\tau_F}T = '+flowradius_str+r'$', **lp.titlestyle)
    lpd.create_folder(lp.outputfolder+"/cont_extr_quality/")


    with open(lp.inputfolder+"/cont_extr_quality/EE_"+flowradius_str+"_cont_quality.txt", 'w') as outfile:
        outfile.write('# continuum extrapolation quality at fixed flowtime for all tauT of Ntau=36 \n')
        outfile.write('# tauT    1/Ntau^2    G/G_norm     err \n')
        
        ### perform cont extr for each tauT at fixed flowtime
        for j,tauT in enumerate(lpd.tauT_imp[36]):
            print("working on ", tauT)
            xdata = numpy.asarray([1/Ntau**2 for k,Ntau in enumerate(Ntaus) if EE_ints[k] is not None])
            
            ### skip if there is a lattice missing or flow time is too small
            if tauT in valid_tauTs and len(xdata) >=4 and flowindex >= flowstarts[j]:
                ### actually perform extr
                ydata = [EE_int[1][j-offset] for EE_int in EE_ints if EE_int is not None]
                edata = [EE_int[2][j-offset] for EE_int in EE_ints if EE_int is not None]
                fitparams, fitparams_err = bootstr.bootstr_from_gauss(fit_sample, data=ydata, data_std_dev=edata, numb_samples = nsamples, sample_size = 1, return_sample=False, args=[xdata,edata] )
                #fitparams, fitparams_err = fit_sample(ydata, xdata, edata)

                results[j][0] = tauT
                results[j][1] = fitparams[1]
                results[j][2] = fitparams_err[1]
            
                ### plot extrapolation
                if mintauTindex is None:
                    mintauTindex = j
                if j > maxtauTindex_plot:
                    maxtauTindex_plot = j
                mycolor = lp.get_color(lpd.tauT_imp[36], j,mintauTindex,maxtauTindex)
                lp.plotstyle_add_point_single.update(dict(fmt=lp.markers[j-18+len(lp.markers)]))
                plots.append(ax.errorbar(xdata, ydata, edata, **lp.plotstyle_add_point_single, color=mycolor, zorder=-100+j, label='{0:.3f}'.format(tauT)))
                ax.errorbar(0, fitparams[1], fitparams_err[1], **lp.plotstyle_add_point_single, color=mycolor, zorder=1)
                x = numpy.linspace(0,0.1,100)
                ax.errorbar(x, extrapolation_ansatz(x, *fitparams), color=mycolor, alpha=1, fmt=':', lw=0.5, zorder=-100)
            else:
                results[j][0] = tauT
                results[j][1] = numpy.nan
                results[j][2] = numpy.nan
                continue
            ### save extrapolation to file
            numpy.savetxt(outfile, numpy.stack(([tauT for dummy in range(len(xdata)+1)], [0, *xdata], [results[j][1], *ydata], [results[j][2], *edata]), axis=-1))
            outfile.write('# \n')
    #lp.verticallinestyle.update(dict())#, dashes=(4,1)) )
    numpy.savetxt(lp.inputfolder+"/cont_extr/EE_"+flowradius_str+"_cont.txt", results, header="tauT    G/G_norm    err")

    ### save continuum extrapolation quality plot for this flow time
    lp.plotstyle_add_point_single.update(dict(fmt='D-'))
    ax.axvline(x=0, ymin=((results[mintauTindex,1]-ymin)/(ymax-ymin)), ymax=((results[maxtauTindex_plot,1]-ymin)/(ymax-ymin)), alpha=1, color='grey', zorder=-1000, lw=0.5, dashes=(5,2))
    lp.legendstyle.update(dict(loc="lower left", bbox_to_anchor=(-0.01,-0.01), columnspacing=0.1, labelspacing=0.25, handletextpad=0, borderpad=0, framealpha=0, handler_map={matplotlib.container.ErrorbarContainer: matplotlib.legend_handler.HandlerErrorbar(xerr_size=0.4)} ))
    ax.legend(handles=plots)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1], title=r'$\tau T=$', **lp.legendstyle, ncol=3) #reverse ordering of legend
    lpd.create_folder(lp.outputfolder+"/cont_extr_quality/")
    matplotlib.pyplot.tight_layout(0)
    fig.savefig(lp.outputfolder+"/cont_extr_quality/EE_"+flowradius_str+"_cont_quality.pdf")
    ax.lines.clear() ; ax.collections.clear() ; plots.clear();

    ### save plot of continuum extrapolation for this flow time
    fig, ax, plots = lp.create_figure(xlims=[0.15,0.51], ylims=[1.4, 4], xlabel=r'$\tau T$', ylabel=r'$'+displaystyle+r'\frac{G_{\tau_F} (\tau)}{G^\mathrm{ norm }_{\tau_F = 0} (\tau)}$', UseTex = UseTex)
    ax.set_title(r'$ \sqrt{8\tau_F}T = $'+flowradius_str)
    ax.axvline(x=_3_spline_interpolate.lower_tauT_limit(flowradius), **lp.verticallinestyle)
    results = numpy.swapaxes(results, 0, 1)
    ax.errorbar(results[0], results[1], results[2], color="black", **lp.plotstyle_add_point)
    lpd.create_folder(lp.outputfolder+"/cont_extr/")
    matplotlib.pyplot.tight_layout(0.1)
    fig.savefig(lp.outputfolder+"/cont_extr/EE_"+flowradius_str+"_cont.pdf")

if __name__ == '__main__':
    main()
