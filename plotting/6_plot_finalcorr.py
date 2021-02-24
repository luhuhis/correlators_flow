#!/usr/local/bin/python

import numpy
import matplotlib
import lib_process_data as lpd

def main():

    qcdtype, fermions, temp, flowtype, corr, add_params = lpd.read_args(skip_conftype=True)

    inputfolder = lpd.get_merged_data_path(qcdtype,corr,"")

    #input
    EE_final = numpy.loadtxt(inputfolder+corr+"_final.txt", unpack=True)
    EE_cont_2015 = numpy.loadtxt("../../data_merged/quenched/multi-level_2015/cont_thomas.dat")
    EE_cont_2015_new = numpy.loadtxt(inputfolder+"/multi-level_2015/EE_2015_new.txt")

    #output
    #   lpd.outputfolder+"/EE_aeq0_teq0.pdf"


    fig, ax, plots = lpd.create_figure(xlims=[0.05,0.51], ylims=[1.5,3.8], xlabel=r'$\tau T$', ylabel=r'$\displaystyle \frac{G^\mathrm{cont }(\tau) }{G^{\substack{ \text{\tiny  norm} \\[-0.5ex] \text{\tiny cont } } } (\tau) }$ ', xlabelpos=(0.95,0.07), ylabelpos=(0.03,0.98)) # 
    lpd.legendstyle.update(dict(loc='lower right', bbox_to_anchor=(1,0.1), labelspacing=0.4, handler_map={matplotlib.container.ErrorbarContainer: matplotlib.legend_handler.HandlerErrorbar(xerr_size=0.5)} ))

    lpd.plotstyle_add_point.update(dict(fmt='s-'))
    plots.append(ax.errorbar(EE_final[0], EE_final[1], EE_final[2], **lpd.plotstyle_add_point, label="Gradient flow method")) #color=matplotlib.cm.gnuplot(0.5)
    #plots.append(ax.errorbar(EE_final_btstrp[0], EE_final_btstrp[1], EE_final_btstrp[2], **lpd.plotstyle_add_point, label="Gradient flow method btstrp"))
    #plots.append(ax.errorbar(lpd.tauT_30_ext[lpd.start:lpd.end], [k[1] for k in EE_final_old], [k[2] for k in EE_final_old], **lpd.plotstyle_add_point, label="Gradient flow method old", color=matplotlib.cm.gnuplot(0.8)))
    lpd.plotstyle_add_point.update(dict(fmt='o-'))
    plots.append(ax.errorbar(EE_cont_2015_new[2:,0], EE_cont_2015_new[2:,1], EE_cont_2015_new[2:,2], **lpd.plotstyle_add_point, label=r'Multi-level method'))#+"\n"+r'w/ $2016$ renorm. update'))#\beta_g (\mathcal{Z}_\mathrm{pert} -1)/6 \approx 0.138
    lpd.plotstyle_add_point.update(dict(fmt='s-'))
    #plots.append(ax.errorbar(EE_cont_2015[:,0], EE_cont_2015[:,1], EE_cont_2015[:,2], **lpd.plotstyle_add_point, label=r'Multi-level method$^1$'+"\n"+r'2015'))#\beta_g (\mathcal{Z}_\mathrm{pert} -1)/6 = 0.079 

    leg = ax.legend(handles=plots, title="", **lpd.legendstyle)
    leg.set_zorder(-1000)
    matplotlib.pyplot.tight_layout(0)
    ### change first xtick label
    fig.canvas.draw()
    ticks=ax.get_xticks().tolist()
    ticks = ['{0:.2f}'.format(x) for x in ticks]
    ticks[-2]=0.5
    #ticks[0]=0
    ax.set_xticklabels(ticks)
    fig.savefig(lpd.outputfolder+"/EE_aeq0_teq0.pdf") 
    print("saved final corr plot", lpd.outputfolder+"/EE_aeq0_teq0.pdf")
    
if __name__ == '__main__':
    main()
