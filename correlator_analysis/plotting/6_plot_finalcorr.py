#!/usr/bin/env python3

import numpy
import lib_process_data as lpd
import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_flow', type=str, required=True)
    parser.add_argument('--input_multilvl', type=str, required=True)
    parser.add_argument('--input_multilvl2', type=str, default="")
    parser.add_argument('--outputfolder', type=str, default="../../plots/")
    args = parser.parse_args()

    EE_final = numpy.loadtxt(args.input_flow, unpack=True)  #
    # EE_cont_2015 = numpy.loadtxt(args.input_multilvl)  # "../../data_merged/quenched/multi-level_2015/cont_thomas.dat"
    EE_cont_2015_new = numpy.loadtxt(args.input_multilvl)  #  inputfolder + "/multi-level_2015/EE_2015_new.txt")
    EE_cont_2015_new_2022 = None
    try:
        EE_cont_2015_new_2022 = numpy.loadtxt(args.input_multilvl2)
    except OSError:
        print("skip multilvl2")

    fig, ax, plots = lpd.create_figure(xlims=[0.05, 0.52], ylims=[1.5, 3.8], xlabel=r'$\tau T$',
                                       ylabel=r'$\displaystyle \frac{G}{G^{\mathrm{norm}}}$')

    ax.errorbar(EE_cont_2015_new[2:, 0], EE_cont_2015_new[2:, 1], EE_cont_2015_new[2:, 2],
                             label=r'Multi-level method', fmt='|')  # +"\n"+r'w/ $2016$ renorm. update'))#\beta_g (\mathcal{Z}_\mathrm{pert} -1)/6 \approx 0.138
    if EE_cont_2015_new_2022 is not None:
        ax.errorbar(EE_cont_2015_new_2022[2:, 0], EE_cont_2015_new_2022[2:, 1], EE_cont_2015_new_2022[2:, 2],
                label=r'Multi-level method', fmt='|')
    ax.errorbar(EE_final[0], EE_final[1], EE_final[2], fmt='|', label="Gradient flow method")
    # plots.append(ax.errorbar(lpd.tauT_30_ext[lpd.start:lpd.end], [k[1] for k in EE_final_old], [k[2] for k in EE_final_old], **lpd.plotstyle_add_point, label="Gradient flow method old", color=matplotlib.cm.gnuplot(0.8)))
    # plots.append(ax.errorbar(EE_cont_2015[:,0], EE_cont_2015[:,1], EE_cont_2015[:,2], 's-', label=r'Multi-level method$^1$'+"\n"+r'2015'))#\beta_g (\mathcal{Z}_\mathrm{pert} -1)/6 = 0.079

    ax.legend(loc='lower right', bbox_to_anchor=(1, 0.1), **lpd.leg_err_size())
    fig.savefig(args.outputfolder + "/EE_flowVSmultilvl.pdf")
    print("saved final corr plot", args.outputfolder + "/EE_flowVSmultilvl.pdf")


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
