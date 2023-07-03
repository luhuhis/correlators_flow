#!/usr/bin/env python3

import lib_process_data as lpd
import numpy
import argparse


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--outputfolder', default="/work/home/altenkort/work/correlators_flow/plots/quenched_1.50Tc_zeuthenFlow/")
    args = parser.parse_args()

    # TODO make this an input
    EE = numpy.loadtxt("/work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/EE/EE_flow_extr_relflow.txt", unpack=True)
    BB = numpy.loadtxt("/work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/BB/BB_flow_extr_relflow.txt", unpack=True)

    plotstyle = dict(fmt='.', markersize=0, fillstyle='none')

    fig, ax, plots = lpd.create_figure(xlims=(0.2, 0.51), ylims=(2.5, 4.1), xlabel=r'$\tau T$', ylabel=r'')

    ax.errorbar(BB[0], BB[1], BB[2], **plotstyle, label=r'$\displaystyle  {(G_B/Z^2_K)}/{G^\mathrm{norm}}$', color="C0")
    ax.errorbar(EE[0], EE[1], EE[2], **plotstyle, label=r'$\displaystyle {G_E}/{G^\mathrm{norm}}$', color="C1")

    ax.legend(loc="upper left", bbox_to_anchor=(0, 1), **lpd.leg_err_size())

    # ax.text(0.85, 0.95, 'preliminary', transform=ax.transAxes,
    #         fontsize=8, color='C1', alpha=1,
    #         ha='center', va='top', rotation='0', zorder=-1000000)

    filename = args.outputfolder + "/EEvsBB.pdf"
    fig.savefig(filename)
    print("saved correlator plot", filename)

    return


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
