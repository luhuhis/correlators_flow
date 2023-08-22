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
    BB_UVLO_IRNLO = numpy.loadtxt("/work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/BB/BB_flow_extr_relflow_UVLO_IRNLO.txt", unpack=True)
    BB_UVNLO_IR_NLO = numpy.loadtxt(
        "/work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/BB/BB_flow_extr_relflow_UVNLO_IRNLO.txt",
        unpack=True)
    BB_UVLO_IRLO = numpy.loadtxt("/work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/BB/BB_flow_extr_relflow_UVLO_IRLO.txt", unpack=True)
    BB_UVNLO_IRLO = numpy.loadtxt(
        "/work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/BB/BB_flow_extr_relflow_UVNLO_IRLO.txt",
        unpack=True)

    plotstyle = dict(fmt='.', markersize=0, fillstyle='none')

    fig, ax, plots = lpd.create_figure(xlims=(0.2, 0.51), ylims=(2.5, 5.1), xlabel=r'$\tau T$', ylabel=r'')

    ax.errorbar(BB_UVLO_IRNLO[0], BB_UVLO_IRNLO[1], BB_UVLO_IRNLO[2], **plotstyle,
                label=r'$  {G_B^\text{phys.}}/{G^\mathrm{norm}}, \bar{\mu}_\text{UV}/\mu_\text{F}={\phantom{0}1.00},\ \bar{\mu}_\text{IR}/T={19.18}$', color="C0")
    ax.errorbar(BB_UVLO_IRLO[0], BB_UVLO_IRLO[1], BB_UVLO_IRLO[2], **plotstyle,
                label=r'$  {G_B^\text{phys.}}/{G^\mathrm{norm}}, \bar{\mu}_\text{UV}/\mu_\text{F}={\phantom{0}1.00},\ \bar{\mu}_\text{IR}/T={\phantom{0}6.28}$',
                color="C2")
    ax.errorbar(BB_UVNLO_IR_NLO[0], BB_UVNLO_IR_NLO[1], BB_UVNLO_IR_NLO[2], **plotstyle,
                label=r'$  {G_B^\text{phys.}}/{G^\mathrm{norm}}, \bar{\mu}_\text{UV}/\mu_\text{F}={10.13},\ \bar{\mu}_\text{IR}/T={19.18}$', color="C1")
    ax.errorbar(BB_UVNLO_IRLO[0], BB_UVNLO_IRLO[1], BB_UVNLO_IRLO[2], **plotstyle,
                label=r'$  {G_B^\text{phys.}}/{G^\mathrm{norm}}, \bar{\mu}_\text{UV}/\mu_\text{F}={10.13},\ \bar{\mu}_\text{IR}/T={\phantom{0}6.28}$',
                color="C3")
    ax.errorbar(EE[0], EE[1], EE[2], **plotstyle, label=r'$\displaystyle {G_E}/{G^\mathrm{norm}}$', color="C4")

    ax.legend(loc="upper left", bbox_to_anchor=(0, 1), **lpd.leg_err_size(), fontsize=6)

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
