#!/usr/local/bin/python3.7 -u
import lib_process_data as lpd
import numpy
import argparse


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--test')
    args = parser.parse_args()

    inputfolder = "/work/home/altenkort/2piTD/"
    outputfolder = "/work/home/altenkort/2piTD/"

    EE = numpy.loadtxt("/work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/EE/EE_flow_extr.txt", unpack=True)
    BB = numpy.loadtxt("/work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/BB/BB_flow_extr.txt", unpack=True)

    plotstyle = dict(fmt='.', linewidth=1, markersize=0, capsize=2, mew=1, fillstyle='none')

    fig, ax, plots = lpd.create_figure(xlims=(0.2,0.51), ylims=(2.5, 4), xlabel=r'$\tau T$', xlabelpos=(0.9,0.08), ylabel=r'', ylabelpos=(0.1,0.9), UseTex=True,
                                       figsize=(2.75,2.5))

    ax.errorbar(EE[0], EE[1], EE[2], **plotstyle, label=r'$\displaystyle {G_E}/{G^\mathrm{norm}}$', color="k")
    ax.errorbar(BB[0], BB[1], BB[2], **plotstyle, label=r'$\displaystyle  {(G_B/Z^2_K)}/{G^\mathrm{norm}}$', color="C1")

    lpd.legendstyle.update(loc="upper left", bbox_to_anchor=(0, 1), fontsize=8, framealpha=0)

    leg = ax.legend(**lpd.legendstyle)

    ax.tick_params(axis='y', direction='out')

    ax.text(0.85, 0.95, 'preliminary', transform=ax.transAxes,
            fontsize=8, color='C1', alpha=1,
            ha='center', va='top', rotation='0', zorder=-1000000)

    filename = outputfolder + "/EEvsBB.pdf"
    fig.savefig(filename)
    print("saved correlator plot", filename)

    return


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
