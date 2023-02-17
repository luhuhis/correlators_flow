#!/usr/bin/env python3
import lib_process_data as lpd
import numpy
import matplotlib.ticker as mtick
import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputfolder', type=str, default="../../data/merged/")
    parser.add_argument('--outputfolder', type=str, default="../../plots/")
    args = parser.parse_args()

    lpd.create_folder(args.outputfolder)
    xdata, I2, I3 = numpy.loadtxt(args.inputfolder+"EE_QED_LPT.dat", unpack=True)
    fig, ax, plots = lpd.create_figure(xlabel=r'$\sqrt{8\tau_\mathrm{F}}/a$', ylabel=r'', xlabelpos=(0.99, 0.07))
    ax.errorbar(numpy.sqrt(xdata*8), I2/I2[0]*100, label=r'operator mixing')  # $-a^4\frac{3}{2}\mathcal{I}_2(\tau_\mathrm{F})$
    ax.errorbar(numpy.sqrt(xdata*8), I3/I3[-1]*100, fmt=':', label=r'flow time renorm.')  # $a^4\frac{1}{4}\mathcal{I}_3(\tau_\mathrm{F})$
    ax.axhline(y=100, **lpd.horizontallinestyle)
    ax.axhline(y=0, **lpd.horizontallinestyle)
    ax.yaxis.set_major_formatter(mtick.PercentFormatter())
    ax.legend(loc="center right", bbox_to_anchor=(1, 0.46), framealpha=0, handlelength=1)  #, title="QED: lattice artifacts \n at NLO")
    fig.savefig(args.outputfolder+"/EE_QED_LPT.pdf")
    print("saved QED LPT plot", args.outputfolder+"/EE_QED_LPT.pdf")


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
