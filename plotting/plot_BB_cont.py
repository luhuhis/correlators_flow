#!/usr/local/bin/python

import lib_process_data as lpd
import numpy
import matplotlib
from matplotlib import pyplot


def main():

    BB05 = numpy.loadtxt("/work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/BB/cont_extr/BB_0.0508_cont.txt", unpack=True)
    BB07 = numpy.loadtxt("/work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/BB/cont_extr/BB_0.0750_cont.txt", unpack=True)
    BB10 = numpy.loadtxt("/work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/BB/cont_extr/BB_0.1000_cont.txt", unpack=True)
    BB12 = numpy.loadtxt("/work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/BB/cont_extr/BB_0.1200_cont.txt", unpack=True)
    BB15 = numpy.loadtxt("/work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/BB/cont_extr/BB_0.1400_cont.txt", unpack=True)
    fig, ax, plots = lpd.create_figure(xlims=[0.2, 0.51], ylims=(2.35, 3.7), xlabel=r'$\tau T$', xlabelpos=(0.95, 0.07), ylabelpos=(0.04, 0.98), UseTex=True,
                                       ylabel=r'$\displaystyle \frac{G^\mathrm{cont}}{ G^{\substack{ \text{\tiny  norm} \\[-0.4ex] \text{\tiny cont } } }_{\tau_\mathrm{F} = 0} }$',)

    n = 5
    colorinput = numpy.linspace(0, 1, n)
    lpd.plotstyle_add_point_single.update(dict(fmt='-'))

    ax.errorbar(BB15[0], BB15[1], BB15[2], **lpd.chmap(lpd.plotstyle_add_point, fmt='-', color=lpd.get_color(colorinput, 0, 0, n-1), label=r'0.140'))
    ax.errorbar(BB12[0], BB12[1], BB12[2], **lpd.chmap(lpd.plotstyle_add_point, fmt='-', color=lpd.get_color(colorinput, 1, 0, n-1), label=r'0.120'))
    ax.errorbar(BB10[0], BB10[1], BB10[2], **lpd.chmap(lpd.plotstyle_add_point, fmt='-', color=lpd.get_color(colorinput, 2, 0, n-1), label=r'0.100'))
    ax.errorbar(BB07[0], BB07[1], BB07[2], **lpd.chmap(lpd.plotstyle_add_point, fmt='-', color=lpd.get_color(colorinput, 3, 0, n-1), label=r'0.075'))
    ax.errorbar(BB05[0], BB05[1], BB05[2], **lpd.chmap(lpd.plotstyle_add_point, fmt='-', color=lpd.get_color(colorinput, 4, 0, n-1), label=r'0.051'))

    lpd.legendstyle.update(dict(loc="lower right", bbox_to_anchor=(1.01, 0.07)))
    ax.legend(title=r'$\sqrt{8\tau_\mathrm{F}}T=$', **lpd.legendstyle)
    # ax.set_title(r'BB', x=0.5, y=0.89)

    fig.savefig("/work/home/altenkort/work/correlators_flow/plots/quenched_1.50Tc_zeuthenFlow/BB_cont_comparison.pdf")


if __name__ == '__main__':
    main()

