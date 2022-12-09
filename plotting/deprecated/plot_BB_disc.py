#!/usr/local/bin/python

import lib_process_data as lpd
import numpy
import matplotlib
from matplotlib import pyplot


def main():

    BB05 = numpy.loadtxt("/work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/BB/cont_extr/BB_0.0508_cont.txt", unpack=True)
    BB10 = numpy.loadtxt("/work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/BB/cont_extr/BB_0.1000_cont.txt", unpack=True)
    BBclov05 = numpy.loadtxt("/work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/BB_clover/cont_extr/BB_clover_0.0508_cont.txt", unpack=True)
    BBclov10 = numpy.loadtxt("/work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/BB_clover/cont_extr/BB_clover_0.1000_cont.txt", unpack=True)
    BBt0 = numpy.loadtxt("/work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/BB/BB_final.txt", unpack=True)
    BBclovt0 = numpy.loadtxt("/work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/BB_clover/BB_clover_final.txt", unpack=True)

    fig, ax, plots = lpd.create_figure(xlims=[0.2, 0.51], ylims=(2.4,3.7), xlabel=r'$\tau T$', xlabelpos=(0.95, 0.07), ylabelpos=(0.08, 0.98), UseTex=True,
                                       ylabel=r'$\displaystyle \frac{G^\mathrm{cont}}{ G^{\substack{ \text{\tiny  norm} \\[-0.4ex] \text{\tiny cont } } }_{\tau_\mathrm{F} = 0} }$',)

    lpd.plotstyle_add_point_single.update(dict(fmt='-'))
    ax.errorbar(BBclov10[0], BBclov10[1], BBclov10[2], color="orange", **lpd.plotstyle_add_point, label=r'clover,\ 0.10')
    ax.errorbar(BB10[0], BB10[1], BB10[2], color="green", **lpd.plotstyle_add_point, label=r'naive,\ \ 0.10')
    lpd.plotstyle_add_point_single.update(dict(fmt='-'))
    ax.errorbar(BBclov05[0], BBclov05[1], BBclov05[2], color="red", **lpd.plotstyle_add_point, label=r'clover,\ 0.05', zorder=-10)
    ax.errorbar(BB05[0], BB05[1], BB05[2], color="violet", **lpd.plotstyle_add_point, label=r'naive,\ \ 0.05', zorder=-10)
    lpd.plotstyle_add_point_single.update(dict(fmt='-'))
    ax.errorbar(BBclovt0[0], BBclovt0[1], BBclovt0[2], color="blue", **lpd.plotstyle_add_point, label=r'clover,\ $\tau_\mathrm{F}\rightarrow$ 0')
    ax.errorbar(BBt0[0], BBt0[1], BBt0[2], color="black", **lpd.plotstyle_add_point, label=r'naive,\ \ $\tau_\mathrm{F}\rightarrow$ 0')

    lpd.legendstyle.update(dict(loc="lower right", bbox_to_anchor=(1.01, 0.07)))
    ax.legend(title=r'\quad discr., $\sqrt{8\tau_\mathrm{F}}T$', **lpd.legendstyle)
    ax.set_title(r'BB', x=0.5, y=0.89)

    matplotlib.pyplot.tight_layout(0.1)
    fig.savefig("/work/home/altenkort/work/correlators_flow/plots/quenched_1.50Tc_zeuthenFlow/BB_disc_comparison.pdf")


if __name__ == '__main__':
    main()
