#!/usr/local/bin/python

import lib_process_data as lpd
import numpy
import matplotlib
from matplotlib import pyplot

a = 0.028057  # lattice spacing in fm from fk scale

a = 1

def convert_dimless_flowtime_to_dimful_flowradius(t):

    return numpy.sqrt(8*t)*a

def main():

    path = "/work/home/altenkort/work/correlators_flow/data/merged/hisq_b8249_zeuthenFlow/EE/"

    EE_32 = numpy.loadtxt(path+"/s096t32_b0824900_m002022_m01011/EE_s096t32_b0824900_m002022_m01011.dat")
    EE_32_err = numpy.loadtxt(path+"/s096t32_b0824900_m002022_m01011/EE_s096t32_b0824900_m002022_m01011.dat")

    EE_64 = numpy.loadtxt(path+"/s064t64_b0824900_m002022_m01011/EE_s064t64_b0824900_m002022_m01011.dat")
    EE_64_err = numpy.loadtxt(path+"/s064t64_b0824900_m002022_m01011/EE_s064t64_b0824900_m002022_m01011.dat")

    flowtimes = numpy.loadtxt(path+"/s064t64_b0824900_m002022_m01011/flowtimes_s064t64_b0824900_m002022_m01011.dat")

    fig, ax, plots = lpd.create_figure(ylims=[4e-6, 5e-5], xlabel=r'$\tau [\mathrm{fm}]$', xlabelpos=(0.95, 0.07), ylabelpos=(0.04, 0.98), UseTex=True,
                                       ylabel=r'$\displaystyle {G^\mathrm{latt}} $',)

    ax.set_yscale('log', nonposy='clip')

    lpd.plotstyle_add_point_single.update(dict(fmt='-', mew=0.5))

    fidxs = (67, 82, 95)

    err_reduce = 0.05
    ax.set_title("!!! Errors scaled down by factor "+str(int(1/err_reduce))+" !!!")

    coloroffset = 2
    n = len(fidxs) + coloroffset
    colorinput = numpy.linspace(0, 1, n)

    tau32 = numpy.asarray(range(1, 17))
    tau64 = numpy.asarray(range(1, 33))

    for i, fidx in enumerate(fidxs):
        flowradius = convert_dimless_flowtime_to_dimful_flowradius(flowtimes[fidx])
        # i *= 2
        ax.errorbar(tau32 * a*1.005, EE_32[fidx, :], EE_32_err[fidx, :]*err_reduce,
                **lpd.chmap(lpd.plotstyle_add_point, fmt='-', color=lpd.get_color(colorinput, i+coloroffset, 0, n-1), label=r'$32;'+'{0:.2f}'.format(flowradius)+r'$'))
        ax.errorbar(tau64 * a*0.995, EE_64[fidx, :], EE_64_err[fidx, :]*err_reduce,
                **lpd.chmap(lpd.plotstyle_add_point, fmt='-', color=lpd.get_color2(colorinput, i+coloroffset, 0, n - 1), label=r'$64;'+'{0:.2f}'.format(flowradius)+r'$'))

    lpd.legendstyle.update(dict(loc="upper right", bbox_to_anchor=(1.01, 0.95)))
    ax.legend(title=r'$N_\tau; \sqrt{8\tau_\mathrm{F}}$', **lpd.legendstyle)
    # ax.set_title(r'BB', x=0.5, y=0.89)

    fig.savefig("/work/home/altenkort/work/correlators_flow/plots/hisq_b8249_zeuthenFlow/EE_hisq_finitevszero.pdf")


if __name__ == '__main__':
    main()

