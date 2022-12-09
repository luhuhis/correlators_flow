#!/usr/bin/env python3
import lib_process_data as lpd
import argparse


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--outputfolder', type=str, required=True)
    args = parser.parse_args()

    plotstyle = dict(fmt='.', markersize=0)

    fig, ax, plots = lpd.create_figure(xlims=(1, 2.5), ylims=(0, 3.5), xlabel=r'$T/T_c$', xlabelpos=(0.9, 0.08), ylabel=r'', ylabelpos=(0.1, 0.9))

    banerjee = [[1.2, 1.8, 0.8], [1.5, 1.55, 0.55], [2.0, 1.2, 0.6]]

    brambilla = [1.52, 1.82, 0.79]

    altenkort = [1.48, 1.07, 0.43]

    ax.errorbar([x[0] for x in banerjee], [y[1] for y in banerjee], [yerr[2] for yerr in banerjee], color='C0', **plotstyle,
                 label=r'$\kappa_B/T^3$, Banerjee et al., \textquotesingle 22')
    ax.errorbar(brambilla[0], brambilla[1], brambilla[2], **plotstyle, color='C2', label=r'$\kappa_B/T^3$, Brambilla et al., \textquotesingle 22')
    ax.errorbar(altenkort[0], altenkort[1], altenkort[2], **plotstyle, color='C1',
                 label=r'$\kappa_B/(T^3Z_K^2)$, \textbf{Altenkort et al.}')

    ax.legend(loc="upper left", bbox_to_anchor=(0, 1))

    # ax.text(0.5/1.5, 0.1, 'preliminary', transform=ax.transAxes,
    #         fontsize=8, color='C1', alpha=1,
    #         ha='center', va='center', rotation='0', zorder=-1000000)

    filename = args.outputfolder + "/kappa_B.pdf"
    fig.savefig(filename)
    print("saved correlator plot", filename)

    return


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
