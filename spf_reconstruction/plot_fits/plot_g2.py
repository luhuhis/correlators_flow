#!/usr/bin/env python3
import lib_process_data as lpd
import numpy
import argparse

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--PhiUV_files', nargs='*')
    parser.add_argument('--labels', nargs='*')
    parser.add_argument('--thresholds', nargs='*', type=float)
    parser.add_argument('--outputpath')
    parser.add_argument('--suffix')
    parser.add_argument('--xlims', nargs=2, type=float)
    parser.add_argument('--ylims', nargs=2, type=float)
    parser.add_argument('--leg_loc', type=str)
    parser.add_argument('--leg_pos', nargs=2, type=float)
    parser.add_argument('--fmt', nargs='*')
    args = parser.parse_args()
    return args


def main():

    args = parse_args()

    fig, ax, plots = lpd.create_figure(xlabel=r'$\omega/T$', ylabel=r'$g^2(\mu)$', xlims=args.xlims, ylims=args.ylims)

    # ax.set_yscale('log')
    ax.set_xscale('log')

    colors = [lpd.get_discrete_color(int(i / 1)) for i in range(len(args.PhiUV_files))]

    Nc=3
    Cf=(Nc ** 2 - 1) / 2 / Nc

    fmts=['.', ':', '-.', '--', '.', '.', ':', '-.']

    dummycounter=0

    handlers = []

    for i, file in enumerate(args.PhiUV_files):
        if file != "dummy":
            data = numpy.load(file)
            xdata = data[:, 0]
            ydata = data[:, 1]
            factor = 1/xdata**3 / Cf *6 * numpy.pi
            # factor = 2 / xdata
            ydata *= factor
            tmp = ax.errorbar(xdata, ydata, fmt=fmts[i], label=args.labels[i], color=colors[i-dummycounter], alpha=1, zorder=-i)
            handlers.append(tmp)
        else:
            tmp = ax.errorbar(1, 1, fmt='.', markersize=0, label=args.labels[i])
            handlers.append(tmp)
            dummycounter+=1

    leg1 = ax.legend(handles=handlers[:4], loc=args.leg_loc, bbox_to_anchor=args.leg_pos, handlelength=1.4, fontsize=8)#, title=r'$N_f, T/T_c, \mu$')  #

    leg2 = ax.legend(handles=handlers[5:], loc="upper right", bbox_to_anchor=(1,1), handlelength=1.4, fontsize=8)#, title=r'$N_f, T/T_c, \mu$')  #

    ax.add_artist(leg1)

    # ax.axvline(x=numpy.pi)
    # ax.axvline(x=2*numpy.pi)


    lpd.create_folder(args.outputpath)
    outfile = args.outputpath + "/UV_spf" + args.suffix + ".pdf"
    fig.savefig(outfile)
    print("saved ", outfile)


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
