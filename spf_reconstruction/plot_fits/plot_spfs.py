#!/usr/bin/env python3
import lib_process_data as lpd
import numpy
import argparse


def load_data(args):

    files = []
    for i, model_id in enumerate(args.model_ids):
        files.append(args.basepath + "/" + model_id + "/spffit.npy")

    # === load data
    xdata = []
    ydata = []
    errorsleft = []
    errorsright = []
    for file in files:
        loadfunc = numpy.load
        try:
            data = loadfunc(file)  # TODO check if this load it correctly. may also need to load new omegaByT   file?
        except OSError:
            print("fail: could not find ", file, "\n")
            xdata.append([numpy.nan,])
            ydata.append([numpy.nan,])
            errorsleft.append([numpy.nan,])
            errorsright.append([numpy.nan,])
            continue

        xdata.append(data[:, 0])
        ydata.append(data[:, 1])
        errorsleft.append(data[:, 2])
        errorsright.append(data[:, 3])

    return len(files), xdata, ydata, errorsleft, errorsright


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--model_ids', type=str, nargs='*')

    parser.add_argument('--xlims', nargs=2, type=float, default=(0.01, 90))
    parser.add_argument('--ylims', nargs=2, type=float, default=(1, 1300))
    parser.add_argument('--labels', nargs='*', type=str)
    parser.add_argument('--colors', nargs='*')

    parser.add_argument('--basepath', type=str)
    parser.add_argument('--outputpath', type=str)
    parser.add_argument('--suffix', type=str)
    parser.add_argument('--corr', type=str, choices=["EE", "BB"])
    parser.add_argument('--leg_pos', nargs=2, default=(0.45, 1))
    parser.add_argument('--leg_loc', default="upper center")
    parser.add_argument('--plot_spf_err', action="store_true")

    args = parser.parse_args()

    if len(args.model_ids) != len(args.labels):
        print("ERROR: need as many labels as model_ids")
        exit(1)

    return args


def main():

    args = parse_args()

    nfiles, xdata, ydata, errorsleft, errorsright = load_data(args)

    fig, ax, plots = lpd.create_figure(xlabel=r'$\omega/T$', ylabel=r'$\displaystyle \frac{2\rho}{ \omega T^2}$', xlims=args.xlims, ylims=args.ylims)

    ax.set_yscale('log')
    ax.set_xscale('log')

    if args.colors is None:
        args.colors = [lpd.get_discrete_color(int(i/2)) for i in range(nfiles)]

    for i in range(nfiles):
        if not numpy.isnan(xdata[i]).any():
            factor = 2/xdata[i]
            y = ydata[i]*factor
            if args.plot_spf_err:
                yerr = [errorsleft[i]*factor, errorsright[i]*factor]
                ax.fill_between(xdata[i][::10], y[::10]-yerr[0][::10], y[::10]+yerr[1][::10], facecolor=args.colors[i], alpha=0.1, zorder=-100-i)

            fmt = '--' if i % 2 == 0 else ':'

            ax.errorbar(xdata[i], y, fmt=fmt, label=args.labels[i], color=args.colors[i], zorder=-i)
            # ax.axvline(x=1, **lpd.verticallinestyle)
            # ax.axvline(x=omegaUV, **lpd.verticallinestyle)

    ax.legend(loc=args.leg_loc, bbox_to_anchor=args.leg_pos, title="model", handlelength=1, fontsize=10)

    lpd.create_folder(args.outputpath)
    outfile = args.outputpath + "/"+args.corr+"_spf" + args.suffix + ".pdf"
    fig.savefig(outfile)
    print("saved ", outfile)


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
