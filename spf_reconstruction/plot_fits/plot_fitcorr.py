#!/usr/bin/env python3
import lib_process_data as lpd
import numpy
import argparse


def load_data(args):

    files = []
    for i, model_id in enumerate(args.model_ids):
        files.append(args.basepath + "/" + model_id + "/corrfit.dat")

    # === load data
    xdata = []
    ydata = []
    errorsleft = []
    errorsright = []
    for file in files:
        loadfunc = numpy.loadtxt
        try:
            data = loadfunc(file)  # TODO check if this load it correctly. may also need to load new omegaByT   file?
        except OSError:
            print("fail: could not find ", file, "\n")
            xdata.append(numpy.nan)
            ydata.append(numpy.nan)
            errorsleft.append(numpy.nan)
            errorsright.append(numpy.nan)
            continue

        xdata.append(data[:, 0])
        ydata.append((data[:, 3] / data[:, 1]))
        errorsleft.append(numpy.abs(ydata[-1]) * numpy.sqrt((data[:, 4] / data[:, 3]) ** 2 + (data[:, 2] / data[:, 1]) ** 2))
        errorsright.append(numpy.abs(ydata[-1]) * numpy.sqrt((data[:, 5] / data[:, 3]) ** 2 + (data[:, 2] / data[:, 1]) ** 2))

    return len(files), xdata, ydata, errorsleft, errorsright


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--model_ids', type=str, nargs='*')

    parser.add_argument('--xlims', nargs=2, type=float, default=(0.23, 0.54))
    parser.add_argument('--ylims', nargs=2, type=float, default=(0.75, 1.15))
    parser.add_argument('--yticks', nargs='*', default="auto")
    parser.add_argument('--xticks', nargs='*', default="auto")
    parser.add_argument('--labels', nargs='*', type=str)
    parser.add_argument('--colors', nargs='*')

    parser.add_argument('--basepath', type=str)
    parser.add_argument('--outputpath', type=str)
    parser.add_argument('--suffix', type=str)
    parser.add_argument('--corr', type=str, choices=["EE", "BB"])

    args = parser.parse_args()

    if len(args.model_ids) != len(args.labels):
        print("ERROR: need as many labels as model_ids")
        exit(1)

    return args


def main():

    args = parse_args()

    nfiles, xdata, ydata, errorsleft, errorsright = load_data(args)

    fig, ax, axtwiny = lpd.create_figure(xlabel=r'$\tau T$', ylabel=r'$\displaystyle \frac{G^\text{model}}{G'+lpd.get_corr_subscript(args.corr)+r'}$', xlims=args.xlims, ylims=args.ylims)

    if args.colors is None and args.corr == "BB":
        args.colors = [lpd.get_discrete_color(int(i)) for i in range(nfiles)]
    elif args.colors is None:
        args.colors = [lpd.get_discrete_color(int(i / 2)) for i in range(nfiles)]

    spacing = 0.002
    if args.corr == "BB":
        spacing = 0.004

    for i in range(nfiles):
        if not numpy.isnan(xdata[i]).any():

            if args.corr != "BB":
                if i % 2 == 0:
                    color = args.colors[i]
                else:
                    color = lpd.lighten_color(args.colors[i], 0.5)
            else:
                color = args.colors[i]
            ax.errorbar(xdata[i] + i * spacing, ydata[i], yerr=[errorsleft[i], errorsright[i]],
                        fmt='|', markersize=0, color=color, label=args.labels[i])
            # ax.errorbar(xdata[i] + i * 0.002, ydata[i],
            #             fmt='-', markersize=0, color=args.colors[i] if i % 2 == 0 else lpd.lighten_color(args.colors[i], 0.5), zorder=-100)
    ax.axhline(y=1, **lpd.horizontallinestyle)

    if args.yticks != "auto":
        ax.set_yticks([float(i) for i in args.yticks])
        axtwiny.set_yticks([float(i) for i in args.yticks])
    if args.xticks != "auto":
        ax.set_xticks([float(i) for i in args.xticks])

    ncols = 2
    if args.corr == "BB":
        ncols = 1

    if args.corr == "BB":
        ax.text(0.98, 0.99, r'\begin{align*} \displaystyle \frac{\bar{\mu}_T}{T}&=19.18 \\[-0.5ex] \displaystyle \frac{\bar{\mu}_{\tau_\mathrm{F}}}{\mu_\mathrm{F}}&=1.50 \end{align*}', transform=ax.transAxes, verticalalignment='top', horizontalalignment='right', fontsize=8)

    ax.legend(title="model", handlelength=1, loc="lower left", bbox_to_anchor=(0, 0), ncols=ncols, **lpd.leg_err_size(1, 0.3))

    lpd.create_folder(args.outputpath)
    outfile = args.outputpath + "/"+args.corr+"_corrfit" + args.suffix + ".pdf"
    fig.savefig(outfile)
    print("saved ", outfile)


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
