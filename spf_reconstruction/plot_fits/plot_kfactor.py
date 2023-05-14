#!/usr/bin/env python3
import lib_process_data as lpd
import numpy
import argparse


def load_data(args):

    files = []
    for i, model_id in enumerate(args.model_ids):
        files.append(args.basepath + "/" + model_id)

    # === load data
    xdata = []
    ydata = []
    errorsleft = []
    errorsright = []
    chisqdof_arr = []
    chisqdof_arr_err = []
    for file in files:
        try:
            data = numpy.loadtxt(file+"/params.dat")
        except OSError as error:
            print(error)
            xdata.append(numpy.nan)
            ydata.append(numpy.nan)
            errorsleft.append(numpy.nan)
            errorsright.append(numpy.nan)
            chisqdof_arr.append(numpy.nan)
            chisqdof_arr_err.append(numpy.nan)
            continue

        idx = 1

        xdata.append(data[idx, 0])

        chisqdof = data[-1, 0]
        chisqdof_arr.append(chisqdof)
        chisqdof_arr_err.append(numpy.fmax(data[-1, 1], data[-1, 2]))

        errorsleft.append(data[idx, 1])
        errorsright.append(data[idx, 2])

    return len(files), xdata, ydata, errorsleft, errorsright, chisqdof_arr, chisqdof_arr_err


def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('--xlims', nargs=2, type=float)
    parser.add_argument('--xticks', nargs="*", default="auto")
    parser.add_argument('--ylims', nargs=2, type=float, default=None)
    parser.add_argument('--model_ids', type=str, nargs='*')
    parser.add_argument('--labels', nargs='*', type=str)
    parser.add_argument('--basepath', type=str)
    parser.add_argument('--fillstyles', nargs='*')
    parser.add_argument('--fmts', nargs='*')
    parser.add_argument('--temperature_on_xaxis', default=False, action="store_true")
    parser.add_argument('--pos', nargs='*', default=None, type=float, help="position of each kfactor value. this could be either just some arbitrary number or the temperature for example.")
    parser.add_argument('--colors', nargs='*')
    parser.add_argument('--outputpath', required=True)
    parser.add_argument('--outputpath_data', required=True)
    parser.add_argument('--suffix', type=str, default="")
    parser.add_argument('--corr', type=str, choices=["EE", "BB"])
    parser.add_argument('--scale_error_by_chisqdof', default=False, action="store_true")
    parser.add_argument('--figsize', default=None, nargs='*')
    parser.add_argument('--Tc_in_GeV', type=float)

    args = parser.parse_args()

    if len(args.model_ids) != len(args.labels):
        print("ERROR: need as many labels as model_ids")
        exit(1)

    if args.ylims is None and args.pos is None:
        args.ylims = (-1, len(args.model_ids)+2.5)

    if args.pos is None:
        args.pos = numpy.flip(numpy.linspace(0.5, len(args.model_ids)+1, len(args.model_ids)))

    return args


def main():

    args = parse_args()

    nfiles, xdata, ydata, errorsleft, errorsright, chisqdof_arr, chisqdof_arr_err = load_data(args)
    labelboxstyle = dict(boxstyle="Round", fc="None", ec="None", alpha=0, pad=0.1, zorder=99)
    fig, ax, plots = lpd.create_figure(figsize=args.figsize, xlims=args.xlims, ylims=args.ylims, ylabelbox=labelboxstyle, xlabelbox=labelboxstyle)

    if args.colors is None:
        args.colors = [lpd.get_discrete_color(int(i/2)) for i in range(nfiles)]

    ax.yaxis.set_label_coords(0.01, 0.99)
    ax.xaxis.set_label_coords(0.99, 0.01)
    ax.set_xlabel(r'$T/T_c$', horizontalalignment='right', verticalalignment='bottom', bbox=lpd.labelboxstyle)
    ax.set_ylabel(r'$K$', horizontalalignment='left', verticalalignment='top', rotation=0, bbox=lpd.labelboxstyle)

    for i in range(nfiles):
        x = args.pos[i]/args.Tc_in_GeV
        y = xdata[i]
        if y is not numpy.nan:
            ax.errorbar(x, y, yerr=[[errorsleft[i]], [errorsright[i]]], color=args.colors[i],
                    label=args.labels[i], markersize=4, lw=0.75, mew=0.75,
                    fmt='.' if args.fmts is None else args.fmts[i],
                    fillstyle=None if args.fillstyles is None else args.fillstyles[i],
            )
        else:
            ax.errorbar(0, 0, label=args.labels[i], markersize=0, alpha=0, lw=0)

    leg = ax.legend(**lpd.leg_err_size(1, 0.5), ncol=2, columnspacing=1, labelspacing=0.5, loc="lower left", bbox_to_anchor=(0.25, 0), framealpha=0)
    for t in leg.texts:
        t.set_multialignment('left')

    ax.text(0.99, 0.99, r'$\rho_\mathrm{model}=\rho_\mathrm{``smax"}$', ha='right', va='top', transform=ax.transAxes, zorder=-1000)
    ax.text(0.23, 0.18, r'$\rho_\mathrm{UV}\! :$', ha='right', va='top', transform=ax.transAxes, zorder=-1000)

    lpd.create_folder(args.outputpath)
    outfile = args.outputpath + "/"+args.corr+"_kfactor" + args.suffix + ".pdf"
    fig.savefig(outfile)
    print("saved ", outfile)


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
