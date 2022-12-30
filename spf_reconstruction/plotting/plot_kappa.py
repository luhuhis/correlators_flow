#!/usr/bin/env python3
import lib_process_data as lpd
import numpy
import argparse


def load_data(args):

    files = []
    for i, model_id in enumerate(args.model_ids):
        files.append(args.basepath + "/" + model_id + "/params.dat")

    # === load data
    xdata = []
    ydata = []
    errorsleft = []
    errorsright = []
    chisqdof_arr = []
    chisqdof_arr_err = []
    for file in files:
        loadfunc = numpy.loadtxt
        try:
            data = loadfunc(file)
        except OSError:
            print("fail: could not find ", file, "\n")
            xdata.append(numpy.nan)
            ydata.append(numpy.nan)
            errorsleft.append(numpy.nan)
            errorsright.append(numpy.nan)
            chisqdof_arr.append(numpy.nan)
            chisqdof_arr_err.append(numpy.nan)
            continue

        idx = 0
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
    parser.add_argument('--temperature_on_xaxis', default=False, action="store_true")
    parser.add_argument('--pos', nargs='*', default=None, type=float, help="position of each kappa value. this could be either just some arbitrary number or the temperature for example.")
    parser.add_argument('--colors', nargs='*')
    parser.add_argument('--outputpath', type=str)
    parser.add_argument('--suffix', type=str)
    parser.add_argument('--corr', type=str, choices=["EE", "BB"])
    parser.add_argument('--scale_error_by_chisqdof', default=False, action="store_true")
    parser.add_argument('--figsize', default=None, nargs='*')
    args = parser.parse_args()

    if len(args.model_ids) != len(args.labels):
        print("ERROR: need as many pos as model_ids")
        exit(1)

    if args.ylims is None and args.pos is None:
        args.ylims = (0, len(args.model_ids)+2.5)

    if args.pos is None:
        args.pos = numpy.flip(numpy.linspace(0.5, len(args.model_ids)+1, len(args.model_ids)))

    return args


def main():

    args = parse_args()

    nfiles, xdata, ydata, errorsleft, errorsright, chisqdof_arr, chisqdof_arr_err = load_data(args)

    fig, ax, plots = lpd.create_figure(figsize=args.figsize, xlims=args.xlims, ylims=args.ylims)

    if args.colors is None:
        args.colors = [lpd.get_discrete_color(int(i/2)) for i in range(nfiles)]

    ax.set_xlabel(r'$\kappa/T^3$', horizontalalignment='left',  verticalalignment='bottom', bbox=lpd.labelboxstyle)
    ax.xaxis.set_label_coords(0.01, 0.01)

    leftmost = 100
    rightmost = 0
    for i in range(nfiles):
        if not numpy.isnan(xdata[i]).any():
            if args.temperature_on_xaxis:
                factor = (args.pos[i] / 1000) ** 3
                x = args.pos[i]
                y = xdata[i] * factor
                ax.errorbar(x, y, yerr=[[errorsleft[i] * factor], [errorsright[i] * factor]], color=args.colors[i], fmt='|')
                thisleft = y-errorsleft[i] * factor
                if leftmost > numpy.min(thisleft):
                    leftmost = thisleft
                thisright = y+errorsright[i] * factor
                if rightmost < numpy.max(thisright):
                    rightmost = thisright
            else:
                x = xdata[i]
                y = args.pos[i]
                scale_error_factor = 1
                if args.scale_error_by_chisqdof:
                    scale_error_factor = numpy.fmax(1, chisqdof_arr[i])
                ax.errorbar(x, y, xerr=[[errorsleft[i]*scale_error_factor], [errorsright[i]*scale_error_factor]],
                            markersize=2, color=args.colors[i], fmt='D', fillstyle='full')
                args.labels[i] += r', $'+lpd.format_float(chisqdof_arr[i], 1)+r'\pm'+lpd.format_float(chisqdof_arr_err[i], 1)+r'$'
                args.labels[i] = r'\begin{flushleft}'+args.labels[i]+r'\end{flushleft}'
                thisleft = x-errorsleft[i] * scale_error_factor
                if leftmost > numpy.min(thisleft):
                    leftmost = thisleft
                thisright = x+errorsright[i] * scale_error_factor
                if rightmost < numpy.max(thisright):
                    rightmost = thisright

    ax.fill_between([leftmost, rightmost], [0, 0], [100, 100], facecolor='k', alpha=0.2, zorder=-1000)
    centralvalue = (leftmost+rightmost)/2
    error = rightmost-centralvalue
    numpy.savetxt("", [centralvalue, error])

    if args.temperature_on_xaxis:
        ax.set_xlabel(r'T\mathrm{[MeV]}')
        ax.set_ylabel(r'\kappa \mathrm{[GeV}^3\mathrm{]}')
    else:
        ax.yaxis.set_label_coords(1.01, 1)
        ax.set_ylabel(r'model, $\chi^2/\mathrm{dof}$', horizontalalignment='left', verticalalignment='top')
        ax.set_yticks(args.pos)
        ax.yaxis.tick_right()
        ax.set_yticklabels(args.labels)
        if args.xticks != "auto":
            ax.set_xticks([float(i) for i in args.xticks])
        ax.yaxis.label.set_size(10)
        ax.tick_params(axis='y', which='major', labelsize=10)

    lpd.create_folder(args.outputpath)
    outfile = args.outputpath + "/"+args.corr+"_kappa" + args.suffix + ".pdf"
    fig.savefig(outfile)
    print("saved ", outfile)


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
