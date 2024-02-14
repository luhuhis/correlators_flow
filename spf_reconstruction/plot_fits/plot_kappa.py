#!/usr/bin/env python3
import lib_process_data as lpd
import numpy as np
import argparse
from matplotlib.ticker import AutoMinorLocator


use_pgf_backend = True
backend = 'pgf'

def load_data(args):

    files = []
    for i, model_id in enumerate(args.model_ids):
        files.append(args.basepath + "/" + model_id)

    kappa_samples = []

    # === load data
    xdata = []
    ydata = []
    errorsleft = []
    errorsright = []
    chisqdof_arr = []
    chisqdof_arr_err = []
    for file in files:
        try:
            data = np.loadtxt(file+"/params.dat")
            fitparam_samples = np.load(file + "/params_samples.npy")
        except OSError as error:
            print(error)
            xdata.append(np.nan)
            ydata.append(np.nan)
            errorsleft.append(np.nan)
            errorsright.append(np.nan)
            chisqdof_arr.append(np.nan)
            chisqdof_arr_err.append(np.nan)
            continue

        idx = 0

        kappa_samples.append(fitparam_samples[:, idx])

        xdata.append(data[idx, 0])

        chisqdof = data[-1, 0]
        chisqdof_arr.append(chisqdof)
        chisqdof_arr_err.append(np.fmax(data[-1, 1], data[-1, 2]))

        errorsleft.append(data[idx, 1])
        errorsright.append(data[idx, 2])

    kappa_samples = np.concatenate([*kappa_samples])
    kappa_median = np.nanmedian(kappa_samples)
    kappa_errs = lpd.dev_by_dist(kappa_samples, return_both_q=True, percentile=68)
    kappa_bounds = [kappa_median-kappa_errs[0], kappa_median+kappa_errs[1]]

    return len(files), xdata, ydata, errorsleft, errorsright, chisqdof_arr, chisqdof_arr_err, kappa_bounds


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
    parser.add_argument('--outputpath', required=True)
    parser.add_argument('--outputpath_data', required=True)
    parser.add_argument('--suffix', type=str, default="")
    parser.add_argument('--corr', type=str, choices=["EE", "BB"])
    parser.add_argument('--scale_error_by_chisqdof', default=False, action="store_true")
    parser.add_argument('--figsize', default=None, nargs='*')
    parser.add_argument('--no_subscript', action="store_true")
    parser.add_argument('--hide_chisq', action="store_true")
    args = parser.parse_args()

    if len(args.model_ids) != len(args.labels):
        print("ERROR: need as many labels as model_ids")
        exit(1)

    if args.ylims is None and args.pos is None:
        args.ylims = (-1, len(args.model_ids)+3)

    if args.pos is None:
        args.pos = np.flip(np.linspace(0, len(args.model_ids)+1, len(args.model_ids)))

    return args


def add_group_labels(ax):
    num_groups = 4
    group_size = 3
    total_points = num_groups * group_size
    y_positions = np.flip(np.linspace(0, total_points + 1, total_points))

    # Assuming you have a list of distinct labels for each group
    group_labels = [r'\parbox[t]{0.85cm}{\raggedleft $2\pi,$}\ \  \parbox[t]{0.5cm}{$1.0$}',
                    r'\parbox[t]{0.85cm}{\raggedleft $2\pi,$}\ \  \parbox[t]{0.5cm}{$1.5$}',
                    r'\parbox[t]{0.85cm}{\raggedleft $19.18,$}\ \ \parbox[t]{0.5cm}{$1.0$}',
                    r'\parbox[t]{0.85cm}{\raggedleft $19.18,$}\ \ \parbox[t]{0.5cm}{$1.5$}']
                    # r'\parbox{2cm}{ \begin{align*}  \displaystyle \bar{\mu}_{\tau_\mathrm{F}}/\mu_\mathrm{F}&=1.5 \\[-1ex] \bar{\mu}_T/T&=19.18 \end{align*}}']

    # Adjust the loop to use these labels and draw separation lines
    for group in range(num_groups):
        # Calculate the middle position of each group for the label
        start_index = group * group_size
        end_index = start_index + group_size
        middle_pos = np.mean(y_positions[start_index:end_index])-0.2

        # Use the distinct label for each group
        label = group_labels[group]
        ax.text(1.27, middle_pos, label, transform=ax.get_yaxis_transform(), verticalalignment='center', horizontalalignment='left', fontsize=8)

        # Draw dashed horizontal lines to separate groups, except after the last group
        if group < num_groups - 1:
            line_pos = (y_positions[end_index - 1] + y_positions[end_index]) / 2
            ax.axhline(y=line_pos, color='k', linestyle=':', lw=0.4) #, xmax=3, clip_on = False)


def main():

    args = parse_args()

    nfiles, xdata, ydata, errorsleft, errorsright, chisqdof_arr, chisqdof_arr_err, kappa_bounds = load_data(args)
    labelboxstyle = dict(boxstyle="Round", fc="None", ec="None", alpha=0, pad=0.1, zorder=99)
    fig, ax, plots = lpd.create_figure(figsize=args.figsize, xlims=args.xlims, ylims=args.ylims, ylabelbox=labelboxstyle, xlabelbox=labelboxstyle, ytwinticks=False, minorticks=False, use_pgf_backend=use_pgf_backend)

    if args.colors is None and args.corr == "BB":
        args.colors = []
        for i in range(int(nfiles/3)):
            args.colors.append(lpd.get_discrete_color(1))
            args.colors.append(lpd.get_discrete_color(2))
            args.colors.append(lpd.get_discrete_color(3))
    elif args.colors is None:
        args.colors = [lpd.get_discrete_color(int(i/2)) for i in range(nfiles)]

    if args.no_subscript:
        corr = ""
    else:
        corr = args.corr

    ax.set_xlabel(r'$\kappa'+lpd.get_corr_subscript(corr)+r'/T^3$', horizontalalignment='left',  verticalalignment='bottom', bbox=None)
    ax.xaxis.set_label_coords(0.01, 0.01)

    leftmost = 100
    rightmost = 0
    for i in range(nfiles):
        if not np.isnan(xdata[i]).any():
            if args.temperature_on_xaxis:
                factor = (args.pos[i] / 1000) ** 3
                x = args.pos[i]
                y = xdata[i] * factor
                ax.errorbar(x, y, yerr=[[errorsleft[i] * factor], [errorsright[i] * factor]], color=args.colors[i], fmt='|')
                thisleft = y-errorsleft[i] * factor
                if leftmost > np.min(thisleft):
                    leftmost = thisleft
                thisright = y+errorsright[i] * factor
                if rightmost < np.max(thisright):
                    rightmost = thisright
            else:
                x = xdata[i]
                y = args.pos[i]
                scale_error_factor = 1
                if args.scale_error_by_chisqdof:
                    scale_error_factor = np.fmax(1, chisqdof_arr[i])
                ax.errorbar(x, y, xerr=[[errorsleft[i]*scale_error_factor], [errorsright[i]*scale_error_factor]],
                            markersize=2, color=args.colors[i], fmt='D', fillstyle='full')
                if not args.hide_chisq:
                    args.labels[i] = r'\parbox{3cm}{' + args.labels[i] + r' \hfill $' + lpd.format_float(chisqdof_arr[i], 1)+r'\pm'+lpd.format_float(chisqdof_arr_err[i], 1)+r'$\null}'
                # else:
                #     r'\parbox{5cm}{' + args.labels[i] + r' \hfill \null}'
                # args.labels[i] = r'\begin{flushleft}'+args.labels[i]+r'\end{flushleft}'
                thisleft = x-errorsleft[i] * scale_error_factor
                if leftmost > np.min(thisleft):
                    leftmost = thisleft
                thisright = x+errorsright[i] * scale_error_factor
                if rightmost < np.max(thisright):
                    rightmost = thisright

    leftmost, rightmost = kappa_bounds
    ax.fill_between([leftmost, rightmost], [-100, -100], [100, 100], facecolor='grey', alpha=0.2, zorder=-1000)
    centralvalue = (leftmost+rightmost)/2
    error = rightmost-centralvalue
    filename = args.outputpath_data + "/" + args.corr + "_kappa" + args.suffix + ".txt"
    with open(filename, 'w') as file:
        file.write('# kappa +- error \n')
        file.write("# $"+lpd.format_float(centralvalue, 1)+r' \pm '+str(lpd.float_ceil(error, 1))+"$\n")
        file.write("# $[" + lpd.format_float(leftmost, 3) + r', ' + lpd.format_float(rightmost, 3) + "]$\n")
        np.savetxt(file, np.c_[centralvalue, error])
    print("saved", filename)

    if args.temperature_on_xaxis:
        ax.set_xlabel(r'T\mathrm{[MeV]}')
        ax.set_ylabel(r'\kappa \mathrm{[GeV}^3\mathrm{]}')
    else:
        ax.yaxis.set_label_coords(1.04, 1.00)
        chisqstr = ""
        if not args.hide_chisq:
            chisqstr = r'$\chi^2/\mathrm{d.o.f.}$'
            ax.set_ylabel(r'\parbox{5.25cm}{$\mathrm{model}$ \hfill ' + chisqstr + r'\null}', horizontalalignment='left', verticalalignment='top',
                          fontsize=5)
        elif args.hide_chisq and corr == "BB":
            ax.set_ylabel(r'\raisebox{0pt}[0pt][0pt]{$\mathrm{model,}$\ \ \, $\scalebox{0.7}{$\frac{\bar{\mu}_T}{T}$},$}\ \scalebox{0.7}{$\frac{\bar{\mu}_{\tau_\mathrm{F}}}{\mu_\mathrm{F}}$}', horizontalalignment='left', verticalalignment='top')
        # ax.set_ylabel(r'\parbox{3cm}{$\mathrm{model}$ \hfill '+chisqstr+r'\null}', horizontalalignment='left', verticalalignment='top', fontsize=5)
        ax.set_yticks(args.pos)
        ax.yaxis.tick_right()
        ax.set_yticklabels(args.labels)
        axtwinx = ax.twiny()
        axtwinx.set_xlim(ax.get_xlim())
        axtwinx.tick_params(axis='x', which='both', direction='in', width=lpd.axeslinewidth, left=False, top=True, right=False, bottom=False, labelleft=False, labeltop=False,
                        labelright=False, labelbottom=False)
        if args.xticks != "auto":
            ax.set_xticks([float(i) for i in args.xticks])
            axtwinx.set_xticks([float(i) for i in args.xticks])
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        axtwinx.xaxis.set_minor_locator(AutoMinorLocator(2))
        ax.yaxis.label.set_size(10)
        ax.tick_params(axis='y', which='major', labelsize=10)

    if corr == "BB":
        add_group_labels(ax)

    lpd.create_folder(args.outputpath)
    outfile = args.outputpath + "/"+args.corr+"_kappa" + args.suffix + ".pdf"
    fig.savefig(outfile, backend=backend)
    print("saved ", outfile)


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
