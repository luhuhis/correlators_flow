#!/usr/bin/env python3
import lib_process_data as lpd
import numpy

import _2_reduce_data as rd

def main():

    parser, requiredNamed = lpd.get_parser()
    requiredNamed.add_argument('--conftype', help='format example: s096t32_b0824900_m002022_m01011', type=str, required=True)
    parser.add_argument('--outputpath', help='where to store the plot', type=str)
    parser.add_argument('--show_id', help="list of indices that shall be plotted. e.g. 0 1 2 3 for the first 4 streams", type=int, nargs='*', default=None)
    parser.add_argument('--vlines', help="where to plot vlines", nargs='*', type=float)
    parser.add_argument('--suffix', help="suffix for output file", default="")
    parser.add_argument('--basepath', type=str, default="")
    parser.add_argument('--basepath_plot', type=str, default="")
    args = parser.parse_args()
    beta, ns, nt, nt_half = lpd.parse_conftype(args.conftype)

    outputfolder = lpd.get_plot_path(args.qcdtype, args.corr, args.conftype, args.basepath_plot) if not args.outputpath else args.outputpath
    lpd.create_folder(outputfolder)

    flow_times, n_flow, n_datafiles, n_streams, n_files_per_stream, XX_data, confnums = rd.load_merged_data(args.qcdtype, args.corr, args.conftype, args.basepath, None)
    polyakov_real = XX_data[1]
    numerator_real = XX_data[0]

    if args.show_id is None:
        args.show_id = range(n_streams)

    xlabel = r'n traj.'
    ylabel = r'$\mathrm{Re}[U_{(\beta,0)}]$'
    ylabelpos = (0.01, 0.99)

    flowidx = n_flow-1
    fig, ax, _ = lpd.create_figure(xlabel=xlabel, ylabel=ylabel, xlabelpos=(0.99, 0.01), ylabelpos=ylabelpos, figsize=(5,2.5))

    ax.set_xlabel(ax.get_xlabel(), fontsize=8, rotation=0, horizontalalignment='right', verticalalignment='bottom')
    ax.set_ylabel(ax.get_ylabel(), fontsize=8, rotation=0, horizontalalignment='left', verticalalignment='top')

    offset = 0
    for k in range(n_streams):
        x = confnums[k]
        y = polyakov_real[offset:offset+n_files_per_stream[k],flowidx,0]
        offset += n_files_per_stream[k]
        if k in args.show_id:
            ax.errorbar(x, y, fmt='.-', lw=0.3, fillstyle='full', markersize=2.5, mew=0, label=str(k), alpha=1)

    legendoptions = dict(edgecolor='none', fancybox=False, facecolor="w", columnspacing=0.1,
                   labelspacing=0.1, borderpad=0.1, handletextpad=0.4, handlelength=1, loc="center left", bbox_to_anchor=[1,0.5])
    ax.legend(**legendoptions)
    if args.vlines:
        for x in args.vlines:
            ax.axvline(x=x, **lpd.verticallinestyle)
    ax.set_title(r'$'+str(ns)+r'^3 \times '+str(nt)+r'$, $\beta='+str(beta)+'$,'+ r'$\sqrt{8\tau_F}/a='+r'{0:.3f}'.format(numpy.sqrt(8*flow_times[-1]))+r', \sqrt{8\tau_F}T \approx '+r'{0:.2f}'.format(numpy.sqrt(8*flow_times[-1])/nt)+r'$'+', ms/5, HISQ', x=0.5, fontsize=8)

    ax.set_xlim(xmin = 0)

    fig.savefig(outputfolder+"/"+"polyakovloop_MCtime"+args.suffix+".pdf")


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
