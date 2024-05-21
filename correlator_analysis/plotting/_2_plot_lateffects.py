#!/usr/bin/env python3

import argparse
import numpy
import matplotlib
import lib_process_data as lpd
from matplotlib import cm
from matplotlib.backends.backend_pdf import PdfPages

import warnings
# # warnings.simplefilter("ignore", OptimizeWarning)
# warnings.simplefilter('error', UserWarning)
warnings.filterwarnings('ignore', r'(.*?)converting a masked element to nan.')
warnings.filterwarnings('ignore', r'(.*?)More than 20 figures have been opened.(.*?)')


def plot_single_flowtime(index, flowradii, args, tauT, XX, XX_err, labels, ylims, xlims):

    nplot = len(XX)
    colors = [*[cm.tab10(nplot-i-1) for i in range(nplot)]]
    # markers = ['s', 'o', 'D', 'H', 'h']
    zorders = numpy.linspace(-50, -1, nplot, dtype=int)

    if args.use_tex:
        displaystyle = r'\displaystyle'
    else:
        displaystyle = r''
    ylabel = r'$ ' + displaystyle + r'\frac{G_E }{G^\mathrm{norm}}$'

    fig, ax, _ = lpd.create_figure(xlims=xlims, ylims=ylims, xlabel=r'$\tau T$', ylabel=ylabel, UseTex=args.use_tex)
    plots = []
    ax.set_xticks((0.0, 0.1, 0.2, 0.3, 0.4, 0.5))

    ax.text(0.97, 0.99, r'$\sqrt{8\tau_\mathrm{F}} T=$ ' + '{0:.3f}'.format(flowradii[index]), ha='right', va='top', transform=ax.transAxes, bbox=lpd.labelboxstyle)

    # loop over all the different data sets and plot settings
    for tauTs, thisXX, thisXX_err, label, color, zorder in zip(tauT, XX, XX_err, labels, colors, zorders):
        thisplot = ax.errorbar(tauTs, thisXX[index], thisXX_err[index], fmt='|', color=color, label=label, zorder=zorder + 2)
        plots.append(thisplot)

    ax.legend(handles=plots)
    handles, labels = ax.get_legend_handles_labels()
    leg = ax.legend(handles[::-1], labels[::-1], title=r'$N_\tau$', loc="lower right", bbox_to_anchor=(1.01, 0.1), **lpd.leg_err_size(1, 0.3))
    leg.set_zorder(100)
    if index != 0:
        lower_limit = lpd.lower_tauT_limit_(flowradii[index])
        ax.axvline(x=lower_limit, **lpd.verticallinestyle)
        # ax.text(lower_limit*1.01, args.lower_limit_text_pos, r'$ 3\sqrt{8\tau_F}T$', ha='left', va='center')

    return fig


def apply_tree_level_imp(XX,nt, flowtimes, corr, flowaction, gaugeaction):
    for i in range(XX.shape[0]):
        flowtime = flowtimes[i]
        for j in range(XX.shape[1]):
            XX[i, j] = XX[i, j] / lpd.G_latt_LO_flow(j, flowtime, corr, nt, flowaction, gaugeaction) * nt**4
    return XX


def load_data(args, gaugeaction, flowaction):
    flowtimesT2 = numpy.loadtxt(args.flowtimesT2)

    # load lattice data and interpolations
    XX = []
    XX_err = []
    tauT = []
    Ntaus = []
    XX_int = []

    for conftype in args.conftypes:

        inputfolder = lpd.get_merged_data_path(args.qcdtype, args.corr, conftype, args.basepath)
        beta, ns, nt, nt_half = lpd.parse_conftype(conftype)
        Ntaus.append(nt)
        tmp = numpy.loadtxt(inputfolder + "/" + args.corr + "_" + conftype + ".dat")
        tmp = apply_tree_level_imp(tmp, nt, flowtimesT2*nt**2, args.corr, flowaction, gaugeaction)
        XX.append(tmp)
        tmp = numpy.loadtxt(inputfolder + "/" + args.corr + "_err_" + conftype + ".dat")
        tmp = numpy.fabs(apply_tree_level_imp(tmp, nt, flowtimesT2 * nt ** 2, args.corr, flowaction, gaugeaction))
        XX_err.append(tmp)

        # TODO load interpolations.
        # tmp = numpy.loadtxt(inputfolder + args.corr + "_" + args.conftype + "_interpolation_mean.npy")
        # XX_int.append(tmp)

        tauT.append(lpd.get_tauTs(nt))
    return tauT, flowtimesT2, XX, XX_err, Ntaus


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--qcdtype', required=True, type=str)
    parser.add_argument('--conftypes', nargs='*', type=str)
    parser.add_argument('--continuum', type=str, help="path to continuum mean")
    parser.add_argument('--continuum_err', type=str, help="path to continuum err")
    parser.add_argument('--flowtimesT2', type=str, help="path to file that contains flowtimes*T^2")
    parser.add_argument('--flow_index_range', default=(0, -1), type=int, nargs=2,
                        help="which flow indices to consider (default considers all). useful if youre only interested in some speficic ones.")
    parser.add_argument('--use_tex', default=True, action="store_false")
    parser.add_argument('--nproc', default=20, type=int, help="number of processes for parallelization over flow times.")
    parser.add_argument('--corr', choices=["EE", "BB"], required=True)
    parser.add_argument('--basepath', type=str)
    parser.add_argument('--ylims', type=float, nargs=2, default=None)
    parser.add_argument('--xlims', type=float, nargs=2, default=[0.15, 0.52])
    parser.add_argument('--outputfolder', type=str, default=None)
    parser.add_argument('--hide_cont', default=False, action="store_true")
    parser.add_argument('--lower_limit_text_pos', type=float)
    args = parser.parse_args()

    fermions, temp, flowtype, gaugeaction, flowaction = lpd.parse_qcdtype(args.qcdtype)

    tauT, flowtimesT2, XX, XX_err, Ntaus = load_data(args, gaugeaction, flowaction)

    labels = [r'$' + str(nt) + r'$' for nt in Ntaus]

    # add cont
    if not args.hide_cont:
        XX.append(numpy.loadtxt(args.continuum))
        XX_err.append(numpy.loadtxt(args.continuum_err))
        tauT.append(lpd.get_tauTs(Ntaus[-1]))
        labels.append(r'cont.')

    if args.flow_index_range == (0, -1):
        indices = range(0, len(flowtimesT2))
    else:
        indices = range(*args.flow_index_range)

    matplotlib.rcParams['figure.max_open_warning'] = 0  # suppress warning due to possibly large number of figures...
    figs = lpd.parallel_function_eval(plot_single_flowtime, indices, args.nproc, numpy.sqrt(8*flowtimesT2), args, tauT, XX, XX_err, labels, args.ylims, args.xlims)

    # save figures to pdf
    lpd.set_rc_params()  # for some reason we need to repeat this here...
    if args.outputfolder is None:
        args.outputfolder = lpd.get_plot_path(args.qcdtype, args.corr, "")
    lpd.create_folder(args.outputfolder)
    filename = args.outputfolder + "/" + args.corr + "_latt_effects.pdf"
    print(f"saving {filename}")
    with PdfPages(filename) as pdf:
        for fig in figs:
            pdf.savefig(fig)
            matplotlib.pyplot.close(fig)


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()