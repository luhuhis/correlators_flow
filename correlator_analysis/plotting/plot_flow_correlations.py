#!/usr/bin/env python3
import numpy
import lib_process_data as lpd
from correlator_analysis.double_extrapolation import _2_reduce_data as rd
from matplotlib.ticker import MaxNLocator
import matplotlib
import cmasher
from matplotlib.backends.backend_pdf import PdfPages

def compute_XX_corr(data):
    """
    function for bootstrap routine that computes an XX correlator (with XX=EE or BB) normalized by the polyakov loop. numerator data (i.e. --X--X--) is first index, polyakov loop is second index of data.
    """
    numerator_mean = numpy.mean(data[0], axis=0)
    print(data[1].shape)
    denominator_mean = numpy.mean(data[1], axis=0)
    print(denominator_mean.shape)
    XX = numerator_mean/denominator_mean
    return XX


def plot_correlation_matrix(index, data, flow_radii, nt, title):
    data = data[index]

    tauT = (index + 1) / nt
    label = '{0:.2f}'.format(tauT)
    title = title + label + r'$'

    data = numpy.corrcoef(data)

    fig, ax, _ = lpd.create_figure(xlims=[0, 0.35], ylims=[0, 0.35], ylabel=r'$\sqrt{8\tau_\mathrm{F}}/\tau$', xlabel=r'$\sqrt{8\tau_\mathrm{F}\prime}/\tau$', constrained_layout=True)
    ax.set_aspect('equal', 'box')

    # cax = ax.pcolormesh(flow_radii, flow_radii, data, cmap=cmasher.fusion, shading='auto', linewidth=0, vmin=-1, vmax=1, rasterized=True)  # 'seismic'

    ax.set_title(title)
    ax.xaxis.set_major_locator(MaxNLocator(4))
    ax.yaxis.set_major_locator(MaxNLocator(4))

    cax = ax.contourf(flow_radii/tauT, flow_radii/tauT, data, numpy.linspace(0, 1, 11), cmap=cmasher.get_sub_cmap(cmasher.fusion, 0.5, 1))
    for c in cax.collections:
        c.set_edgecolor("face")
    cbar = fig.colorbar(cax, ax=ax, shrink=0.7)
    cbar.set_ticks(numpy.arange(0, 1.01, 0.2))
    ax.axvline(x=1.5 / nt / tauT, color='grey', alpha=0.8, zorder=100, lw=1)
    ax.axhline(y=1.5 / nt / tauT, color='grey', alpha=0.8, zorder=100, lw=1)

    return fig


def main():

    # parse cmd line arguments
    parser, requiredNamed = lpd.get_parser()
    requiredNamed.add_argument('--conftype', help="format: s096t20_b0824900 for quenched or s096t20_b0824900_m002022_m01011 for hisq", required=True)
    parser.add_argument('--basepath', type=str, default="../../data/merged/")
    parser.add_argument('--outputfolder', default="./")
    parser.add_argument('--nproc', type=int, default=20)
    args = parser.parse_args()

    beta, ns, nt, nt_half = lpd.parse_conftype(args.conftype)

    flow_times, n_flow, n_datafiles, n_streams, n_files_per_stream, XX_data, confnums = rd.load_merged_data(args.qcdtype, args.corr, args.conftype, args.basepath, None, only_metadata=True)
    flow_radii = numpy.sqrt(8*flow_times)/nt

    # # combined matrix for EE correlator
    # XX_std_dev = numpy.sqrt((numerator_std_dev / denominator_mean) ** 2 + (
    #             numerator_mean * denominator_std_dev / (denominator_mean ** 2)) ** 2)

    # EE_numerator = numpy.asarray(XX_data[0])
    # polyakov_real = numpy.asarray(XX_data[1])

    # print(EE_numerator.shape)
    # print(polyakov_real.shape)

    # data = numpy.concatenate((EE_numerator, polyakov_real), axis=2)

    XX_samples = numpy.load(lpd.get_merged_data_path(args.qcdtype, args.corr, args.conftype, args.basepath) + args.corr + "_" + args.conftype + "_samples.npy")
    print(XX_samples.shape)

    # XX_samples, _, _ = bootstr.bootstr(est.compute_XX_corr, data, numb_samples=10000, sample_size=n_datafiles, conf_axis=0, return_sample=True,
    #                                          same_rand_for_obs=True, parallelize=True, nproc=args.nproc, seed=0, err_by_dist=True)

    # XX_samples = numpy.asarray(XX_samples)
    XX_samples = numpy.swapaxes(XX_samples, 0, 2)

    # bring data in correct shape for numpy.cov
    # polyakov_real = numpy.copy(polyakov_real[:, :, 0])
    # polyakov_real = numpy.swapaxes(polyakov_real, 0, 1)
    # print(polyakov_real.shape)
    # plot_correlation_matrix(polyakov_real, flow_radii, r'$\mathrm{corr}[X(\tau_\mathrm{F}), X(\tau_\mathrm{F}\prime)], X= U(\beta, 0) $', "poly")

    title = r'$\mathrm{corr}[G(\tau_\mathrm{F}), G(\tau_\mathrm{F}\prime)], \tau T ='
    figs = lpd.parallel_function_eval(plot_correlation_matrix, range(nt_half), args.nproc, XX_samples, flow_radii, nt, title)
    # data = numpy.copy(EE_numerator[:, :, i])
    # data = numpy.swapaxes(data, 0, 1)
    # plot_correlation_matrix(data, flow_radii/tauT, r'$\mathrm{corr}[X(\tau_\mathrm{F}), X(\tau_\mathrm{F}\prime)], \newline X= U(\beta, \tau) E(\tau) U(\tau, 0) E(0), \tau T ='+label+r'$', args.outputfolder+"EE_tauT"+label)

    # save figures to pdf
    print("save figures...")
    lpd.set_rc_params()  # for some reason we need to repeat this here...
    with PdfPages(args.outputfolder+"/"+args.conftype+"/EE_"+args.conftype+"_correlation.pdf") as pdf:
        for fig in figs:
            pdf.savefig(fig)
            matplotlib.pyplot.close(fig)


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
