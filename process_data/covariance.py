#!python3
import numpy
import lib_process_data as lpd
import _2_reduce_data as rd
from matplotlib.ticker import MaxNLocator
from matplotlib import pyplot as plt
from latqcdtools.statistics import bootstr


def compute_XX_corr(data):
    """
    function for bootstrap routine that computes an XX correlator (with XX=EE or BB) normalized by the polyakov loop. numerator data (i.e. --X--X--) is first index, polyakov loop is second index of data.
    """
    numerator_mean = numpy.mean(data[0], axis=0)
    denominator_mean = numpy.mean(data[1], axis=0)
    XX = numerator_mean/denominator_mean
    return XX


def plot_correlation_matrix(data, flow_radii, title, prefix):
    data = numpy.corrcoef(data)

    fig, ax, _ = lpd.create_figure(xlims=[0, 0.3], ylims=[0, 0.3], xlabel=r'$\sqrt{8\tau_F}T$', ylabel=None, xlabelpos=(0.6, -0.2), ylabelpos=(-0.2, 1.15),
                                   tickpad=1,
                                   figsize=(3, 2.5), UseTex=False, fig=None, subplot=111, no_ax=False)
    plt.rc('axes', unicode_minus=False)

    cax = ax.pcolormesh(flow_radii, flow_radii, data, cmap='seismic', shading='auto', linewidth=0, vmin=-1, vmax=1, rasterized=True)
    fig.colorbar(cax)

    ax.set_title(title, fontsize=6)

    ax.set_aspect('equal', 'box')
    ax.xaxis.set_major_locator(MaxNLocator(5))
    ax.yaxis.set_major_locator(MaxNLocator(5))

    fig.savefig(prefix + "_correlation.pdf")
    plt.close(fig)


def main():

    # parse cmd line arguments
    parser, requiredNamed = lpd.get_parser()
    requiredNamed.add_argument('--conftype', help="format: s096t20_b0824900 for quenched or s096t20_b0824900_m002022_m01011 for hisq", required=True)
    args = parser.parse_args()

    beta, ns, nt, nt_half = lpd.parse_conftype(args.conftype)

    flow_times, n_flow, n_datafiles, n_streams, XX_data = rd.load_merged_data(args.qcdtype, args.corr, args.conftype)

    flow_radii = numpy.sqrt(8*flow_times)/nt

    # # combined matrix for EE correlator
    # XX_std_dev = numpy.sqrt((numerator_std_dev / denominator_mean) ** 2 + (
    #             numerator_mean * denominator_std_dev / (denominator_mean ** 2)) ** 2)

    EE_numerator = numpy.asarray(XX_data[0])
    polyakov_real = numpy.asarray(XX_data[1])

    XX_samples, XX, XX_err = bootstr.bootstr(compute_XX_corr, XX_data, numb_samples=1000, sample_size=n_datafiles, conf_axis=1, return_sample=True,
                                             same_rand_for_obs=False, parallelize=True, nproc=20, seed=0, err_by_dist=True)

    print(numpy.asarray(XX))
    print(numpy.asarray(XX_err))

    XX_samples = numpy.asarray(XX_samples)
    XX_samples = numpy.swapaxes(XX_samples, 0, 2)

    # bring data in correct shape for numpy.cov
    polyakov_real = numpy.copy(polyakov_real[:, :, 0])
    polyakov_real = numpy.swapaxes(polyakov_real, 0, 1)
    print(polyakov_real.shape)
    plot_correlation_matrix(polyakov_real, flow_radii, r'$\mathrm{corr}[X(\tau_\mathrm{F}), X(\tau_\mathrm{F}\prime)], X= U(\beta, 0) $', "poly")

    for i in range(nt_half):
        label = '{0:.2f}'.format((i+1)/nt)
        data = numpy.copy(EE_numerator[:, :, i])
        data = numpy.swapaxes(data, 0, 1)
        plot_correlation_matrix(data, flow_radii, r'$\mathrm{corr}[X(\tau_\mathrm{F}), X(\tau_\mathrm{F}\prime)], X= U(\beta, \tau) E(\tau) U(\tau, 0) E(0), \tau T ='+label+r'$', "EE_tauT"+label)
        plot_correlation_matrix(XX_samples[i], flow_radii, r'$\mathrm{corr}[G(\tau_\mathrm{F}), G(\tau_\mathrm{F}\prime)], \tau T ='+label+r'$', "EE_corr_tauT"+label)


if __name__ == '__main__':
    main()
    lpd.save_script_call()
