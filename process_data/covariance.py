#!/usr/local/bin/python3.7m -u
import numpy
import lib_process_data as lpd
import _2_reduce_data as rd
from matplotlib.ticker import MaxNLocator
from matplotlib import pyplot as plt


def plot_correlation_matrix(data, idx, flow_radii, title, prefix):
    # bring data in correct shape for numpy.cov
    data = numpy.copy(data[:, :, idx])
    data = numpy.swapaxes(data, 0, 1)
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

    plot_correlation_matrix(polyakov_real, nt_half, flow_radii, r'$\mathrm{corr}[X(\tau_\mathrm{F}), X(\tau_\mathrm{F}\prime)], X= U(\beta, 0) $', "poly")

    for i in range(nt_half):
        label = '{0:.2f}'.format((i+1)/nt)
        plot_correlation_matrix(EE_numerator, i, flow_radii, r'$\mathrm{corr}[X(\tau_\mathrm{F}), X(\tau_\mathrm{F}\prime)], X= U(\beta, \tau) E(\tau) U(\tau, 0) E(0), \tau T ='+label+r'$', "EE_tauT"+label)


if __name__ == '__main__':
    main()
    lpd.save_script_call()