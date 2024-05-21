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


def plot_correlation_matrix(index, data, flow_radii, nt):
    data = data[index]

    tauT = (index + 1) / nt
    label = '{0:.2f}'.format(tauT)

    title = r'\begin{center}$\mathrm{corr}[G_E(\tau_\mathrm{F}), G_E(\tau_\mathrm{F}\prime)]$ \\ $\tau T ='
    title = title + label + r'$ \end{center}'

    data = numpy.corrcoef(data)

    xydata = flow_radii/tauT
    
    overwrite_with_perturbative_formula = False
    if overwrite_with_perturbative_formula:
        for i in range(len(xydata)):
            for j in range(len(xydata)):
                r2 = (xydata[i]/xydata[j])**2
                data[i,j] = (4 * r2 / (1+r2)**2)

    fig, ax, _ = lpd.create_figure(xlims=[0, 0.35], ylims=[0, 0.35], constrained_layout=True,
                                   ytwinticks=False, xlabelbox=None, ylabelbox=None)
    ax.set_aspect('equal', 'box')

    ax.xaxis.tick_top()
    ax.invert_yaxis()
    ax.set_xlabel(r'$\sqrt{8\tau_\mathrm{F}\prime}/\tau$', ha='left', va='center', rotation=90, zorder=1000)
    ax.xaxis.set_label_coords(0.02, 0.5)
    ax.set_ylabel(r'$\sqrt{8\tau_\mathrm{F}}/\tau$', ha='center', va='top', rotation=0, zorder=1000)
    ax.yaxis.set_label_coords(0.5, 0.98)

    # cax = ax.pcolormesh(flow_radii, flow_radii, data, cmap=cmasher.fusion, shading='auto', linewidth=0, vmin=-1, vmax=1, rasterized=True)  # 'seismic'

    ax.text(0.5, -0.02, title, ha='center', va='top', transform=ax.transAxes)

    ax.xaxis.set_major_locator(MaxNLocator(4))
    ax.yaxis.set_major_locator(MaxNLocator(4))

    cax = ax.contourf(xydata, xydata, data, numpy.linspace(0, 1, 21), cmap=cmasher.get_sub_cmap(cmasher.fusion, 0.5, 1))
    for c in cax.collections:
        c.set_edgecolor("face")
    cbar = fig.colorbar(cax, ax=ax, shrink=0.74)
    cbar.set_ticks(numpy.arange(0, 1.01, 0.2))
    ax.vlines(x=0.25, ymin=0.2495, ymax=0.3005, color='C1', alpha=1, zorder=100, lw=0.5, linestyles='dashed')
    ax.vlines(x=0.3, ymin=0.2495, ymax=0.3005, color='C1', alpha=1, zorder=100, lw=0.5, linestyles='dashed')
    ax.hlines(y=0.25, xmin=0.2495, xmax=0.3005, color='C1', alpha=1, zorder=100, lw=0.5, linestyles='dashed')
    ax.hlines(y=0.3, xmin=0.2495, xmax=0.3005, color='C1', alpha=1, zorder=100, lw=0.5, linestyles='dashed')

    onespacingflow = 1 / nt / tauT
    ax.hlines(y=onespacingflow, xmin=0, xmax=onespacingflow, color='C3', alpha=1, zorder=100, lw=0.5, linestyles='dashdot')
    ax.vlines(x=onespacingflow, ymin=0, ymax=onespacingflow, color='C3', alpha=1, zorder=100, lw=0.5, linestyles='dashdot')

    return fig


def main():

    # parse cmd line arguments
    parser, requiredNamed = lpd.get_parser()
    requiredNamed.add_argument('--conftype', help="format: s096t20_b0824900 for quenched or s096t20_b0824900_m002022_m01011 for hisq", required=True)
    parser.add_argument('--basepath', type=str)
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

    figs = lpd.parallel_function_eval(plot_correlation_matrix, range(nt_half), args.nproc, XX_samples, flow_radii, nt)
    # figs = lpd.serial_function_eval(plot_correlation_matrix, range(nt_half), XX_samples, flow_radii, nt)
    # data = numpy.copy(EE_numerator[:, :, i])
    # data = numpy.swapaxes(data, 0, 1)
    # plot_correlation_matrix(data, flow_radii/tauT, r'$\mathrm{corr}[X(\tau_\mathrm{F}), X(\tau_\mathrm{F}\prime)], \newline X= U(\beta, \tau) E(\tau) U(\tau, 0) E(0), \tau T ='+label+r'$', args.outputfolder+"EE_tauT"+label)

    # save figures to pdf
    filename = args.outputfolder+"/"+args.conftype+"/EE_"+args.conftype+"_correlation.pdf"
    print(f"saving {filename}")
    lpd.set_rc_params()  # for some reason we need to repeat this here...
    with PdfPages(filename) as pdf:
        for fig in figs:
            pdf.savefig(fig)
            matplotlib.pyplot.close(fig)


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
