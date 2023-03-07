#!/usr/bin/env python3
import numpy

# https://stackoverflow.com/questions/44810437/problems-with-curve-fit-fitting-highly-correlated-data

import lib_process_data as lpd
import numpy as np
from scipy.optimize import curve_fit
from correlator_analysis.double_extrapolation import _2_reduce_data as rd
import argparse
np.set_printoptions(linewidth=np.inf, precision=3, threshold=np.inf, floatmode="fixed")


def f(x, a, b):
    return a + b*x


def pert_corr(tauf1, tauf2):
    r2 = (tauf1 / tauf2) ** 2
    return 4 * r2 / (1 + r2) ** 2



def get_data_corrmatrix(indices, tau_index):
    parser = argparse.ArgumentParser()

    parser.add_argument('--qcdtype', default="quenched_1.50Tc_zeuthenFlow")
    parser.add_argument('--corr', default="EE")
    parser.add_argument('--conftype', default="s144t36_b0754400")
    parser.add_argument('--basepath', type=str, default="../../../data/merged/")
    parser.add_argument('--outputfolder', default="./")
    args = parser.parse_args()

    # relflows = 8*flow_times/nt /0.5
    basepath = lpd.get_merged_data_path(args.qcdtype, args.corr, args.conftype, args.basepath)
    x = np.loadtxt(basepath + "EE_s144t36_b0754400_relflows.txt")**2
    XX_samples = np.load(basepath + args.corr + "_" + args.conftype + "_interpolation_relflow_samples.npy")
    XX_samples = np.swapaxes(XX_samples, 0, 2)[tau_index]
    # XX_samples = np.swapaxes(XX_samples, 0, 1)

    print(XX_samples.shape)
    y = np.nanmedian(XX_samples, axis=0)
    e = lpd.dev_by_dist(XX_samples, axis=0)

    tmpcov = np.cov(XX_samples, rowvar=False)
    tmpcorr = np.corrcoef(XX_samples, rowvar=False)

    # filter out data points outside of extr window
    x = x[indices]
    y = y[indices]
    e = e[indices]

    ndata = len(x)
    data_cov = np.empty((ndata, ndata))
    data_corr = np.empty((ndata, ndata))
    for i in range(ndata):
        for j in range(ndata):
            data_cov[i, j] = tmpcov[indices[i], indices[j]]
            data_corr[i, j] = tmpcorr[indices[i], indices[j]]

    print("data correlation matrix:\n", data_corr)

    # get bootstrap estimate
    nsamples = XX_samples.shape[0]
    fitparams_arr = []
    for m in range(nsamples):
        fitparams, _ = curve_fit(f, x, XX_samples[m][indices], sigma=e)
        fitparams_arr.append(fitparams)
    fitparams = np.nanmedian(fitparams_arr, axis=0)
    fitparams_dev = lpd.dev_by_dist(fitparams_arr, axis=0)
    bootstrap = numpy.column_stack((fitparams, fitparams_dev))

    return x, y, e, data_cov, data_corr, bootstrap


def filter_covariance_matrix(cov_mat, num_largest_eigenvalues, eigenvalue_cutoff):
    # Compute the eigendecomposition of the covariance matrix
    eigenvalues, eigenvectors = np.linalg.eig(cov_mat)

    print(eigenvalues)

    # Find the indices of the largest eigenvalues
    largest_eigenvalue_indices = np.argsort(eigenvalues)[-num_largest_eigenvalues:]

    # Set all eigenvalues to zero except the largest ones
    eigenvalues_filtered = np.zeros_like(eigenvalues)
    eigenvalues_filtered[largest_eigenvalue_indices] = eigenvalues[largest_eigenvalue_indices]

    # Find the smallest filtered eigenvalue and use it as the regularization constant
    regularization_constant = np.max([eigenvalue_cutoff, np.min(eigenvalues_filtered[eigenvalues_filtered > 0])])

    # Add the regularization constant to any filtered out eigenvalues smaller than it
    eigenvalues_filtered[eigenvalues_filtered < regularization_constant] = regularization_constant
    print(eigenvalues_filtered)

    # Reconstruct the covariance matrix using the filtered eigenvectors and eigenvalues
    cov_mat_filtered = eigenvectors @ np.diag(eigenvalues_filtered) @ eigenvectors.T

    return cov_mat_filtered


def cov_to_corr(cov):
    # Calculate the standard deviations of each variable
    stds = np.sqrt(np.diag(cov))

    # Create a diagonal matrix of the standard deviations
    stds_matrix = np.diag(1 / stds)

    # Calculate the correlation matrix
    corr = stds_matrix @ cov @ stds_matrix

    return corr


def main():

    tau_index = -1

    # indices = range(11, 21, 2)
    indices = [11, 16, 21]
    # indices = range(11, 21, 1)

    x, y, e, data_cov, data_corr, bootstrap = get_data_corrmatrix(indices, tau_index)

    data_cov = filter_covariance_matrix(data_cov, 2, 0)

    print("filtered correlation matrix\n", cov_to_corr(data_cov))

    # creat covariance matrix using data and correlation matrix from pert theory.
    # ndata = len(x)
    # cov = np.empty((ndata, ndata))
    # for i in range(ndata):
    #     for j in range(ndata):
    #         if i == j:
    #             cov[i, i] = e[i]**2
    #         else:
    #             cov[i, j] = pert_corr(x[i], x[j]) * e[i]*e[j]

    # print correlation matrix
    # corr = np.empty(cov.shape)
    # for i in range(cov.shape[0]):
    #     for j in range(cov.shape[1]):
    #         corr[i,j] = cov[i,j] / np.sqrt(cov[i,i]*cov[j,j])
    # print("pert. correlation matrix\n", corr)

    params_naive, pcov_naive = curve_fit(f, x, y, sigma=e)
    fit_err_naive = np.sqrt(np.diag(pcov_naive))
    params, pcov = curve_fit(f, x, y, sigma=data_cov, p0=params_naive) #[i, j]
    fit_err = np.sqrt(np.diag(pcov))

    print("fit errors with cov matrix", params, "+-", fit_err)
    print("fit errors w/o  cov matrix", params_naive, "+-", fit_err_naive)
    print("fit error bootstrap       ", bootstrap[:,0], "+-", bootstrap[:,1])

    color_naive = 'C1'
    color = 'C2'
    color_boot = 'C3'

    fig, ax, _ = lpd.create_figure(xlims=[-0.005, 0.1], ylims=[3.5, 3.9],
                                   xlabel=r'$8\tau_\mathrm{F}/\tau^2$', ylabel=r'$\displaystyle\frac{G}{G^\mathrm{norm}}$')

    ax.errorbar(x, y, e, label=r'data ($N_\tau=36, \tau T =0.5$)', fmt='|', color='C0')
    ax.errorbar([0, *x], f(np.array([0, *x]), *params_naive), fmt=':', linewidth=1, label="uncorrelated fit", color=color_naive, zorder=-1)
    ax.errorbar([0, *x], f(np.array([0, *x]), *params), linewidth=1, fmt='--', label="correlated fit\nwith eigenvalue cutoff", color=color, zorder=-2)
    ax.errorbar([0, *x], f(np.array([0, *x]), bootstrap[0, 0], bootstrap[1, 0]), linewidth=1, fmt='-', label="bootstrap", color=color_boot, zorder=-3)
    ax.errorbar(-0.00025, bootstrap[0, 0], bootstrap[0, 1], color=color_boot, fmt='|', zorder=0)
    ax.errorbar(0.00025, params[0], fit_err[0], color=color, zorder=1)
    ax.errorbar(0, params_naive[0], fit_err_naive[0], color=color_naive, zorder=2)

    # ax.grid(zorder=-100, lw=0.5, color='grey')

    # ax.text(0.01, 0.01, r'$\tau T =0.5$', ha='left', va='bottom', transform=ax.transAxes)
    ax.legend(framealpha=1, handlelength=1.5, fontsize=8)

    fig.savefig("correlated_fit.pdf")

    return










if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()