#!/usr/local/bin/python3.7m -u

import lib_process_data as lpd
import numpy as np
import math
import scipy.integrate
import scipy.optimize
import scipy.interpolate
from typing import NamedTuple
import latqcdtools.bootstr

global args


def Gnorm(tauT):
    norm = np.pi ** 2 * (np.cos(np.pi * tauT) ** 2 / np.sin(np.pi * tauT) ** 4 + 1 / (3 * np.sin(np.pi * tauT) ** 2))
    return norm

# ==============================================================


def Kernel(OmegaByT, tauT):
    return np.cosh(OmegaByT / 2 - OmegaByT * tauT) / np.sinh(OmegaByT / 2)


def En(n, OmegaByT, mu):
    x = np.log(1 + OmegaByT / np.pi)
    y = x / (1 + x)
    if mu == "alpha":
        return np.sin(np.pi * n * y)
    elif mu == "beta":
        return np.sin(np.pi * y) * np.sin(np.pi * n * y)


class SpfArgs(NamedTuple):
    model: int
    mu: str
    constrain: bool
    PhiuvByT3: scipy.interpolate.InterpolatedUnivariateSpline
    n_max: int


# ========= Spf models =========

def SpfByT3(OmegaByT, MaxOmegaByT, spfargs, *fit_params_0):
    if args.model == 3:
        return np.maximum(0.5*fit_params_0[0][0]*OmegaByT, spfargs.PhiuvByT3(OmegaByT)*fit_params_0[0][1])
    if spfargs.constrain:
        coef_tmp = 1
        for i in range(1, spfargs.n_max + 1):
            coef_tmp += fit_params_0[0][i] * En(i, MaxOmegaByT, spfargs.mu)
        c_nmax = (spfargs.PhiuvByT3(MaxOmegaByT) / np.sqrt((0.5 * fit_params_0[0][0] * MaxOmegaByT) ** 2 + (spfargs.PhiuvByT3(MaxOmegaByT)) ** 2) - coef_tmp) / En(
                spfargs.n_max + 1, MaxOmegaByT, spfargs.mu)
    coef = 1
    for i in range(1, spfargs.n_max+1):
        coef += fit_params_0[0][i]*En(i, OmegaByT, spfargs.mu)

    if spfargs.constrain:
        coef += c_nmax * En(spfargs.n_max + 1, OmegaByT, spfargs.mu)

    if coef < 0:
        if args.verbose:
            print("negative SPF")
        return np.inf
        # return np.inf  # this results in infinite chisq whenever the spf becomes zero.
    return np.sqrt((0.5*fit_params_0[0][0]*OmegaByT)**2+(spfargs.PhiuvByT3(OmegaByT))**2)*coef


# ==============================================================


def Integrand(OmegaByT, tauT, MaxOmegaByT, spfargs, *fit_params_0):
    return 1. / np.pi * Kernel(OmegaByT, tauT) * SpfByT3(OmegaByT, MaxOmegaByT, spfargs, *fit_params_0)


def TargetCorr(tauT, MaxOmegaByT, spfargs, *fit_params_0):
    CorrTrial = []
    for i in range(len(tauT)):
        try:
            CorrTrial.append(scipy.integrate.quad(lambda OmegaByT: Integrand(OmegaByT, tauT[i], MaxOmegaByT, spfargs, *fit_params_0), 0, MaxOmegaByT)[0])
        except OverflowError as e:
            print(str(e) + " for integration appears at tauT=" + str(tauT[i]))
            return np.inf
    return CorrTrial


def chisq_dof(fit_params_0, xdata, ydata_sample, edata, MaxOmegaByT, spfargs, verbose=False):
    res = (ydata_sample - TargetCorr(xdata, MaxOmegaByT, spfargs, fit_params_0)) / edata
    chisqdof = np.sum(res ** 2) / (len(xdata) - len(fit_params_0))
    if verbose:
        print(['{0:.7f} '.format(i) for i in fit_params_0], '{0:.4f}'.format(chisqdof))
    return chisqdof


def get_results(ydata_sample, fit_params_0, xdata, edata, MaxOmegaByT, spfargs, PhiUV, MaxIter, NtauT, corr, model, verbose=False):
    if verbose:
        print("current correlator sample:", ydata_sample)

    # fun: The objective function to be minimized.
    # x0: Initial guess.
    # args: Extra arguments passed to the objective function
    fit_res = scipy.optimize.minimize(fun=chisq_dof, x0=fit_params_0, args=(xdata, ydata_sample, edata, MaxOmegaByT, spfargs, verbose), method='Nelder-Mead',
                                      options={'maxiter': MaxIter, 'disp': True, 'xatol': args.tol, 'fatol': args.tol}, callback=None)  # bounds=[[0, np.inf], *([[-1, 1]]*spfargs.n_max)])

    # now use the fit results for the parameters to compute the fitted spectral function and correlator

    # spectral function
    Spf = []
    for i in range(len(PhiUV)):
        Spf.append(SpfByT3(PhiUV[i][0], MaxOmegaByT, spfargs, fit_res.x))
    return_nan = False
    if any(n < 0 for n in Spf):
        print("negative spf for this sample. returning nan.")
        return_nan = True

    # because in the fit we take the square of PhiIR, kappa can be negative. Since it is a square, the sign doesn't matter, so just take the absolute:
    if model == 2:
        if fit_res.x[0] < 0 and return_nan is False:
            print("kappa negative but spf positive. using absolute kappa value.")
            fit_res.x[0] = math.fabs(fit_res.x[0])

    # correlator
    fit_corr = TargetCorr(xdata, MaxOmegaByT, spfargs, fit_res.x)
    for i in range(NtauT):
        fit_corr[i] /= Gnorm(corr[i][0])

    # chisq
    chisqdof = chisq_dof(fit_res.x, xdata, ydata_sample, edata, MaxOmegaByT, spfargs, verbose)  # Note: this chisq_dof is appended to the end of Record!

    # stack the result into one long array, because the bootstrap is limited to 1D arrays. we'll need to accordingly extract this again later.
    result = np.hstack((fit_res.x, Spf, fit_corr, chisqdof))

    if return_nan:
        nans = np.empty(result.shape)
        nans[:] = np.nan
        return nans
    else:
        return result


def main():
    parser, requiredNamed = lpd.get_parser()

    requiredNamed.add_argument('--model', help='which model to use', choices=[2, 3], type=int, required=True)
    requiredNamed.add_argument('--PathPhiUV', help='the full path of the input phiuv in omega/T phi/T^3', type=str, required=True)
    requiredNamed.add_argument('--PhiUVtype', help='specify it this is LO (a) or NLO (b)', type=str, choices=["a", "b"], required=True)
    parser.add_argument('--mu', help='which "en" function to use', choices=["alpha", "beta"], type=str)
    parser.add_argument('--nmax', help='what nmax to use. valid only for model 1,2.', type=int, choices=[1, 2, 3, 4, 5, 6, 7])
    parser.add_argument('--constrain', help='force the spf to reach the UV limit at large omega', action="store_true")

    parser.add_argument('--nsamples', help='number of bootstrap samples to draw from the gaussian distribution', default=200, type=int)

    parser.add_argument('--PathOutputFolder', help='the path of output folder like /a/b', type=str)
    parser.add_argument('--PathInputFolder', help='the path of input folder like /a/b', type=str)

    parser.add_argument('--maxiter', help='maximum number of nelder-mead iterations', default=1000, type=int)
    parser.add_argument('--tol', help='subsequent iterations have to differ less than this for the function and param values', type=float, default=0.00001)
    parser.add_argument('--nproc', help='number of processes for the parallel bootstrap', default=1, type=int)
    parser.add_argument('--verbose', help='output current fit parameters at each iteration', action="store_true")
    parser.add_argument('--seed', help='seed for gaussian bootstrap sample drawings', default=None, type=int)
    parser.add_argument('--start_from_mean_fit', help='fit the mean first and use its fit params as the initial guess for all other fits', action="store_true")
    parser.add_argument('--asym_err', help='get asymmetric errors from the bootstrap to better reflect the distribution', action="store_true")
    parser.add_argument('--min_tauT', help='ignore corr data below this tauT', type=float, default=0)

    parser.add_argument('--add_suffix', help='add an extra suffix to the output files in order to not overwrite previous ones with similar parameters on a '
                                             'different data set', type=str, default="")

    global args
    args = parser.parse_args()

    if args.model == 2 and (not args.mu or not args.nmax):
        print("ERROR: Need mu and nmax for model 2.")
        return 1

    constrainstr = "s1" if not args.constrain else "s2"  # s1 = dont constrain, s2 = constrain
    startstr = "_d" if not args.start_from_mean_fit else "_m"
    if args.model == 2:
        modelidentifier = str(args.model)+"_"+args.PhiUVtype+"_"+constrainstr+"_"+str(args.mu)+"_"+str(args.nmax)
    elif args.model == 3:
        modelidentifier = str(args.model)+"_"+args.PhiUVtype
    fileidentifier = modelidentifier+"_"+str(args.nsamples)+"_"+'{:.0e}'.format(args.tol)+startstr+"_"+str(args.min_tauT) + args.add_suffix

    # read in the normalized correlator:
    inputfolder = "../"+lpd.get_merged_data_path(args.qcdtype, args.corr, "") if not args.PathInputFolder else args.PathInputFolder
    corr = np.loadtxt(inputfolder+args.corr+"_flow_extr.txt")

    # remove lines with NaN's:
    corr = corr[~np.isnan(corr).any(axis=1)]
    # filter out below min_tauT:
    corr = corr[corr[:, 0] >= args.min_tauT, :]

    # beta, ns, nt, nt_half = lpd.parse_conftype(args.conftype)
    NtauT = len(corr)

    # read in the phiuv. columns in order: OmegaByT, PhiUVByT3, err
    PhiUV = np.loadtxt(args.PathPhiUV)
    PhiUV = PhiUV[:, 0:2]

    # interpolate the UV spf for the integration
    # spline order: 1 linear, 2 quadratic, 3 cubic ...
    order = 3
    PhiuvByT3 = scipy.interpolate.InterpolatedUnivariateSpline(PhiUV[:, 0], PhiUV[:, 1], k=order)
    MaxOmegaByT = PhiUV[-1][0]

    # get rid of the pert. normalization in the correlator data
    CorrByT4 = np.copy(corr)
    for i in range(NtauT):
        CorrByT4[i][2] = corr[i][2] * Gnorm(corr[i][0])
        CorrByT4[i][1] = corr[i][1] * Gnorm(corr[i][0])

    xdata = CorrByT4[:, 0]
    ydata = CorrByT4[:, 1]
    edata = CorrByT4[:, 2]

    # set up initial guess for the fitted parameters. the initial guess for kappa is 1, and for the c_n is 0. for model 3 we only have one other fit parameter
    # (the overall coefficient for the UV part), whose initial guess is also 1.
    if args.model == 3:
        fit_params_0 = [1, 1]
    else:
        fit_params_0 = [1.]
        if args.constrain:
            args.nmax -= 1
        for i in range(args.nmax):
            fit_params_0.append(0.)

    nparam = len(fit_params_0)

    # constant parameters only used by the function SpfByT3
    spfargs = SpfArgs(args.model, args.mu, args.constrain, PhiuvByT3, nparam-1)

    if args.start_from_mean_fit:
        print("finding initial guess for fit params...")
        mean_result = get_results(ydata, fit_params_0, xdata, edata, MaxOmegaByT, spfargs, PhiUV, args.maxiter, NtauT, corr, args.model, verbose=True)

    nomega = len(PhiUV[:, 0])
    structure = np.asarray((
        (0, nparam),
        (nparam, nparam+nomega),
        (nparam+nomega, nparam+nomega+NtauT),
        (nparam+nomega+NtauT, nparam+nomega+NtauT+1)
    ), dtype=int)

    if args.start_from_mean_fit and not np.isnan(mean_result).any():
        fit_params_0 = np.asarray(mean_result[structure[0, 0]:structure[0, 1]])
    print("Initial guess for fit params:", fit_params_0)

    samples, results, error = \
        latqcdtools.bootstr.bootstr_from_gauss(get_results, ydata, edata, args.nsamples, sample_size=1, return_sample=True, seed=args.seed, err_by_dist=True,
                                               useCovariance=False, parallelize=True, nproc=args.nproc, asym_err=args.asym_err,
                                               args=(fit_params_0, xdata, edata, MaxOmegaByT, spfargs, PhiUV, args.maxiter, NtauT, corr, args.model, args.verbose))

    # make handling of the left and right 68-quantiles easier
    error = np.asarray(error)
    error = np.swapaxes(error, 0, 1)

    if args.PathOutputFolder:
        outputfolder = args.PathOutputFolder+"/"+fileidentifier+"/"
    else:
        outputfolder = inputfolder+"/spf/"+fileidentifier+"/"
    lpd.create_folder(outputfolder)

    np.savetxt(outputfolder+"samples_structure_"+fileidentifier+".dat", structure, header='This file contains pairs (a,b) of indeces with which to split the array '
                                                                                          'in the samples.npy in the correct way. Example: Spf=samples(a:b). \n rows for (a,b) in order: fit_resx, Spf, '
                                                                                          'fit_corr, chisqdof', fmt='%i')
    np.save(outputfolder+"samples_"+fileidentifier, samples)

    # extract the various quantities that have been stacked together due to bootstrap being limited to 1D arrays
    fit_resx = results[structure[0, 0]:structure[0, 1]]
    fit_resx_err = error[structure[0, 0]:structure[0, 1]]
    Spf = results[structure[1, 0]:structure[1, 1]]
    Spf_err = error[structure[1, 0]:structure[1, 1]]
    fit_corr = results[structure[2, 0]:structure[2, 1]]
    fit_corr_err = error[structure[2, 0]:structure[2, 1]]
    chisqdof = results[structure[3, 0]]
    chisqdof_err = error[structure[3, 0]]

    # combine fit params and chisqdof into one object for file storage
    fit_resx_data = np.column_stack((fit_resx, fit_resx_err))
    chisqdof_data = np.asarray(((chisqdof, *chisqdof_err),))
    fitparams_chisqdof = np.concatenate((fit_resx_data, chisqdof_data), axis=0)

    print("\nThe first line contains kappa/T^3. The last line contains chisq/dof. In between are c_i.\n"
          "param                  error_left    error_right:")
    print(fitparams_chisqdof)

    # make it say 1e-4 instead of 0.0001

    # save reconstructed fit spf
    np.savetxt(outputfolder + "/spffit_" + fileidentifier + ".dat", np.column_stack((PhiUV[:, 0], Spf, Spf_err)), fmt='%16.15e',
               header='omegaByT            SpfByT3               err(-/+)')

    # save fitparams and chisqdof
    np.savetxt(outputfolder + "/params_" + fileidentifier + ".dat", fitparams_chisqdof, fmt='%22.15e',
               header="The first line contains kappa/T^3. The last line contains chisq/dof. In between are c_i.\n"
                      "param                  error_left                 error_right:")

    # save reconstructed correlator
    np.savetxt(outputfolder + "/corrfit_" + fileidentifier + ".dat",
               np.column_stack((xdata, corr[:, 1], corr[:, 2], fit_corr, fit_corr_err)),
               fmt='%16.15e', header='tauT                corr(orig)            err(orig)             corr(fit)             err(-/+)')


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
