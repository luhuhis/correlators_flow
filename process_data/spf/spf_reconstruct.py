import lib_process_data as lpd
import numpy as np
import math
import sys
import scipy.integrate
import scipy.optimize
import scipy.interpolate
from typing import NamedTuple


def Gnorm(tauT):
    norm = np.pi ** 2 * (np.cos(np.pi * tauT) ** 2 / np.sin(np.pi * tauT) ** 4 + 1 / (3 * np.sin(np.pi * tauT) ** 2))
    return norm


def bootstr_from_gauss(data, data_std_dev, Seed):
    data = np.asarray(data)
    data_std_dev = np.asarray(data_std_dev)
    np.random.seed(Seed)
    numb_observe = len(data)
    sampleval = []
    for k in range(numb_observe):
        sampleval.append(np.random.normal(data[k], data_std_dev[k]))
    sampleval = np.array(sampleval)
    return sampleval

# ==============================================================


def Kernel(OmegaByT, tauT):
    kernel = np.cosh(OmegaByT / 2 - OmegaByT * tauT) / np.sinh(OmegaByT / 2)
    return kernel


def En(n, OmegaByT, mu):
    x = np.log(1 + OmegaByT / np.pi)
    y = x / (1 + x)
    if mu == "alpha":
        return np.sin(np.pi * n * y)
    elif mu == "beta":
        return np.sin(np.pi * y) * np.sin(np.pi * n * y)
    else:
        print("wrong mu!!!")
        sys.exit(2)


class SpfArgs(NamedTuple):
    model: int
    mu: str
    constrain: bool
    PhiuvByT3: scipy.interpolate.InterpolatedUnivariateSpline


def SpfByT3(OmegaByT, MaxOmegaByT, spfargs, *fit_params_0):
    if spfargs.model == 3:
        return np.maximum(0.5 * fit_params_0[0][0] * OmegaByT, spfargs.PhiuvByT3(OmegaByT) * fit_params_0[0][1])
    else:
        n_max = len(*fit_params_0) - 1
        if spfargs.constrain:
            coef_tmp = 1
            for i in range(1, n_max + 1):
                coef_tmp += fit_params_0[0][i] * En(i, MaxOmegaByT, spfargs.mu)
            if spfargs.model == 2:
                c_nmax = (spfargs.PhiuvByT3(MaxOmegaByT) / np.sqrt((0.5 * fit_params_0[0][0] * MaxOmegaByT) ** 2 + (spfargs.PhiuvByT3(MaxOmegaByT)) ** 2) - coef_tmp) / En(
                    n_max + 1, MaxOmegaByT, spfargs.mu)
            else:
                c_nmax = (spfargs.PhiuvByT3(MaxOmegaByT) / (0.5 * fit_params_0[0][0] * MaxOmegaByT + spfargs.PhiuvByT3(MaxOmegaByT)) - coef_tmp) / En(n_max + 1, MaxOmegaByT, spfargs.mu)

        coef = 1
        for i in range(1, n_max+1):
            coef += fit_params_0[0][i] * En(i, OmegaByT, spfargs.mu)

        if spfargs.constrain:
            coef += c_nmax * En(n_max + 1, OmegaByT, spfargs.mu)

        if coef < 0:
            print("negative spf")
            return 1e20
        if spfargs.model == 2:
            return np.sqrt((0.5 * fit_params_0[0][0] * OmegaByT) ** 2 + (spfargs.PhiuvByT3(OmegaByT)) ** 2) * coef
        else:
            return (0.5 * fit_params_0[0][0] * OmegaByT + spfargs.PhiuvByT3(OmegaByT)) * coef


def Integrand(OmegaByT, tauT, MaxOmegaByT, spfargs, *fit_params_0):
    return 1. / np.pi * Kernel(OmegaByT, tauT) * SpfByT3(OmegaByT, MaxOmegaByT, spfargs, *fit_params_0)


def TargetCorr(tauT, MaxOmegaByT, spfargs, *fit_params_0):
    CorrTrial = []
    for i in range(len(tauT)):
        try:
            CorrTrial.append(scipy.integrate.quad(lambda OmegaByT: Integrand(OmegaByT, tauT[i], MaxOmegaByT, spfargs, *fit_params_0), 0, MaxOmegaByT)[0])
        except OverflowError as e:
            print(str(e) + " for integration appears at tauT=" + str(tauT[i]))
    return CorrTrial


def chisq_dof(fit_params_0, xdata, ydata, edata, MaxOmegaByT, spfargs, Record):
    res = (ydata - TargetCorr(xdata, MaxOmegaByT, spfargs, fit_params_0)) / edata
    chisqdof = np.sum(res ** 2) / (len(xdata) - len(fit_params_0))
    record = list(fit_params_0)
    record.append(chisqdof)
    Record.append(record)
    print(record)
    return chisqdof


def main():
    """this script currently only does one bootstrap sample. so instead, wrap everything in bootstrap"""
    parser, requiredNamed = lpd.get_parser()

    requiredNamed.add_argument('--model', help='which model to use', choices=[1, 2, 3], type=int)
    requiredNamed.add_argument('--PathPhiUV', help='the full path of the input phiuv in omega/T phi/T^3', type=str)
    requiredNamed.add_argument('--PathOutputFolder', help='the path of output folder like /a/b', type=str)
    requiredNamed.add_argument('--ID', help='number we name this run, and also serve as seed', type=int)
    requiredNamed.add_argument('--mu', help='which "en" function to use', choices=["alpha", "beta"], type=str)
    requiredNamed.add_argument('--nmax', help='what nmax to use. valid only for model 1,2.', choices=[1, 2, 3, 4, 5, 6], type=int)
    requiredNamed.add_argument('--constrain', help='force the spf to reach the UV limit at large omega', action="store_true")

    args = parser.parse_args()

    MaxIter = 2000

    # read in the normalized correlator
    inputfolder = "../"+lpd.get_merged_data_path(args.qcdtype, args.corr, "")
    corr = np.genfromtxt(inputfolder+args.corr+"_final.txt", missing_values=None, usemask=True, invalid_raise=False)

    # remove lines with NaN's
    corr = corr[~np.isnan(corr).any(axis=1)]

    # beta, ns, nt, nt_half = lpd.parse_conftype(args.conftype)
    NtauT = len(corr)

    # read in the phiuv. columns in order: OmegaByT, PhiUVByT3, err
    PhiUV = np.loadtxt(args.PathPhiUV)
    PhiUV = PhiUV[:, 0:2]
    # interpolate and extrapolate the spf to everywhere for the integration in the next
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

    # resampling using bootstrp
    corr_sample = bootstr_from_gauss(ydata, edata, args.ID)

    # set up initial guess for the fitted parameters. the initial guess for kappa is 1, and for the c_n is 0. for model 3 we only have one other fit parameter
    # (the overall coefficient for the UV part), whose initial guess is also 1.
    if args.model == 3:
        fit_params_0 = [1, 1]
    else:
        fit_params_0 = [1]
        if args.constrain:
            args.nmax -= 1
        for i in range(args.nmax):
            fit_params_0.append(0.)

    Record = []
    spfargs = SpfArgs(args.model, args.mu, args.constrain, PhiuvByT3)  # constant parameters only used by the function SpfByT3

    # on the next line the fit happens
    # fun: The objective function to be minimized.
    # x0: Initial guess.
    # args: Extra arguments passed to the objective function
    res = scipy.optimize.minimize(fun=chisq_dof, x0=fit_params_0, args=(xdata, corr_sample, edata, MaxOmegaByT, spfargs, Record), method='Nelder-Mead',
                                  options={'maxiter': MaxIter, 'disp': True}, callback=None)

    # because in the fit we take the square of PhiIR, kappa can be negative. Since it is a square, the sign doesn't matter, so just take the absolute:
    if args.model == 2:
        res.x[0] = math.fabs(res.x[0])

    # use the fit results for the parameters to compute the fitted spectral function and correlator
    try:
        Spf = []
        for i in range(len(PhiUV)):
            Spf.append(SpfByT3(PhiUV[i][0], args.model, args.mu, PhiuvByT3, res.x))
        if any(n < 0 for n in Spf):
            print("negative spf")
            sys.exit(4)
        np.savetxt(args.PathOutputFolder + "/Spf_fitted_" + str(args.ID) + ".dat", np.column_stack((PhiUV[:, 0], Spf)), fmt='%16.15e')

        chisqdof = chisq_dof(res.x, xdata, corr_sample, edata, MaxOmegaByT, args.model, args.mu, PhiuvByT3, Record)  # Note: this chisq_dof is appended to the end of Record!

        Corr_fit = TargetCorr(xdata, args.model, args.mu, MaxOmegaByT, PhiuvByT3, res.x)
        for i in range(NtauT):
            Corr_fit[i] /= Gnorm(corr[i][0])
        resx = " ".join(map(str, res.x))
        write_final = open(args.PathOutputFolder + "/Final_collection_" + str(args.ID) + ".dat", "w+")
        write_final.write("#Fitted params:\n")
        write_final.write("%s\n" % resx)
        write_final.write("#Fit success:\n")
        write_final.write("%s\n" % res.success)
        write_final.write("#Fit status:\n")
        write_final.write("%s\n" % res.status)
        write_final.write("#Chisq_dof:\n")
        write_final.write("%e\n" % chisqdof)
        write_final.close()

        np.savetxt(args.PathOutputFolder + "/Record_" + str(args.ID) + ".dat", Record, fmt='%16.15e')

        np.savetxt(args.PathOutputFolder + "/Corr_original_fitted_" + str(args.ID) + ".dat", np.column_stack((xdata, corr[:, 1], corr[:, 2], Corr_fit)),
                   fmt='%16.15e')
    except:
        print("fit fails")  # this error message may be a bit misleading
        sys.exit(3)


if __name__ == '__main__':
    main()
    lpd.save_script_call()
