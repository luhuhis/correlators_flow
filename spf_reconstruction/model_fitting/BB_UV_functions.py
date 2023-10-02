from scipy import integrate
from scipy.interpolate import CubicSpline
import rundec
import numpy as np
import lib_process_data as lpd

def nb(x):
    if x > 200:
        return 0
    else:
        return 1. / (np.exp(x) - 1.)


def nf(x):
    if x > 200:
        return 0
    else:
        return 1. / (np.exp(x) + 1.)


def ColorPart(OmegaByT, qByT):
    X = (qByT ** 2 + 2 * OmegaByT ** 2) * np.log(abs((qByT + OmegaByT) / (qByT - OmegaByT)))
    Y = qByT * OmegaByT * (np.log(abs(qByT ** 2 - OmegaByT ** 2) / OmegaByT ** 2) + OmegaByT ** 2 / (
                OmegaByT ** 2 - qByT ** 2) - 2)
    Z = OmegaByT ** 4 / qByT * (
                1. / (OmegaByT + qByT) * np.log((OmegaByT + qByT) / OmegaByT) + 1 / (qByT - OmegaByT) * np.log(
            OmegaByT / abs(qByT - OmegaByT)))
    return (X + Y + Z) * nb(qByT)


def FlavorPart(OmegaByT, qByT):
    X = (qByT ** 2 + OmegaByT ** 2 / 2) * np.log(abs((qByT + OmegaByT) / (qByT - OmegaByT)))
    Y = qByT * OmegaByT * (np.log(abs(qByT ** 2 - OmegaByT ** 2) / OmegaByT ** 2) - 1)
    return (X + Y) * nf(qByT)


def get_mu_free(OmegaByT, Nc, Nf):
    if OmegaByT <= 1e-6:
        OmegaByT = 1e-6
    colorPart = \
    integrate.quad(lambda x: ColorPart(OmegaByT, x) + ColorPart(OmegaByT, 2 * OmegaByT - x), OmegaByT, 2 * OmegaByT,
                   epsabs=0)[0] + integrate.quad(lambda x: ColorPart(OmegaByT, x), 2 * OmegaByT, 200, epsabs=0)[0]
    if Nf == 0:
        return 0 if Nc * colorPart > 0 else Nc * colorPart
    else:
        flavorPart = \
        integrate.quad(lambda x: FlavorPart(OmegaByT, x) + FlavorPart(OmegaByT, 2 * OmegaByT - x), OmegaByT,
                       2 * OmegaByT, epsabs=0)[0] + \
        integrate.quad(lambda x: FlavorPart(OmegaByT, x), 2 * OmegaByT, 200, epsabs=0)[0]
        return 0 if Nc * colorPart + Nf * flavorPart > 0 else Nc * colorPart + Nf * flavorPart


def get_cBsq(ir_scale, T_in_GeV, Npoints, gamma_0, Lambda_MSbar, Nf, Nloop):
    crd = rundec.CRunDec()
    g2_arr = []
    for OmegaByT in np.logspace(-6, 3, Npoints, base=10):
        mu = np.sqrt(4 * (OmegaByT * T_in_GeV) ** 2. + ir_scale ** 2)
        Alphas = crd.AlphasLam(Lambda_MSbar, mu, Nf, Nloop)
        g2 = 4. * np.pi * Alphas
        g2_arr.append([mu, g2])
    g2_arr = np.array(g2_arr)

    threshold = 1e-10
    diff = np.empty(len(g2_arr))
    diff[0] = np.inf  # always retain the 1st element
    diff[1:] = np.diff(g2_arr[:, 0])
    mask = diff > threshold
    new_arr = g2_arr[mask]

    spl = CubicSpline(new_arr[:, 0], new_arr[:, 1], axis=0, bc_type=((2, 0), (2, 0)), extrapolate=True)
    cBsq = np.asarray(lpd.parallel_function_eval(calc, range(len(g2_arr)), 128, g2_arr, spl, gamma_0))
    return cBsq


def calc(index, g2_arr, spl, gamma_0):
    this_x = g2_arr[:index, 0]
    this_y = 2. / np.array(this_x) * (gamma_0 * spl(this_x))
    integral = integrate.trapz(this_y, this_x)
    return np.exp(integral)