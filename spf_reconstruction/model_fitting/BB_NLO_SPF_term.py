from scipy import integrate
import numpy as np


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



