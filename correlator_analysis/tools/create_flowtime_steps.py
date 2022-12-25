#!/usr/env/bin python3
import numpy as np
import lib_process_data as lpd
from matplotlib import pyplot


def conv_tF_to_tauFT2(tF, Ntau):
    tF = np.asarray(tF)
    return tF / Ntau ** 2


def conv_tauFT2_to_sqrt8tauFT(tauFT2):
    tauFT2 = np.asarray(tauFT2)
    return np.sqrt(8 * tauFT2)


def conv_tF_to_sqrt8tauFT(tF, Ntau):
    return conv_tauFT2_to_sqrt8tauFT(conv_tF_to_tauFT2(tF, Ntau))


def conv_tF_to_sqrt8tF(tF):
    return np.sqrt(8*tF)


def conv_sqrt8tauFT_to_tauFT2(sqrt8tauFT):
    sqrt8tauFT = np.asarray(sqrt8tauFT)
    return sqrt8tauFT ** 2 / 8


def conv_tauFT2_to_tF(tauFT2, Ntau):
    tauFT2 = np.asarray(tauFT2)
    return tauFT2 * Ntau ** 2


def conv_sqrt8tauFT_to_tF(sqrt8tauFT, Ntau):
    sqrt8tauFT = np.asarray(sqrt8tauFT)
    return conv_tauFT2_to_tF(conv_sqrt8tauFT_to_tauFT2(sqrt8tauFT), Ntau)


def conv_sqrt8tF_to_tF(sqrt8tF):
    return sqrt8tF**2/8


def stepsizes(values):
    return np.asarray([values[i + 1] - val for i, val in enumerate(values) if values[i] != values[-1]])


def get_flowtimes(Ntau, max_sqrt8tF=None, first_step_in_sqrt8tF=0.1, stepsize_increment_rate=1.03, max_stepsize_in_tF=0.15):
    # this creates flowtime steps that are similar (but with slighlty smaller stepsizes) to what the adaptive stepsize algorithm would
    # produce with accuracy=1e-05 and initial stepsize=1e-04
    if max_sqrt8tF is None:
        max_sqrt8tF = Ntau/2

    proposed_tF = [0,  0.0001, conv_sqrt8tF_to_tF(first_step_in_sqrt8tF)]
    current_tF = proposed_tF[-1]
    current_stepsize_in_tF = conv_sqrt8tF_to_tF(first_step_in_sqrt8tF)
    while current_tF < conv_sqrt8tF_to_tF(max_sqrt8tF):
        # update stepsize
        if current_stepsize_in_tF < max_stepsize_in_tF:
            current_stepsize_in_tF = conv_sqrt8tF_to_tF(conv_tF_to_sqrt8tF(current_stepsize_in_tF)*stepsize_increment_rate)
            if current_stepsize_in_tF > max_stepsize_in_tF:
                current_stepsize_in_tF = max_stepsize_in_tF
        current_tF += current_stepsize_in_tF
        proposed_tF.append(float(np.format_float_positional(current_tF, precision=4, unique=False, fractional=False, trim='k')))  # round to 4 significant digits.

    # remove lastly added tF, which is too large
    tF = np.asarray([*proposed_tF[:-1], conv_sqrt8tF_to_tF(max_sqrt8tF)])

    return tF


Ntau = 20
tF = get_flowtimes(20)
fig, ax, _ = lpd.create_figure()

ax.set_yscale('log')
ax.errorbar(range(len(tF)), tF, fmt='x')
pyplot.show()


np.savetxt("flowtimes_Nt"+str(Ntau)+".dat", tF, fmt='%.8f')