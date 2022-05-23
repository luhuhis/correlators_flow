#!/usr/local/bin/python
import numpy as np


# forwards
def conv_tF_to_tauFT2(tF, Ntau):
    tF = np.asarray(tF)
    return tF / Ntau ** 2


def conv_tauFT2_to_sqrt8tauFT(tauFT2):
    tauFT2 = np.asarray(tauFT2)
    return np.sqrt(8 * tauFT2)


def conv_tF_to_sqrt8tauFT(tF, Ntau):
    return conv_tauFT2_to_sqrt8tauFT(conv_tF_to_tauFT2(tF, Ntau))


# backwards
def conv_sqrt8tauFT_to_tauFT2(sqrt8tauFT):
    sqrt8tauFT = np.asarray(sqrt8tauFT)
    return sqrt8tauFT ** 2 / 8


def conv_tauFT2_to_tF(tauFT2, Ntau):
    tauFT2 = np.asarray(tauFT2)
    return tauFT2 * Ntau ** 2


def conv_sqrt8tauFT_to_tF(sqrt8tauFT, Ntau):
    sqrt8tauFT = np.asarray(sqrt8tauFT)
    return conv_tauFT2_to_tF(conv_sqrt8tauFT_to_tauFT2(sqrt8tauFT), Ntau)


def stepsizes(values):
    return np.asarray([values[i + 1] - val for i, val in enumerate(values) if values[i] != values[-1]])


# proposed stepsize

Ntau = 64

max_stepsize = 0.075  # 0.15
max_sqrt8tauFT = 0.5  # 0.3

under_150 = 0
# proposed_sqrt8tauFT = [0,conv_tF_to_sqrt8tauFT(0.001,Ntau)]
proposed_sqrt8tauFT = [0, 0.0001]  # 0.001
current_sqrt8tauFT = proposed_sqrt8tauFT[-1]  # inital step
current_stepsize = 0.0004  # 0.001.  start stepsize (after first step)

# while conv_sqrt8tauFT_to_tF(proposed_sqrt8tauFT[-1],Ntau) - conv_sqrt8tauFT_to_tF(proposed_sqrt8tauFT[-2],Ntau) < max_stepsize and current_sqrt8tauFT < max_sqrt8tauFT:
while current_sqrt8tauFT <= max_sqrt8tauFT:  # use this and fixed_stepsize=0.15 to get the output that you can put into the gradientFlow
    current_sqrt8tauFT += current_stepsize
    proposed_sqrt8tauFT.append(current_sqrt8tauFT)
    current_stepsize *= 1.01
    if current_sqrt8tauFT < 0.150:
        under_150 += 1

proposed_sqrt8tauFT = np.asarray(proposed_sqrt8tauFT[:-1])
proposed_sqrt8tauFT = np.append(proposed_sqrt8tauFT, max_sqrt8tauFT)
proposed_tF = conv_sqrt8tauFT_to_tF(proposed_sqrt8tauFT, Ntau)

while proposed_tF[-1] < conv_sqrt8tauFT_to_tF(max_sqrt8tauFT, Ntau):
    print(proposed_tF[-1], conv_sqrt8tauFT_to_tF(0.5,Ntau))
    proposed_tF = np.append(proposed_tF, proposed_tF[-1]+max_stepsize)


proposed_sqrt8tauFT = conv_tF_to_sqrt8tauFT(proposed_tF[:-1],Ntau)
proposed_sqrt8tauFT = np.append(proposed_sqrt8tauFT, max_sqrt8tauFT)


print(proposed_sqrt8tauFT)
print(proposed_tF)

np.savetxt("flowradii_fine.dat", proposed_sqrt8tauFT, fmt='%.8f')
np.savetxt("flowtimes_fine.dat", proposed_tF, fmt='%.8f')  # newline=' ',

