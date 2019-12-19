#!/usr/bin/python3
import math
import numpy
import scipy.optimize as opt

def free_corr( tauT ):
    return math.pi**2 * ( math.cos(math.pi*tauT)**2 / math.sin(math.pi*tauT)**4 + 1/(3*math.sin(math.pi*tauT)**2))

def free_corr_flow( tauT, flowtime, Nt ):
    result = 0
    for n in range(-50,50):
        x = tauT*Nt + n*Nt
        xi_sq = x**2 / (8 * flowtime)
        result += 1/x**4  * ((xi_sq**2 + xi_sq + 1)*math.exp(-xi_sq) - 1)
    result *= - 1 / math.pi**2 * Nt**4
    return result

Nt = 16
tauT_imp_Nt16 = (6.883194e-02, 1.085070e-01, 1.616768e-01, 2.247303e-01, 2.914720e-01, 3.571978e-01, 4.178695e-01, 4.527297e-01)
flowtimes = (.00080, .00320, .00720, .01280, .02000, .02880, .03920, .05120, .06480, .08000, .09680, .11520, .13520, .15680, .18000, .20480, .23120, .25920, .28880, .32000, .50000, .72000, 1.28000, 2.00000, 2.88000, 3.92000, 5.12000, 6.48000, 8.00000)

print("for unitless flow time t = 0.0008 (--> tau_F = t * a^2):")
print("   tauT      G_free G_free_flow")
for i in range(0,8):
    print('{0: >.5f} {1: >11.5f} {2: >11.5f}'.format(tauT_imp_Nt16[i], free_corr(tauT_imp_Nt16[i]), free_corr_flow( tauT_imp_Nt16[i], flowtimes[0], Nt ) ))


print(opt.root_scalar(free_corr, bracket=[-0.5,0.5]))
