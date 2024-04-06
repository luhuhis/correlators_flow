#!/usr/bin/env python3

from numpy.random import normal
import numpy as np
from lmfit import minimize, Parameters, fit_report
import sys
import matplotlib.ticker as mticker
import lib_process_data as lpd
from matplotlib.ticker import LogLocator
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('inputpath', default="/work/home/altenkort/work/correlators_flow/plots/quenched_1.50Tc_zeuthenFlow/coupling/")
parser.add_argument('outputpath', default="/work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/coupling/")

args = parser.parse_args()


def bootstr_from_gauss(x, data, data_err, Seed):
    data = np.asarray(data)
    data_err = np.asarray(data_err)
    np.random.seed(Seed)
    numb_observe = len(data)
    sampleval = []
    for k in range(numb_observe):
        sampleval.append(normal(data[k], data_err[k]))
    sampleval = np.array(sampleval)
    return sampleval


def fun(a, g2, flag):
    if flag == 1:
        return a * g2 ** 2
    elif flag == 2:
        return a * g2
    else:
        print("no such model.exit")
        sys.exit(1)


def residual(fit_params_0, g2, y, yerr, flag):
    return (fun(fit_params_0['a'], g2, flag) - y) / yerr


def fit(g2, sample, yerr, flag):
    fit_params_0 = Parameters()
    fit_params_0.add('a', min=0, max=10, value=0.2)
    out = minimize(residual, fit_params_0, method="leastsq", kws={"g2": g2, "y": sample, "yerr": yerr, "flag": flag}, scale_covar=False)
    # print(fit_report(out))
    fits = []
    for name, param in out.params.items():
        fits.append(param.value)
    fits.append(out.redchi)
    return fits


def fit_direct(g2, y, yerr, flag):
    fit_params_0 = Parameters()
    fit_params_0.add('a', min=0, max=10, value=0.2)
    out = minimize(residual, fit_params_0, method="leastsq", kws={"g2": g2, "y": y, "yerr": yerr, "flag": flag}, scale_covar=False)
    print(fit_report(out))
    for name, param in out.params.items():
        value = param.value
        error = param.stderr
    return value, error


Data = []
file = args.inputpath + "/quenched_kappa_B.dat"
print("load", file)
quenched_kappa_B = np.loadtxt(file)
for ii in range(len(quenched_kappa_B)):
    Data.append(
        [quenched_kappa_B[ii, 1], (quenched_kappa_B[ii, 2] + quenched_kappa_B[ii, 3]) / 2, (quenched_kappa_B[ii, 3] - quenched_kappa_B[ii, 2]) / 2.])

file = args.inputpath + "/quenched_kappa_E.dat"
print("load", file)
quenched_kappa_E = np.loadtxt(file)
for ii in range(len(quenched_kappa_E)):
    Data.append(
        [quenched_kappa_E[ii, 1], (quenched_kappa_E[ii, 2] + quenched_kappa_E[ii, 3]) / 2, (quenched_kappa_E[ii, 3] - quenched_kappa_E[ii, 2]) / 2.])

file = args.inputpath + "/g2_qcd.dat"
print("load", file)
g2_qcd = np.loadtxt(file)

file = args.inputpath + "/kappa_qcd_E.dat"
print("load", file)
kappa_qcd_E = np.loadtxt(file)
for ii in range(len(g2_qcd)):
    Data.append([g2_qcd[ii, 1], kappa_qcd_E[ii, 1], kappa_qcd_E[ii, 2]])

file = args.inputpath + "/kappa_qcd_B.dat"
print("load", file)
kappa_qcd_B = np.loadtxt(file)
for ii in range(len(g2_qcd)):
    Data.append([g2_qcd[ii, 1], kappa_qcd_B[ii, 1], (kappa_qcd_B[ii, 2] + kappa_qcd_B[ii, 3]) / 2.])
Data = np.array(Data)
data = Data[Data[:, 0].argsort()]
fit_data = []
g2 = data[:, 0]
y = data[:, 1]
yerr = data[:, 2]

Nsamples = 200
Samples = []
for s in range(Nsamples):
    Samples.append(np.array(bootstr_from_gauss(g2, y, yerr, s)))
Samples = np.array(Samples)  # Nsamples x Nd

FitSample = []
for s in range(Nsamples):
    FitSample.append(fit(g2, Samples[s], yerr, 1))
FitSample = np.array(FitSample)
value = np.mean(FitSample, axis=0)
error = np.std(FitSample, axis=0)
print("[mean, error] from gaussian bootstrap: ", value[0], error[0])

value, error = fit_direct(g2, y, yerr, 1)
print("[mean, error] from direct fit: ", value, error)

value2, error2 = fit_direct(g2, y, yerr, 2)

# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'], 'size':14})
# plt.rc('text', usetex=True)


fig, ax, ax_twiny = lpd.create_figure(xlims=[0.5, 6], ylims=[4e-2, 2e1], xlabel=r'$g^2(\mu=2\pi T)$')

ax.set_xscale("log")
ax.set_yscale("log")
ax_twiny.set_yscale("log")

ax.xaxis.set_minor_formatter(mticker.ScalarFormatter())
ax.xaxis.set_minor_locator(LogLocator(base=2))
ax.xaxis.set_major_formatter(mticker.StrMethodFormatter("${x:.2f}$"))
# ax.set_xticks = [0.5, 1, 2, 4, 6]

redcolor = 'tab:red'
bluecolor = 'tab:blue'
magentacolor = 'tab:pink'
greencolor = 'tab:green'

redcolor = 'C0'
bluecolor = 'C2'
magentacolor = 'C1'
greencolor = 'C3'

ms = 3.5
lw = 1
cs = 1.5
mew = 0.75

ax.errorbar(quenched_kappa_E[:, 1], (quenched_kappa_E[:, 2] + quenched_kappa_E[:, 3]) / 2., (quenched_kappa_E[:, 3] - quenched_kappa_E[:, 2]) / 2.,
            color=redcolor, fmt='o', fillstyle='none', capsize=cs, linewidth=lw, markersize=ms, label=r'$\kappa_E/T^3,\ N_f=0$', zorder=2, mew=mew)
ax.errorbar(quenched_kappa_B[:, 1], (quenched_kappa_B[:, 2] + quenched_kappa_B[:, 3]) / 2., (quenched_kappa_B[:, 3] - quenched_kappa_B[:, 2]) / 2.,
            color=bluecolor, fmt='o', fillstyle='none', capsize=cs, linewidth=lw, markersize=ms, label=r'$\kappa_B/T^3,\ N_f=0$', zorder=3, mew=mew)

ax.errorbar(g2_qcd[:, 1], kappa_qcd_E[:, 1], kappa_qcd_E[:, 2], color=redcolor, fmt='o', capsize=cs, linewidth=lw, markersize=ms, ecolor=redcolor,
            label=r'$\kappa_E/T^3,\ N_f=2+1$', zorder=4, mew=mew)
ax.errorbar(np.array(g2_qcd[:, 1]), kappa_qcd_B[:, 1], (kappa_qcd_B[:, 2] + kappa_qcd_B[:, 3]) / 2., color=bluecolor, fmt='o', capsize=cs,
            linewidth=lw, markersize=ms, ecolor=bluecolor, label=r'$\kappa_B/T^3,\ N_f=2+1$', zorder=5, mew=mew)

x = np.linspace(0.5, 6, 1000)
ax.errorbar(x, value * x ** 2, linewidth=lw, linestyle='solid', color=magentacolor, label='', zorder=-10)
ax.fill_between(x, (value - error) * x ** 2, (value + error) * x ** 2, facecolor=magentacolor, edgecolor="none", alpha=0.4, label=r'', zorder=-11)

ax.errorbar(x, value2 * x, linewidth=lw, linestyle='solid', color=greencolor, label='', zorder=-12)
ax.fill_between(x, (value2 - error2) * x, (value2 + error2) * x, facecolor=greencolor, edgecolor="none", alpha=0.4, label=r'', zorder=-13)

handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, loc='upper left', bbox_to_anchor=(0, 1.0), fontsize=9)

file = args.outputpath + "/compare_kappa_g2.pdf"
print("save", file)
fig.savefig(file)
