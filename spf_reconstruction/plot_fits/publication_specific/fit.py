import concurrent.futures
from functools import partial
from numpy.random import randint
from numpy.random import normal
from scipy.interpolate import CubicSpline
from scipy import integrate
import numpy as np
from lmfit import minimize, Parameters, fit_report
from scipy.optimize import minimize as fastfit
import sys
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib import cm
from numpy import linspace
import matplotlib.ticker as mticker

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
    if flag==1:
        return a*g2**2
    elif flag==2:
        return a*g2
    else:
        print("no such model.exit")
        sys.exit(1)

def residual(fit_params_0, g2, y, yerr, flag):
    return (fun(fit_params_0['a'], g2, flag)-y)/yerr

def fit(g2, sample, yerr, flag):
    fit_params_0 = Parameters()
    fit_params_0.add('a', min=0, max=10,value=0.2)
    out = minimize(residual, fit_params_0, method="leastsq", kws={"g2": g2, "y": sample, "yerr": yerr, "flag": flag}, scale_covar=False)
    #print(fit_report(out))
    fits=[]
    for name, param in out.params.items():
        fits.append(param.value) 
    fits.append(out.redchi)
    return fits

def fit_direct(g2, y, yerr, flag):
    fit_params_0 = Parameters()
    fit_params_0.add('a', min=0, max=10,value=0.2)
    out = minimize(residual, fit_params_0, method="leastsq", kws={"g2": g2, "y": y, "yerr": yerr, "flag": flag}, scale_covar=False)
    #print(fit_report(out))
    for name, param in out.params.items():
        value=param.value
        error=param.stderr
    return value, error

Data=[]
quenched_kappa_B=np.loadtxt("quenched_kappa_B.dat")
for ii in range(len(quenched_kappa_B)):
    Data.append([quenched_kappa_B[ii,1], (quenched_kappa_B[ii,2]+quenched_kappa_B[ii,3])/2, (quenched_kappa_B[ii,3]-quenched_kappa_B[ii,2])/2.])

quenched_kappa_E=np.loadtxt("quenched_kappa_E.dat")
for ii in range(len(quenched_kappa_E)):
    Data.append([quenched_kappa_E[ii,1], (quenched_kappa_E[ii,2]+quenched_kappa_E[ii,3])/2, (quenched_kappa_E[ii,3]-quenched_kappa_E[ii,2])/2.])


g2_qcd = np.loadtxt("g2_qcd.dat")
kappa_qcd_E=np.loadtxt("kappa_qcd_E.dat")
for ii in range(len(g2_qcd)):
    Data.append([g2_qcd[ii, 1], kappa_qcd_E[ii,1], kappa_qcd_E[ii,2]])

kappa_qcd_B=np.loadtxt("kappa_qcd_B.dat")
for ii in range(len(g2_qcd)):
    Data.append([g2_qcd[ii, 1], kappa_qcd_B[ii,1], (kappa_qcd_B[ii,2]+kappa_qcd_B[ii,3])/2.])
Data=np.array(Data)
data = Data[Data[:, 0].argsort()]
fit_data=[]
g2=data[:,0]
y=data[:,1]
yerr=data[:,2]

Nsamples = 200
Samples = []
for s in range(Nsamples):
    Samples.append(np.array(bootstr_from_gauss(g2, y, yerr, s)))
Samples=np.array(Samples) #Nsamples x Nd

FitSample=[]
for s in range(Nsamples):
    FitSample.append(fit(g2, Samples[s], yerr, 1))
FitSample = np.array(FitSample)
value=np.mean(FitSample,axis=0)
error=np.std(FitSample,axis=0)
print("[mean, error] from gaussian bootstrap: ", value[0],error[0])

value, error = fit_direct(g2, y, yerr, 1)
print("[mean, error] from direct fit: ", value,error)

value2, error2 = fit_direct(g2, y, yerr, 2)

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'], 'size':14})
plt.rc('text', usetex=True)

plt.figure(figsize=(4,4))
plt.xscale("log")
plt.yscale("log")

axes = plt.gca()
axes.set_xlim([1,6])
axes.set_ylim([2e-1,2e1])
axes.xaxis.set_minor_formatter(mticker.ScalarFormatter())


axes.errorbar(quenched_kappa_E[:,1], (quenched_kappa_E[:,2]+quenched_kappa_E[:,3])/2., (quenched_kappa_E[:,3]-quenched_kappa_E[:,2])/2., color='red', fmt='o', fillstyle='none',capsize=2.0,linewidth=1.5, markersize=5, label=r'$\kappa_E/T^3,\ N_f=0$')
axes.errorbar(quenched_kappa_B[:,1], (quenched_kappa_B[:,2]+quenched_kappa_B[:,3])/2., (quenched_kappa_B[:,3]-quenched_kappa_B[:,2])/2., color='blue', fmt='o', fillstyle='none',capsize=2.0,linewidth=1.5, markersize=5, label=r'$\kappa_B/T^3,\ N_f=0$')

axes.errorbar(g2_qcd[:,1], kappa_qcd_E[:,1], kappa_qcd_E[:,2], color='red', fmt='o', capsize=2.0,linewidth=1.5, markersize=6, ecolor="red",label=r'$\kappa_E/T^3,\ N_f=2+1$')
axes.errorbar(np.array(g2_qcd[:,1])+0.06, kappa_qcd_B[:,1], (kappa_qcd_B[:,2]+kappa_qcd_B[:,3])/2., color='blue', fmt='o', capsize=2.0,linewidth=1.5, markersize=6, ecolor="blue", label=r'$\kappa_B/T^3,\ N_f=2+1$')

x = np.linspace(1,6,1000)
plt.plot(x,value*x**2, linewidth=1, linestyle='solid', color='magenta', label='')
plt.fill_between(x, (value-error)*x**2, (value+error)*x**2,  facecolor = 'magenta', edgecolor="none", alpha = 0.4, label=r'')

plt.plot(x,value2*x, linewidth=1, linestyle='solid', color='green', label='')
plt.fill_between(x, (value2-error2)*x, (value2+error2)*x,  facecolor = 'green', edgecolor="none", alpha = 0.4, label=r'')

plt.xlabel(r'$g^2(\mu=2\pi T)$', fontsize=12)

handles, labels = axes.get_legend_handles_labels()
lgd=plt.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.25, 1.0), prop={'size': 9}, frameon=True, handletextpad=0.1)
plt.savefig("compare_kappa_g2.pdf",  bbox_extra_artists=(lgd,), bbox_inches='tight')
