#!/usr/local/bin/python
import math
import numpy
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from matplotlib.legend_handler import HandlerErrorbar
from os import listdir
import sys
from latqcdtools import scales_quenched

qcdtype="quenched"
outputfolder="/home/altenkort/master/work/plots/"+qcdtype+"/"
inputfolder="/home/altenkort/master/work/data_merged/"+qcdtype+"/"

EE_16 = numpy.loadtxt(inputfolder+"/s096t16_b0719200/EE_s096t16_b0719200.dat")
EE_24 = numpy.loadtxt(inputfolder+"/s096t24_b0719200/EE_s096t24_b0719200.dat")
EE_28 = numpy.loadtxt(inputfolder+"/s096t28_b0719200/EE_s096t28_b0719200.dat")
EE_32 = numpy.loadtxt(inputfolder+"/s096t32_b0719200/EE_s096t32_b0719200.dat")
EE_err_16 = numpy.loadtxt(inputfolder+"/s096t16_b0719200/EE_err_s096t16_b0719200.dat")
EE_err_24 = numpy.loadtxt(inputfolder+"/s096t24_b0719200/EE_err_s096t24_b0719200.dat")
EE_err_28 = numpy.loadtxt(inputfolder+"/s096t28_b0719200/EE_err_s096t28_b0719200.dat")
EE_err_32 = numpy.loadtxt(inputfolder+"/s096t32_b0719200/EE_err_s096t32_b0719200.dat")

flow_radius_16 = list(numpy.around(numpy.loadtxt(inputfolder+"/s096t16_b0719200/flowradius_s096t16_b0719200.dat"), 6)) #need to change this to flow -times-
flow_times = [(16*i)**2 /8 for i in flow_radius_16] #they are the same for all these lattices

tauT_16 = list(numpy.around(numpy.loadtxt(inputfolder+"/s096t16_b0719200/tauT_imps096t16_b0719200.dat"), 4))
tauT_24 = list(numpy.around(numpy.loadtxt(inputfolder+"/s096t24_b0719200/tauT_imps096t24_b0719200.dat"), 4))
tauT_28 = list(numpy.around(numpy.loadtxt(inputfolder+"/s096t28_b0719200/tauT_imps096t28_b0719200.dat"), 4))
tauT_32 = list(numpy.around(numpy.loadtxt(inputfolder+"/s096t32_b0719200/tauT_imps096t32_b0719200.dat"), 4))


def free_corr( tauT ):
    norm = math.pi**2 * ( math.cos(math.pi*tauT)**2 / math.sin(math.pi*tauT)**4 + 1/(3*math.sin(math.pi*tauT)**2))
    return norm

def T_Tc(Ntau):
    return (1/0.7457 * numpy.exp(scales_quenched.ln_r0(7.192)) * 1/Ntau)
    #return (1/0.7457 * scales_quenched.a_r0_invGeV(7.192) * 1/Ntau)

def normalize_EE( data, data_err, tauT_imp):
    datashape = data.shape
    Ntau = datashape[1]*2
    for i in range(0, datashape[0]):
        for j in range(0, datashape[1]):
            data[i,j] = data[i,j] * free_corr(tauT_imp[j]) * T_Tc(Ntau)**4
            #if data[i,j] < 0:
                #data[i,j] = None
            data_err[i,j] = data_err[i,j] * free_corr(tauT_imp[j]) * T_Tc(Ntau)**4
            
def remove_negatives( data ):
    datashape = data.shape
    for i in range(0, datashape[0]):
        for j in range(0, datashape[1]):
            if data[i,j] < 0:
                data[i,j] = None

#/ Ntau**4 * 35.6345**4

#this actually removes the normalization factor and multiplies T/T_C
normalize_EE(EE_16, EE_err_16, tauT_16)
normalize_EE(EE_24, EE_err_24, tauT_24)
normalize_EE(EE_28, EE_err_28, tauT_28)
normalize_EE(EE_32, EE_err_32, tauT_32)
EE_16_bak = numpy.copy(EE_16)
EE_24_bak = numpy.copy(EE_24)
EE_28_bak = numpy.copy(EE_28)
EE_32_bak = numpy.copy(EE_32)

tauT_16 = numpy.arange(1/16,9/16,1/16)
tauT_24 = numpy.arange(1/24,13/24,1/24)
tauT_28 = numpy.arange(1/28,15/28,1/28)
tauT_32 = numpy.arange(1/32,17/32,1/32)
tauTc_16 = [numpy.around(i/T_Tc(16), 3) for i in tauT_16]
tauTc_24 = [numpy.around(i/T_Tc(24), 3) for i in tauT_24]
tauTc_28 = [numpy.around(i/T_Tc(28), 3) for i in tauT_28]
tauTc_32 = [numpy.around(i/T_Tc(32), 3) for i in tauT_32]
print(tauTc_16)
print(tauTc_24)
print(tauTc_28)
print(tauTc_32)

#-------------------
#set plot parameters
#-------------------
figurestyle        = dict(bbox_inches="tight")
plotstyle_points   = dict(fmt='D-', linewidth=1, markersize=4, capsize=2, mew=0.5, fillstyle='none')
legendstyle        = dict(loc="upper right", frameon=True, framealpha=0, title="", labelspacing=0.1, borderpad=0.1, handletextpad=0.4)
figurestyle        = dict(bbox_inches="tight", pad_inches=0.05)
labelboxstyle      = dict(boxstyle="square", fc="w", ec='none', alpha=0.7, pad=0.15, zorder=999999)
plotstyle_fill     = dict(linewidth=0.5)
xlabelstyle        = dict(horizontalalignment='right', verticalalignment='bottom',  bbox=labelboxstyle, zorder=999998)
ylabelstyle        = dict(horizontalalignment='left', verticalalignment='top', rotation=0, bbox=labelboxstyle, zorder=999999)



plt.rc('text', usetex=True)
plt.rc('text.latex')
plt.rc('font', family='serif', size=16)
fig = plt.figure() 
ax = fig.add_subplot(1,1,1) 
ax.xaxis.set_label_coords(0.975,0.025)
ax.yaxis.set_label_coords(0.01,0.97)
ax.set_xlabel(r'$\tau T_C$', **xlabelstyle)
ax.set_ylabel(r'$\displaystyle \frac{G_{\tau_F}(\tau T_C)}{T_C^4}$', **ylabelstyle)


ax.set_xlim([0,0.525])
ax.set_yscale('log')
ax.set_ylim([10,1000000])
plots = []

offset=0.02
flowstart=10
flowend=21

remove_negatives(EE_16)
remove_negatives(EE_24)
remove_negatives(EE_28)
remove_negatives(EE_32)


#PLOT LATTICE EFFECTS FOR ANIMATION
for i in range(0,len(flow_times)):
#for i in range(0,20):
#for i in (0,6,14,19):
    ax.set_title(r'$\tau_F =$ '+'{0:.3f}'.format(flow_times[i]), x=0.5, y=0.9, bbox=labelboxstyle, zorder=999999)
    plots.append(ax.errorbar(tauTc_32[:16], EE_32[i,:], EE_err_32[i,:], label=r'$T=1.11 T_C$', **plotstyle_points, color=cm.gnuplot(0.90), zorder=0))
    plots.append(ax.errorbar(tauTc_32[:14], EE_28[i,:], EE_err_28[i,:], label=r'$T=1.27 T_C$', **plotstyle_points, color=cm.gnuplot(0.65), zorder=-1))
    plots.append(ax.errorbar(tauTc_32[:12], EE_24[i,:], EE_err_24[i,:], label=r'$T-1.49 T_C$', **plotstyle_points, color=cm.gnuplot(0.40), zorder=-2))
    plots.append(ax.errorbar(tauTc_32[:8], EE_16[i,:], EE_err_16[i,:], label=r'$T=2.23 T_C$', **plotstyle_points, color=cm.gnuplot(0.03), zorder=-3))
    ax.legend(handles=plots, labels=[r'$T/T_C=1.11$', r'$T/T_C=1.27$', r'$T/T_C=1.49$', r'$T/T_C=2.23$'], **legendstyle, handler_map={type(plots[0]): HandlerErrorbar(xerr_size=0.4)}, handlelength=1.5)
    fig.savefig(outputfolder+"/single_flow_mass_shift/EE_flow_mass_shift_"+'{0:.4f}'.format(flow_times[i])+".pdf", **figurestyle) 
    ax.lines.clear() ; ax.collections.clear() ; plots.clear()

EE_16 = EE_16_bak
EE_24 = EE_24_bak
EE_28 = EE_28_bak
EE_32 = EE_32_bak


plt.rc('text', usetex=True)
plt.rc('text.latex')
plt.rc('font', family='serif', size=16)

fig = plt.figure() 
ax = fig.add_subplot(1,1,1) 
ax.xaxis.set_label_coords(0.975,0.025)
ax.yaxis.set_label_coords(0.01,0.97)
ax.set_xlabel(r'$\tau T_C$', **xlabelstyle)
#ax.set_ylabel(r'$\displaystyle \frac{ \left(G_{\tau_F}^{2.23\,T_C} - G_{\tau_F}^{1.12\,T_C}\right) (\tau T_C)}{T_C^4} $', **ylabelstyle)
ax.set_ylabel(r'')
ax.set_ylim([0,150])
#ax.set_yscale('log')
#ax.set_ylim([-1000,-2000])
for j in (0,6,14,19):
    ax.set_title(r'$\tau_F =$ '+'{0:.3f}'.format(flow_times[j]), x=0.98, y=0.955, bbox=labelboxstyle, horizontalalignment='right', verticalalignment='top', zorder=999999)
    plots.append(ax.errorbar(tauTc_32[:8], EE_16[j,:8]-EE_32[j,:8], (EE_err_16[j,:8]+EE_err_32[j,:8]), **plotstyle_points, zorder=-10, color='black'))
    #plots.append(ax.errorbar(tauTc_32[:8], EE_16[j,:8], EE_err_16[j,:8], **plotstyle_points, zorder=-10, color='black'))
    #plots.append(ax.errorbar(tauTc_32[:8], EE_32[j,:8], EE_err_32[j,:8], **plotstyle_points, zorder=-10, color='blue'))
    #ax.axvline(x=(0.5*1/T_Tc(16)), ymin=0, ymax=1, color='grey', alpha=0.8, zorder=-100, dashes=(4,4), lw=0.5)
    plots.append(ax.errorbar(tauTc_32[8:16]-tauTc_32[8]+tauTc_32[0], EE_32[j,8:16], EE_err_32[j,8:16], **plotstyle_points, zorder=-10, color='red'))
    ax.legend(handles=plots, labels=[r'$T_C^{-4} G_\Delta(\tau T_C)$', r'$T_C^{-4} G_1(\tau T_C+\Delta \tau T_C)$'], loc="upper left", bbox_to_anchor=(0.005,0.985), frameon=True, framealpha=0.5, edgecolor='none', fancybox=False, facecolor="w", title="", labelspacing=0.1, borderpad=0.1, handletextpad=0.4, handler_map={type(plots[0]): HandlerErrorbar(xerr_size=0.4)})
    fig.savefig(outputfolder+"/temp/EE_flow_mass_shift_diff_"+'{0:.4f}'.format(flow_times[j])+".pdf", **figurestyle) 
    ax.lines.clear() ; ax.collections.clear() ; plots.clear()

#NAIVE INTEGRAL
for j in (14,19): #the two flowtimes
    #Gdelta=tauTc_32[0] * (EE_16[j,0] - EE_32[j,0])
    Gdelta=0
    G1=0
    Gdelta_sum=0
    G1_sum=0
    for i in (0,1,2,3,4,5,6):
    #for i in (0,1,2,3,4,5,6,7):
        #Gdelta_sum+=(EE_16[j,i] - EE_32[j,i]) #* (tauTc_32[1]-tauTc_32[0])
        #G1_sum += EE_32[j,8+i] #* (tauTc_32[1]-tauTc_32[0])
        Gdelta += (tauTc_32[  i+1]-tauTc_32[  i]) * (EE_16[j,  i]-EE_32[j,  i  ] + EE_16[j,i+1]-EE_32[j,i+1]) / 2
        G1 +=     (tauTc_32[8+i+1]-tauTc_32[8+i]) * (EE_32[j,8+i]+EE_32[j,8+i+1])                             / 2
    #print("flow time = ", flow_times[j])
    #print("Gdelta = ", Gdelta, "G1 = ", G1)
    print("integral over -2*Gdelta + integral over 2*G1 = ", -2*Gdelta+2*G1)
    #print(Gdelta_sum, G1_sum)
    #print("naive = ", -2*Gdelta_sum+2*G1_sum)
#ax.set_xlim([0])
#ax.set_yscale('log')
#plots.append(ax.errorbar(flow_times, EE_16[:,i]-EE_32[:,i], (EE_err_16[:, i] + EE_err_32[:,i]), **plotstyle_points))
#print("tauTc=", tauT_16[i]/T_Tc(16))
#print("t_F      G(tauTc, t_F)             rel. err in %")
#for i in range(0,4):
    #print('{0:.4f}'.format(tauT_16[i]), EE_16[j, :] - EE_32[j,:], (EE_err_16[j, :] + EE_err_32[j,:])/(EE_16[j, i] - EE_32[j,i])*100)
#print("\n")
