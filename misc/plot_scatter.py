#!/usr/bin/python3
import math
import numpy
from os import listdir
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from latqcdtools import jackknife

def get_color(myarray, i, start, end):
    return cm.gnuplot((myarray[i]-myarray[start])/(myarray[end-1]-myarray[start]))


try:
    qcdtype = sys.argv[1]
    conftype = sys.argv[2]
    beta = float(sys.argv[3])
    ns = int(sys.argv[4])
    nt = int(sys.argv[5])
except IndexError:
    exit("Invalid Arguments: merge_raw_data.py <qcdtype> <conftype> <beta> <ns> <nt>")

print(qcdtype, conftype, beta, ns, nt)
outputfolder="../../data_merged/"+qcdtype+"/"+conftype+"/"
inputfolder="../../data_raw/"+qcdtype+"/"+conftype+"/"

#---------
#load data
#---------
EE_numerator_real, EE_numerator_imag, polyakov_real, polyakov_imag, flow_times = ([] for i in range(5))
n_datafiles, n_streams = (0,0)

for stream_folder in listdir(inputfolder):
    if stream_folder.startswith(conftype+"_"):
        n_streams += 1
        for datafile in listdir(inputfolder+"/"+stream_folder):
            if datafile.startswith("zeuthenFlow"): #_acc0.000010_sts0.001000_ColElecCorrTimeSlices_"):#+conftype+"_U"):
                path = inputfolder+"/"+stream_folder+"/"+datafile
                tmp = numpy.loadtxt(path) # read in values as a numpy.ndarray
                if n_datafiles == 0:
                    flow_times = tmp[:,0]
                polyakov_real.append(tmp[:,1])
                polyakov_imag.append(tmp[:,2])
                EE_numerator_real.append(tmp[:,3:int((3+int(nt)/2))])
                EE_numerator_imag.append(tmp[:,int((4+int(nt)/2)):])
                n_datafiles += 1


plt.rc('text', usetex=True)
plt.rc('text.latex')
plt.rc('font', family='serif')

fig = plt.figure() 
ax = fig.add_subplot(1,1,1) 
labelboxstyle      = dict(boxstyle="round", fc="w", ec='none', alpha=0.7)
figurestyle        = dict(bbox_inches="tight")
ax.set_title(r'pure SU(3) $|$ $T\approx 1.5 T_C$ $|$ $'+str(ns)+'^3 \\times $'+str(nt)+' $|$ $\#_\mathrm{conf}=\,$'+str(n_datafiles)+' $|$ imp. dist.')
ax.xaxis.set_label_coords(0.975,0.050)
ax.yaxis.set_label_coords(0.025,0.975)
ax.set_ylabel(r'$\displaystyle \mathrm{Im}(\langle EE \rangle)$', 
              bbox=labelboxstyle, horizontalalignment='left', verticalalignment='top', rotation=0)
ax.set_xlabel(r'$\mathrm{Re}(\langle EE \rangle)$', horizontalalignment='right', bbox=labelboxstyle)
ax.axhline()
ax.axvline()

flowend=10
#poly_mean_real=numpy.empty(30)
#poly_mean_imag=numpy.empty(30)
#for i in range(0,n_datafiles):
    #for j in range(0,30):
        #poly_mean_real[j] += polyakov_real[i][j]
        #poly_mean_imag[j] += polyakov_imag[i][j]
#for j in range(0,30):
    #poly_mean_real[j] /= n_datafiles
    #poly_mean_imag[j] /= n_datafiles
#for j in range(0,flowend):
    #plt.scatter(poly_mean_imag[j], poly_mean_real[j], s=2, c=get_color(flow_times, j, 0, flowend))
        
#int(nt/2-)
for i in range(0,flowend):
    real=list(j[20,i] for j in EE_numerator_real[0:n_datafiles:10])
    imag=list(j[20,i] for j in EE_numerator_imag[0:n_datafiles:10])
    #real=list(j[i] for j in polyakov_real[0:n_datafiles:10])
    #imag=list(j[i] for j in polyakov_imag[0:n_datafiles:10])
    plt.scatter(real, imag,  s=0.05, label=str(flow_times[i]), c=get_color(flow_times, i, 0, flowend))

#ax.set_xlim([-0.00001,0.0001]); ax.set_ylim([-0.00005,0.0002])
#ax.set_xlim([-0.0001,0.001]); ax.set_ylim([-0.0005,0.002])
ax.legend()
fig.savefig("scatter.pdf", **figurestyle) 
