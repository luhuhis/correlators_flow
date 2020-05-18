import numpy
import matplotlib
from matplotlib import pyplot, cm, container, legend_handler
import lib_process_data as lpd

def get_color(myarray, i, start, end):
    return matplotlib.cm.gnuplot((myarray[i]-myarray[start])/(myarray[end-1]-myarray[start]))

#styles
mylinewidth        = 0.5
figurestyle        = dict(bbox_inches="tight", pad_inches=0)
plotstyle_points   = dict(fmt='D-', linewidth=1, markersize=4, capsize=2, mew=0.5, fillstyle='none')
plotstyle_lines    = dict(fmt='-', linewidth=0.5, mew=0.25, mec='grey', markersize=2, capsize=3)
plotstyle_add_point= dict(fmt='D-', linewidth=0.5, mew=0.25, mec='grey', markersize=2, capsize=3)
plotstyle_add_point_single= dict(fmt='D', linewidth=1, mew=0.5, markersize=3, capsize=3, fillstyle='none')
labelboxstyle      = dict(boxstyle="Round", fc="w", ec="None", alpha=0.8, pad=0.1, zorder=99)
legendstyle        = dict(loc="center left", bbox_to_anchor=(1,0.5), frameon=True, framealpha=0.8, edgecolor='none', fancybox=False, facecolor="w", labelspacing=0.1, borderpad=0.1, handletextpad=0.4, handlelength=1, handler_map={matplotlib.container.ErrorbarContainer: matplotlib.legend_handler.HandlerErrorbar(xerr_size=1)})
plotstyle_fill     = dict(linewidth=0.5)
xlabelstyle        = dict(bbox=labelboxstyle, zorder=99)#horizontalalignment='right', verticalalignment='bottom', bbox=labelboxstyle)
ylabelstyle        = dict(bbox=labelboxstyle, horizontalalignment='right', rotation=0, zorder=99)#horizontalalignment='right', verticalalignment='top', rotation=0, bbox=labelboxstyle)
titlestyle         = dict(x=0.5, y=0.9, bbox=labelboxstyle, verticalalignment='top', zorder=99)
verticallinestyle  = dict(ymin=0, ymax=1, color='grey', alpha=0.8, zorder=-100, dashes=(4,4), lw=0.5)
flowlimitplotstyle = dict(fmt='|', mew=0.7, markersize=15)

"""some parameters"""
n_skiprows=15 #skip the discrete points at the start of the continuum limit file
offset=0.02 #for flow limits
#index for tauT
start=2
end=14
#index for tau_F
flowstart=10
flowend=21

def create_figure(xlims=None, ylims=None, xlabel="", ylabel="", UseTex = True):
    if UseTex == True:
        matplotlib.pyplot.rc('text', usetex=True)
        matplotlib.pyplot.rc('font', family='serif', size=14)
    fig = matplotlib.pyplot.figure() 
    ax = fig.add_subplot(1,1,1) 
    #ax.xaxis.set_label_coords(0.99,0.01)
    #ax.yaxis.set_label_coords(0.01,0.97)
    if xlims is not None:
        ax.set_xlim(xlims)
    if ylims is not None:
        ax.set_ylim(ylims)
    ax.set_xlabel(xlabel, **xlabelstyle)
    ax.set_ylabel(ylabel, **ylabelstyle)
    plots = []
    return fig, ax, plots


"""load data"""
qcdtype="quenched"
outputfolder="../../plots/"+qcdtype+"/"
inputfolder="../../data_merged/"+qcdtype+"/"

EE_16 = numpy.loadtxt(inputfolder+"/s064t16_b0687361/EE_s064t16_b0687361.dat")
EE_20 = numpy.loadtxt(inputfolder+"/s080t20_b0703500/EE_s080t20_b0703500.dat")
EE_24 = numpy.loadtxt(inputfolder+"/s096t24_b0719200/EE_s096t24_b0719200.dat")
EE_30 = numpy.loadtxt(inputfolder+"/s120t30_b0739400/EE_s120t30_b0739400.dat")
EE_36 = numpy.loadtxt(inputfolder+"/s144t36_b0754400/EE_s144t36_b0754400.dat")
EE = [EE_16, EE_20, EE_24, EE_30, EE_36]
EE_err_16 = numpy.loadtxt(inputfolder+"/s064t16_b0687361/EE_err_s064t16_b0687361.dat")
EE_err_20 = numpy.loadtxt(inputfolder+"/s080t20_b0703500/EE_err_s080t20_b0703500.dat")
EE_err_24 = numpy.loadtxt(inputfolder+"/s096t24_b0719200/EE_err_s096t24_b0719200.dat")
EE_err_30 = numpy.loadtxt(inputfolder+"/s120t30_b0739400/EE_err_s120t30_b0739400.dat")
EE_err_36 = numpy.loadtxt(inputfolder+"/s144t36_b0754400/EE_err_s144t36_b0754400.dat")
EE_err = [EE_err_16, EE_err_20, EE_err_24, EE_err_30, EE_err_36]
tauT_16 = lpd.tauT_imp[16]
tauT_20 = lpd.tauT_imp[20]
tauT_24 = lpd.tauT_imp[24]
tauT_30 = lpd.tauT_imp[30]
tauT_36 = lpd.tauT_imp[36]
tauT = [tauT_16, tauT_20, tauT_24, tauT_30, tauT_36]
#tauT_30_ext = (*tauT_30,0.5)
tauT_30_ext = (0.229167, 0.250000, 0.270833, 0.291667, 0.312500, 0.333333, 0.354167, 0.375000, 0.395833, 0.416667, 0.437500, 0.458333, 0.479167, 0.500000)
flow_radius = numpy.loadtxt(inputfolder+"/s064t16_b0687361/flowradii_s064t16_b0687361.dat")
