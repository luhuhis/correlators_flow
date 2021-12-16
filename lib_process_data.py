import numpy
import re
import math
import argparse
from collections import ChainMap
import matplotlib
from matplotlib import pyplot, cm, container, legend_handler
import sys


def save_script_call(add_folder=None):
    """ save full script call in logfile """
    import datetime
    now = datetime.datetime.now()
    dt_string = now.strftime("%Y/%m/%d %H:%M:%S")
    import os
    file_name = os.path.basename(sys.argv[0])
    folderlist = ("./logs/",)
    if add_folder is not None:
        folderlist = (*folderlist, add_folder)
    for folder in folderlist:
        with open(folder + "/" + file_name + ".log", 'a') as logfile:
            logfile.write(dt_string + "\n" + " ".join(sys.argv) + '\n')
    return


def print_script_call():
    """ save full script call in logfile """
    import datetime
    now = datetime.datetime.now()
    dt_string = now.strftime("%Y/%m/%d %H:%M:%S")
    print(dt_string + "\n" + " ".join(sys.argv) + '\n')
    return


# remove everything to the left of (and including) the first delimiter:
def remove_left_of_first(delimiter, label):
    return re.sub(r'(^.*?)' + delimiter, '', label)


# remove everything to the left of (and including) the last delimiter:
def remove_left_of_last(delimiter, label):
    return re.sub(r'(^.*)' + delimiter, '', label)


# remove everything to the right of (and including) the last delimiter:
# TODO FIX THIS ONE
def remove_right_of_last(delimiter, label):
    print("ERROR THIS FUNCTION DOESNT WORK")
    exit()
    return re.sub(r'(^.*?)' + delimiter + r'(.*?)', r'\1', label)


# remove everything to the right of (and including) the first delimiter:
def remove_right_of_first(delimiter, label):
    return re.sub(r'(^.*?)' + delimiter + r'(.*)', r'\1', label)


# === extract information from qcdtype and conftype ===


def parse_conftype(conftype):
    beta = remove_left_of_first('_b', conftype)
    beta = remove_right_of_first('_', beta)
    beta = float(beta) / 100000
    ns = remove_left_of_first('s', conftype)
    ns = int(remove_right_of_first('t', ns))
    nt = remove_left_of_first('t', conftype)
    nt = int(remove_right_of_first('_b', nt))
    nt_half = int(nt / 2)
    return beta, ns, nt, nt_half


def parse_qcdtype(qcdtype):
    fermions = remove_right_of_first('_', qcdtype)
    flowtype = remove_left_of_last('_', qcdtype)
    temp = remove_left_of_first('_', qcdtype)
    temp = remove_right_of_first('_', temp)
    return fermions, temp, flowtype


# === read and parse cmd line arguments ===


def get_parser():
    parser = argparse.ArgumentParser()
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('--qcdtype', help="format doesnt matter, only used for finding data. example: quenched_1.50Tc_zeuthenFlow", required=True)
    requiredNamed.add_argument('--corr', choices=['EE', 'EE_clover', 'BB', 'BB_clover'], help="choose from EE, EE_clover, BB, BB_clover", required=True)

    # optional
    # parser.add_argument('--conftype', help="format: s064t64_b0824900_m002022_m01011")

    # return qcdtype, fermions, temp, flowtype, corr, add_params
    # return onftype, beta, ns, nt, nt_half, qcdtype, fermions, temp, flowtype, corr, add_params
    return parser, requiredNamed


def get_merged_data_path(qcdtype, corr, conftype):
    return "../../data/merged/" + qcdtype + "/" + corr + "/" + conftype + "/"


def get_raw_data_path(qcdtype, conftype):
    return "../../data/raw/" + qcdtype + "/" + conftype + "/"


def get_plot_path(qcdtype, corr, conftype):
    return "../../plots/" + qcdtype + "/" + corr + "/" + conftype + "/"


def create_folder(*paths):
    from pathlib import Path
    for path in paths:
        Path(path).mkdir(parents=True, exist_ok=True)
    return


def write_flow_times(outfile, flow_times):
    outfile.write('# flow times: ')
    for i in flow_times:
        outfile.write(str(i) + " ")
    outfile.write('\n')


def chmap(mydict, **kwargs):
    c = ChainMap(kwargs, mydict)
    return c


# === functions related to the EE correlator ===

# improved tauTs from mathematica
# tauT_improved = {16:(6.883194e-02, 1.085070e-01, 1.616768e-01, 2.247303e-01, 2.914720e-01, 3.571978e-01, 4.178695e-01, 4.527297e-01),
# 20:(5.506555e-02, 8.680590e-02, 1.293465e-01, 1.798279e-01, 2.334154e-01, 2.867992e-01, 3.389318e-01, 3.894339e-01, 4.362603e-01, 4.632426e-01),
# 24:(4.588796e-02, 7.233829e-02, 1.077897e-01, 1.498636e-01, 1.945462e-01, 2.391203e-01, 2.828289e-01, 3.257494e-01, 3.679753e-01, 4.092103e-01, 4.476142e-01, 4.697792e-01),
# 28:(0.039332543108858704, 0.06200429713148807, 0.09239146163313365, 0.12845583966778362, 0.16676106589862263, 0.20498393407330573, 0.2424922712244153, 0.2793926631132959, 0.31587145378525494, 0.3520157792617203, 0.38777053007225404, 0.4227840462554289, 0.45543361809657334, 0.4742844242530861),
# 30:(3.671037e-02, 5.787068e-02, 8.623207e-02, 1.198924e-01, 1.556450e-01, 1.913227e-01, 2.263377e-01, 2.607949e-01, 2.948799e-01, 3.287019e-01, 3.622833e-01, 3.955382e-01, 4.281210e-01, 4.585112e-01, 4.760585e-01),
# 36:(0.030591978346158538, 0.048225570112454395, 0.07186010884057005, 0.0999106894769821, 0.12970558033931176, 0.15943980536845198, 0.1886256825235924, 0.2173546715249435, 0.24579006824642974, 0.2740405351256582, 0.302163599109069, 0.3301833888944985, 0.3580983972194321, 0.385875234263269, 0.41341447770969214, 0.44041363885624496, 0.46560294274945346, 0.480147821801901)}


def get_tauTs(nt):
    return numpy.arange(1 / nt, 0.5001, 1 / nt)  # default tauT
    # return tauT_improved[nt] #tree-level improved tauT for XX correlators


# only this called in the other scripts
# tauT_imp = {}
# for key in tauT_unimproved:
# tauT_imp[key] = [*tauT_improved[key][0:5], *tauT_unimproved[key][5:]]

def EE_cont_LO(tauT):
    return math.pi ** 2 * (math.cos(math.pi * tauT) ** 2 / math.sin(math.pi * tauT) ** 4 + 1 / (3 * math.sin(math.pi * tauT) ** 2))


EE_latt_LO_flow = {}
for Ntau in (16, 20, 24, 30, 32, 36, 44, 48, 56, 64):
    EE_latt_LO_flow[Ntau] = numpy.loadtxt("/home/altenkort/work/correlators_flow/data/merged/quenched_pertLO_wilsonFlow/EE/EE_latt_flow_" + str(Ntau) + '.dat')


def improve_corr_factor(tauTindex, nt, flowindex, improve=True, improve_with_flow=False):
    nt_half = int(nt / 2)
    if improve:
        if improve_with_flow:
            return nt ** 4 / EE_latt_LO_flow[nt][tauTindex + nt_half * flowindex][2]  # here, "G_norm" is pert latt corr at zero flow time
        else:
            return nt ** 4 / EE_latt_LO_flow[nt][tauTindex][2]  # here, "G_norm" is pert latt corr at zero flow time
    else:
        return nt ** 4 / EE_cont_LO(tauT_imp[nt][tauTindex])  # here, "G_norm" is pert cont corr at zero flow time


def lower_tauT_limit(flowradius, coarsest_Ntau=20):
    return flowradius/numpy.sqrt(8*0.014) + 1/coarsest_Ntau


def upper_flow_limit(tauT, coarsest_Ntau=20):
    return (tauT - 1 / coarsest_Ntau) * numpy.sqrt(8 * 0.014)


# === some data that does not really belong to extra files ===
ToverTc = {16: 1.510424, 20: 1.473469, 24: 1.484769, 30: 1.511761, 36: 1.504171}


# === plotting ===

def get_color(myarray, i, start=0, end=-1):
    if myarray[end] == myarray[start]:
        return matplotlib.cm.gnuplot(0)
    # i = end - 1 - i + start
    return matplotlib.cm.gnuplot((myarray[i] - myarray[start]) / (myarray[end] - myarray[start]) * 0.9)


def get_color2(myarray, i, start=0, end=-1):
    if myarray[end] == myarray[start]:
        return matplotlib.cm.viridis(0)
    # i = end - 1 - i + start
    return matplotlib.cm.viridis((myarray[i] - myarray[start]) / (myarray[end] - myarray[start]) * 0.9)



# styles
fontsize = 10
mylinewidth = 0.5
figurestyle = dict(bbox_inches="tight", pad_inches=0)
plotstyle_points = dict(fmt='D-', linewidth=1, markersize=4, capsize=2, mew=0.5, fillstyle='none')
plotstyle_lines = dict(fmt='-', linewidth=0.5, mew=0.25, mec='grey', markersize=2, capsize=3)
plotstyle_add_point = dict(fmt='D-', fillstyle='none', markersize=5, mew=0.25, lw=0.5, elinewidth=0.25, capsize=1.2)
plotstyle_add_point_single = plotstyle_add_point
plotstyle_add_point_single.update(dict(fmt='D'))
labelboxstyle = dict(boxstyle="Round", fc="w", ec="None", alpha=0.9, pad=0.1, zorder=99)
legendstyle = dict(loc="center left", bbox_to_anchor=(1, 0.5), frameon=True, framealpha=0.8, edgecolor='none', fancybox=False, facecolor="w", columnspacing=0.1,
                   labelspacing=0.1, borderpad=0.1, handletextpad=0.4, handlelength=1,
                   handler_map={matplotlib.container.ErrorbarContainer: matplotlib.legend_handler.HandlerErrorbar(xerr_size=1, yerr_size=0.3)})
plotstyle_fill = dict(linewidth=0.5)
xlabelstyle = dict(bbox=labelboxstyle, zorder=-1)  # horizontalalignment='right', verticalalignment='bottom', bbox=labelboxstyle)
ylabelstyle = dict(bbox=labelboxstyle, horizontalalignment='left', verticalalignment='top', rotation=0,
                   zorder=-1)  # horizontalalignment='right', , rotation=0, bbox=labelboxstyle)
titlestyle = dict(x=0.5, y=0.9, bbox=labelboxstyle, verticalalignment='top', zorder=99, fontsize=10)
verticallinestyle = dict(ymin=0, ymax=1, color='grey', alpha=0.8, zorder=-10000, dashes=(4, 4), lw=0.5)
horizontallinestyle = dict(xmin=0, xmax=1, color='grey', alpha=0.8, zorder=-10000, dashes=(4, 4), lw=0.5)
flowlimitplotstyle = dict(fmt='|', mew=0.7, markersize=5)

markers = ['.', '+', 'x', 'P', '*', 'X', 'o', 'v', 's', 'H', '8', 'd', 'p', '^', 'h', 'D', '<', '>']


def create_figure(xlims=None, ylims=None, xlabel="", ylabel="", xlabelpos=(0.99, 0.01), ylabelpos=(0.01, 0.97), tickpad=2,
                  figsize=(3+3/8, 3+3/8 - 1/2.54), UseTex=True):
    if UseTex:
        matplotlib.pyplot.rc('text', usetex=True)
        matplotlib.pyplot.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{mathtools}')
    matplotlib.pyplot.rc('font', family='serif', size=fontsize)
    linewidth = 0.5
    matplotlib.rcParams['axes.linewidth'] = linewidth
    fig = matplotlib.pyplot.figure(figsize=figsize, constrained_layout=True)
    ax = fig.add_subplot(1, 1, 1)
    if xlabelpos is not None:
        ax.xaxis.set_label_coords(*xlabelpos)
    if ylabelpos is not None:
        ax.yaxis.set_label_coords(*ylabelpos)
    ax.minorticks_off()
    ax.tick_params(direction='in', pad=tickpad, width=linewidth)
    if xlims is not None:
        ax.set_xlim(xlims)
    if ylims is not None:
        ax.set_ylim(ylims)
    ax.set_xlabel(xlabel, **xlabelstyle)
    ax.set_ylabel(ylabel, **ylabelstyle)
    plots = []
    matplotlib.rc('image', cmap='Set1')
    # fig.set_tight_layout(dict(pad=0.4))
    return fig, ax, plots


# === some parameters ===
n_skiprows = 15  # skip the discrete points at the start of the continuum limit file
offset = 0.02  # for flow limits
# index for tauT
start = 2
end = 14
# index for tau_F
flowstart = 10
flowend = 21
