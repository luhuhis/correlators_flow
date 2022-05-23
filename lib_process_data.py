import numpy
import re
import math
import argparse
from collections import ChainMap
import matplotlib
from matplotlib import pyplot, cm, container, legend_handler
import sys
import scipy.interpolate


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
        create_folder(folder)
        with open(folder + "/" + file_name + ".log", 'a') as logfile:
            logfile.write(dt_string + "\n" + " ".join(sys.argv) + '\n')
    return


def print_script_call():
    """ print full script call to std out """
    import datetime
    now = datetime.datetime.now()
    dt_string = now.strftime("%Y/%m/%d %H:%M:%S")
    print("\n" + dt_string + "\n" + " ".join(sys.argv) + '\n')
    return


def remove_left_of_first(delimiter, label):
    """ remove everything to the left of (and including) the first delimiter """
    return re.sub(r'(^.*?)' + delimiter, '', label)


def remove_left_of_last(delimiter, label):
    """remove everything to the left of (and including) the last delimiter"""
    return re.sub(r'(^.*)' + delimiter, '', label)


def remove_right_of_last(delimiter, label):
    """remove everything to the right of (and including) the last delimiter"""
    return re.sub(r'(^.*)' + delimiter + r'(.*)', r'\1', label)


def remove_right_of_first(delimiter, label):
    """remove everything to the right of (and including) the first delimiter"""
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

def get_tauTs(nt):
    return numpy.arange(1 / nt, 0.5001, 1 / nt)  # default tauT
    # return tauT_improved[nt] #tree-level improved tauT for XX correlators


def EE_cont_LO(tauT):
    return math.pi ** 2 * (math.cos(math.pi * tauT) ** 2 / math.sin(math.pi * tauT) ** 4 + 1 / (3 * math.sin(math.pi * tauT) ** 2))


G_latt_LO_flow_storage = {}


def G_latt_LO_flow(tau: int, flowtime, corr: str, Nt: int, flowaction: str, gaugeaction: str):
    """This returns the perturbative LO EE or BB correlator under flow for various gauge and/or flow actions.
    Data from numerical calculations is interpolated in flowtime to return any input flow time(s) (numpy arrays of flowtimes are supported).
    Once one dataset is loaded and interpolated it is stored in a global dictonary for future access."""

    identifier = str(Nt)+flowaction+gaugeaction

    # store for future access
    global G_latt_LO_flow_storage
    if identifier not in G_latt_LO_flow_storage.keys():
        #TODO make these paths an input instead of hard code
        file = "/work/home/altenkort/work/correlators_flow/data/merged/pert_LO/"+corr+"_pert_latt_"+flowaction+"_"+gaugeaction+"_Nt"+str(Nt)+".dat"
        flowtimes = numpy.loadtxt("/work/home/altenkort/work/correlators_flow/data/merged/pert_LO/flowtimes.dat")
        try:
            tmp = numpy.loadtxt(file)
        except OSError:
            print("Error in latt_LO_flow: could not load file ", file)
            exit(1)
        interpolations = []
        for i in range(tmp.shape[1]):
            interpolations.append(scipy.interpolate.InterpolatedUnivariateSpline(flowtimes, tmp[:, i], k=3, ext=2))
        G_latt_LO_flow_storage[identifier] = interpolations

    # return the corresponding flowtime-interpolated correlator evaluated at the given flowtime(s)
    this_int = G_latt_LO_flow_storage[identifier][tau]
    return this_int(flowtime)


def lower_tauT_limit_(flowradius, max_FlowradiusBytauT=numpy.sqrt(8*0.014), tauT_offset=1/20):  # note: sqrt(8*0.014) ~= 0.33
    return flowradius/max_FlowradiusBytauT + tauT_offset


def upper_flowradius_limit_(tauT, max_FlowradiusBytauT=numpy.sqrt(8*0.014), tauT_offset=1/20):  # note: sqrt(8*0.014) ~= 0.33
    return (tauT - tauT_offset) * max_FlowradiusBytauT


# === some data that does not really belong to extra files ===
ToverTc = {16: 1.510424, 20: 1.473469, 24: 1.484769, 30: 1.511761, 36: 1.504171}


# === plotting ===

def get_color(myarray, i, start=0, end=-1, scale_factor=0.9):
    if myarray[end] == myarray[start]:
        return matplotlib.cm.gnuplot(0)
    # i = end - 1 - i + start
    return matplotlib.cm.gnuplot((myarray[i] - myarray[start]) / (myarray[end] - myarray[start]) * scale_factor)


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
xlabelstyle = dict(bbox=labelboxstyle, zorder=-1, fontsize=fontsize)  # horizontalalignment='right', verticalalignment='bottom', bbox=labelboxstyle)
ylabelstyle = dict(bbox=labelboxstyle, horizontalalignment='left', verticalalignment='top', rotation=0,
                   zorder=-1, fontsize=fontsize)  # horizontalalignment='right', , rotation=0, bbox=labelboxstyle)
titlestyle = dict(x=0.5, y=0.9, bbox=labelboxstyle, verticalalignment='top', zorder=99, fontsize=10)
verticallinestyle = dict(ymin=0, ymax=1, color='grey', alpha=0.8, zorder=-10000, dashes=(4, 4), lw=0.5)
horizontallinestyle = dict(xmin=0, xmax=1, color='grey', alpha=0.8, zorder=-10000, dashes=(4, 4), lw=0.5)
flowlimitplotstyle = dict(fmt='|', mew=0.7, markersize=5)

markers = ['.', '+', 'x', 'P', '*', 'X', 'o', 'v', 's', 'H', '8', 'd', 'p', '^', 'h', 'D', '<', '>', '.', '+', 'x', 'P', '*', 'X', 'o', 'v', 's', 'H', '8',
           'd', 'p', '^', 'h', 'D', '<', '>', '.', '+', 'x', 'P', '*', 'X', 'o', 'v', 's', 'H', '8', 'd', 'p', '^', 'h', 'D', '<', '>', '.', '+', 'x', 'P',
           '*', 'X', 'o', 'v', 's', 'H', '8', 'd', 'p', '^', 'h', 'D', '<', '>']


def create_figure(xlims=None, ylims=None, xlabel="", ylabel="", xlabelpos=(0.99, 0.03), ylabelpos=(0.01, 0.97), tickpad=2,
                  figsize=(3+3/8, 3+3/8 - 1/2.54), UseTex=True, fig=None, subplot=111, no_ax=False):
    if UseTex:
        matplotlib.pyplot.rc('text', usetex=True)
        matplotlib.pyplot.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{mathtools}')
    matplotlib.pyplot.rc('font', family='cmr10', size=fontsize)
    matplotlib.rcParams['mathtext.fontset'] = 'cm'
    linewidth = 0.5
    matplotlib.rcParams['axes.linewidth'] = linewidth
    if fig is None:
        fig = matplotlib.pyplot.figure(figsize=figsize, constrained_layout=True)
    if not no_ax:
        ax = fig.add_subplot(subplot)
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
    else:
        ax = None
    plots = []
    matplotlib.rc('image', cmap='Set1')
    # fig.set_tight_layout(dict(pad=0.4))
    return fig, ax, plots
