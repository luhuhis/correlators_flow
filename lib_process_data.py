import numpy
import re
import math
import argparse
from collections import ChainMap
import matplotlib
from matplotlib import pyplot, container, legend_handler
import sys
import scipy.interpolate
import concurrent.futures


def format_float(number, digits=3):
    thisformat = '{0:.'+str(int(digits))+'f}'
    return thisformat.format(number)


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

    if fermions == "hisq":
        gaugeaction = "LW"
    elif fermions == "quenched":
        gaugeaction = "Wilson"
    else:
        gaugeaction = None
    if flowtype == "zeuthenFlow":
        flowaction = "Zeuthen"
    elif flowtype == "wilsonFlow":
        flowaction = "Wilson"
    else:
        flowaction = None

    return fermions, temp, flowtype, gaugeaction, flowaction


# === read and parse cmd line arguments ===
def get_parser():
    parser = argparse.ArgumentParser()
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('--qcdtype', help="format doesnt matter, only used for finding data. example: quenched_1.50Tc_zeuthenFlow", required=True)
    requiredNamed.add_argument('--corr', choices=['EE', 'EE_clover', 'BB', 'BB_clover'], help="choose from EE, EE_clover, BB, BB_clover", required=True)
    return parser, requiredNamed


def get_merged_data_path(qcdtype, corr, conftype, basepath="../../data/merged/"):
    return basepath + "/" + qcdtype + "/" + corr + "/" + conftype + "/"


def get_raw_data_path(qcdtype, conftype, basepath="../../data/raw/"):
    return basepath + qcdtype + "/" + conftype + "/"


def get_plot_path(qcdtype, corr, conftype, basepath="../../plots/"):
    return basepath + qcdtype + "/" + corr + "/" + conftype + "/"


def create_folder(*paths):
    from pathlib import Path
    for path in paths:
        if not Path(path).is_dir():
            print("Creating path ", path)
            Path(path).mkdir(parents=True)  # exist_ok=True
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


def EE_cont_LO(tauT, flowradius=0):
    EE_cont_zeroflow = math.pi ** 2 * (math.cos(math.pi * tauT) ** 2 / math.sin(math.pi * tauT) ** 4 + 1 / (3 * math.sin(math.pi * tauT) ** 2))
    if numpy.isscalar(flowradius) and flowradius == 0:
        return EE_cont_zeroflow
    else:
        flowradius = numpy.asarray(flowradius)
        tmpsum = 0
        with numpy.errstate(divide='ignore', invalid='ignore'):  # supress divide by 0 warnings.
            for m in range(-4, 5):
                tmpsum += 1/(tauT+m)**4 * ((tauT+m)**4/flowradius**4 + (tauT+m)**2/flowradius**2 + 1) * numpy.exp(-((tauT+m)**2/flowradius**2))
            answer = numpy.asarray(-(1/numpy.pi**2) * tmpsum + EE_cont_zeroflow)
        # fix flowradius=0 entries
        indices = numpy.where(flowradius == 0)
        answer[indices] = EE_cont_zeroflow
        return answer


G_latt_LO_flow_storage = {}


def G_latt_LO_flow(tau_index: int, flowtime, corr: str, Nt: int, flowaction: str, gaugeaction: str):
    """This returns the perturbative LO EE or BB correlator under flow for various gauge and/or flow actions.
    Data from numerical calculations is interpolated in flowtime to return any input flow time(s) (numpy arrays of flowtimes are supported).
    Once one dataset is loaded and interpolated it is stored in a global dictonary for future access.

    param tau_index
        Min = 0 (--> tau/a=1), Max = Nt/2-1 (--> tau/a=Nt/2).

    """

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
    this_int = G_latt_LO_flow_storage[identifier][tau_index]
    return this_int(flowtime)


def lower_tauT_limit_(flowradius, max_FlowradiusBytauT=numpy.sqrt(8*0.014), tauT_offset=0):  # note: sqrt(8*0.014) ~= 0.33
    return flowradius/max_FlowradiusBytauT + tauT_offset


def upper_flowradius_limit_(tauT, max_FlowradiusBytauT=numpy.sqrt(8*0.014), tauT_offset=0):  # note: sqrt(8*0.014) ~= 0.33
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


def get_discrete_color(i):
    i = i % 10
    return "C"+str(i)


def get_marker(i):
    mymarkers = ['o', 's', 'D', 'H', 'p']
    return mymarkers[i]


def lighten_color(color, amount=0.5):
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])


def leg_err_size(x=1, y=0.3):
    return dict(handler_map={matplotlib.container.ErrorbarContainer: matplotlib.legend_handler.HandlerErrorbar(xerr_size=x, yerr_size=y)})


labelboxstyle = dict(boxstyle="Round", fc="w", ec="None", alpha=0.9, pad=0.1, zorder=99)
verticallinestyle = dict(color='grey', alpha=0.8, zorder=-10000, dashes=(4, 4))  # ymin=0, ymax=1,
horizontallinestyle = dict(color='grey', alpha=0.8, zorder=-10000, dashes=(4, 4))  # xmin=0, xmax=1,z
markers = ['.', '+', 'x', 'P', '*', 'X', 'o', 'v', 's', 'H', '8', 'd', 'p', '^', 'h', 'D', '<', '>', '.', '+', 'x', 'P', '*', 'X', 'o', 'v', 's', 'H', '8',
           'd', 'p', '^', 'h', 'D', '<', '>', '.', '+', 'x', 'P', '*', 'X', 'o', 'v', 's', 'H', '8', 'd', 'p', '^', 'h', 'D', '<', '>', '.', '+', 'x', 'P',
           '*', 'X', 'o', 'v', 's', 'H', '8', 'd', 'p', '^', 'h', 'D', '<', '>']
linewidth = 1.5
axeslinewidth = 0.5
plot_fontsize = 11


def set_rc_params():
    matplotlib.pyplot.rc('font', family='cmr10', size=plot_fontsize)
    matplotlib.pyplot.rc('axes', unicode_minus=False)
    matplotlib.rcParams['mathtext.fontset'] = 'cm'
    matplotlib.rcParams['figure.constrained_layout.use'] = True
    matplotlib.rcParams['axes.titlesize'] = plot_fontsize
    matplotlib.rcParams['axes.linewidth'] = axeslinewidth
    matplotlib.rcParams['lines.linewidth'] = linewidth
    matplotlib.rcParams['lines.markeredgewidth'] = linewidth
    matplotlib.rcParams['errorbar.capsize'] = 2*linewidth
    matplotlib.rcParams['legend.handlelength'] = 0.2
    matplotlib.rcParams['legend.frameon'] = True
    matplotlib.rcParams['legend.framealpha'] = 0.8
    matplotlib.rcParams['legend.edgecolor'] = 'none'
    matplotlib.rcParams['legend.fancybox'] = False
    matplotlib.rcParams['legend.facecolor'] = 'w'
    matplotlib.rcParams['legend.labelspacing'] = 0.15
    matplotlib.rcParams['legend.columnspacing'] = 0.1
    matplotlib.rcParams['legend.borderpad'] = 0.1
    matplotlib.rcParams['legend.handletextpad'] = 0.4
    matplotlib.rcParams['legend.fontsize'] = plot_fontsize
    matplotlib.rcParams['legend.title_fontsize'] = plot_fontsize


def create_figure(xlims=None, ylims=None, xlabel="", ylabel="", xlabelpos=(0.98, 0.01), ylabelpos=(0.01, 0.98), tickpad=1,
                  figsize=None, UseTex=True, fig=None, subplot=111, no_ax=False, constrained_layout=True, xlabelbox=labelboxstyle, ylabelbox=labelboxstyle):

    set_rc_params()
    if UseTex:
        matplotlib.pyplot.rc('text', usetex=True)
        matplotlib.pyplot.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{mathtools}')

    if figsize == "fullwidth":
        figsize = (15, 7)
    elif figsize == "fullwidth_slim":
        figsize = (15, 5)
    elif figsize == "wide" or figsize == ['wide']:
        figsize = (10.875, 7)  # 3/4 of the width
    elif figsize is None:
        figsize = (7, 7)
    elif len(figsize) == 2:
        figsize = (float(figsize[0]), float(figsize[1]))

    if fig is None:
        fig = matplotlib.pyplot.figure(figsize=[val*1/2.54 for val in figsize], constrained_layout=constrained_layout)
        fig.set_constrained_layout_pads(w_pad=0.01*1/2.54, h_pad=0.01*1/2.54)  # set 0.1mm (instead of 0) padding to not clip any axes lines.
    if not no_ax:
        ax = fig.add_subplot(subplot)
        if xlabelpos is not None:
            ax.xaxis.set_label_coords(*xlabelpos)
        if ylabelpos is not None:
            ax.yaxis.set_label_coords(*ylabelpos)
        ax.minorticks_off()
        ax.tick_params(pad=tickpad, width=axeslinewidth, direction="out")
        if xlims is not None:
            ax.set_xlim(xlims)
        if ylims is not None:
            ax.set_ylim(ylims)
        ax.set_xlabel(xlabel, horizontalalignment='right', verticalalignment='bottom', bbox=xlabelbox)
        ax.set_ylabel(ylabel, horizontalalignment='left', verticalalignment='top', rotation=0, bbox=ylabelbox)
    else:
        ax = None
    plots = []
    matplotlib.rc('image', cmap='Set1')
    return fig, ax, plots


# class that is used in parallel_function_eval down below
class ComputationClass:
    def __init__(self, function, input_array, nproc, *add_param):
        self._nproc = nproc  # number of processes
        self._input_array = input_array
        self._function = function

        # additional arguments for actual_computation
        self._add_param = add_param

        # compute the result when class is initialized
        self._result = self.parallelization_wrapper()

    def parallelization_wrapper(self):
        results = []
        single_return_value = False
        with concurrent.futures.ProcessPoolExecutor(max_workers=self._nproc) as executor:
            for result in executor.map(self.pass_argument_wrapper, self._input_array):
                if type(result) is tuple:
                    results.append(list(result))
                else:
                    results.append(result)
                    single_return_value = True
        if single_return_value:
            return results
        else:
            return list(map(list, zip(*results)))  # "transpose" the list to allow for multiple return values like a normal function.

    def pass_argument_wrapper(self, single_input):
        return self._function(single_input, *self._add_param)

    def getResult(self):
        return self._result


# in parallel, compute function(input_array)
def parallel_function_eval(function, input_array, nproc, *add_param):
    computer = ComputationClass(function, input_array, nproc, *add_param)
    return computer.getResult()


def serial_function_eval(function, input_array, *add_param):
    results = []
    single_return_value = False
    for i in input_array:
        result = function(i, *add_param)
        if type(result) is tuple:
            results.append(list(result))
        else:
            results.append(result)
            single_return_value = True
    if single_return_value:
        return results
    else:
        return list(map(list, zip(*results)))


def dev_by_dist(data, axis=0, return_both_q=False, percentile=68):
    """Calculate the distance between the median and 68% quantiles. Returns the larger of the two distances. This
    method is used sometimes to estimate error, for example in the bootstrap."""
    data = numpy.asarray(data)
    median = numpy.nanmedian(data, axis)
    numb_data = data.shape[axis]
    idx_dn = max(int(numpy.floor((numb_data-1) / 2 - percentile/100 * (numb_data-1) / 2)), 0)
    idx_up = min(int(numpy.ceil((numb_data-1) / 2 + percentile/100 * (numb_data-1) / 2)), numb_data-1)
    sorted_data = numpy.sort(data - numpy.expand_dims(median, axis), axis=axis)
    q_l = numpy.take(sorted_data, idx_dn, axis)
    q_r = numpy.take(sorted_data, idx_up, axis)
    if return_both_q:
        return numpy.abs(q_l), numpy.abs(q_r)
    else:
        return numpy.max(numpy.stack((numpy.abs(q_l), numpy.abs(q_r)), axis=0), axis=0)


def print_var(prefix, var):
    # print(prefix, var)
    return var