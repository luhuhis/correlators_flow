#!/usr/bin/env python3

import lib_process_data as lpd
import numpy
import scipy.integrate
import scipy.interpolate
import scipy.optimize
import argparse
import itertools
import matplotlib
from matplotlib import gridspec
from spf_reconstruction.model_fitting.spf_reconstruct import SpfArgs, Gnorm, Integrand
from matplotlib.backends.backend_pdf import PdfPages



# TODO change phiUV to npy file

def set_title(args, ax):

    def get_scale_string(min_scale_str, run_scale_str):
        return ', $\\:\\:\\mu=\\sqrt{[' + min_scale_str + ']^2+' + run_scale_str + '^2}$'

    def get_model_str(model, PhiUV_label, min_scale_str, run_scale_str):
        if model == "max":
            model_str = r'Fit: $\:\rho(\omega)= \mathrm{max}[\kappa \omega /2T, \, c\: \phi_\mathrm{UV}^\mathrm{' + PhiUV_label + r'}(\mu)]$' + get_scale_string(
                min_scale_str, run_scale_str)

        elif model == "pnorm":
            model_str = r'Dashed lines: fit with model $\rho(\omega)= \sqrt{[\kappa \omega /2T\,]^2 + [c\: \phi_\mathrm{UV}^\mathrm{' + PhiUV_label + r'}(\mu)]^2}$' + get_scale_string(
                min_scale_str, run_scale_str)
        elif model == "plaw":
            model_str = r'Fit: $\:\rho(\omega)= \rho_{IR}(\kappa, \omega) + \rho_{line}(\omega) + \rho^\mathrm{' + PhiUV_label + r'}_{UV}(\omega, c, \mu)$' + get_scale_string(
                min_scale_str, run_scale_str)
        else:
            # TODO add line model etc
            model_str = ""
        return model_str

    # set title
    model_str = get_model_str(args.model, args.PhiUV_label, args.min_scale, args.run_scale)
    title = model_str+args.title
    if title != "":
        ax.text(0.5, 1, title, transform=ax.transAxes, ha='center', va='bottom')


def set_text(args, ax, ax_kappa):
    if args.plot_quenched_systematics:
        ax_kappa.text((3.40024856+1.4370101)/2/4, 0.85, "total \n systematic\n uncertainty \n from model \n choice at \n $a\\rightarrow 0$,\n$\\tau_F\\rightarrow 0$", transform=ax_kappa.transAxes, fontsize=6, color='k', alpha=1, ha='center', va='center', rotation='0')
        ax_kappa.fill_between([1.4370101, 3.40024856], [0, 0], [10, 10], zorder=-10000, color='k', alpha=0.15, lw=0)

    if args.custom_text is not None:
        x = float(args.custom_text[0])
        y = float(args.custom_text[1])
        text = args.custom_text[2]
        ha = args.custom_text[3]
        va = args.custom_text[4]
        transform = None
        if args.custom_text[5] == "rel":
            transform = ax.transAxes

        ax.text(x, y, text, ha=ha, va=va, transform=transform)
        # ax.text(0.35/0.5, 0.4, "Fit $\\tau T \\geq 0.35$", transform=ax.transAxes, fontsize=6,
        #         color='k', alpha=1, ha='center', va='center', rotation='0', zorder=10, bbox=dict(facecolor='white', edgecolor='white'))


def set_watermark(args, ax, ax_kappa):
    # set watermark
    if args.show_watermark and not args.no_kappa_plot:
        ax.text(0.67, 0.2, 'HotQCD preliminary', transform=ax.transAxes,
                fontsize=16, color='gray', alpha=0.5,
                ha='center', va='center', rotation='20', zorder=-1000000)

        ax_kappa.text(0.51, 0.2, 'HotQCD preliminary', transform=ax_kappa.transAxes,
                      fontsize=8, color='gray', alpha=0.5,
                      ha='center', va='center', rotation='20', zorder=-1000000)


def set_legend(args, ax):

    # Legend
    for i in range(args.leg_n_dummies):
        ax.errorbar(0, 0, label=' ', markersize=0, alpha=0, lw=0)

    ax.legend()
    handles, labels = ax.get_legend_handles_labels()
    n = len(handles)
    n2 = int(n/2)

    if args.leg_interleave:
        index = numpy.asarray([(i, i+n2) for i in range(0, n2)], dtype=int).flatten()
    else:
        index = range(n)

    legend = ax.legend([handles[i] for i in index], [labels[i] for i in index],
                       loc=args.leg_loc, bbox_to_anchor=args.leg_pos, ncol=args.leg_n_col, **lpd.leg_err_size(), columnspacing=0.5, framealpha=args.leg_framealpha)

    legend.set_title(args.leg_title)
    # texts = legend.get_texts()
    # for text in texts:
    #     text.set_linespacing(0.85)
    # legend._handler_map
    # legend._legend_box.sep = 5  # separation between legend title and entries
    # legend._legend_box.align = "left"


def plot_error_interpolation(xdata, ydata, edata, y_int, e_int, tauT):

    x_int = numpy.linspace(xdata[0], xdata[-1], 1000)

    fig, ax, _ = lpd.create_figure()
    ax.set_xlim(0, tauT/2)
    ax.set_ylim(0, 12)
    ax.errorbar(xdata, ydata, edata, color='k', zorder=10)
    ax.fill_between(x_int, y_int(x_int)-e_int(x_int), y_int(x_int)+e_int(x_int), zorder=9)

    ax.text(0.99, 0.99, lpd.format_float(tauT), ha='right', va='top', transform=ax.transAxes)

    return fig


def plot_relative_flow(args, ax, plot_offset, spfargs_arr, params_ansatz_arr):
    # TODO put this part under correlator_analysis, and just plot the corresponding corrs here?
    # TODO add option to only calculate these corrs, don't plot anything.
    # find correlator at various fixed flowradiusBytauT
    if args.qcdtype is not None and args.corr is not None and args.conftype is not None:
        n = max(len(args.flowradiusBytauT), len(args.qcdtype), len(args.corr), len(args.conftype))
        args.color_data = args.color_data[plot_offset:plot_offset + n]
        args.markers = args.markers[plot_offset:plot_offset + n]
        args.fillstyle = args.fillstyle[plot_offset:plot_offset + n]
        args.leg_labels = args.leg_labels[plot_offset:plot_offset + n]
        for flowradiusBytauT, qcdtype, corr, conftype, color, marker, fillstyle, label_suffix, spfargs, params, a in \
                zip(args.flowradiusBytauT if len(args.flowradiusBytauT) > 1 else itertools.cycle(args.flowradiusBytauT),
                    args.qcdtype if len(args.qcdtype) > 1 else itertools.cycle(args.qcdtype),
                    args.corr if len(args.corr) > 1 else itertools.cycle(args.corr),
                    args.conftype if len(args.conftype) > 1 else itertools.cycle(args.conftype),
                    args.color_data,
                    args.markers if len(args.markers) > 1 else itertools.cycle(args.markers),
                    args.fillstyle if len(args.fillstyle) > 1 else itertools.cycle(args.fillstyle),
                    args.leg_labels if args.leg_labels is not None else ["" for _ in range(n)],
                    spfargs_arr if len(spfargs_arr) > 0 else [None for _ in range(n)],
                    params_ansatz_arr if len(params_ansatz_arr) > 0 else [None for _ in range(n)],
                    args.plot_in_fm if args.plot_in_fm[0] is not None else [None for _ in range(n)]):
            if flowradiusBytauT is not None:
                inputfolder = lpd.get_merged_data_path(qcdtype, corr, conftype, basepath="../../../../data/merged/")
                these_flowtimes = numpy.loadtxt(inputfolder + "flowtimes_" + conftype + ".dat")
                XX = numpy.loadtxt(inputfolder + corr + "_" + conftype + ".dat")
                XX_err = numpy.loadtxt(inputfolder + corr + "_err_" + conftype + ".dat")
                beta, ns, nt, nt_half = lpd.parse_conftype(conftype)
                these_flowradii = numpy.sqrt(8 * these_flowtimes) / nt
                these_tauT = numpy.arange(1 / nt, 0.501, 1 / nt)

                fermions, temp, flowtype, gaugeaction, flowaction = lpd.parse_qcdtype(qcdtype)

                # interpolate between flow times
                y_int = []
                e_int = []
                # figs = []
                for i in range(len(these_tauT)):
                    G_latt_LO_flow = lpd.G_latt_LO_flow(i, these_flowtimes, corr, nt, flowaction, gaugeaction)
                    ydata = numpy.asarray(XX[:, i]) * nt ** 4 / G_latt_LO_flow
                    edata = numpy.asarray(XX_err[:, i]) * nt ** 4 / numpy.fabs(G_latt_LO_flow)
                    xdata = these_flowradii
                    min_index = numpy.fabs(these_flowradii-args.min_flowradius).argmin()
                    y_int.append(scipy.interpolate.InterpolatedUnivariateSpline(xdata[min_index:], ydata[min_index:], k=3, ext=2, check_finite=True))
                    e_int.append(scipy.interpolate.InterpolatedUnivariateSpline(xdata[min_index:], edata[min_index:], k=3, ext=2, check_finite=True))

                    # figs.append(plot_error_interpolation(xdata[min_index:], ydata[min_index:], edata[min_index:], y_int[i], e_int[i], (i+1)/nt))

                # save figures to pdf
                # print("save figures...")
                # lpd.set_rc_params()
                # with PdfPages(args.output_path + "EE_" + conftype + "_relflow_int.pdf") as pdf:
                #     for fig in figs:
                #         pdf.savefig(fig)
                #         matplotlib.pyplot.close(fig)

                find_and_plot_and_save_relflow(args, ax, flowradiusBytauT, these_tauT, y_int, e_int,
                                               inputfolder, color, marker, fillstyle, spfargs, params,
                                               int(-1000 * flowradiusBytauT),
                                               nt if args.plot_in_lattice_units else None, a, label_suffix)


def find_and_plot_and_save_relflow(args, ax, flowradiusBytauT, possible_tauTs, y_int_list, e_int_list, out_file_path, color, marker, fillstyle,
                                   spfargs, params, zorder=-10, nt=None, lattice_spacing=None, label=""):

    relflow_edata = []
    relflow_ydata = []
    relflow_xdata = []
    for i, tauT in enumerate(possible_tauTs):
        wanted_flowradius = flowradiusBytauT * tauT
        if y_int_list[i] is not None:
            try:
                relflow_edata.append(e_int_list[i](wanted_flowradius))
                relflow_ydata.append(y_int_list[i](wanted_flowradius))
                relflow_xdata.append(tauT)
            except ValueError:
                print("Warn: ValueError for flow radius", wanted_flowradius)
                pass
    relflow_xdata = numpy.asarray(relflow_xdata)
    relflow_ydata = numpy.asarray(relflow_ydata)
    relflow_edata = numpy.asarray(relflow_edata)
    flowstr = r'{0:.2f}'.format(flowradiusBytauT)
    lpd.create_folder(out_file_path + "/rel_flow/")
    numpy.savetxt(out_file_path + "/rel_flow/"+args.corr[0]+"_relflow_" + flowstr + ".dat",
                  numpy.column_stack((relflow_xdata, relflow_ydata, relflow_edata)),
                  header="tauT      G/Gnorm_sqrt(8tauF)/tau=" + flowstr + "       err")

    if args.norm_by_Gansatz:
        relflow_ydata = [val*Gnorm(relflow_xdata[i])/G_ansatz(spfargs, params, relflow_xdata[i]) for i, val in enumerate(relflow_ydata)]
        relflow_edata = [val*Gnorm(relflow_xdata[i])/G_ansatz(spfargs, params, relflow_xdata[i]) for i, val in enumerate(relflow_edata)]

    if lattice_spacing is None:
        if nt is None:
            xdata = relflow_xdata
        else:
            xdata = numpy.asarray(relflow_xdata) * nt
    else:
        xdata = numpy.asarray(relflow_xdata) * nt * lattice_spacing
    ax.errorbar(xdata, relflow_ydata, relflow_edata, zorder=zorder, fmt=marker, color=color, fillstyle=fillstyle,
                label=label, markeredgewidth=0.75, elinewidth=0.75, markersize=4)  # TODO remove
    # TODO add custom labels
    if not args.no_connection:
        ax.errorbar(xdata, relflow_ydata, zorder=-100 * zorder,  markersize=0, alpha=0.5, fmt='-', color=color)


def plot_extrapolated_data(args, ax, plot_offset):
    if args.plot_flow_extr is not None:
        for i, path in enumerate(args.plot_flow_extr):
            try:
                XX_flow_extr = numpy.loadtxt(path, unpack=True)
                ax.errorbar(XX_flow_extr[0]*args.flow_extr_custom_units[i], XX_flow_extr[1], XX_flow_extr[2],
                            zorder=-10000, color=args.color_data[i], fillstyle=args.fillstyle[i], fmt=args.markers[i], label=args.leg_labels[i], markersize=4, markeredgewidth=0.75, elinewidth=0.75)
                if not args.no_connection:
                    ax.errorbar(XX_flow_extr[0], XX_flow_extr[1], zorder=-100000, markersize=0, color=args.color_data[i], fmt='-', alpha=0.5)
                plot_offset += 1
            except OSError:
                print("WARN: skipping tf->0 extr., could not find ", path)
    return plot_offset


# TODO fix this
def plot_relflow_continuum_data(args, ax, plot_offset):

    for datapath, flowpath in zip(args.plot_cont_extr, args.cont_flowtimes):
        plot_offset += 1
        conftype = args.conftype[0]  # TODO fix conftype

        # TODO update this for new sample-based format
        # load continuum quenched data and interpolate it
        EE_cont_arr = []
        flowtimes = numpy.loadtxt(flowpath)
        flowradii = numpy.sqrt(8*flowtimes)
        min_idx = 0
        max_idx = 0
        stop = False

        # TODO replace this loop over flowradius with one line that loads the npy file containing all flow and all samples.
        # TODO then compute the mean+std and use that further down the line
        for flowradius in flowradii:
            flowradius_str = r'{0:.4f}'.format(flowradius)
            fitparam_file = path+"/cont_extr/EE_" + flowradius_str + "_cont.txt"
            try:
                EE_cont_arr.append(numpy.loadtxt(fitparam_file, unpack=True))
                stop = True
                max_idx += 1
            except OSError:
                # print("could not find ", fitparam_file)
                if not stop:
                    min_idx += 1
        max_idx += min_idx
        flowradii = flowradii[min_idx:max_idx]  # filter which flowradii got actually read in
        possible_tauTs = EE_cont_arr[0][0]  # doesn't matter which one we choose here, they are all the same

        # interpolate between flowtimes for each tauT. for that we first need to remove the nan's
        cont_int = []
        cont_err_int = []
        for i in range(len(possible_tauTs)):
            ydata = numpy.asarray([EE_cont[1][i] for EE_cont in EE_cont_arr])
            edata = numpy.asarray([EE_cont[2][i] for EE_cont in EE_cont_arr])
            mask = ~numpy.isnan(ydata)
            xdata = flowradii[mask]
            ydata = ydata[mask]
            edata = edata[mask]
            if len(xdata) >= 4:
                cont_int.append(scipy.interpolate.InterpolatedUnivariateSpline(xdata, ydata, k=3, ext=2, check_finite=True))
                cont_err_int.append(scipy.interpolate.InterpolatedUnivariateSpline(xdata, edata, k=3, ext=2, check_finite=True))
            else:
                cont_int.append(None)
                cont_err_int.append(None)
        outfolder = path+"/cont_rel_flow/"
        lpd.create_folder(outfolder)

        find_and_plot_and_save_relflow(args, ax, args.flowradiusBytauT[0], possible_tauTs, cont_int, cont_err_int, outfolder, args.color_data[1], r'a\rightarrow 0',
                                       None, None, zorder=-9999)

    return plot_offset


def plot_fit_corrs(args, ax, xpoints_arr, model_means, just_UVs, fitparams, nts):
    zorders = range(len(xpoints_arr))
    lattice_spacings = args.plot_in_fm

    # plot reconstructed correlators
    for xpoints, model_mean, just_UV, fitparam, color, zorder, nt, a in itertools.zip_longest(xpoints_arr, model_means, just_UVs, fitparams, args.color_fit,
                                                                                              zorders, nts, lattice_spacings):
        if xpoints is not None:
            if args.plot_in_lattice_units:
                xpoints = xpoints * nt
            if a is not None:
                xpoints = xpoints * nt * a

            if not args.no_label:
                label = "$ \\kappa/T^3=" + "{0:.2f}".format(fitparam[0][0]) + "\\pm " + "{0:.2f}".format(numpy.fmax(fitparam[0][1], fitparam[0][2])) \
                      + "$\n $ c=" + "{0:.2f}".format(fitparam[1][0]) + "\\pm " + "{0:.2f}".format(numpy.fmax(fitparam[1][1], fitparam[1][2])) + "$"
                # "$ G_\\mathrm{fit}^{\\rho =\\mathrm{max(IR,c\\,UV_{a})}}$ \n "
                # + "$\\chi^2/\\mathrm{dof}=" + "{0:.1f}".format(fitparam[2][0]) + " \\pm " + "{0:.1f}".format(numpy.fmax(fitparam[2][1], fitparam[2][2])) +"$",
            else:
                label = ""
            ax.errorbar(xpoints, model_mean, color=color, label=label, fmt='--', zorder=-10000+zorder)

            if not args.no_just_UV:
                ax.errorbar(xpoints, just_UV, color=color, alpha=0.25, fmt='-.', zorder=-10000+zorder)


def G_ansatz(spfargs, params, tauT):
    return scipy.integrate.quad(lambda OmegaByT: Integrand(OmegaByT, tauT, spfargs, params), spfargs.MinOmegaByT, spfargs.MaxOmegaByT)[0]


def plot_fits(args, ax, ax_kappa):

    def get_nts(conftypes):
        nts = []
        if conftypes is not None:
            for conftype in conftypes:
                _, _, nt, _ = lpd.parse_conftype(conftype)
                nts.append(nt)
        return nts

    nts = get_nts(args.conftype)

    xpoints_arr = []
    Gmodel_by_norm_arr = []
    fitparams = []
    just_UVs = []
    spfargs_arr = []
    params_ansatz_arr = []

    if args.kappa_ypos and args.fitparam_files:
        pos = numpy.asarray(list(args.kappa_ypos))
        pos = pos[:len(args.fitparam_files)]

    posidx = 0
    # calculate fit corrs for all given fitparam and PhiUV files, and plot kappas
    if args.fitparam_files and args.model:
        if args.norm_by_Gansatz:
            # fit kappa
            xvalues = []
            yvalues = []
            evalues = []

            for fitparam_file, temp in \
                    zip(args.fitparam_files,
                        args.kappa_temps if args.kappa_temps is not None else [None for _ in range(len(args.fitparam_files))]):

                if args.deduce_fitparam_files:
                    fitparam_file = args.fitparam_basepath + "/" + fitparam_file + "/params.dat"
                else:
                    fitparam_file = args.fitparam_basepath + "/" + fitparam_file

                params, err_left, err_right = numpy.loadtxt(fitparam_file, unpack=True)
                factor = (temp / 1000) ** 3
                if not numpy.isnan(params[0]):
                    xvalues.append(temp)
                    yvalues.append(params[0] * factor)
                    evalues.append(numpy.fmax(err_left[0], err_right[0]) * factor)

            result = scipy.optimize.curve_fit(lambda _x, _m, _b: _m * _x + _b, xdata=xvalues, ydata=yvalues, sigma=evalues, p0=(0.001, -0.1))[0]
            m = result[0]
            b = result[1]

            xline = numpy.linspace(0, 500, 5)
            ax_kappa.errorbar(xline, m * xline + b, fmt=':', lw=0.5, color='gray', alpha=0.5)

        for fitparam_file, PhiUV_file, temp in \
                zip(args.fitparam_files,
                    args.PhiUV_files if len(args.PhiUV_files) > 1 else itertools.cycle(args.PhiUV_files),
                    args.kappa_temps if args.kappa_temps is not None else [None for _ in range(len(args.fitparam_files))]):

            # kappa plot
            if args.deduce_fitparam_files:
                fitparam_file = args.fitparam_basepath + "/" + fitparam_file + "/params.dat"
            else:
                fitparam_file = args.fitparam_basepath + "/" + fitparam_file

            # find out lower tauT based on file path
            # tmp = lpd.remove_left_of_last('/', fitparam_file)
            # tmp = lpd.remove_right_of_last('.dat', tmp)
            # tmp = "0" + lpd.remove_left_of_first('_0', tmp)
            # this_mintauT = float(lpd.remove_right_of_first('_', tmp))
            # TODO, if min tauT is "detected" and it is the same as other min tauT, then display it in the title!?
            # TODO, also do this for other properties of the SPF fit!!! look into plot_fit.py or so!

            this_mintauT = 0.05

            print("read in ", fitparam_file)
            params, err_left, err_right = numpy.loadtxt(fitparam_file, unpack=True)
            fitparams.append(numpy.column_stack((params, err_left, err_right)))

            if args.norm_by_Gansatz:
                params_ansatz = numpy.copy(params)
                params_ansatz[0] = (m * temp + b) / (temp / 1000) ** 3
                params_ansatz_arr.append(params_ansatz)

            if not args.no_kappa_plot:
                if not args.kappa_xaxis_temp:
                    ax_kappa.errorbar(params[0], pos[posidx], xerr=[[err_left[0]], [err_right[0]]], color=args.color_fit[posidx], fmt='x-', fillstyle='none',
                                      markersize=0, zorder=-10)
                else:
                    factor = 1 if not args.kappa_in_GeV else (temp / 1000) ** 3
                    ax_kappa.errorbar(temp, params[0] * factor, yerr=[[err_left[0] * factor], [err_right[0] * factor]], color=args.color_fit[posidx], fmt='x-',
                                      fillstyle='none',
                                      markersize=0, zorder=-10)
                    # print(temp, r'{0:.2f}'.format(params[0]), r'{0:.2f}'.format(err_left[0]), r'{0:.2f}'.format(err_right[0]))

                # if posidx == 0 and args.plot_extr:
                # ax_kappa.axvline(params[0], **lpd.chmap(lpd.verticallinestyle, dashes=(2, 2), color='k', alpha=0.5))

            posidx += 1

            # load PhiUV and calculate model correlators
            PhiUV_file = args.PhiUV_basepath + PhiUV_file
            print("read in ", PhiUV_file)
            PhiUV = numpy.loadtxt(PhiUV_file, unpack=True)[0:3]
            PhiUVByT3 = scipy.interpolate.InterpolatedUnivariateSpline(PhiUV[0], PhiUV[1], k=3, ext=2, check_finite=True)
            MinOmegaByT = PhiUV[0][0]
            MaxOmegaByT = PhiUV[0][-1]

            Gmodel_by_norm = []
            just_UV = []

            spfargs = SpfArgs(args.model, "", False, PhiUVByT3, -1, args.OmegaByT_IR, args.OmegaByT_UV, args.p, MinOmegaByT, MaxOmegaByT)
            spfargs_arr.append(spfargs)

            xpoints_arr.append(numpy.linspace(this_mintauT, 0.5, args.npoints))
            for tauT in xpoints_arr[-1]:

                Gmodel = scipy.integrate.quad(lambda OmegaByT: Integrand(OmegaByT, tauT, spfargs, params), MinOmegaByT, MaxOmegaByT)[0]

                normalization = Gnorm(tauT)
                if args.norm_by_Gansatz:
                    normalization = G_ansatz(spfargs, params_ansatz, tauT)

                Gmodel_by_norm.append(Gmodel / normalization)

                UV_params = numpy.copy(params)
                UV_params[0] = 0  # set kappa to zero
                GmodelUV = scipy.integrate.quad(lambda OmegaByT: Integrand(OmegaByT, tauT, spfargs, UV_params), MinOmegaByT, MaxOmegaByT)[0]
                just_UV.append(GmodelUV / normalization)

            Gmodel_by_norm_arr.append(Gmodel_by_norm)
            just_UVs.append(just_UV)

        # kappa plot fit range
        # TODO add check in beginning for error when mintauT, kappa_pos and fitparam files have different number of arguments
        if not args.no_kappa_plot and not args.kappa_xaxis_temp:
            labels_plot = []
            for mintauT in args.min_tauT:
                labels_plot.append("$" + str(mintauT) + "$")
                ax_kappa.text(0.02, args.kappa_labelheight, "Fit $\\tau T \\geq$", transform=ax_kappa.transAxes, fontsize=6, color='k', alpha=1, ha='left',
                              va='center', rotation='0')
                for label in ax_kappa.yaxis.get_majorticklabels():
                    label.set_ha("left")
                ax_kappa.tick_params(axis='y', length=0, direction='in', labelsize=7, pad=-2)
                from matplotlib.ticker import MaxNLocator
                ax_kappa.xaxis.set_major_locator(MaxNLocator(4))
            ax_kappa.set_yticks(pos)
            ax_kappa.set_yticklabels(labels_plot)

        # kappa plot temp
        if not args.no_kappa_plot and args.kappa_xaxis_temp:
            ax_kappa.set_xlim((0, 400))
            ax_kappa.set_ylim(args.kappa_xlims)
            ax_kappa.set_xlabel("$T\\, \\mathrm{[MeV]}$") #fontsize=fontsize
            if args.kappa_in_GeV:
                ax_kappa.set_ylabel("$\\kappa \\, [\\mathrm{GeV}^3]$")  #, fontsize=fontsize
            else:
                ax_kappa.set_ylabel("$\\kappa /T^3$")  #, fontsize=fontsize # [\\mathrm{GeV}^3]
            ax_kappa.xaxis.set_label_coords(0.8, 0.06)
            ax_kappa.yaxis.set_label_coords(0.15, 0.99)  # 0.07
            xpoints = numpy.linspace(0, 500, 100)
            if args.kappa_in_GeV:
                ax_kappa.errorbar(xpoints, 4 * numpy.pi / 0.9 * (xpoints / 1000) ** 3, fmt='--', color='grey', alpha=0.8)
                ax_kappa.text(140, 0.25, 'AdS/CFT estimate', fontsize=5, color='grey', alpha=0.8, ha='left')
                ax_kappa.text(270, 0.20, 'Linear fit', fontsize=5, color='grey', alpha=0.8)
            else:
                ax_kappa.axhline(4 * numpy.pi / 0.9, **lpd.horizontallinestyle)  # 4*numpy.pi/0.9 ~= 14
                ax_kappa.text(390, 4 * numpy.pi / 0.9 * 0.99, 'AdS/CFT estimate', fontsize=5, color='grey', alpha=0.8, ha='right', va='top')

        plot_fit_corrs(args, ax, xpoints_arr, Gmodel_by_norm_arr, just_UVs, fitparams, nts)

    return spfargs_arr, params_ansatz_arr


def prepare_plot_canvas(args):

    # set xlabel
    if args.plot_in_lattice_units:
        xlabel = r'$\tau/a$'
    elif args.plot_in_fm[0] is not None:
        xlabel = r'$\tau \mathrm{[fm]}$'
    else:
        xlabel = r'$\tau T$'

    if args.norm_by_Gansatz:
        ylabel = "$\\frac{G^\\mathrm{imp}}{G^\\mathrm{ansatz}}$"
    else:
        ylabel = "$\\dfrac{G}{G^\\mathrm{norm}}$"

    ax_kappa = None
    if not args.no_kappa_plot:
        fig, _, _ = lpd.create_figure(figsize=args.figsize, UseTex=args.usetex, subplot=None, no_ax=True)
        spec = gridspec.GridSpec(figure=fig, ncols=2, nrows=1, width_ratios=[3, 1], wspace=0)
        _, ax, _ = lpd.create_figure(xlims=args.xlims, ylims=args.ylims, xlabel=xlabel, ylabel=ylabel,
                                     UseTex=args.usetex, fig=fig, subplot=spec[0])
        ax.set_xlabel(xlabel)
        _, ax_kappa, _ = lpd.create_figure(ylims=(0, 2.5), UseTex=args.usetex, fig=fig, subplot=spec[1])
        kappa_ylabel = "$\\kappa /T^3$"
        ax_kappa.set_ylabel(kappa_ylabel, x=0.7, y=0.1)
        ax_kappa.set_xlim(args.kappa_xlims)
        if args.kappa_ylims:
            ax_kappa.set_ylim(args.kappa_ylims)

    else:
        fig, ax, _ = lpd.create_figure(figsize=args.figsize, UseTex=args.usetex, xlims=args.xlims, xlabel=xlabel, ylims=args.ylims, ylabel=ylabel)

    if args.norm_by_Gansatz:
        ax.axhline(y=1, **lpd.horizontallinestyle)

    # plot dashed vlines
    if not args.no_vlines and args.tauT_vlines is not None and args.color_vlines is not None:
        for tauT, color in zip(args.tauT_vlines, args.color_vlines):
            ax.axvline(tauT, **lpd.chmap(lpd.verticallinestyle, color=color, alpha=1))

    if args.xticks != "auto":
        ax.set_xticks([float(i) for i in args.xticks])

    return fig, ax, ax_kappa


def parse_args():

    def set_output_path(outputpath, qcdtype, corr):
        if not outputpath:
            if qcdtype:
                outputpath = lpd.get_plot_path(qcdtype[0], "", "") + "/rel_flow/"
            if qcdtype and corr:
                outputpath = lpd.get_plot_path(qcdtype[0], corr[0], "") + "/rel_flow/"
        lpd.create_folder(outputpath)
        return outputpath

    parser = argparse.ArgumentParser()

    parser.add_argument('--plot_flow_extr', help='paths to files containing flow-extr data.', default=None, type=str, nargs='*')
    parser.add_argument('--flow_extr_custom_units', nargs='*', type=float)
    parser.add_argument('--plot_cont_extr', help='paths to files containing continuum-extr data.', default=None, type=str, nargs='*')
    parser.add_argument('--cont_flowtimes', help='paths to files containing continuum flowtimes in temperature units', default=None, type=str, nargs='*')

    parser.add_argument('--plot_quenched_systematics', action="store_true")  # TODO soften hard code

    parser.add_argument('--flowradiusBytauT', nargs='*', type=float, help='which sqrt(8tauF)/tauT to plot')
    parser.add_argument('--min_flowradius', default=0.075, type=float, help="set to at least one lattice spacing of the coarsest lattice.")

    # fitparam options
    parser.add_argument('--fitparam_files', nargs='*', help='plot model 3 spf for these fit params')
    parser.add_argument('--fitparam_basepath', default="", type=str, help='shared basepath of all fitparam files, prepended to --fitparam_files.')
    parser.add_argument('--deduce_fitparam_files', help='with this option it suffices to only specify the fitparam file folder. it only works for a specific '
                                                        'folder structure/naming convention', action="store_true")

    # general plot options
    parser.add_argument('--usetex', action="store_true")
    parser.add_argument('--ylims', help='custom ylims', nargs=2, type=float, default=(-0.03, 4))
    parser.add_argument('--xlims', help='custom xlims', nargs=2, type=float, default=(-0.01, 0.55))
    parser.add_argument('--xticks', default="auto", nargs="*")
    parser.add_argument('--title', help='title prefix of plot', default="", type=str)
    parser.add_argument('--custom_text', help="show some custom text. modify this script to customize it. (i dont feel like adding 10 more parameters "
                                              "just for some text)", nargs='*')

    parser.add_argument('--plot_in_lattice_units', help="only works if the lattice spacing is constant across different conftypes", action="store_true")
    parser.add_argument('--plot_in_fm', nargs='*', type=float, default=[None, ], help="plot tau in fm. arguments here are the lattice spacings in order in fm.")
    parser.add_argument('--min_tauT', type=float, help='lower limits from which on the corr is fitted. appear in the kappa plot.', nargs='*')
    # parser.add_argument('--min_tauT_plot', type=float, help='lower limit from which tauT on the reconstructed corr is valid/should be plotted. '
    #                                                         'default is to parse it from fitparam file.', nargs='*', default=[0.05,])
    parser.add_argument('--tauT_vlines', type=float, help='where to plot vertical dashed lines', nargs='*')
    parser.add_argument('--npoints', help='how many tauT to plot between 0 and 0.5 for the integrated spfs (model correlators)', default=100, type=int)
    parser.add_argument('--norm_by_Gansatz', help="normalize by G_ansatz instead of Gnorm", action="store_true")

    # legend options
    parser.add_argument('--leg_pos', nargs=2, default=(0.02, 1), type=float, help="where to put the upper left corner of the legend")
    parser.add_argument('--leg_loc', type=str, default='lower right', help="anchor point of the legend")
    parser.add_argument('--leg_n_dummies', help="how many dummy entries to put in the legend in order to align things nicely", default=0, type=int)
    parser.add_argument('--leg_n_col', help="how many columns the legend should have", default=1, type=int)
    parser.add_argument('--leg_title', help="Put this at end of legend title", type=str, default="")
    parser.add_argument('--leg_labels', help="additional strings to put after the Ntau in the labels in legend", nargs='*', type=str)
    parser.add_argument('--leg_interleave', help="whether to interleave the given data in the legend. this is useful if you want to pair extrapolated data with nonextrapolated data.",
                        default=False, action="store_true")
    parser.add_argument('--leg_framealpha', default=0.8, type=int)

    # plot styling
    parser.add_argument('--figsize', help="size of the figure", default=None, nargs='*')
    parser.add_argument('--markers', nargs='*', default='|', type=str, help="fmt keywords for the plots")
    parser.add_argument('--fillstyle', nargs='*', default='none', type=str, help="fillstyle keywords for the plots")
    parser.add_argument('--color_data', help="specify colors in order that the data should have", type=str, nargs='*',
                        default=('C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9'))
    parser.add_argument('--color_fit', help="specify colors in order that the fits should have", type=str, nargs='*',
                        default=('C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9'))
    parser.add_argument('--color_vlines', help="specify colors in order that the vlines should have", type=str, nargs='*',
                        default=('C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9'))

    # hide certain plot elements
    parser.add_argument('--no_vlines', help="hide vertical lines that indicate lower fit limit", action="store_true")
    parser.add_argument('--no_connection', help="hide connections between data points", action="store_true")
    parser.add_argument('--no_label', help="hide labels of fit corrs", action="store_true")
    parser.add_argument('--no_just_UV', action="store_true", help="do not plot dash-dotted UV parts")

    # SPF options
    parser.add_argument('--model', type=str, help="which spf model to use. choices: see spf_reconstruct.py")
    # TODO get spf args from spf_reconstruct.py instead of repeating them here
    parser.add_argument('--OmegaByT_IR', type=float)
    parser.add_argument('--OmegaByT_UV', type=float)
    parser.add_argument('--p', type=float)
    parser.add_argument('--PhiUV_basepath', default="", type=str, help='shared basepath of all PhiUV files, prepended to --PhiUV_files.')
    parser.add_argument('--PhiUV_files', type=str, nargs='*', help='path to PhiUV file(s), space seperated. these should match the fitparam_files. '
                                                                   'if you specify a single PhiUV_file but multiple fitparam_files, than the single one is used '
                                                                   'for all of them. file should contain three columns: omega/T, PhiUVByT3, err')
    parser.add_argument('--PhiUV_label', type=str, help="will appear in the title of the plot. e.g. LO or NLO.", default="")
    parser.add_argument('--min_scale', default="", type=str, help='will appear in title of plot')
    parser.add_argument('--run_scale', default="", type=str, help='will appear in title of plot')

    # options for kappa plot (kappa on x-axis, fit range on y-axis)
    parser.add_argument('--no_kappa_plot', action="store_true", help="hide kappa plot")
    parser.add_argument('--show_watermark', action="store_true")
    parser.add_argument('--kappa_labelheight', help="how high the label that says Fit tauT is in the plot", default=0.75, type=float)
    parser.add_argument('--kappa_ypos', nargs='*', type=float, help="ypos of kappas in kappa plot, in order")
    parser.add_argument('--kappa_xlims', help="custom xlims for kappa plot", nargs=2, type=float, default=(0, 3))
    parser.add_argument('--kappa_ylims', help="custom xlims for kappa plot", nargs=2, type=float, default=(0, 3))

    # options for other version of kappa plot (kappa on y-axis, temp on x-axis)
    parser.add_argument('--kappa_xaxis_temp', help="plot xaxis:temp, yaxis:kappa in kappa plot, instead of xaxis:kappa yaxis:fitrange", action="store_true")
    parser.add_argument('--kappa_temps', help="temps of kappas of fit files in order", nargs='*', type=float)
    parser.add_argument('--kappa_in_GeV', help="use in together with --kappa_xaxis_temp. shows y-axis of kappa in GeV", action="store_true")

    # when you want to load finite lattice spacing
    parser.add_argument('--qcdtype', help="list of qcdtypes. format doesnt matter, only used for finding data. example: quenched_1.50Tc_zeuthenFlow", nargs='*',
                        type=str)
    parser.add_argument('--corr', help="list of corrs. choose from EE, EE_clover, BB, BB_clover", type=str, nargs='*')
    parser.add_argument('--conftype', nargs='*', type=str, help="list of conftypes, e.g. s080t20_b0703500 s096t24_b0719200 s120t30_b0739400 s144t36_b0754400")

    # file output options
    parser.add_argument('--output_path', help="where to save the plot.")
    parser.add_argument('--output_suffix', type=str, help="suffix to add to plot file output name")

    args = parser.parse_args()

    # check if arguments make sense
    if args.plot_cont_extr is not None:
        if not all(element == args.flowradiusBytauT[0] for element in args.flowradiusBytauT):
            # TODO remove this limitation
            print("ERROR: when plotting the quenched cont/flow extr, all flowradiusBytauT have to be identical")
            exit(1)
    if args.fitparam_files and not args.model:
        print("ERROR no model was specified although fit parameters were given")
        exit(1)
    if not args.no_kappa_plot and not args.kappa_ypos and not args.kappa_xaxis_temp:
        print("ERROR: need kappa_ypos")
        exit(1)

    if args.fitparam_files is not None and args.min_tauT is not None and (len(args.min_tauT) != len(args.fitparam_files)):
        print("ERROR need the same number of arguments for min_tauT as for fitparam_files")
        exit(1)
    if args.plot_in_lattice_units and args.plot_in_fm[0] is not None:
        print("ERROR choose either plot_in_lattice_units and plot_in_fm")
        exit(1)

    if args.conftype is not None and args.flowradiusBytauT is None:
        print("ERROR: need flowradiusBytauT")
        exit(1)
    args.output_path = set_output_path(args.output_path, args.qcdtype, args.corr)

    return args


def main():

    # TODO load flowtime instead of flowradius file

    args = parse_args()

    fig, ax, ax_kappa = prepare_plot_canvas(args)

    spfargs_array, params_ansatz_array = plot_fits(args, ax, ax_kappa)
    # we return some fit information, in case we want to normalize the relative flow correlators by some fit correlator.

    plot_offset = plot_extrapolated_data(args, ax, plot_offset=0)
    # we return an offset, so that plot_relative_flow knows that something else has already been plotted, so that it can skip already used colors.

    # plot_offset = plot_relflow_continuum_data(args, ax, plot_offset)

    plot_relative_flow(args, ax, plot_offset, spfargs_array, params_ansatz_array)

    set_legend(args, ax)
    set_watermark(args, ax, ax_kappa)
    set_text(args, ax, ax_kappa)
    set_title(args, ax)

    # save figure
    file = args.output_path + "/"+args.corr[0]+"_relflow" + args.output_suffix + ".pdf"
    print("saving ", file)
    fig.savefig(file)
    matplotlib.pyplot.close(fig)


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
