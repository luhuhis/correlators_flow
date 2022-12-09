#!/usr/bin/env python3
import numpy
from latqcdtools.statistics import bootstr
import lib_process_data as lpd

import process_data._4_continuum_extr as ce
import scipy


def identity(input):
    return input


def cont_extr(Ntaus, tauT_data, corr_data, corr_err_data, output_tauTs, nsamples):
    output_tauT_len = len(output_tauTs)

    # result variables
    cont_corr = numpy.empty(output_tauT_len)
    cont_corr[:] = numpy.nan
    cont_corr_err = numpy.empty(output_tauT_len)
    cont_corr_err[:] = numpy.nan

    # interpolate all coarser lattices to tauT of the finest lattice
    interpolations = numpy.ndarray((len(Ntaus), 2, output_tauT_len))
    interpolations[:] = numpy.nan
    for i, Ntau in enumerate(Ntaus):
        xdata = tauT_data[i]
        ydata = corr_data[i]
        edata = corr_err_data[i]

        # generate mock bootstrap samples
        samples, _, _ = bootstr.bootstr_from_gauss(identity, ydata, edata, nsamples, return_sample=True)

        # result variables
        theoutputdata = []
        # perform spline interpolation for each sample
        for m in range(nsamples):
            ydata = samples[m]
            # thefitparams[m] = _3_spline_interpolate.interpolate_EE_flow(xdata, ydata, None, theknots, order, constraints, output_xdata=None, return_params=True)
            spline = scipy.interpolate.CubicSpline(xdata, ydata, bc_type=((2, 0.0), (1, 0.0)))

            theoutputdata.append(spline(output_tauTs))

        interpolations[i][0] = numpy.mean(theoutputdata, axis=0)
        interpolations[i][1] = numpy.std(theoutputdata, axis=0)
        print("done interpolation Ntau", Ntau)

    # continuum extrapolation
    xdata = numpy.asarray([1 / Ntau ** 2 for k, Ntau in enumerate(Ntaus)])
    for j in range(output_tauT_len):
        ydata = [interpolation[0][j] for interpolation in interpolations]
        edata = [interpolation[1][j] for interpolation in interpolations]
        fitparams, fitparams_err = bootstr.bootstr_from_gauss(ce.fit_sample, data=ydata, data_std_dev=edata, numb_samples=nsamples,
                                                              sample_size=1, return_sample=False, args=[xdata, edata, [0, 3]])

        cont_corr[j] = fitparams[1]
        cont_corr_err[j] = fitparams_err[1]
    return cont_corr, cont_corr_err


def main():

    basepath = "/work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/EE/multi-level_2015/"

    nsamples = 10000
    tauT_data = []
    corr_data = []
    corr_err_data = []
    # lattices = ["192_48", "144_36", "96_24", "80_20", "64_16"]  # ,
    Ntaus = [48, 36, 24, 20, 16]
    for Ntau in Ntaus:
        # for latt,Ntau in zip(lattices,Ntaus):
        tauT, corr, err = numpy.loadtxt(basepath+"EE_2015_Nt" + str(Ntau) + ".dat", unpack=True)
        # tmp = numpy.loadtxt("./diffusion_"+latt+".dat", unpack=True)
        nt_half = int(Ntau / 2)
        tauT_data.append(tauT[:nt_half])
        corr_data.append(corr[:nt_half])
        corr_err_data.append(err[:nt_half])

    # output_tauTs = numpy.asarray(lpd.tauT_imp[36]) #usually these should be the ones of the finest lattice
    output_tauTs = tauT_data[1]  # usually these should be the ones of the finest lattice
    cont_corr, cont_corr_err = cont_extr(Ntaus, tauT_data, corr_data, corr_err_data, output_tauTs, nsamples)

    # save corr estimate in file and plot it
    numpy.savetxt(basepath+"./EE_2015_new_2022.txt", numpy.stack((output_tauTs, cont_corr, cont_corr_err), axis=-1), header="#tauT    corr    corr_err")
    fig, ax, plots = lpd.create_figure(xlims=[0, 0.51], ylims=[0, 4], xlabel=r'$\tau T$', ylabel=r'$ \frac{G (\tau)}{G^\mathrm{ norm } (\tau)}$')
    ax.errorbar(output_tauTs, cont_corr, cont_corr_err, color='black')
    fig.savefig(basepath+"./EE_2015_new_cont_2022.pdf")


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
