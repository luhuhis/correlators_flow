#!/usr/bin/env python3

import lib_process_data as lpd
import numpy
import scipy.signal
from latqcdtools.statistics import statistics as stat
from latqcdtools.statistics import bootstr

def reduce_function(data):
    return numpy.mean(data)


def formatfloat(number):
    return '{0:.5f}'.format(number)


def main():

    n = 10000
    x = range(n)

    noise = numpy.random.normal(0, 1, n)
    truetauint = 30
    kernel = [1/truetauint for _ in range(truetauint)]  # [1, 1, 1, 1, 1 .... 1]
    conv = scipy.signal.fftconvolve(noise, kernel, mode='valid')

    fig, ax, plots = lpd.create_figure()

    tpickmax = 50
    nblocks = int(len(conv)/tpickmax)
    tau_int, tau_inte, tau_intbias, itpick = stat.getTauInt(conv, nblocks, tpickmax, acoutfileName='acor.d', showPlot=True)
    print(tau_int, tau_inte, tau_intbias, itpick)

    nplot = 1000
    ax.plot(x[:nplot], conv[:nplot]/len(conv), label="correlated signal")

    blockeddata = []
    blocklength = truetauint*2
    nblocks = int(len(conv)/blocklength)
    for i in range(nblocks):
        blockeddata.append(numpy.mean(conv[i*blocklength:(i+1)*blocklength]))

    blockeddata = numpy.asarray(blockeddata)

    print("true values:                                       ", formatfloat(0), formatfloat(1/numpy.sqrt(n)))
    # print("direct calculation (mean and std error):           ", formatfloat(numpy.mean(noise)), formatfloat(numpy.std(noise, ddof=1) / numpy.sqrt(numpy.size(noise))))
    print("direct calculation (mean and std error):           ", formatfloat(numpy.mean(conv)),
          formatfloat(numpy.std(conv, ddof=1) / numpy.sqrt(numpy.size(conv))))
    XX_samples, XX, XX_err = bootstr.bootstr(reduce_function, conv, numb_samples=10000, sample_size=len(conv), conf_axis=0,
                                             return_sample=True, parallelize=True, nproc=4, seed=0, err_by_dist=True)
    print("independent bootstrap, median and 68th percentiles:", formatfloat(XX), formatfloat(XX_err))

    # ax.hist(XX_samples, bins="fd")
    # matplotlib.pyplot.show()

    XX, XX_err = bootstr.bootstr(reduce_function, conv, numb_samples=10000, sample_size=len(conv), conf_axis=0,
                                 same_rand_for_obs=True, parallelize=True, nproc=4, seed=0, err_by_dist=False)
    print("independent bootstrap, mean and standard deviation:", formatfloat(XX), formatfloat(XX_err))

    XX, XX_err = bootstr.bootstr(reduce_function, blockeddata, numb_samples=10000, sample_size=len(blockeddata), conf_axis=0,
                                             same_rand_for_obs=True, parallelize=True, nproc=4, seed=0, err_by_dist=False)

    print("blocked bootstrap,     median and 68th percentiles:", formatfloat(XX), formatfloat(XX_err))
    XX, XX_err = bootstr.bootstr(reduce_function, blockeddata, numb_samples=10000, sample_size=len(blockeddata), conf_axis=0,
                                             same_rand_for_obs=True, parallelize=True, nproc=4, seed=0, err_by_dist=True)
    print("blocked bootstrap:     mean and standard deviation:", formatfloat(XX), formatfloat(XX_err))

    fig.savefig("example2.pdf")

    return


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
