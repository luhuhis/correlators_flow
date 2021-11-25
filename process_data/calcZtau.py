#!/usr/local/bin/python3.7
import numpy
from scipy import integrate
import lib_process_data as lpd


def find_index(array, value):
    array = numpy.asarray(array)
    return (numpy.abs(array - value)).argmin()


def main():
    lowtfT2 = 0.002719
    path = "/work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/coupling/"
    confdict = {"16": "s064t64_b0687361", "20": "s080t80_b0703500", "24": "s096t96_b0719200", "30": "s096t120_b0739400", "36": "s096t144_b0754400"}
    Nts = (16, 20, 24, 30, 36)
    tfT2_arr = []
    EbyT4_arr = []
    EbyT4_err_arr = []
    integrand_arr = []
    integrand_err_arr = []
    for Nt in Nts:
        tbya2, t2E, t2E_err = numpy.loadtxt(path + "flow_t2E_" + confdict[str(Nt)] + ".dat", unpack=True)
        tfT2_arr.append(tbya2/Nt**2)
        EbyT4_arr.append(t2E / tbya2**2 * Nt**4)
        EbyT4_err_arr.append(t2E_err / tbya2**2 * Nt**4)
        integrand_arr.append(t2E / tbya2**2 * Nt**4 * tbya2/Nt**2)
        integrand_err_arr.append(t2E_err / tbya2**2 * Nt**4 * tbya2/Nt**2)

    fig, ax, plots = lpd.create_figure(xlabel=r'$\tau_F T^2$', xlabelpos=(1.1, 0.05), ylabelpos=(0.1, 0.9), xlims=(0, 0.0032))
    # ax.set_yscale('log')

    ax.set_title(r'$ \frac{1}{\tau_F} \left(\frac{-3g^2(\tau_F)}{8\pi^2}\right)$', y=0.6)
    #  r'$\tau_F T^2 \frac{\langle E \rangle}{T^4} $'

    for i, Nt in enumerate(Nts):
        ax.errorbar(tfT2_arr[i], -2*EbyT4_arr[i]*tfT2_arr[i], -2*EbyT4_err_arr[i]*tfT2_arr[i], fmt='x-', lw=0.5, markersize=2, mew=0.3, capsize=1, label=str(Nt), zorder=-Nt)
    ax.axvline(x=lowtfT2, **lpd.verticallinestyle)
    ax.legend(**lpd.legendstyle, title=r'$N_\tau$')
    fig.savefig(path+"EbyT4.pdf")

    Z2_arr = []
    x_arr = []

    for i, Nt in enumerate(Nts):
        lowindex = find_index(tfT2_arr[i], lowtfT2)

        temp_Z2 = []
        temp_x = []
        for tfT2 in tfT2_arr[i]:
            highindex = find_index(tfT2_arr[i], tfT2)

            if highindex < lowindex:
                y = integrand_arr[i][highindex:lowindex+1]
                x = tfT2_arr[i][highindex:lowindex+1]
                sign = -1 * -1
            else:
                y = integrand_arr[i][lowindex:highindex+1]
                x = tfT2_arr[i][lowindex:highindex+1]
                sign = -1
            # TODO: DO BOOTSTRAP FOR THE ERRORS ?
            integral = integrate.trapz(y, x)
            temp_Z2.append(numpy.exp(sign * 2 * integral))

            temp_x.append(tfT2)

        Z2_arr.append(temp_Z2)
        x_arr.append(temp_x)

        numpy.savetxt(path + "Z2_" + str(Nt) + ".dat", numpy.stack((x_arr[i], Z2_arr[i]), axis=-1), header="tf T^2        Z^2")

    fig2, ax2, plots2 = lpd.create_figure(xlabel=r'$\tau_F T^2$', xlabelpos=(1.1, 0.05), ylabelpos=(0.1, 0.9), xlims=(0, 0.0032))
    # ax2.set_yscale('log')

    ax2.set_title(r'$Z_\tau^2 $',  y=0.85)

    for i, Nt in enumerate(Nts):
        ax2.errorbar(x_arr[i], Z2_arr[i], fmt='x-', lw=0.5, markersize=2, mew=0.3, capsize=1, label=str(Nt), zorder=-Nt)
    ax2.axvline(x=lowtfT2, **lpd.verticallinestyle)
    ax2.axhline(y=1, **lpd.horizontallinestyle)
    ax2.legend(**lpd.legendstyle, title=r'$N_\tau$')
    fig2.savefig(path+"Z2.pdf")


if __name__ == '__main__':
    main()
    lpd.save_script_call()
