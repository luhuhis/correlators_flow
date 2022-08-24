#!/usr/local/bin/python3.7 -u
import numpy
from scipy import integrate
import lib_process_data as lpd


def find_index(array, value):
    array = numpy.asarray(array)
    return (numpy.abs(array - value)).argmin()


def main():
    lowtauFT2 = 0.002719  # this comes from choosing overline{mu} such that Z_B =1
    inputpath = "/work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/coupling/"
    outputpath = "/work/home/altenkort/work/correlators_flow/plots/quenched_1.50Tc_zeuthenFlow/coupling/"
    lpd.create_folder(outputpath)
    confdict = {"16": "s064t64_b0687361", "20": "s080t80_b0703500", "24": "s096t96_b0719200", "30": "s096t120_b0739400", "36": "s096t144_b0754400"}
    # confdict = {"64": "s064t64_b0687361", "80": "s080t80_b0703500", "96": "s096t96_b0719200", "120": "s096t120_b0739400", "144": "s096t144_b0754400"}
    Nts = (20, 24, 30, 36)
    # Nts = (64, 80, 96, 120, 144)
    tauFT2_arr = []
    prefactor = 128 * numpy.pi ** 2 / 24
    g2_arr = []
    g2_err_arr = []
    integrand_arr = []
    integrand_err_arr = []
    for Nt in Nts:
        tauFbya2, tauF2E, tauF2E_err = numpy.loadtxt(inputpath + "flow_t2E_" + confdict[str(Nt)] + ".dat", unpack=True)
        tauFT2_arr.append(tauFbya2/Nt**2)
        g2_arr.append(prefactor * tauF2E)
        g2_err_arr.append(prefactor * tauF2E_err)
        integrand_arr.append((-3*g2_arr[-1])/(8*numpy.pi**2) / tauFT2_arr[-1])
        integrand_err_arr.append((-3*g2_err_arr[-1])/(8*numpy.pi**2) / tauFT2_arr[-1])

    fig, ax, _ = lpd.create_figure(xlabel=r'$\tau_F T^2$', xlabelpos=(1.1, 0.05), ylabelpos=(0.1, 0.9))  # xlims=(0, 0.0032), xlims=(0, 0.0032))
    # ax.set_yscale('log')

    ax.set_title(r'$ \frac{1}{\tau_FT^2} \left(\frac{-3g^2(\tau_F)}{8\pi^2}\right), \quad T=1.5T_c$', y=0.6)
    #  r'$\tau_F T^2 \frac{\langle E \rangle}{T^4} $'

    for i, Nt in enumerate(Nts):
        ax.errorbar(tauFT2_arr[i], integrand_arr[i], integrand_err_arr[i], fmt='x-', lw=0.5, markersize=2, mew=0.3, capsize=1, label=str(Nt), zorder=-Nt)
    ax.axvline(x=lowtauFT2, **lpd.verticallinestyle)
    ax.legend(**lpd.legendstyle, title=r'$N_\tau$')
    fig.savefig(outputpath+"integrand.pdf")

    Z2_arr = []
    x_arr = []

    for i, Nt in enumerate(Nts):
        lowindex = find_index(tauFT2_arr[i], lowtauFT2)

        temp_Z2 = []
        temp_x = []
        for tfT2 in tauFT2_arr[i]:
            highindex = find_index(tauFT2_arr[i], tfT2)

            if highindex < lowindex:
                y = integrand_arr[i][highindex:lowindex+1]
                x = tauFT2_arr[i][highindex:lowindex+1]
                sign = -1
            else:
                y = integrand_arr[i][lowindex:highindex+1]
                x = tauFT2_arr[i][lowindex:highindex+1]
                sign = 1
            # TODO: DO BOOTSTRAP FOR THE ERRORS ?
            integral = integrate.trapz(y, x)
            temp_Z2.append(numpy.exp(sign * integral))

            temp_x.append(tfT2)

        Z2_arr.append(temp_Z2)
        x_arr.append(temp_x)

        numpy.savetxt(inputpath + "Z2_" + str(Nt) + ".dat", numpy.stack((x_arr[i], Z2_arr[i]), axis=-1), header="tf T^2        Z^2")

    fig2, ax2, plots2 = lpd.create_figure(xlabel=r'$\tau_F T^2$', xlabelpos=(1.1, 0.05), ylabelpos=(0.1, 0.9))  # xlims=(0, 0.0032)
    # ax2.set_yscale('log')

    ax2.set_title(r'$Z_f^2$',  y=0.85, x=0.12)
    labeldict = {20: "7.035", 24: "7.192", 30: "7.394", 36: "7.544"}


    for i, Nt in enumerate(Nts):
        ax2.errorbar(x_arr[i], Z2_arr[i], fmt='x-', lw=0.5, markersize=2, mew=0.3, capsize=1, label=labeldict[Nt], zorder=-Nt)
    ax2.axvline(x=lowtauFT2, **lpd.verticallinestyle)
    ax2.axhline(y=1, **lpd.horizontallinestyle)
    ax2.legend(**lpd.legendstyle, title=r'$\beta$')
    fig2.savefig(outputpath+"Zf2.pdf")


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
