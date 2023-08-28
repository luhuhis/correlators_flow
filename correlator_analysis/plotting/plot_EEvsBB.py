#!/usr/bin/env python3

import lib_process_data as lpd
import numpy
import argparse
from correlator_analysis.double_extrapolation.BB_renormalization.compute_Zf2 import get_scale_choices
import matplotlib.pyplot


class Plotter:
    def __init__(self, args):
        self.filename = args.outputfolder + "/EEvsBB.pdf"
        self.inputfolder = args.inputfolder
        self.plotstyle = dict(fmt='.', markersize=0, fillstyle='none')

    def _setup_plot(self):
        self.fig, self.ax, self.plots = lpd.create_figure(xlims=(0.235, 0.51), ylims=(2.5, 4.4), xlabel=r'$\tau T$',
                                                          ylabel=r'$\displaystyle \frac{G_X}{G^\mathrm{norm}}$')
        self.ax.errorbar(0, 0, label=r'$X,\, \frac{\bar{\mu}_\text{UV}}{\mu_\text{F}},\ \ \frac{\bar{\mu}_\text{IR}}{T}$',
                         markersize=0, alpha=0, lw=0)
        self.ax.set_prop_cycle(None)

    def _plot_BB(self):
        scale_choices = get_scale_choices()
        for scale_choice in scale_choices:
            BB = numpy.loadtxt(self.inputfolder + "/BB/BB_flow_extr_relflow_" + scale_choice.choice_label + ".txt", unpack=True)
            self.ax.errorbar(BB[0], BB[1], BB[2], **self.plotstyle,
                             label=r'$B, '+scale_choice.choice_label_for_plot_short + r'$')

    def _plot_EE(self):
        EE = numpy.loadtxt(self.inputfolder + "/EE/EE_flow_extr_relflow.txt", unpack=True)
        self.ax.errorbar(EE[0], EE[1], EE[2], **self.plotstyle, label=r'$E$')

    def _plot(self):
        self._setup_plot()
        self._plot_BB()
        self._plot_EE()
        self._finalize_plot()

    def _finalize_plot(self):
        self.ax.legend(**lpd.leg_err_size(), fontsize=8, title_fontsize=8, loc="lower right",
                       bbox_to_anchor=(1, 0.05), framealpha=0)

        self.fig.savefig(self.filename)
        print("saved correlator plot", self.filename)
        matplotlib.pyplot.close(self.fig)

    @classmethod
    def plot(cls, *args):
        # This method only exists such that you don't have to manually instantiate the class.
        instance = cls(*args)
        return instance._plot()


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--inputfolder', default="/work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/")
    parser.add_argument('--outputfolder', default="/work/home/altenkort/work/correlators_flow/plots/quenched_1.50Tc_zeuthenFlow/")
    args = parser.parse_args()

    Plotter.plot(args)

    return


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
