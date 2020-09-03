#!/bin/bash

#plot all figures from the paper
./plot_pert_correlators.py 20
./1_plot_correlator.py quenched s144t36_b0754400
./2_plot_lateffects.py 60
./2_plot_lateffects.py 100
../process_data/_4_continuum_extr.py 85
../process_data/_5_flowtime_extr.py
./6_plot_finalcorr.py
