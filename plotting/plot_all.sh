#!/bin/bash

#plot all figures from the paper (15)

./0_plot_pert_correlators.py 20
#EE_pert_contvslatt_flow.pdf

./1_plot_correlator.py quenched s144t36_b0754400
#

./2_plot_lateffects.py 60
#

./2_plot_lateffects.py 100
#

../process_data/_4_continuum_extr.py 85
#

../process_data/_5_flowtime_extr.py
#

./6_plot_finalcorr.py
#

../process_data/_flowtime_extr_coarsest.py
#
