#!/bin/bash

# plot perturbative corrs
./plotting/00_plot_QED_LPT.py --inputfolder ../data/merged/ --outputfolder ../plots/
./plotting/0_plot_pert_correlators.py --Ntau 24 --inputfolder "../data/merged/quenched_pertLO_wilsonFlow/EE/"
./process_data/plot_pert_latt.py --Nt 20 --corr EE --flowtime_file ../data/merged/pert_LO/flowtimes.dat --outputpath ../plots/pertLO/ --inputpath ../data/merged/pert_LO/ --tau 5
./process_data/plot_pert_latt.py --Nt 20 --corr EE --flowtime_file ../data/merged/pert_LO/flowtimes.dat --outputpath ../plots/pertLO/ --inputpath ../data/merged/pert_LO/ --tau 5 --exp


qcdtypes=("quenched_1.5Tc_zeuthenFlow" "hisq_ms5_zeuthenFlow")
for qcdtype in "${qcdtypes[@]}"; do
# TODO add BB
./process_data/example_usage/1_merge_data.sh $qcdtype EE
./process_data/example_usage/1_xestimate_autocorrelations.sh $qcdtype EE
./process_data/example_usage/3_spline_interpolate.sh $qcdtype
./process_data/example_usage/4_continuum_extr.sh $qcdtype
./process_data/example_usage/5_flowtime_extr.sh $qcdtype

done