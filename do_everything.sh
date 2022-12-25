#!/bin/bash

# plot perturbative corrs
./plotting/00_plot_QED_LPT.py --inputfolder ../data/merged/ --outputfolder ../plots/
./plotting/0_plot_pert_correlators.py --Ntau 24 --inputfolder "../data/merged/quenched_pertLO_wilsonFlow/EE/"
./process_data/plot_pert_latt.py --Nt 20 --corr EE --flowtime_file ../data/merged/pert_LO/flowtimes.dat --outputpath ../plots/pertLO/ --inputpath ../data/merged/pert_LO/ --tau 5
./process_data/plot_pert_latt.py --Nt 20 --corr EE --flowtime_file ../data/merged/pert_LO/flowtimes.dat --outputpath ../plots/pertLO/ --inputpath ../data/merged/pert_LO/ --tau 5 --exp

# TODO change input phiUV to npy file
./plotting/plot_integrand.py --outputpath ../plots/ --PathPhiUV ../data/merged/quenched_1.50Tc_zeuthenFlow/EE/spf/max_LO_T0.472_min2piT_w1_500_0.0_exp0_quenched_cont_f0/phiUV.dat

./plotting/plot_EEvsBB.py
./process_data/covariance.py --qcdtype quenched_1.50Tc_zeuthenFlow --conftype s144t36_b0754400 --corr EE --basepath ../data/merged/ --outputfolder ../plots/quenched_1.50Tc_zeuthenFlow/EE/
./plotting/plot_kappaB.py --outputfolder ../plots/quenched_1.50Tc_zeuthenFlow/BB/


qcdtypes=("quenched_1.5Tc_zeuthenFlow" "hisq_ms5_zeuthenFlow")
corrs=("EE" "BB")
for corr in "${corrs[@]}" ; do
    for qcdtype in "${qcdtypes[@]}"; do
        ./process_data/example_usage/1_merge_data.sh $qcdtype $corr
        ./process_data/example_usage/1_xestimate_autocorrelations.sh $qcdtype EE
        ./process_data/example_usage/3_spline_interpolate.sh $qcdtype
        ./process_data/example_usage/4_continuum_extr.sh $qcdtype
        ./process_data/example_usage/5_flowtime_extr.sh $qcdtype
    done
done

./plotting/6_plot_finalcorr.py --outputfolder ../plots/quenched_1.50Tc_zeuthenFlow/EE/ --input_flow ../data/merged/quenched_1.50Tc_zeuthenFlow/EE/EE_flow_extr.txt --input_multilvl ../data/merged/quenched_1.50Tc_zeuthenFlow/EE/multi-level_2015/EE_2015_new.txt