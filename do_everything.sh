#!/bin/bash

# TODO update all the paths here!

# plot perturbative corrs
#./perturbative_corr/plot_QED_LPT.py --inputfolder ../data/merged/ --outputfolder ../plots/
#./plotting/0_plot_pert_correlators.py --Ntau 24 --inputfolder "../data/merged/quenched_pertLO_wilsonFlow/EE/"

#for tau in 1 2 3 4 5 6 7 8 9 10 ; do
#    ./perturbative_corr/plot_tree_level_imp.py --Nt 30 --corr EE --flowtime_file ../data/merged/pert_LO/flowtimes.dat --outputpath ../plots/pertLO/ --inputpath ../data/merged/pert_LO/ --tau $tau
#done

# TODO change input phiUV to npy file
#./plotting/plot_integrand.py --outputpath ../plots/ --PathPhiUV ../data/merged/quenched_1.50Tc_zeuthenFlow/EE/spf/max_LO_T0.472_min2piT_w1_500_0.0_exp0_quenched_cont_f0/phiUV.dat

#./correlator_analysis/plotting/plot_EEvsBB.py
#./process_data/covariance.py --qcdtype quenched_1.50Tc_zeuthenFlow --conftype s144t36_b0754400 --corr EE --basepath ../data/merged/ --outputfolder ../plots/quenched_1.50Tc_zeuthenFlow/EE/
#./plotting/plot_kappaB.py --outputfolder ../plots/quenched_1.50Tc_zeuthenFlow/BB/


#./correlator_analysis/double_extrapolation/BB_renormalization/example_usage/extrapolate_coupling.sh
#
#
#qcdtypes=("quenched_1.5Tc_zeuthenFlow" "hisq_ms5_zeuthenFlow")
#corrs=("EE" "BB")
#for corr in "${corrs[@]}" ; do
#    for qcdtype in "${qcdtypes[@]}"; do
#        if [ "$qcdtype" == "hisq_ms5_zeuthenFlow" ] && [ "$corr" == "BB" ] ; then
#            continue
#        fi
#        ./process_data/example_usage/1_merge_data.sh $qcdtype $corr
#        ./process_data/example_usage/1_xestimate_autocorrelations.sh $qcdtype $corr
#        ./process_data/example_usage/3_spline_interpolate.sh $qcdtype $corr
#        ./process_data/example_usage/4_continuum_extr.sh $qcdtype $corr
#        ./process_data/example_usage/5_flowtime_extr.sh $qcdtype $corr
#    done
#done

#./correlator_analysis/plotting/plot_flow_correlations.py --qcdtype quenched_1.50Tc_zeuthenFlow --corr EE --conftype s144t36_b0754400 --basepath "../data/merged/" --outputfolder "/work/home/altenkort/work/correlators_flow/plots/quenched_1.50Tc_zeuthenFlow/EE/"


#./correlator_analysis/plotting/example_usage/2_plot_lateffects.sh quenched_1.50Tc_zeuthenFlow EE
#./correlator_analysis/plotting/example_usage/2_plot_lateffects.sh quenched_1.50Tc_zeuthenFlow BB
#./correlator_analysis/plotting/example_usage/2_plot_lateffects.sh hisq_ms5_zeuthenFlow EE
#
./correlator_analysis/plotting/6_plot_finalcorr.py --outputfolder ../plots/quenched_1.50Tc_zeuthenFlow/EE/ --input_flow ../data/merged/quenched_1.50Tc_zeuthenFlow/EE/EE_flow_extr_relflow.txt --input_multilvl ../data/merged/quenched_1.50Tc_zeuthenFlow/EE/multi-level_2015/EE_2015_new.txt
#
#./correlator_analysis/plotting/example_usage/plot_flow_dep.sh
#
#
#./spf_reconstruction/plot_fits/example_usage/plot_UV_spfs.sh
#./spf_reconstruction/plot_fits/example_usage/plot_fits_quenched.sh
#./spf_reconstruction/plot_fits/example_usage/plot_fits_hisq.sh
#./spf_reconstruction/plot_fits/example_usage/plot_final_kappas.sh
