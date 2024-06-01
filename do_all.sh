#!/bin/bash

# tar -xzf correlators_flow.tar.gz   
# OR
git clone --branch thesis https://github.com/luhuhis/correlators_flow.git
chmod -R +x ./correlators_flow
# tar -xzf AnalysisToolbox.tar.gz 
# OR
git clone https://github.com/LatticeQCD/AnalysisToolbox.git
( cd AnalysisToolbox && git checkout f9eee73d45d7b981153db75cfaf2efa2b4cefa9c )
tar -xzf data.tar.gz
tar -xzf figures.tar.gz
tar -xzf output_data.tar.gz
export PYTHONPATH=$(pwd)/correlators_flow:$(pwd)/AnalysisToolbox:${PYTHONPATH} BASEPATH_RAW_DATA=$(pwd)/input BASEPATH_WORK_DATA=$(pwd)/output_data BASEPATH_PLOT=$(pwd)/figures
export G_PERT_LO_DIR=${BASEPATH_RAW_DATA}/quenched_1.50Tc_zeuthenFlow/pert_LO
export NPROC=20
poetry install
./correlators_flow/perturbative_corr/plot_QED_LPT.py --inputfolder ${BASEPATH_RAW_DATA} --outputfolder ${BASEPATH_PLOT}
./correlators_flow/perturbative_corr/plot_pert_correlators.py --Ntau 24 --inputfolder ${BASEPATH_RAW_DATA}/quenched_1.50Tc_zeuthenFlow/pert_LO/ --outputfolder ${BASEPATH_PLOT}/pertLO
./correlators_flow/perturbative_corr/plot_tree_level_imp.py --Nt 30 --corr EE --flowtime_file ${BASEPATH_RAW_DATA}/quenched_1.50Tc_zeuthenFlow/pert_LO/flowtimes.dat --outputpath ${BASEPATH_PLOT}/pertLO/ --inputpath ${BASEPATH_RAW_DATA}/quenched_1.50Tc_zeuthenFlow/pert_LO/ --tau 5
./correlators_flow/perturbative_corr/plot_tree_level_imp.py --Nt 30 --corr EE --flowtime_file ${BASEPATH_RAW_DATA}/quenched_1.50Tc_zeuthenFlow/pert_LO/flowtimes.dat --outputpath ${BASEPATH_PLOT}/pertLO/ --inputpath ${BASEPATH_RAW_DATA}/quenched_1.50Tc_zeuthenFlow/pert_LO/ --tau 10
./correlators_flow/correlator_analysis/double_extrapolation/example_usage/1_merge_data.sh quenched_1.50Tc_zeuthenFlow EE ${BASEPATH_RAW_DATA} ${BASEPATH_WORK_DATA}
./correlators_flow/correlator_analysis/double_extrapolation/example_usage/1_merge_data.sh quenched_1.50Tc_zeuthenFlow BB ${BASEPATH_RAW_DATA} ${BASEPATH_WORK_DATA}
./correlators_flow/correlator_analysis/double_extrapolation/example_usage/1_merge_data.sh hisq_ms5_zeuthenFlow EE ${BASEPATH_RAW_DATA} ${BASEPATH_WORK_DATA} ${BASEPATH_RAW_DATA}/hisq_ms5_zeuthenFlow/reference_flowtimes
./correlators_flow/correlator_analysis/double_extrapolation/example_usage/2_reduce_data.sh quenched_1.50Tc_zeuthenFlow EE ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT} ${NPROC}
./correlators_flow/correlator_analysis/double_extrapolation/example_usage/2_reduce_data.sh quenched_1.50Tc_zeuthenFlow BB ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT} ${NPROC}
./correlators_flow/correlator_analysis/double_extrapolation/example_usage/2_reduce_data.sh hisq_ms5_zeuthenFlow EE ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT} ${NPROC}
./correlators_flow/correlator_analysis/plotting/example_usage/2_plot_lateffects.sh quenched_1.50Tc_zeuthenFlow EE ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT} ${NPROC}
./correlators_flow/correlator_analysis/plotting/example_usage/2_plot_lateffects.sh hisq_ms5_zeuthenFlow EE ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT} ${NPROC}
./correlators_flow/correlator_analysis/plotting/example_usage/plot_flow_dep.sh quenched_1.50Tc_zeuthenFlow EE ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT}
./correlators_flow/correlator_analysis/plotting/example_usage/plot_flow_dep.sh quenched_1.50Tc_zeuthenFlow BB ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT}
./correlators_flow/correlator_analysis/plotting/example_usage/plot_flow_dep.sh hisq_ms5_zeuthenFlow EE ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT}
./correlators_flow/correlator_analysis/double_extrapolation/example_usage/3_spline_interpolate.sh quenched_1.50Tc_zeuthenFlow EE ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT} ${NPROC}
./correlators_flow/correlator_analysis/double_extrapolation/example_usage/3_spline_interpolate.sh quenched_1.50Tc_zeuthenFlow BB ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT} ${NPROC}
./correlators_flow/correlator_analysis/double_extrapolation/example_usage/3_spline_interpolate.sh hisq_ms5_zeuthenFlow EE ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT} ${NPROC}
./correlators_flow/correlator_analysis/double_extrapolation/example_usage/4_continuum_extr.sh quenched_1.50Tc_zeuthenFlow EE ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT} ${NPROC}
./correlators_flow/correlator_analysis/double_extrapolation/example_usage/4_continuum_extr.sh quenched_1.50Tc_zeuthenFlow BB ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT} ${NPROC}
./correlators_flow/correlator_analysis/double_extrapolation/example_usage/4_continuum_extr.sh hisq_ms5_zeuthenFlow EE ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT} ${NPROC}
./correlators_flow/correlator_analysis/plotting/plot_flow_correlations.py --qcdtype quenched_1.50Tc_zeuthenFlow --corr EE --conftype s144t36_b0754400 --basepath ${BASEPATH_WORK_DATA} --outputfolder ${BASEPATH_PLOT}/quenched_1.50Tc_zeuthenFlow/EE/ --nproc ${NPROC}
./correlators_flow/correlator_analysis/double_extrapolation/example_usage/5_flowtime_extr.sh quenched_1.50Tc_zeuthenFlow EE ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT} ${NPROC}
./correlators_flow/correlator_analysis/double_extrapolation/example_usage/5_flowtime_extr.sh hisq_ms5_zeuthenFlow EE ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT} ${NPROC}
./correlators_flow/multi-level/cont_extr_new.py --basepath ${BASEPATH_RAW_DATA}
./correlators_flow/correlator_analysis/plotting/6_plot_finalcorr.py --outputfolder ${BASEPATH_PLOT}/quenched_1.50Tc_zeuthenFlow/EE/ --input_flow ${BASEPATH_WORK_DATA}/quenched_1.50Tc_zeuthenFlow/EE/EE_flow_extr_relflow.txt --input_multilvl ${BASEPATH_RAW_DATA}/multi-level_2015/EE_2015_new.txt
./correlators_flow/correlator_analysis/double_extrapolation/BB_renormalization/example_usage/extrapolate_coupling.sh ${BASEPATH_RAW_DATA} ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT} 6.40
./correlators_flow/correlator_analysis/double_extrapolation/BB_renormalization/example_usage/compute_Z.sh ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT}
./correlators_flow/correlator_analysis/double_extrapolation/example_usage/5_flowtime_extr.sh quenched_1.50Tc_zeuthenFlow BB ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT} ${NPROC}
./correlators_flow/correlator_analysis/plotting/plot_EEvsBB.py --inputfolder ${BASEPATH_WORK_DATA}/quenched_1.50Tc_zeuthenFlow/ --outputfolder ${BASEPATH_PLOT}/quenched_1.50Tc_zeuthenFlow/
./correlators_flow/correlator_analysis/relative_flow/example_usage/Nf0.sh ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT}
./correlators_flow/correlator_analysis/relative_flow/example_usage/Nf3TemperatureComparison_Paper.sh ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT}
./correlators_flow/spf_reconstruction/plot_fits/example_usage/plot_g2.sh ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT}
./correlators_flow/spf_reconstruction/plotting/plot_integrand.py --outputpath ${BASEPATH_PLOT} --Nf 0 --min_scale eff --T_in_GeV 0.472 --omega_prefactor "1" --order LO --corr EE --mu_IR_by_T 1
./correlators_flow/spf_reconstruction/model_fitting/example_usage/spf_reconstruct.sh quenched_1.50Tc_zeuthenFlow EE  ${BASEPATH_WORK_DATA} NO ${NPROC}
./correlators_flow/spf_reconstruction/model_fitting/example_usage/spf_reconstruct.sh quenched_1.50Tc_zeuthenFlow BB  ${BASEPATH_WORK_DATA} NO ${NPROC}
./correlators_flow/spf_reconstruction/model_fitting/example_usage/spf_reconstruct.sh hisq_ms5_zeuthenFlow EE  ${BASEPATH_WORK_DATA} NO ${NPROC}
./correlators_flow/spf_reconstruction/plot_fits/example_usage/plot_fits_quenched.sh EE ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT} yes
./correlators_flow/spf_reconstruction/plot_fits/example_usage/plot_fits_quenched.sh BB ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT} yes
./correlators_flow/spf_reconstruction/plot_fits/example_usage/plot_fits_quenched.sh BB ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT} no
./correlators_flow/spf_reconstruction/plot_fits/example_usage/plot_fits_hisq.sh ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT}
./correlators_flow/spf_reconstruction/plot_fits/example_usage/plot_kfactors.sh ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT}
./correlators_flow/spf_reconstruction/plot_fits/example_usage/plot_final_kappas.sh ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT} EE
./correlators_flow/spf_reconstruction/plot_fits/example_usage/plot_final_kappas.sh ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT} BB
./correlators_flow/spf_reconstruction/plot_fits/example_usage/plot_final_kappas.sh ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT} hisq_thesis
./correlators_flow/spf_reconstruction/plot_fits/publication_specific/2024-BB-paper/fit_kappa_to_g2_g4.py --outputpath ${BASEPATH_PLOT}/quenched_1.50Tc_zeuthenFlow/
./correlators_flow/spf_reconstruction/plotting/plot_2piTD.py --outputfolder ${BASEPATH_PLOT}
