#!/bin/bash

# quenched
# general quenched comparison
params=(" \
--output_suffix _quenched_0 \
--fitparam_files \
5_a_300_0.24_exp0_quenched_1.5Tc_cont_tauF0_relflow_0.20_fit \
--no_kappa_plot \
" " \
--output_suffix _quenched_1 \
--fitparam_files \
5_a_300_0.24_exp0_quenched_1.5Tc_cont_tauF0_relflow_0.20_fit \
5_a_300_1e-04_d_0.2_relflow_0.20_piT_1.0 \
--tauT_vlines 0.25 \
" " \
--output_suffix _quenched_2 \
--fitparam_files \
5_a_300_0.24_exp0_quenched_1.5Tc_cont_tauF0_relflow_0.20_fit \
5_a_300_1e-04_d_0.2_relflow_0.20_piT_1.0 \
5_a_300_1e-04_d_0.4_relflow_0.20_piT_1.0 \
--tauT_vlines 0.25 0.416 \
" " \
--output_suffix _quenched_3 \
--fitparam_files \
5_a_300_0.24_exp0_quenched_1.5Tc_cont_tauF0_relflow_0.20_fit \
5_a_300_1e-04_d_0.2_relflow_0.20_piT_1.0 \
5_a_300_1e-04_d_0.4_relflow_0.20_piT_1.0 \
5_a_300_1e-04_d_0.4_quenched_1.5Tc_Nt36_relflow_0.20_piT_1.0 \
--tauT_vlines 0.25 0.416 0.4167 \
--conftype s144t36_b0754400 \
" " \
--output_suffix _quenched_4 \
--fitparam_files \
5_a_300_0.24_exp0_quenched_1.5Tc_cont_tauF0_relflow_0.20_fit \
5_a_300_1e-04_d_0.2_relflow_0.20_piT_1.0 \
5_a_300_1e-04_d_0.4_relflow_0.20_piT_1.0 \
5_a_300_1e-04_d_0.4_quenched_1.5Tc_Nt36_relflow_0.20_piT_1.0 \
5_a_300_1e-04_d_0.4_quenched_1.5Tc_Nt30_relflow_0.20_piT_1.0 \
--conftype s144t36_b0754400 s120t30_b0739400 \
--tauT_vlines 0.25 0.416 0.4167 0.433 \
" " \
--output_suffix _quenched_5 \
--fitparam_files \
5_a_300_0.24_exp0_quenched_1.5Tc_cont_tauF0_relflow_0.20_fit \
5_a_300_1e-04_d_0.2_relflow_0.20_piT_1.0 \
5_a_300_1e-04_d_0.4_relflow_0.20_piT_1.0 \
5_a_300_1e-04_d_0.4_quenched_1.5Tc_Nt36_relflow_0.20_piT_1.0 \
5_a_300_1e-04_d_0.4_quenched_1.5Tc_Nt30_relflow_0.20_piT_1.0 \
5_a_500_0.38_exp2_quenched_1.5Tc_Nt20_relflow0.20_piT_1.0w \
--conftype s144t36_b0754400 s120t30_b0739400 s080t20_b0703500 \
--tauT_vlines 0.25 0.416 0.4167 0.433 \
")

for param in "${params[@]}" ; do
./plot_rec_corr_fixFlowBytauT.py  --npoints 100 --ylims 2.5 3.75 \
$param \
--xlims 0.2 0.512 --kappa_xlims -0.1 2.5 \
--title ", $\\:\\: \\mu=\\sqrt{[\\pi T\\,]^2+\\omega^2}$" \
--no_connection \
--model 5 \
--plot_quenched_extr \
--flowradiusBytauT 0.2 \
--qcdtype quenched_1.50Tc_zeuthenFlow \
--corr EE \
--deduce_fitparam_files --fitparam_basepath ~/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/EE/spf/ \
--PhiUV_basepath ~/work/correlators_flow/data/merged/spf_coupling/ \
--PhiUV_files quenched_1.5Tc_SPF_LO_Nf0_0.472_piT_1.0_smax.dat \
--output_path ~/work/correlators_flow/plots/quenched_1.50Tc_zeuthenFlow/EE/ \
--color_data k C0 C1 C2 C3 \
--color_fit k C0 C9 C1 C2 C3 \
--color_vlines k C9 C1 C2 \
--kappa_ypos 1.75 1.45 1.3 1 0.85 0.6 \
--leg_n_dummies 0 \
--min_tauT 0.250 0.250 0.416 0.416 0.433 \
--min_tauT_plot 0.05 \
--no_label \
--no_just_UV \
--leg_pos 0.15 1 \
&
done


./plot_rec_corr_fixFlowBytauT.py  --npoints 100 --ylims 2.5 3.75 \
--output_suffix _quenched_5_IMP \
--conftype s144t36_b0754400 s080t20_b0703500 \
--tauT_vlines 0.25 0.416 0.4167 0.433 \
--xlims 0.2 0.512 --kappa_xlims -0.1 2.5 \
--title ", $\\:\\: \\mu=\\sqrt{[\\pi T\\,]^2+\\omega^2}$" \
--no_connection \
--model 5 \
--plot_quenched_extr \
--flowradiusBytauT 0.2 \
--qcdtype quenched_1.50Tc_zeuthenFlow \
--corr EE \
--deduce_fitparam_files --fitparam_basepath ~/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/EE/spf/ \
--PhiUV_basepath ~/work/correlators_flow/data/merged/spf_coupling/ \
--PhiUV_files quenched_1.5Tc_SPF_LO_Nf0_0.472_piT_1.0_smax.dat \
--output_path ~/work/correlators_flow/plots/quenched_1.50Tc_zeuthenFlow/EE/ \
--color_data k C0 C1 C3 \
--color_fit k C0 C9 C1 C3 \
--color_vlines k C9 C1 C2 \
--kappa_ypos 1.75 1.45 1.3 1 0.6 \
--leg_n_dummies 0 \
--min_tauT 0.250 0.250 0.416 0.433 \
--min_tauT_plot 0.05 \
--no_label \
--no_just_UV \
--leg_pos 0.15 1 \
