#!/bin/bash

# quenched
# general quenched comparison

params=(" \
--suffix _quenched_0 \
--fitparam_files \
5_a_300_0.24_exp0_quenched_1.5Tc_cont_tauF0_relflow_0.20_fit \
--no_kappa_plot \
" " \
--suffix _quenched_1 \
--fitparam_files \
5_a_300_0.24_exp0_quenched_1.5Tc_cont_tauF0_relflow_0.20_fit \
5_a_300_1e-04_d_0.2_relflow_0.20_piT_1.0 \
--tauT_vlines 0.25 \
" " \
--suffix _quenched_2 \
--fitparam_files \
5_a_300_0.24_exp0_quenched_1.5Tc_cont_tauF0_relflow_0.20_fit \
5_a_300_1e-04_d_0.2_relflow_0.20_piT_1.0 \
5_a_300_1e-04_d_0.4_relflow_0.20_piT_1.0 \
--tauT_vlines 0.25 0.416 \
" " \
--suffix _quenched_3 \
--fitparam_files \
5_a_300_0.24_exp0_quenched_1.5Tc_cont_tauF0_relflow_0.20_fit \
5_a_300_1e-04_d_0.2_relflow_0.20_piT_1.0 \
5_a_300_1e-04_d_0.4_relflow_0.20_piT_1.0 \
5_a_300_1e-04_d_0.4_quenched_1.5Tc_Nt36_relflow_0.20_piT_1.0 \
--tauT_vlines 0.25 0.416 0.4167 \
--conftype s144t36_b0754400 \
" " \
--suffix _quenched_4 \
--fitparam_files \
5_a_300_0.24_exp0_quenched_1.5Tc_cont_tauF0_relflow_0.20_fit \
5_a_300_1e-04_d_0.2_relflow_0.20_piT_1.0 \
5_a_300_1e-04_d_0.4_relflow_0.20_piT_1.0 \
5_a_300_1e-04_d_0.4_quenched_1.5Tc_Nt36_relflow_0.20_piT_1.0 \
5_a_300_1e-04_d_0.4_quenched_1.5Tc_Nt30_relflow_0.20_piT_1.0 \
--conftype s144t36_b0754400 s120t30_b0739400 \
--tauT_vlines 0.25 0.416 0.4167 0.433 \
" " \
--suffix _quenched_5 \
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
--xlims 0.2 0.512 --xlims_kappa -0.1 2.5 \
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
--outputpath ~/work/correlators_flow/plots/quenched_1.50Tc_zeuthenFlow/EE/ \
--color_data k C0 C1 C2 C3 \
--color_fit k C0 C9 C1 C2 C3 \
--color_vlines k C9 C1 C2 \
--kappa_ypos 1.75 1.45 1.3 1 0.85 0.6 \
--n_dummies 0 \
--min_tauT 0.250 0.250 0.416 0.416 0.433 \
--min_tauT_plot 0.05 \
--no_label \
--no_just_UV \
--leg_pos 0.15 1 \
&
done



# ================
# ================ HISQ
# ================
./plot_rec_corr_fixFlowBytauT.py --suffix _hisq --npoints 100 --ylims 1.8 9 \
--xlims 0.14 0.512 --xlims_kappa 0 12 \
--title ", $\\:\\:\\mu=\\sqrt{[2\\pi T\\,]^2+\\omega^2}$" \
--no_connection \
--model 5 \
--flowradiusBytauT 0.3 0.3 0.3 0.2 \
--qcdtype hisq_b8249_zeuthenFlow \
--corr EE \
--conftype \
s096t36_b0824900_m002022_m01011 \
s096t32_b0824900_m002022_m01011 \
s096t28_b0824900_m002022_m01011 \
s096t20_b0824900_m002022_m01011 \
--deduce_fitparam_files --fitparam_basepath ~/work/correlators_flow/data/merged/hisq_b8249_zeuthenFlow/EE/spf/ \
--fitparam_files \
5_a_500_0.4_exp2_hisq_nt36_relflow0.30_2piT_1w \
5_a_500_0.43_exp2_hisq_nt32_relflow0.30_2piT_1w \
5_a_500_0.42_exp2_hisq_nt28_relflow0.30_2piT_1w \
5_a_500_0.39_exp2_hisq_nt20_relflow0.20_2piT_1w \
--color_vlines k C0 C1 C2 \
--PhiUV_basepath ~/work/correlators_flow/data/merged/spf_coupling/ \
--PhiUV_files \
hisq_nt36_SPF_LO_Nf3_0.196_2piT_1.0_smax.dat \
hisq_nt32_SPF_LO_Nf3_0.220_2piT_1.0_smax.dat \
hisq_nt28_SPF_LO_Nf3_0.251_2piT_1.0_smax.dat \
hisq_nt20_SPF_LO_Nf3_0.352_2piT_1.0_smax.dat \
--outputpath ~/work/correlators_flow/plots/hisq_b8249_zeuthenFlow/EE/ \
--color_data k C0 C1 C2 \
--color_fit k C0 C1 C2 \
--kappa_ypos 1.45 1.3 1.15 1 0.85 \
--n_dummies 0 \
--tauT_vlines 0.4167 0.4375 0.4285 0.4 \
--min_tauT 0.4167 0.4375 0.4285 0.4 \
--min_tauT_plot 0.05 \
--no_label \
--kappa_labelheight 0.65 \
--leg_title_suffix ",\\: T\\," \
--leg_pos 0.15 1 \
--data_label_suffix ",\\: 196 \\,\\mathrm{MeV}" ",\\: 220\\,\\mathrm{MeV} " ",\\: 251\\,\\mathrm{MeV}" ",\\: 352\\,\\mathrm{MeV}" \
--plot_kappa_temp \
--temps 196 220 251 352 \
&


if 'false' ; then

# comparison of different temps in lattice units
./plot_rec_corr_fixFlowBytauT.py --suffix temps_latt_units_0.2 --npoints 100 --ylims 0 10 --xlims 0 19 --title 'HISQ, $ms/5$' --min_tauT 0.01 --plot_in_lattice_units \
--flowradiusBytauT 0.2 0.2 0.2 \
--qcdtype hisq_b8249_zeuthenFlow \
--corr EE \
--conftype s096t36_b0824900_m002022_m01011 s096t32_b0824900_m002022_m01011 s096t20_b0824900_m002022_m01011 \
--model 5 \
&

./plot_rec_corr_fixFlowBytauT.py --suffix temps_latt_units_0.3 --npoints 100 --ylims 0 10 --xlims 0 19 --title 'HISQ, $ms/5$' --min_tauT 0.01 --plot_in_lattice_units \
--flowradiusBytauT 0.3 0.3 0.3 \
--qcdtype hisq_b8249_zeuthenFlow \
--corr EE \
--conftype s096t36_b0824900_m002022_m01011 s096t32_b0824900_m002022_m01011 s096t20_b0824900_m002022_m01011 \
--model 5 --min_tauT_plot 0.05 \
&


# compare different scales
./plot_rec_corr_fixFlowBytauT.py --suffix scale_comp --npoints 100 --ylims 0 6 --xlims 0 0.505 --title 'HISQ, $ms/5$' --min_tauT 0.01 \
--flowradiusBytauT 0.3 0.3 0.3 \
--qcdtype hisq_b8249_zeuthenFlow \
--corr EE \
--conftype s096t20_b0824900_m002022_m01011 \
--deduce_fitparam_files --fitparam_basepath ~/work/correlators_flow/data/merged/hisq_b8249_zeuthenFlow/EE/spf/ \
--fitparam_files 3_a_260_1e-04_d_0.4_frankstein_0.300_nt20_piT_1.0 3_a_260_1e-04_d_0.4_frankstein_0.300_nt20_2piT_1.0 3_a_260_1e-04_d_0.4_frankstein_0.300_nt20_2piT_opt \
--PhiUV_basepath /home/altenkort/work/correlators_flow/data/merged/spf_coupling/ \
--PhiUV_files hisq_nt20_SPF_LO_Nf3_0.352_piT_1.0_smax.dat hisq_nt20_SPF_LO_Nf3_0.352_2piT_1.0_smax.dat hisq_nt20_SPF_LO_Nf3_0.352_2piT_opt_smax.dat \
--model 5 --min_tauT_plot 0.05 \
&

# compare temps at two rel flowtimes
./plot_rec_corr_fixFlowBytauT.py --suffix temp_comp_0.2_0.3 --npoints 100 --ylims 0 10 --xlims 0 0.505 --title 'HISQ, $ms/5$' --min_tauT 0.01 \
--flowradiusBytauT 0.3 0.3 0.3 \
--qcdtype hisq_b8249_zeuthenFlow \
--corr EE \
--conftype s096t20_b0824900_m002022_m01011 \
--deduce_fitparam_files --fitparam_basepath ~/work/correlators_flow/data/merged/hisq_b8249_zeuthenFlow/EE/spf/ \
--fitparam_files 3_a_260_1e-04_d_0.4_frankstein_0.300_nt20_piT_1.0 3_a_260_1e-04_d_0.4_frankstein_0.300_nt20_2piT_1.0 3_a_260_1e-04_d_0.4_frankstein_0.300_nt20_2piT_opt \
--PhiUV_basepath /home/altenkort/work/correlators_flow/data/merged/spf_coupling/ \
--PhiUV_files hisq_nt20_SPF_LO_Nf3_0.352_piT_1.0_smax.dat hisq_nt20_SPF_LO_Nf3_0.352_2piT_1.0_smax.dat hisq_nt20_SPF_LO_Nf3_0.352_2piT_opt_smax.dat \
--model 5 --min_tauT_plot 0.05 \
&

fi

wait