#!/bin/bash
#without any fits

#s120t30_b0739400 s096t24_b0719200

#./plot_rec_corr_fixFlowBytauT.py  \
#--ylims 2.5 3.75 \
#--output_suffix _quenched_corrs \
#--conftype s144t36_b0754400 s080t20_b0703500 \
#--xlims 0.2 0.512 \
#--plot_quenched_extr \
#--flowradiusBytauT 0.25 \
#--qcdtype quenched_1.50Tc_zeuthenFlow \
#--corr EE \
#--output_path ~/work/correlators_flow/plots/quenched_1.50Tc_zeuthenFlow/EE/ \
#--color_data k C0 C1 C2 C3 C4 C5 \
#--no_kappa_plot \
#--leg_pos 0.15 1


./plot_rec_corr_fixFlowBytauT.py  \
--output_suffix _quenched_corrs_withFits \
--conftype s144t36_b0754400 s080t20_b0703500 \
--xlims 0.23 0.512 \
--ylims 2.6 3.8 \
--kappa_ylims 0 5 \
--kappa_xlims 0 4 \
--kappa_labelheight 0.57 \
--leg_pos 0.75 0.65 \
--kappa_ypos 2.65 2.2 1.93 1.65 \
--plot_quenched_extr --plot_quenched_systematics \
--flowradiusBytauT 0.20 \
--qcdtype quenched_1.50Tc_zeuthenFlow \
--corr EE \
--output_path ~/work/correlators_flow/plots/quenched_1.50Tc_zeuthenFlow/EE/ \
--color_data k C0 C1 C2 C3 C4 C5 \
--color_fit k C0 C1 C2 \
--deduce_fitparam_files --fitparam_basepath ~/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/EE/ \
--fitparam_files \
/spf/pnorm2.0_LO_T0.472_min2piT_w1_500_0.0_exp0_quenched_cont_f0 \
/cont_rel_flow/rel_flow/spf/pnorm2.0_LO_T0.472_min2piT_w1_500_0.35_exp0_quenched_cont_f0.20 \
/s144t36_b0754400/rel_flow/spf/pnorm2.0_LO_T0.472_min2piT_w1_500_0.35_exp0_quenched_s144t36_b0754400_f0.20 \
/s080t20_b0703500/rel_flow/spf/pnorm2.0_LO_T0.472_min2piT_w1_500_0.35_exp0_quenched_s080t20_b0703500_f0.20 \
--PhiUV_basepath ~/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/EE/ \
--PhiUV_files \
/spf/pnorm2.0_LO_T0.472_min2piT_w1_500_0.0_exp0_quenched_cont_f0/phiUV.dat \
/cont_rel_flow/rel_flow/spf/pnorm2.0_LO_T0.472_min2piT_w1_500_0.35_exp0_quenched_cont_f0.20/phiUV.dat \
/s144t36_b0754400/rel_flow/spf/pnorm2.0_LO_T0.472_min2piT_w1_500_0.35_exp0_quenched_s144t36_b0754400_f0.20/phiUV.dat \
/s080t20_b0703500/rel_flow/spf/pnorm2.0_LO_T0.472_min2piT_w1_500_0.35_exp0_quenched_s080t20_b0703500_f0.20/phiUV.dat \
--min_tauT 0.25 0.35 0.35 0.35 \
--leg_n_col 1 \
--model pnorm --p 2 --min_scale "2\pi T\," --run_scale "\omega" \
--no_just_UV --no_connection --no_label \

#--tauT_vlines 0.35 \

#  --OmegaByT_IR 1.0 --OmegaByT_UV 6.28
# /s144t36_b0754400/rel_flow/spf/pnorm2.0_LO_T0.472_min2piT_w1_500_0.35_exp0_quenched_s144t36_b0754400_f0.25/phiUV.dat \
# /s144t36_b0754400/rel_flow/spf/pnorm2.0_LO_T0.472_min2piT_w1_500_0.35_exp0_quenched_s144t36_b0754400_f0.25 \
#s144t36_b0754400