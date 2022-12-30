#!/bin/bash


# EE quenched 1.5 Tc

finalsamples=10000
highsamples=1000
lowsamples=100

labels_and_models=(
'$\text{max}_\text{NLO}$'                   "max_NLO_Nf0_T0.472_minpiT_wopt_${finalsamples}smpls_tauTgtr0.24_final"
'$\text{max}_\text{LO\hphantom{N}}$'        "max_LO_Nf0_T0.472_minpiT_w1_${finalsamples}smpls_tauTgtr0.24_final"
'$\text{smax}_\text{NLO}$'                  "smax_NLO_Nf0_T0.472_minpiT_wopt_${finalsamples}smpls_tauTgtr0.24_final"
'$\text{smax}_\text{LO\hphantom{N}}$'       "smax_LO_Nf0_T0.472_minpiT_w1_${finalsamples}smpls_tauTgtr0.24_final"
'$\text{plaw}_\text{NLO}$'                  "plaw_wIR1.0_wUV3.14_NLO_Nf0_T0.472_minpiT_wopt_${highsamples}smpls_tauTgtr0.24_final"
'$\text{plaw}_\text{LO\hphantom{N}}$'       "plaw_wIR1.0_wUV3.14_LO_Nf0_T0.472_minpiT_w1_${highsamples}smpls_tauTgtr0.24_final"
'$\text{trig1}_\text{NLO}$'                 "fourier_NLO_s1_beta_1_Nf0_T0.472_minpiT_wopt_${finalsamples}smpls_tauTgtr0.24_final"
'$\text{trig1}_\text{LO\hphantom{N}}$'      "fourier_LO_s1_alpha_1_Nf0_T0.472_minpiT_w1_${finalsamples}smpls_tauTgtr0.24_final"
'$\text{trig2}_\text{NLO}$'                 "fourier_NLO_s1_beta_2_Nf0_T0.472_minpiT_wopt_${finalsamples}smpls_tauTgtr0.24_final"
'$\text{trig2}_\text{LO\hphantom{N}}$'      "fourier_LO_s1_alpha_2_Nf0_T0.472_minpiT_w1_${highsamples}smpls_tauTgtr0.24_final"
#'$\text{trig3}_\text{NLO}$'                 "fourier_NLO_s1_beta_3_Nf0_T0.472_minpiT_wopt_${highsamples}smpls_tauTgtr0.24_final"
#'$\text{trig3}_\text{LO\hphantom{N}}$'      "fourier_LO_s1_alpha_3_Nf0_T0.472_minpiT_w1_${highsamples}smpls_tauTgtr0.24_final"
#'$\text{trig4}_\text{NLO}$'                 "fourier_NLO_s1_beta_4_Nf0_T0.472_minpiT_wopt_${lowsamples}smpls_tauTgtr0.24_final"
#'$\text{trig4}_\text{LO\hphantom{N}}$'      "fourier_LO_s1_alpha_4_Nf0_T0.472_minpiT_w1_${highsamples}smpls_tauTgtr0.24_final"
#'$\text{trig5}_\text{NLO}$'                 "fourier_NLO_s1_beta_5_Nf0_T0.472_minpiT_wopt_${lowsamples}smpls_tauTgtr0.24_final"
#'$\text{trig5}_\text{LO\hphantom{N}}$'      "fourier_LO_s1_alpha_5_Nf0_T0.472_minpiT_w1_${highsamples}smpls_tauTgtr0.24_final"
)

models=()
labels=()
for i in "${!labels_and_models[@]}" ; do
    if [ $((i % 2)) -eq 0 ] ; then
        labels+=("${labels_and_models[i]}")
    else
        models+=("${labels_and_models[i]}")
    fi
done


plot_kappa() {
    ../plot_kappa.py \
--xlims -0.15 4.15 \
--xticks 0 1 2 3 4 \
--model_ids "${models[@]}" \
--labels "${labels[@]}" \
--basepath ~/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/EE/spf/ \
--outputpath ~/work/correlators_flow/plots/quenched_1.50Tc_zeuthenFlow/EE/ \
--suffix _quenched_1.5Tc \
--corr EE \
--scale_error_by_chisqdof \
--figsize 7 5
}

plot_fitcorr() {
    ../plot_fitcorr.py \
--model_ids "${models[@]}" \
--labels "${labels[@]}" \
--basepath ~/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/EE/spf/ \
--outputpath ~/work/correlators_flow/plots/quenched_1.50Tc_zeuthenFlow/EE/ \
--suffix _quenched_1.5Tc \
--corr EE
}

plot_spfs(){
    ../plot_spfs.py \
--model_ids "${models[@]}" \
--labels "${labels[@]}" \
--basepath ~/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/EE/spf/ \
--outputpath ~/work/correlators_flow/plots/quenched_1.50Tc_zeuthenFlow/EE/ \
--suffix _quenched_1.5Tc \
--corr EE \
--plot_spf_err \
--ylims 1 150
}

plot_spfs
plot_kappa
plot_fitcorr



