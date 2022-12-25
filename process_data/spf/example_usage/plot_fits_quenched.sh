#!/bin/bash


# EE quenched 1.5 Tc

#nsamples=1000
highsamples=1000
#lowsamples=100

labels_and_models=(
"max NLO"         "max_NLO_Nf0_T0.472_minpiT_wopt_${highsamples}smpls_tauTgtr0.24_final"
"max LO"          "max_LO_Nf0_T0.472_minpiT_w1_${highsamples}smpls_tauTgtr0.24_final"
"smax NLO"        "smax_NLO_Nf0_T0.472_minpiT_wopt_${highsamples}smpls_tauTgtr0.24_final"
"smax LO"         "smax_LO_Nf0_T0.472_minpiT_w1_${highsamples}smpls_tauTgtr0.24_final"
#"line NLO"        "line_wIR1.0_wUV3.14_NLO_Nf0_T0.472_minpiT_wopt_${highsamples}smpls_tauTgtr0.24_final"
#"line LO"         "line_wIR1.0_wUV3.14_LO_Nf0_T0.472_minpiT_w1_${highsamples}smpls_tauTgtr0.24_final"
#"fourier NLO 1"   "fourier_NLO_s1_beta_1_Nf0_T0.472_minpiT_wopt_${highsamples}smpls_tauTgtr0.24_final"
#"fourier LO 1"    "fourier_LO_s1_alpha_1_Nf0_T0.472_minpiT_w1_${highsamples}smpls_tauTgtr0.24_final"
"fourier NLO 2"  "fourier_NLO_s1_beta_1_Nf0_T0.472_minpiT_wopt_${highsamples}smpls_tauTgtr0.24_final"
"fourier LO 2"   "fourier_LO_s1_alpha_1_Nf0_T0.472_minpiT_w1_${highsamples}smpls_tauTgtr0.24_final"
"fourier NLO 3"   "fourier_NLO_s1_beta_2_Nf0_T0.472_minpiT_wopt_${highsamples}smpls_tauTgtr0.24_final"
"fourier LO 3"    "fourier_LO_s1_alpha_2_Nf0_T0.472_minpiT_w1_${highsamples}smpls_tauTgtr0.24_final"
#"fourier NLO 4"   fourier_NLO_s1_beta_3_Nf0_T0.472_minpiT_wopt_${lowsamples}smpls_tauTgtr0.24_final
#"fourier LO 4"    fourier_LO_s1_alpha_3_Nf0_T0.472_minpiT_w1_${lowsamples}smpls_tauTgtr0.24_final
#"fourier NLO 5"   fourier_NLO_s1_beta_4_Nf0_T0.472_minpiT_wopt_${lowsamples}smpls_tauTgtr0.24_final
#"fourier LO 5"    fourier_LO_s1_alpha_4_Nf0_T0.472_minpiT_w1_${lowsamples}smpls_tauTgtr0.24_final
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
--xlims 0 4.2 \
--model_ids "${models[@]}" \
--labels "${labels[@]}" \
--basepath ~/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/EE/spf/ \
--outputpath ~/work/correlators_flow/plots/quenched_1.50Tc_zeuthenFlow/EE/ \
--suffix _quenched_1.5Tc \
--corr EE \
--scale_error_by_chisqdof
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
--ylims 0.8 150
}

plot_spfs
plot_kappa
plot_fitcorr



