#!/bin/bash

plot_quenched(){

    basepath="/work/home/altenkort/work/correlators_flow/data/merged/UV_spfs"

    labels_and_models=(
#        "EE_LO_piT_1" "${basepath}/SPF_LO_Nf0_0.472_piT_1_hmax.npy"
        '$N_f=0, T=1.5T_c$' "dummy"
        '$\mu=\mathrm{max}(\mu_\mathrm{eff}, 1)$' "${basepath}/SPF_LO_Nf0_0.472_eff_1_hmax.npy"
#        "LO_2piT_1" "${basepath}/SPF_LO_Nf0_0.472_2piT_1_hmax.npy"
#        "NLO_piT_1" "${basepath}/SPF_NLO_Nf0_0.472_piT_1_hmax.npy"
#        "NLO_opt_1" "${basepath}/SPF_NLO_Nf0_0.472_eff_1_hmax.npy"
#        "NLO_2piT_1" "${basepath}/SPF_NLO_Nf0_0.472_2piT_1_hmax.npy"
#        "LO_piT_opt" "${basepath}/SPF_LO_Nf0_0.472_piT_opt_hmax.npy"
#        "LO_opt_opt" "${basepath}/SPF_LO_Nf0_0.472_eff_opt_hmax.npy"
#        "LO_2piT_opt" "${basepath}/SPF_LO_Nf0_0.472_2piT_opt_hmax.npy"
#        "EE_NLO_piT_opt" "${basepath}/SPF_NLO_Nf0_0.472_piT_opt_hmax.npy"
#        "EE_LO_eff_opt" "${basepath}/SPF_LO_Nf0_0.472_eff_opt_hmax.npy"
        '$\mu=\mathrm{max}(\mu_\mathrm{eff}, \mu_\mathrm{opt})$' "${basepath}/SPF_NLO_Nf0_0.472_eff_opt_hmax.npy"
        '$\mu=\mathrm{max}(\mu_\mathrm{eff}, \mu_\mathrm{opt,B})$' "${basepath}/SPF_LO_Nf0_0.472_eff_optBB_hmax.npy"
#        "BB_(N)LO_piT" "${basepath}/SPF_LO_Nf0_0.472_piT_optBB_hmax.npy"
#        "BB_(N)LO_piT" "${basepath}/SPF_LO_Nf0_0.472_piT_optBB_hmax.npy"
#        "BB_(N)LO_eff_piT" "${basepath}/SPF_LO_Nf0_0.472_eff_optBBpiT_hmax.npy"
#        "BB_(N)LO_piT_piT" "${basepath}/SPF_LO_Nf0_0.472_piT_optBBpiT_hmax.npy"
#        "NLO_2piT_opt" "${basepath}/SPF_NLO_Nf0_0.472_2piT_opt_hmax.npy"
        ' ' "dummy"
        '$N_f=3, T=1.6T_c$' "dummy"
        '$\mu=\mathrm{max}(\mu_\mathrm{eff}, 1)$' "${basepath}/SPF_LO_Nf3_0.251_eff_1_hmax.npy"
        '$\mu=\mathrm{max}(\mu_\mathrm{eff}, \mu_\mathrm{opt})$' "${basepath}/SPF_NLO_Nf3_0.251_eff_opt_hmax.npy"
    )

    models=()
    labels=()
    for i in "${!labels_and_models[@]}"; do
        if [ $((i % 2)) -eq 0 ]; then
            labels+=("${labels_and_models[i]}")
        else
            models+=("${labels_and_models[i]}")
        fi
    done

    ../plot_UV_spf.py \
        --xlims 0.3 900 \
        --ylims 0 5.5 \
        --PhiUV_files "${models[@]}" \
        --labels "${labels[@]}" \
        --outputpath /work/home/altenkort/work/correlators_flow/plots/quenched_1.50Tc_zeuthenFlow/coupling/ \
        --suffix _EE_quenched_1.5Tc \
        --leg_pos 1 1 \
        --leg_loc "upper right" \
        --fmt
}

plot_hisq(){

    basepath="/work/home/altenkort/work/correlators_flow/data/merged/UV_spfs"

    labels_and_models=(
#        "LO_piT_1" "${basepath}/SPF_LO_Nf3_0.251_piT_1_hmax.npy"
        "LO_eff_1" "${basepath}/SPF_LO_Nf3_0.251_eff_1_hmax.npy"
#        "LO_2piT_1" "${basepath}/SPF_LO_Nf3_0.251_2piT_1_hmax.npy"
#        "NLO_piT_1" "${basepath}/SPF_NLO_Nf3_0.251_piT_1_hmax.npy"
#        "NLO_opt_1" "${basepath}/SPF_NLO_Nf3_0.251_eff_1_hmax.npy"
#        "NLO_2piT_1" "${basepath}/SPF_NLO_Nf3_0.251_2piT_1_hmax.npy"
#        "LO_piT_opt" "${basepath}/SPF_LO_Nf3_0.251_piT_opt_hmax.npy"
#        "LO_opt_opt" "${basepath}/SPF_LO_Nf3_0.251_eff_opt_hmax.npy"
#        "LO_2piT_opt" "${basepath}/SPF_LO_Nf3_0.251_2piT_opt_hmax.npy"
#        "NLO_piT_opt" "${basepath}/SPF_NLO_Nf3_0.251_piT_opt_hmax.npy"
        "NLO_eff_opt" "${basepath}/SPF_NLO_Nf3_0.251_eff_opt_hmax.npy"
#        "NLO_2piT_opt" "${basepath}/SPF_NLO_Nf3_0.251_2piT_opt_hmax.npy"

    )

    models=()
    labels=()
    for i in "${!labels_and_models[@]}"; do
        if [ $((i % 2)) -eq 0 ]; then
            labels+=("${labels_and_models[i]}")
        else
            models+=("${labels_and_models[i]}")
        fi
    done

    ../plot_UV_spf.py \
        --xlims 0.1 1000 \
        --ylims 0 6 \
        --T_in_GeV 0.251 --Nf 3 \
        --PhiUV_files "${models[@]}" \
        --labels "${labels[@]}" \
        --outputpath /work/home/altenkort/work/correlators_flow/plots/quenched_1.50Tc_zeuthenFlow/coupling/ \
        --suffix _EE_hisq_1.6Tc \
        --leg_pos 1 0.5 \
        --leg_loc "center right"

}

plot_quenched
#plot_hisq &
wait