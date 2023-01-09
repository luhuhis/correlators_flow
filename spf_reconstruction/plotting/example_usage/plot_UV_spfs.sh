#!/bin/bash

plot_quenched(){

    basepath="/work/home/altenkort/work/correlators_flow/data/merged/UV_spfs"

    labels_and_models=(
#        "LO_piT_1" "${basepath}/SPF_LO_Nf0_0.472_piT_1_hmax.npy"
        "LO_opt_1" "${basepath}/SPF_LO_Nf0_0.472_opt_1_hmax.npy"
#        "LO_2piT_1" "${basepath}/SPF_LO_Nf0_0.472_2piT_1_hmax.npy"
#        "NLO_piT_1" "${basepath}/SPF_NLO_Nf0_0.472_piT_1_hmax.npy"
#        "NLO_opt_1" "${basepath}/SPF_NLO_Nf0_0.472_opt_1_hmax.npy"
#        "NLO_2piT_1" "${basepath}/SPF_NLO_Nf0_0.472_2piT_1_hmax.npy"
#        "LO_piT_opt" "${basepath}/SPF_LO_Nf0_0.472_piT_opt_hmax.npy"
#        "LO_opt_opt" "${basepath}/SPF_LO_Nf0_0.472_opt_opt_hmax.npy"
#        "LO_2piT_opt" "${basepath}/SPF_LO_Nf0_0.472_2piT_opt_hmax.npy"
#        "NLO_piT_opt" "${basepath}/SPF_NLO_Nf0_0.472_piT_opt_hmax.npy"
        "NLO_opt_opt" "${basepath}/SPF_NLO_Nf0_0.472_opt_opt_hmax.npy"
#        "NLO_2piT_opt" "${basepath}/SPF_NLO_Nf0_0.472_2piT_opt_hmax.npy"
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
        --xlims 0.01 1000 \
        --ylims 0 10 \
        --T_in_GeV 0.472 --Nf 0 \
        --PhiUV_files "${models[@]}" \
        --labels "${labels[@]}" \
        --outputpath /work/home/altenkort/work/correlators_flow/plots/quenched_1.50Tc_zeuthenFlow/ \
        --suffix _EE_quenched_1.5Tc \
        --leg_pos 1 0.5 \
        --leg_loc "center right"
}

plot_hisq(){

    basepath="/work/home/altenkort/work/correlators_flow/data/merged/UV_spfs"

    labels_and_models=(
#        "LO_piT_1" "${basepath}/SPF_LO_Nf3_0.251_piT_1_hmax.npy"
        "LO_opt_1" "${basepath}/SPF_LO_Nf3_0.251_opt_1_hmax.npy"
#        "LO_2piT_1" "${basepath}/SPF_LO_Nf3_0.251_2piT_1_hmax.npy"
#        "NLO_piT_1" "${basepath}/SPF_NLO_Nf3_0.251_piT_1_hmax.npy"
#        "NLO_opt_1" "${basepath}/SPF_NLO_Nf3_0.251_opt_1_hmax.npy"
#        "NLO_2piT_1" "${basepath}/SPF_NLO_Nf3_0.251_2piT_1_hmax.npy"
#        "LO_piT_opt" "${basepath}/SPF_LO_Nf3_0.251_piT_opt_hmax.npy"
#        "LO_opt_opt" "${basepath}/SPF_LO_Nf3_0.251_opt_opt_hmax.npy"
#        "LO_2piT_opt" "${basepath}/SPF_LO_Nf3_0.251_2piT_opt_hmax.npy"
#        "NLO_piT_opt" "${basepath}/SPF_NLO_Nf3_0.251_piT_opt_hmax.npy"
        "NLO_opt_opt" "${basepath}/SPF_NLO_Nf3_0.251_opt_opt_hmax.npy"
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
        --xlims 0.01 1000 \
        --ylims 0 10 \
        --T_in_GeV 0.251 --Nf 3 \
        --PhiUV_files "${models[@]}" \
        --labels "${labels[@]}" \
        --outputpath /work/home/altenkort/work/correlators_flow/plots/quenched_1.50Tc_zeuthenFlow/ \
        --suffix _EE_hisq_1.6Tc \
        --leg_pos 1 0.5 \
        --leg_loc "center right"

}

plot_quenched &
plot_hisq &
wait