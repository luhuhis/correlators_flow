#!/bin/bash

basepath=$1
outputpath=$2

if [ -z "$basepath" ] || [ -z "$outputpath" ] ; then
    echo "Usage: $0 basepath outputpath"
    exit
fi

basepath="${basepath}/UV_spfs/"


plot_quenched(){

    labels_and_models=(
#        "EE_LO_piT_1" "${basepath}/SPF_LO_Nf0_0.472_piT_1_hmax.npy"
        '$N_f=0,\ T=1.5\,T_c$' "dummy"
        '$\mu=\mathrm{max}(\mu_\mathrm{DR}, \omega)$' "${basepath}/SPF_EE_LO_Nf0_0.472_eff_1_hmax.npy"
#        "LO_2piT_1" "${basepath}/SPF_LO_Nf0_0.472_2piT_1_hmax.npy"
#        "NLO_piT_1" "${basepath}/SPF_NLO_Nf0_0.472_piT_1_hmax.npy"
#        "NLO_opt_1" "${basepath}/SPF_NLO_Nf0_0.472_eff_1_hmax.npy"
#        "NLO_2piT_1" "${basepath}/SPF_NLO_Nf0_0.472_2piT_1_hmax.npy"
#        "LO_piT_opt" "${basepath}/SPF_LO_Nf0_0.472_piT_opt_hmax.npy"
#        "LO_opt_opt" "${basepath}/SPF_LO_Nf0_0.472_eff_opt_hmax.npy"
#        "LO_2piT_opt" "${basepath}/SPF_LO_Nf0_0.472_2piT_opt_hmax.npy"
#        "EE_NLO_piT_opt" "${basepath}/SPF_NLO_Nf0_0.472_piT_opt_hmax.npy"
#        "EE_LO_eff_opt" "${basepath}/SPF_LO_Nf0_0.472_eff_opt_hmax.npy"
        '$\mu=\mathrm{max}(\mu_\mathrm{DR}, \mu_\mathrm{opt})$' "${basepath}/SPF_EE_LO_Nf0_0.472_eff_opt_hmax.npy"

#        '$\mu=\mathrm{max}(\mu_\mathrm{eff}, \mu_\mathrm{opt,B})$' "${basepath}/SPF_LO_Nf0_0.472_eff_optBB_hmax.npy"

#        "BB_(N)LO_piT" "${basepath}/SPF_LO_Nf0_0.472_piT_optBB_hmax.npy"
#        "BB_(N)LO_piT" "${basepath}/SPF_LO_Nf0_0.472_piT_optBB_hmax.npy"
#        "BB_(N)LO_eff_piT" "${basepath}/SPF_LO_Nf0_0.472_eff_optBBpiT_hmax.npy"
#        "BB_(N)LO_piT_piT" "${basepath}/SPF_LO_Nf0_0.472_piT_optBBpiT_hmax.npy"
#        "NLO_2piT_opt" "${basepath}/SPF_NLO_Nf0_0.472_2piT_opt_hmax.npy"
        '\begin{center}$N_f=3,\ T\approx 1.4\, T_c,$ \\[-0.5ex] $m_l=m_s/5$\end{center}' "dummy"
        '$\mu=\mathrm{max}(\mu_\mathrm{DR}, \omega)$' "${basepath}/SPF_EE_LO_Nf3_0.251_eff_1_hmax.npy"
        '$\mu=\mathrm{max}(\mu_\mathrm{DR}, \mu_\mathrm{opt})$' "${basepath}/SPF_EE_LO_Nf3_0.251_eff_opt_hmax.npy"
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

    (
        cd "$(dirname "$0")" || exit
        ../plot_g2.py \
            --xlims 0.3 900 \
            --ylims 0 4.5 \
            --PhiUV_files "${models[@]}" \
            --labels "${labels[@]}" \
            --outputpath ${outputpath}/quenched_1.50Tc_zeuthenFlow/coupling/ \
            --suffix _EE_quenched_1.5Tc \
            --leg_pos 0 0 \
            --leg_loc "lower left" \
            --fmt
    )
}


plot_quenched
