#!/usr/bin/env bash

basepath_work_data="${1}"
basepath_plot="${2}"
selector="${3:-"all"}"

plot_hisq_thesis(){

        basepath="$basepath_work_data/hisq_ms5_zeuthenFlow/EE/"

        suffix=""

        ../plot_final_kappas.py \
            --input_kappa_files \
            $basepath_work_data/quenched_1.50Tc_zeuthenFlow/EE/EE_kappa_quenched_1.5Tc.txt \
            ${basepath}/T195/EE_kappa_T195${suffix}.txt \
            ${basepath}/T220/EE_kappa_T220${suffix}.txt \
            ${basepath}/T251/EE_kappa_T251${suffix}.txt \
            ${basepath}/T293/EE_kappa_T293${suffix}.txt \
            --labels \
            "\textbf{this work}* (flow)" '\textbf{this work} (flow)' "" "" "" "" \
            --colors \
            k m m m m \
            --outputpath $basepath_plot/hisq_ms5_zeuthenFlow/EE/ \
            --suffix "hisq_thesis" \
            --temps_in_GeV \
            0.270 0.195 0.220 0.251 0.293 \
            --Tc_in_GeV 0.180 \
            --leg_ncol 1 \
            --xlims 1 3.1 \
            --ylims 0 16.5 \
            --corr EE \
            --plot_EE_quenched_lit --add_leg_titles --xlabelpos 0.93 0.01 \
            --leg_fontsize 9 

}

plot_hisq_paper(){


        addargs="--no_subscript"


        basepath="$basepath_work_data/hisq_ms5_zeuthenFlow/EE/"

        ../plot_final_kappas.py \
            --input_kappa_files \
            ${basepath}/T195/EE_kappa_T195${suffix}.txt \
            ${basepath}/T220/EE_kappa_T220${suffix}.txt \
            ${basepath}/T251/EE_kappa_T251${suffix}.txt \
            ${basepath}/T293/EE_kappa_T293${suffix}.txt \
            ${basepath}/s096t36_b0824900_m002022_m01011/EE_kappa_T195_finiteflow_paper.txt \
            ${basepath}/s096t32_b0824900_m002022_m01011/EE_kappa_T220_finiteflow_paper.txt \
            ${basepath}/s096t28_b0824900_m002022_m01011/EE_kappa_T251_finiteflow_paper.txt \
            ${basepath}/s096t24_b0824900_m002022_m01011/EE_kappa_T293_finiteflow_paper.txt \
            ${basepath}/s096t20_b0824900_m002022_m01011/EE_kappa_T352_finiteflow_paper.txt \
            --labels \
            '$ a\rightarrow 0, \tau_\mathrm{F} \rightarrow 0$' "" "" "" '\begin{flushleft}$ a^{-1}=7.036\,\mathrm{GeV}$, \newline $\sqrt{8\tau_\mathrm{F}}/\tau=0.3$ \end{flushleft}' "" "" "" "" \
            --fillstyles \
            none none none none full full full full full \
            --fmts \
            s s s s o o o o o \
            --zorders \
            1 1 1 1 0 0 0 0 0 \
            --colors \
            C0 C0 C0 C0 C1 C1 C1 C1 C1 \
            --markersize 4 \
            --outputpath $basepath_plot/hisq_ms5_zeuthenFlow/EE/ \
            --suffix "hisq${suffix}" \
            --temps_in_GeV \
            0.195 0.220 0.251 0.293 0.198 0.222 0.253 0.295 0.352 \
            --Tc_in_GeV 0.180 \
            --leg_ncol 1 \
            --xlims 0.9 2.4 \
            --ylims 0 16.5 \
            --corr EE \
            $addargs

}

plot_quenched_EE(){

    basepath="$basepath_work_data/quenched_1.50Tc_zeuthenFlow/EE/"

    ../plot_final_kappas.py \
        --input_kappa_files \
        ${basepath}/EE_kappa_quenched_1.5Tc.txt \
        --labels \
        '\textbf{this work}* (flow)' \
        --outputpath "$basepath_plot/quenched_1.50Tc_zeuthenFlow/EE/" \
        --suffix "EE_quenched_literature" \
        --temps_in_GeV \
        0.4725 \
        --Tc_in_GeV 0.315 \
        --xlims 1 3.1 \
        --ylims 0 7.5 \
        --plot_EE_quenched_lit --corr EE

}

plot_quenched_BB(){

    basepath="$basepath_work_data/quenched_1.50Tc_zeuthenFlow/BB/"

    ../plot_final_kappas.py \
        --input_kappa_files \
        ${basepath}/BB_kappa_quenched_1.5Tc.txt \
        --labels \
        '\textbf{this work} (flow)' \
        --outputpath "$basepath_plot/quenched_1.50Tc_zeuthenFlow/BB/" \
        --suffix "BB_quenched_literature" \
        --temps_in_GeV \
        0.482 \
        --Tc_in_GeV 0.315 \
        --xlims 1 3.1 \
        --ylims 0 3.4 \
        --plot_BB_quenched_lit --corr BB

}

(
    cd "$(dirname $0)" || exit

    if [ "$selector" == "all" ]; then
        plot_hisq_paper &
        plot_hisq_thesis &
        plot_quenched_EE &
        plot_quenched_BB &
    elif [ "$selector" == "hisq" ]; then
        plot_hisq_paper &
    elif [ "$selector" == "hisq_thesis" ]; then
        plot_hisq_thesis &
    elif [ "$selector" == "EE" ]; then
        plot_quenched_EE &
    elif [ "$selector" == "BB" ]; then
        plot_quenched_BB &
    else
        echo "Error: unknown selector"
    fi
    wait
)