#!/usr/bin/env bash



plot_hisq(){

    for suffix in "" "_paper" ; do

        if [ "${suffix}" == "_paper" ] ; then
            addargs="--no_subscript"
        else
            addargs="--plot_analytical_results"
        fi

        basepath="/work/home/altenkort/work/correlators_flow/data/merged/hisq_ms5_zeuthenFlow/EE/"

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
            '$ a\rightarrow 0, \tau_\mathrm{F} \rightarrow 0$' "" "" "" '\begin{flushleft}$ a^{-1}=7.06\,\mathrm{GeV}$, \newline $\sqrt{8\tau_\mathrm{F}}/\tau=0.3$ \end{flushleft}' "" "" "" "" \
            --fillstyles \
            none none none none full full full full full \
            --fmts \
            s s s s o o o o o \
            --zorders \
            1 1 1 1 0 0 0 0 0 \
            --colors \
            C0 C0 C0 C0 C1 C1 C1 C1 C1 \
            --markersize 4 \
            --outputpath /work/home/altenkort/work/correlators_flow/plots/hisq_ms5_zeuthenFlow/EE/ \
            --suffix "hisq${suffix}" \
            --temps_in_GeV \
            0.195 0.220 0.251 0.293 0.198 0.222 0.253 0.295 0.352 \
            --Tc_in_GeV 0.180 \
            --leg_ncol 1 \
            --xlims 0.9 2.4 \
            --ylims 0 16.5 \
            --corr EE \
            $addargs

    done

}

plot_quenched_EE(){

    basepath="/work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/EE/"

    ../plot_final_kappas.py \
        --input_kappa_files \
        ${basepath}/EE_kappa_quenched_1.5Tc.txt \
        --labels \
        '\textbf{this work} (flow)' \
        --outputpath "/work/home/altenkort/work/correlators_flow/plots/quenched_1.50Tc_zeuthenFlow/EE/" \
        --suffix "EE_quenched_literature" \
        --temps_in_GeV \
        0.493 \
        --Tc_in_GeV 0.315 \
        --xlims 0.9 3.1 \
        --ylims 0 7.5 \
        --plot_EE_quenched_lit --corr EE

}

plot_quenched_BB(){

    basepath="/work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/BB/"

    ../plot_final_kappas.py \
        --input_kappa_files \
        ${basepath}/BB_kappa_quenched_1.5Tc.txt \
        --labels \
        '\textbf{this work} (flow)' \
        --outputpath "/work/home/altenkort/work/correlators_flow/plots/quenched_1.50Tc_zeuthenFlow/BB/" \
        --suffix "BB_quenched_literature" \
        --temps_in_GeV \
        0.482 \
        --Tc_in_GeV 0.315 \
        --xlims 0.9 3.1 \
        --ylims 0 4.5 \
        --plot_BB_quenched_lit --corr BB

}

plot_hisq &
plot_quenched_EE &
plot_quenched_BB &

wait