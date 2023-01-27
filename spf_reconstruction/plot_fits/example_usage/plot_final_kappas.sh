#!/usr/bin/env bash

basepath="/work/home/altenkort/work/correlators_flow/data/merged/hisq_ms5_zeuthenFlow/EE/"

suffix="" #_paper

plot_hisq(){

    ../plot_final_kappas.py \
        --input_kappa_files \
        ${basepath}/T196/EE_kappa_T196${suffix}.txt \
        ${basepath}/T220/EE_kappa_T220${suffix}.txt \
        ${basepath}/T251/EE_kappa_T251${suffix}.txt \
        ${basepath}/T296/EE_kappa_T296${suffix}.txt \
        \
        --labels \
        196 \
        220 \
        251 \
        296 \
        --outputpath /work/home/altenkort/work/correlators_flow/plots/hisq_ms5_zeuthenFlow/EE/ \
        --suffix "hisq" \
        --temps_in_GeV \
        0.196 \
        0.220 \
        0.251 \
        0.296 \
        --Tc_in_GeV 0.180 \
        --leg_hide \
        --xlims 0.9 2.5 \
        --ylims 0 16.5 \
        --plot_analytical_results --corr EE

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