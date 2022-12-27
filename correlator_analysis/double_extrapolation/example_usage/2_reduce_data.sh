#!/bin/bash


qcdtype=$1
corr=$2
if [ -z "$qcdtype" ] || [ -z "$corr" ] ; then
    echo "Usage: $0 qcdtype corr"
    echo "choices for qcdtype: quenched_1.50Tc_zeuthenFlow hisq_ms5_zeuthenFlow"
    echo "choices for corr: EE BB EE_clover BB_clover"
    exit
fi

if [ "$qcdtype" == quenched_1.50Tc_zeuthenFlow ] ; then
    arr_conftypes=("s080t20_b0703500" "s096t24_b0719200" "s120t30_b0739400" "s144t36_b0754400")
        min_trajs=(
        "0 0 0 0 0 0 0 0 0 0"
        "0 0 0 0 0 0 0 0 0 0"
        "0 0 0 0 0 0 0 0 0 0"
        "0 0 0 0 0 0 0 0 0 0"
    )
    MC_stepsize=500
    add_args="--skip_binning --already_equally_spaced"
elif  [ "$qcdtype" == hisq_ms5_zeuthenFlow ] ; then
    arr_conftypes=(
    "s064t20_b0803600" "s064t22_b0814700" "s096t24_b0824900_m002022_m01011"
    "s064t20_b0785700" "s064t24_b0806800" "s096t28_b0824900_m002022_m01011"
    "s064t20_b0770400" "s064t24_b0791300" "s096t32_b0824900_m002022_m01011"
    "s064t20_b0757000" "s064t24_b0777700" "s096t36_b0824900_m002022_m01011"
    )
#    add_args="--include_bias"  # only do this once to get a sense for thermalization, then adjust manually.
    min_trajs=(
        "0 100 100 100"                 "  0   0   0   0 0 0 0 0 0 0 0 0" "  0   0   0"
        "0   0   0   0"                 "  0   0   0   0 0 0 0 0 0 0 0 0" "200 200 200 200"
        "0   0   0   0 0 0 0 0 0 0 0 0" "700 100 400 400"                 "175 175"
        "0   0   0   0 0 0 0 0 0 0 0 0" "200 300 100 100"                 "  0   0   0   0"
    )
    MC_stepsize=10
    # TODO switch back to using these, simply because bias estimate is not reliable.
#    min_trajs=(
#        "0 0 0 0" "0 0 0 0 0 0" "0 0 0"
#        "0 0 0 0" "200 0 130 100 0 0" "175 175 175 175"
#        "0 0 0 0 0 0" "700 400 400 400" "150 150"
#        "0 0 0 0 0 0" "200 200 200 200" "0 0 0 0"
#    )
fi




for idx in "${!arr_conftypes[@]}"; do
    conftype=${arr_conftypes[idx]}
    min_conf="--min_conf ${min_trajs[idx]}"
    ../_2_reduce_data.py $add_args --MC_stepsize $MC_stepsize $min_conf --qcdtype $qcdtype --conftype $conftype --corr $corr --basepath ../../../../data/merged/ --basepath_plot ../../../../plots/
done
wait