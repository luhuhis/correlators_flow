#!/bin/bash


qcdtype=$1
corr=$2
if [ -z "$qcdtype" ] ; then
    echo "Usage: $0 qcdtype corr"
    echo "choices for qcdtype: quenched_1.50Tc_zeuthenFlow hisq_ms5_zeuthenFlow"
    echo "choices for corr: EE BB EE_clover BB_clover"
    exit
fi


if [ "$qcdtype" == hisq_ms5_zeuthenFlow ] ; then
    arr_conftypes=(
    "s064t20_b0803600" "s064t22_b0814700" "s096t24_b0824900_m002022_m01011"
    "s064t20_b0785700" "s064t24_b0806800" "s096t28_b0824900_m002022_m01011"
    "s064t20_b0770400" "s064t24_b0791300" "s096t32_b0824900_m002022_m01011"
    "s064t20_b0757000" "s064t24_b0777700" "s096t36_b0824900_m002022_m01011"
)
fi

#    min_trajs=(
#        "0 0 0 0" "0 0 0 0 0 0" "0 0 0"
#        "0 0 0 0" "200 0 130 100 0 0" "175 175 175 175"
#        "0 0 0 0 0 0" "700 400 400 400" "150 150"
#        "0 0 0 0 0 0" "200 200 200 200" "0 0 0 0"  # TODO remove <2000 data files because of gaps for nt=36
#    )
min_trajs=(
        "0 0 0 0" "0 0 0 0 0 0" "0 0 0"
        "0 0 0 0" "0 0 0 0 0 0" "0 0 0 0"
        "0 0 0 0 0 0" "0 0 0 0" "0 0"
        "0 0 0 0 0 0" "0 0 0 0" "0 0 0 0"  # TODO remove <2000 data files because of gaps for nt=36
    )


for idx in "${!arr_conftypes[@]}"; do
    conftype=${arr_conftypes[idx]}
    min_conf="--min_conf ${min_trajs[idx]}"
    ../_1_xestimate_autocorrelations.py --suffix _test2 $min_conf --qcdtype $qcdtype --conftype $conftype --corr $corr --basepath ../../../data/merged/ --basepath_plot ../../../plots/ &
done
wait