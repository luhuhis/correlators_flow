#!/bin/bash

qcdtype=$1
if [ -z "$qcdtype" ] ; then
    echo "Usage: $0 qcdtype"
    echo "choices for qcdtype: quenched_1.50Tc_zeuthenFlow hisq_ms5_zeuthenFlow"
    exit
fi

if [ "$qcdtype" == quenched_1.50Tc_zeuthenFlow ] ; then
    arr_conftypes=("s080t20_b0703500" "s096t24_b0719200" "s120t30_b0739400" "s144t36_b0754400")
elif [ "$qcdtype" == hisq_ms5_zeuthenFlow ] ; then
    arr_conftypes=(
    "s064t20_b0803600" "s064t22_b0814700" "s096t24_b0824900_m002022_m01011"
    "s064t20_b0785700" "s064t24_b0806800" "s096t28_b0824900_m002022_m01011"
    "s064t20_b0770400" "s064t24_b0791300" "s096t32_b0824900_m002022_m01011"
    "s064t20_b0757000" "s064t24_b0777700" "s096t36_b0824900_m002022_m01011"
)
    arr_discards=(
        "0 0 0 0" "0 0 0 0 0 0" "0 0 0"
        "0 0 0 0" "0 0 0 0 0 0" "400 400 400 400"
        "0 0 0 0 0 0" "0 200 200 200" "400 400"
        "0 0 0 0 0 0" "0 100 100 100" "400 400 400 400"
    )
fi


for idx in "${!arr_conftypes[@]}" ; do
    conftype=${arr_conftypes[idx]}
    if [ "$qcdtype" == hisq_ms5_zeuthenFlow ] ; then
        discard="--n_discard ${arr_discards[idx]}"
    elif [ "$qcdtype" == quenched_1.50Tc_zeuthenFlow ] ; then
        discard="--n_discard 0 0 0 0 0 0 0 0 0 0"
    fi
        ../_2_reduce_data.py  $discard --qcdtype "$qcdtype" --corr EE --conftype "$conftype" --basepath ../../../data/merged/ --nproc 5 &

done
wait