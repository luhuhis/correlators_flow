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
    arr_conftypes=("s064t16_b0687361" "s080t20_b0703500" "s096t24_b0719200" "s120t30_b0739400" "s144t36_b0754400")
    if [ "$corr" == "EE" ] ; then
        acc_sts="--acc_sts acc0.000010_sts0.000010"
        add_args="--legacy --basepath ../../../../data/raw/"
    elif [ "$corr" == "BB" ] ; then
        acc_sts="--acc_sts sts0.150000"
        add_args="--basepath ../../../../data/raw/"
    fi

elif [ "$qcdtype" == hisq_ms5_zeuthenFlow ] ; then
    arr_conftypes=(
    "s064t20_b0803600" "s064t22_b0814700" "s096t24_b0824900_m002022_m01011"
    "s064t20_b0785700" "s064t24_b0806800" "s096t28_b0824900_m002022_m01011"
    "s064t20_b0770400" "s064t24_b0791300" "s096t32_b0824900_m002022_m01011"
    "s064t20_b0757000" "s064t24_b0777700" "s096t36_b0824900_m002022_m01011"
)
    temps=(
        "296"
        "251"
        "220"
        "196"
        )
    add_args="--basepath /work/data/altenkort/gradientFlow"
    flowradiusbasepath="/home/altenkort/work/correlators_flow/data/merged/hisq_ms5_zeuthenFlow/EE/"
    flowradii_refs=()
    for temp in "${temps[@]}"; do
        flowradii_refs+=("--reference_flowradii $flowradiusbasepath/T${temp}/flowradii_T${temp}.dat")
    done
fi


for idx in "${!arr_conftypes[@]}" ; do
    conftype="${arr_conftypes[idx]}"

    if [ "$qcdtype" == hisq_ms5_zeuthenFlow ] ; then
        if [ "$conftype" == "s096t24_b0824900_m002022_m01011" ] ; then
            acc_sts=""
        else
            acc_sts="--acc_sts sts0.150000"
        fi
        if [ "$conftype" == "s096t36_b0824900_m002022_m01011" ] ; then
            even_more_args="--min_conf_nr 2000"
        else
            even_more_args=""
        fi
    fi

    ../_1_merge_data.py \
    ${flowradii_refs[idx/3]} \
    --output_basepath ../../../../data/merged/ \
    $even_more_args $add_args \
    --qcdtype $qcdtype --corr $corr $acc_sts --conftype $conftype
done
