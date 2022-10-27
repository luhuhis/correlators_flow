#!/bin/bash

qcdtype=$1

if [ -z "$qcdtype" ] ; then
    echo "Usage: $0 qcdtype"
    echo "choices for qcdtype: quenched_1.50Tc_zeuthenFlow hisq_ms5_zeuthenFlow"
    exit
fi

ncpu=20

if [ "$qcdtype" == quenched_1.50Tc_zeuthenFlow ] ; then
    arr_conftypes=("s080t20_b0703500 s096t24_b0719200 s120t30_b0739400 s144t36_b0754400")
    arr_output_suffix=("" )
    ylims="1.5 4"
    max_flow_idx=170
    min_flow_idx=0
    corrs="EE" # BB_clover EE BB
    add_args=""
elif [ "$qcdtype" == hisq_ms5_zeuthenFlow ] ; then
    arr_conftypes=(\
    "s064t20_b0770400 s064t24_b0791300 s096t32_b0824900_m002022_m01011"
    "s064t20_b0785700 s064t24_b0806800 s096t28_b0824900_m002022_m01011"
    "s064t20_b0803600 s064t22_b0814700 s096t24_b0824900_m002022_m01011"
    "s064t20_b0757000 s064t24_b0777700 s096t36_b0824900_m002022_m01011")
    arr_output_suffix=("T220" "T251" "T296" "T196")
    ylims="0 12"
    max_flow_idx=170
    min_flow_idx=0
    corrs="EE" # BB_clover EE BB
fi

for idx in "${!arr_conftypes[@]}" ; do

    for corr in $corrs; do
        rm -f /home/altenkort/work/correlators_flow/plots/$qcdtype/$corr/${arr_output_suffix[idx]}/cont_extr*/*.pdf
        rm -f /home/altenkort/work/correlators_flow//data/merged/$qcdtype/$corr/${arr_output_suffix[idx]}/cont_extr*/*.txt
    done

    if [ "${arr_output_suffix[idx]}" ] ; then
        sufargs="--output_suffix"
    fi

    args="$add_args --max_FlowradiusBytauT_offset 0 --max_FlowradiusBytauT 0.5  $sufargs ${arr_output_suffix[idx]} --qcdtype $qcdtype --conftypes ${arr_conftypes[idx]} --corr $corr --custom_ylims $ylims"

    for ((i=min_flow_idx; i < max_flow_idx; i+=20)) ; do
        for corr in $corrs; do
            for ((j=0; j < ncpu; j+=1)) ; do
                python3 -u _4_continuum_extr.py $args --flow_index $((i+j))   &
            done
            wait
        done
    done
done