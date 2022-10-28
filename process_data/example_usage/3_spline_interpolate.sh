#!/bin/bash

qcdtype=$1

if [ -z "$qcdtype" ] ; then
    echo "Usage: $0 qcdtype"
    echo "choices for qcdtype: quenched_1.50Tc_zeuthenFlow hisq_ms5_zeuthenFlow"
    exit
fi

ncpu=20
nsamples=1000
min_flow_idx=0

if [ "$qcdtype" == quenched_1.50Tc_zeuthenFlow ] ; then
    arr_conftypes=("s080t20_b0703500" "s096t24_b0719200" "s120t30_b0739400")
    arr_int_Nt=(36 36 36)
    max_flow_idx=170 #100
    corrs="EE"
    ylims="0 5"
    add_args="--max_FlowradiusBytauT 0.5 --max_FlowradiusBytauT_offset=0"
elif [ "$qcdtype" == hisq_ms5_zeuthenFlow ] ; then
    arr_conftypes=( "s064t20_b0785700" "s064t24_b0806800" "s064t20_b0803600" "s064t22_b0814700" "s064t20_b0770400" "s064t24_b0791300" "s064t20_b0757000" "s064t24_b0777700")
    arr_int_Nt=(28 28 24 24 32 32 36 36)
    max_flow_idx=100
    corrs="EE"
    ylims="0 12"
    add_args="--max_FlowradiusBytauT 0.5 --max_FlowradiusBytauT_offset=0"
fi

for idx in "${!arr_conftypes[@]}"; do
    conftype="${arr_conftypes[idx]}"
    for corr in "${corrs[@]}"; do
        rm -f /home/altenkort/work/correlators_flow/plots/$qcdtype/$corr/$conftype/interpolations/*interpolation.pdf
        rm -f /home/altenkort/work/correlators_flow//data/merged/$qcdtype/$corr/$conftype/interpolations/*interpolation.txt

        args="$add_args --ylims $ylims --qcdtype $qcdtype --conftype $conftype --corr $corr --int_Nt ${arr_int_Nt[idx]} --nsamples $nsamples"

        ../_3_spline_interpolate.py "$args"

    done
done

