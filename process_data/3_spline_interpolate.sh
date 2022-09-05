#!/bin/bash

# "s080t20_b0703500 s096t24_b0719200 s120t30_b0739400 s144t36_b0754400" #""
arr_conftypes=("s064t20_b0785700" "s064t20_b0803600" "s064t24_b0777700" "s064t24_b0791300")
arr_int_Nt=(28 24 36 32)

for idx in "${!arr_conftypes[@]}" ; do

    conftypes=${arr_conftypes[idx]}
    int_Nt=${arr_int_Nt[idx]}
    max_flow_idx=100 #170
    min_flow_idx=0
    qcdtype=hisq_ms5_zeuthenFlow  # quenched_1.50Tc_zeuthenFlow
    corrs="EE" # BB BB_clover EE
    nsamples=1000
    max_FlowradiusBytauT=""  # "--max_FlowradiusBytauT 0.33"


    for conftype in $conftypes; do
        for corr in $corrs; do
            rm -f /home/altenkort/work/correlators_flow/plots/$qcdtype/$corr/$conftype/interpolations/*interpolation.pdf
            rm -f /home/altenkort/work/correlators_flow//data/merged/$qcdtype/$corr/$conftype/interpolations/*interpolation.txt
        done
    done

    args="--ylims 0 12 $max_FlowradiusBytauT --qcdtype $qcdtype --conftype $conftype --corr $corr --int_Nt $int_Nt --nsamples $nsamples"

    for corr in $corrs; do
        for conftype in $conftypes ; do
            for ((i=min_flow_idx; i < max_flow_idx; i+=10)) ; do
                python3 -u _3_spline_interpolate.py $args --flow_index $((i))   &
                python3 -u _3_spline_interpolate.py $args --flow_index $((i+1)) &
                python3 -u _3_spline_interpolate.py $args --flow_index $((i+2)) &
                python3 -u _3_spline_interpolate.py $args --flow_index $((i+3)) &
                python3 -u _3_spline_interpolate.py $args --flow_index $((i+4)) &
                python3 -u _3_spline_interpolate.py $args --flow_index $((i+5)) &
                python3 -u _3_spline_interpolate.py $args --flow_index $((i+6)) &
                python3 -u _3_spline_interpolate.py $args --flow_index $((i+7)) &
                python3 -u _3_spline_interpolate.py $args --flow_index $((i+8)) &
                python3 -u _3_spline_interpolate.py $args --flow_index $((i+9)) &
                wait
            done
        done
    done

done