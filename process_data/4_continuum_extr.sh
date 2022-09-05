#!/bin/bash

#conftypes="s080t20_b0703500 s096t24_b0719200 s120t30_b0739400 s144t36_b0754400"

arr_conftypes=(\
    "s064t20_b0785700 s096t28_b0824900_m002022_m01011"\
    "s064t20_b0803600 s096t24_b0824900_m002022_m01011"
    "s064t24_b0777700 s096t36_b0824900_m002022_m01011"
    "s064t24_b0791300 s096t32_b0824900_m002022_m01011")

arr_output_suffix=("T251" "T296" "T196" "T220")

#parameters
max_flow_idx=100 # 170 #135
min_flow_idx=20 #50
qcdtype=hisq_ms5_zeuthenFlow  # quenched_1.50Tc_zeuthenFlow

ylims="0 12" # $"1.5 4"
corrs="EE" # BB_clover EE BB
int="" #--int_only --int_Nt 64"
tex="" #"--use_tex"
max_FlowradiusBytauT=""  # "--max_FlowradiusBytauT 0.5"

for idx in "${!arr_conftypes[@]}" ; do

    for corr in $corrs; do
        rm -f /home/altenkort/work/correlators_flow/plots/$qcdtype/$corr/${arr_output_suffix[idx]}/cont_extr*/*.pdf
        rm -f /home/altenkort/work/correlators_flow//data/merged/$qcdtype/$corr/${arr_output_suffix[idx]}/cont_extr*/*.txt
    done

    args="--output_suffix ${arr_output_suffix[idx]} --qcdtype $qcdtype --conftypes ${arr_conftypes[idx]@Q} --corr $corr $tex --custom_ylims $ylims $max_FlowradiusBytauT $int"

    for ((i=min_flow_idx; i < max_flow_idx; i+=10)) ; do
      for corr in $corrs; do
        python3 -u _4_continuum_extr.py $args --flow_index $((i))   &
        python3 -u _4_continuum_extr.py $args --flow_index $((i+1)) &
        python3 -u _4_continuum_extr.py $args --flow_index $((i+2)) &
        python3 -u _4_continuum_extr.py $args --flow_index $((i+3)) &
        python3 -u _4_continuum_extr.py $args --flow_index $((i+4)) &
        python3 -u _4_continuum_extr.py $args --flow_index $((i+5)) &
        python3 -u _4_continuum_extr.py $args --flow_index $((i+6)) &
        python3 -u _4_continuum_extr.py $args --flow_index $((i+7)) &
        python3 -u _4_continuum_extr.py $args --flow_index $((i+8)) &
        python3 -u _4_continuum_extr.py $args --flow_index $((i+9)) &
        wait
      done
    done
done