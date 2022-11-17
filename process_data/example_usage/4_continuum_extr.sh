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
    corrs="EE" # BB_clover EE BB
    add_args=""
    min_flowradius=0.055
elif [ "$qcdtype" == hisq_ms5_zeuthenFlow ] ; then
    arr_conftypes=(\
    "s064t20_b0770400 s064t24_b0791300 s096t32_b0824900_m002022_m01011"
    "s064t20_b0785700 s064t24_b0806800 s096t28_b0824900_m002022_m01011"
    "s064t20_b0803600 s064t22_b0814700 s096t24_b0824900_m002022_m01011"
    "s064t20_b0757000 s064t24_b0777700 s096t36_b0824900_m002022_m01011")
    arr_output_suffix=("T220" "T251" "T296" "T196")
    ylims="0 12"
    min_flowradius=0.055
    corrs="EE" # BB_clover EE BB
fi

for idx in "${!arr_conftypes[@]}" ; do

#    for corr in $corrs; do
#        rm -f /home/altenkort/work/correlators_flow/plots/$qcdtype/$corr/${arr_output_suffix[idx]}/cont_extr*/*.pdf
#        rm -f /home/altenkort/work/correlators_flow//data/merged/$qcdtype/$corr/${arr_output_suffix[idx]}/cont_extr*/*.txt
#    done

    if [ "${arr_output_suffix[idx]}" ] ; then
        sufargs="--output_suffix"
    fi
    files=""
    for conftype in ${arr_conftypes[idx]} ; do
        files="$files $conftype/flowradii_$conftype.dat"
    done

    for corr in $corrs; do
#        ../find_common_flowtimes.py --basepath ../../../data/merged/$qcdtype/$corr/ --files $files --output ../../../data/merged/$qcdtype/$corr/${arr_output_suffix[idx]}/flowradii_${arr_output_suffix[idx]}.dat
        args="$add_args --nproc 35 --min_flowradius $min_flowradius --basepath ../../../data/merged/ --basepath_plot ../../../plots/ --max_FlowradiusBytauT_offset 0 --max_FlowradiusBytauT 0.4  $sufargs ${arr_output_suffix[idx]} --qcdtype $qcdtype --conftypes ${arr_conftypes[idx]} --corr $corr --custom_ylims $ylims"
        ../_4_continuum_extr.py $args
    done
done