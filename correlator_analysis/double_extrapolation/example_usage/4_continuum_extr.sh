#!/bin/bash

qcdtype=$1
corr=$2
addsuf=$3
if [ -z "$qcdtype" ] || [ -z "$corr" ] ; then
    echo "Usage: $0 qcdtype corr"
    echo "choices for qcdtype: quenched_1.50Tc_zeuthenFlow hisq_ms5_zeuthenFlow"
    echo "choices for corr: EE BB EE_clover BB_clover"
    exit
fi


ncpu=20

if [ "$qcdtype" == quenched_1.50Tc_zeuthenFlow ] ; then
    arr_conftypes=("s080t20_b0703500 s096t24_b0719200 s120t30_b0739400 s144t36_b0754400")
    arr_output_suffix=("" )
    ylims="2.4 3.9"
    add_args="--relflow --nsamples 1000"
    min_flowradius=0.05
elif [ "$qcdtype" == hisq_ms5_zeuthenFlow ] ; then
    arr_conftypes=(
    "s064t20_b0757000 s064t24_b0777700 s096t36_b0824900_m002022_m01011"
    "s064t20_b0770400 s064t24_b0791300 s096t32_b0824900_m002022_m01011"
    "s064t20_b0785700 s064t24_b0806800 s096t28_b0824900_m002022_m01011"
    "s064t20_b0803600 s064t22_b0814700 s096t24_b0824900_m002022_m01011"
    )
    arr_output_suffix=("T195" "T220" "T251" "T293")
    arr_ylims=("4 9.9" "3.9 8.9"  "3.5 7.9" "3.5 6.9")
    min_flowradius=0.05
    add_args="--ansatz custom --relflow --nsamples 10000"
fi

for idx in "${!arr_conftypes[@]}" ; do

    if [ "$qcdtype" == hisq_ms5_zeuthenFlow ] ; then
        ylims=${arr_ylims[idx]}
    fi

    if [ "${arr_output_suffix[idx]}" ] ; then
        sufargs="--output_suffix"
    fi
    files=""
    for conftype in ${arr_conftypes[idx]} ; do
        files="$files $conftype/flowradii_$conftype.dat"
    done

#        ../find_common_flowtimes.py --basepath ../../../../data/merged/$qcdtype/$corr/ --files $files --output ../../../../data/merged/$qcdtype/$corr/${arr_output_suffix[idx]}/flowradii_${arr_output_suffix[idx]}.dat
        args="$add_args --use_tex --nproc $ncpu --min_flowradius $min_flowradius --basepath ../../../../data/merged/
        --basepath_plot ../../../../plots/ --max_FlowradiusBytauT_offset 0 --max_FlowradiusBytauT 0.31  --min_FlowradiusBytauT 0.2  $sufargs ${arr_output_suffix[idx]}${addsuf}
        --qcdtype $qcdtype --conftypes ${arr_conftypes[idx]} --corr $corr --custom_ylims $ylims"
        ../_4_continuum_extr.py $args

done
wait