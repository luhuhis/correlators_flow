#!/bin/bash

qcdtype=$1

# TODO add BB

if [ -z "$qcdtype" ] ; then
    echo "Usage: $0 qcdtype"
    echo "choices for qcdtype: quenched_1.50Tc_zeuthenFlow hisq_ms5_zeuthenFlow"
    exit
fi

ncpu=20

if [ "$qcdtype" == quenched_1.50Tc_zeuthenFlow ] ; then
    arr_conftypes=("s080t20_b0703500 s096t24_b0719200 s120t30_b0739400 s144t36_b0754400")
    continuum="../../../../data/merged/quenched_1.50Tc_zeuthenFlow/EE/cont_extr/EE_cont.dat"
    continuum_err="../../../../data/merged/quenched_1.50Tc_zeuthenFlow/EE/cont_extr/EE_cont_err.dat"
    flowtimesT2="../../../../data/merged/quenched_1.50Tc_zeuthenFlow/EE/cont_extr/EE_cont_flowtimesT2.dat"
    ylims="2 4.1"
    add_args="--hide_cont --xlims 0 0.52 --lower_limit_text_pos 3.5"
elif [ "$qcdtype" == hisq_ms5_zeuthenFlow ] ; then

    # TODO make plots for hisq
    arr_conftypes=(\
    "s064t20_b0803600 s064t22_b0814700 s096t24_b0824900_m002022_m01011"
    "s064t20_b0785700 s064t24_b0806800 s096t28_b0824900_m002022_m01011"
    "s064t20_b0770400 s064t24_b0791300 s096t32_b0824900_m002022_m01011"
    "s064t20_b0757000 s064t24_b0777700 s096t36_b0824900_m002022_m01011")
    arr_output_suffix=("T296" "T251" "T220" "T196")
    ylims_arr=("3.5 10.5" "3.5 10.5" "3.5 10.5" "3.5 10.5")
    basepath="../../../../data/merged/hisq_ms5_zeuthenFlow/EE"
    add_args="--hide_cont --lower_limit_text_pos 8"
fi


for idx in "${!arr_conftypes[@]}" ; do

    if [ "$qcdtype" == hisq_ms5_zeuthenFlow ] ; then
        temp=${arr_output_suffix[idx]}
        flowtimesT2="${basepath}/${temp}/cont_extr/EE_cont_flowtimesT2.dat"
        continuum="${basepath}/${temp}/cont_extr/EE_cont.dat"
        continuum_err="${basepath}/${temp}/cont_extr/EE_cont_err.dat"
        outputfolder="--outputfolder ../../../../plots/$qcdtype/EE/$temp/"
        ylims=${ylims_arr[idx]}
    fi

../2_plot_lateffects.py \
--corr EE --qcdtype $qcdtype \
--conftypes ${arr_conftypes[idx]} \
--continuum ${continuum} \
--continuum_err ${continuum_err} \
--flowtimesT2 ${flowtimesT2} \
--nproc $ncpu \
--ylims $ylims \
$outputfolder \
$add_args

done


