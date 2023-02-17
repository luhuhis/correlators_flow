#!/bin/bash

qcdtype=$1
corr=$2
if [ -z "$qcdtype" ] || [ -z "$corr" ] ; then
    echo "Usage: $0 qcdtype corr"
    echo "choices for qcdtype: quenched_1.50Tc_zeuthenFlow hisq_ms5_zeuthenFlow"
    echo "choices for corr: EE BB EE_clover BB_clover"
    exit
fi

ncpu=20

if [ "$qcdtype" == quenched_1.50Tc_zeuthenFlow ] ; then
    arr_conftypes=("s080t20_b0703500 s096t24_b0719200 s120t30_b0739400 s144t36_b0754400")
    continuum="../../../../data/merged/quenched_1.50Tc_zeuthenFlow/${corr}/cont_extr/${corr}_cont.dat"
    continuum_err="../../../../data/merged/quenched_1.50Tc_zeuthenFlow/${corr}/cont_extr/${corr}_cont_err.dat"
    flowtimesT2="../../../../data/merged/quenched_1.50Tc_zeuthenFlow/${corr}/cont_extr/${corr}_cont_flowtimesT2.dat"
    ylims="1 4.1"
    outputfolder="--outputfolder /work/home/altenkort/work/correlators_flow/plots/$qcdtype/${corr}"
    add_args="--hide_cont --xlims 0 0.52 --lower_limit_text_pos 3.5"
elif [ "$qcdtype" == hisq_ms5_zeuthenFlow ] ; then

    arr_conftypes=(\
    "s064t20_b0803600 s064t22_b0814700 s096t24_b0824900_m002022_m01011"
    "s064t20_b0785700 s064t24_b0806800 s096t28_b0824900_m002022_m01011"
    "s064t20_b0770400 s064t24_b0791300 s096t32_b0824900_m002022_m01011"
    "s064t20_b0757000 s064t24_b0777700 s096t36_b0824900_m002022_m01011")
    arr_output_suffix=("T293" "T251" "T220" "T195")
    ylims_arr=("3 10.5" "3 10.5" "3 10.5" "3 10.5")
    basepath="../../../../data/merged/hisq_ms5_zeuthenFlow/${corr}"
    add_args="--hide_cont --lower_limit_text_pos 8"
fi


for idx in "${!arr_conftypes[@]}" ; do

    if [ "$qcdtype" == hisq_ms5_zeuthenFlow ] ; then
        temp=${arr_output_suffix[idx]}
        flowtimesT2="${basepath}/${temp}/cont_extr/${corr}_cont_flowtimesT2.dat"
        continuum="${basepath}/${temp}/cont_extr/${corr}_cont.dat"
        continuum_err="${basepath}/${temp}/cont_extr/${corr}_cont_err.dat"
        outputfolder="--outputfolder /work/home/altenkort/work/correlators_flow/plots/$qcdtype/${corr}/$temp/"
        ylims=${ylims_arr[idx]}
    fi

../_2_plot_lateffects.py \
--corr ${corr} --qcdtype $qcdtype \
--conftypes ${arr_conftypes[idx]} \
--continuum ${continuum} \
--continuum_err ${continuum_err} \
--flowtimesT2 ${flowtimesT2} \
--nproc $ncpu \
--ylims $ylims \
$outputfolder $add_args

done


wait