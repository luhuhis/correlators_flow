#!/bin/bash

qcdtype=$1
corr=$2
basepath_work_data=$3
basepath_plot=$4
if [ -z "$qcdtype" ] || [ -z "$corr" ] || [ -z "$basepath_work_data" ] || [ -z "$basepath_plot" ] ; then
    echo "Usage: $0 qcdtype corr basepath_work_data basepath_plot"
    echo "choices for qcdtype: quenched_1.50Tc_zeuthenFlow hisq_ms5_zeuthenFlow"
    echo "choices for corr: EE BB EE_clover BB_clover"
    echo "Example: $0 hisq_ms5_zeuthenFlow EE ../../../../data/merged/ ../../../../plots/"
    exit
fi

if [ "$qcdtype" == "quenched_1.50Tc_zeuthenFlow" ] ; then

    conftypes_quenched=(
    #"s064t16_b0687361"
    "s080t20_b0703500"
    "s096t24_b0719200"
    "s120t30_b0739400"
    "s144t36_b0754400"
    )

    for conftype in "${conftypes_quenched[@]}" ; do

        (
            cd "$(dirname "$0")" || exit

            ../plot_flow_dependency.py \
         --qcdtype quenched_1.50Tc_zeuthenFlow \
         --corr $corr \
         --conftype $conftype \
         --basepath "$basepath_work_data" \
         --basepath_plot "$basepath_plot" \
         --xlims 0 0.27 \
         --ylims 1 4.4 \
         --ticklocations 0.1 0.25 0.33 0.4 0.45 0.5 \
         --leg_pos 1 0.5 --leg_ncol 1 --leg_lw 0 --leg_pad 0.5 \
         --leg_loc "center left"

         ../plot_flow_dependency.py \
         --qcdtype quenched_1.50Tc_zeuthenFlow \
         --corr $corr \
         --conftype $conftype \
         --basepath "$basepath_work_data" \
         --basepath_plot "$basepath_plot" \
         --xlims 0 0.11 \
         --ylims 1 4.4 \
         --ticklocations 0.1 0.2 0.25 0.3 \
         --leg_pos 1 0.5 --leg_ncol 1 --leg_lw 0 --leg_pad 0.5 \
         --leg_loc "center left" \
         --suffix "zoom"

        )



done
fi


if [ "$qcdtype" == "hisq_ms5_zeuthenFlow" ] ; then


    conftype=s096t24_b0824900_m002022_m01011

        (
            cd "$(dirname "$0")" || exit

            ../plot_flow_dependency.py \
                --qcdtype $qcdtype \
                --corr $corr \
                --conftype $conftype \
                --basepath "$basepath_work_data" \
                --basepath_plot "$basepath_plot" \
                --xlims 0 0.27 \
                --ylims 1 8 \
                --ticklocations 0.1 0.25 0.33 0.4 0.45 0.5 \
                --leg_pos 1 0.5 --leg_ncol 1 --leg_lw 0 --leg_pad 0.5 \
                --leg_loc "center left"

            ../plot_flow_dependency.py \
                --qcdtype $qcdtype \
                --corr $corr \
                --conftype $conftype \
                --basepath "$basepath_work_data" \
                --basepath_plot "$basepath_plot" \
                --xlims 0 0.11 \
                --ylims 1 8 \
                --ticklocations 0.1 0.2 0.25 0.3 \
                --leg_pos 1 0.5 --leg_ncol 1 --leg_lw 0 --leg_pad 0.5 \
                --leg_loc "center left" \
                --suffix "zoom"
        )

fi