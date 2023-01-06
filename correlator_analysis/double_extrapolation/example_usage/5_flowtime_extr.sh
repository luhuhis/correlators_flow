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

    quenched_extr(){
        ../_5_flowtime_extr.py \
        --corr $corr \
        --qcdtype quenched_1.50Tc_zeuthenFlow \
        --custom_ylims $ylims \
        --coarsest_Nt 20 --finest_Nt 36 \
        --flowtimes_finest /work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/$corr/s144t36_b0754400/flowtimes_s144t36_b0754400.dat \
        --plot_all_flowtimes \
        --max_FlowradiusBytauT 0.3 \
        --min_FlowradiusBytauT 0.25 \
        --basepath ../../../../data/merged/ --basepath_plot ../../../../plots/ \
        --use_tex --min_tauT_plot 0.25 \
        $Zf2file $noextr $output_suffix
    }

    Zf2file=""
    ylims="2.2 4"

    if [ "$corr" == "BB" ] ; then
        ylims="2.5 3.8"
        noextr="--no_extr"
        output_suffix="--output_suffix _no_extr"
        quenched_extr
        Zf2file="--Zf2_file /work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/coupling/Z2_cont.dat"
        ylims="2.9 4"
    fi
    output_suffix=""
    noextr=""
    quenched_extr

elif [ "$qcdtype" == hisq_ms5_zeuthenFlow ] ; then
    coarse_Nt=20
    fine_Nts=(36 32 28 24)
    temps=(196 220 251 296)
    ylims=("4.5 10.2" "4 9" "4 8" "4 7.2")
    for idx in "${!temps[@]}" ; do
        ../_5_flowtime_extr.py --min_tauT_plot 0.25 --corr $corr --qcdtype $qcdtype --custom_ylims ${ylims[idx]} --coarsest_Nt $coarse_Nt --finest_Nt ${fine_Nts[idx]} \
        --flowtimes_finest /work/home/altenkort/work/correlators_flow/data/merged/hisq_ms5_zeuthenFlow/EE/s096t${fine_Nts[idx]}_b0824900_m002022_m01011/flowtimes_s096t${fine_Nts[idx]}_b0824900_m002022_m01011.dat \
        --plot_all_flowtimes --temp_subfolder T${temps[idx]} --max_FlowradiusBytauT_offset 0 --max_FlowradiusBytauT 0.3 --min_FlowradiusBytauT 0.25 \
        --basepath ../../../../data/merged/ --basepath_plot ../../../../plots/ --use_tex --slope_bounds -4 -1
    done

fi
