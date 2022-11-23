#!/bin/bash

qcdtype=$1

if [ -z "$qcdtype" ] ; then
    echo "Usage: $0 qcdtype"
    echo "choices for qcdtype: quenched_1.50Tc_zeuthenFlow hisq_ms5_zeuthenFlow"
    exit
fi

if [ "$qcdtype" == quenched_1.50Tc_zeuthenFlow ] ; then
../_5_flowtime_extr.py --corr EE --qcdtype quenched_1.50Tc_zeuthenFlow --custom_ylims 2.2 4 --coarsest_Nt 20 --finest_Nt 36 \
--flowtimes_finest /work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/EE/s144t36_b0754400/flowtimes_s144t36_b0754400.dat \
--plot_all_flowtimes --max_FlowradiusBytauT 0.3 --min_FlowradiusBytauT 0.25 --basepath ../../../data/merged/ --basepath_plot ../../../plots/
elif [ "$qcdtype" == hisq_ms5_zeuthenFlow ] ; then
    coarse_Nt=20
    fine_Nts=(36 32 28 24)
    temps=(196 220 251 296)
    ylims=("4 11" "4 10" "3 9" "3 7.5")
    for idx in "${!temps[@]}" ; do
        ../_5_flowtime_extr.py --corr EE --qcdtype $qcdtype --custom_ylims ${ylims[idx]} --coarsest_Nt $coarse_Nt --finest_Nt ${fine_Nts[idx]} \
        --flowtimes_finest /work/home/altenkort/work/correlators_flow/data/merged/hisq_ms5_zeuthenFlow/EE/s096t${fine_Nts[idx]}_b0824900_m002022_m01011/flowtimes_s096t${fine_Nts[idx]}_b0824900_m002022_m01011.dat \
        --plot_all_flowtimes --temp_subfolder T${temps[idx]} --max_FlowradiusBytauT_offset 0 --max_FlowradiusBytauT 0.33 --min_FlowradiusBytauT 0.25 \
        --basepath ../../../data/merged/ --basepath_plot ../../../plots/
    done

fi

