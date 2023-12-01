#!/bin/bash

# TODO convert all paths to inputs of the script instead of hard-coding them

basepath_raw_data=$1
basepath_work_data=$2
basepath_plot=$3

# this path contains "reference" flow times which are useful if there are
# multiple lattice spacings for each temperature that have different flow times
flowradiusbasepath=$5

if [ -z "$basepath_raw_data" ] || [ -z "$basepath_work_data" ] ; then
    echo "Usage: $0 input_basepath output_basepath"
    echo "Example usage: $0 /work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/coupling/ /work/home/altenkort/work/correlators_flow/plots/quenched_1.50Tc_zeuthenFlow/coupling/ /work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/coupling/"
    exit
fi


extrapolate_coupling(){
    ../extrapolate_coupling.py \
        --calc_cont \
        --input_basepath "$basepath_raw_data/quenched_1.50Tc_zeuthenFlow/coupling/" \
        --input_files \
        flow_t2E_s096t144_b0754400.dat \
        flow_t2E_s096t120_b0739400.dat \
        flow_t2E_s096t96_b0719200.dat \
        --outputpath_plot "$basepath_work_data/quenched_1.50Tc_zeuthenFlow/coupling/" \
        --outputpath_data "$basepath_plot/quenched_1.50Tc_zeuthenFlow/coupling/" \
        --Nts 144 120 96 \
        --betas 7.544 7.394 7.192

#flow_t2E_s080t80_b0703500.dat \
# 80
# 7.035
}


(
    cd "$(dirname "$0")" || exit

    mycmd="extrapolate_coupling"
    eval $mycmd

)



