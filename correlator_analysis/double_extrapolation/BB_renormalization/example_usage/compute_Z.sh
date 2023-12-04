#!/bin/bash

basepath_work_data=$1
basepath_plot=$2


compute_Z(){
    ../compute_Z.py \
    --g2_file "$basepath_work_data/quenched_1.50Tc_zeuthenFlow/coupling/g2_MSBAR_runFromMu_4.0.txt" \
    --outputpath_plot "$basepath_plot/quenched_1.50Tc_zeuthenFlow/coupling/" \
    --outputpath_data "$basepath_work_data/quenched_1.50Tc_zeuthenFlow/coupling/"
}


(
    cd "$(dirname "$0")" || exit
    mycmd="compute_Z"
    eval $mycmd
)



