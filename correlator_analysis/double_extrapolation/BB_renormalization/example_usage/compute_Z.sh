#!/bin/bash


# TODO add plots for all the different possibilities of Z_total

compute_Zf2(){
    ../compute_Zf2.py \
    --g2_file \
    "/work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/coupling/g2_MSBAR_runFromMu_4.0.txt" \
    --outputpath_plot "/work/home/altenkort/work/correlators_flow/plots/quenched_1.50Tc_zeuthenFlow/coupling/" \
    --outputpath_data "/work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/coupling/"
}

compute_Zf2
