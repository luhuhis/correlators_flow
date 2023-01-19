#!/bin/bash


extrapolate_coupling(){
    ../extrapolate_coupling.py \
        --ref_scale r0 \
        --calc_cont \
        --input_basepath "/work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/coupling/" \
        --input_files \
        flow_t2E_s096t144_b0754400.dat \
        flow_t2E_s096t120_b0739400.dat \
        flow_t2E_s096t96_b0719200.dat \
        --outputpath_plot "/work/home/altenkort/work/correlators_flow/plots/quenched_1.50Tc_zeuthenFlow/coupling/" \
        --outputpath_data "/work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/coupling/" \
        --Nts 144 120 96 \
        --betas 7.544 7.394 7.192

#flow_t2E_s080t80_b0703500.dat \
# 80
# 7.035
}



compute_Zf2(){
    ../compute_Zf2.py \
    --T_by_Tc 1.5 \
    --g2_file "/work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/coupling/g2_r0_cont_extr.txt" \
    --outputpath_plot "/work/home/altenkort/work/correlators_flow/plots/quenched_1.50Tc_zeuthenFlow/coupling/" \
    --outputpath_data "/work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/coupling/"
}


extrapolate_coupling
compute_Zf2

# flow_t2E_s064t64_b0687361.dat \
# 6.87361
# 64

#    confdict = {"64": "s064t64_b0687361", "80": "s080t80_b0703500", "96": "s096t96_b0719200", "120": "s096t120_b0739400", "144": "s096t144_b0754400"}
