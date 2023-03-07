#!/bin/bash

qcdtype="quenched_1.50Tc_zeuthenFlow"

for flowradiusBytauT in "0.2" "0.25" "0.30" ; do
../plot_rec_corr_fixFlowBytauT.py \
    --output_suffix _quenched_$flowradiusBytauT --xlims 0.15 0.52 --xticks 0.2 0.3 0.4 0.5 \
    --min_flowradius 0.05 \
    --output_path ~/work/correlators_flow/plots/$qcdtype/EE/ \
    --ylims 2.1 3.7 \
    --flowradiusBytauT $flowradiusBytauT \
    --qcdtype $qcdtype --corr EE \
    --conftype \
    s144t36_b0754400 \
    s120t30_b0739400 \
    s096t24_b0719200 \
    s080t20_b0703500 \
    --color_data \
    C0 C1 C2 C3 \
    --markers o s D H \
    --fillstyle none none none none \
    --no_label --no_connection --usetex --figsize 7 7 \
    --leg_title '$N_\tau$' --leg_pos 1 0.05 --leg_loc "lower right" --leg_n_col 1 --leg_n_dummies 0 \
    --custom_text 0.5 0.99 '$\sqrt{8\tau_\mathrm{F}}/\tau='${flowradiusBytauT}'$' center top rel \
    --leg_labels \
    "36" "32" "24" "20"

done