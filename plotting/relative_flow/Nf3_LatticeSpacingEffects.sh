#!/bin/bash

flows=("0.20" "0.25" "0.30")
conftypes=("s096t36_b0824900_m002022_m01011 s064t24_b0777700" "s096t32_b0824900_m002022_m01011 s064t24_b0791300"
"s096t28_b0824900_m002022_m01011 s064t20_b0785700" "s096t24_b0824900_m002022_m01011 s064t20_b0803600")
labels_temp=("196" "220" "251" "296")
a=("0.028 0.0421" "0.028 0.0374" "0.028 0.0392" "0.028 0.0336")
for idx in "${!conftypes[@]}" ; do
    for flow in "${flows[@]}" ; do

        ./plot_rec_corr_fixFlowBytauT.py \
        --output_suffix _hisq_lat_effects_${labels_temp[idx]}_$flow  \
        --output_path ~/work/correlators_flow/plots/hisq_ms5_zeuthenFlow/EE/ \
        --ylims 0 10 \
        --xlims 0 0.51 \
        --title " " \
        --flowradiusBytauT "$flow" \
        --qcdtype hisq_ms5_zeuthenFlow --corr EE \
        --conftype \
        ${conftypes[idx]} \
        --no_kappa_plot \
        --no_label \
        --plot_extr /home/altenkort/work/correlators_flow/data/merged/hisq_ms5_zeuthenFlow/EE/T${labels_temp[idx]} \
        --leg_title_suffix ",\\:\\: T\\," --leg_pos 0.9 0.1 --leg_n_dummies 0 --leg_label_showNtinsteadofa \
        --leg_label_suffix  \
        ",\\: ${labels_temp[idx]}\\,\\mathrm{MeV}" \
        ",\\: ${labels_temp[idx]}\\,\\mathrm{MeV}" \


    done
done
#--no_connection
#",\\: 251\\,\\mathrm{MeV}" \


#* nt * lattice_spacing