#!/bin/bash

# ================
# ================ HISQ EE
# ================

qcdtype=hisq_ms5_zeuthenFlow

suffix=("" "_phys")
xlims=("--xlims 0 0.51" "--plot_in_lattice_units --xlims 0 17")



# plot comparing different temperatures of doube-extrapolated correlators
for idz in {0..1} ; do
                ./plot_rec_corr_fixFlowBytauT.py \
                --output_suffix _hisq_final${suffix[idz]} \
                --output_path ~/work/correlators_flow/plots/$qcdtype/EE/rel_flow/ \
                ${xlims[idz]} \
                --no_just_UV --no_kappa_plot --no_connection \
                --ylims -0.5 11 \
                --usetex \
                --title " " \
                --plot_flow_extr \
                /home/altenkort/work/correlators_flow/data/merged/$qcdtype/EE/T196/EE_flow_extr.txt \
                /home/altenkort/work/correlators_flow/data/merged/$qcdtype/EE/T220/EE_flow_extr.txt \
                /home/altenkort/work/correlators_flow/data/merged/$qcdtype/EE/T251/EE_flow_extr.txt \
                /home/altenkort/work/correlators_flow/data/merged/$qcdtype/EE/T296/EE_flow_extr.txt \
                --qcdtype $qcdtype --corr EE \
                --leg_title_suffix "T\\," --leg_loc "center left" --leg_pos 1 0.5 --leg_n_dummies 0 --leg_n_col 1 --leg_label_showNtinsteadofa \
                --leg_label_suffix  \
                "196,\\: " "220,\\: " "251,\\: " "296,\\: "
    #TODO add option for figsize here
done
wait


# plot comparing continuum + flow time to zero with finite lattice spacing and flowtime
for idz in {0..1} ; do
                ./plot_rec_corr_fixFlowBytauT.py \
                --output_suffix _hisq_final_vs_latt${suffix[idz]} \
                --output_path ~/work/correlators_flow/plots/$qcdtype/EE/rel_flow/ \
                ${xlims[idz]} \
                --no_just_UV --no_kappa_plot --no_connection \
                --ylims -0.5 11 \
                --usetex \
                --title " " \
                --plot_flow_extr \
                /home/altenkort/work/correlators_flow/data/merged/$qcdtype/EE/T196/EE_flow_extr.txt \
                /home/altenkort/work/correlators_flow/data/merged/$qcdtype/EE/T220/EE_flow_extr.txt \
                /home/altenkort/work/correlators_flow/data/merged/$qcdtype/EE/T251/EE_flow_extr.txt \
                /home/altenkort/work/correlators_flow/data/merged/$qcdtype/EE/T296/EE_flow_extr.txt \
                --flowradiusBytauT 0.25 \
                --conftype \
                s096t36_b0824900_m002022_m01011 \
                s096t32_b0824900_m002022_m01011 \
                s096t28_b0824900_m002022_m01011 \
                s096t24_b0824900_m002022_m01011 \
                --qcdtype $qcdtype --corr EE \
                --leg_title_suffix "T\\," --leg_loc "center left" --leg_pos 1 0.5 --leg_n_dummies 0 --leg_n_col 1 --leg_label_showNtinsteadofa \
                --leg_label_suffix  \
                "196,\\: " "220,\\: " "251,\\: " "296,\\: " \
                "196,\\: " "220,\\: " "251,\\: " "296,\\: "  \

done
wait


