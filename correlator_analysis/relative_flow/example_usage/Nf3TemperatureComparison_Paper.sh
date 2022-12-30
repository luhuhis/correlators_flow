#!/bin/bash

qcdtype="hisq_ms5_zeuthenFlow"

# plot comparing different temperatures of doube-extrapolated correlators

plot_comparison_of_different_temperatures_of_doublextrapolated_correlators() {
    ../plot_rec_corr_fixFlowBytauT.py \
        --output_suffix _hisq_final \
        --output_path ~/work/correlators_flow/plots/$qcdtype/EE/ \
        --xlims 0.24 0.52 \
        --ylims 3.5 10.5 \
        --no_kappa_plot \
        --usetex \
        --no_connection \
        --conftype \
        s096t20_b0824900_m002022_m01011 \
        --plot_flow_extr \
        /home/altenkort/work/correlators_flow/data/merged/$qcdtype/EE/T196/EE_flow_extr.txt \
        /home/altenkort/work/correlators_flow/data/merged/$qcdtype/EE/T220/EE_flow_extr.txt \
        /home/altenkort/work/correlators_flow/data/merged/$qcdtype/EE/T251/EE_flow_extr.txt \
        /home/altenkort/work/correlators_flow/data/merged/$qcdtype/EE/T296/EE_flow_extr.txt \
        --qcdtype $qcdtype --corr EE \
        --figsize 8.5 7 \
        --leg_title '$T [\mathrm{MeV}]$' --leg_loc "upper center" --leg_pos 0.4 1 --leg_n_dummies 0 --leg_n_col 2  \
        --leg_labels \
        "196" "220" "251" "296" "352" "196" "220" "251" "296" \
        --markers o s D H p o s D H \
        --color_data C0 C1 C2 C3 C4 C0 C1 C2 C3 \
        --fillstyle full full full full none none none none \
        --flowradiusBytauT 0.3 \
        --flow_extr_custom_units 1 1 1 1 #36 32 28 24
}

# --leg_interleave
#        --conftype \
#        s096t36_b0824900_m002022_m01011 \
#        s096t32_b0824900_m002022_m01011 \
#        s096t28_b0824900_m002022_m01011 \
#        s096t24_b0824900_m002022_m01011 \
# --xlims 0 19 \
#        --plot_in_lattice_units \

plot_comparison_of_lattice_spacing_and_temperature_effects() {

    for flowradiusBytauT in "0.25" "0.33"; do

        ../plot_rec_corr_fixFlowBytauT.py \
            --min_flowradius 0.05 \
            --output_suffix _hisq_$flowradiusBytauT \
            --output_path ~/work/correlators_flow/plots/$qcdtype/EE/ \
            --xlims 0.24 0.52 \
            --xticks 0.25 0.3 0.35 0.4 0.45 0.5 \
            --no_kappa_plot \
            --no_just_UV \
            --ylims 3.5 9.5 \
            --flowradiusBytauT $flowradiusBytauT \
            --qcdtype $qcdtype --corr EE \
            --conftype \
            s096t36_b0824900_m002022_m01011 \
            s064t20_b0757000 \
            s096t32_b0824900_m002022_m01011 \
            s064t20_b0770400 \
            s096t28_b0824900_m002022_m01011 \
            s064t20_b0785700 \
            s096t24_b0824900_m002022_m01011 \
            s064t20_b0803600 \
            --color_data \
            C0 C0 C1 C1 C2 C2 C3 C3 \
            --markers o o s s D D H H \
            --fillstyle full none full none full none full none \
            --no_label --no_connection --usetex --figsize 8.5 7 \
            --leg_title '$T [\mathrm{MeV}], N_\tau$' --leg_pos 1 0.5 --leg_loc "center left" --leg_n_dummies 0 \
            --custom_text 0.5 0.99 '$\sqrt{8\tau_\mathrm{F}}/\tau='${flowradiusBytauT}'$' center top rel \
            --leg_labels \
            "196, 36" \
            "196, 20" \
            "220, 32" \
            "220, 20" \
            "251, 28" \
            "251, 20" \
            "296, 24" \
            "296, 20"

    done

}

plot_comparison_of_different_temperatures_of_doublextrapolated_correlators
plot_comparison_of_lattice_spacing_and_temperature_effects
