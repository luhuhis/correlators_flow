#!/bin/bash

qcdtype="hisq_ms5_zeuthenFlow"

basepath=${1:-"/work/home/altenkort/work/correlators_flow/data/merged/"}
basepath_work_data="$basepath/$qcdtype/EE/"
basepath_plot="${2:-"/work/home/altenkort/work/correlators_flow/plots/"}/${qcdtype}/EE/"

suffix="" #_paper"


plot_comparison_of_different_temperatures_of_doublextrapolated_correlators_thesis() {

    params=(
#        "--xlims 0 19 --xticks 0 3 6 9 12 15 18 --hide_fits --x_scales 36 32 28 24 20 36 32 28 24 --xlabel $\tau/a$ --output_suffix _hisq_final_phys --flow_extr_custom_units 36 32 28 24"
        "--xlims 0.24 0.52 --xticks 0.25 0.3 0.35 0.4 0.45 0.5 --output_suffix _hisq_final --flow_extr_custom_units 1 1 1 1"
    )

    fit_args=(
        "--npoints 50 --show_UV_corrs"
        "--fit_basepath $basepath_work_data"
        "--fit_folders"
        "T195/spf/smax_NLO_Nf3_T0.195_mineff_wopt_1000smpls_tauTgtr0.24_23-02-16_relflow"
        "T220/spf/smax_NLO_Nf3_T0.220_mineff_wopt_1000smpls_tauTgtr0.24_23-02-16_relflow"
        "T251/spf/smax_NLO_Nf3_T0.251_mineff_wopt_1000smpls_tauTgtr0.24_23-02-16_relflow"
        "T293/spf/smax_NLO_Nf3_T0.293_mineff_wopt_1000smpls_tauTgtr0.24_23-02-16_relflow"
    )

    for x_units in "${params[@]}" ; do

        # TODO revert ylims 4 9.5, leg_pos to 0 0.6 and leg_loc to "center left"
        ../plot_rec_corr_fixFlowBytauT.py \
            --output_path ${basepath_plot} --basepath ${basepath}\
            ${x_units} \
            --ylims 3 9.5 \
            --usetex \
            ${fit_args[@]} \
            --no_connection \
            --plot_flow_extr \
            $basepath_work_data/T195${suffix}/EE_flow_extr_relflow.txt \
            $basepath_work_data/T220${suffix}/EE_flow_extr_relflow.txt \
            $basepath_work_data/T251${suffix}/EE_flow_extr_relflow.txt \
            $basepath_work_data/T293${suffix}/EE_flow_extr_relflow.txt \
            --qcdtype $qcdtype --corr EE \
            --figsize 7 7 \
            --leg_title '$T [\mathrm{MeV}]$' --leg_pos 1 0 --leg_loc "lower right" --leg_n_dummies 3 --leg_n_col 2 --leg_framealpha 0 \
            --leg_labels \
            "195" "220" "251" "293" "195" "220" "251" "293" \
            --markers o s D H o s D H \
            --color_data C0 C1 C2 C3 C0 C1 C2 C3 \
            --fillstyle none none none none full full full full  \
            --flowradiusBytauT 0.3 --min_flowradius 0.09 \

    done

}

plot_comparison_of_different_temperatures_of_doublextrapolated_correlators_paper() {

    params=(
#        "--xlims 0 19 --xticks 0 3 6 9 12 15 18 --hide_fits --x_scales 36 32 28 24 20 36 32 28 24 --xlabel $\tau/a$ --output_suffix _hisq_final_phys --flow_extr_custom_units 36 32 28 24"
        "--xlims 0.24 0.52 --xticks 0.25 0.3 0.35 0.4 0.45 0.5 --output_suffix _hisq_final --flow_extr_custom_units 1 1 1 1"
    )

    fit_args=(
        "--npoints 50 --show_UV_corrs"
        "--fit_basepath $basepath_work_data"
        "--fit_folders"
        "T195/spf/smax_NLO_Nf3_T0.195_mineff_wopt_1000smpls_tauTgtr0.24_23-02-16_relflow"
        "T220/spf/smax_NLO_Nf3_T0.220_mineff_wopt_1000smpls_tauTgtr0.24_23-02-16_relflow"
        "T251/spf/smax_NLO_Nf3_T0.251_mineff_wopt_1000smpls_tauTgtr0.24_23-02-16_relflow"
        "T293/spf/smax_NLO_Nf3_T0.293_mineff_wopt_1000smpls_tauTgtr0.24_23-02-16_relflow"
        "s096t20_b0824900_m002022_m01011/spf/smax_NLO_Nf3_T0.352_mineff_wopt_1000smpls_tauTgtr0.26_23-02-16_0.30" #TODO remove this and conftype below later for thesis
    )

    for x_units in "${params[@]}" ; do

        # TODO revert ylims 4 9.5, leg_pos to 0 0.6 and leg_loc to "center left"
        ../plot_rec_corr_fixFlowBytauT.py \
            --output_path $basepath_plot --basepath $basepath\
            ${x_units} \
            --ylims 3 9.5 \
            --usetex \
            ${fit_args[@]} \
            --no_connection \
            --conftype s096t20_b0824900_m002022_m01011 \
            --plot_flow_extr \
            $basepath_work_data/T195${suffix}/EE_flow_extr_relflow.txt \
            $basepath_work_data/T220${suffix}/EE_flow_extr_relflow.txt \
            $basepath_work_data/T251${suffix}/EE_flow_extr_relflow.txt \
            $basepath_work_data/T293${suffix}/EE_flow_extr_relflow.txt \
            --qcdtype $qcdtype --corr EE \
            --figsize 7 7 \
            --leg_title '$T [\mathrm{MeV}]$' --leg_pos 1 0 --leg_loc "lower right" --leg_n_dummies 3 --leg_n_col 2 --leg_framealpha 0 \
            --leg_labels \
            "195" "220" "251" "293" "352*" "195" "220" "251" "293" \
            --markers o s D H p o s D H \
            --color_data C0 C1 C2 C3 C4 C0 C1 C2 C3 \
            --fillstyle none none none none full full full full  \
            --flowradiusBytauT 0.3 --min_flowradius 0.09 \

    done

}

plot_comparison_of_lattice_spacing_and_temperature_effects() {

    for flowradiusBytauT in "0.25" "0.30"; do
        xparams=(
            "--output_suffix _hisq_${flowradiusBytauT} --xlims 0.24 0.52 --xticks 0.25 0.3 0.35 0.4 0.45 0.5"
        )
        for xparam in "${xparams[@]}" ; do

        ../plot_rec_corr_fixFlowBytauT.py \
            ${xparam} \
            --min_flowradius 0.05 \
            --output_path $basepath_plot --basepath $basepath\
            --ylims 3 9.25 \
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
            --fillstyle none full none full none full none full \
            --no_label --no_connection --usetex --figsize 7 7 \
            --leg_title '$T [\mathrm{MeV}], N_\tau$' --leg_pos 1 0.05 --leg_loc "lower right" --leg_n_col 2 --leg_n_dummies 0 \
            --custom_text 0.5 0.99 '$\sqrt{8\tau_\mathrm{F}}/\tau='${flowradiusBytauT}'$' center top rel \
            --leg_labels \
            "195, 36" \
            "195, 20" \
            "220, 32" \
            "220, 20" \
            "251, 28" \
            "251, 20" \
            "293, 24" \
            "293, 20" \

        done
    done

}

plot_comparison_of_temperature_effects_phys() {

    for flowradiusBytauT in "0.25" "0.30"; do
        xparams=(
            "--xlims 0 19 --xticks 0 3 6 9 12 15 18 --x_scales 36 32 28 24 20  --xlabel $\tau/a$ --output_suffix _hisq_${flowradiusBytauT}_phys"
        )
        for xparam in "${xparams[@]}" ; do

        ../plot_rec_corr_fixFlowBytauT.py \
            ${xparam} \
            --min_flowradius 0.05 \
            --output_path $basepath_plot --basepath $basepath\
            --ylims 0 9.25 \
            --flowradiusBytauT $flowradiusBytauT \
            --qcdtype $qcdtype --corr EE \
            --conftype \
            s096t36_b0824900_m002022_m01011 \
            s096t32_b0824900_m002022_m01011 \
            s096t28_b0824900_m002022_m01011 \
            s096t24_b0824900_m002022_m01011 \
            s096t20_b0824900_m002022_m01011 \
            --color_data \
            C0 C1 C2 C3 C4 \
            --markers o s D H p \
            --fillstyle none none none none none \
            --no_label --no_connection --usetex --figsize 7 7 \
            --leg_title '$T [\mathrm{MeV}], N_\tau$' --leg_pos 1 0.05 --leg_loc "lower right" --leg_n_col 2 --leg_n_dummies 0 \
            --custom_text 0.5 0.99 '$\sqrt{8\tau_\mathrm{F}}/\tau='${flowradiusBytauT}'$' center top rel \
            --leg_labels \
            "195, 36" \
            "220, 32" \
            "251, 28" \
            "293, 24" \
            "352, 20" \

        done
    done

}

(
    cd "$(dirname $0)" || exit

    plot_comparison_of_different_temperatures_of_doublextrapolated_correlators_thesis
    # plot_comparison_of_different_temperatures_of_doublextrapolated_correlators_paper
    plot_comparison_of_lattice_spacing_and_temperature_effects
    # plot_comparison_of_temperature_effects_phys


)
