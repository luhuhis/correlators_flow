#!/bin/bash

# TODO change phiUV to npy file

# ================
# ================ HISQ EE
# ================
qcdtype=hisq_ms5_zeuthenFlow

wUV=6.28
models_suffix=( maxLO smaxLO plawLO maxNLO smaxNLO plawNLO ) # lineNLO lineLO
models_file=( max pnorm2.0 "plaw_wIR1.0_wUV${wUV}")  # "line_wIR1.0_wUV${wUV}"
models=( max pnorm plaw ) # "line_wIR1.0_wUV${wUV}"
modelparams=( "" "--p 2" "--OmegaByT_IR 1 --OmegaByT_UV ${wUV}" )
order=("LO" "NLO")
scale=("min2piT_w1" "minopt_wopt")
minscale=('2\pi T\,' '\mathrm{opt}')
runscale=('\omega' '\mathrm{opt}')
# --xlims 0.14 0.512  \
#
flows=("0.25" "0.30")  # "0.15" "0.20"
suffix=("" "_phys")
# 0 1 2 ,   0 1
#--kappa_in_GeV  --kappa_xlims 0 0.5 \
#--kappa_ypos 1.45 1.3 1.15 1 0.85
xlims=("--xlims 0 0.505" "--plot_in_lattice_units --xlims 0 17")
#0.8 1.2
#--ylims 0.6 1.4 \
#--norm_by_Gansatz \
for idz in {0..1} ; do
    for idx in "${!models_file[@]}"; do
        for idy in "${!order[@]}" ; do
            for idf in "${!flows[@]}" ; do
                ./plot_rec_corr_fixFlowBytauT.py \
                --output_suffix _hisq_${models_suffix[idx+idy*${#models_file[@]}]}_f${flows[idf]}${suffix[idz]} --npoints 100  \
                --output_path ~/work/correlators_flow/plots/$qcdtype/EE/ \
                ${xlims[idz]} \
                --kappa_xaxis_temp --kappa_xlims 0 15 --kappa_labelheight 0.65 \
                --kappa_temps 196 220 251 296 352 \
                --no_just_UV \
                --ylims -0.25 10 \
                --title " " \
                --model ${models[idx]} ${modelparams[idx]} \
                --flowradiusBytauT "${flows[idf]}" \
                --qcdtype $qcdtype --corr EE \
                --conftype \
                s096t36_b0824900_m002022_m01011 \
                s096t32_b0824900_m002022_m01011 \
                s096t28_b0824900_m002022_m01011 \
                s096t24_b0824900_m002022_m01011 \
                s096t20_b0824900_m002022_m01011 \
                --deduce_fitparam_files --fitparam_basepath ~/work/correlators_flow/data/merged/$qcdtype/EE/ \
                --fitparam_files \
                s096t36_b0824900_m002022_m01011/spf/${models_file[idx]}_${order[idy]}_Nf3_T0.196_${scale[idy]}_500_0.35_exp0_hisq_nt36_f${flows[idf]} \
                s096t32_b0824900_m002022_m01011/spf/${models_file[idx]}_${order[idy]}_Nf3_T0.220_${scale[idy]}_500_0.35_exp0_hisq_nt32_f${flows[idf]} \
                s096t28_b0824900_m002022_m01011/spf/${models_file[idx]}_${order[idy]}_Nf3_T0.251_${scale[idy]}_500_0.35_exp0_hisq_nt28_f${flows[idf]} \
                s096t24_b0824900_m002022_m01011/spf/${models_file[idx]}_${order[idy]}_Nf3_T0.296_${scale[idy]}_500_0.35_exp0_hisq_nt24_f${flows[idf]} \
                s096t20_b0824900_m002022_m01011/spf/${models_file[idx]}_${order[idy]}_Nf3_T0.352_${scale[idy]}_500_0.35_exp0_hisq_nt20_f${flows[idf]} \
                --PhiUV_basepath ~/work/correlators_flow/data/merged/$qcdtype/EE/ \
                --PhiUV_files \
                s096t36_b0824900_m002022_m01011/spf/${models_file[idx]}_${order[idy]}_Nf3_T0.196_${scale[idy]}_500_0.35_exp0_hisq_nt36_f${flows[idf]}/phiUV.dat \
                s096t32_b0824900_m002022_m01011/spf/${models_file[idx]}_${order[idy]}_Nf3_T0.220_${scale[idy]}_500_0.35_exp0_hisq_nt32_f${flows[idf]}/phiUV.dat \
                s096t28_b0824900_m002022_m01011/spf/${models_file[idx]}_${order[idy]}_Nf3_T0.251_${scale[idy]}_500_0.35_exp0_hisq_nt28_f${flows[idf]}/phiUV.dat \
                s096t24_b0824900_m002022_m01011/spf/${models_file[idx]}_${order[idy]}_Nf3_T0.296_${scale[idy]}_500_0.35_exp0_hisq_nt24_f${flows[idf]}/phiUV.dat \
                s096t20_b0824900_m002022_m01011/spf/${models_file[idx]}_${order[idy]}_Nf3_T0.352_${scale[idy]}_500_0.35_exp0_hisq_nt20_f${flows[idf]}/phiUV.dat \
                --min_scale "${minscale[idy]}" --run_scale "${runscale[idy]}" \
                --tauT_vlines 0.35 \
                --min_tauT 0.35 0.35 0.35 0.35 0.35 \
                --no_label --no_connection --show_watermark --custom_text \
                --leg_title_suffix ",\\:\\: T\\," --leg_pos 0.15 1 --leg_n_dummies 0 --leg_label_showNtinsteadofa \
                --leg_label_suffix  ",\\: 196 \\,\\mathrm{MeV}" ",\\: 220\\,\\mathrm{MeV} " ",\\: 251\\,\\mathrm{MeV}" ",\\: 296\\,\\mathrm{MeV}" ",\\: 352\\,\\mathrm{MeV}" \
                &
            done
        done
    done
done
wait