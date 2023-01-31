#!/bin/bash

input_suffix="23-01-19"
minscale="eff"

plot_kappa() {
        ../plot_kappa.py \
    --xlims -0.5 16.5 \
    --xticks 0 5 10 15 \
    --model_ids "${models[@]}" \
    --labels "${labels[@]}" \
    --basepath "${basepath}" \
    --outputpath "${outputpath}" \
    --suffix "${suffix}" \
    --corr EE \
    --figsize 7 7 \
    --colors C0 C1 C2 C3 C4 C5 C6 C7 C8 C9 \
    --outputpath_data $outputpath_data
}

plot_fitcorr() {
        ../plot_fitcorr.py \
    --model_ids "${models[@]}" \
    --labels "${labels[@]}" \
    --basepath "${basepath}" \
    --outputpath "${outputpath}" \
    --suffix "${suffix}" \
    --corr EE \
    --yticks 0.9 0.95 1 1.05 1.1 \
    --ylims 0.75 1.25 \
    --colors C0 C1 C2 C3 C4 C5 C6 C7 C8 C9
}

plot_spfs(){
        ../plot_spfs.py \
    --model_ids "${models[@]}" \
    --labels "${labels[@]}" \
    --basepath "${basepath}" \
    --outputpath "${outputpath}" \
    --suffix "${suffix}" \
    --corr EE \
    --plot_spf_err \
    --ylims 1 1500 \
    --xlims 0.1 150 \
    --colors C0 C1 C2 C3 C4 C5 C6 C7 C8 C9
}


set_labels_and_models(){

if [ "${selector}" == "paper" ] ; then
labels_and_models=(
    '$\text{max}_\text{NLO}$'                     "max_NLO_Nf3_T0.${temp}_min${minscale}_wopt_${nsamples}smpls_tauTgtr0.24_${input_suffix}"
    '$\text{max}_\text{LO\hphantom{N}}$'          "max_LO_Nf3_T0.${temp}_min${minscale}_w1_${nsamples}smpls_tauTgtr0.24_${input_suffix}"
    '$\text{smax}_\text{NLO}$'                    "smax_NLO_Nf3_T0.${temp}_min${minscale}_wopt_${nsamples}smpls_tauTgtr0.24_${input_suffix}"
    '$\text{smax}_\text{LO\hphantom{N}}$'         "smax_LO_Nf3_T0.${temp}_min${minscale}_w1_${nsamples}smpls_tauTgtr0.24_${input_suffix}"
    '$\text{plaw}_\text{NLO}$'                    "plaw_wIR1.0_wUV6.2832_NLO_Nf3_T0.${temp}_min${minscale}_wopt_${nsamples}smpls_tauTgtr0.24_${input_suffix}"
    '$\text{plaw}_\text{LO\hphantom{N}}$'         "plaw_wIR1.0_wUV6.2832_LO_Nf3_T0.${temp}_min${minscale}_w1_${nsamples}smpls_tauTgtr0.24_${input_suffix}"
#    '$\text{trig1}_\text{NLO}$'                 "trig_NLO__beta_1_Nf3_T0.${temp}_min${minscale}_wopt_${nsamples}smpls_tauTgtr0.24_${input_suffix}"
    #'$\text{trig1}_\text{LO\hphantom{N}}$'        "fourier_LO_s1_alpha_1_Nf3_T0.${temp}_min${minscale}_w1_${nsamples}smpls_tauTgtr0.24_${input_suffix}"
    #'$\text{trig2}_\text{NLO}$'                   "fourier_NLO_s1_beta_2_Nf3_T0.${temp}_min${minscale}_wopt_${nsamples}smpls_tauTgtr0.24_${input_suffix}"
    #'$\text{trig2}_\text{LO\hphantom{N}}$'        "fourier_LO_s1_alpha_2_Nf3_T0.${temp}_min${minscale}_w1_${nsamples}smpls_tauTgtr0.24_${input_suffix}"
    )
    add_suffix="_paper"
else
    mintauT=0.24
    labels_and_models=(
    '$\text{max}_\text{NLO}$'                     "max_NLO_Nf3_T0.${temp}_min${minscale}_wopt_${nsamples}smpls_tauTgtr${mintauT}_${input_suffix}"
    '$\text{max}_\text{LO\hphantom{N}}$'          "max_LO_Nf3_T0.${temp}_min${minscale}_w1_${nsamples}smpls_tauTgtr${mintauT}_${input_suffix}"
    '$\text{smax}_\text{NLO}$'                    "smax_NLO_Nf3_T0.${temp}_min${minscale}_wopt_${nsamples}smpls_tauTgtr${mintauT}_${input_suffix}"
    '$\text{smax}_\text{LO\hphantom{N}}$'         "smax_LO_Nf3_T0.${temp}_min${minscale}_w1_${nsamples}smpls_tauTgtr${mintauT}_${input_suffix}"
    '$\text{plaw}_\text{NLO}$'                    "plaw_any_wUV6.2832_NLO_Nf3_T0.${temp}_min${minscale}_wopt_${nsamples}smpls_tauTgtr${mintauT}_${input_suffix}"
    '$\text{plaw}_\text{LO\hphantom{N}}$'         "plaw_any_wUV6.2832_LO_Nf3_T0.${temp}_min${minscale}_w1_${nsamples}smpls_tauTgtr${mintauT}_${input_suffix}"
#    '$\text{plaw}_\text{NLO}$'                    "plaw_wIR1.0_wUV6.2832_NLO_Nf3_T0.${temp}_min${minscale}_wopt_${nsamples}smpls_tauTgtr${mintauT}_${input_suffix}"
#    '$\text{plaw}_\text{LO\hphantom{N}}$'         "plaw_wIR1.0_wUV6.2832_LO_Nf3_T0.${temp}_min${minscale}_w1_${nsamples}smpls_tauTgtr${mintauT}_${input_suffix}"
    '$\text{trig1}_\text{NLO}$'                   "trig_NLO__beta_1_Nf3_T0.${temp}_min${minscale}_wopt_${nsamples}smpls_tauTgtr${mintauT}_${input_suffix}"
    '$\text{trig1}_\text{LO\hphantom{N}}$'        "trig_LO__alpha_1_Nf3_T0.${temp}_min${minscale}_w1_${nsamples}smpls_tauTgtr${mintauT}_${input_suffix}"
#    '$\text{trig2}_\text{NLO}$'                 "trig_NLO__beta_2_Nf3_T0.${temp}_min${minscale}_wopt_${nsamples}smpls_tauTgtr${mintauT}_${input_suffix}"
#    '$\text{trig2}_\text{LO\hphantom{N}}$'      "trig_LO__alpha_2_Nf3_T0.${temp}_min${minscale}_w1_${nsamples}smpls_tauTgtr${mintauT}_${input_suffix}"
    )
    add_suffix=""
fi
}

conftypes=(s096t36_b0824900_m002022_m01011 s096t32_b0824900_m002022_m01011 s096t28_b0824900_m002022_m01011 s096t24_b0824900_m002022_m01011 s096t20_b0824900_m002022_m01011)
temps=(196 220 251 296 352)
for selector2 in "zeroflow" "finiteflow" ; do
    for selector in "thesis" "paper" "finiteflow" ; do
        for j in "${!temps[@]}" ; do

            temp=${temps[j]}
            nsamples=500


            set_labels_and_models

            if [ "${selector2}" == "finiteflow" ] ; then
                conftype=${conftypes[j]}
                suffix="_T${temp}_finiteflow${add_suffix}"
                basepath=/work/home/altenkort/work/correlators_flow/data/merged/hisq_ms5_zeuthenFlow/EE//${conftype}/spf/
                outputpath=/work/home/altenkort/work/correlators_flow/plots/hisq_ms5_zeuthenFlow/EE//${conftype}/
                outputpath_data=/work/home/altenkort/work/correlators_flow/data/merged/hisq_ms5_zeuthenFlow/EE/${conftype}/
            else
                if [ "${temp}" == "352" ]; then
                    continue
                fi
                suffix="_T${temp}${add_suffix}"
                basepath=/work/home/altenkort/work/correlators_flow/data/merged/hisq_ms5_zeuthenFlow/EE//T${temp}/spf/
                outputpath=/work/home/altenkort/work/correlators_flow/plots/hisq_ms5_zeuthenFlow/EE//T${temp}/
                outputpath_data=/work/home/altenkort/work/correlators_flow/data/merged/hisq_ms5_zeuthenFlow/EE/T${temp}
            fi
            models=()
            labels=()
            for i in "${!labels_and_models[@]}" ; do
                if [ $((i % 2)) -eq 0 ] ; then
                    labels+=("${labels_and_models[i]}")
                else
                    models+=("${labels_and_models[i]}")
                fi
            done

            plot_spfs &
            plot_kappa &
            plot_fitcorr &

        done
    done
done
wait