#!/bin/bash


plot_kappa() {
        ../plot_kappa.py \
    --xlims -0.5 15 \
    --xticks 0 5 10 15 \
    --model_ids "${models[@]}" \
    --labels "${labels[@]}" \
    --basepath "${basepath}" \
    --outputpath "${outputpath}" \
    --suffix "${suffix}" \
    --corr EE \
    --scale_error_by_chisqdof \
    --figsize 7 7 \
    --colors C0 C1 C2 C3 C4 C5 C6 C7 C8 C9
}

plot_fitcorr() {
        ../plot_fitcorr.py \
    --model_ids "${models[@]}" \
    --labels "${labels[@]}" \
    --basepath "${basepath}" \
    --outputpath "${outputpath}" \
    --suffix "${suffix}" \
    --corr EE \
    --yticks 0.98 0.99 1 1.01 1.02 \
    --ylims 0.94 1.04 \
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
    --ylims 1 1000 \
    --colors C0 C1 C2 C3 C4 C5 C6 C7 C8 C9
}


set_labels_and_models(){

if [ "${selector}" == "paper" ] ; then
labels_and_models=(
    '$\text{max}_\text{NLO}$'                     "max_NLO_Nf3_T0.${temp}_min2piT_wopt_${lowsamples}smpls_tauTgtr0.24_paper"
    #'$\text{max}_\text{LO\hphantom{N}}$'          "max_LO_Nf3_T0.${temp}_min2piT_w1_${lowsamples}smpls_tauTgtr0.24_paper"
    '$\text{smax}_\text{NLO}$'                    "smax_NLO_Nf3_T0.${temp}_min2piT_wopt_${lowsamples}smpls_tauTgtr0.24_paper"
    #'$\text{smax}_\text{LO\hphantom{N}}$'         "smax_LO_Nf3_T0.${temp}_min2piT_w1_${lowsamples}smpls_tauTgtr0.24_paper"
    '$\text{plaw}_\text{NLO}$'                    "plaw_wIR1.0_wUV6.28_NLO_Nf3_T0.${temp}_min2piT_wopt_${lowsamples}smpls_tauTgtr0.24_paper"
    #'$\text{plaw}_\text{LO\hphantom{N}}$'         "plaw_wIR1.0_wUV6.28_LO_Nf3_T0.${temp}_min2piT_w1_${lowsamples}smpls_tauTgtr0.24_paper"
    '$\text{trig1}_\text{NLO}$'                 "trig_NLO__beta_1_Nf3_T0.${temp}_min2piT_wopt_1000smpls_tauTgtr0.24_paper"
    #'$\text{trig1}_\text{LO\hphantom{N}}$'        "fourier_LO_s1_alpha_1_Nf3_T0.${temp}_min2piT_w1_${lowsamples}smpls_tauTgtr0.24_paper"
    #'$\text{trig2}_\text{NLO}$'                   "fourier_NLO_s1_beta_2_Nf3_T0.${temp}_min2piT_wopt_${lowsamples}smpls_tauTgtr0.24_paper"
    #'$\text{trig2}_\text{LO\hphantom{N}}$'        "fourier_LO_s1_alpha_2_Nf3_T0.${temp}_min2piT_w1_${lowsamples}smpls_tauTgtr0.24_paper"
    )
    add_suffix="_paper"
else
    labels_and_models=(
    '$\text{max}_\text{NLO}$'                     "max_NLO_Nf3_T0.${temp}_min2piT_wopt_${lowsamples}smpls_tauTgtr0.24_paper"
    '$\text{max}_\text{LO\hphantom{N}}$'          "max_LO_Nf3_T0.${temp}_min2piT_w1_${lowsamples}smpls_tauTgtr0.24_paper"
    '$\text{smax}_\text{NLO}$'                    "smax_NLO_Nf3_T0.${temp}_min2piT_wopt_${lowsamples}smpls_tauTgtr0.24_paper"
    '$\text{smax}_\text{LO\hphantom{N}}$'         "smax_LO_Nf3_T0.${temp}_min2piT_w1_${lowsamples}smpls_tauTgtr0.24_paper"
    '$\text{plaw}_\text{NLO}$'                    "plaw_wIR1.0_wUV6.28_NLO_Nf3_T0.${temp}_min2piT_wopt_${lowsamples}smpls_tauTgtr0.24_paper"
    '$\text{plaw}_\text{LO\hphantom{N}}$'         "plaw_wIR1.0_wUV6.28_LO_Nf3_T0.${temp}_min2piT_w1_${lowsamples}smpls_tauTgtr0.24_paper"
#    '$\text{trig1}_\text{NLO}$'                 "trig_NLO__beta_1_Nf3_T0.${temp}_min2piT_wopt_1000smpls_tauTgtr0.24_paper"
#    '$\text{trig1}_\text{LO\hphantom{N}}$'      "trig_LO__alpha_1_Nf3_T0.${temp}_min2piT_w1_1000smpls_tauTgtr0.24_paper"
#    '$\text{trig2}_\text{NLO}$'                 "trig_NLO__beta_2_Nf3_T0.${temp}_min2piT_wopt_1000smpls_tauTgtr0.24_paper"
#    '$\text{trig2}_\text{LO\hphantom{N}}$'      "trig_LO__alpha_2_Nf3_T0.${temp}_min2piT_w1_1000smpls_tauTgtr0.24_paper"
    )
    add_suffix=""
fi
}


temps=( 196 220 251 296)
for selector in "thesis" "paper" ; do
    for j in "${!temps[@]}" ; do

        temp=${temps[j]}
        lowsamples=1000

        set_labels_and_models

        suffix="_T${temp}${add_suffix}"
        basepath=/work/home/altenkort/work/correlators_flow/data/merged/hisq_ms5_zeuthenFlow/EE//T${temp}/spf/
        outputpath=/work/home/altenkort/work/correlators_flow/plots/hisq_ms5_zeuthenFlow/EE//T${temp}/


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
wait