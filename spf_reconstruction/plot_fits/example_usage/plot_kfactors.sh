#!/usr/bin/env bash

basepath_work_data="${1:-"/work/home/altenkort/work/correlators_flow/data/merged/"}/hisq_ms5_zeuthenFlow/EE/"
basepath_plot="${2:-"/work/home/altenkort/work/correlators_flow/plots"}/hisq_ms5_zeuthenFlow/EE/"

# TODO check whether this is even correct
input_suffix="23-02-16_relflow" # "23-01-19"
minscale="eff"
nsamples=1000 #500

# TODO fix nt=20

suffix=""

plot_kfactors() {

    ../plot_kfactor.py \
    --temperature_on_xaxis \
    --xlims 1 2.1 --ylims 0 2.9 \
    --xticks 1 1.5 2.0 \
    --model_ids "${models[@]}" \
    --labels '\begin{flushright}$\scriptstyle a\rightarrow 0,$ \\[-1ex] $\scriptstyle\tau_\mathrm{F}\rightarrow 0$\hphantom{,} \end{flushright}' "NLO" "LO" "" "" "" "" "" "" \
             '\begin{flushleft}{\scriptsize nonzero} \\[-1.5ex] $\scriptstyle a$ {\scriptsize and} $\scriptstyle\tau_\mathrm{F}$\end{flushleft}' "NLO" "LO" "" "" "" "" "" "" "" "" \
    --basepath "${basepath_work_data}" \
    --outputpath "${basepath_plot}" \
    --suffix "${suffix}" \
    --corr EE \
    --figsize 7 7 \
    --colors k C0 C1 C0 C1 C0 C1 C0 C1 \
             k C2 C3 C2 C3 C2 C3 C2 C3 C2 C3 \
    --fillstyles none none none none none none none none none \
                 none full full full full full full full full full full \
    --fmts '.' 's' 'D' 's' 'D' 's' 'D' 's' 'D'\
           '.' 'o' 'p' 'o' 'p' 'o' 'p' 'o' 'p' 'o' 'p' \
    --outputpath_data $basepath_work_data \
    --pos 0 0.195 0.195 0.220 0.220 0.251 0.251 0.293 0.293 \
          0 0.195 0.195 0.220 0.220 0.251 0.251 0.293 0.293 0.352 0.352 \
    --Tc_in_GeV 0.180
}


set_labels_and_models(){
    labels_and_models=()
    labels_and_models+=("test" "dummy")
    temps=(195 220 251 293) #352
    for j in "${!temps[@]}" ; do
        temp=${temps[j]}

        labels_and_models+=(
    '$\text{smax}_\text{NLO}$'                    "T${temp}/spf/smax_NLO_Nf3_T0.${temp}_min${minscale}_wopt_${nsamples}smpls_tauTgtr0.24_${input_suffix}"
    '$\text{smax}_\text{LO\hphantom{N}}$'         "T${temp}/spf/smax_LO_Nf3_T0.${temp}_min${minscale}_w1_${nsamples}smpls_tauTgtr0.24_${input_suffix}"
        )
    done


    conftypes=(s096t36_b0824900_m002022_m01011 s096t32_b0824900_m002022_m01011 s096t28_b0824900_m002022_m01011 s096t24_b0824900_m002022_m01011 s096t20_b0824900_m002022_m01011)
    temps=(195 220 251 293 352)
    labels_and_models+=("test2" "dummy")
    for j in "${!conftypes[@]}" ; do
        conftype=${conftypes[j]}
        temp=${temps[j]}

        if [ "${conftype}" == s096t20_b0824900_m002022_m01011 ] ; then
            mintauT=0.26
            input_suffix="23-02-16_0.30"
            nsamples=1000
        else
            mintauT=0.24
            input_suffix="23-02-16_0.30"
            nsamples=1000
        fi

        labels_and_models+=(
    '$\text{smax}_\text{NLO}$'                    "${conftype}/spf/smax_NLO_Nf3_T0.${temp}_min${minscale}_wopt_${nsamples}smpls_tauTgtr${mintauT}_${input_suffix}"
    '$\text{smax}_\text{LO\hphantom{N}}$'         "${conftype}/spf/smax_LO_Nf3_T0.${temp}_min${minscale}_w1_${nsamples}smpls_tauTgtr${mintauT}_${input_suffix}"
        )
    done

}

set_labels_and_models
models=()
labels=()
for i in "${!labels_and_models[@]}" ; do
    if [ $((i % 2)) -eq 0 ] ; then
        labels+=("${labels_and_models[i]}")
    else
        models+=("${labels_and_models[i]}")
    fi
done


(
    cd "$(dirname $0)" || exit
    plot_kfactors
)