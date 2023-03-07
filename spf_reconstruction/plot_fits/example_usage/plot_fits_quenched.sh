#!/bin/bash


# quenched 1.5 Tc

for corr in "EE" "BB" ; do # "EE"

    if [ "${corr}" == "EE" ] ; then
        input_suffix="23-02-26"
        nsamples=1000
        minscale=eff
        ylims="1 200"
        add_suffix=""

        labels_and_models=(
            '$\text{max}_\text{NLO}$'                   "max_NLO_Nf0_T0.472_min${minscale}_wopt_${nsamples}smpls_tauTgtr0.24_${input_suffix}"
            '$\text{max}_\text{LO}$'                    "max_LO_Nf0_T0.472_min${minscale}_w1_${nsamples}smpls_tauTgtr0.24_${input_suffix}"
            '$\text{smax}_\text{NLO}$'                  "smax_NLO_Nf0_T0.472_min${minscale}_wopt_${nsamples}smpls_tauTgtr0.24_${input_suffix}"
            '$\text{smax}_\text{LO}$'                   "smax_LO_Nf0_T0.472_min${minscale}_w1_${nsamples}smpls_tauTgtr0.24_${input_suffix}"
            '$\text{plaw}_\text{NLO}$'                  "plaw_wIR1.0_wUV6.2832_NLO_Nf0_T0.472_min${minscale}_wopt_${nsamples}smpls_tauTgtr0.24_${input_suffix}"
            '$\text{plaw}_\text{LO}$'                   "plaw_wIR1.0_wUV6.2832_LO_Nf0_T0.472_min${minscale}_w1_${nsamples}smpls_tauTgtr0.24_${input_suffix}"
            '$\text{trig1}_\text{NLO}$'                 "trig_NLO__beta_1_Nf0_T0.472_min${minscale}_wopt_${nsamples}smpls_tauTgtr0.24_${input_suffix}"
            '$\text{trig1}_\text{LO}$'                  "trig_LO__alpha_1_Nf0_T0.472_min${minscale}_w1_${nsamples}smpls_tauTgtr0.24_${input_suffix}"
            '$\text{trig2}_\text{NLO}$'                 "trig_NLO__beta_2_Nf0_T0.472_min${minscale}_wopt_${nsamples}smpls_tauTgtr0.24_${input_suffix}"
            '$\text{trig2}_\text{LO}$'                  "trig_LO__alpha_2_Nf0_T0.472_min${minscale}_w1_${nsamples}smpls_tauTgtr0.24_${input_suffix}"
#            '$\text{trig1}_\text{NLO}$'                 "trig_NLO__alpha_1_Nf0_T0.472_min${minscale}_wopt_${nsamples}smpls_tauTgtr0.24_${input_suffix}"
#            '$\text{trig1}_\text{LO}$'      "trig_LO__beta_1_Nf0_T0.472_min${minscale}_w1_${nsamples}smpls_tauTgtr0.24_${input_suffix}"
#            '$\text{trig2}_\text{NLO}$'                 "trig_NLO__alpha_2_Nf0_T0.472_min${minscale}_wopt_${nsamples}smpls_tauTgtr0.24_${input_suffix}"
#            '$\text{trig2}_\text{LO}$'      "trig_LO__beta_2_Nf0_T0.472_min${minscale}_w1_${nsamples}smpls_tauTgtr0.24_${input_suffix}"
        )

    elif [ "${corr}" == "BB" ] ; then
        input_suffix="23-01-19"
        nsamples=500
        minscale=eff
        ylims="0.1 5000"
        add_suffix=""

        labels_and_models=(
            '$\text{max}_\text{NLO}$'                   "max_LO_Nf0_T0.472_min${minscale}_woptBB_${nsamples}smpls_tauTgtr0.24_${input_suffix}"
            '$\text{max}_\text{LO}$'                    "max_LO_Nf0_T0.472_min${minscale}_w1_${nsamples}smpls_tauTgtr0.24_${input_suffix}"
            '$\text{smax}_\text{NLO}$'                  "smax_LO_Nf0_T0.472_min${minscale}_woptBB_${nsamples}smpls_tauTgtr0.24_${input_suffix}"
            '$\text{smax}_\text{LO}$'                   "smax_LO_Nf0_T0.472_min${minscale}_w1_${nsamples}smpls_tauTgtr0.24_${input_suffix}"
            '$\text{plaw}_\text{NLO}$'                  "plaw_any_wUV3.1416_LO_Nf0_T0.472_min${minscale}_woptBB_${nsamples}smpls_tauTgtr0.24_${input_suffix}"
            '$\text{plaw}_\text{LO}$'                   "plaw_any_wUV3.1416_LO_Nf0_T0.472_min${minscale}_w1_${nsamples}smpls_tauTgtr0.24_${input_suffix}"
            '$\text{trig1}_\text{NLO}$'                 "trig_LO__beta_1_Nf0_T0.472_min${minscale}_woptBB_${nsamples}smpls_tauTgtr0.24_${input_suffix}"
            '$\text{trig1}_\text{LO}$'                  "trig_LO__alpha_1_Nf0_T0.472_min${minscale}_w1_${nsamples}smpls_tauTgtr0.24_${input_suffix}"
            '$\text{trig2}_\text{NLO}$'                 "trig_LO__beta_2_Nf0_T0.472_min${minscale}_woptBB_${nsamples}smpls_tauTgtr0.24_${input_suffix}"
            '$\text{trig2}_\text{LO}$'                  "trig_LO__alpha_2_Nf0_T0.472_min${minscale}_w1_${nsamples}smpls_tauTgtr0.24_${input_suffix}"
        )

    fi

#    lowsamples=100



    models=()
    labels=()
    for i in "${!labels_and_models[@]}" ; do
        if [ $((i % 2)) -eq 0 ] ; then
            labels+=("${labels_and_models[i]}")
        else
            models+=("${labels_and_models[i]}")
        fi
    done


    plot_kappa() {
        ../plot_kappa.py \
    --xlims -0.15 4.15 \
    --xticks 0 1 2 3 4 \
    --model_ids "${models[@]}" \
    --labels "${labels[@]}" \
    --basepath ~/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/${corr}/spf/ \
    --outputpath ~/work/correlators_flow/plots/quenched_1.50Tc_zeuthenFlow/${corr}/ \
    --outputpath_data ~/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/${corr}/ \
    --suffix _quenched_1.5Tc${add_suffix} \
    --corr ${corr} \
    --figsize 7 7
    }

    plot_fitcorr() {
        ../plot_fitcorr.py \
    --model_ids "${models[@]}" \
    --labels "${labels[@]}" \
    --basepath ~/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/${corr}/spf/ \
    --outputpath ~/work/correlators_flow/plots/quenched_1.50Tc_zeuthenFlow/${corr}/ \
    --suffix _quenched_1.5Tc${add_suffix} \
    --corr ${corr} \
    --ylims 0.94 1.04 \
    --xticks 0.25 0.30 0.35 0.40 0.45 0.50 \
    --yticks 0.98 0.99 1 1.01 1.02
    }

    plot_spfs(){
        ../plot_spfs.py \
    --model_ids "${models[@]}" \
    --labels "${labels[@]}" \
    --basepath ~/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/${corr}/spf/ \
    --outputpath ~/work/correlators_flow/plots/quenched_1.50Tc_zeuthenFlow/${corr}/ \
    --suffix _quenched_1.5Tc${add_suffix} \
    --corr ${corr} \
    --ylims ${ylims}
    }

    plot_spfs &
    plot_kappa &
    plot_fitcorr &

done
wait
