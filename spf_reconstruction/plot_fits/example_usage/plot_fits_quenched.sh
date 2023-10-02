#!/bin/bash


# quenched 1.5 Tc

for corr in "BB" ; do # "EE"

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
#            '$\text{trig1}_\text{NLO}$'                 "trig_NLO__beta_1_Nf0_T0.472_min${minscale}_wopt_${nsamples}smpls_tauTgtr0.24_${input_suffix}"
#            '$\text{trig1}_\text{LO}$'                  "trig_LO__alpha_1_Nf0_T0.472_min${minscale}_w1_${nsamples}smpls_tauTgtr0.24_${input_suffix}"
#            '$\text{trig2}_\text{NLO}$'                 "trig_NLO__beta_2_Nf0_T0.472_min${minscale}_wopt_${nsamples}smpls_tauTgtr0.24_${input_suffix}"
#            '$\text{trig2}_\text{LO}$'                  "trig_LO__alpha_2_Nf0_T0.472_min${minscale}_w1_${nsamples}smpls_tauTgtr0.24_${input_suffix}"
        )

        for i in "${!labels_and_models[@]}" ; do
          if [ $((i % 2)) -eq 0 ] ; then
              labels+=("${labels_and_models[i]}")
          else
              models+=("${labels_and_models[i]}")
          fi
       done

    fi



    if [ "${corr}" == "BB" ] ; then

    models=()
    labels=()

      # "23-01-19"
    nsamples=100
    minscale=eff
    ylims="0.1 5000"
    add_suffix=""

    input_corr_suffixes=(
    "ref4.0_UVLO_IRNLO"
    "ref4.0_UVLO_IRLO"
    "ref4.0_UVNLO_IRNLO"
    "ref4.0_UVNLO_IRLO"
    )


    for input_corr_suffix in "${input_corr_suffixes[@]}" ; do

        # 08-27
        # 10-01
        input_suffix="23-10-01-${input_corr_suffix}"

        labels_and_models=(                              #max_LO_Nf0_T0.472_mineff_w1_250smpls_tauTgtr0.24_23-07-03
            '$\text{max}_\text{NLO}$'                   "max_NLO_Nf0_T0.472_min${minscale}_woptBB_${nsamples}smpls_tauTgtr0.24_${input_suffix}"
            '$\text{max}_\text{LO}$'                    "max_LO_Nf0_T0.472_min${minscale}_w1_${nsamples}smpls_tauTgtr0.24_${input_suffix}"
            '$\text{smax}_\text{NLO}$'                  "smax_NLO_Nf0_T0.472_min${minscale}_woptBB_${nsamples}smpls_tauTgtr0.24_${input_suffix}"
            '$\text{smax}_\text{LO}$'                   "smax_LO_Nf0_T0.472_min${minscale}_w1_${nsamples}smpls_tauTgtr0.24_${input_suffix}"
            '$\text{plaw}_\text{NLO}$'                  "plaw_wIR1.0_wUV6.2832_NLO_Nf0_T0.472_min${minscale}_woptBB_${nsamples}smpls_tauTgtr0.24_${input_suffix}"
            '$\text{plaw}_\text{LO}$'                   "plaw_wIR1.0_wUV6.2832_LO_Nf0_T0.472_min${minscale}_w1_${nsamples}smpls_tauTgtr0.24_${input_suffix}"
#            '$\text{trig1}_\text{NLO}$'                 "trig_NLO__beta_1_Nf0_T0.472_min${minscale}_woptBB_${nsamples}smpls_tauTgtr0.24_${input_suffix}"
#            '$\text{trig1}_\text{LO}$'                  "trig_LO__alpha_1_Nf0_T0.472_min${minscale}_w1_${nsamples}smpls_tauTgtr0.24_${input_suffix}"
#            '$\text{trig2}_\text{NLO}$'                 "trig_NLO__beta_2_Nf0_T0.472_min${minscale}_woptBB_${nsamples}smpls_tauTgtr0.24_${input_suffix}"
#            '$\text{trig2}_\text{LO}$'                  "trig_LO__alpha_2_Nf0_T0.472_min${minscale}_w1_${nsamples}smpls_tauTgtr0.24_${input_suffix}"
        )

      for i in "${!labels_and_models[@]}" ; do
          if [ $((i % 2)) -eq 0 ] ; then
              labels+=("${labels_and_models[i]}")
          else
              models+=("${labels_and_models[i]}")
          fi
      done
    done

    fi

#    lowsamples=100






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
    --figsize 9 9
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

#    plot_spfs
    plot_kappa
#    plot_fitcorr

done
wait
