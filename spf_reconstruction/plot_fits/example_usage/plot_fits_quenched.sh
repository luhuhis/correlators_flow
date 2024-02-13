#!/bin/bash

corr=$1
basepath_work_data=$2
basepath_plot=$3


if [ -z "$corr" ] || [ -z "$basepath_work_data" ] || [ -z "$basepath_plot" ]; then
    echo "Usage: $0 corr basepath_work_data basepath_plot"
    echo "Example usage: $0 BB ~/work/correlators_flow/data/merged/ ~/work/correlators_flow/plots/"
    exit
fi


if [ "${corr}" == "EE" ] ; then
    input_suffix="23-02-26"
    nsamples=1000
    minscale=eff
    ylims="1 200"
    add_suffix=""
    xlims="--xlims -0.15 4.15 --xticks 0 1 2 3 4"
    figsize="9 9"

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

elif [ "${corr}" == "BB" ] ; then

    models=()
    labels=()

      # "23-01-19"
    nsamples=1000
    minscale=eff
    ylims="0.1 5000"
    add_suffix=""
    xlims="--xlims -0.15 2.5    --xticks 0 0.5 1 1.5 2 2.5 --hide_chisq"
    figsize="7 5"

    input_corr_suffixes=(
    "ref6.28_UVLO_IRLO"
    "ref6.28_UVNLO_IRLO"
    "ref6.28_UVLO_IRNLO"
    "ref6.28_UVNLO_IRNLO"
    )


    for input_corr_suffix in "${input_corr_suffixes[@]}" ; do

        # 08-27
        # 10-01
        # 10-23
        input_suffix="24-02-08-${input_corr_suffix}"

        # smax_NLO_Nf0_T0.472_minmu_IR_NLO_w2_100smpls_tauTgtr0.24_23-10-23-ref4.0_UVLO_IRNLO

        # Adjusting the minscale based on the suffix
    #    if [[ $input_corr_suffix == *"_IRLO" ]]; then
    #        minscale="2piT" # replace with the desired value
    #    elif [[ $input_corr_suffix == *"_IRNLO" ]]; then
    #        minscale="mu_IR_NLO" # replace with the desired value
    #    fi                                                                                       mineff_wBB_NLO_full

        minscale="eff"

        prefactor="1"

        labels_and_models=(                              #max_LO_Nf0_T0.472_mineff_w1_250smpls_tauTgtr0.24_23-07-03
            '$\mathrm{max}$'                   "max_NLO_Nf0_T0.472_min${minscale}_w${prefactor}_${nsamples}smpls_tauTgtr0.24_${input_suffix}"
#            '$\text{max}_\text{LO}$'                    "max_LO_Nf0_T0.472_min${minscale}_w1_${nsamples}smpls_tauTgtr0.24_${input_suffix}"
            '$\mathrm{smax}$'                  "smax_NLO_Nf0_T0.472_min${minscale}_w${prefactor}_${nsamples}smpls_tauTgtr0.24_${input_suffix}"
                                                       # smax_NLO_Nf0_T0.472_minmu_IR_NLO_w2_100smpls_tauTgtr0.24_23-10-23-ref4.0_UVLO_IRNLO
#            '$\text{smax}_\text{LO}$'                   "smax_LO_Nf0_T0.472_min${minscale}_w1_${nsamples}smpls_tauTgtr0.24_${input_suffix}"
            '$\mathrm{plaw}$'                  "plaw_wIR1.0_wUV6.2832_NLO_Nf0_T0.472_min${minscale}_w${prefactor}_${nsamples}smpls_tauTgtr0.24_${input_suffix}"
#            '$\text{plaw}_\text{LO}$'                   "plaw_wIR1.0_wUV6.2832_LO_Nf0_T0.472_min${minscale}_w1_${nsamples}smpls_tauTgtr0.24_${input_suffix}"
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


plot_kappa() {

    ../plot_kappa.py \
$xlims \
--model_ids "${models[@]}" \
--labels "${labels[@]}" \
--basepath $basepath_work_data/quenched_1.50Tc_zeuthenFlow/${corr}/spf/ \
--outputpath $basepath_plot/quenched_1.50Tc_zeuthenFlow/${corr}/ \
--outputpath_data $basepath_work_data/quenched_1.50Tc_zeuthenFlow/${corr}/ \
--suffix _quenched_1.5Tc${add_suffix} \
--corr ${corr} \
--figsize $figsize
}

plot_fitcorr() {
    ../plot_fitcorr.py \
--model_ids "${models[@]}" \
--labels "${labels[@]}" \
--basepath $basepath_work_data/quenched_1.50Tc_zeuthenFlow/${corr}/spf/ \
--outputpath $basepath_plot/quenched_1.50Tc_zeuthenFlow/${corr}/ \
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
--basepath $basepath_work_data/quenched_1.50Tc_zeuthenFlow/${corr}/spf/ \
--outputpath $basepath_plot/quenched_1.50Tc_zeuthenFlow/${corr}/ \
--suffix _quenched_1.5Tc${add_suffix} \
--corr ${corr} \
--ylims ${ylims}
}


(
    cd "$(dirname $0)" || exit
    plot_spfs
    plot_kappa
    plot_fitcorr
)



