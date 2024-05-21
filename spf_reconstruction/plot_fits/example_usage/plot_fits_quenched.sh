#!/bin/bash

corr=$1
basepath_work_data=$2
basepath_plot=$3
plot_all_renorm_points=$4

if [ -z "$corr" ] || [ -z "$basepath_work_data" ] || [ -z "$basepath_plot" ] || [ -z "$plot_all_renorm_points" ] ; then
    echo "Usage: $0 corr basepath_work_data basepath_plot"
    exit
fi


if [ "${corr}" == "EE" ] ; then
    input_suffix="23-02-26_2piT"
    nsamples=1000
    minscale=eff
    ylims="1 200"
    add_suffix=""
    xlims="--xlims -0.15 4.15 --xticks 0 1 2 3 4"
    figsize="7 7"

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

    nsamples=1000
    minscale=eff
    ylims="0.1 5000"
    add_suffix=""
    xlims="--xlims -0.15 2.5    --xticks 0 0.5 1 1.5 2 2.5 --hide_chisq"
    figsize="7 5"

    if [ "${plot_all_renorm_points}" == "yes" ] ; then
        input_corr_suffixes=(
        "ref6.28_UVLO_IRLO"
        "ref6.28_UVNLO_IRLO"
        "ref6.28_UVLO_IRNLO"
        "ref6.28_UVNLO_IRNLO"
        )
    else
        input_corr_suffixes=(
        "ref6.28_UVNLO_IRNLO"
        )
    fi


    for input_corr_suffix in "${input_corr_suffixes[@]}" ; do

        input_suffix="24-02-08-${input_corr_suffix}"
        minscale="eff"
        prefactor="1"

        labels_and_models=(
            '$\mathrm{max}$'                   "max_NLO_Nf0_T0.472_min${minscale}_w${prefactor}_${nsamples}smpls_tauTgtr0.24_${input_suffix}"
            '$\mathrm{smax}$'                  "smax_NLO_Nf0_T0.472_min${minscale}_w${prefactor}_${nsamples}smpls_tauTgtr0.24_${input_suffix}"
            '$\mathrm{plaw}$'                  "plaw_wIR1.0_wUV6.2832_NLO_Nf0_T0.472_min${minscale}_w${prefactor}_${nsamples}smpls_tauTgtr0.24_${input_suffix}"
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

    if [ "${corr}" == "BB" ] && [ "${plot_all_renorm_points}" == "yes" ] ; then
        plot_kappa
    elif [ "${corr}" == "BB" ] && [ "${plot_all_renorm_points}" != "yes" ] ; then
        plot_spfs
        plot_fitcorr
    else
        plot_spfs
        plot_kappa
        plot_fitcorr
    fi
)



