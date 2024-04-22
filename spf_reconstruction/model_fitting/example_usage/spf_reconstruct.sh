#!/bin/bash

# TODO add flag to choose between EE BB, hisq quenched etc

basepath_work_data="${1:-"/work/home/altenkort/work/correlators_flow/data/merged/"}"
use_cluster="${2:-"yes"}"
nproc="${3:-"1"}"

if [ "$use_cluster" == "yes" ]; then
    nproc=128
    spfbatch() {
        sbatch --qos=debug --partition=cpu_compute --time=10-00:00:00 --ntasks=1 --cpus-per-task=128 --nodes=1 "$@"
    }
else
    spfbatch() {
        export LD_LIBRARY_PATH=$HOME/libffi/lib64:$LD_LIBRARY_PATH
        "$@"
    }
fi

minscale=eff

set_UV_params_EE() {
    LO="--order LO --omega_prefactor 1 --min_scale $minscale"
    NLO="--order NLO --omega_prefactor opt --min_scale $minscale"
}

set_models_EE() {
    models=(
        "--model max $LO"
        "--model max $NLO"
        "--model smax $LO"
        "--model smax $NLO"
        "--model plaw_any $LO --OmegaByT_UV $OmegaByT_UV"
        "--model plaw_any $NLO --OmegaByT_UV $OmegaByT_UV"
        "--model plaw $LO --OmegaByT_IR 1 --OmegaByT_UV 6.2832"
        "--model plaw $NLO --OmegaByT_IR 1 --OmegaByT_UV 6.2832"
        "--model trig $LO --mu alpha --nmax 1 "
        "--model trig $NLO --mu beta --nmax 1 "
        "--model trig $LO --mu alpha --nmax 2 "
        "--model trig $NLO --mu beta --nmax 2 "
    )
}

submit_quenched_EE() {

    OmegaByT_UV=3.1416

    set_UV_params_EE
    set_models_EE

    nsamples=1000
    for i in "${!models[@]}"; do
        spfbatch ../spf_reconstruct.py \
            --output_path $basepath_work_data/quenched_1.50Tc_zeuthenFlow/EE/spf/ \
            --add_suffix 23-02-26_2piT \
            --input_corr $basepath_work_data/quenched_1.50Tc_zeuthenFlow/EE/EE_flow_extr_relflow.npy \
            --min_tauT 0.24 \
            --nproc $nproc \
            --T_in_GeV 0.472 \
            --corr EE \
            --Nf 0 \
            --nsamples $nsamples \
            ${models[i]}
    done
}

submit_quenched_BB() {

    #    OmegaByT_UV=3.1416
    input_corr_suffixes=(
        "ref6.28_UVLO_IRLO"
        "ref6.28_UVLO_IRNLO"
        "ref6.28_UVNLO_IRLO"
        "ref6.28_UVNLO_IRNLO"
    )

    nsamples=1000
    for input_corr_suffix in "${input_corr_suffixes[@]}"; do

        #Adjusting the minscale based on the suffix
        if [[ $input_corr_suffix == *"_IRLO" ]]; then
            mu_IR_by_T="LO" # replace with the desired value
        elif [[ $input_corr_suffix == *"_IRNLO" ]]; then
            mu_IR_by_T="NLO" # replace with the desired value
        fi

        minscale="eff"

        NLO_naive_soft="--order NLO --omega_prefactor 1 --min_scale $minscale --mu_IR_by_T $mu_IR_by_T --max_type smooth"

        models=(
            "--model max $NLO_naive_soft"

            "--model smax $NLO_naive_soft"

            "--model plaw $NLO_naive_soft --OmegaByT_IR 1 --OmegaByT_UV 6.2832 "
        )

        for i in "${!models[@]}"; do
            echo "============================================================================================"
            echo "${models[i]}"
            spfbatch ../spf_reconstruct.py \
                --output_path $basepath_work_data/quenched_1.50Tc_zeuthenFlow/BB/spf/ \
                --add_suffix 24-02-08-${input_corr_suffix} \
                --input_corr $basepath_work_data/quenched_1.50Tc_zeuthenFlow/BB/BB_flow_extr_relflow_${input_corr_suffix}.npy \
                --min_tauT 0.24 \
                --nproc $nproc \
                --T_in_GeV 0.472 \
                --corr BB \
                --Nf 0 \
                --nsamples $nsamples \
                ${models[i]}
        done
    done
}

submit_hisq() {

    OmegaByT_UV=6.2832

    temps=(195 220 251 293)
    int_nts=(36 32 28 24)
    relflowsuffix="_relflow"
    set_UV_params_EE

    models=(
        "--model max $LO"
        "--model max $NLO"
        "--model smax $LO"
        "--model smax $NLO"
        "--model plaw $LO --OmegaByT_IR 1 --OmegaByT_UV $OmegaByT_UV"
        "--model plaw $NLO --OmegaByT_IR 1 --OmegaByT_UV $OmegaByT_UV"
    )

    nsamples=1000

    for j in "${!temps[@]}"; do
        #        if [ $j -eq 0 ] ; then
        temp=${temps[j]}
        int_nt=${int_nts[j]}
        for i in "${!models[@]}"; do
            #            if [ $i -eq 0 ] ; then
            spfbatch ../spf_reconstruct.py \
                --output_path $basepath_work_data/hisq_ms5_zeuthenFlow/EE//T${temp}/spf/ \
                --add_suffix 23-02-16_relflow \
                --input_corr $basepath_work_data/hisq_ms5_zeuthenFlow/EE//T${temp}/EE_flow_extr${relflowsuffix}.npy \
                --min_tauT 0.24 \
                --nproc $nproc \
                --T_in_GeV 0.${temp} --corr_from_combined_fit_nt $int_nt \
                --corr EE \
                --Nf 3 \
                --nsamples $nsamples \
                ${models[i]}
            #            fi
        done
        #        fi
    done
}

submit_hisq_finite_a_and_tf() {

    OmegaByT_UV=6.2832

    Nts=(36 32 28 24 20)        #
    temps=(195 220 251 293 352) #
    set_UV_params_EE

    models=(
        "--model max $LO"
        "--model max $NLO"
        "--model smax $LO"
        "--model smax $NLO"
        "--model plaw $LO --OmegaByT_IR 1 --OmegaByT_UV $OmegaByT_UV"
        "--model plaw $NLO --OmegaByT_IR 1 --OmegaByT_UV $OmegaByT_UV"
    )

    nsamples=1000

    for j in "${!Nts[@]}"; do
        Nt=${Nts[j]}
        temp=${temps[j]}
        if [ "${temp}" == "352" ]; then
            mintauT=0.26
        else
            mintauT=0.24
        fi
        for i in "${!models[@]}"; do
            spfbatch ../spf_reconstruct.py \
                --output_path $basepath_work_data/hisq_ms5_zeuthenFlow/EE//s096t${Nt}_b0824900_m002022_m01011/spf/ \
                --add_suffix 23-02-16_0.30 --relflow 0.30 \
                --input_corr $basepath_work_data/hisq_ms5_zeuthenFlow/EE//s096t${Nt}_b0824900_m002022_m01011/EE_s096t${Nt}_b0824900_m002022_m01011_interpolation_relflow_samples.npy \
                --min_tauT $mintauT \
                --relflow_file $basepath_work_data/hisq_ms5_zeuthenFlow/EE//s096t${Nt}_b0824900_m002022_m01011/EE_s096t${Nt}_b0824900_m002022_m01011_relflows.txt \
                --nproc $nproc \
                --corr EE \
                --T_in_GeV 0.${temp} \
                --Nf 3 \
                --nsamples $nsamples \
                ${models[i]}
        done
    done
}

(
    cd "$(dirname $0)" || exit

    #    submit_quenched_EE
    submit_quenched_BB
    #    submit_hisq
    #    submit_hisq_finite_a_and_tf

)
