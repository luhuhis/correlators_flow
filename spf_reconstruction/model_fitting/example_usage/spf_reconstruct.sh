#!/bin/bash
source ~/.bashrc



minscale=eff

set_UV_params_EE(){
    LO="--PhiUV_order LO --omega_prefactor 1 --min_scale $minscale"
    NLO="--PhiUV_order NLO --omega_prefactor opt --min_scale $minscale"
}

set_models_EE(){
    models=(
        "--model max $LO"
        "--model max $NLO"
        "--model smax $LO"
        "--model smax $NLO"
        "--model plaw_any $LO --OmegaByT_UV $OmegaByT_UV"
        "--model plaw_any $NLO --OmegaByT_UV $OmegaByT_UV"
        "--model trig $LO --mu alpha --nmax 1 "
        "--model trig $LO --mu beta --nmax 1 "
        "--model trig $NLO --mu alpha --nmax 1 "
        "--model trig $NLO --mu beta --nmax 1 "
        "--model trig $LO --mu alpha --nmax 2 "
        "--model trig $LO --mu beta --nmax 2 "
        "--model trig $NLO --mu alpha --nmax 2 "
        "--model trig $NLO --mu beta --nmax 2 "
    )
}

submit_quenched_EE(){

    OmegaByT_UV=3.1416

    set_UV_params_EE
    set_models_EE

    nsamples=64
    for i in "${!models[@]}" ; do
         spfbatch ../spf_reconstruct.py \
            --output_path /work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/EE/spf/ \
            --add_suffix 23-01-30 \
            --input_corr /work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/EE/EE_flow_extr.npy \
            --min_tauT 0.24 \
            --nproc 64 \
            --T_in_GeV 0.472 \
            --Nf 0 \
            --nsamples $nsamples \
            ${models[i]}
    done
}

submit_quenched_BB(){

    OmegaByT_UV=3.1416

    LO="--PhiUV_order LO --omega_prefactor 1 --min_scale $minscale"
    NLO="--PhiUV_order LO --omega_prefactor optBB --min_scale $minscale"

    models=(
        "--model max $LO"
        "--model max $NLO"
        "--model smax $LO"
        "--model smax $NLO"
        "--model plaw_any $LO --OmegaByT_UV $OmegaByT_UV"
        "--model plaw_any $NLO --OmegaByT_UV $OmegaByT_UV"
        "--model trig $LO --mu alpha --nmax 1 "
        "--model trig $LO --mu beta --nmax 1 "
        "--model trig $NLO --mu alpha --nmax 1 "
        "--model trig $NLO --mu beta --nmax 1 "
        "--model trig $LO --mu alpha --nmax 2 "
        "--model trig $LO --mu beta --nmax 2 "
        "--model trig $NLO --mu alpha --nmax 2 "
        "--model trig $NLO --mu beta --nmax 2 "
    )

    nsamples=10000
    for i in "${!models[@]}" ; do
         spfbatch ../spf_reconstruct.py \
            --output_path /work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/BB/spf/ \
            --add_suffix 23-01-30 \
            --input_corr /work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/BB/BB_flow_extr.npy \
            --min_tauT 0.24 \
            --nproc 128 \
            --T_in_GeV 0.472 \
            --Nf 0 \
            --nsamples $nsamples \
            ${models[i]}
    done
}


submit_hisq(){

    OmegaByT_UV=6.2832

    temps=(195 220 251 293)
    int_nts=(36 32 28 24)
    relflowsuffix="_relflow"
    set_UV_params_EE

    models=(
#        "--model max $LO"
#        "--model max $NLO"
#        "--model smax $LO"
#        "--model smax $NLO"
        "--model plaw_any $LO --OmegaByT_UV $OmegaByT_UV" # --prevent_overfitting 1
        "--model plaw_any $NLO --OmegaByT_UV $OmegaByT_UV" #--prevent_overfitting 1
#        "--model plaw $LO --OmegaByT_IR 1 --OmegaByT_UV $OmegaByT_UV"
#        "--model plaw $NLO --OmegaByT_IR 1 --OmegaByT_UV $OmegaByT_UV"
        "--model trig $LO --mu alpha --nmax 1 --prevent_overfitting 2"
        "--model trig $NLO --mu beta --nmax 1 --prevent_overfitting 2"
        "--model trig $LO --mu alpha --nmax 2 --prevent_overfitting 2"
        "--model trig $NLO --mu beta --nmax 2 --prevent_overfitting 2"
    )

    nsamples=1000

    for j in "${!temps[@]}" ; do
#        if [ $j -eq 0 ] ; then
        temp=${temps[j]}
        int_nt=${int_nts[j]}
        for i in "${!models[@]}" ; do
#            if [ $i -eq 0 ] ; then
             spfbatch ../spf_reconstruct.py \
                --output_path /work/home/altenkort/work/correlators_flow/data/merged/hisq_ms5_zeuthenFlow/EE//T${temp}/spf/ \
                --add_suffix 23-02-16_relflow \
                --input_corr /work/home/altenkort/work/correlators_flow/data/merged/hisq_ms5_zeuthenFlow/EE//T${temp}/EE_flow_extr${relflowsuffix}.npy \
                --min_tauT 0.24 \
                --nproc 128 \
                --T_in_GeV 0.${temp} --corr_from_combined_fit_nt $int_nt \
                --Nf 3 \
                --nsamples $nsamples \
                ${models[i]}
#            fi
        done
#        fi
    done
}

submit_hisq_finite_a_and_tf(){

    OmegaByT_UV=6.2832

    Nts=(36 32 28 24 20)  #
    temps=(195 220 251 293 352)  #
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

    for j in "${!Nts[@]}" ; do
        Nt=${Nts[j]}
        temp=${temps[j]}
        for i in "${!models[@]}" ; do
             spfbatch ../spf_reconstruct.py \
                --output_path /work/home/altenkort/work/correlators_flow/data/merged/hisq_ms5_zeuthenFlow/EE//s096t${Nt}_b0824900_m002022_m01011/spf/ \
                --add_suffix 23-02-16_0.30 --relflow 0.30 \
                --input_corr /work/home/altenkort/work/correlators_flow/data/merged/hisq_ms5_zeuthenFlow/EE//s096t${Nt}_b0824900_m002022_m01011/EE_s096t${Nt}_b0824900_m002022_m01011_interpolation_relflow_samples.npy \
                --min_tauT 0.24 \
                --relflow_file /work/home/altenkort/work/correlators_flow/data/merged/hisq_ms5_zeuthenFlow/EE//s096t${Nt}_b0824900_m002022_m01011/EE_s096t${Nt}_b0824900_m002022_m01011_relflows.txt \
                --nproc 128 \
                --T_in_GeV 0.${temp} \
                --Nf 3 \
                --nsamples $nsamples \
                ${models[i]}
        done
    done
}


#submit_quenched_EE
#submit_quenched_BB
submit_hisq
#submit_hisq_finite_a_and_tf