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

    nsamples=500
    for i in "${!models[@]}" ; do
         spfbatch ../spf_reconstruct.py \
            --output_path /work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/EE/spf/ \
            --add_suffix 23-01-15 \
            --input_corr /work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/EE/EE_flow_extr.npy \
            --min_tauT 0.24 \
            --nproc 4 \
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

    nsamples=500
    for i in "${!models[@]}" ; do
         spfbatch ../spf_reconstruct.py \
            --output_path /work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/BB/spf/ \
            --add_suffix 23-01-19 \
            --input_corr /work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/BB/BB_flow_extr.npy \
            --min_tauT 0.24 \
            --nproc 4 \
            --T_in_GeV 0.472 \
            --Nf 0 \
            --nsamples $nsamples \
            ${models[i]}
    done
}


submit_hisq(){

    OmegaByT_UV=6.2832

    temps=(196 220 251 296)
    set_UV_params_EE

    models=(
        "--model max $LO"
        "--model max $NLO"
        "--model smax $LO"
        "--model smax $NLO"
        "--model plaw_any $LO --OmegaByT_UV $OmegaByT_UV"
        "--model plaw_any $NLO --OmegaByT_UV $OmegaByT_UV"
        "--model plaw $LO --OmegaByT_IR 1 --OmegaByT_UV $OmegaByT_UV"
        "--model plaw $NLO --OmegaByT_IR 1 --OmegaByT_UV $OmegaByT_UV"
        "--model trig $LO --mu alpha --nmax 1 --prevent_overfitting 1"
        "--model trig $NLO --mu beta --nmax 1 --prevent_overfitting 1"
        "--model trig $LO --mu alpha --nmax 2 --prevent_overfitting 1"
        "--model trig $NLO --mu beta --nmax 2 --prevent_overfitting 1"
    )

    nsamples=500

    for j in "${!temps[@]}" ; do
        temp=${temps[j]}
        for i in "${!models[@]}" ; do
             spfbatch ../spf_reconstruct.py \
                --output_path /work/home/altenkort/work/correlators_flow/data/merged/hisq_ms5_zeuthenFlow/EE//T${temp}/spf/ \
                --add_suffix 23-01-19 \
                --input_corr /work/home/altenkort/work/correlators_flow/data/merged/hisq_ms5_zeuthenFlow/EE//T${temp}/EE_flow_extr.npy \
                --min_tauT 0.24 \
                --nproc 4 \
                --T_in_GeV 0.${temp} \
                --Nf 3 \
                --nsamples $nsamples \
                ${models[i]}
        done
    done
}

submit_hisq_finite_a_and_tf(){

    OmegaByT_UV=6.2832

    Nts=(36 32 28 24 20)
    temps=(196 220 251 296 352)
    set_UV_params_EE

    models=(
        "--model max $LO"
        "--model max $NLO"
        "--model smax $LO"
        "--model smax $NLO"
        "--model plaw_any $LO --OmegaByT_UV $OmegaByT_UV"
        "--model plaw_any $NLO --OmegaByT_UV $OmegaByT_UV"
        "--model plaw $LO --OmegaByT_IR 1 --OmegaByT_UV $OmegaByT_UV"
        "--model plaw $NLO --OmegaByT_IR 1 --OmegaByT_UV $OmegaByT_UV"
        "--model trig $LO --mu alpha --nmax 1 --prevent_overfitting 1"
        "--model trig $NLO --mu beta --nmax 1 --prevent_overfitting 1"
        "--model trig $LO --mu alpha --nmax 2 --prevent_overfitting 1"
        "--model trig $NLO --mu beta --nmax 2 --prevent_overfitting 1"
    )

    nsamples=500

    for j in "${!Nts[@]}" ; do
        Nt=${Nts[j]}
        temp=${temps[j]}
        for i in "${!models[@]}" ; do
             spfbatch ../spf_reconstruct.py \
                --output_path /work/home/altenkort/work/correlators_flow/data/merged/hisq_ms5_zeuthenFlow/EE//s096t${Nt}_b0824900_m002022_m01011/spf/ \
                --add_suffix 23-01-19 \
                --input_corr /work/home/altenkort/work/correlators_flow/data/merged/hisq_ms5_zeuthenFlow/EE//s096t${Nt}_b0824900_m002022_m01011/rel_flow/EE_relflow_0.30.dat \
                --min_tauT 0.24 \
                --mock_bootstrap \
                --nproc 4 \
                --T_in_GeV 0.${temp} \
                --Nf 3 \
                --nsamples $nsamples \
                ${models[i]}
        done
    done
}


#submit_quenched_EE
submit_quenched_BB
#submit_hisq
#submit_hisq_finite_a_and_tf