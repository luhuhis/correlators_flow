#!/bin/bash
source ~/.bashrc


set_UV_params(){
    LO="--PhiUV_order LO --omega_prefactor 1 --min_scale $minscale --max_type hard"
    NLO="--PhiUV_order NLO --omega_prefactor opt --min_scale $minscale --max_type hard"
}

set_models(){
    models=(
#        "--model max $LO"
#        "--model max $NLO"
#        "--model smax $LO"
#        "--model smax $NLO"
#        "--model plaw $LO --OmegaByT_IR 1 --OmegaByT_UV $OmegaByT_UV"
#        "--model plaw $NLO --OmegaByT_IR 1 --OmegaByT_UV $OmegaByT_UV"
        "--model trig $LO --mu alpha --nmax 1 --prevent_overfitting "
        "--model trig $NLO --mu beta --nmax 1 --prevent_overfitting"
        "--model trig $LO --mu alpha --nmax 2 --prevent_overfitting"
        "--model trig $NLO --mu beta --nmax 2 --prevent_overfitting"
#        "--model trig $LO --mu alpha --nmax 3 --prevent_overfitting"
#        "--model trig $NLO --mu beta --nmax 3 --prevent_overfitting"
#        "--model fourier $LO --mu alpha --nmax 4 --prevent_overfitting"
#        "--model fourier $NLO --mu beta --nmax 4 --prevent_overfitting"
#        "--model fourier $LO --mu alpha --nmax 5 --prevent_overfitting"
#        "--model fourier $NLO --mu beta --nmax 5 --prevent_overfitting"
    )
}

submit_quenched_EE(){
    minscale=piT
    OmegaByT_UV=3.14
    set_UV_params
    set_models

    nsamples=1000
    for i in "${!models[@]}" ; do
         spfbatch ../spf_reconstruct.py \
            --output_path /work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/EE/spf/ \
            --add_suffix final \
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
    minscale=piT
    OmegaByT_UV=3.14
    set_UV_params
    set_models

    nsamples=100
    for i in "${!models[@]}" ; do
         spfbatch ../spf_reconstruct.py \
            --output_path /work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/BB/spf/ \
            --add_suffix final \
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
    temps=( 196 220 251 296)
    minscale=2piT
    OmegaByT_UV=6.28
    set_UV_params
    set_models

    nsamples=1000

    for j in "${!temps[@]}" ; do
        temp=${temps[j]}
        for i in "${!models[@]}" ; do
             spfbatch ../spf_reconstruct.py \
                --output_path /work/home/altenkort/work/correlators_flow/data/merged/hisq_ms5_zeuthenFlow/EE//T${temp}/spf/ \
                --add_suffix paper \
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

#submit_quenched_EE
#submit_quenched_BB
submit_hisq