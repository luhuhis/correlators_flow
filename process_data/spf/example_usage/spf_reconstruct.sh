#!/bin/bash
source ~/.bashrc

#--mu \
#--nmax \
#--OmegaByT_IR \
#--OmegaByT_UV \
#     --p \

LO="--PhiUV_order LO --omega_prefactor 1 --min_scale piT --max_type hard"
NLO="--PhiUV_order NLO --omega_prefactor opt --min_scale piT --max_type hard"

models=(
    "--model max $LO"
#    "--model max $NLO"
#    "--model smax $LO"
#    "--model smax $NLO"
#    "--model line $LO --OmegaByT_IR 1 --OmegaByT_UV 3.14"
#    "--model line $NLO --OmegaByT_IR 1 --OmegaByT_UV 3.14"
#    "--model fourier $LO --mu alpha --nmax 1"
#    "--model fourier $NLO --mu beta --nmax 1"
#    "--model fourier $LO --mu alpha --nmax 2"
#    "--model fourier $NLO --mu beta --nmax 2"
#    "--model fourier $LO --mu alpha --nmax 3"
#    "--model fourier $NLO --mu beta --nmax 3"
#    "--model fourier $LO --mu alpha --nmax 4 --prevent_overfitting"
#    "--model fourier $NLO --mu beta --nmax 4 --prevent_overfitting"
#    "--model fourier $LO --mu alpha --nmax 5 --prevent_overfitting"
#    "--model fourier $NLO --mu beta --nmax 5 --prevent_overfitting"
)

for i in "${!models[@]}" ; do
#spfbatch
     spfbatch ../spf_reconstruct.py \
        --output_path /work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/EE/spf/ \
        --add_suffix final \
        --input_corr /work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/EE/EE_flow_extr.npy \
        --min_tauT 0.24 \
        --nproc 4 \
        --T_in_GeV 0.472 \
        --Nf 0 \
        --nsamples 100 \
        ${models[i]}
done