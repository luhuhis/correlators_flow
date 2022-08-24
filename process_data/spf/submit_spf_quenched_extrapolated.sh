#!/bin/bash
if [ 1 -eq 2 ] ; then
for model in "--model pnorm --p 2" "--model max" "--model sum" "--model plaw --OmegaByT_IR 1.0 --OmegaByT_UV 6.28" ; do
    for PhiUV in "--PhiUVtype LO  --min_scale 2piT --omega_prefactor 1" "--PhiUVtype NLO --min_scale opt --omega_prefactor opt" ; do  #
        ./spf_reconstruct.py --add_suffix quenched_cont_f0 --nsamples 500 --nproc 20 \
                    --output_path /home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/EE/ \
                    --input_corr /home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/EE/EE_flow_extr.txt \
                    --min_tauT 0 $model \
                    --max_type smooth --T_in_GeV 0.472 --Nf 0 $PhiUV
    done
done
fi

source ~/.bashrc
for model in "--model pnorm --p 2" "--model max" "--model sum" "--model plaw --OmegaByT_IR 1.0 --OmegaByT_UV 6.28" ; do
    for PhiUV in "--PhiUVtype LO  --min_scale 2piT --omega_prefactor 1" "--PhiUVtype NLO --min_scale opt --omega_prefactor opt" ; do  #
        spfbatch ./spf_reconstruct.py --add_suffix quenched_cont_f0 --nsamples 500 --nproc 20 \
                    --output_path /home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/EE/cont_rel_flow/rel_flow/ \
                    --input_corr /home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/EE/cont_rel_flow/rel_flow/EE_relflow_0.25.dat \
                    --min_tauT 0 $model \
                    --max_type smooth --T_in_GeV 0.472 --Nf 0 $PhiUV
    done
done