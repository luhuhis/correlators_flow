#!/bin/bash
source ~/.bashrc

###### SYSTEMATIC TESTS FOR HISQ

# LINE model


#
#spfbatch ./spf_reconstruct.py --add_suffix quenched_1.5Tc_conttauF0_piT_1w_wIR0.01 \
#--nsamples 200 --nproc 20 --output_path ~/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/EE/ \
#--input_corr ~/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/EE/EE_flow_extr.txt \
#--min_tauT 0.24 --error_exponent 0 \
#--model line --OmegaByT_IR 0.01 \
#--PhiUVtype LO \
#--T_in_GeV 0.472 --Nf 0 --min_scale piT --omega_prefactor 1 --max_type smooth


#--PathPhiUV ~/work/correlators_flow/data/merged/spf_coupling/quenched_1.5Tc_SPF_LO_Nf0_0.472_piT_1.0_smax.dat \
#quenched_1.5Tc_cont_tauF0_line0.01_piT_1.0w


# Nt=24 at 0.2 is done

Nts=( 36 32 28 24 20 )
temps=( 196 220 251 296 352 )
flows=( 0.25 0.20 )  #0.20 0.15
#OmegaByT_UVs=( 7.143 6.364 5.578 4.728 3.977 )

if [ 1 -eq 2 ] ; then
#    if [ ${Nts[idx]} == 24 ] ; then
for flow in "${flows[@]}" ; do
    for idx in "${!Nts[@]}"; do
        # ${OmegaByT_UVs[idx]}
        for model in "--model plaw --OmegaByT_IR 0.4 --OmegaByT_UV 6.28" ; do  #  "--model pnorm --p 2" "--model max"
            for PhiUV in "--PhiUVtype LO  --min_scale 2piT --omega_prefactor 1" "--PhiUVtype NLO --min_scale opt --omega_prefactor opt" ; do  #
                spfbatch ./spf_reconstruct.py --add_suffix hisq_nt${Nts[idx]}_f${flow} --nsamples 500 --nproc 4 \
                --output_path /home/altenkort/work/correlators_flow/data/merged/hisq_b8249_zeuthenFlow/EE/s096t${Nts[idx]}_b0824900_m002022_m01011/ \
                --input_corr /home/altenkort/work/correlators_flow/data/merged/hisq_b8249_zeuthenFlow/EE/s096t${Nts[idx]}_b0824900_m002022_m01011/rel_flow/EE_relflow_${flow}.dat \
                --min_tauT 0.35 $model \
                --max_type smooth --T_in_GeV 0.${temps[idx]} --Nf 3 $PhiUV
                done
            done
    #    fi
    done
done
fi

# continuum fits
#for flow in "0.20" ; do
#    for model in "--model pnorm --p 2" ; do #"--model max" "--model sum" "--model plaw --OmegaByT_IR 1.0 --OmegaByT_UV 6.28" ; do
#        for PhiUV in "--PhiUVtype LO  --min_scale 2piT --omega_prefactor 1" ; do # "--PhiUVtype NLO --min_scale opt --omega_prefactor opt" ; do  #
#            ./spf_reconstruct.py --add_suffix quenched_cont_f$flow --nsamples 500 --nproc 34 \
#                    --output_path /home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/EE/cont_rel_flow/rel_flow/ \
#                    --input_corr /home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/EE/cont_rel_flow/rel_flow/EE_relflow_${flow}.dat \
#                    --min_tauT 0.35 $model \
#                    --max_type smooth --T_in_GeV 0.472 --Nf 0 $PhiUV
#        done
#    done
#done
#fi

# finite a and tf fits
if [ 1 -eq 2 ] ; then
conftypes=(s080t20_b0703500 s144t36_b0754400)  # 24 30
Nts=(20 36)  # 24 30
flow="0.20"
for idx in "${!conftypes[@]}" ; do
    for model in "--model pnorm --p 2" "--model max" "--model sum" "--model plaw --OmegaByT_IR 1.0 --OmegaByT_UV 6.28" ; do
        for PhiUV in "--PhiUVtype LO  --min_scale 2piT --omega_prefactor 1" "--PhiUVtype NLO --min_scale opt --omega_prefactor opt" ; do  #
            ./spf_reconstruct.py --add_suffix quenched_${conftypes[idx]}_f${flow} --nsamples 500 --nproc 34 \
                        --output_path /home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/EE/${conftypes[idx]}/rel_flow/ \
                        --input_corr /home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/EE/${conftypes[idx]}/rel_flow/EE_relflow_${flow}.dat \
                        --min_tauT 0.35 $model \
                        --max_type smooth --T_in_GeV 0.472 --Nf 0 $PhiUV
        done
    done
done
fi

for model in "--model pnorm --p 2" "--model max" "--model sum" "--model plaw --OmegaByT_IR 1.0 --OmegaByT_UV 6.28" ; do
    for PhiUV in "--PhiUVtype LO  --min_scale 2piT --omega_prefactor 1" "--PhiUVtype NLO --min_scale opt --omega_prefactor opt" ; do  #
        ./spf_reconstruct.py --add_suffix quenched_BB --nsamples 500 --nproc 34 \
            --output_path /home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/BB/ \
            --input_corr /home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/BB/BB_flow_extr.txt \
            --min_tauT 0 $model \
            --max_type smooth --T_in_GeV 0.472 --Nf 0 $PhiUV
    done
done
