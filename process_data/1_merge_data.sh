#!/bin/bash
#conftypes_quenched=${1:-s064t16_b0687361 s080t20_b0703500 s096t24_b0719200} #s120t30_b0739400 s096t16_b0719200 s096t28_b0719200 s096t32_b0719200
conftypes_hisq=(s096t36_b0824900_m002022_m01011 s096t28_b0824900_m002022_m01011 s064t64_b0824900_m002022_m01011 s096t20_b0824900_m002022_m01011 s096t24_b0824900_m002022_m01011 s096t32_b0824900_m002022_m01011 s096t40_b0824900_m002022_m01011 s096t44_b0824900_m002022_m01011 s096t48_b0824900_m002022_m01011 s096t56_b0824900_m002022_m01011)

for conftype in "${conftypes_hisq[@]}" ; do
    ./_1_merge_data.py --n_discard 400 --qcdtype hisq_b8249_zeuthenFlow --corr EE --acc_sts sts0.150000 --conftype $conftype --basepath /work/data/altenkort/gradientFlow
done



# conftypes_hisq=${2:-s064t16_b07188_m00113_m0306 s064t16_b07010_m00132_m0357 s064t16_b07095_m00124_m0334 s064t16_b07130_m00119_m0322 s064t16_b07054_m00129_m0348 s064t16_b07156_m00116_m0314}
# for conftype in $conftypes_hisq ; do
#     #rm -rf ../data_merged/hisq/$conftype
#     #rm -rf ../data_merged/hisq/$conftype/single_flow
#     mkdir -p ../data_merged/hisq/$conftype
#     beta=${conftype#*_b}; beta=${beta%%_*}; beta=`bc <<< "scale=5;$beta/1000"`
#     ns=${conftype#s}; ns=${ns%%t*}
#     nt=${conftype#*t}; nt=${nt%%_b*}
#     python3 merge_raw_data.py hisq $conftype $beta $ns $nt &
# done

wait
