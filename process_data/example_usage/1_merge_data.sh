#!/bin/bash


#!/bin/bash

qcdtype=$1
corr=$2
if [ -z "$qcdtype" ] ; then
    echo "Usage: $0 qcdtype corr"
    echo "choices for qcdtype: quenched_1.50Tc_zeuthenFlow hisq_ms5_zeuthenFlow"
    echo "choices for corr: EE BB EE_clover BB_clover"
    exit
fi


if [ "$qcdtype" == quenched_1.50Tc_zeuthenFlow ] ; then
    arr_conftypes=("s064t16_b0687361 s080t20_b0703500 s096t24_b0719200 s120t30_b0739400 s144t36_b0754400")
    acc_sts="acc0.000010_sts0.000010"
elif [ "$qcdtype" == hisq_ms5_zeuthenFlow ] ; then
    arr_conftypes=(
    "s064t20_b0803600" "s064t22_b0814700" "s096t24_b0824900_m002022_m01011"
    "s064t20_b0785700" "s064t24_b0806800" "s096t28_b0824900_m002022_m01011"
    "s064t20_b0770400" "s064t24_b0791300" "s096t32_b0824900_m002022_m01011"
    "s064t20_b0757000" "s064t24_b0777700" "s096t36_b0824900_m002022_m01011"
)
    arr_discards=(
        "0 100 100 100" "0 0 0 0 0 0" "0 0 0"
        "0 100 100 100" "0 0 0 0 0 0" "400 400 400 400"
        "0 0 0 0 0 0" "0 200 200 200" "400 400"
        "0 0 0 0 0 0" "0 100 100 100" "400 400 400 400"
    )
    temps=(
        "296"
        "251"
        "220"
        "196"
        )
    basepath="/home/altenkort/work/correlators_flow/data/merged/hisq_ms5_zeuthenFlow/EE/"
    add_args="--basepath /work/data/altenkort/gradientFlow"
    flowradii_refs=()
    for temp in "${temps[@]}"; do
        flowradii_refs+=("--reference_flowradii $basepath/T${temp}/flowradii_T${temp}.dat")
    done
fi


for idx in "${!arr_conftypes[@]}" ; do
    conftype="${arr_conftypes[idx]}"

    if [ "$qcdtype" == quenched_1.50Tc_zeuthenFlow ] ; then
        discard="0 0 0 0 0 0 0 0 0 0"
    elif [ "$qcdtype" == hisq_ms5_zeuthenFlow ] ; then
        discard=${arr_discards[idx]}
        if [ "$conftype" == "s096t24_b0824900_m002022_m01011" ] ; then
            acc_sts=""
        else
            acc_sts="--acc_sts sts0.150000"
        fi
        if [ "$conftype" == "s096t36_b0824900_m002022_m01011" ] ; then
            excess="--excess_workaround"
        else
            excess=""
        fi
    fi


    ../_1_merge_data.py ${flowradii_refs[idx/3]} --output_basepath ../../../data/merged/ --n_discard $discard $excess $add_args --qcdtype $qcdtype --corr $corr $acc_sts --conftype $conftype
done
wait


#--n_discard 400
#hisq_b8249_zeuthenFlow
#
#
#




#TODO update this and incorporate intersecting flowtimes

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

#wait
