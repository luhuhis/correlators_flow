#!/bin/bash

qcdtype=$1
if [ -z "$qcdtype" ] ; then
    echo "Usage: $0 qcdtype"
    echo "choices for qcdtype: quenched_1.50Tc_zeuthenFlow hisq_ms5_zeuthenFlow"
    exit
fi


if [ "$qcdtype" == quenched_1.50Tc_zeuthenFlow ] ; then
    conftypes=(s080t20_b0703500 s096t24_b0719200 s120t30_b0739400 s144t36_b0754400)
elif [ "$qcdtype" == hisq_ms5_zeuthenFlow ] ; then
    conftypes=(s064t16_b07188_m00113_m0306 s064t16_b07010_m00132_m0357 s064t16_b07095_m00124_m0334 s064t16_b07130_m00119_m0322 s064t16_b07054_m00129_m0348 s064t16_b07156_m00116_m0314)
fi


for conftype in "${conftypes[@]}" ; do
    ../_2_reduce_data.py --qcdtype "$qcdtype" --corr EE --conftype "$conftype" &
done
wait