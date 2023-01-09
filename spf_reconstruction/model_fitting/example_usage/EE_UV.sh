#!/bin/bash

for T in "--T_in_GeV 0.472 --Nf 0" "--T_in_GeV 0.251 --Nf 3" ; do
    for min_scale in "piT" "2piT" "opt" ; do
        for omega_prefactor in "1" "opt" ; do
            ../EE_UV_spf.py \
            $T \
            --min_scale $min_scale \
            --omega_prefactor $omega_prefactor \
            --outputpath "/work/home/altenkort/work/correlators_flow/data/merged/UV_spfs/" \

        done
    done
done
wait


