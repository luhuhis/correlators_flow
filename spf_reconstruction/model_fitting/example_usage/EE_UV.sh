#!/bin/bash

outputpath=$1

if [ -z "$outputpath" ] ; then
    echo "Usage: $0 outputpath"
    exit
fi

outputpath="${outputpath}/UV_spfs/"
mkdir -p "${outputpath}"

min_scale="eff"
for T in "--T_in_GeV 0.472 --Nf 0" "--T_in_GeV 0.251 --Nf 3" ; do
    for omega_prefactor in "1" "opt" ; do
        (
            cd "$(dirname $0)" || exit
            ../compute_UV_spf.py \
            ${T} \
            --min_scale $min_scale \
            --omega_prefactor $omega_prefactor \
            --outputpath "${outputpath}" \
            --corr EE \
            --order LO  # dummy value
        )
    done
done


