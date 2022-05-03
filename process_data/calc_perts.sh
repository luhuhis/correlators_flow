#!/bin/bash
export OMP_NUM_THREADS=1

# TODO FIX THE FLOW TIME FILES
corrs=("BB")  # "EE"
Ntaus=(20)  #  24 28 32
flow_actions=("Wilson" "Zeuthen" "LW")
gauge_actions=("Wilson" "LW")

for corr in "${corrs[@]}" ; do
    for Ntau in "${Ntaus[@]}"; do
        for flow_action in "${flow_actions[@]}"  ; do
            for gauge_action in "${gauge_actions[@]}" ; do
                ./calc_pert_latt_corr_flow.py \
                --printprogress --nproc 38 \
                --flowtimes_file /home/altenkort/work/slurm_scripts/run_gradientFlow/flowtimes/old/nt"${Ntau}".txt \
                --corr "$corr" --Nt "$Ntau" --flow_action "$flow_action" --gauge_action "$gauge_action" \
                --outputpath .  --Nspace 16
            done
        done
    done
done