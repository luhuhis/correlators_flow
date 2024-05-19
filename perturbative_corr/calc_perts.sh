#!/bin/bash

#SBATCH --partition=cpu_compute
#SBATCH --array=0-299
#SBATCH --qos=regular
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=128
#SBATCH --ntasks=1

export OMP_NUM_THREADS=1

flow_actions=("Zeuthen" "Wilson" "LW")
gauge_actions=("Wilson" "LW")
corrs=("EE" "BB")  #
Ntaus=(30 20 24 28 32 36 16 18 22 26 34 38 40 42 44 46 48 50 52 54 56 58 60 62 64)
param_list=()

for flow_action in "${flow_actions[@]}"  ; do
    for gauge_action in "${gauge_actions[@]}" ; do
        for corr in "${corrs[@]}" ; do
            for Ntau in "${Ntaus[@]}"; do
               param_list+=("--corr $corr --Nt $Ntau --flow_action $flow_action --gauge_action $gauge_action ")
            done
        done
    done
done

# TODO remove hard coded paths

srun -n1 -u /home/altenkort/work/correlators_flow/scripts/perturbative_corr/calc_pert_latt_corr_flow.py ${param_list[$((SLURM_ARRAY_TASK_ID))]} \
    --nproc 128 --printprogress \
    --flowtimes_file ~/work/correlators_flow/data/merged/pert_LO/flowtimes.dat \
    --outputpath ~/work/correlators_flow/data/merged/pert_LO/
