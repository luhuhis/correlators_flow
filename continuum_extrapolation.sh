#!/bin/bash
#SBATCH --job-name=EE_continuum_extrapolation
#SBATCH --output=../data_merged/quenched/continuum_limit/logs/%x_%A_%a.out
#SBATCH --error=../data_merged/quenched/continuum_limit/logs/%x_%A_%a.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=altenkort@physik.uni-bielefeld.de

#SBATCH --partition=volta
#SBATCH --array=0-15
#SBATCH --ntasks=1
echo "${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} ${SLURM_JOB_NAME}"; date; hostname; pwd;
mkdir -p ../data_merged/quenched/continuum_limit/logs/
flowradii=(0.0 0.005 0.01 0.015 0.02 0.025 0.03 0.035 0.04 0.045 0.05 0.055 0.06 0.065 0.07 0.075) #0.08 0.085 0.09 0.095 0.1 0.125 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 ; do --xmin 0.045 --xmax 0.47
srun python -u ~/code/analysistoolbox/bin/extrapolate.py --data-input direct --method gauss_btstr --nknots 1 --order 2 --randomization-factor 0 --nsamples 1000 --no-tex --outname=EE_${flowradii[$SLURM_ARRAY_TASK_ID]} --folder=../data_merged/quenched/continuum_limit/ ../data_merged/quenched/s*t*_b*/continuum_limit/EE_${flowradii[$SLURM_ARRAY_TASK_ID]}_Nt*.dat --xmin=`bc <<< "scale=5; 1*${flowradii[$SLURM_ARRAY_TASK_ID]}"` --xmax=`bc <<< "scale=5; 0.5-(1*(0.075-${flowradii[$SLURM_ARRAY_TASK_ID]}))"`
date

# --constraints 0 1 0 0.5 1 0
