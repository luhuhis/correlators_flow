#!/bin/bash

#Job parameters
#SBATCH --job-name=spf_reconstruct
#SBATCH --output=./logs/%x_%A_%a.out
#SBATCH --error=./logs/%x_%A_%a.err

#SBATCH --time=10-00:00:00
#SBATCH --partition=volta
#SBATCH --qos=compute_cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --array=1-16

nsamples=998

params=()

for corr in EE BB ; do
  for constrain in "--constrain" "" ; do
    for mu in alpha beta; do
      for nmax in 4 5 ; do
         params+=("$constrain --nmax ${nmax} --mu ${mu} --corr ${corr}")
      done
    done
  done
done

#Print some information
echo -e "Start $(date +"%F %T") | $SLURM_JOB_ID $SLURM_JOB_NAME | $(hostname) | $(pwd) \n"

#Job steps
echo "${params[$((SLURM_ARRAY_TASK_ID-1))]}"
srun ./spf_reconstruct.py ${params[$((SLURM_ARRAY_TASK_ID-1))]} --seed 123456 --ignore_nans --guess_from_mean --qcdtype quenched_1.50Tc_zeuthenFlow --model 2 --PathPhiUV /work/data/htshu/ee_spf/PHIUV_a.dat --nsamples $nsamples --nproc 4 --PhiUVtype a --maxiter 2000


echo -e "End $(date +"%F %T") | $SLURM_JOB_ID $SLURM_JOB_NAME | $(hostname) | $(pwd) \n"