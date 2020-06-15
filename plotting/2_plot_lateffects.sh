#! /bin/bash
#SBATCH --job-name=plot_interpolations
#SBATCH --output=./%x.out
#SBATCH --error=./%x.err
#SBATCH --time=01:00:00
#SBATCH --partition=compute_cpu
#SBATCH --ntasks=40
#SBATCH --cores-per-socket=1

echo -e "Start `date +"%F %T"` | $SLURM_JOB_ID $SLURM_JOB_NAME | `hostname` | `pwd` \n" #log start time etc.
rm -f /work/home/altenkort/work/EE_correlator/plots/quenched/lattice_effects/*
for i in {40..150} ; do
    python 2_plot_lateffects.py $i &    
done
wait
echo -e "Start `date +"%F %T"` | $SLURM_JOB_ID $SLURM_JOB_NAME | `hostname` | `pwd` \n" #log end time etc.
 
