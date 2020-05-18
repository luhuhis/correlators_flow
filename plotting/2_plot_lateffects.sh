#! /bin/bash
#SBATCH --job-name=plot_interpolations
#SBATCH --output=./%x_%j.out
#SBATCH --error=./%x_%j.err
#SBATCH --time=01:00:00
#SBATCH --partition=compute_cpu
#SBATCH --ntasks=50
#SBATCH --nodes=5
#SBATCH --cores-per-socket=5
echo -e "Start `date +"%F %T"` | $SLURM_JOB_ID $SLURM_JOB_NAME | `hostname` | `pwd` \n" #log start time etc.
rm -f /work/home/altenkort/work/EE_correlator/plots/quenched/lattice_effects/*
for i in {0..220} ; do
    srun --ntasks=1 --nodes=1 --exclusive python3 2_plot_lateffects.py $i &    
done
wait
echo -e "Start `date +"%F %T"` | $SLURM_JOB_ID $SLURM_JOB_NAME | `hostname` | `pwd` \n" #log end time etc.
 
