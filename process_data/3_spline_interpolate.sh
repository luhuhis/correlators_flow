#! /bin/bash
#SBATCH --job-name=spline_interpolate
#SBATCH --output=./%x.out
#SBATCH --error=./%x.err
#SBATCH --time=23:00:00
#SBATCH --partition=compute_cpu
#SBATCH --cores-per-socket=1
#SBATCH --ntasks=100

echo -e "Start `date +"%F %T"` | $SLURM_JOB_ID $SLURM_JOB_NAME | `hostname` | `pwd` \n" #log start time etc.
for conftype in s064t16_b0687361 s080t20_b0703500 s120t30_b0739400 s096t24_b0719200 ; do #s144t36_b0754400 #
      rm -f /work/home/altenkort/work/EE_correlator/plots/quenched/$conftype/interpolations/*
      rm -f /work/home/altenkort/work/EE_correlator/data_merged/quenched/$conftype/interpolations/*
    for i in {42..150} ; do
        srun --ntasks=1 --nodes=1 --exclusive python3.6 3_spline_interpolate.py quenched $conftype $i &    
    done
done
wait
echo -e "Start `date +"%F %T"` | $SLURM_JOB_ID $SLURM_JOB_NAME | `hostname` | `pwd` \n" #log end time etc.
