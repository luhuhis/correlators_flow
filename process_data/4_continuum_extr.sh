#! /bin/bash
#SBATCH --job-name=cont_extr
#SBATCH --output=./%x.out
#SBATCH --error=./%x.err
#SBATCH --time=24:00:00
#SBATCH --partition=compute_cpu
#SBATCH --ntasks=83
#SBATCH --cores-per-socket=1

echo -e "Start `date +"%F %T"` | $SLURM_JOB_ID $SLURM_JOB_NAME | `hostname` | `pwd` \n" #log start time etc.
rm -f /work/home/altenkort/work/EE_correlator/data_merged/quenched/cont_extr/*
rm -f /work/home/altenkort/work/EE_correlator/data_merged/quenched/cont_extr_quality/*
rm -f /work/home/altenkort/work/EE_correlator/plots/quenched/cont_extr/*
rm -f /work/home/altenkort/work/EE_correlator/plots/quenched/cont_extr_quality/*
#srun --ntasks=1 --nodes=1 --exclusive
for i in {50..150} ; do
    srun --ntasks=1 --nodes=1 --exclusive python 4_continuum_extr.py $i &
done
wait
echo -e "End `date +"%F %T"` | $SLURM_JOB_ID $SLURM_JOB_NAME | `hostname` | `pwd` \n" #log end time etc.
