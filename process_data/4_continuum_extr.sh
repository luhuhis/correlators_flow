#! /bin/bash
#SBATCH --job-name=cont_extr
#SBATCH --output=./%x_%j.out
#SBATCH --error=./%x_%j.err
#SBATCH --time=10:00:00
#SBATCH --partition=compute_cpu
#SBATCH --ntasks-per-socket=4
#SBATCH --cores-per-socket=4
#SBATCH --sockets-per-node=2
#SBATCH --nodes=10
echo -e "Start `date +"%F %T"` | $SLURM_JOB_ID $SLURM_JOB_NAME | `hostname` | `pwd` \n" #log start time etc.
for i in {0..220} ; do
    srun --ntasks=1 --nodes=1 --exclusive python3 4_continuum_extr.py $i &    
done
wait
echo -e "Start `date +"%F %T"` | $SLURM_JOB_ID $SLURM_JOB_NAME | `hostname` | `pwd` \n" #log end time etc.
