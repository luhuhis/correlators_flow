#!/bin/bash
rm -f ../../data_merged/quenched/final_corr_samples/*
rm -f ../../plots/quenched/final_corr_samples/*

sbatchscript=$(cat <<EOF
#!/bin/bash
#SBATCH --job-name=EE_btstrp
#SBATCH --output=./%x.out
#SBATCH --error=./%x.err
#SBATCH --time=00:05:00
#SBATCH --partition=devel_cpu
#SBATCH --ntasks=1
#SBATCH --cores-per-socket=1

#SBATCH --array=0-9999

srun python complete_analysis_on_a_sample.py \$SLURM_ARRAY_TASK_ID

EOF
)
sbatch <(cat <<< "$sbatchscript")
