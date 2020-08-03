#!/bin/bash
rm -f /work/home/altenkort/work/EE_correlator/plots/quenched/lattice_effects/*
rm -f /home/altenkort/.cache/matplotlib/tex.cache/*matplotlib-lock

sbatchscript=$(cat <<EOF
#! /bin/bash
#SBATCH --job-name=plot_interpolations
#SBATCH --output=./%x.out
#SBATCH --error=./%x.err
#SBATCH --time=24:00:00
#SBATCH --partition=devel_cpu
#SBATCH --ntasks=1
##SBATCH --array=0-133

export PATH=\$PATH:/usr/local/texlive/2020/bin/x86_64-linux
parameters=()
for i in {0..133} ; do
    srun --exclusive python -u 2_plot_lateffects.py \$i
done
#parameters+=("\$i")
#\${parameters[\$SLURM_ARRAY_TASK_ID]}
EOF
)
sbatch <(cat <<< "$sbatchscript")
