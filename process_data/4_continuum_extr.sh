rm -f /work/home/altenkort/work/EE_correlator/data_merged/quenched/cont_extr/*
rm -f /work/home/altenkort/work/EE_correlator/data_merged/quenched/cont_extr_quality/*
rm -f /work/home/altenkort/work/EE_correlator/plots/quenched/cont_extr/*
rm -f /work/home/altenkort/work/EE_correlator/plots/quenched/cont_extr_quality/*

sbatchscript=$(cat <<EOF
#! /bin/bash
#SBATCH --job-name=cont_extr
#SBATCH --output=./%x.out
#SBATCH --error=./%x.err
#SBATCH --time=01:00:00
#SBATCH --partition=devel_cpu
#SBATCH --ntasks=1
#SBATCH --array=0-83
parameters=()
for i in {50..133} ; do
    parameters+=("\$i")
done

export PATH=\$PATH:/usr/local/texlive/2020/bin/x86_64-linux
srun python -u _4_continuum_extr.py \${parameters[\$SLURM_ARRAY_TASK_ID]}

EOF
)
sbatch <(cat <<< "$sbatchscript")
