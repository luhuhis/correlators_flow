#!/bin/bash

#parameters
qcdtype=quenched_1.50Tc_zeuthenFlow
conftypes="s064t16_b0687361 s080t20_b0703500 s096t24_b0719200 s120t30_b0739400"
corrs="BB" # BB_clover"

for corr in $corrs; do
    rm -f /home/altenkort/work/correlators_flow/plots/$qcdtype/$corr/cont_extr*/*.pdf
    rm -f /home/altenkort/work/correlators_flow//data/merged/$qcdtype/$corr/$conftype/cont_extr*/*.txt
done


for corr in $corrs; do
    for i in {0..60} ; do
            python _4_continuum_extr.py --qcdtype $qcdtype --conftypes $conftypes --corr $corr --flow_index $i 
    done
done



#deprecated:
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
#sbatch <(cat <<< "$sbatchscript")
