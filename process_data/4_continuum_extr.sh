#!/bin/bash

#parameters
max_flow_idx=135
min_flow_idx=50
qcdtype=quenched_1.50Tc_zeuthenFlow
conftypes="s080t20_b0703500 s096t24_b0719200 s120t30_b0739400 s144t36_b0754400"
corrs="EE" #BB BB_clover

for corr in $corrs; do
    rm -f /home/altenkort/work/correlators_flow/plots/$qcdtype/$corr/cont_extr*/*.pdf
    rm -f /home/altenkort/work/correlators_flow//data/merged/$qcdtype/$corr/cont_extr*/*.txt
done

for ((i=min_flow_idx; i < max_flow_idx; i++)) ; do
  for corr in $corrs; do
    /usr/local/bin/python3.7m -u _4_continuum_extr.py --qcdtype $qcdtype --conftypes $conftypes --corr $corr --flow_index $i --use_tex --custom_ylims 1.5 4 --int_only --int_Nt 64
  done
done
#wait


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
