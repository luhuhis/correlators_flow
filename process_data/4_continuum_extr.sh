#!/bin/bash

#parameters
max_flow_idx=170 #135
min_flow_idx=0 #50
qcdtype=quenched_1.50Tc_zeuthenFlow
conftypes="s080t20_b0703500 s096t24_b0719200 s120t30_b0739400 s144t36_b0754400"
corrs="EE" # BB_clover EE BB
int="" #--int_only --int_Nt 64"
tex="" #"--use_tex"
max_FlowradiusBytauT=""  # "--max_FlowradiusBytauT 0.5"

for corr in $corrs; do
    rm -f /home/altenkort/work/correlators_flow/plots/$qcdtype/$corr/cont_extr*/*.pdf
    rm -f /home/altenkort/work/correlators_flow//data/merged/$qcdtype/$corr/cont_extr*/*.txt
done

for ((i=min_flow_idx; i < max_flow_idx; i+=10)) ; do
  for corr in $corrs; do
    /usr/local/bin/python3.7m -u _4_continuum_extr.py  --qcdtype $qcdtype --conftypes $conftypes --corr $corr --flow_index $((i)) $tex --custom_ylims 1.5 4   $max_FlowradiusBytauT $int &
    /usr/local/bin/python3.7m -u _4_continuum_extr.py  --qcdtype $qcdtype --conftypes $conftypes --corr $corr --flow_index $((i+1)) $tex --custom_ylims 1.5 4 $max_FlowradiusBytauT $int &
    /usr/local/bin/python3.7m -u _4_continuum_extr.py  --qcdtype $qcdtype --conftypes $conftypes --corr $corr --flow_index $((i+2)) $tex --custom_ylims 1.5 4 $max_FlowradiusBytauT $int &
    /usr/local/bin/python3.7m -u _4_continuum_extr.py  --qcdtype $qcdtype --conftypes $conftypes --corr $corr --flow_index $((i+3)) $tex --custom_ylims 1.5 4 $max_FlowradiusBytauT $int &
    /usr/local/bin/python3.7m -u _4_continuum_extr.py  --qcdtype $qcdtype --conftypes $conftypes --corr $corr --flow_index $((i+4)) $tex --custom_ylims 1.5 4 $max_FlowradiusBytauT $int &
    /usr/local/bin/python3.7m -u _4_continuum_extr.py  --qcdtype $qcdtype --conftypes $conftypes --corr $corr --flow_index $((i+5)) $tex --custom_ylims 1.5 4 $max_FlowradiusBytauT $int &
    /usr/local/bin/python3.7m -u _4_continuum_extr.py  --qcdtype $qcdtype --conftypes $conftypes --corr $corr --flow_index $((i+6)) $tex --custom_ylims 1.5 4 $max_FlowradiusBytauT $int &
    /usr/local/bin/python3.7m -u _4_continuum_extr.py  --qcdtype $qcdtype --conftypes $conftypes --corr $corr --flow_index $((i+7)) $tex --custom_ylims 1.5 4 $max_FlowradiusBytauT $int &
    /usr/local/bin/python3.7m -u _4_continuum_extr.py  --qcdtype $qcdtype --conftypes $conftypes --corr $corr --flow_index $((i+8)) $tex --custom_ylims 1.5 4 $max_FlowradiusBytauT $int &
    /usr/local/bin/python3.7m -u _4_continuum_extr.py  --qcdtype $qcdtype --conftypes $conftypes --corr $corr --flow_index $((i+9)) $tex --custom_ylims 1.5 4 $max_FlowradiusBytauT $int &
    wait
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
