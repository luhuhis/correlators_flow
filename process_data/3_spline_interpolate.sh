#!/bin/bash

max_flow_idx=100
qcdtype=quenched_1.50Tc_zeuthenFlow
conftypes="s064t16_b0687361 s080t20_b0703500 s096t24_b0719200 " #  s144t36_b0754400"
corrs="BB BB_clover"

for conftype in $conftypes; do  
    for corr in $corrs; do
        rm -f /home/altenkort/work/correlators_flow/plots/$qcdtype/$corr/$conftype/interpolations/*interpolation.pdf
        rm -f /home/altenkort/work/correlators_flow//data/merged/$qcdtype/$corr/$conftype/interpolations/*interpolation.txt
    done
done


for corr in $corrs; do
    for conftype in $conftypes ; do  
        for i in {0..60} ; do
            python _3_spline_interpolate.py --qcdtype $qcdtype --conftype $conftype --corr $corr --flow_index $i --int_Nt 30 
        done
    done
done


sbatchscript=$(cat <<EOF
#!/bin/bash
#SBATCH --job-name=spline_interpolate
#SBATCH --output=./%x.out
#SBATCH --error=./%x.err
#SBATCH --time=00:30:00
#SBATCH --partition=volta
#SBATCH --qos=compute_cpu
#SBATCH --ntasks=1
#SBATCH --array=0-600
#SBATCH --open-mode=append

export PATH=\$PATH:/usr/local/texlive/2020/bin/x86_64-linux
parameters=()
for corr in $corrs; do
    for conftype in $conftypes ; do  
        for i in {0..$max_flow_idx} ; do
            parameters+=("--qcdtype $qcdtype --conftype \$conftype --corr \$corr --flow_index \$i")
        done
    done
done
echo "\${parameters[\$SLURM_ARRAY_TASK_ID]}"
srun python -u _3_spline_interpolate.py  \${parameters[\$SLURM_ARRAY_TASK_ID]}

EOF
)
#sbatch <(cat <<< "$sbatchscript")
