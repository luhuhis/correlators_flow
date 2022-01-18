#!/bin/bash

max_flow_idx=135
min_flow_idx=50
qcdtype=quenched_1.50Tc_zeuthenFlow
conftypes="s144t36_b0754400" #"s080t20_b0703500 s096t24_b0719200 s120t30_b0739400"
corrs="EE" # BB BB_clover
int_Nt=64
nsamples=1000

for conftype in $conftypes; do
    for corr in $corrs; do
        rm -f /home/altenkort/work/correlators_flow/plots/$qcdtype/$corr/$conftype/interpolations/*interpolation.pdf
        rm -f /home/altenkort/work/correlators_flow//data/merged/$qcdtype/$corr/$conftype/interpolations/*interpolation.txt
    done
done


for corr in $corrs; do
    for conftype in $conftypes ; do
        for ((i=min_flow_idx; i < max_flow_idx; i++)) ; do
            /usr/local/bin/python3.7m -u _3_spline_interpolate.py --qcdtype $qcdtype --conftype $conftype --corr $corr --flow_index $i --min_Nt 20 --int_Nt $int_Nt --nsamples $nsamples
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
