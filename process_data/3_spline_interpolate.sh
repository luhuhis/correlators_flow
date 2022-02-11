#!/bin/bash

max_flow_idx=170
min_flow_idx=0
qcdtype=quenched_1.50Tc_zeuthenFlow
conftypes="s080t20_b0703500 s096t24_b0719200 s120t30_b0739400 s144t36_b0754400" #""
corrs="EE" # BB BB_clover EE
int_Nt=36
nsamples=1000
max_FlowradiusBytauT=""  # "--max_FlowradiusBytauT 0.33"


for conftype in $conftypes; do
    for corr in $corrs; do
        rm -f /home/altenkort/work/correlators_flow/plots/$qcdtype/$corr/$conftype/interpolations/*interpolation.pdf
        rm -f /home/altenkort/work/correlators_flow//data/merged/$qcdtype/$corr/$conftype/interpolations/*interpolation.txt
    done
done


for corr in $corrs; do
    for conftype in $conftypes ; do
        for ((i=min_flow_idx; i < max_flow_idx; i+=10)) ; do
            /usr/local/bin/python3.7m -u _3_spline_interpolate.py $max_FlowradiusBytauT --qcdtype $qcdtype --conftype $conftype --corr $corr --flow_index $((i))   --int_Nt $int_Nt --nsamples $nsamples &
            /usr/local/bin/python3.7m -u _3_spline_interpolate.py $max_FlowradiusBytauT --qcdtype $qcdtype --conftype $conftype --corr $corr --flow_index $((i+1)) --int_Nt $int_Nt --nsamples $nsamples &
            /usr/local/bin/python3.7m -u _3_spline_interpolate.py $max_FlowradiusBytauT --qcdtype $qcdtype --conftype $conftype --corr $corr --flow_index $((i+2)) --int_Nt $int_Nt --nsamples $nsamples &
            /usr/local/bin/python3.7m -u _3_spline_interpolate.py $max_FlowradiusBytauT --qcdtype $qcdtype --conftype $conftype --corr $corr --flow_index $((i+3)) --int_Nt $int_Nt --nsamples $nsamples &
            /usr/local/bin/python3.7m -u _3_spline_interpolate.py $max_FlowradiusBytauT --qcdtype $qcdtype --conftype $conftype --corr $corr --flow_index $((i+4)) --int_Nt $int_Nt --nsamples $nsamples &
            /usr/local/bin/python3.7m -u _3_spline_interpolate.py $max_FlowradiusBytauT --qcdtype $qcdtype --conftype $conftype --corr $corr --flow_index $((i+5)) --int_Nt $int_Nt --nsamples $nsamples &
            /usr/local/bin/python3.7m -u _3_spline_interpolate.py $max_FlowradiusBytauT --qcdtype $qcdtype --conftype $conftype --corr $corr --flow_index $((i+6)) --int_Nt $int_Nt --nsamples $nsamples &
            /usr/local/bin/python3.7m -u _3_spline_interpolate.py $max_FlowradiusBytauT --qcdtype $qcdtype --conftype $conftype --corr $corr --flow_index $((i+7)) --int_Nt $int_Nt --nsamples $nsamples &
            /usr/local/bin/python3.7m -u _3_spline_interpolate.py $max_FlowradiusBytauT --qcdtype $qcdtype --conftype $conftype --corr $corr --flow_index $((i+8)) --int_Nt $int_Nt --nsamples $nsamples &
            /usr/local/bin/python3.7m -u _3_spline_interpolate.py $max_FlowradiusBytauT --qcdtype $qcdtype --conftype $conftype --corr $corr --flow_index $((i+9)) --int_Nt $int_Nt --nsamples $nsamples &
            wait
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
