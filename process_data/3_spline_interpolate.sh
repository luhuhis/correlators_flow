#!/bin/bash
conftypes="s080t20_b0703500 s096t24_b0719200 s120t30_b0739400 s144t36_b0754400"
for conftype in $conftypes; do  
    rm -f /work/home/altenkort/work/EE_correlator/plots/quenched/$conftype/interpolations/*
    rm -f /work/home/altenkort/work/EE_correlator/data_merged/quenched/$conftype/interpolations/*
done

sbatchscript=$(cat <<EOF
#!/bin/bash
#SBATCH --job-name=spline_interpolate
#SBATCH --output=./%x.out
#SBATCH --error=./%x.err
#SBATCH --time=01:00:00
#SBATCH --partition=devel_cpu
#SBATCH --ntasks=1
#SBATCH --array=0-536
#SBATCH --open-mode=append

export PATH=\$PATH:/usr/local/texlive/2020/bin/x86_64-linux
parameters=()
for conftype in $conftypes ; do  
    for i in {0..133} ; do
        parameters+=("\$conftype \$i")
    done
done
echo "\${parameters[\$SLURM_ARRAY_TASK_ID]}"
srun python -u _3_spline_interpolate.py quenched \${parameters[\$SLURM_ARRAY_TASK_ID]}

EOF
)
sbatch <(cat <<< "$sbatchscript")
