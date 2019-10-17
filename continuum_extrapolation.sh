#!/bin/bash
# rm -rf ../data_merged/quenched/continuum_limit/
mkdir -p ../data_merged/quenched/continuum_limit/
mkdir -p ../data_merged/quenched/continuum_limit/logs/
extrapolate="srun -n 1 -N 1 --exclusive python -u /home/altenkort/code/analysistoolbox/bin/extrapolate.py"
path1=../data_merged/quenched
args="--nknots 1 2 3 --order 2 --constraints 0 2 0 0.5 1 0 --randomization-factor 0 --nsamples 100 --base-point=0.31 --xmax=0.5 --data-input sample --method from_sample  --no-tex --plot-xmin 0.134 --plot-xmax 0.5 --folder=$path1/continuum_limit/ --log-level INFO"
#--nknots 3
chi_dof_file=../data_merged/quenched/continuum_limit/logs/chi_dof.txt

sbatchscript=$(cat <<EOF
#!/bin/bash
#SBATCH --job-name=EE_continuum_extrapolation
#SBATCH --open-mode=append
#SBATCH --output=../data_merged/quenched/continuum_limit/logs/%x_.out
#SBATCH --error=../data_merged/quenched/continuum_limit/logs/%x_.err
#SBATCH --mail-type=NONE
#SBATCH --mail-user=luis.a@hotmail.de


#SBATCH --partition=volta
#SBATCH --nodes=2
#SBATCH --tasks-per-node=12

echo "\${SLURM_JOB_ID} \${SLURM_JOB_NAME} \`date\` \`hostname\` \`pwd\`"

flowradii=(0.0500 0.0550 0.0600 0.0650 0.0700 0.0750 0.0800 0.0850 0.0900 0.0950 0.1000)

for i in "\${flowradii[@]}" ; do  
$extrapolate --xmin=\`bc <<< "scale=5; \${i}/sqrt(8*0.014)+0.02"\` $args --outname=EE_\${i} $path1/s*t*_b*/single_flow/EE_\${i}_Nt*_btstrp_samples.dat >> $path1/continuum_limit/logs/\${SLURM_JOB_NAME}_\${i}.out &
done
wait

#save chi_dofs
rm -f $chi_dof_file
echo "#flowradius chi2_dof \`date\`" >> $chi_dof_file
for i in "\${flowradii[@]}" ; do   
    echo "\$i \`tail -n 3 ../data_merged/quenched/continuum_limit/logs/EE_continuum_extrapolation_\${i}.out | grep 'Chi^2/d.o.f. of first full fit:' | grep -Eo '[+-]?[0-9]+([.][0-9]+)?' | tail -n 1\`" >> $chi_dof_file
done

srun -n 1 -N 1 --exclusive python -u ./misc/av_chi_dof.py

echo "done \${SLURM_JOB_ID} \`date\`"
EOF
)
sbatch <(cat <<< "$sbatchscript")
