#!/bin/bash
for samples in 1000 ; do
  for obs in kappa chisqdof ; do
    for corr in EE BB ; do
      for constrain in "--constrain" "" ; do
        for nmax in 4 5 ; do
          for mu in alpha beta; do
             ./plot_dist.py --model 2 --qcdtype quenched_1.50Tc_zeuthenFlow --corr $corr $constrain --nmax ${nmax} --mu ${mu} --nsamples $samples --PhiUVtype a --obs $obs --PathOutputFolder ./plots/
          done
        done
      done
    done
    wait
  done
done
wait