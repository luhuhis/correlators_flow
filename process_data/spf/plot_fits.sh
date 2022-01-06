#!/bin/bash
for obs in "kappa --xlims 1.5 12.5" "chisqdof --xlims 0 9" ; do
for tol in "1e-04" "1e-05" "1e-06" ; do
  for start in "" "--start_from_mean_fit" ; do
    ./plot_fits.py --qcdtype quenched_1.50Tc_zeuthenFlow --corr EE --nsamples 500 --PathOutputFolder . --tol $tol $start --obs $obs --usetex
  done
done
done
