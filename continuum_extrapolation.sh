#!/bin/bash
for flowradius in 0.0 0.005 0.01 0.015 0.02 0.025 0.03 0.035 0.04 0.045 0.05 0.055 0.06 0.065 0.07 0.075 ; do #0.08 0.085 0.09 0.095 0.1 0.125 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 ; do
    ../../code/analysistoolbox/bin/extrapolate.py --method gauss_btstr --nknots 2 3 4 5 --order 3 --randomization-factor 0.9 --xmax 0.5 --nsamples 1000 --constraints 0.5 1 0 --plot-xmax 0.5 --plot-xmin 0 --outname=EE_$flowradius.dat --folder=../data_merged/quenched/continuum_limit/ ../data_merged/quenched/s*t*_b*/continuum_limit/EE_${flowradius}_Nt*.dat
done

#--log-level PROGRESS
