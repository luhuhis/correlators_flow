#!/bin/bash

# Data publication for "Heavy Quark Diffusion from 2+1 Flavor Lattice QCD with 320 MeV Pion Mass", 2023
# Authors:
#   Luis Altenkort, altenkort@physik.uni-bielefeld.de, Universität Bielefeld#%
#   Olaf Kaczmarek, Universität Bielefeld
#   Rasmus Larsen, University of Stavanger
#   Swagato Mukherjee, Brookhaven National Laboratory
#   Peter Petreczky, Brookhaven National Laboratory
#   Hai-Tao Shu, Universität Regensburg
#   Simon Stendebach, TU Darmstadt
#   (HotQCD collaboration)

# TODO Requirements?
# pip3 install ... matplotplot numpy scipy rundec
# gnuplot

git clone git@github.com:luhuhis/correlators_flow.git

PYTHONPATH=$PYTHONPATH:$(pwd)/correlators_flow
export PYTHONPATH

# these paths should be absolute paths. (relative paths will be relative to the corresponding example usage folders.)
export BASEPATH_RAW_DATA=your-path
export BASEPATH_WORK_DATA=your-path
export BASEPATH_PLOT=your-path
export G_PERT_LO_DIR # TODO

# EXTRACT DATA
# tar -xzf data.tar.gz $INPUT_BASEPATH

# The data and plots are already contained in here, but can in principle be created by calling these script on the source data file.

# merge individual measurement text files (output from SIMULATeQCD) into a small number of larger numpy files. meta data is saved to text files. this can take multiple hours, mostly depending on file system speed.
scriptpath1="./correlators_flow/correlator_analysis/double_extrapolation/example_usage"
./$scriptpath1/1_merge_data.sh hisq_ms5_zeuthenFlow EE $BASEPATH_RAW_DATA $BASEPATH_WORK_DATA

# TODO write which files have now been created, and what they contain roughly?

# load the merged data files, "clean" the time series to make it equally spaced, then plot the time history of the Polyakovloop, then bin configurations according to the integrated autocorrelation time, then perform bootstrap resampling and save the samples to numpy files.
./$scriptpath1/2_reduce_data.sh hisq_ms5_zeuthenFlow EE $BASEPATH_WORK_DATA $BASEPATH_PLOT

# interpolate in Euclidean time and in flow time
./$scriptpath1/3_spline_interpolate.sh hisq_ms5_zeuthenFlow EE $BASEPATH_WORK_DATA $BASEPATH_PLOT

# continuum extrapolation
./$scriptpath1/4_continuum_extr.sh hisq_ms5_zeuthenFlow EE $BASEPATH_WORK_DATA $BASEPATH_PLOT

# flow time extrapolation
./$scriptpath1/5_flowtime_extr.sh hisq_ms5_zeuthenFlow EE $BASEPATH_WORK_DATA $BASEPATH_PLOT

# create relative flow time plots
./correlators_flow/correlator_analysis/relative_flow/example_usage/Nf3TemperatureComparison_Paper.sh $BASEPATH_WORK_DATA $BASEPATH_PLOT

# plot 2piTD figure
(
    cd 2piTD_plot || exit
    gnuplot 2piTDs.plt
)

# spf fits. The data line is optional
./correlators_flow/spf_reconstruction/model_fitting/example_usage/spf_reconstruct.sh $BASEPATH_WORK_DATA NO

# plot spf fit results
./correlators_flow/spf_reconstruction/plot_fits/example_usage/plot_fits_hisq.sh $BASEPATH_WORK_DATA $BASEPATH_PLOT

# k factor
./correlators_flow/spf_reconstruction/plot_fits/example_usage/plot_kfactors.sh $BASEPATH_WORK_DATA $BASEPATH_PLOT

# plot final kappa plot
./correlators_flow/spf_reconstruction/plot_fits/example_usage/plot_final_kappas.sh $BASEPATH_WORK_DATA $BASEPATH_PLOT

# TODO: make files executable locally and push them
# TODO include pertLO files in data pub
# TODO why are there two interpolation plots?
