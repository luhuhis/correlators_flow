#!/bin/bash

# Data publication instructions for TODO


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

# tar -xzf data.tar.gz $INPUT_BASEPATH

# TODO: make files executable locally and push them

scriptpath1="./correlators_flow/correlator_analysis/double_extrapolation/example_usage"


# merge individual measurement text files (output from SIMULATeQCD) into a small number of larger numpy files. meta data is saved to text files. this can take multiple hours, mostly depending on file system speed.
./$scriptpath1/1_merge_data.sh hisq_ms5_zeuthenFlow EE $BASEPATH_RAW_DATA $BASEPATH_WORK_DATA

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

# spf fits



# TODO probably some other files are missing which are hardcoded in lib process data. test this with the new environment on the other laptop.