#!/bin/bash

# Data publication instructions for TODO


# TODO Requirements?
# pip3 install ... matplotplot numpy scipy rundec


git clone git@github.com:luhuhis/correlators_flow.git

export PYTHONPATH=$PYTHONPATH:$(pwd)/correlators_flow

export BASEPATH_RAW_DATA=<your-path>
export BASEPATH_WORK_DATA=<your-path>
export BASEPATH_PLOT=<your-path>

# tar -xzf data.tar.gz $INPUT_BASEPATH

# TODO: make files executable locally and push them

# merge individual measurement text files (output from SIMULATeQCD) into a small number of larger numpy files. meta data is saved to text files. this can take multiple hours, mostly depending on file system speed.
./correlators_flow/correlator_analysis/double_extrapolation/example_usage/1_merge_data.sh hisq_ms5_zeuthenFlow EE $BASEPATH_RAW_DATA $BASEPATH_WORK_DATA

# load the merged data files, "clean" the time series to make it equally spaced, then plot the time history of the Polyakovloop, then bin configurations according to the integrated autocorrelation time, then perform bootstrap resampling and save the samples to numpy files.
./correlators_flow/correlator_analysis/double_extrapolation/example_usage/2_reduce_data.sh hisq_ms5_zeuthenFlow EE $BASEPATH_WORK_DATA $BASEPATH_PLOT

# interpolate in Euclidean time and in flow time
./correlators_flow/correlator_analysis/double_extrapolation/example_usage/3_spline_interpolate.sh hisq_ms5_zeuthenFlow EE $BASEPATH_WORK_DATA $BASEPATH_PLOT


#
#
#4_continuum double_extrapolation
#5_flowtime double_extrapolation


# spf reconstruction
