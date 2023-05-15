#!/bin/bash

# ---
# Data publication for "Heavy Quark Diffusion from 2+1 Flavor Lattice QCD with 320 MeV Pion Mass", 2023, arXiv 2302.08501
# Authors:
#   Luis Altenkort, altenkort@physik.uni-bielefeld.de, Universität Bielefeld#%
#   Olaf Kaczmarek, Universität Bielefeld
#   Rasmus Larsen, University of Stavanger
#   Swagato Mukherjee, Brookhaven National Laboratory
#   Peter Petreczky, Brookhaven National Laboratory
#   Hai-Tao Shu, Universität Regensburg
#   Simon Stendebach, TU Darmstadt
#   (HotQCD collaboration)
# ---

# ---
# Requirements
# - python 3.11
#   - packages: matplotlib numpy scipy rundec natsort
# - gnuplot 5
#---

# ---
# INSTRUCTIONS

# 1. SET UP THE ENVIRONMENT

# Unzip the analysis scripts. This will create a new folder "correlators_flow" in the current directory.
tar -xzf correlators_flow.tar.gz
# OR
git clone https://github.com/luhuhis/correlators_flow.git # TODO fix release tag

# Add the scripts folder to the PYTHONPATH.
PYTHONPATH=$PYTHONPATH:$(pwd)/correlators_flow
export PYTHONPATH

# Unzip the AnalysisToolbox, which is a dependency of "correlators_flow". This will create a folder "AnalysisToobox"
tar -xzf analysistoolbox.tar.gz
# OR
git clone https://github.com/LatticeQCD/AnalysisToolbox.git # TODO fix commit tag

# Add the AnalysisToolbox folder to the PYTHONPATH.
PYTHONPATH=$PYTHONPATH:$(pwd)/AnalysisToolbox
export PYTHONPATH

# Unzip raw data (gradientFlow output from SIMULATeQCD). This will create a new folder "hisq_ms5_zeuthenFlow".
tar -xzf hisq_ms5_zeuthenFlow.tar.gz

# Add the *parent* directory of the new folder (i.e., the directory where you called the "tar" command) to the environment. This should usually be $(pwd).
# For example, if you called tar in the folder /home/data, then a new folder /home/data/hisq_ms5_zeuthenFlow has been created. The parent of that directory is then /home/data.
BASEPATH_RAW_DATA=$(pwd)
export BASEPATH_RAW_DATA

# Unzip perturbative correlator data. This will create a new folder "pertLO". We want to add this new folder to the environment, e.g. /home/data/pertLO.
tar -xzf pertLO
G_PERT_LO_DIR=$(pwd)/pertLO
export G_PERT_LO_DIR

# Create two directories where you want to store the data of the intermediate steps (resampling, interpolation, extrapolations, ...) as well as final output data and plots.
export BASEPATH_WORK_DATA=#<YOUR-PATH-HERE>
export BASEPATH_PLOT=#<YOUR-PATH-HERE>


# 2. RUN SCRIPTS

# 2.1 Double extrapolation
# TODO The data and plots are already contained in here, but can in principle be created by calling these script on the source data file.

# Merge individual EE correlator measurement text files (output from SIMULATeQCD) into a small number of larger numpy files (binary format).
# Meta data is saved to text files. This can take multiple hours, mostly depending on file system speed.
tmppath="./correlators_flow/correlator_analysis/double_extrapolation/example_usage"
./$tmppath/1_merge_data.sh hisq_ms5_zeuthenFlow EE $BASEPATH_RAW_DATA $BASEPATH_WORK_DATA $BASEPATH_RAW_DATA/hisq_ms5_zeuthenFlow/reference_flowtimes

# TODO write which files have now been created, and what they contain roughly?

# Load the merged data files, extract an equally spaced MCMC time series, then plot the time history of the Polyakovloop at a large flow time.
# Then bin configurations according to the integrated autocorrelation time, then perform bootstrap resampling of the uncorrelated blocks and save the samples to numpy files (binary format).
./$tmppath/2_reduce_data.sh hisq_ms5_zeuthenFlow EE $BASEPATH_WORK_DATA $BASEPATH_PLOT

# interpolate in Euclidean time and in flow time
./$tmppath/3_spline_interpolate.sh hisq_ms5_zeuthenFlow EE $BASEPATH_WORK_DATA $BASEPATH_PLOT

# continuum extrapolation
./$tmppath/4_continuum_extr.sh hisq_ms5_zeuthenFlow EE $BASEPATH_WORK_DATA $BASEPATH_PLOT

# flow time extrapolation
./$tmppath/5_flowtime_extr.sh hisq_ms5_zeuthenFlow EE $BASEPATH_WORK_DATA $BASEPATH_PLOT

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
