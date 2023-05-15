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
# - bash 5.0 (Linux shell)
# - python 3.11
#   - packages: matplotlib numpy scipy rundec natsort colorama
# - working latex installation (e.g. texlive-full) for matplotlib tex rendering
# - gnuplot 5
#---

# ---
# INSTRUCTIONS
# These instructions can be copied into a bash 5.0 shell.
# If the files are not executable by default, you can make them executable via chmod +x <filename>.
# If the correct bash interpreter is not found automatically, you can execute the scripts with bash "manually" via
# "bash <path/to/script>"

# 1. SET UP THE ENVIRONMENT

# Unzip the analysis scripts. This will create a new folder "correlators_flow" in the current directory.
tar -xzf correlators_flow.tar.gz
# OR
# git clone https://github.com/luhuhis/correlators_flow.git # TODO fix release tag

# Add the scripts folder to the PYTHONPATH.
PYTHONPATH=$PYTHONPATH:$(pwd)/correlators_flow
export PYTHONPATH

# Unzip the AnalysisToolbox, which is a dependency of "correlators_flow". This will create a folder "AnalysisToobox"
tar -xzf AnalysisToolbox.tar.gz
# OR
# git clone https://github.com/LatticeQCD/AnalysisToolbox.git # TODO fix commit tag

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

# Specify how many processes can be created for parallelization.
export NPROC=1

# 2. RUN SCRIPTS

# 2.1 Double extrapolation
# TODO The data and plots are already contained in here, but can in principle be created by calling these script on the source data file.

# 2.1.1 Merge individual small text files into large binary (npy) files
# Merge individual EE correlator measurement text files (output from SIMULATeQCD) into a small number of larger numpy files (binary format).
# Meta data is saved to text files. This can take some time, mostly depending on file system speed (with conventional HDDs it may take hours).
tmppath="./correlators_flow/correlator_analysis/double_extrapolation/example_usage"
./$tmppath/1_merge_data.sh hisq_ms5_zeuthenFlow EE $BASEPATH_RAW_DATA $BASEPATH_WORK_DATA $BASEPATH_RAW_DATA/hisq_ms5_zeuthenFlow/reference_flowtimes

# Afterwards, the following files have been created inside
# $BASEPATH_WORK_DATA/hisq_ms5_zeuthenFlow/EE/<conftype>/

# flowtimes_<conftype>.dat              | flowtimes that were measured, corresponding to the data inside the .npy files
# n_datafiles_<conftype>.dat            | metadata (number of files per stream, MCMC trajectory number, etc.)
# EE_imag_<conftype>_merged.npy         | merged raw data, imaginary part of EE correlator
# EE_real_<conftype>_merged.npy         | merged raw data, real part of EE correlator
# polyakov_imag_<conftype>_merged.npy   | merged raw data, imaginary part of polyakovloop
# polyakov_real_<conftype>_merged.npy   | merged raw data, real part of polyakovloop

# where <conftype> may be, for example,
# s096t36_b0824900_m002022_m01011 (meaning Ns=96, Nt=36, beta=8.249, m_l=0.002022, m_s=0.01011). m_l and m_s are optional.

# 2.1.2 Resampling
# Load the merged data files, extract an equally spaced MCMC time series, then plot the time history of the Polyakov loop at a large flow time.
# Then bin configurations according to the integrated autocorrelation time, then perform bootstrap resampling of the uncorrelated blocks and save the samples to numpy files (binary format).
./$tmppath/2_reduce_data.sh hisq_ms5_zeuthenFlow EE $BASEPATH_WORK_DATA $BASEPATH_PLOT $NPROC

# Afterwards, the following files have been created:

# $BASEPATH_WORK_DATA/hisq_ms5_zeuthenFlow/EE/<conftype>/
# EE_<conftype>_samples.npy   | bootstrap samples of EE correlator
# EE_flow_cov_<conftype>.npy  | flow time correlation matrix (based on bootstrap samples)
# EE_<conftype>.dat           | mean EE correlator (useful for checking the data / plotting)
# EE_err_<conftype>.dat       | std_dev of EE correlator (useful for checking the data / plotting)

# $BASEPATH_PLOT/hisq_ms5_zeuthenFlow/EE/<conftype>
# polyakovloop_MCtime.pdf     | shows the MCMC time series of the polyakovloop at a large flow time

# interpolate in Euclidean time and in flow time
./$tmppath/3_spline_interpolate.sh hisq_ms5_zeuthenFlow EE $BASEPATH_WORK_DATA $BASEPATH_PLOT $NPROC

# continuum extrapolation
./$tmppath/4_continuum_extr.sh hisq_ms5_zeuthenFlow EE $BASEPATH_WORK_DATA $BASEPATH_PLOT $NPROC

# flow time extrapolation
./$tmppath/5_flowtime_extr.sh hisq_ms5_zeuthenFlow EE $BASEPATH_WORK_DATA $BASEPATH_PLOT $NPROC

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


# TODO fix nproc in spline and cont extr
