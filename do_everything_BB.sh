#!/bin/bash

# ---
# Data publication for "Lattice B-field correlators for heavy quarks", 2023, arXiv # TODO
# Authors:
#   Luis Altenkort, altenkort@physik.uni-bielefeld.de, Universität Bielefeld
#   David de la Cruz, TU Darmstadt
#   Olaf Kaczmarek, Universität Bielefeld
#   Guy D. Moore, TU Darmstadt
#   Hai-Tao Shu, Brookhaven National Laboratory
# ---

# ---
# Requirements
# - bash 5.0 (Linux shell)
# - python 3.11
#   - packages: matplotlib numpy scipy rundec natsort colorama
# - working latex installation with packages amsmath and mathtools for matplotlib tex rendering
# - gnuplot 5
#---

# ---
# INSTRUCTIONS
# These instructions can be copied into a bash 5.0 shell.

# 1. SET UP THE ENVIRONMENT

# TODO add a summary / one-liner of setting the env

export NPROC=100 && export PYTHONPATH=$(pwd)/correlators_flow:$(pwd)/AnalysisToolbox:$PYTHONPATH && export BASEPATH_RAW_DATA=$(pwd) && export BASEPATH_WORK_DATA=$(pwd)/output_data && export BASEPATH_PLOT=$(pwd)/figures && export G_PERT_LO_DIR=$BASEPATH_RAW_DATA/quenched_1.50Tc_zeuthenFlow/pert_LO && export tmppath="./correlators_flow/correlator_analysis/double_extrapolation/example_usage"


# Unzip the analysis scripts. This will create a new folder "correlators_flow" in the current directory.
tar -xzf correlators_flow.tar.gz
# OR
# git clone --branch v0.1.0 https://github.com/luhuhis/correlators_flow.git
# TODO update git

# Add the scripts folder to the PYTHONPATH.
export PYTHONPATH=$PYTHONPATH:$(pwd)/correlators_flow

# Unzip the AnalysisToolbox, which is a dependency of "correlators_flow". This will create a folder "AnalysisToobox"
tar -xzf AnalysisToolbox.tar.gz
# OR
# git clone https://github.com/LatticeQCD/AnalysisToolbox.git
# cd AnalysisToolbox
# git checkout f9eee73d45d7b981153db75cfaf2efa2b4cefa9c
# cd ..

# Add the AnalysisToolbox folder to the PYTHONPATH.
export PYTHONPATH=$PYTHONPATH:$(pwd)/AnalysisToolbox


# Create new directory for input data and change into it
mkdir input && cd input

# Unzip raw data (gradientFlow output from SIMULATeQCD). This will create a new folder "hisq_ms5_zeuthenFlow".
tar -xzf ../input_data.tar.gz

# Add the *parent* directory of the new folder (i.e., the "input" directory where you called the "tar" command) to the
# environment. This should usually be $(pwd).
# For example, if you called tar in the folder /home/data/input, then a new folder
# /home/data/input/hisq_ms5_zeuthenFlow has been created. The parent of that directory is then /home/data/input.
export BASEPATH_RAW_DATA=$(pwd)

# We also need to add another subdirectory to the environment
export G_PERT_LO_DIR=$BASEPATH_RAW_DATA/quenched_1.50Tc_zeuthenFlow/pert_LO

# Now, create two directories where you want to store the data of the intermediate steps
# (resampling, interpolation, extrapolations, ...) as well as final output data and figures.

cd ..
mkdir output_data
mkdir figures

export BASEPATH_WORK_DATA=$(pwd)/output_data
export BASEPATH_PLOT=$(pwd)/figures

# Specify how many processes can be created for parallelization.
export NPROC=1

# Note that the finished figures are also contained in "figures.tar.gz" and can just be extracted for convenience,
# skipping the whole data processing steps.
# tar -xzf figures.tar.gz # TODO

# Note that the finished output data is also contained in "output_data.tar.gz" and can just be extracted for
# convenience, skipping the whole processing steps.
# tar -xzf output_data.tar.gz # TODO


# 2. RUN SCRIPTS

# 2.1 Double extrapolation
#
# Files are either saved in plain text (.txt, .dat) or in numpy binary format (.npy).
# Some steps output median and standard deviation over all bootstrap samples as plain text files, and then the
# actual underlying bootstrap samples in numpy format (binary).

# IMPORTANT: If the files are not executable by default, just make every file in ./correlators_flow/
# executable using chmod +x. # TODO make this fix part of the readme

# 2.1.1 Merge individual small text files into large binary (npy) files
# Merge individual BB correlator measurement text files (output from SIMULATeQCD) into a small number of larger numpy files (binary format).
# Metadata is saved to text files. This can take some time, mostly depending on file system speed (with conventional HDDs it may take hours).
export tmppath="./correlators_flow/correlator_analysis/double_extrapolation/example_usage"
./$tmppath/1_merge_data.sh quenched_1.50Tc_zeuthenFlow BB $BASEPATH_RAW_DATA $BASEPATH_WORK_DATA

# Afterwards, the following files have been created in
# $BASEPATH_WORK_DATA/quenched_1.50Tc_zeuthenFlow/BB/<conftype>/

# flowtimes_<conftype>.dat              | flowtimes that were measured, corresponding to the data in the .npy files
# n_datafiles_<conftype>.dat            | metadata (number of files per stream, MCMC trajectory number, etc.)
# BB_imag_<conftype>_merged.npy         | merged raw data, imaginary part of BB correlator
# BB_real_<conftype>_merged.npy         | merged raw data, real part of BB correlator
# polyakov_imag_<conftype>_merged.npy   | merged raw data, imaginary part of polyakovloop
# polyakov_real_<conftype>_merged.npy   | merged raw data, real part of polyakovloop

# where <conftype> may be, for example,
# s144t36_b0754400 (meaning Ns=144, Nt=36, beta=7.544)

# 2.1.2 Resampling
# The double-extrapolation of the correlator data as well as the spectral reconstruction fits are later performed on
# each individual bootstrap sample.
#
# Load the merged data files, extract an equally spaced MCMC time series, then plot the time history of the Polyakov loop at a large flow time.
# Then bin configurations according to the integrated autocorrelation time, then perform bootstrap resampling of the uncorrelated blocks and save the samples to numpy files (binary format).
./$tmppath/2_reduce_data.sh quenched_1.50Tc_zeuthenFlow BB $BASEPATH_WORK_DATA $BASEPATH_PLOT $NPROC

# Afterward, the following files have been created in
# $BASEPATH_WORK_DATA/quenched_1.50Tc_zeuthenFlow/BB/<conftype>/

# BB_<conftype>_samples.npy   | bootstrap samples of BB correlator
# BB_flow_cov_<conftype>.npy  | flow time correlation matrix (based on bootstrap samples)
# BB_<conftype>.dat           | median BB correlator (useful for checking the data / plotting)
# BB_err_<conftype>.dat       | std_dev of BB correlator (useful for checking the data / plotting)

# and in
# $BASEPATH_PLOT/quenched_1.50Tc_zeuthenFlow/BB/<conftype>/

# polyakovloop_MCtime.pdf     | shows the MCMC time series of the polyakovloop at a large flow time

# 2.1.3 Interpolation
# Interpolate the BB correlator in Euclidean time and in flow time, such that a common set of normalized flow times
# is available across all lattices and temperatures.
./$tmppath/3_spline_interpolate.sh quenched_1.50Tc_zeuthenFlow BB $BASEPATH_WORK_DATA $BASEPATH_PLOT $NPROC

# Afterward, the following files have been created in
# $BASEPATH_WORK_DATA/quenched_1.50Tc_zeuthenFlow/BB/<conftype>/

# BB_<conftype>_relflows.txt                         | Normalized flow times (\sqrt{8_\tau_F}/\tau)
# BB_<conftype>_interpolation_relflow_mean.npy       | Median of the interpolated BB correlator for the corresponding normalized flow times (binary format, npy)
# BB_<conftype>_interpolation_relflow_samples.npy    | Interpolations of each individual bootstrap sample of the BB correlator for the corresponding normalized flow times (binary format, npy)

# and in
# $BASEPATH_PLOT/quenched_1.50Tc_zeuthenFlow/BB/<conftype>/

# BB_interpolation_relflow.pdf                       | Multi-page PDF containing plots of interpolation in Eucl. time at different normalized flow times
# BB_interpolation_relflow_combined.pdf              | Plot of interpolation in Eucl. time for two normalized flow times

# 2.1.4 Continuum extrapolation
# Take the continuum limit of the BB correlator using a combined fit on each sample
./$tmppath/4_continuum_extr.sh quenched_1.50Tc_zeuthenFlow BB $BASEPATH_WORK_DATA $BASEPATH_PLOT $NPROC

# Afterward, the following files have been created in
# $BASEPATH_WORK_DATA/quenched_1.50Tc_zeuthenFlow/BB/cont_extr/

# BB_cont_relflows.dat        | Normalized flow times at which the continuum extrapolation was done
# BB_cont_relflow.dat         | Median of continuum-extrapolated BB correlator at corresponding normalized flow times
# BB_cont_relflow_err.dat     | Std dev of continuum-extrapolated BB correlator at corresponding normalized flow times
# BB_cont_relflow_samples.npy | Continuum extrapolations of the BB correlator on each individual bootstrap sample

# and in
# $BASEPATH_PLOT/quenched_1.50Tc_zeuthenFlow/BB/

# BB_cont_quality_relflow.pdf | FIG 4 in the paper. Multipage PDF containing plots of continuum extrapolation of
# BB correlator for the corresponding normalized flow times

# 2.1.5 Coupling calculations
# Continuum extrapolation of the flow-scheme coupling measured on zero temperature lattices, TODO is this already reflected in the paper?
# and conversion from flow scheme to the MSBAR coupling at one scale, then perturbative 5-loop running to other relevant scales.

./correlators_flow/correlator_analysis/double_extrapolation/BB_renormalization/example_usage/extrapolate_coupling.sh $BASEPATH_RAW_DATA $BASEPATH_WORK_DATA $BASEPATH_PLOT

# Afterward, the following files have been created in
# $BASEPATH_WORK_DATA/quenched_1.50Tc_zeuthenFlow/coupling/

# g2_muF_cont_extr.txt        | Continuum-extrapolated flow scheme coupling (g^2) as function of flow scale mu_F=1/sqrt(8 tau_F) in units of temperature T
# g2_MSBAR_runFromMu_4.0.txt  | Continuum MS-BAR coupling (g^2) as a function of MS-BAR scale mu in units of temperature T (matched at mu_F/T=4 and then perturbatively run)

# and in
# $BASEPATH_PLOT/quenched_1.50Tc_zeuthenFlow/coupling/

# g2.pdf                      | Flow- and MSBAR-scheme couplings (g^2) as a function of scale in temperature units # TODO change this plot, decide what to show and what not
# g2_cont_extr.pdf            | Flow-scheme coupling (g^2) as a function of squared lattice spacing (= inverse Nt^2 at fixed temperature)


# 2.1.6 Carry out renormalization of G_B by computing \Zmatch(\mubarir,\mubaruv,\muflow)

./correlators_flow/correlator_analysis/double_extrapolation/BB_renormalization/example_usage/compute_Z.sh $BASEPATH_WORK_DATA $BASEPATH_PLOT

# TODO use new input_data.tar.gz

# Afterward, the following files have been created in
# $BASEPATH_WORK_DATA/quenched_1.50Tc_zeuthenFlow/coupling/

# Z_match_ref4.0_UVLO_IRLO.dat   | Z_match as a function of scale, with mu=... # TODO
# Z_match_ref4.0_UVNLO_IRLO.dat  | Z_match as a function of scale, with mu=... # TODO

# $BASEPATH_PLOT/quenched_1.50Tc_zeuthenFlow/coupling/

# integrand.pdf          | integrand plot # TODO remove
# Zref4.0_UVLO_IRLO.pdf  | individual Z factors for chice  # TODO remove
# Zref4.0_UVNLO_IRLO.pdf | individual Z factors  # TODO remove
# Z_total.pdf            | All considered versions of Z_match as a function of flow scale
# Z_total_flowtime.pdf   | All considered versions of Z_match as a function of flow radius 1/(8 tau_F) T

# 2.1.7 Flow-time-to-zero extrapolation
# Take the flow-time-to-zero limit of the BB correlator using a combined fit on each sample
./$tmppath/5_flowtime_extr.sh quenched_1.50Tc_zeuthenFlow BB $BASEPATH_WORK_DATA $BASEPATH_PLOT $NPROC

# Afterward, the following files have been created in
# $BASEPATH_WORK_DATA/quenched_1.50Tc_zeuthenFlow/BB/

# BB_flow_extr_relflow.npy | Flow-time-to-zero extrapolated continuum BB correlator for each bootstrap sample
# BB_flow_extr_relflow.txt | Median and std dev of flow-time-to-zero extrapolated continuum BB correlator

# and in
# $BASEPATH_PLOT/quenched_1.50Tc_zeuthenFlow/BB/

# BB_flow_extr_quality_relflow.pdf | FIG 5 in the paper.

# 2.1.8 Comparison plot: EE vs BB correlator after continuum and flow-time-to-zero extrapolations
./correlators_flow/correlator_analysis/plotting/plot_EEvsBB.py --inputfolder $BASEPATH_WORK_DATA/quenched_1.50Tc_zeuthenFlow/ --outputfolder $BASEPATH_PLOT/quenched_1.50Tc_zeuthenFlow/

# Afterward, the following files have been created in
# $BASEPATH_PLOT/quenched_1.50Tc_zeuthenFlow

# EEvsBB.pdf | Continuum- and flow-time-extrapolated color-magnetic and -electric correlators as a function of tau T at 1.5 Tc


# TODO: add EE/EE_flow_extr_relflow.txt to output data somehow...


# 2.2 Spectral reconstruction [OPTIONAL, this takes a lot of computing time, so the output files are already included]
./correlators_flow/spf_reconstruction/model_fitting/example_usage/spf_reconstruct.sh $BASEPATH_WORK_DATA NO $NPROC

# Afterward, the following files have been created in
# $BASEPATH_WORK_DATA/quenched_1.50Tc_zeuthenFlow/BB/spf/<model>_<rhoUV>_Nf3_T<Temp-in-MeV>_<min_scale>_<scale>_tauTgtr<min_tauT>_<suffix>
# and in
# $BASEPATH_WORK_DATA/quenched_1.50Tc_zeuthenFlow/BB/<conftype-with-smallest-lattice-spacing>/spf
# which reflect fits to the double-extrapolated BB correlator as well as finite lattice spacing and flow time fits,
# respectively.

# corrfit.dat            | Median input BB correlator and fitted model correlator
# params_samples.npy     | Model spectral function fit parameters for each bootstrap sample, as well as chisq/dof
# params.dat             | Median spectral function fit parameters and 34th percentiles
# phIUV.npy              | UV part of fitted model spectral function as function of omega/T (binary numpy format)
# samples.npy            | Copy of the input correlator bootstrap samples but multiplied by Gnorm (= actual fit input)
# spffit.npy             | Median spectral function with left/right 34th percentiles as function of omega/T


# 2.3 Create correlator plots at fixed normalized flow time
./correlators_flow/correlator_analysis/relative_flow/example_usage/Nf3TemperatureComparison_Paper.sh $BASEPATH_WORK_DATA $BASEPATH_PLOT

# Afterwards, the following files have been created in
# $BASEPATH_PLOT/quenched_1.50Tc_zeuthenFlow/BB/

# BB_relflow_hisq_final.pdf | Fig. 2 in the paper
# BB_relflow_hisq_0.25.pdf  | Top panel of Fig. 1 in the paper
# BB_relflow_hisq_0.30.pdf  | Bottom panel of Fig. 1 in the paper

# and in
# $BASEPATH_WORK_DATA/quenched_1.50Tc_zeuthenFlow/BB/<conftype>/relflow/

# BB_relflow_0.25.dat       | BB correlator at normalized flow time 0.25 (see Fig. 1)
# BB_relflow_0.30.dat       | BB correlator at normalized flow time 0.30 (see Fig. 1)


# 2.4.1 Plot spectral reconstruction fit results #1
./correlators_flow/spf_reconstruction/plot_fits/example_usage/plot_fits_hisq.sh $BASEPATH_WORK_DATA $BASEPATH_PLOT

# Afterwards, the following files have been created in
# $BASEPATH_PLOT/quenched_1.50Tc_zeuthenFlow/BB/T<T-in-MeV>/

# BB_spf_T<T-in-MeV>_paper.pdf      | Fig. 6 in the paper
# BB_corrfit_T<T-in-MeV>_paper.pdf  | Fig. 7 in the paper
# BB_kappa_T<T-in-MeV>_paper.pdf    | Fig. 9 in the paper

# and in
# $BASEPATH_WORK_DATA/quenched_1.50Tc_zeuthenFlow/BB/T<Temp-in-MeV>/

# BB_kappa_T<T-in_MeV>_paper.txt | Final results for kappa, see Table III in the paper

# and in
# $BASEPATH_WORK_DATA/quenched_1.50Tc_zeuthenFlow/BB/<conftype>/

# BB_kappa_T<T-in-MeV>_finiteflow_paper.txt | Final results for kappa from finite a and tau_F data, see Table III in the paper

# 2.4.2 Plot spectral reconstruction fit results #2
./correlators_flow/spf_reconstruction/plot_fits/example_usage/plot_kfactors.sh $BASEPATH_WORK_DATA $BASEPATH_PLOT

# Afterwards, the following files have been created in
# $BASEPATH_PLOT/quenched_1.50Tc_zeuthenFlow/BB/

# BB_kfactor.pdf | Fig. 8 in the paper


# 2.4.3 Plot spectral reconstruction fit results #3
./correlators_flow/spf_reconstruction/plot_fits/example_usage/plot_final_kappas.sh $BASEPATH_WORK_DATA $BASEPATH_PLOT

# Afterwards, the following files have been created in
# $BASEPATH_PLOT/quenched_1.50Tc_zeuthenFlow/BB/

# kappa_hisq_paper.pdf | Fig. 10 in the paper


# 2.4.4 Plot Figure 3 (2piTD)
tar -xzf figure_3_2piTD.tar.gz
# Please refer to the README file: ./figure_3_2piTD/README
cd ./figure_3_2piTD/
sh 2piTDs.sh
# Now, ./figure_3_2piTD/2piTDs.pdf exists.


# --------------------------------------------------------------------------------------------------------------------------------------------------------------

./1_merge_data.sh quenched_1.50Tc_zeuthenFlow BB /work/data/altenkort/gradientFlow ../../../../data/merged/

