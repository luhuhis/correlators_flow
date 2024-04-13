#!/bin/bash

## Data publication for "Lattice B-field correlators for heavy quarks", 2023, arXiv:2402.09337
## Authors:
##   Luis Altenkort, altenkort@physik.uni-bielefeld.de, Universität Bielefeld
##   David de la Cruz, TU Darmstadt
##   Olaf Kaczmarek, Universität Bielefeld
##   Guy D. Moore, TU Darmstadt
##   Hai-Tao Shu, Brookhaven National Laboratory

## Requirements
## - bash 5.0 (Linux shell)
## - python 3.11
##   - packages: matplotlib numpy scipy rundec natsort colorama
## - working latex installation with packages amsmath and mathtools for matplotlib tex rendering
## - gnuplot 5

## INSTRUCTIONS
## These instructions can be copied into a bash 5.0 shell.

## 1. Setup
export NPROC=20  # Change this accordingly to number of processors for parallelization
#tar -xzf correlators_flow.tar.gz  # Unzip the analysis scripts. This will create a new folder "correlators_flow" in the current directory.
## OR
# git clone --branch v0.1.0 https://github.com/luhuhis/correlators_flow.git  # TODO update git
git clone https://github.com/luhuhis/correlators_flow.git  # TODO update git
chmod -R +x ./correlators_flow  # Make all scripts executable

tar -xzf AnalysisToolbox.tar.gz # Unzip the AnalysisToolbox (python package), which is a dependency of "correlators_flow". This will create a folder "AnalysisToobox"
## OR
# git clone https://github.com/LatticeQCD/AnalysisToolbox.git
# ( cd AnalysisToolbox && git checkout f9eee73d45d7b981153db75cfaf2efa2b4cefa9c )

tar -xzf data.tar.gz # Unzip data (gradientFlow output from SIMULATeQCD) and create various folders.

## Create environment variables for paths
export PYTHONPATH=$(pwd)/correlators_flow:$(pwd)/AnalysisToolbox:${PYTHONPATH} BASEPATH_RAW_DATA=$(pwd)/input BASEPATH_WORK_DATA=$(pwd)/output_data BASEPATH_PLOT=$(pwd)/figures
export G_PERT_LO_DIR=${BASEPATH_RAW_DATA}/quenched_1.50Tc_zeuthenFlow/pert_LO

## Note that the finished figures are also contained in "figures.tar.gz" and can just be extracted for convenience,
## skipping the whole data processing steps.
# tar -xzf figures.tar.gz # TODO

## Note that the finished output data is also contained in "output_data.tar.gz" and can just be extracted for
## convenience, skipping the whole processing steps.
# tar -xzf output_data.tar.gz # TODO


## 2. RUN SCRIPTS
## 2.1 Double extrapolation

# Files are either saved in plain text (.txt, .dat) or in numpy binary format (.npy).
# Some steps output median and standard deviation over all bootstrap samples as plain text files, and then the
# actual underlying bootstrap samples in numpy format (binary).

## 2.1.1 Merge individual small text files into large binary (npy) files
# Merge individual BB correlator measurement text files (output from SIMULATeQCD) into a small number of larger numpy files (binary format).
# Metadata is saved to text files. This can take some time, mostly depending on file system speed (with conventional HDDs it may take hours).
./correlators_flow/correlator_analysis/double_extrapolation/example_usage/1_merge_data.sh quenched_1.50Tc_zeuthenFlow BB ${BASEPATH_RAW_DATA} ${BASEPATH_WORK_DATA}

# Afterward, the following files have been created in:
# $BASEPATH_WORK_DATA/quenched_1.50Tc_zeuthenFlow/BB/<conftype>/,
# where <conftype> may be, for example:
# s144t36_b0754400 (meaning Ns=144, Nt=36, beta=7.544).

# flowtimes_<conftype>.dat              | flowtimes that were measured, corresponding to the data in the .npy files
# n_datafiles_<conftype>.dat            | metadata (number of files per stream, MCMC trajectory number, etc.)
# BB_imag_<conftype>_merged.npy         | merged raw data, imaginary part of BB correlator
# BB_real_<conftype>_merged.npy         | merged raw data, real part of BB correlator
# polyakov_imag_<conftype>_merged.npy   | merged raw data, imaginary part of polyakovloop
# polyakov_real_<conftype>_merged.npy   | merged raw data, real part of polyakovloop

## 2.1.2 Resampling
# The double-extrapolation of the correlator data as well as the spectral reconstruction fits are later performed on
# each individual bootstrap sample.

# Load the merged data files, extract an equally spaced MCMC time series, then plot the time history of the Polyakov loop at a large flow time.
# Then bin configurations according to the integrated autocorrelation time, then perform bootstrap resampling of the uncorrelated blocks and save the samples to numpy files (binary format).
./correlators_flow/correlator_analysis/double_extrapolation/example_usage/2_reduce_data.sh quenched_1.50Tc_zeuthenFlow BB ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT} ${NPROC}

# Afterward, the following files have been created in
# $BASEPATH_WORK_DATA/quenched_1.50Tc_zeuthenFlow/BB/<conftype>/

# BB_<conftype>_samples.npy   | bootstrap samples of BB correlator
# BB_flow_cov_<conftype>.npy  | flow time correlation matrix (based on bootstrap samples)
# BB_<conftype>.dat           | median BB correlator (useful for checking the data / plotting)
# BB_err_<conftype>.dat       | std_dev of BB correlator (useful for checking the data / plotting)

# and in
# $BASEPATH_PLOT/quenched_1.50Tc_zeuthenFlow/BB/<conftype>/

# polyakovloop_MCtime.pdf     | shows the MCMC time series of the polyakovloop at a large flow time

## 2.1.3 Interpolation
# Interpolate the BB correlator in Euclidean time and in flow time, such that a common set of normalized flow times
# is available across all lattices and temperatures.
./correlators_flow/correlator_analysis/double_extrapolation/example_usage/3_spline_interpolate.sh quenched_1.50Tc_zeuthenFlow BB ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT} ${NPROC}

# Afterward, the following files have been created in
# $BASEPATH_WORK_DATA/quenched_1.50Tc_zeuthenFlow/BB/<conftype>/

# BB_<conftype>_relflows.txt                         | Normalized flow times (\sqrt{8_\tau_F}/\tau)
# BB_<conftype>_interpolation_relflow_mean.npy       | Median of the interpolated BB correlator for the corresponding normalized flow times (binary format, npy)
# BB_<conftype>_interpolation_relflow_samples.npy    | Interpolations of each individual bootstrap sample of the BB correlator for the corresponding normalized flow times (binary format, npy)

# and in
# $BASEPATH_PLOT/quenched_1.50Tc_zeuthenFlow/BB/<conftype>/

# BB_interpolation_relflow.pdf                       | Multi-page PDF containing plots of interpolation in Eucl. time at different normalized flow times
# BB_interpolation_relflow_combined.pdf              | Plot of interpolation in Eucl. time for two normalized flow times

## 2.1.4 Continuum extrapolation
# Take the continuum limit of the BB correlator using a fit on each sample
./correlators_flow/correlator_analysis/double_extrapolation/example_usage/4_continuum_extr.sh quenched_1.50Tc_zeuthenFlow BB ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT} ${NPROC}

# Afterward, the following files have been created in
# $BASEPATH_WORK_DATA/quenched_1.50Tc_zeuthenFlow/BB/cont_extr/

# BB_cont_relflows.dat        | Normalized flow times at which the continuum extrapolation was done
# BB_cont_relflow.dat         | Median of continuum-extrapolated BB correlator at corresponding normalized flow times
# BB_cont_relflow_err.dat     | Std dev of continuum-extrapolated BB correlator at corresponding normalized flow times
# BB_cont_relflow_samples.npy | Continuum extrapolations of the BB correlator on each individual bootstrap sample

# and in
# $BASEPATH_PLOT/quenched_1.50Tc_zeuthenFlow/BB/

# BB_cont_quality_relflow.pdf | Contains FIG 3 in the paper. Multipage PDF containing plots of continuum extrapolation of BB correlator for the corresponding normalized flow times.

## 2.1.5 Coupling calculations
# Continuum extrapolation of the flow-scheme coupling measured on zero temperature lattices,
# and conversion from flow scheme to the MSBAR scheme coupling at one scale, then perturbative 5-loop running to other relevant scales.

./correlators_flow/correlator_analysis/double_extrapolation/BB_renormalization/example_usage/extrapolate_coupling.sh ${BASEPATH_RAW_DATA} ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT}

# Afterward, the following files have been created in
# $BASEPATH_WORK_DATA/quenched_1.50Tc_zeuthenFlow/coupling/

# g2_muF_cont_extr.txt         | Continuum-extrapolated flow scheme coupling (g^2) as function of flow scale mu_F=1/sqrt(8 tau_F) in units of temperature T
# g2_MSBAR_runFromMu_6.28.txt  | Continuum MS-BAR coupling (g^2) as a function of MS-BAR scale mu in units of temperature T (matched at mu_F/T=2pi and then perturbatively run)

# and in
# $BASEPATH_PLOT/quenched_1.50Tc_zeuthenFlow/coupling/

# g2.pdf                      | Right panel of FIG 1 in the paper. Flow- and MSBAR-scheme couplings (g^2) as a function of scale in temperature units.
# g2_cont_extr.pdf            | Left panel of FIG 1 in the paper. Flow-scheme coupling (g^2) as a function of squared lattice spacing (= inverse Nt^2 at fixed temperature)

## 2.1.6 Carry out renormalization of G_B by computing Z_match

./correlators_flow/correlator_analysis/double_extrapolation/BB_renormalization/example_usage/compute_Z.sh ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT}

# Afterward, the following files have been created in
# $BASEPATH_WORK_DATA/quenched_1.50Tc_zeuthenFlow/coupling/

# Z_match_ref6.28_UVLO_IRNLO.dat    | Z_match as a function of scale, with \bar{\mu}_T/T=19.18, \bar{\mu}_{\tau_F}/mu_F=1
# Z_match_ref6.28_UVNLO_IRNLO.dat   | Z_match as a function of scale, with \bar{\mu}_T/T=19.18, \bar{\mu}_{\tau_F}/mu_F=1.50
# Z_match_ref6.28_UVLO_IRLO.dat     | Z_match as a function of scale, with \bar{\mu}_T/T=2piT, \bar{\mu}_{\tau_F}/mu_F=1
# Z_match_ref6.28_UVNLO_IRLO.dat    | Z_match as a function of scale, with \bar{\mu}_T/T=2piT, \bar{\mu}_{\tau_F}/mu_F=1.50

# and in
# $BASEPATH_PLOT/quenched_1.50Tc_zeuthenFlow/coupling/

# Z_total.pdf            | Left panel of FIG 2 in the paper. All considered versions of Z_match as a function of flow scale
# Z_total_flowtime.pdf   | Right panel of FIG 2 in the paper. All considered versions of Z_match as a function of flow radius 1/(8 tau_F) T

## 2.1.7 Flow-time-to-zero extrapolation
# Take the flow-time-to-zero limit of the BB correlator using a combined fit on each sample
./correlators_flow/correlator_analysis/double_extrapolation/example_usage/5_flowtime_extr.sh quenched_1.50Tc_zeuthenFlow BB ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT} ${NPROC}

# Afterward, the following files have been created in
# $BASEPATH_WORK_DATA/quenched_1.50Tc_zeuthenFlow/BB/

# BB_flow_extr_relflow_ref6.28_<scale-choices>.npy | Flow-time-to-zero extrapolated continuum BB correlator for each bootstrap sample
# BB_flow_extr_relflow_ref6.28_<scale-choices>.txt | Median and std dev of flow-time-to-zero extrapolated continuum BB correlator

# using the following scales choices:
#                                               | \bar{\mu}_T/T,  \bar{\mu}_{\tau_F}/mu_F
# BB_flow_extr_relflow_ref6.28_UVLO_IRLO.npy    | 2\pi            1
# BB_flow_extr_relflow_ref6.28_UVLO_IRLO.txt    |
#                                               |
# BB_flow_extr_relflow_ref6.28_UVLO_IRNLO.npy   | 19.18           1
# BB_flow_extr_relflow_ref6.28_UVLO_IRNLO.txt   |
#                                               |
# BB_flow_extr_relflow_ref6.28_UVNLO_IRLO.npy   | 2\pi            1.50
# BB_flow_extr_relflow_ref6.28_UVNLO_IRLO.txt   |
#                                               |
# BB_flow_extr_relflow_ref6.28_UVNLO_IRNLO.npy  | 19.18           1.50
# BB_flow_extr_relflow_ref6.28_UVNLO_IRNLO.txt  |


# and in
# $BASEPATH_PLOT/quenched_1.50Tc_zeuthenFlow/BB/

# BB_flow_extr_quality_no_extr.pdf | Left panel of FIG 4 in the paper. Bare continuum BB correlator as a function of flow time.
# BB_flow_extr_quality_relflow.pdf | Right panel of FIG 4 in the paper. Renormalized continuum BB correlator as a function of flow time with flow time extrapolation.

## 2.1.8 Comparison plot: EE vs BB correlator after continuum and flow-time-to-zero extrapolations
./correlators_flow/correlator_analysis/plotting/plot_EEvsBB.py --inputfolder ${BASEPATH_WORK_DATA}/quenched_1.50Tc_zeuthenFlow/ --outputfolder ${BASEPATH_PLOT}/quenched_1.50Tc_zeuthenFlow/

# Afterward, the following file has been created in
# $BASEPATH_PLOT/quenched_1.50Tc_zeuthenFlow

# EEvsBB.pdf | FIG 5 in the paper. Continuum- and flow-time-extrapolated color-magnetic and -electric correlators

## 2.2 Spectral reconstruction [OPTIONAL, this takes a lot of computing time, so the output files are already included]
./correlators_flow/spf_reconstruction/model_fitting/example_usage/spf_reconstruct.sh ${BASEPATH_WORK_DATA} NO ${NPROC}
# TODO create separate zip file for all output files

# Afterward, the following files have been created in
# $BASEPATH_WORK_DATA/quenched_1.50Tc_zeuthenFlow/BB/spf/<model>_<rho-UV-order>_Nf0_T0.472_<min_scale>_<running-scale-coefficient>_tauTgtr0.24_<suffix>

# corrfit.dat            | Median input BB correlator and fitted model correlator
# params_samples.npy     | Model spectral function fit parameters for each bootstrap sample, as well as chisq/dof
# params.dat             | Median spectral function fit parameters and 34th percentiles
# phIUV.npy              | UV part of fitted model spectral function as function of omega/T (binary numpy format)
# samples.npy            | Copy of the input correlator bootstrap samples but multiplied by Gnorm (= actual fit input)
# spffit.npy             | Median spectral function with left/right 34th percentiles as function of omega/T


## 2.2.1 Plot spectral reconstruction fit results for kappa_B
./correlators_flow/spf_reconstruction/plot_fits/example_usage/plot_fits_quenched.sh BB ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT} yes

# Afterward, the following files have been created in
# $BASEPATH_PLOT/quenched_1.50Tc_zeuthenFlow/BB/

# TODO: additional text file: /work/data/altenkort/gradientFlow/BB_paper_datapub_test/output_data/quenched_1.50Tc_zeuthenFlow/BB//BB_kappa_quenched_1.5Tc.txt

# BB_kappa_quenched_1.5Tc.pdf    | FIG 7 in the paper

## 2.2.2 Plot comparison to literature (this depends on the previous call to `plot_fits_quenched.sh`!)
./correlators_flow/spf_reconstruction/plot_fits/example_usage/plot_final_kappas.sh ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT} BB
# TODO /work/data/altenkort/gradientFlow/BB_paper_datapub_test/figures/quenched_1.50Tc_zeuthenFlow/BB//kappa_BB_quenched_literature.pdf
# Afterward, the following files have been created in
# $BASEPATH_PLOT/quenched_1.50Tc_zeuthenFlow/BB/
# kappa_BB_quenched_literature.pdf      | FIG 8 in the paper


## 2.2.3 Plot spectral reconstruction fit results for spectral function model shapes and model correlators
# Note that this script (`plot_fits_quenched.sh`) needs to be run again now with last argument being "no" instead of "yes":
./correlators_flow/spf_reconstruction/plot_fits/example_usage/plot_fits_quenched.sh BB ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT} no

# Afterward, the following files have been created in
# $BASEPATH_PLOT/quenched_1.50Tc_zeuthenFlow/BB/
# BB_spf_quenched_1.5Tc.pdf      | Left panel of FIG 6 in the paper
# BB_corrfit_quenched_1.5Tc.pdf  | Right panel of FIG 6 in the paper


## 4. Plot fit to g^2 and g^4
./correlators_flow/spf_reconstruction/plot_fits/publication_specific/fit.py ${BASEPATH_RAW_DATA}/quenched_1.50Tc_zeuthenFlow/literature/ ${BASEPATH_PLOT}/quenched_1.50Tc_zeuthenFlow/
# Afterward, the following files have been created in
# $BASEPATH_PLOT/quenched_1.50Tc_zeuthenFlow/BB/
# compare_kappa_g2.pdf                  | FIG 9 in the paper
