# Data publication for "The diffusion of heavy quarks from lattice QCD", PhD thesis by Luis Altenkort, 2024, Bielefeld University

## TODO

- complete all zip files
- correct the figure numbers, add appendix numbers
- add descriptions for EE files
- add poetry install instructions
- create new release once everything else is done

## Requirements

- Bash 5.0 (Linux shell)
- Python poetry, or manually install packages according to `pyproject.toml`
- LaTeX installation with packages `amsmath` and `mathtools`
- gnuplot 5

## Instructions

Download all files from the data publication, move them to a new folder and navigate into that folder.
The following instructions can then copied into a bash 5.0 shell **in the same order as they appear below**. For convenience, all commands in this document are consolidated in `do_everything_thesis.sh`, which can be run to do everything with one script call.

## Unzip or clone analysis scripts

This will create a new folder `correlators_flow` in the current directory:

``` shell
# tar -xzf correlators_flow.tar.gz   
# OR
git clone --branch thesis https://github.com/luhuhis/correlators_flow.git  
```

Make all scripts executable:s

``` shell
chmod -R +x ./correlators_flow 
```

Unzip or clone the AnalysisToolbox (custom python package), which is a dependency of "correlators_flow". This will create a folder "AnalysisToobox":

``` shell
# tar -xzf AnalysisToolbox.tar.gz 
# OR
git clone https://github.com/LatticeQCD/AnalysisToolbox.git
( cd AnalysisToolbox && git checkout f9eee73d45d7b981153db75cfaf2efa2b4cefa9c )
```

## Unzip data

This extracts the gradientFlow output from SIMULATeQCD and creates various folders.

``` shell
tar -xzf data.tar.gz 
```

Note that the finished figures are also contained in "figures.tar.gz" and can just be extracted for convenience,
skipping the whole data processing steps:

``` shell
tar -xzf figures.tar.gz
```

Note that the finished output data is also contained in "output_data.tar.gz" and can just be extracted for
convenience, skipping the whole processing steps:

``` shell
tar -xzf output_data.tar.gz
```

## Setup python and shell environment

Set environment variables:

``` shell
export PYTHONPATH=$(pwd)/correlators_flow:$(pwd)/AnalysisToolbox:${PYTHONPATH} BASEPATH_RAW_DATA=$(pwd)/input BASEPATH_WORK_DATA=$(pwd)/output_data BASEPATH_PLOT=$(pwd)/figures
export G_PERT_LO_DIR=${BASEPATH_RAW_DATA}/quenched_1.50Tc_zeuthenFlow/pert_LO
```

Change this accordingly to number of processors for parallelization:

``` shell
export NPROC=20  
```

Install the python environment:
``` shell
poetry install
```

## Note on output files

Files are either saved in plain text (.txt, .dat) or in numpy binary format (.npy).
Some steps output median and standard deviation over all bootstrap samples as plain text files, and then the
actual underlying bootstrap samples in numpy format.

In the following, these variables are often used in file names:

- `<qcdtype>` is either `quenched_1.50Tc_zeuthenFlow` or `hisq_ms5_zeuthenFlow`
- `<conftype>` may be, for example, `s144t36_b0754400` (meaning Ns=144, Nt=36, beta=7.544)
- `<corr>` is either `EE` or `BB`

## Create figures of perturbative correlators

Create **Figure 4.1** at `${BASEPATH_PLOT}/EE_QED_LPT.pdf`:

```shell
./correlators_flow/perturbative_corr/plot_QED_LPT.py --inputfolder ${BASEPATH_RAW_DATA} --outputfolder ${BASEPATH_PLOT}
```

Create **Figure 5.2** at `${BASEPATH_PLOT}/pertLO/EE_pert_contvslatt_flow.pdf`:

```shell
./correlators_flow/perturbative_corr/plot_pert_correlators.py --Ntau 24 --inputfolder ${BASEPATH_RAW_DATA}/quenched_1.50Tc_zeuthenFlow/pert_LO/ --outputfolder ${BASEPATH_PLOT}/pertLO
```

Create **Figure 5.3** at `${BASEPATH_PLOT}/pertLO//pert_latt_comparison_EE_Nt30_<tau>.pdf` with `tau=5` and `tau=10`:

```shell
./correlators_flow/perturbative_corr/plot_tree_level_imp.py --Nt 30 --corr EE --flowtime_file ${BASEPATH_RAW_DATA}/quenched_1.50Tc_zeuthenFlow/pert_LO/flowtimes.dat --outputpath ${BASEPATH_PLOT}/pertLO/ --inputpath ${BASEPATH_RAW_DATA}/quenched_1.50Tc_zeuthenFlow/pert_LO/ --tau 5
./correlators_flow/perturbative_corr/plot_tree_level_imp.py --Nt 30 --corr EE --flowtime_file ${BASEPATH_RAW_DATA}/quenched_1.50Tc_zeuthenFlow/pert_LO/flowtimes.dat --outputpath ${BASEPATH_PLOT}/pertLO/ --inputpath ${BASEPATH_RAW_DATA}/quenched_1.50Tc_zeuthenFlow/pert_LO/ --tau 10
```

## Merge correlator measurement files

Merge individual small text files into larger binary numpy files.

Merge individual correlator measurement text files (output from SIMULATeQCD) into a small number of larger numpy files (binary format).
Metadata is saved to text files. This can take some time, mostly depending on file system speed (with slow HDDs it may take hours).

```shell
./correlators_flow/correlator_analysis/double_extrapolation/example_usage/1_merge_data.sh quenched_1.50Tc_zeuthenFlow EE ${BASEPATH_RAW_DATA} ${BASEPATH_WORK_DATA}
./correlators_flow/correlator_analysis/double_extrapolation/example_usage/1_merge_data.sh quenched_1.50Tc_zeuthenFlow BB ${BASEPATH_RAW_DATA} ${BASEPATH_WORK_DATA}
./correlators_flow/correlator_analysis/double_extrapolation/example_usage/1_merge_data.sh hisq_ms5_zeuthenFlow EE ${BASEPATH_RAW_DATA} ${BASEPATH_WORK_DATA} ${BASEPATH_RAW_DATA}/hisq_ms5_zeuthenFlow/reference_flowtimes
```

Afterward, the following files have been created in
`$BASEPATH_WORK_DATA/<qcdtype>/<corr>/<conftype>/`:

| File | Comment |
| --- | --- |
| `flowtimes_<conftype>.dat`              | flowtimes that were measured, corresponding to the data in the .files |
| `n_datafiles_<conftype>.dat`           | metadata (number of files per stream, MCMC trajectory number, etc.) |
| `<corr>_imag_<conftype>_merged.npy`     | merged raw data, imaginary part of correlator |
| `<corr>_real_<conftype>_merged.npy`     | merged raw data, real part of correlator |
| `polyakov_imag_<conftype>_merged.npy`  | merged raw data, imaginary part of polyakovloop |
| `polyakov_real_<conftype>_merged.npy`   | merged raw data, real part of polyakovloop |

## Bootstrap resampling of uncorrelated blocks

Reminder: The double-extrapolation of the correlator data as well as the spectral reconstruction fits are later performed on each individual bootstrap sample.

Load the merged data files, extract an equally spaced MCMC time series, then plot the time history of the Polyakov loop at a large flow time.
Then bin configurations according to the integrated autocorrelation time, then perform bootstrap resampling of the uncorrelated blocks and save the samples to numpy files (binary format).

```shell
./correlators_flow/correlator_analysis/double_extrapolation/example_usage/2_reduce_data.sh quenched_1.50Tc_zeuthenFlow EE ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT} ${NPROC}
./correlators_flow/correlator_analysis/double_extrapolation/example_usage/2_reduce_data.sh quenched_1.50Tc_zeuthenFlow BB ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT} ${NPROC}
./correlators_flow/correlator_analysis/double_extrapolation/example_usage/2_reduce_data.sh hisq_ms5_zeuthenFlow EE ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT} ${NPROC}
```

Afterward, the following files have been created in
`$BASEPATH_WORK_DATA/<qcdtype>/<corr>/<conftype>/`

| File | Comment |
| --- | --- |
| `<corr>_<conftype>_samples.npy`   | Bootstrap samples of correlator |
| `<corr>_flow_cov_<conftype>.npy`  | Flow time correlation matrix (based on bootstrap samples) |
| `<corr>_<conftype>.dat`           | Median correlator (useful for checking the data / plotting) |
| `<corr>_err_<conftype>.dat`       | Std dev of correlator (useful for checking the data / plotting) |

and, in `$BASEPATH_PLOT/<qcdtype>/<corr>/<conftype>/`

| File | Comment |
| --- | --- |
| `polyakovloop_MCtime.pdf`     | **Figure 6.1 and 7.1.** Shows the MCMC time series of the polyakovloop at a large flow time. |

## Plot lattice spacing effects

Create **Figure 6.2** at `${BASEPATH_PLOT}/quenched_1.50Tc_zeuthenFlow/EE/EE_latt_effects.pdf`:

```shell
./correlators_flow/correlator_analysis/plotting/example_usage/2_plot_lateffects.sh quenched_1.50Tc_zeuthenFlow EE ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT} ${NPROC}
```

Optional: create the same figures for the 2+1-flavor cases at `${BASEPATH_PLOT}/hisq_ms5_zeuthenFlow/EE/T<T-in-MeV>/EE_latt_effects.pdf`:

```shell
./correlators_flow/correlator_analysis/plotting/example_usage/2_plot_lateffects.sh hisq_ms5_zeuthenFlow EE ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT} ${NPROC}
```

## Plot flow time dependency

Create **Figure 6.3**, **Figure 6.13**, and **Figure 7.3**:

```shell
./correlators_flow/correlator_analysis/plotting/example_usage/plot_flow_dep.sh quenched_1.50Tc_zeuthenFlow EE ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT}
./correlators_flow/correlator_analysis/plotting/example_usage/plot_flow_dep.sh quenched_1.50Tc_zeuthenFlow BB ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT}
./correlators_flow/correlator_analysis/plotting/example_usage/plot_flow_dep.sh hisq_ms5_zeuthenFlow EE ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT}
```

The figures can be found at:

```shell
${BASEPATH_PLOT}/quenched_1.50Tc_zeuthenFlow/EE/s080t20_b0703500/EE_s080t20_b0703500_flow_dep.pdf
${BASEPATH_PLOT}/quenched_1.50Tc_zeuthenFlow/EE/s080t20_b0703500/EE_s080t20_b0703500_flow_depzoom.pdf
${BASEPATH_PLOT}/quenched_1.50Tc_zeuthenFlow/BB/s080t20_b0703500/BB_s080t20_b0703500_flow_dep.pdf
${BASEPATH_PLOT}/quenched_1.50Tc_zeuthenFlow/BB/s080t20_b0703500/BB_s080t20_b0703500_flow_depzoom.pdf
${BASEPATH_PLOT}/hisq_ms5_zeuthenFlow/EE/s096t24_b0824900_m002022_m01011/EE_s096t24_b0824900_m002022_m01011_flow_dep.pdf
${BASEPATH_PLOT}/hisq_ms5_zeuthenFlow/EE/s096t24_b0824900_m002022_m01011/EE_s096t24_b0824900_m002022_m01011_flow_depzoom.pdf
```

## Interpolation

Interpolate the correlator in Euclidean time and in flow time, such that a common set of normalized flow times
is available across all lattices and temperatures.

```shell
./correlators_flow/correlator_analysis/double_extrapolation/example_usage/3_spline_interpolate.sh quenched_1.50Tc_zeuthenFlow EE ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT} ${NPROC}
./correlators_flow/correlator_analysis/double_extrapolation/example_usage/3_spline_interpolate.sh quenched_1.50Tc_zeuthenFlow BB ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT} ${NPROC}
./correlators_flow/correlator_analysis/double_extrapolation/example_usage/3_spline_interpolate.sh hisq_ms5_zeuthenFlow EE ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT} ${NPROC}
```

Afterward, the following files have been created in
`$BASEPATH_WORK_DATA/<qcdtype>/<corr>/<conftype>/`:

| File | Comment |
| --- | --- |
| `<corr>_<conftype>_relflows.txt`                         | Normalized flow times $\sqrt{8_{\tau_F}}/\tau$ |
| `<corr>_<conftype>_interpolation_relflow_mean.npy`       | Median of the interpolated correlator for the corresponding normalized flow times (binary format, npy) |
| `<corr>_<conftype>_interpolation_relflow_samples.npy`    | Interpolations of each individual bootstrap sample of the correlator for the corresponding normalized flow times (binary format, npy) |

and, in `$BASEPATH_PLOT/<qcdtype>/<corr>/<conftype>/`:

| File | Comment |
| --- | --- |
| `<corr>_interpolation_relflow.pdf`                       | **Figure 6.4**. Multi-page PDF containing plots of interpolation in Eucl. time at different normalized flow times. |
| `<corr>_interpolation_relflow_combined.pdf`              | **Figure 6.5**. Plot of interpolation in Eucl. time for two normalized flow times. |

## Continuum extrapolation

Take the continuum limit of the correlators using a fit on each sample

```shell
./correlators_flow/correlator_analysis/double_extrapolation/example_usage/4_continuum_extr.sh quenched_1.50Tc_zeuthenFlow EE ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT} ${NPROC}
./correlators_flow/correlator_analysis/double_extrapolation/example_usage/4_continuum_extr.sh quenched_1.50Tc_zeuthenFlow BB ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT} ${NPROC}
./correlators_flow/correlator_analysis/double_extrapolation/example_usage/4_continuum_extr.sh hisq_ms5_zeuthenFlow EE ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT} ${NPROC}
```

Afterward, the following files have been created in
`$BASEPATH_WORK_DATA/<qcdtype>/<corr>/cont_extr/`:

| File | Comment |
| --- | --- |
| `<corr>_cont_relflows.dat`        | Normalized flow times at which the continuum extrapolation was done |
| `<corr>_cont_relflow.dat`         | Median of continuum-extrapolated correlator at corresponding normalized flow times |
| `<corr>_cont_relflow_err.dat`     | Std dev of continuum-extrapolated correlator at corresponding normalized flow times |
| `<corr>_cont_relflow_samples.npy` | Continuum extrapolations of the correlator on each individual bootstrap sample |

and in `$BASEPATH_PLOT/<qcdtype>/<corr>/`:

| File | Comment |
| --- | --- |
| `<corr>_cont_quality_relflow.pdf` | **Figure 6.6**, **Figure 7.4**. This is a multipage PDF containing plots of continuum extrapolation of correlator for the corresponding normalized flow times. |

## Plot flow time correlations

Create **Figure 6.7** at 

```shell
./correlators_flow/correlator_analysis/plotting/plot_flow_correlations.py --qcdtype quenched_1.50Tc_zeuthenFlow --corr EE --conftype s144t36_b0754400 --basepath ${BASEPATH_WORK_DATA} --outputfolder ${BASEPATH_PLOT}/quenched_1.50Tc_zeuthenFlow/EE/ --nproc ${NPROC}
```

## Coupling calculations

Continuum extrapolation of the flow-scheme coupling measured on zero temperature lattices,
and conversion from flow scheme to the MSBAR scheme coupling at one scale, then perturbative 5-loop running to other relevant scales.

```shell
./correlators_flow/correlator_analysis/double_extrapolation/BB_renormalization/example_usage/extrapolate_coupling.sh ${BASEPATH_RAW_DATA} ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT} 6.40
```

Afterward, the following files have been created in
`$BASEPATH_WORK_DATA/quenched_1.50Tc_zeuthenFlow/coupling/`

| File | Comment |
| --- | --- |
| `g2_muF_cont_extr.txt`         | Continuum-extrapolated flow scheme coupling (g^2) as function of flow scale mu_F=1/sqrt(8 tau_F) in units of temperature T |
| `g2_MSBAR_runFromMu_6.28.txt`  | Continuum MS-BAR coupling (g^2) as a function of MS-BAR scale mu in units of temperature T (matched at mu_F/T=2pi and then perturbatively run) |

and, in `$BASEPATH_PLOT/quenched_1.50Tc_zeuthenFlow/coupling/`

| File | Comment |
| --- | --- |
| `g2.pdf`                      | Right panel of FIG 1 in the paper. Flow- and MSBAR-scheme couplings (g^2) as a function of scale in temperature units. |
| `g2_cont_extr.pdf`            | Left panel of FIG 1 in the paper. Flow-scheme coupling (g^2) as a function of squared lattice spacing (= inverse Nt^2 at fixed temperature) |


## Carry out renormalization of $G_B$ by computing $Z$

```shell
./correlators_flow/correlator_analysis/double_extrapolation/BB_renormalization/example_usage/compute_Z.sh ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT}
```

Afterward, the following files have been created in
`$BASEPATH_WORK_DATA/quenched_1.50Tc_zeuthenFlow/coupling/`

| File | Comment |
| --- | --- |
| `Z_match_ref6.28_UVLO_IRNLO.dat`    | Z_match as a function of scale, with \bar{\mu}_T/T=19.18, \bar{\mu}_{\tau_F}/mu_F=1 |
| `Z_match_ref6.28_UVNLO_IRNLO.dat`   | Z_match as a function of scale, with \bar{\mu}_T/T=19.18, \bar{\mu}_{\tau_F}/mu_F=1.50 |
| `Z_match_ref6.28_UVLO_IRLO.dat`     | Z_match as a function of scale, with \bar{\mu}_T/T=2piT, \bar{\mu}_{\tau_F}/mu_F=1 |
| `Z_match_ref6.28_UVNLO_IRLO.dat`    | Z_match as a function of scale, with \bar{\mu}_T/T=2piT, \bar{\mu}_{\tau_F}/mu_F=1.50 |

and, in `$BASEPATH_PLOT/quenched_1.50Tc_zeuthenFlow/coupling/`

| File | Comment |
| --- | --- |
| `Z_total.pdf`            | Left panel of FIG 2 in the paper. All considered versions of Z_match as a function of flow scale |
| `Z_total_flowtime.pdf`   | Right panel of FIG 2 in the paper. All considered versions of Z_match as a function of flow radius 1/(8 tau_F) T |

## Flow-time-to-zero extrapolation of $G_E$

```shell
./correlators_flow/correlator_analysis/double_extrapolation/example_usage/5_flowtime_extr.sh quenched_1.50Tc_zeuthenFlow EE ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT} ${NPROC}
```

Take the flow-time-to-zero limit of the EE correlator using a combined fit on each sample:

```shell
./correlators_flow/correlator_analysis/double_extrapolation/example_usage/5_flowtime_extr.sh hisq_ms5_zeuthenFlow EE ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT} ${NPROC}
```

## Flow-time-to-zero extrapolation of renormalized $G_B$

```shell
./correlators_flow/correlator_analysis/double_extrapolation/example_usage/5_flowtime_extr.sh quenched_1.50Tc_zeuthenFlow BB ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT} ${NPROC}
```

Afterward, the following files have been created in
`$BASEPATH_WORK_DATA/quenched_1.50Tc_zeuthenFlow/BB/`

| File | Comment | Scale |
| --- | --- | --- |
| ` BB_flow_extr_relflow_ref6.28_<scale-choices>.npy` | Flow-time-to-zero extrapolated continuum BB correlator for each bootstrap sample |  |
| `BB_flow_extr_relflow_ref6.28_<scale-choices>.txt` | Median and std dev of flow-time-to-zero extrapolated continuum BB correlator | |

using the following scales choices:

| File |$\bar{\mu}_T/T$ |  $\bar{\mu}_{\tau_F}/mu_F$ |
| --- | --- | --- |
| `BB_flow_extr_relflow_ref6.28_UVLO_IRLO.npy`    | $2\pi$         |   1 |
| `BB_flow_extr_relflow_ref6.28_UVLO_IRLO.txt`    | $2\pi$         |   1 |
   | |
| `BB_flow_extr_relflow_ref6.28_UVLO_IRNLO.npy`   | 19.18 |          1 |
| `BB_flow_extr_relflow_ref6.28_UVLO_IRNLO.txt`   | 19.18 |          1 |
|                                                | | |
| `BB_flow_extr_relflow_ref6.28_UVNLO_IRLO.npy`   | $2\pi$ |           1.50 |
| `BB_flow_extr_relflow_ref6.28_UVNLO_IRLO.txt`   | $2\pi$ |           1.50 |
|                                               | | |
| `BB_flow_extr_relflow_ref6.28_UVNLO_IRNLO.npy`  | 19.18 |          1.50 |
| `BB_flow_extr_relflow_ref6.28_UVNLO_IRNLO.txt`  | 19.18 |          1.50 |

and, in `$BASEPATH_PLOT/quenched_1.50Tc_zeuthenFlow/BB/`

| File | Comment |
| --- | --- |
| `BB_flow_extr_quality_no_extr.pdf` | Left panel of FIG 4 in the paper. Bare continuum BB correlator as a function of flow time. |
| `BB_flow_extr_quality_relflow.pdf` | Right panel of FIG 4 in the paper. Renormalized continuum BB correlator as a function of flow time with flow time extrapolation. |

## Comparison with multi-level results

```shell
./correlators_flow/multi-level/cont_extr_new.py --basepath ${BASEPATH_RAW_DATA}
./correlators_flow/correlator_analysis/plotting/6_plot_finalcorr.py --outputfolder ${BASEPATH_PLOT}/quenched_1.50Tc_zeuthenFlow/EE/ --input_flow ${BASEPATH_WORK_DATA}/quenched_1.50Tc_zeuthenFlow/EE/EE_flow_extr_relflow.txt --input_multilvl ${BASEPATH_RAW_DATA}/multi-level_2015/EE_2015_new.txt
```

TODO add compare EE vs Multilevel

## Compare final $G_E$ and $G_B$

```shell
./correlators_flow/correlator_analysis/plotting/plot_EEvsBB.py --inputfolder ${BASEPATH_WORK_DATA}/quenched_1.50Tc_zeuthenFlow/ --outputfolder ${BASEPATH_PLOT}/quenched_1.50Tc_zeuthenFlow/
```

Afterward, the following file has been created in
`$BASEPATH_PLOT/quenched_1.50Tc_zeuthenFlow`

| File | Comment |
| --- | --- |
| `EEvsBB.pdf` | FIG 5 in the paper. Continuum- and flow-time-extrapolated color-magnetic and -electric correlators |

## Create figures that illustrate spectral function models and reconstruction process

```shell
./correlators_flow/spf_reconstruction/plot_fits/example_usage/plot_g2.sh ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT}
```

```shell
./correlators_flow/spf_reconstruction/plotting/plot_integrand.py --outputpath ${BASEPATH_PLOT} --Nf 0 --min_scale eff --T_in_GeV 0.472 --omega_prefactor "1" --order LO --corr EE --mu_IR_by_T 1
```

## Spectral reconstruction 

Note: this takes a lot of computing time, so the output files are already included.

```shell
./correlators_flow/spf_reconstruction/model_fitting/example_usage/spf_reconstruct.sh quenched_1.50Tc_zeuthenFlow EE  ${BASEPATH_WORK_DATA} NO ${NPROC}
./correlators_flow/spf_reconstruction/model_fitting/example_usage/spf_reconstruct.sh quenched_1.50Tc_zeuthenFlow BB  ${BASEPATH_WORK_DATA} NO ${NPROC}
./correlators_flow/spf_reconstruction/model_fitting/example_usage/spf_reconstruct.sh hisq_ms5_zeuthenFlow EE  ${BASEPATH_WORK_DATA} NO ${NPROC}
```

Afterward, the following files have been created in
`$BASEPATH_WORK_DATA/<qcdtype>/<corr>/spf/<model>_<rho-UV-order>_Nf<nf>_T<T-in-MeV>_<min_scale>_<running-scale-coefficient>_tauTgtr0.24_<suffix>`

| File | Comment |
| --- | --- |
| `corrfit.dat`            | Median input BB correlator and fitted model correlator
| `params_samples.npy`     | Model spectral function fit parameters for each bootstrap sample, as well as chisq/dof |
| `params.dat`             | Median spectral function fit parameters and 34th percentiles
| `phIUV.npy`              | UV part of fitted model spectral function as function of omega/T (binary numpy format) |
| `samples.npy`            | Copy of the input correlator bootstrap samples but multiplied by Gnorm (= actual fit input) |
| `spffit.npy`             | Median spectral function with left/right 34th percentiles as function of omega/T |

## Plot spectral reconstruction fit results for quenched QCD

```shell
./correlators_flow/spf_reconstruction/plot_fits/example_usage/plot_fits_quenched.sh EE ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT} yes
./correlators_flow/spf_reconstruction/plot_fits/example_usage/plot_fits_quenched.sh BB ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT} yes
```

Afterward, the following files have been created in
`$BASEPATH_PLOT/quenched_1.50Tc_zeuthenFlow/<corr>/`

| File | Comment |
| --- | --- |
| `<corr>_kappa_quenched_1.5Tc.pdf`    | FIG 7 in the paper |

and, in `$BASEPATH_WORK_DATA/quenched_1.50Tc_zeuthenFlow/<corr>`

| File | Comment |
| --- | --- |
| `<corr>_kappa_quenched_1.5Tc.txt`  |  Final kappa value |

## Plot comparison to literature

Note: this especially depends on the previous call to `plot_fits_quenched.sh`

```shell
./correlators_flow/spf_reconstruction/plot_fits/example_usage/plot_final_kappas.sh ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT} EE
./correlators_flow/spf_reconstruction/plot_fits/example_usage/plot_final_kappas.sh ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT} BB
./correlators_flow/spf_reconstruction/plot_fits/example_usage/plot_final_kappas.sh ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT} hisq_thesis
```

Afterward, the following files have been created in
`$BASEPATH_PLOT/quenched_1.50Tc_zeuthenFlow/BB/`

| File | Comment |
| --- | --- |
| kappa_BB_quenched_literature.pdf      | FIG 8 in the paper |

## Plot spectral reconstruction fit results for spectral function model shapes and model correlators

Note that this script (`plot_fits_quenched.sh`) needs to be run again now with last argument being "no" instead of "yes":

```shell
./correlators_flow/spf_reconstruction/plot_fits/example_usage/plot_fits_quenched.sh BB ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT} no
```

Afterward, the following files have been created in
`$BASEPATH_PLOT/quenched_1.50Tc_zeuthenFlow/BB/`

| File | Comment |
| --- | --- |
| BB_spf_quenched_1.5Tc.pdf      | Left panel of FIG 6 in the paper |
| BB_corrfit_quenched_1.5Tc.pdf  | Right panel of FIG 6 in the paper |

## Plot fit to $g^2$ and $g^4$

```shell
./correlators_flow/spf_reconstruction/plot_fits/publication_specific/2024-BB-paper/fit_kappa_to_g2_g4.py --outputpath ${BASEPATH_PLOT}/quenched_1.50Tc_zeuthenFlow/
```

Afterward, the following files have been created in
`$BASEPATH_PLOT/quenched_1.50Tc_zeuthenFlow/BB/`

| File | Comment |
| --- | --- |
| `compare_kappa_g2.pdf`           | FIG 9 in the paper |

# HISQ only

## Create correlator plots at fixed normalized flow time

```shell
./correlators_flow/correlator_analysis/relative_flow/example_usage/Nf3TemperatureComparison_Paper.sh ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT}
```

Afterwards, the following files have been created inside
`$BASEPATH_PLOT/hisq_ms5_zeuthenFlow/EE/`

| File | Comment |
| --- | --- |
| EE_relflow_hisq_final.pdf | Fig. 2 in the paper
| EE_relflow_hisq_0.25.pdf  | Top panel of Fig. 1 in the paper
|  EE_relflow_hisq_0.30.pdf  | Bottom panel of Fig. 1 in the paper

and, in `$BASEPATH_WORK_DATA/hisq_ms5_zeuthenFlow/EE/<conftype>/relflow/`

| File | Comment |
| --- | --- |
| EE_relflow_0.25.dat       | EE correlator at normalized flow time 0.25 (see Fig. 1) |
| EE_relflow_0.30.dat       | EE correlator at normalized flow time 0.30 (see Fig. 1) |

```shell
./correlators_flow/spf_reconstruction/plot_fits/example_usage/plot_fits_hisq.sh ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT}
```

```shell
./correlators_flow/spf_reconstruction/plot_fits/example_usage/plot_kfactors.sh ${BASEPATH_WORK_DATA} ${BASEPATH_PLOT}
```

```shell
./correlators_flow/spf_reconstruction/plotting/plot_2piTD.py --outputfolder ${BASEPATH_PLOT}
```
