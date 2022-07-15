#!/bin/bash

# quenched
# general quenched comparison
if 'false' ; then
params=(" \
--output_suffix _quenched_0 \
--fitparam_files \
5_a_300_0.24_exp0_quenched_1.5Tc_cont_tauF0_relflow_0.20_fit \
--no_kappa_plot \
" " \
--output_suffix _quenched_1 \
--fitparam_files \
5_a_300_0.24_exp0_quenched_1.5Tc_cont_tauF0_relflow_0.20_fit \
5_a_300_1e-04_d_0.2_relflow_0.20_piT_1.0 \
--tauT_vlines 0.25 \
" " \
--output_suffix _quenched_2 \
--fitparam_files \
5_a_300_0.24_exp0_quenched_1.5Tc_cont_tauF0_relflow_0.20_fit \
5_a_300_1e-04_d_0.2_relflow_0.20_piT_1.0 \
5_a_300_1e-04_d_0.4_relflow_0.20_piT_1.0 \
--tauT_vlines 0.25 0.416 \
" " \
--output_suffix _quenched_3 \
--fitparam_files \
5_a_300_0.24_exp0_quenched_1.5Tc_cont_tauF0_relflow_0.20_fit \
5_a_300_1e-04_d_0.2_relflow_0.20_piT_1.0 \
5_a_300_1e-04_d_0.4_relflow_0.20_piT_1.0 \
5_a_300_1e-04_d_0.4_quenched_1.5Tc_Nt36_relflow_0.20_piT_1.0 \
--tauT_vlines 0.25 0.416 0.4167 \
--conftype s144t36_b0754400 \
" " \
--output_suffix _quenched_4 \
--fitparam_files \
5_a_300_0.24_exp0_quenched_1.5Tc_cont_tauF0_relflow_0.20_fit \
5_a_300_1e-04_d_0.2_relflow_0.20_piT_1.0 \
5_a_300_1e-04_d_0.4_relflow_0.20_piT_1.0 \
5_a_300_1e-04_d_0.4_quenched_1.5Tc_Nt36_relflow_0.20_piT_1.0 \
5_a_300_1e-04_d_0.4_quenched_1.5Tc_Nt30_relflow_0.20_piT_1.0 \
--conftype s144t36_b0754400 s120t30_b0739400 \
--tauT_vlines 0.25 0.416 0.4167 0.433 \
" " \
--output_suffix _quenched_5 \
--fitparam_files \
5_a_300_0.24_exp0_quenched_1.5Tc_cont_tauF0_relflow_0.20_fit \
5_a_300_1e-04_d_0.2_relflow_0.20_piT_1.0 \
5_a_300_1e-04_d_0.4_relflow_0.20_piT_1.0 \
5_a_300_1e-04_d_0.4_quenched_1.5Tc_Nt36_relflow_0.20_piT_1.0 \
5_a_300_1e-04_d_0.4_quenched_1.5Tc_Nt30_relflow_0.20_piT_1.0 \
5_a_500_0.38_exp2_quenched_1.5Tc_Nt20_relflow0.20_piT_1.0w \
--conftype s144t36_b0754400 s120t30_b0739400 s080t20_b0703500 \
--tauT_vlines 0.25 0.416 0.4167 0.433 \
")

for param in "${params[@]}" ; do
./plot_rec_corr_fixFlowBytauT.py  --npoints 100 --ylims 2.5 3.75 \
$param \
--xlims 0.2 0.512 --kappa_xlims -0.1 2.5 \
--title ", $\\:\\: \\mu=\\sqrt{[\\pi T\\,]^2+\\omega^2}$" \
--no_connection \
--model 5 \
--plot_quenched_extr \
--flowradiusBytauT 0.2 \
--qcdtype quenched_1.50Tc_zeuthenFlow \
--corr EE \
--deduce_fitparam_files --fitparam_basepath ~/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/EE/spf/ \
--PhiUV_basepath ~/work/correlators_flow/data/merged/spf_coupling/ \
--PhiUV_files quenched_1.5Tc_SPF_LO_Nf0_0.472_piT_1.0_smax.dat \
--output_path ~/work/correlators_flow/plots/quenched_1.50Tc_zeuthenFlow/EE/ \
--color_data k C0 C1 C2 C3 \
--color_fit k C0 C9 C1 C2 C3 \
--color_vlines k C9 C1 C2 \
--kappa_ypos 1.75 1.45 1.3 1 0.85 0.6 \
--leg_n_dummies 0 \
--min_tauT 0.250 0.250 0.416 0.416 0.433 \
--min_tauT_plot 0.05 \
--no_label \
--no_just_UV \
--leg_pos 0.15 1 \
&
done


./plot_rec_corr_fixFlowBytauT.py  --npoints 100 --ylims 2.5 3.75 \
--output_suffix _quenched_5_IMP \
--conftype s144t36_b0754400 s080t20_b0703500 \
--tauT_vlines 0.25 0.416 0.4167 0.433 \
--xlims 0.2 0.512 --kappa_xlims -0.1 2.5 \
--title ", $\\:\\: \\mu=\\sqrt{[\\pi T\\,]^2+\\omega^2}$" \
--no_connection \
--model 5 \
--plot_quenched_extr \
--flowradiusBytauT 0.2 \
--qcdtype quenched_1.50Tc_zeuthenFlow \
--corr EE \
--deduce_fitparam_files --fitparam_basepath ~/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/EE/spf/ \
--PhiUV_basepath ~/work/correlators_flow/data/merged/spf_coupling/ \
--PhiUV_files quenched_1.5Tc_SPF_LO_Nf0_0.472_piT_1.0_smax.dat \
--output_path ~/work/correlators_flow/plots/quenched_1.50Tc_zeuthenFlow/EE/ \
--color_data k C0 C1 C3 \
--color_fit k C0 C9 C1 C3 \
--color_vlines k C9 C1 C2 \
--kappa_ypos 1.75 1.45 1.3 1 0.6 \
--leg_n_dummies 0 \
--min_tauT 0.250 0.250 0.416 0.433 \
--min_tauT_plot 0.05 \
--no_label \
--no_just_UV \
--leg_pos 0.15 1 \

fi

# ================
# ================ HISQ EE
# ================
wUV=3.14
models_suffix=( maxLO smaxLO lineLO maxNLO smaxNLO lineNLO )
models_file=( max pnorm2.0 line_wIR1.0_wUV${wUV} )
models=( max pnorm line )
modelparams=( "" "--p 2" "--OmegaByT_IR 1 --OmegaByT_UV ${wUV}" )
order=("LO" "NLO")
scale=("min2piT_w1" "minopt_wopt")
minscale=('2\pi T\,' '\mathrm{opt}')
runscale=('\omega' '\mathrm{opt}')
# --xlims 0.14 0.512  \
#
flows=("0.25" "0.30")  # "0.15" "0.20"
suffix=("" "_phys")
# 0 1 2 ,   0 1
# --ylims 0 9 \
#--kappa_ypos 1.45 1.3 1.15 1 0.85
xlims=("--xlims 0 0.505" "--plot_in_lattice_units --xlims 0 17")
#0.8 1.2
for idz in {0..1} ; do
    for idx in "${!models_file[@]}"; do
        for idy in "${!order[@]}" ; do
            for idf in "${!flows[@]}" ; do
                ./plot_rec_corr_fixFlowBytauT.py \
                --output_suffix _hisq_${models_suffix[idx+idy*3]}_f${flows[idf]}${suffix[idz]} --npoints 100  \
                --output_path ~/work/correlators_flow/plots/hisq_b8249_zeuthenFlow/EE/ \
                ${xlims[idz]} \
                --kappa_in_GeV \
                --ylims 0.6 1.4 \
                --norm_by_Gansatz \
                --kappa_xaxis_temp --kappa_xlims 0 0.5 --kappa_labelheight 0.65 \
                --kappa_temps 196 220 251 296 352 \
                --title " " \
                --model ${models[idx]} ${modelparams[idx]} \
                --flowradiusBytauT "${flows[idf]}" \
                --qcdtype hisq_b8249_zeuthenFlow --corr EE \
                --conftype \
                s096t36_b0824900_m002022_m01011 \
                s096t32_b0824900_m002022_m01011 \
                s096t28_b0824900_m002022_m01011 \
                s096t24_b0824900_m002022_m01011 \
                s096t20_b0824900_m002022_m01011 \
                --deduce_fitparam_files --fitparam_basepath ~/work/correlators_flow/data/merged/hisq_b8249_zeuthenFlow/EE/ \
                --fitparam_files \
                s096t36_b0824900_m002022_m01011/spf/${models_file[idx]}_${order[idy]}_Nf3_T0.196_${scale[idy]}_500_0.35_exp0_hisq_nt36_f${flows[idf]} \
                s096t32_b0824900_m002022_m01011/spf/${models_file[idx]}_${order[idy]}_Nf3_T0.220_${scale[idy]}_500_0.35_exp0_hisq_nt32_f${flows[idf]} \
                s096t28_b0824900_m002022_m01011/spf/${models_file[idx]}_${order[idy]}_Nf3_T0.251_${scale[idy]}_500_0.35_exp0_hisq_nt28_f${flows[idf]} \
                s096t24_b0824900_m002022_m01011/spf/${models_file[idx]}_${order[idy]}_Nf3_T0.296_${scale[idy]}_500_0.35_exp0_hisq_nt24_f${flows[idf]} \
                s096t20_b0824900_m002022_m01011/spf/${models_file[idx]}_${order[idy]}_Nf3_T0.352_${scale[idy]}_500_0.35_exp0_hisq_nt20_f${flows[idf]} \
                --PhiUV_basepath ~/work/correlators_flow/data/merged/hisq_b8249_zeuthenFlow/EE/ \
                --PhiUV_files \
                s096t36_b0824900_m002022_m01011/spf/${models_file[idx]}_${order[idy]}_Nf3_T0.196_${scale[idy]}_500_0.35_exp0_hisq_nt36_f${flows[idf]}/phiUV.dat \
                s096t32_b0824900_m002022_m01011/spf/${models_file[idx]}_${order[idy]}_Nf3_T0.220_${scale[idy]}_500_0.35_exp0_hisq_nt32_f${flows[idf]}/phiUV.dat \
                s096t28_b0824900_m002022_m01011/spf/${models_file[idx]}_${order[idy]}_Nf3_T0.251_${scale[idy]}_500_0.35_exp0_hisq_nt28_f${flows[idf]}/phiUV.dat \
                s096t24_b0824900_m002022_m01011/spf/${models_file[idx]}_${order[idy]}_Nf3_T0.296_${scale[idy]}_500_0.35_exp0_hisq_nt24_f${flows[idf]}/phiUV.dat \
                s096t20_b0824900_m002022_m01011/spf/${models_file[idx]}_${order[idy]}_Nf3_T0.352_${scale[idy]}_500_0.35_exp0_hisq_nt20_f${flows[idf]}/phiUV.dat \
                --min_scale "${minscale[idy]}" --run_scale "${runscale[idy]}" \
                --tauT_vlines 0.35 \
                --min_tauT 0.35 0.35 0.35 0.35 0.35 --min_tauT_plot 0.05 \
                --no_label --no_connection \
                --leg_title_suffix ",\\: T\\," --leg_pos 0.15 1 --leg_n_dummies 0 \
                --leg_label_suffix  ",\\: 196 \\,\\mathrm{MeV}" ",\\: 220\\,\\mathrm{MeV} " ",\\: 251\\,\\mathrm{MeV}" ",\\: 296\\,\\mathrm{MeV}" ",\\: 352\\,\\mathrm{MeV}" \
                &
            done
        done
    done
done

if [ 1 -eq 2 ] ; then
# ================
# ================ HISQ BB
# ================
flows=( 0.20 0.25 0.30 )
xlims=("--xlims 0 0.505" "--plot_in_lattice_units --xlims 0 19")
suffix=("" "_phys")
corrs=("EE" "BB")
# todo loop over flows here
# \
#0.8 1.2
for idz in {0..2} ; do
    for idy in {0..1} ; do
        for idx in {0..1}; do
            ./plot_rec_corr_fixFlowBytauT.py \
            --output_suffix _hisq${suffix[idy]}_f${flows[idz]} \
            --output_path ~/work/correlators_flow/plots/hisq_b8249_zeuthenFlow/BB/ \
            ${xlims[idy]} \
            --ylims 0 10  \
            --title " " \
            --no_kappa_plot \
            --flowradiusBytauT "${flows[idz]}" \
            --qcdtype hisq_b8249_zeuthenFlow --corr ${corrs[idx]} \
            --conftype \
            s096t36_b0824900_m002022_m01011 \
            s096t32_b0824900_m002022_m01011 \
            s096t28_b0824900_m002022_m01011 \
            s096t24_b0824900_m002022_m01011 \
            s096t20_b0824900_m002022_m01011 \
            --tauT_vlines 0.35 \
            --min_tauT 0.35 0.35 0.35 0.35 0.35 --min_tauT_plot 0.05 \
            --no_label --no_connection \
            --leg_title_suffix ",\\: T\\," --leg_pos 0.15 1 --leg_n_dummies 0 \
            --leg_label_suffix  ",\\: 196 \\,\\mathrm{MeV}" ",\\: 220\\,\\mathrm{MeV} " ",\\: 251\\,\\mathrm{MeV}" ",\\: 296\\,\\mathrm{MeV}" ",\\: 352\\,\\mathrm{MeV}"
        done
    done
done
fi
#
# \
#
#
#
#

if 'false' ; then

# comparison of different temps in lattice units
./plot_rec_corr_fixFlowBytauT.py --output_suffix temps_latt_units_0.2 --npoints 100 --ylims 0 10 --xlims 0 19 --title 'HISQ, $ms/5$' --min_tauT 0.01 --plot_in_lattice_units \
--flowradiusBytauT 0.2 0.2 0.2 \
--qcdtype hisq_b8249_zeuthenFlow \
--corr EE \
--conftype s096t36_b0824900_m002022_m01011 s096t32_b0824900_m002022_m01011 s096t20_b0824900_m002022_m01011 \
--model 5 \
&

./plot_rec_corr_fixFlowBytauT.py --output_suffix temps_latt_units_0.3 --npoints 100 --ylims 0 10 --xlims 0 19 --title 'HISQ, $ms/5$' --min_tauT 0.01 --plot_in_lattice_units \
--flowradiusBytauT 0.3 0.3 0.3 \
--qcdtype hisq_b8249_zeuthenFlow \
--corr EE \
--conftype s096t36_b0824900_m002022_m01011 s096t32_b0824900_m002022_m01011 s096t20_b0824900_m002022_m01011 \
--model 5 --min_tauT_plot 0.05 \
&


# compare different scales
./plot_rec_corr_fixFlowBytauT.py --output_suffix scale_comp --npoints 100 --ylims 0 6 --xlims 0 0.505 --title 'HISQ, $ms/5$' --min_tauT 0.01 \
--flowradiusBytauT 0.3 0.3 0.3 \
--qcdtype hisq_b8249_zeuthenFlow \
--corr EE \
--conftype s096t20_b0824900_m002022_m01011 \
--deduce_fitparam_files --fitparam_basepath ~/work/correlators_flow/data/merged/hisq_b8249_zeuthenFlow/EE/spf/ \
--fitparam_files 3_a_260_1e-04_d_0.4_frankstein_0.300_nt20_piT_1.0 3_a_260_1e-04_d_0.4_frankstein_0.300_nt20_2piT_1.0 3_a_260_1e-04_d_0.4_frankstein_0.300_nt20_2piT_opt \
--PhiUV_basepath /home/altenkort/work/correlators_flow/data/merged/spf_coupling/ \
--PhiUV_files hisq_nt20_SPF_LO_Nf3_0.352_piT_1.0_smax.dat hisq_nt20_SPF_LO_Nf3_0.352_2piT_1.0_smax.dat hisq_nt20_SPF_LO_Nf3_0.352_2piT_opt_smax.dat \
--model 5 --min_tauT_plot 0.05 \
&

# compare temps at two rel flowtimes
./plot_rec_corr_fixFlowBytauT.py --output_suffix temp_comp_0.2_0.3 --npoints 100 --ylims 0 10 --xlims 0 0.505 --title 'HISQ, $ms/5$' --min_tauT 0.01 \
--flowradiusBytauT 0.3 0.3 0.3 \
--qcdtype hisq_b8249_zeuthenFlow \
--corr EE \
--conftype s096t20_b0824900_m002022_m01011 \
--deduce_fitparam_files --fitparam_basepath ~/work/correlators_flow/data/merged/hisq_b8249_zeuthenFlow/EE/spf/ \
--fitparam_files 3_a_260_1e-04_d_0.4_frankstein_0.300_nt20_piT_1.0 3_a_260_1e-04_d_0.4_frankstein_0.300_nt20_2piT_1.0 3_a_260_1e-04_d_0.4_frankstein_0.300_nt20_2piT_opt \
--PhiUV_basepath /home/altenkort/work/correlators_flow/data/merged/spf_coupling/ \
--PhiUV_files hisq_nt20_SPF_LO_Nf3_0.352_piT_1.0_smax.dat hisq_nt20_SPF_LO_Nf3_0.352_2piT_1.0_smax.dat hisq_nt20_SPF_LO_Nf3_0.352_2piT_opt_smax.dat \
--model 5 --min_tauT_plot 0.05 \
&

fi

wait