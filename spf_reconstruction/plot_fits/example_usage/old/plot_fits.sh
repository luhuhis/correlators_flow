#!/bin/bash
if false; then
    ./plot_fits.py --outputpath ~/work/correlators_flow/plots/quenched_1.50Tc_zeuthenFlow/EE/spf/ \
        --file_basepath ~/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/EE/spf/ \
        --files \
        smax_NLO_20_0.24_exp0_quenched_1.5Tc_conttauF0_0.5opt_w0.5opt/spffit.dat \
        smax_NLO_20_0.24_exp0_quenched_1.5Tc_conttauF0_opt_wopt/spffit.dat \
        smax_NLO_20_0.24_exp0_quenched_1.5Tc_conttauF0_2opt_w2opt/spffit.dat \
        max_NLO_20_0.24_exp0_quenched_1.5Tc_conttauF0_0.5opt_w0.5opt/spffit.dat \
        max_NLO_20_0.24_exp0_quenched_1.5Tc_conttauF0_opt_wopt/spffit.dat \
        max_NLO_20_0.24_exp0_quenched_1.5Tc_conttauF0_2opt_w2opt/spffit.dat \
        --labels "smax 0.5opt" "smax opt" "smax 2opt" "max 0.5opt" "max opt" "max 2opt" \
        --obs spf --suffix comparison_scales --title 'NLO, $\mu = \mu_\mathrm{opt},  \omega_\mathrm{UV} = 2.2T$' &
fi

if false; then
    ./plot_fits.py --outputpath ~/work/correlators_flow/plots/quenched_1.50Tc_zeuthenFlow/EE/spf/ \
        --file_basepath ~/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/EE/spf/ \
        --files \
        max_NLO_20_0.24_exp0_quenched_1.5Tc_conttauF0_opt_wopt/spffit.dat \
        pnorm5_NLO_20_0.24_exp0_quenched_1.5Tc_conttauF0_opt_wopt/spffit.dat \
        pnorm3_NLO_20_0.24_exp0_quenched_1.5Tc_conttauF0_opt_wopt/spffit.dat \
        smax_NLO_20_0.24_exp0_quenched_1.5Tc_conttauF0_opt_wopt/spffit.dat \
        pnorm1.5_NLO_20_0.24_exp0_quenched_1.5Tc_conttauF0_opt_wopt/spffit.dat \
        pnorm1.0_NLO_20_0.24_exp0_quenched_1.5Tc_conttauF0_opt_wopt/spffit.dat \
        --labels "pnorminf (max)" "pnorm5" "pnorm3" "pnorm2 (smax)" "pnorm1.5" "pnorm1 (sum)" \
        --obs spf --suffix comparison_norm --title 'NLO, $\mu = \mu_\mathrm{opt}$' &
fi

if false ; then
for obs in spf corr; do
    ./plot_fits.py --outputpath ~/work/correlators_flow/plots/quenched_1.50Tc_zeuthenFlow/EE/spf/ \
        --file_basepath ~/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/EE/spf/ \
        --files \
        max_NLO_20_0.24_exp0_quenched_1.5Tc_conttauF0_opt_wopt/${obs}fit.dat \
        pnorm1.0_NLO_20_0.24_exp0_quenched_1.5Tc_conttauF0_opt_wopt/${obs}fit.dat \
        max_LO_20_0.24_exp0_quenched_1.5Tc_conttauF0_0.5opt_w0.5opt/${obs}fit.dat \
        max_LO_20_0.24_exp0_quenched_1.5Tc_conttauF0_opt_wopt/${obs}fit.dat \
        max_LO_T0.472_min2opt_w2opt_20_0.24_exp0_quenched_1.5Tc_conttauF0/${obs}fit.dat \
        pnorm1.0_LO_20_0.24_exp0_quenched_1.5Tc_conttauF0_0.5opt_w0.5opt/${obs}fit.dat \
        pnorm1.0_LO_20_0.24_exp0_quenched_1.5Tc_conttauF0_opt_wopt/${obs}fit.dat \
        pnorm1.0_LO_20_0.24_exp0_quenched_1.5Tc_conttauF0_2opt_w2opt/${obs}fit.dat \
        line_NLO_20_0.24_exp0_quenched_1.5Tc_conttauF0_opt_wopt_wIR0.4/${obs}fit.dat \
        line_NLO_20_0.24_exp0_quenched_1.5Tc_conttauF0_opt_wopt_wIR1/${obs}fit.dat \
        pnorm2.0_LO_T0.472_minopt_wopt_20_0.24_exp0_quenched_1.5Tc_conttauF0/${obs}fit.dat \
        line_wIR0.4_wUV3.0_NLO_T0.472_minopt_wopt_20_0.24_exp0_quenched_1.5Tc_conttauF0/${obs}fit.dat \
        line_wIR0.4_wUV3.5_NLO_T0.472_minopt_wopt_20_0.24_exp0_quenched_1.5Tc_conttauF0/${obs}fit.dat \
        --labels "NLO max" "NLO sum" "LO max, 0.5opt" "LO max, opt" "LO max, 2opt" "LO sum 0.5opt" "LO sum opt" "LO sum 2opt" "NLO line opt wIR=0.4" "NLO line opt wIR=1" "LO smax opt" "NLO line opt wIR0.4 wUV3" "NLO line opt wIR0.4 wUV3.5"\
        --obs ${obs} --suffix comparison_LO_NLO --title 'LO using different scales vs NLO at $\mu = \mu_\mathrm{opt}$' &
done
fi

# ======= SPF CORR KAPPA

if [ 1 -eq 2 ] ; then
for OmegaByT_UV in 3.14 6.28 ; do

# for output file names
models_str=( maxLO smaxLO "plawLO_wUV${OmegaByT_UV}" maxNLO smaxNLO "plawNLO_wUV${OmegaByT_UV}" )

# for input file names
models_file=( "max" "pnorm2.0" "plaw_wIR1.0_wUV${OmegaByT_UV}" )
order=("LO" "NLO")
scale=("min2piT_w1" "minopt_wopt")
filenames=("spffit" "spffit" "corrfit" "params" "params")
obs=("spf" "spf_phys" "corr" "kappa"  "chisqdof")
flows=("0.20" "0.25" "0.30")  #"0.15" "0.20"

# for plot title
scale_title=()

# --xlims 0.14 0.512 --ylims 1.8 9 \
#--ylims 0.001 0.1 --xlims 0.001 1
for idf in "${!flows[@]}"; do
    for idx in "${!models_file[@]}"; do
        for idy in "${!order[@]}" ; do
            for idz in "${!obs[@]}"; do
                ./plot_fits.py --outputpath ~/work/correlators_flow/plots/hisq_b8249_zeuthenFlow/EE/spf/ \
                --file_basepath ~/work/correlators_flow/data/merged/hisq_b8249_zeuthenFlow/EE/ \
                --files \
                s096t36_b0824900_m002022_m01011/spf/${models_file[idx]}_${order[idy]}_Nf3_T0.195_${scale[idy]}_500_0.35_exp0_hisq_nt36_f${flows[idf]}/${filenames[idz]}.dat \
                s096t32_b0824900_m002022_m01011/spf/${models_file[idx]}_${order[idy]}_Nf3_T0.220_${scale[idy]}_500_0.35_exp0_hisq_nt32_f${flows[idf]}/${filenames[idz]}.dat \
                s096t28_b0824900_m002022_m01011/spf/${models_file[idx]}_${order[idy]}_Nf3_T0.251_${scale[idy]}_500_0.35_exp0_hisq_nt28_f${flows[idf]}/${filenames[idz]}.dat \
                s096t24_b0824900_m002022_m01011/spf/${models_file[idx]}_${order[idy]}_Nf3_T0.293_${scale[idy]}_500_0.35_exp0_hisq_nt24_f${flows[idf]}/${filenames[idz]}.dat \
                s096t20_b0824900_m002022_m01011/spf/${models_file[idx]}_${order[idy]}_Nf3_T0.352_${scale[idy]}_500_0.35_exp0_hisq_nt20_f${flows[idf]}/${filenames[idz]}.dat \
                --Nts 36 32 28 24 20 \
                --pos 195 220 251 293 352 \
                --labels "0.195" "0.220" "0.251" "0.293" "0.352" --legtitle "\$T\$ in GeV"\
                --OmegaByT_UV ${OmegaByT_UV} \
                --obs ${obs[idz]} --suffix ${models_str[idx+idy*3]}_f${flows[idf]} --title "${models_str[idx+idy*3]} ${scale[idy]} f${flows[idf]}" \
                &
            done
        done
    done
done
done
wait
fi


# ======= SYSTEMATICS

if [ 1 -eq 1 ] ; then
models_str=( plawLO_2piT smaxLO maxLO plawNLO_2piT smaxNLO maxNLO )
models_file=( plaw_wIR1.0_wUV6.28 pnorm2.0 max )

order=("LO" "NLO")
scale=("min2piT_w1" "minopt_wopt")
# todo add labels here
model_label=('``line"' '``smooth"' '``step"')
order_label=("\\ \\ LO" "NLO")

filenames=("spffit" "corrfit" "params")
obs=("spf" "corr" "kappa")

flows=("0.20" "0.25")
colors=(C0 C1 C2 C3)

#"0.15" "0.20"
# --xlims 0.14 0.512 --ylims 1.8 9 \
#--ylims 0.001 0.1 --xlims 0.001 1
Nts=(20 24 28 32 36)  # ${models_str[idx+idy*3]} ${scale[idy]} f${flows[idf]}
temp=(352 293 251 220 195)
for Nt_idx in "${!Nts[@]}"; do
    filearray=()
    labelarray=()
    counter=1
    posarray=()
    Ntsarray=()
    colorarray=()
    for idx in "${!models_file[@]}"; do
        for idy in "${!order[@]}" ; do
            for idf in "${!flows[@]}"; do
                filearray+=("s096t${Nts[Nt_idx]}_b0824900_m002022_m01011/spf/${models_file[idx]}_${order[idy]}_Nf3_T0.${temp[Nt_idx]}_${scale[idy]}_500_0.35_exp0_hisq_nt${Nts[Nt_idx]}_f${flows[idf]}/params.dat")
                labelarray+=("${model_label[idx]}, ${order_label[idy]}, \\ \\ $ ${flows[idf]}$")
                colorarray+=(${colors[idx]})
                Ntsarray+=(${Nts[Nt_idx]})
                posarray+=($counter)
                counter=$((counter+1))

            done
        done
    done
    ./plot_fits.py --outputpath ~/work/correlators_flow/plots/hisq_b8249_zeuthenFlow/EE/spf/ \
                    --file_basepath ~/work/correlators_flow/data/merged/hisq_b8249_zeuthenFlow/EE/ \
                    --files ${filearray[@]} \
                    --Nts ${Ntsarray[@]} \
                    --pos ${posarray[@]} \
                    --labels "${labelarray[@]}" \
                    --colors ${colorarray[@]} \
                    --kappa_swap_axes \
                    --ylims 0 13 \
                    --xlims 0 15 \
                    --usetex \
                    --obs kappa --suffix systematics_Nt${Nts[Nt_idx]} --title "\$T=${temp[Nt_idx]}\, \mathrm{MeV}\$" --xlabel "\$\\kappa/T^3 \$" --ylabel " "
done
wait
fi

#2_NLO_s1_alpha_1_20_0.24_exp0_quenched_1.5Tc_conttauF0_opt_wopt/spffit.dat \
#2_NLO_s1_beta_1_20_0.24_exp0_quenched_1.5Tc_conttauF0_opt_wopt/spffit.dat \
#

# TRY smooth max with x^4

# "step" 'line, $\omega_\mathrm{IR}=1T$'

#step_NLO_20_0.24_exp0_quenched_1.5Tc_conttauF0_opt_wopt/spffit.dat \
#line_NLO_20_0.24_exp0_quenched_1.5Tc_conttauF0_opt_wopt_wIR1/spffit.dat \

#--PhiUV_file smax_NLO_20_0.24_exp0_quenched_1.5Tc_conttauF0_opt_wopt/phiUV.dat \

# HISQ data

#./plot_fits.py --outputpath ~/work/correlators_flow/plots/quenched_1.50Tc_zeuthenFlow/EE/spf/ \
#--file_basepath ~/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/EE/spf/ \
#--files \
#smax_NLO_20_0.24_exp0_quenched_1.5Tc_conttauF0_0.25opt_w0.25opt/spffit.dat \
#smax_NLO_20_0.24_exp0_quenched_1.5Tc_conttauF0_0.5opt_w0.5opt/spffit.dat \
#smax_NLO_20_0.24_exp0_quenched_1.5Tc_conttauF0_opt_wopt/spffit.dat \
#smax_NLO_20_0.24_exp0_quenched_1.5Tc_conttauF0_2opt_w2opt/spffit.dat \
#max_NLO_20_0.24_exp0_quenched_1.5Tc_conttauF0_0.25opt_w0.25opt/spffit.dat \
#max_NLO_20_0.24_exp0_quenched_1.5Tc_conttauF0_0.5opt_w0.5opt/spffit.dat \
#max_NLO_20_0.24_exp0_quenched_1.5Tc_conttauF0_opt_wopt/spffit.dat \
#max_NLO_20_0.24_exp0_quenched_1.5Tc_conttauF0_2opt_w2opt/spffit.dat \
#2_NLO_s1_alpha_1_20_0.24_exp0_quenched_1.5Tc_conttauF0_opt_wopt/spffit.dat \
#2_NLO_s1_beta_1_20_0.24_exp0_quenched_1.5Tc_conttauF0_opt_wopt/spffit.dat \
#2_NLO_s1_alpha_2_20_0.24_exp0_quenched_1.5Tc_conttauF0_opt_wopt/spffit.dat \
#2_NLO_s1_beta_2_20_0.24_exp0_quenched_1.5Tc_conttauF0_opt_wopt/spffit.dat \
#--labels "smax 0.25opt" "smax 0.5opt" "smax opt" "smax 2opt" "max 0.25opt" "max 0.5opt" "max opt" "max 2opt" "1, sin" "1, sin2" "2, sin" "2, sin2" \
#--obs spf --suffix comparison --title '$\mu = \mu_\mathrm{opt},  \omega_\mathrm{UV} = 2.2T$'
