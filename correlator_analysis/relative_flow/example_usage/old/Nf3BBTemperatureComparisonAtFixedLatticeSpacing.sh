#!/bin/bash

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