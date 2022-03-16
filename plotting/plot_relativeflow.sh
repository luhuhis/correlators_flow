#!/bin/bash

# quenched
#TODO quenched plot!!!!




#HISQ

# comparison of different temps in lattice units
./plot_rec_corr_fixFlowBytauT.py --suffix temps_latt_units --npoints 100 --ylims 0 10 --xlims 0 19 --title 'HISQ, $ms/5$' --min_tauT 0.01 --plot_in_lattice_units \
--flowradiusBytauT 0.3 0.3 0.3 \
--qcdtype hisq_b8249_zeuthenFlow \
--corr EE EE EE \
--conftype s096t36_b0824900_m002022_m01011 s096t32_b0824900_m002022_m01011 s096t20_b0824900_m002022_m01011 \



# compare different scales
./plot_rec_corr_fixFlowBytauT.py --suffix scale_comp --npoints 100 --ylims 0 6 --xlims 0 0.505 --title 'HISQ, $ms/5$' --min_tauT 0.01 \
--flowradiusBytauT 0.3 0.3 0.3 \
--qcdtype hisq_b8249_zeuthenFlow \
--corr EE \
--conftype s096t20_b0824900_m002022_m01011 \
--deduce_fitparam_files --fitparam_basepath ~/work/correlators_flow/data/merged/hisq_b8249_zeuthenFlow/EE/spf/ \
--fitparam_files 3_a_260_1e-04_d_0.4_frankstein_0.300_nt20_piT_1.0 3_a_260_1e-04_d_0.4_frankstein_0.300_nt20_2piT_1.0 3_a_260_1e-04_d_0.4_frankstein_0.300_nt20_2piT_opt \
--PhiUV_basepath /home/altenkort/work/correlators_flow/data/merged/spf_coupling/ \
--PhiUV_files hisq_nt20_SPF_LO_Nf3_0.352_piT_1.0_smax.dat hisq_nt20_SPF_LO_Nf3_0.352_2piT_1.0_smax.dat hisq_nt20_SPF_LO_Nf3_0.352_2piT_opt_smax.dat \


