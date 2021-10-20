#!/bin/bash
for i in 1 2 ; do
  for argument in \
  "--conftype s096t32_b0824900_m002022_m01011" \
  "--conftype s096t32_b0824900_m002022_m01011 --part_obs numerator" \
  "--conftype s096t32_b0824900_m002022_m01011 --part_obs polyakovloop --custom_xlims3 -0.1 11.25" \
  "--conftype s064t64_b0824900_m002022_m01011" \
  "--conftype s064t64_b0824900_m002022_m01011 --part_obs numerator --rel_err_limit 0.5 --custom_ylims -35 10" \
  "--conftype s064t64_b0824900_m002022_m01011 --reconstruct " \
  "--conftype s064t64_b0824900_m002022_m01011 --reconstruct --conftype_2 s096t32_b0824900_m002022_m01011 --custom_ylims 0.7 1.5" \
  "--conftype s064t64_b0824900_m002022_m01011 --reconstruct --part_obs numerator --custom_ylims -0.1 0.32" \
  "--conftype s064t64_b0824900_m002022_m01011 --part_obs polyakovloop --custom_xlims3 -0.1 11.25 --rel_err_limit 2" \
  ; do
  ./1_plot_correlator.py --qcdtype hisq_b8249_zeuthenFlow --corr EE --wideaspect --only_plot_no 2 3 --rel_err_limit 0.2 --ignore_pert_lims --show_TauByA --show_flowtime_instead_of_flowradii --flowselectionfile ./supplements/validflowtimes_ms5_Nt32and64.dat --custom_xlims3 0 11.25 --tau_selection 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 $argument &
  done
  wait
done