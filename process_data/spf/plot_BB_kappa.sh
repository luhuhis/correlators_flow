#!/bin/bash

# ======= SYSTEMATICS

if [ 1 -eq 1 ] ; then
models_file=( plaw_wIR1.0_wUV6.28 pnorm2.0 max )

order="LO"
scale="min2piT_w1"
model_label=('``line"' '``smooth"' '``step"')
order_label="\\ \\ LO"

colors=(C0 C1 C2 C3)

filearray=()
labelarray=()
counter=3
posarray=()
colorarray=()
for idx in "${!models_file[@]}"; do
    filearray+=("spf/${models_file[idx]}_${order}_T0.472_${scale}_500_0.0_exp0_quenched_BB/params.dat")
    labelarray+=("${model_label[idx]}")
    colorarray+=(${colors[idx]})
    posarray+=($counter)
    counter=$((counter+1))
done
./plot_fits.py --outputpath ~/work/correlators_flow/plots/quenched_1.50Tc_zeuthenFlow/BB/spf/ \
                --file_basepath /work/home/altenkort/work/correlators_flow/data/merged/quenched_1.50Tc_zeuthenFlow/BB/ \
                --files ${filearray[@]} \
                --pos ${posarray[@]} \
                --labels "${labelarray[@]}" \
                --colors ${colorarray[@]} \
                --kappa_swap_axes \
                --ylims 0 6 \
                --xlims 0 2 \
                --obs kappa --suffix systematics --title "\$T=1.5\,T_c\$" --xlabel "\$\\kappa_B/T^3 \$" --ylabel " "
wait
fi

#--usetex \