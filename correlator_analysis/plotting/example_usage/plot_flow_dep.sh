#!/bin/bash
conftypes_quenched=(
#"s064t16_b0687361"
"s080t20_b0703500"
"s096t24_b0719200"
"s120t30_b0739400"
"s144t36_b0754400"
)

# TODO add option to only plot EE or BB and decdie between quenched/hist


for conftype in "${conftypes_quenched[@]}" ; do
    for corr in "EE" "BB" ; do

        ../plot_flow_dependency.py \
         --qcdtype quenched_1.50Tc_zeuthenFlow \
         --corr $corr \
         --conftype $conftype \
         --basepath "/work/home/altenkort/work/correlators_flow/data/merged/" \
         --basepath_plot "/work/home/altenkort/work/correlators_flow/plots/" \
         --xlims 0 0.27 \
         --ylims 1 4.4 \
         --ticklocations 0.1 0.25 0.33 0.4 0.45 0.5 \
         --leg_pos 1 0.5 --leg_ncol 1 --leg_lw 0 --leg_pad 0.5 \
         --leg_loc "center left" \
         &

             ../plot_flow_dependency.py \
         --qcdtype quenched_1.50Tc_zeuthenFlow \
         --corr $corr \
         --conftype $conftype \
         --basepath "/work/home/altenkort/work/correlators_flow/data/merged/" \
         --basepath_plot "/work/home/altenkort/work/correlators_flow/plots/" \
         --xlims 0 0.11 \
         --ylims 1 4.4 \
         --ticklocations 0.1 0.2 0.25 0.3 \
         --leg_pos 1 0.5 --leg_ncol 1 --leg_lw 0 --leg_pad 0.5 \
         --leg_loc "center left" \
         --suffix "zoom" \
        &

    done
done


conftype=s096t24_b0824900_m002022_m01011
qcdtype=hisq_ms5_zeuthenFlow

../plot_flow_dependency.py \
 --qcdtype $qcdtype \
 --corr EE \
 --conftype $conftype \
 --basepath "/work/home/altenkort/work/correlators_flow/data/merged/" \
 --basepath_plot "/work/home/altenkort/work/correlators_flow/plots/" \
 --xlims 0 0.27 \
 --ylims 1 8 \
 --ticklocations 0.1 0.25 0.33 0.4 0.45 0.5 \
 --leg_pos 1 0.5 --leg_ncol 1 --leg_lw 0 --leg_pad 0.5 \
 --leg_loc "center left" \
&


     ../plot_flow_dependency.py \
 --qcdtype $qcdtype \
 --corr EE \
 --conftype $conftype \
 --basepath "/work/home/altenkort/work/correlators_flow/data/merged/" \
 --basepath_plot "/work/home/altenkort/work/correlators_flow/plots/" \
 --xlims 0 0.11 \
 --ylims 1 8 \
 --ticklocations 0.1 0.2 0.25 0.3 \
 --leg_pos 1 0.5 --leg_ncol 1 --leg_lw 0 --leg_pad 0.5 \
 --leg_loc "center left" \
 --suffix "zoom" \
&

wait
