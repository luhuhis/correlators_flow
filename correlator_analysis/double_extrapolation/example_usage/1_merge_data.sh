#!/bin/bash

qcdtype=$1
corr=$2
basepath_raw_data=$3
basepath_work_data=$4

# this path contains "reference" flow times which are useful if there are
# multiple lattice spacings for each temperature that have different flow times
flowradiusbasepath=$5

if [ -z "$qcdtype" ] || [ -z "$corr" ] || [ -z "$basepath_raw_data" ] || [ -z "$basepath_work_data" ] ; then
    echo "Usage: $0 qcdtype corr input_basepath output_basepath [flowradiusbasepath]"
    echo "choices for qcdtype: quenched_1.50Tc_zeuthenFlow hisq_ms5_zeuthenFlow"
    echo "choices for corr: EE BB EE_clover BB_clover"
    echo "Usage example: $0 hisq_ms5_zeuthenFlow EE /work/data/altenkort/gradientFlow ../../../../data/merged/ /home/altenkort/work/correlators_flow/data/merged/hisq_ms5_zeuthenFlow/EE/"
    exit
fi


if [ "$qcdtype" == quenched_1.50Tc_zeuthenFlow ] ; then
    arr_conftypes=("s064t16_b0687361" "s080t20_b0703500" "s096t24_b0719200" "s120t30_b0739400" "s144t36_b0754400")

    # here we need to adjust the file names via --acc_sts due to historical reasons
    if [ "$corr" == "EE" ] ; then
        acc_sts="--acc_sts acc0.000010_sts0.000010"
        add_args="--legacy --basepath ../../../../data/raw/"
    elif [ "$corr" == "BB" ] ; then
        acc_sts="--acc_sts sts0.150000"
        add_args="--basepath ../../../../data/raw/"
    fi

elif [ "$qcdtype" == hisq_ms5_zeuthenFlow ] ; then
    arr_conftypes=(
    "s064t20_b0803600" "s064t22_b0814700" "s096t24_b0824900_m002022_m01011"
    "s064t20_b0785700" "s064t24_b0806800" "s096t28_b0824900_m002022_m01011"
    "s064t20_b0770400" "s064t24_b0791300" "s096t32_b0824900_m002022_m01011"
    "s064t20_b0757000" "s064t24_b0777700" "s096t36_b0824900_m002022_m01011"
    "s096t20_b0824900_m002022_m01011"
)
    temps=(
        "293"
        "251"
        "220"
        "195"
        )
    add_args="--basepath $basepath_raw_data"
    # TODO check whether these files are used somewhere else or only here. maybe just assume they are not used anywhere, then test the entire pipeline
    flowradii_refs=()
    for temp in "${temps[@]}"; do
        # the file name for these is hardcoded, since this is anyway just an example
        reffile="$flowradiusbasepath/T${temp}/flowradii_T${temp}.dat"
        if [ ! -f "$reffile" ] ; then
            echo "WARNING: $reffile does not exist"
        fi
        flowradii_refs+=("--reference_flowradii $reffile")
    done
fi


for idx in "${!arr_conftypes[@]}" ; do
    conftype="${arr_conftypes[idx]}"

    if [ "$qcdtype" == hisq_ms5_zeuthenFlow ] ; then

        # for historical reasons the file names of the Nt=24 lattice are different
        if [ "$conftype" == "s096t24_b0824900_m002022_m01011" ] ; then
            acc_sts=""
        else
            acc_sts="--acc_sts sts0.150000"
        fi

        # throw away the first 2000 traj for Nt=36 due to thermalization.
        # we do that here instead of in the next step because for some reason many of the files below 2000 are broken...
        if [ "$conftype" == "s096t36_b0824900_m002022_m01011" ] ; then
            even_more_args="--min_conf_nr 2000"
        else
            even_more_args=""
        fi

        # for this lattice we don't have other lattice spacings at the same temperature, which is why we also don't need any common flowtime values.
        if [ "$conftype" == "s096t20_b0824900_m002022_m01011" ] ; then
            flowradiiref=""
        else
            flowradiiref=${flowradii_refs[idx / 3]}
        fi

    fi

    (
    echo "work dir: $(dirname $0)" && cd "$(dirname $0)" || exit

    mycmd="
    ../_1_merge_data.py \
        $flowradiiref \
        --output_basepath $basepath_work_data \
        $even_more_args $add_args \
        --qcdtype $qcdtype --corr $corr $acc_sts --conftype $conftype
        "

    # uncomment these lines to confirm each script call
#    echo "$mycmd"
#    echo -en "\n y/n? "
#    read -r input
#    if [ "$input" == "y" ]; then
        eval $mycmd
#    fi

    )
done
