#!/bin/bash
for file in \
../plots/quenched/s*t*_b*/*.pdf \
../plots/quenched/*.pdf \
../plots/quenched/single_flow/*.pdf \
../plots/quenched/single_flow_mass_shift/*.pdf \
../data_merged/quenched/continuum_limit/*.pdf
do
    pngfolder=${file%/*}/png
    mkdir -p $pngfolder
    outname=$(basename $file); outname=${outname%.pdf}
    pdftoppm $file $pngfolder/$outname -png -r 300 -singlefile &
done
wait

# for file in *.pdf ; do outname=$(basename $file); outname=${outname%.pdf} ; pdftoppm $file ./png/$outname -png -r 300 -singlefile ; done 



