#!/bin/bash
for file in /home/altenkort/master/work/data_analysis/plots/quenched/s*t*_b*/*.pdf ; do
    pngfolder=${file%/*}/png
    rm -rf $pngfolder
    mkdir -p $pngfolder
    outname=$(basename $file); outname=${outname%.pdf}
    pdftoppm $file $pngfolder/$outname -png -r 300 -singlefile &
done


# for file in *.pdf ; do outname=$(basename $file); outname=${outname%.pdf} ; pdftoppm $file ./png/$outname -png -r 300 -singlefile ; done 



