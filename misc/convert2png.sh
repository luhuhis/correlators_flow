#!/bin/bash
for file in ../plots/*/*/*.pdf ; do
    pngfolder=${file%/*}/png
    mkdir -p $pngfolder
    outname=$(basename $file); outname=${outname%.pdf}
    pdftoppm $file $pngfolder/$outname -png -r 300 -singlefile 
done
