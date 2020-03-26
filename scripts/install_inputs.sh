#!/bin/bash

rootdir=$1; shift

if [ -z "$rootdir" ]; then
    echo "Specify the location of the inputs"
    exit 1
fi

makelink() {
    base=$1; shift
    ext=$1; shift
    dir=$1; shift

    for stage in 6ad 8ad 7ad; do
        ln -rsf $rootdir/${base}_${stage}*.$ext ShapeFit/$dir/${base}_$stage.$ext
    done
}

makelink Theta13-inputs_P17B_inclusive txt Inputs
makelink accidental_eprompt_shapes root Spectra
makelink ibd_eprompt_shapes root Spectra
