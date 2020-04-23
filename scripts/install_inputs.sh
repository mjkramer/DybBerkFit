#!/bin/bash

rootdir=$1; shift

if [ -z "$rootdir" ]; then
    echo "Specify the location of the inputs"
    exit 1
fi

copy() {
    base=$1; shift
    ext=$1; shift
    dir=$1; shift

    for stage in 6ad 8ad 7ad; do
        cp --remove-destination $rootdir/${base}_${stage}*.$ext ShapeFit/$dir/${base}_$stage.$ext
    done
}

copy Theta13-inputs_P17B_inclusive txt Inputs
copy accidental_eprompt_shapes root Spectra
copy ibd_eprompt_shapes root Spectra
