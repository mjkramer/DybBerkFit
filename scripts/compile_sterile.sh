#!/usr/bin/env bash


BASE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"/..
DBG=${LBNL_FIT_DEBUG:+g}

cd $BASE/ShapeFit
root -b -q LoadClasses.C -e ".L genOscProbTable.C+$DBG" \
    -e ".L fit_shape_3d.C+$DBG"
root -b -q LoadClasses.C -e ".L CLs/fit_shape_3d_CLs.C+$DBG" # \
    # -e ".L CLs/make_data_contours_CLs.C+$DBG"

cd $BASE/toySpectra
root -b -q LoadClasses.C -e ".L CLs/genToySpectraTree_parscans_CLs.C+$DBG"
