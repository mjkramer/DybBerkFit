#!/bin/bash

BASE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"/..
DBG=${LBNL_FIT_DEBUG:+g}

cd $BASE/ShapeFit
root -b -q LoadClasses.C -e ".L make_evis_to_enu_matrix_fine_P17B.C+$DBG" \
    -e ".L build_covmatrix.C+$DBG" \
    -e ".L fit_shape_2d_P17B.C+$DBG"

cd $BASE/toySpectra
root -b -q LoadClasses.C -e ".L genToySpectraTree.C+$DBG" \
     -e ".L genEvisToEnuMatrix.C+$DBG" \
     -e ".L genSuperHistograms.C+$DBG" \
     -e ".L genPredictedIBD.C+$DBG"
