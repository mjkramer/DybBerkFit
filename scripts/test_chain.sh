#!/bin/bash

# NOTE: Run prep_dirs.sh and install_inputs.sh first

step=$1; shift
step=${step:-all}

# RECOMPILE=1

BASE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"/..

ROOT=$(which root)
root() {
    test $RECOMPILE && $BASE/toySpectra/clean.sh
    $ROOT "$@"
}

# Compile stuff in advance to avoid race conditions when we parallelize
precompile() {
    cd $BASE/ShapeFit
    root -b -q LoadClasses.C -e '.L build_covmatrix.C+'

    cd $BASE/toySpectra
    root -b -q LoadClasses.C -e '.L genToySpectraTree.C+'
}

precompile

# echo "Using IHEP fast-n spectrum (see Config.h)"




# -------------------------- Generate reactor spectra --------------------------
genReactor() {
    cd $BASE/ReactorPowerCalculator
    root -b -q 'Produce_Isotope_SpectraP17B_unblinded.C(1)'
    cd isotope_spectra_by_Beda
    root -b -q make_combined_spectra_P17B_unblinded.C
}

# ------------------------ Generate ToyMC config files -------------------------
genToyConf() {
    cd $BASE/toySpectra/data_file
    ./generate_data_file.py
}

# --------------------------- Generate ToyMC samples ---------------------------
genToys() {
    cd $BASE/toySpectra
    ## sigsys:
    root -b -q LoadClasses.C -e '.L genToySpectraTree.C+' 'rungenToySpectraTree.C(2)' &
    # sleep 60
    ## bgsys
    root -b -q LoadClasses.C -e '.L genToySpectraTree.C+' 'rungenToySpectraTree.C(3)' &
    wait
}

# -------------------------- Generate evis/enu matrix --------------------------
genEvisEnu() {
    cd $BASE/toySpectra
    root -b -q LoadClasses.C genEvisToEnuMatrix.C+
    cd ../ShapeFit
    root -b -q LoadClasses.C make_evis_to_enu_matrix_fine_P17B.C+
}

# ------------------------- Generate super histograms --------------------------
genSuperHists() {
    cd $BASE/toySpectra
    root -b -q LoadClasses.C genSuperHistograms.C+ #&
    # sleep 20
}

# --------------------------- Generate PredictedIBD ----------------------------
genPredIBD() {
    cd $BASE/toySpectra
    root -b -q LoadClasses.C genPredictedIBD.C+ #&
    # wait
}

# ------------------------ Generate covariance matrices ------------------------
genCovMat() {
    cd $BASE/ShapeFit
    ## sigsys
    root -b -q LoadClasses.C -e '.L build_covmatrix.C+' 'run_build_covmatrix.C(9)'  &
    # sleep 60
    ## bgsys
    root -b -q LoadClasses.C -e '.L build_covmatrix.C+' 'run_build_covmatrix.C(21)' &
    wait
}

# ------------------------------------ Fit! ------------------------------------
shapeFit() {
    cd $BASE/ShapeFit
    root -b -q LoadClasses.C fit_shape_2d_P17B.C+
}

all() {
    genReactor
    genToyConf
    genToys
    genEvisEnu
    genSuperHists
    genPredIBD
    genCovMat
    shapeFit

}

eval $step
