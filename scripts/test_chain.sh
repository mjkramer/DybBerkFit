#!/bin/bash

# NOTE: Run prep_dirs.sh and install_inputs.sh first

# To log timing data:
# (time scripts/test_chain.sh) 2>&1 | tee output.time.log

step=$1; shift
step=${step:-all}

# RECOMPILE=1

BASE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"/..
DBG=${LBNL_FIT_DEBUG:+g}

ROOT=$(which root)
root() {
    test $RECOMPILE && $BASE/toySpectra/clean.sh
    time $ROOT "$@"
}

# Compile stuff in advance to avoid race conditions when we parallelize
time $BASE/scripts/compile.sh

# echo "Using IHEP fast-n spectrum (see Config.h)"

# -------------------------- Generate reactor spectra --------------------------
genReactor() {
    cd $BASE/ReactorPowerCalculator
    root -b -q "Produce_Isotope_SpectraP17B_unblinded.C(1)"
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
    export OMP_NUM_THREADS=10
    cd $BASE/toySpectra
    ## sigsys:
    root -b -q LoadClasses.C -e ".L genToySpectraTree.C+$DBG" "rungenToySpectraTree.C(2)" &
    ## bgsys:
    root -b -q LoadClasses.C -e ".L genToySpectraTree.C+$DBG" "rungenToySpectraTree.C(3)" &
    wait
}

# -------------------------- Generate evis/enu matrix --------------------------
genEvisEnu() {
    export OMP_NUM_THREADS=30
    cd $BASE/toySpectra
    root -b -q LoadClasses.C genEvisToEnuMatrix.C+$DBG
    cd ../ShapeFit
    root -b -q LoadClasses.C make_evis_to_enu_matrix_fine_P17B.C+$DBG
}

# ------------------------- Generate super histograms --------------------------
genSuperHists() {
    cd $BASE/toySpectra
    root -b -q LoadClasses.C genSuperHistograms.C+$DBG
}

# --------------------------- Generate PredictedIBD ----------------------------
genPredIBD() {
    cd $BASE/toySpectra
    root -b -q LoadClasses.C genPredictedIBD.C+$DBG
}

# ------------------------ Generate covariance matrices ------------------------
genCovMat() {
    export OMP_NUM_THREADS=8
    cd $BASE/ShapeFit
    ## sigsys:
    root -b -q LoadClasses.C -e ".L build_covmatrix.C+$DBG" "run_build_covmatrix.C(9)" &
    ## bgsys:
    root -b -q LoadClasses.C -e ".L build_covmatrix.C+$DBG" "run_build_covmatrix.C(21)" &
    wait
}

# ------------------------------------ Fit! ------------------------------------
shapeFit() {
    export OMP_NUM_THREADS=12
    cd $BASE/ShapeFit
    root -b -q LoadClasses.C fit_shape_2d_P17B.C+$DBG
}

all() {
    genToyConf &
    genReactor &
    wait
    genToys &
    genEvisEnu &
    genSuperHists &
    genPredIBD &
    wait
    genCovMat
    shapeFit
}

eval $step
