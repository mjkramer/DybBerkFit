#!/bin/bash

# NOTE: First run init.sh and prepare input data (see README.org)

# To log timing data:
# (time scripts/test_chain.sh) 2>&1 | tee output.time.log

step=$1; shift
step=${step:-all}

if [ -z $LBNL_FIT_OUTDIR ]; then
    export LBNL_FIT_OUTDIR=$LBNL_FIT_INDIR
fi

# RECOMPILE=1

BASE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"/..
DBG=${LBNL_FIT_DEBUG:+g}

ROOT=$(which root)
root() {
    test $RECOMPILE && $BASE/toySpectra/clean.sh
    time $ROOT "$@"
}

set_threads() {
    local threads=$1
    if [ -n "$LBNL_FIT_MAXTHREADS" ]; then
        if ((threads > LBNL_FIT_MAXTHREADS)); then
            threads=$LBNL_FIT_MAXTHREADS
        fi
    fi
    export OMP_NUM_THREADS=$threads
}

# Compile stuff in advance to avoid race conditions when we parallelize
time $BASE/scripts/compile.sh

[ -n "$LBNL_FIT_OUTDIR" ] && mkdir -p "$LBNL_FIT_OUTDIR"

# echo "Using IHEP fast-n spectrum (see Config.h)"

# --------------------------- Generate ToyMC samples ---------------------------
genToys() {
    set_threads 10
    cd $BASE/toySpectra
    root -b -q LoadClasses.C -e ".L genToySpectraTree.C+$DBG" "rungenToySpectraTree.C(\"sigsys\")" &
    root -b -q LoadClasses.C -e ".L genToySpectraTree.C+$DBG" "rungenToySpectraTree.C(\"bgsys\")" &
    root -b -q LoadClasses.C -e ".L genToySpectraTree.C+$DBG" "rungenToySpectraTree.C(\"dm2ee\")" &
    wait
}

# Generate Asimov toys (for sterile analysis)
# Run this in an interactive salloc session
genToys_parscans() {
    unset OMP_NUM_THREADS
    root -b -q LoadClasses.C -e ".L genToySpectraTree_parscans.C+$DBG" rungenToySpectraTree_parscans.C
}

# -------------------------- Generate evis/enu matrix --------------------------
genEvisEnu() {
    set_threads 30
    cd $BASE/toySpectra
    root -b -q LoadClasses.C genEvisToEnuMatrix.C+$DBG
    cd ../ShapeFit
    root -b -q LoadClasses.C make_evis_to_enu_matrix_fine_P17B.C+$DBG
}

# ------------------------- Generate super histograms --------------------------
# The super histograms are the nominal xsec-weighted spectra produced by each
# core in each stage. The units are essentially arbitrary, since ratios are
# always taken when using these to extrapolate from near to far. Although the
# file contains a separate histogram for each AD, there is actually no AD
# dependence.
genSuperHists() {
    cd $BASE/toySpectra
    root -b -q LoadClasses.C genSuperHistograms.C+$DBG
}

# --------------------------- Generate PredictedIBD ----------------------------
# The PredictedIBD file which contains the background-free no-oscillation IBD
# spectra of each detector. As far as I can tell, this is only used in order to
# calculate a `summed' covariance matrix in which the matrices of the three
# stages (6, 8, 7AD) are combined, with the weighting determined by the
# PredictedIBD counts. In turn, the summed matrix just seems to be a diagnostic
# that is not actually used during the fit.-- mk
genPredIBD() {
    cd $BASE/toySpectra
    root -b -q LoadClasses.C genPredictedIBD.C+$DBG
}

# Generate OscProbTables for faster near->far extrapolation
genOscProb() {
    cd $BASE/toySpectra
    root -b -q LoadClasses.C genOscProbTable.C+$DBG
}

# ------------------------ Generate covariance matrices ------------------------
genCovMat() {
    set_threads 8
    cd $BASE/ShapeFit
    root -b -q LoadClasses.C -e ".L build_covmatrix.C+$DBG" "run_build_covmatrix.C(\"sigsys\", 0)" &
    root -b -q LoadClasses.C -e ".L build_covmatrix.C+$DBG" "run_build_covmatrix.C(\"bgsys\", 1)" &
    root -b -q LoadClasses.C -e ".L build_covmatrix.C+$DBG" "run_build_covmatrix.C(\"dm2ee\", 0)" &
    wait
}

# ------------------------------------ Fit! ------------------------------------
shapeFit() {
    local period=${1:--1}       # default = -1 (6+8+7 AD)
    set_threads 12
    cd $BASE/ShapeFit
    root -b -q LoadClasses.C "fit_shape_sterile.C+$DBG(${period})"
}

# Generate Asimove delta-chi2 (4nu vs 3nu)
genAsimovDChi2() {
    set_threads 12
    cd $BASE/ShapeFit
    root -b -q LoadClasses.C fit_shape_3d_CLs.C+$DBG run_fit_shape_3d_CLs_3n.C &
    root -b -q LoadClasses.C fit_shape_3d_CLs.C+$DBG run_fit_shape_3d_CLs_4n.C &
    wait
}

# Generate contours
contours() {
    cd $BASE/ShapeFit
    root -b -q LoadClasses.C "make_data_contours_cls.C+$DBG(${period})"
    # root -b -q LoadClasses.C "make_data_contours_comparison.C+$DBG(${period})"
}

all() {
    genToys &
    genEvisEnu &
    genSuperHists &
    genPredIBD &
    wait
    genOscProb &
    genCovMat &
    wait
    shapeFit &
    genAsimovDChi2 &
    wait
    contours
}

eval $step "$@"
