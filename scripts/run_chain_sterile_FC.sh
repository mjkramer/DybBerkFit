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
    # set_threads 10
    set_threads $(( $(nproc) / 3 ))
    cd $BASE/toySpectra
    root -b -q LoadClasses.C -e ".L genToySpectraTree.C+$DBG" "rungenToySpectraTree.C(\"sigsys\")" &
    root -b -q LoadClasses.C -e ".L genToySpectraTree.C+$DBG" "rungenToySpectraTree.C(\"bgsys\")" &
    root -b -q LoadClasses.C -e ".L genToySpectraTree.C+$DBG" "rungenToySpectraTree.C(\"dm2ee\")" &
    wait
    unset OMP_NUM_THREADS
}

# Generate toys with syst/stat variations for getting Dchi2 distributions
# run in salloc
genToys_parscans() {
    # # unset OMP_NUM_THREADS
    # export OMP_NUM_THREADS=4
    # cd $BASE/toySpectra
    # srun -n 16 root -b -q LoadClasses.C -e ".L genToySpectraTree.C+$DBG" FC/rungenToySpectraTree_parscans_FC.C+
    # unset OMP_NUM_THREADS
    ntasks=16
    set_threads $(($(nproc) / $ntasks))
    cd $BASE/ShapeFit/FC/grid_job
    ./prep_parscans_job_FC.sh
    ./submit_parscans_job_FC.sh -N 10 -t 4:00:00 --ntasks-per-node $ntasks
    unset OMP_NUM_THREADS
}

# -------------------------- Generate evis/enu matrix --------------------------
genEvisEnu() {
    # set_threads 30
    # leave 2 threads for super hists, predicted IBD
    set_threads $(( $(nproc) - 2 ))
    cd $BASE/toySpectra
    root -b -q LoadClasses.C genEvisToEnuMatrix.C+$DBG
    cd ../ShapeFit
    root -b -q LoadClasses.C make_evis_to_enu_matrix_fine_P17B.C+$DBG
    unset OMP_NUM_THREADS
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
    unset OMP_NUM_THREADS
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
    unset OMP_NUM_THREADS
}

# Generate OscProbTables for faster near->far extrapolation
# TODO: Multithread?
genOscProb() {
    # cd $BASE/ShapeFit/FC/grid_job
    # ./prep_fit_job_FC.sh
    cd $BASE/ShapeFit
    root -b -q LoadClasses.C genOscProbTable.C+$DBG
    unset OMP_NUM_THREADS
}

# ------------------------ Generate covariance matrices ------------------------
genCovMat() {
    # due to rounding down, should have at least 1 thread available for osc prob table
    # set_threads $(( $(nproc) / 3 ))
    set_threads 16
    cd $BASE/ShapeFit
    root -b -q LoadClasses.C -e ".L build_covmatrix.C+$DBG" "run_build_covmatrix.C(\"sigsys\", 0)" &
    root -b -q LoadClasses.C -e ".L build_covmatrix.C+$DBG" "run_build_covmatrix.C(\"bgsys\", 1)" &
    root -b -q LoadClasses.C -e ".L build_covmatrix.C+$DBG" "run_build_covmatrix.C(\"dm2ee\", 0)" &
    wait
    unset OMP_NUM_THREADS
}

# ------------------------------------ Fit! ------------------------------------
# run in salloc
shapeFit() {
    # local period=${1:--1}       # default = -1 (6+8+7 AD)
    set_threads 12
    cd $BASE/ShapeFit
    # root -b -q LoadClasses.C "fit_shape_3d.C+$DBG(${period})"
    root -b -q LoadClasses.C "fit_shape_3d.C+$DBG"
    # nominal:
    root -b -q LoadClasses.C "fit_shape_3d.C+$DBG(true)"
    unset OMP_NUM_THREADS
}

# Generate delta-chi2 (point - bestfit)
genDChi2() {
    ntasks=16
    set_threads $(($(nproc) / $ntasks))
    cd $BASE/ShapeFit/FC/grid_job
    ./prep_fit_job_FC.sh
    ./submit_fit_job_FC.sh -N 10 -t 4:00:00 --ntasks-per-node $ntasks
    unset OMP_NUM_THREADS
}

# Generate contours
contours() {
    cd $BASE/ShapeFit
    root -b -q LoadClasses.C "FC/genContours_FC.C+$DBG"
    root -b -q LoadClasses.C "FC/make_data_contours_FC.C+$DBG"
    # expected (nominal):
    # root -b -q LoadClasses.C "FC/make_data_contours_FC.C+$DBG(true)"
    # root -b -q LoadClasses.C "FC/make_mc_contours_FC.C+$DBG(true)"

    # root -b -q LoadClasses.C "make_data_contours_comparison.C+$DBG"
    unset OMP_NUM_THREADS
}

# NOTE: all() doesn't actually work on login nodes since we use srun in genToys_parscans
all() {
    # genToys

    # genToys_parscans

    genEvisEnu &
    genSuperHists &
    genPredIBD &
    wait

    genOscProb &
    genCovMat &
    wait

    # shapeFit

    # genDChi2

    # contours
}

eval $step "$@"
