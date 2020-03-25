#!/bin/bash

# NOTE: Run prep_dirs.sh and install_example.sh first

RECOMPILE=0

BASE=$(pwd)

ROOT=$(which root)
root() {
    test $RECOMPILE && $BASE/toySpectra/clean.sh
    $ROOT "$@"
}

# Compile stuff in advance to avoid race conditions when we parallelize
precompile() {
    dummyScript=$(mktemp --suffix=.C)

    pushd ShapeFit
    cat > $dummyScript <<EOF
    {
    gROOT->ProcessLine(".L build_covmatrix.C+");
    }
EOF
    root -b -q LoadClasses.C $dummyScript
    popd

    pushd toySpectra
    cat > $dummyScript <<EOF
    {
    gROOT->ProcessLine(".L genToySpectraTree.C+");
    }
EOF
    root -b -q LoadClasses.C $dummyScript
    popd

    rm $dummyScript
}

# precompile

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
    root -b -q 'rungenToySpectraTree.C(2)' #& # sigsys
    # sleep 60
    root -b -q 'rungenToySpectraTree.C(3)' #& # bgsys
    # wait
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
    root -b -q 'run_build_covmatrix.C(9)'  #& # sigsys
    # sleep 60
    root -b -q 'run_build_covmatrix.C(21)' #& # bgsys
    # wait
}

# ------------------------------------ Fit! ------------------------------------
shapeFit() {
    cd $BASE/ShapeFit
    root -b -q LoadClasses.C fit_shape_2d_P17B.C+

}

# The following only must be run once
genReactor
genToyConf

# The following must(?) be repeated for each selection
genToys
genEvisEnu
genSuperHists
genPredIBD
genCovMat
shapeFit
