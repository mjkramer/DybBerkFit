#!/bin/bash

# NOTE: Run prep_dirs.sh and install_example.sh first

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

pushd ReactorPowerCalculator
root -b -q 'Produce_Isotope_SpectraP17B_unblinded.C(1)'
cd isotope_spectra_by_Beda
root -b -q make_combined_spectra_P17B_unblinded.C
popd

pushd toySpectra/data_file
./generate_data_file.py
popd

pushd toySpectra
./clean.sh
root -b -q 'rungenToySpectraTree.C(2)' #& # sigsys
# sleep 60
./clean.sh
root -b -q 'rungenToySpectraTree.C(3)' #& # bgsys
# wait
popd

pushd toySpectra
./clean.sh
root -b -q LoadClasses.C genEvisToEnuMatrix.C+
cd ../ShapeFit
root -b -q LoadClasses.C make_evis_to_enu_matrix_fine_P17B.C+
popd

pushd toySpectra
./clean.sh
root -b -q LoadClasses.C genSuperHistograms.C+ #&
# sleep 20
./clean.sh
root -b -q LoadClasses.C genPredictedIBD.C+ #&
# wait
popd

pushd ShapeFit
../toySpectra/clean.sh
root -b -q 'run_build_covmatrix.C(9)'  #& # sigsys
# sleep 60
../toySpectra/clean.sh
root -b -q 'run_build_covmatrix.C(21)' #& # bgsys
# wait
popd

pushd ShapeFit
../toySpectra/clean.sh
root -b -q LoadClasses.C fit_shape_2d_P17B.C+
popd
