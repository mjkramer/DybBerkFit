#!/bin/bash

BASE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"/..

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

genReactor &
genToyConf &
wait
