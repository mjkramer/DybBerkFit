#!/bin/bash

BASE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"/..

# -------------------------- Generate reactor spectra --------------------------
genReactor() {
    cd $BASE/ReactorPowerCalculator
    # root -b -q "Produce_Isotope_SpectraP17B_unblinded.C(1)"
    # cd isotope_spectra_by_Beda
    if [[ -n $LBNL_FIT_BLINDED && $LBNL_FIT_BLINDED != "0" ]]; then
        tag=blinded
        power_file=WeeklyAvg_v4v5v3v1_blinded.txt
        pushd WeeklyAvg
        ./make_WeeklyAvg_v4v5v3v1_blinded.py
        popd
    else
        tag=unblinded
        power_file=weekly_power_fulldata_release.dat
    fi
    ./prod_isotope_spec.py --power-file WeeklyAvg/$power_file --outdir isotope_spectra_v4v5v3v1_$tag
    cd isotope_spectra_v4v5v3v1_$tag
    root -b -q "make_combined_spectra_v4v5v3v1_$tag.C(true, true)"
    # Uncomment below to make spectra w/o SNF and/or noneq
    # (also uncomment in toySpectra/genPredForCompare.C)
    # root -b -q "make_combined_spectra_v4v5v3v1_$tag.C(true, false)"
    # root -b -q "make_combined_spectra_v4v5v3v1_$tag.C(false, true)"
    # root -b -q "make_combined_spectra_v4v5v3v1_$tag.C(false, false)"
}

# ------------------------ Generate ToyMC config files -------------------------
genToyConf() {
    cd $BASE/toySpectra/data_file
    ./generate_data_file.py
}

genReactor &
genToyConf &
wait
