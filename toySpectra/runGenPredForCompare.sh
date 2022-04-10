#!/usr/bin/env bash

root -b -q LoadClasses.C "genPredForCompare.C+(true, true)"

# Uncomment below to make spectra w/o SNF and/or noneq
# (also uncomment in scripts/init.sh)

# root -b -q LoadClasses.C "genPredForCompare.C+(true, false)"
# root -b -q LoadClasses.C "genPredForCompare.C+(false, true)"
# root -b -q LoadClasses.C "genPredForCompare.C+(false, false)"
