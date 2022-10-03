#!/usr/bin/env bash

cd $LBNL_FIT_HOME/toySpectra

root -b -q LoadClasses.C "FC/dump_parscans_FC.C(\"$LBNL_FIT_OUTDIR/parscans_input.txt\")"
