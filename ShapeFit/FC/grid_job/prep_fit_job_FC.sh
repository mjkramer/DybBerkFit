#!/usr/bin/env bash

listdir=$LBNL_FIT_OUTDIR/FC_input
mkdir -p $listdir

if [[ -e $listdir/input.list ]]; then
    echo "Input list already exists, exiting"
    exit 1
fi

find $LBNL_FIT_OUTDIR/toys_parscans -name '*.root' | sort > $listdir/input.list
