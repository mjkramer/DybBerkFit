#!/bin/bash

for i in `seq 0 1`; do
    root -q -b "run_build_covmatrix.C($i)"
done

#for i in `seq 0 18`; do
#    root -q -b "run_build_covmatrix.C($i)"
#done
