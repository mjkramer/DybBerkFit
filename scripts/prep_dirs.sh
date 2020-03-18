#!/bin/bash

# Run me from ../

mkdir ShapeFit/Spectra
mkdir ShapeFit/Inputs

mkdir outputs
mkdir ShapeFit/Flux
mkdir ShapeFit/PredictedIBD
mkdir ShapeFit/covariance_matrices
mkdir ShapeFit/fit_result_files

# mkdir DayaBaySpectrum/p15a_reactor/output
# ln -s DayaBaySpectrum/p15a_reactor/output reactor_covmatrix/p15a/isotope_spectra

ln -s ../ReactorPowerCalculator/isotope_spectra_by_Beda reactor_covmatrix/p17b_unblinded

mkdir ShapeFit/Distances
ln -s ../../toySpectra/unblinded_baseline.txt ShapeFit/Distances/
