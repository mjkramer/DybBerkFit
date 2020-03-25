#!/bin/bash

while getopts "m" opt; do
    case $opt in
        m)
            cmd=mv
            ;;
        *)
            echo "Unknown option -$opt"
            ;;
    esac
done
shift $((OPTIND-1))

cmd=${cmd:-cp -R}

bakdir=$1; shift

if [ -z "$bakdir" ]; then
    echo "Specify a backup dir"
    exit 1
fi

mkdir -p $bakdir
mkdir -p $bakdir/outputs

mkdir -p $bakdir/isotope_spectra
$cmd ReactorPowerCalculator/isotope_spectra_by_Beda/reactor_*.txt $bakdir/isotope_spectra

mkdir -p $bakdir/data_file
$cmd toySpectra/data_file/dyb_data_*.txt $bakdir/data_file
cp $bakdir/data_file/dyb_data_v1_nominal.txt $bakdir/data_file/dyb_data_v1_nominal_noosc.txt toySpectra/data_file

mkdir -p $bakdir/Flux
$cmd ShapeFit/Flux/SuperHistograms_*.root $bakdir/Flux

mkdir -p $bakdir/PredictedIBD
$cmd ShapeFit/PredictedIBD/PredictedIBD_*.root $bakdir/PredictedIBD

mkdir -p $bakdir/covariance_matrices
$cmd ShapeFit/covariance_matrices/matrix_*.txt $bakdir/covariance_matrices

$cmd ShapeFit/final_covmatrix_*.txt $bakdir

mkdir -p $bakdir/fit_result_files
$cmd ShapeFit/fit_result_files/fit_shape_2d_*.root $bakdir/fit_result_files

$cmd ShapeFit/matrix_evis_to_enu*.txt $bakdir

$cmd outputs/evis_to_enu_*.root $bakdir/outputs
$cmd outputs/toySpectra_*.root $bakdir/outputs
