#!/usr/bin/env bash

name=$(basename $LBNL_FIT_OUTDIR)
logdir=$SCRATCH/logs/FC_fits_band/$name
mkdir -p $logdir

sbatch -o $logdir/slurm-%j.txt $@ fit_job_FC_band.sh
