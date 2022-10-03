#!/usr/bin/env bash

name=$(basename $LBNL_FIT_OUTDIR)
logdir=$SCRATCH/logs/FC_parscans/$name
mkdir -p $logdir

sbatch -o $logdir/slurm-%j.txt $@ parscans_job_FC.sh
