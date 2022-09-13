#!/usr/bin/env bash

name=$(basename $LBNL_FIT_OUTDIR)
logdir=$SCRATCH/logs/FC_fits/$name
mkdir -p $logdir

listdir=$LBNL_FIT_OUTDIR/FC_input

sbatch -o $logdir/slurm-%j.txt $@ fit_job_FC.sh $listdir
