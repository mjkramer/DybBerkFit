#!/usr/bin/env bash
#SBATCH -N 1 -C haswell -L project
#SBATCH -t 06:00:00 -A dayabay -q regular

echo THREADS $OMP_NUM_THREADS

name=$(basename $LBNL_FIT_OUTDIR)
logdir=$SCRATCH/logs/FC_parscans/$name/$SLURM_JOBID
mkdir -p $logdir

mkdir -p $LBNL_FIT_OUTDIR/toys_parscans

srun -o "$logdir"/output-%j.%t.txt --cpus-per-task $OMP_NUM_THREADS --kill-on-bad-exit=0 --no-kill -- ./parscans_worker_FC.py $LBNL_FIT_OUTDIR/parscans_input.txt
