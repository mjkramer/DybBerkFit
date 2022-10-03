#!/usr/bin/env bash
#SBATCH -N 1 -C haswell -L project
#SBATCH -t 06:00:00 -A dayabay -q regular

listdir=$1; shift

echo THREADS $OMP_NUM_THREADS

name=$(basename $LBNL_FIT_OUTDIR)
logdir=$SCRATCH/logs/FC_fits/$name/$SLURM_JOBID
mkdir -p $logdir

srun -o "$logdir"/output-%j.%t.txt --cpus-per-task $OMP_NUM_THREADS --kill-on-bad-exit=0 --no-kill -- ./fit_worker_FC.py $listdir
