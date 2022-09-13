#!/usr/bin/env bash

#SBATCH -N 1 -C haswell -L project
#SBATCH -t 06:00:00 -A dayabay -q regular

listdir=$1; shift

srun -l -n $SLURM_NTASKS --kill-on-bad-exit=0 --no-kill -- ./fit_worker_FC.py $listdir
