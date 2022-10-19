#!/usr/bin/env bash
#SBATCH -N 1 -C haswell -L project
#SBATCH -t 03:00:00 -A dayabay -q regular

# name=$(basename $LBNL_FIT_OUTDIR)
# logdir=$SCRATCH/logs/FC_fits_band/$name
# while true; do
#     top -u $UID -n1 -b
#     echo
#     sleep 60
# done > $logdir/monitor-${SLURM_JOBID}.txt &

name=$(basename $LBNL_FIT_OUTDIR)
logdir=$SCRATCH/logs/FC_fits_band/$name

srun -o $logdir/output-%j.%t.txt --ntasks=$SLURM_NNODES --ntasks-per-node=1 --kill-on-bad-exit=0 --no-kill fit_task_FC_band.sh

# cd $LBNL_FIT_HOME/ShapeFit

# toyfile=toys_parscans/toySpectra_allsys_w_dm2ee_and_stat_s2t13_0.0850_dm2ee_0.00250_s2t14_0.0000_dm214_0.00000.root
# # root -b -q LoadClasses.C "fit_shape_3d.C+$DBG(false, \"$toyfile\", 500)"
# root -b -q LoadClasses.C "fit_shape_3d_quick.C+$DBG(false, \"$toyfile\", 500)"
