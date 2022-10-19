#!/usr/bin/env bash

cd $LBNL_FIT_HOME/ShapeFit

toyfile=toys_parscans/toySpectra_allsys_w_dm2ee_and_stat_s2t13_0.0850_dm2ee_0.00250_s2t14_0.0000_dm214_0.00000.root
# root -b -q LoadClasses.C "fit_shape_3d.C+$DBG(false, \"$toyfile\", 500)"
root -b -q LoadClasses.C "fit_shape_3d_quick.C+$DBG(false, \"$toyfile\", 500, $SLURM_NODEID, $SLURM_NNODES)"
