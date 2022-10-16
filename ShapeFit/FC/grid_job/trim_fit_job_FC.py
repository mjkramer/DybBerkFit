#!/usr/bin/env python3

# XXX include null (3nu) hypothesis (first line in input.list.full)?

import argparse
import os
import shutil as sh

import numpy as np

import ROOT as R
oldpwd = os.getcwd()
os.chdir('../..')
R.gInterpreter.ExecuteMacro('LoadClasses.C')
os.chdir(oldpwd)


def read_paths(inpfile):
    result = {}
    for line in open(inpfile):
        path = line.strip()
        assert path.endswith('.root')
        parts = os.path.basename(path)[:-5].split('_')
        assert parts[-4] == 's2t14' and parts[-2] == 'dm214'
        s2t14, dm214 = parts[-3], parts[-1] # note: these are strings
        result[(s2t14, dm214)] = path
    return result


def get_s2t14_steps():
    return [np.exp(np.log(R.Config.s22t14start) + R.Config.log_s22t14_step * i)
            for i in range(R.Config.nsteps_s22t14)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cls_file', help='ROOT file containing CLs contour')
    # ap.add_argument('--cls-graph', default='dm241_vs_sin22theta14_cls_95cl_0',
    #                 help='Name of the CLs contour TGraph')
    # ap.add_argument('--fraction', default=0.2,
    #                 help='Fraction of grid points to include')
    ap.add_argument('--npoints-right', default=5,
                    help='Number of points (in s2t14) to the right of the first point with CLs < 0.05 (including the point))')
    ap.add_argument('--npoints-left', default=5,
                    help='Number of points (in s2t14) to the left of the first point with CLs < 0.05 (not including the point))')
    args = ap.parse_args()

    direc = f'{os.getenv("LBNL_FIT_OUTDIR")}/FC_input'
    if not os.path.exists(f'{direc}/input.list.full'):
        sh.copy(f'{direc}/input.list', f'{direc}/input.list.full')

    pathdict = read_paths(f'{direc}/input.list.full')

    ranger_dm41 = R.Ranger()
    ranger_dm41_2 = R.Ranger()
    ranger_dm41.nsteps = R.Config.nsteps_dm214
    ranger_dm41.min = R.Config.dm214start
    ranger_dm41.max = R.Config.dm214end
    ranger_dm41.setLogScale()
    ranger_dm41_2.nsteps = R.Config.nsteps_dm214_2
    ranger_dm41_2.min = R.Config.dm214start_2
    ranger_dm41_2.max = R.Config.dm214end_2

    s2t14_steps = get_s2t14_steps()

    f = R.TFile(args.cls_file)
    h = f.Get('h_cls')
    assert h.GetNbinsX() == len(s2t14_steps)
    assert h.GetNbinsY() == R.Config.nsteps_dm214_all

    outf = open(f'{direc}/input.list', 'w')

    for dm2_bin in range(1, h.GetNbinsY()+1):
        idm2 = dm2_bin - 1
        if idm2 < R.Config.nsteps_dm214:
            dm214 = ranger_dm41.returnVal(idm2)
        else:
            dm214 = ranger_dm41_2.returnVal(idm2 - R.Config.nsteps_dm214)

        cls_vals = [h.GetBinContent(s2t_bin, dm2_bin)
                    for s2t_bin in range(1, len(s2t14_steps)+1)]
        icrit = [i for (i, x) in enumerate(cls_vals) if x < 0.05][0]
        start, end = icrit - args.npoints_left, icrit + args.npoints_right
        if end > len(s2t14_steps):
            start -= end - len(s2t14_steps)
            end = len(s2t14_steps)
        if start < 0:
            end += -start
            start = 0

        for i in range(start, end):
            try:
                path = pathdict[('%4.4f' % s2t14_steps[i], '%5.5f' % dm214)]
            except IndexError:
                print('WTF', i)
                continue
            outf.write(path + '\n')

    outf.close()


if __name__ == '__main__':
    main()
