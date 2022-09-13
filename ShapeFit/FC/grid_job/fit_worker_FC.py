#!/usr/bin/env python3

import argparse
import functools
import os

from listfile import LockfileListReader, LockfileListWriter


@functools.lru_cache
def get_outdir():
    path = f'{os.getenv("LBNL_FIT_OUTDIR")}/FC_fits'
    os.system(f'mkdir -p {path}')
    return path


def get_outpath(toypath: str):
    toyfile = os.path.basename(toypath)
    x = toyfile.rfind('s2t13_')
    toyconfig = toyfile[len('toySpectra_'):x-1]
    params_and_ext = toyfile[x:]
    return f'{get_outdir()}/FC_{toyconfig}_{params_and_ext}'


def process(toypath: str):
    outpath = get_outpath(toypath)
    # Need to restore cwd for LockfileListReader/Writer
    oldcwd = os.getcwd()
    os.chdir(f'{os.getenv("LBNL_FIT_HOME")}/ShapeFit')
    cmd = f'root -b -q LoadClasses.C FC/fit_shape_3d_FC.C+("{toypath}", "{outpath}")'
    os.system(cmd)
    os.chdir(oldcwd)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('listdir', help='Location of input.list etc.')
    args = ap.parse_args()

    reader = LockfileListReader(f'{args.listdir}/input.list')
    writer = LockfileListWriter(f'{args.listdir}/input.list.done')

    with writer:
        for line in reader:
            toypath = line.strip()
            process(toypath)
            writer.log(line)


if __name__ == '__main__':
    # unbuf_stdout()
    main()
