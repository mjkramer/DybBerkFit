#!/usr/bin/env python3

import argparse
import os

from listfile import LockfileListReader, LockfileListWriter

TOYCONFIG = 'allsys_w_dm2ee_and_stat'

def process(line):
    s2t13, dm2ee, s2t14, dm214 = line.strip().split()
    # Need to restore cwd for LockfileListReader/Writer
    oldcwd = os.getcwd()
    os.chdir(f'{os.getenv("LBNL_FIT_HOME")}/toySpectra')
    outdir = f'{os.getenv("LBNL_FIT_OUTDIR")}/toys_parscans'
    outpath = f'{outdir}/toySpectra_{TOYCONFIG}_s2t13_{s2t13}_dm2ee_{dm2ee}_s2t14_{s2t14}_dm214_{dm214}.root'
    toyconfigfile = f'{os.getenv("LBNL_FIT_HOME")}/toySpectra/data_file/dyb_data_v1_{TOYCONFIG}.txt'
    cmd = f'root -b -q LoadClasses.C \'genToySpectraTree.C+("{toyconfigfile}", "{outpath}", {s2t13}, {dm2ee}, {s2t14}, {dm214})\''
    os.system(cmd)
    os.chdir(oldcwd)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('listfile', help='parscans_input.txt')
    args = ap.parse_args()

    reader = LockfileListReader(f'{args.listfile}')
    writer = LockfileListWriter(f'{args.listfile}.done')

    with writer:
        for line in reader:
            process(line)
            writer.log(line)


if __name__ == '__main__':
    # unbuf_stdout()
    main()
