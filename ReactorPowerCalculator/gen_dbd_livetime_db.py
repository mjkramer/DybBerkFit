#!/usr/bin/env python3

# This is less precise than gen_dbd_livetime_files.py

import argparse
import datetime as D
import gzip
import os
import sys

import pandas as pd

DAY_ZERO = D.datetime(2011, 12, 24)


def gzopen(path):
    return gzip.open(path) if path.endswith(".gz") else open(path)

def parse_path(path):
    "Returns (runno, fileno, site)"
    base = path.strip().split("/")[-1]
    parts = base.split(".")
    return int(parts[2]), int(parts[6][1:]), int(parts[4][2])

def read_grl(grlpath):
    return pd.DataFrame((parse_path(line) for line in gzopen(grlpath)),
                        columns=["runno", "fileno", "site"]) \
             .set_index(["runno", "fileno"])

def read_db(firstrun, lastrun):
    passwd = os.environ["DBPASS"]
    con = f"mysql://dayabay:{passwd}@dybdb1.ihep.ac.cn/offline_db"
    query = f"""SELECT runno, fileno, timestart, timeend FROM DaqRawDataFileInfo
                NATURAL JOIN DaqRawDataFileInfoVld
                WHERE stream LIKE 'EH_-Merged' AND streamtype = 'Physics'
                AND runno between {firstrun} and {lastrun}"""
    df = pd.read_sql(query, con).set_index(['runno', 'fileno'])
    df["day"] = (df.timeend - DAY_ZERO).map(lambda _: _.days) \
                                       .fillna(0).astype(int)  # for null stamps
    df["len_s"] = (df.timeend - df.timestart).map(lambda _: _.seconds)
    return df

def get_data(grlpath):
    grl = read_grl(grlpath)
    firstrun, _ = grl.index.min()
    lastrun, _ = grl.index.max()
    db = read_db(firstrun, lastrun)
    return db.join(grl, how="right")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("listfile",
                    help="List of input DAQ files (e.g. for stage1)")
    args = ap.parse_args()

    if not os.environ["DBPASS"]:
        print("You must set DBPASS!")
        return

    df = get_data(args.listfile)
    out_df = df.groupby(['day', 'site'])["len_s"].sum()
    out_df.to_csv(sys.stdout, sep="\t", header=False, float_format="%.3f")

if __name__ == '__main__':
    main()
