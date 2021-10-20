#!/usr/bin/env python3

import argparse
import os

import pandas as pd

EMIN = 1.85
EMAX = 12.75
STRIDE = 5                      # Spacing of the samples we take from the model

SCALE = 1e-18
NOM_POWER_GW = 2.9
MEV_PER_GW_S = 6.2415e21

STAGE_NDET = [6, 8, 7]
STAGE_START_WEEK = [0, 42, 264]
STAGE_END_WEEK = [33, 261, 297]

DEFAULT_LIVETIME_FILE = "dbd_livetime_P17B.txt"
DEFAULT_POWER_FILE = "WeeklyAvg/WeeklyAvg_P17B_by_Beda.txt"
DEFAULT_ISOTOPE_FILE = "fissionIsotopeTable_v1.txt"
DEFAULT_MODEL_OPTION = 1

ISOTOPES = ["U235", "U238", "Pu239", "Pu241"]


# takes the mean across the three halls
def weekly_livetime(path=DEFAULT_LIVETIME_FILE):
    dbd = pd.read_csv(path, sep=r"\s+",
                      names=["day", "hall", "livetime_s"])
    dbd["week"] = dbd["day"] // 7
    return dbd.groupby("day").mean() \
                             .groupby("week").sum()["livetime_s"]


def weekly_power(path=DEFAULT_POWER_FILE):
    return pd.read_csv(path, sep=r"\s+",
                       names=["week", "core", "start", "end", "frac_pow",
                              "DUMMY",
                              "frac_U235", "frac_U238",
                              "frac_Pu239", "frac_Pu241"]) \
             .drop(columns=["start", "end", "DUMMY"]) \
             .set_index(["week", "core"])


def weekly_data(livetime_file=DEFAULT_LIVETIME_FILE,
                power_file=DEFAULT_POWER_FILE):
    lt = weekly_livetime(livetime_file)
    power = weekly_power(power_file)
    return power.join(lt)


def e_per_fis(path=DEFAULT_ISOTOPE_FILE):
    table = pd.read_csv(path, sep=r"\s+", comment="#",
                        names=["isotope", "EperFis", "NomFisFrac"]) \
              .set_index("isotope")["EperFis"].to_dict()
    return {"U235": table[1], "U238": table[2],
            "Pu239": table[3], "Pu241": table[4]}


def model_spec_path(option=DEFAULT_MODEL_OPTION):
    return {
        0: "BCW/fissionIsotopeSpectra_Huber_Linear.txt",
        1: "BCW/fissionIsotopeSpectra_Huber_Fit.txt",
        2: "BCW/fissionIsotopeSpectra_ILL_Muller_Linear.txt",
        3: "BCW/fissionIsotopeSpectra_ILL_Muller_Fit.txt",
        4: "LBNL_Spectra/fissionIsotopeSpectra_Huber_v0.txt"
    }[option]


def read_model_spec(option=DEFAULT_MODEL_OPTION):
    "nuebar/MeV/fission"
    path = model_spec_path(option)
    data = pd.read_csv(path, sep=r"\s+", comment="#",
                       names=["Enu", "spec_U235", "spec_U238",
                              "spec_Pu239", "spec_Pu241"]) \
             .set_index("Enu")
    data = data.query(f"{EMIN} <= Enu <= {EMAX}")
    return data.sort_index().iloc[::STRIDE]

def cross(df1, df2):
    return df1.reset_index().assign(dummy=0) \
                            .merge(df2.reset_index().assign(dummy=0),
                                   on="dummy") \
                            .drop(columns="dummy")

def full_data(livetime_file=DEFAULT_LIVETIME_FILE,
              power_file=DEFAULT_POWER_FILE,
              isotope_file=DEFAULT_ISOTOPE_FILE,
              model_option=DEFAULT_MODEL_OPTION):

    weekly = weekly_data(livetime_file, power_file)
    modelspec = read_model_spec(model_option)
    df = cross(weekly, modelspec).set_index(["week", "core", "Enu"])

    fisE = e_per_fis(isotope_file)

    # fission energy for each week/core
    df["MeVPerFis"] = sum(df[f"frac_{iso}"] * fisE[iso] for iso in ISOTOPES)
    # number of fissions for each week/core
    df["fissions"] = df["frac_pow"] * NOM_POWER_GW * MEV_PER_GW_S \
        * df["livetime_s"] / df["MeVPerFis"]

    for iso in ISOTOPES:
        # number of nuebar (per MeV) for each week/core/isotope, as a fn of Enu
        df[f"nNu_{iso}"] = df[f"frac_{iso}"] * df["fissions"] * df[f"spec_{iso}"]

    return df.dropna()


def dump_spectra(full_df, istage, outdir):
    ndet = STAGE_NDET[istage]
    start = STAGE_START_WEEK[istage]
    end = STAGE_END_WEEK[istage]

    df = full_df.loc[start:end].copy()
    gb = df.groupby(["core", "Enu"])
    # now we are summing over weeks, indexing over (core, Enu)
    livetime_s = gb["livetime_s"].sum()

    for iso in ISOTOPES:
        spec = gb[f"nNu_{iso}"].sum() / livetime_s * SCALE
        table = spec.reset_index().pivot(columns="core", index="Enu") \
                                  .reset_index()
        # HACK since we want only two decimals here, six for other cols
        table["Enu"] = table["Enu"].map(lambda x: "%.2f" % x)

        path = f"{outdir}/reactor_{ndet}AD_{iso}.txt"
        table.to_csv(path, sep="\t", header=None, index=False,
                     float_format="%.6f")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="isotope_spectra")
    ap.add_argument("--livetime-file", default=DEFAULT_LIVETIME_FILE)
    ap.add_argument("--power-file", default=DEFAULT_POWER_FILE)
    ap.add_argument("--isotope-file", default=DEFAULT_ISOTOPE_FILE)
    ap.add_argument("--model-option", type=int, default=DEFAULT_MODEL_OPTION)
    args = ap.parse_args()

    os.system(f"mkdir -p {args.outdir}")

    df = full_data(args.livetime_file, args.power_file,
                   args.isotope_file, args.model_option)

    for istage in range(3):
        dump_spectra(df, istage, args.outdir)


if __name__ == '__main__':
    main()