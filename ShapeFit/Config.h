#pragma once

#include "TROOT.h"

namespace Config {

const int Ncores = 6;
const int Ndetectors = 8;
const int Nhalls = 3;
const int Nstage = 3;

const double lowBinInflation = 0; // was 0.1 originally

const int MaxPredictions = 16; // 4 near AD predicting 4 far AD: 4 x 4 = 16
const int detConfigEH[Ndetectors] = {1, 1, 2, 2, 3, 3, 3, 3}; // EH
const int detConfigAD[Ndetectors] = {1, 2, 1, 2, 1, 2, 3, 4}; // AD

// This part needs to be modified when the detector config changes
const int NdetectorsConfig[3][Ndetectors] = {{1, 1, 1, 0, 1, 1, 1, 0},
                                             {1, 1, 1, 1, 1, 1, 1, 1},
                                             {0, 1, 1, 1, 1, 1, 1, 1}};

// chi2 map binning
// const Int_t nsteps = 101;
// const Double_t s22t13start=0.06;
// const Double_t s22t13end=0.11;
// const Int_t nsteps_dm2 = 101;
// const Double_t dm2eestart=2.1e-3;
// const Double_t dm2eeend=2.9e-3;

const Int_t nsteps = 51;
Double_t s22t13start = 0.06;
Double_t s22t13end = 0.11;

const Int_t nsteps_dm2 = 31;
Double_t dm2eestart = 2.1e-3;
Double_t dm2eeend = 2.8e-3;

} // namespace Config
