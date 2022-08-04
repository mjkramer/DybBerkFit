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

const Int_t nsteps = 31;
Double_t s22t13start = 0.07;
Double_t s22t13end = 0.10;

const Int_t nsteps_dm2 = 31;
Double_t dm2eestart = 2.1e-3;
Double_t dm2eeend = 2.8e-3;

const double S22T13 = 0.084;
const double DM2EE = 2.49e-3;

const int nsteps_s22t14 = 31;
const double s22t14start=1.0e-3;
const double s22t14end=1.0e-0;
const double log_s22t14_step = log(s22t14end/s22t14start) / (double)(nsteps_s22t14-1);

//grid in dm241 has to be same as that of toy MC generation -----------
const int nsteps_dm214 = 120;
const int nsteps_dm214_2 = 100;
const int nsteps_dm214_all = nsteps_dm214 + nsteps_dm214_2;

const double old_log_dm214_step = log(0.05/0.0005) / 99.0; //Set log step first
const double dm214start = exp(log(0.0005) + old_log_dm214_step*(-20));
const double dm214end = 0.05;

const double dm214start_2 = 0.05 + 0.25/100;
const double dm214end_2 = 0.3;

const double log_dm214_step = log(dm214end/dm214start) / 119.0;
const double lin_dm214_step = (dm214end_2 - dm214start_2) / 99.0;

} // namespace Config
