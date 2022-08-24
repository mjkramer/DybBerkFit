#include "PredSet.h"
PredSet::PredSet()
{
  for (int istage = 0; istage < Nstage; ++istage) {
    for (int idet1 = 0; idet1 < Ndetectors; ++idet1) {
      for (int idet2 = 0; idet2 < Ndetectors; ++idet2) {
        PredEvts[istage][idet1][idet2] = new TH1F();
        // HACK: Prevent TH1::Copy (in SetPred) from calling gDirectory->Delete
        // (which is behind a lock that the threads apparently can fight over)
        PredEvts[istage][idet1][idet2]->SetDirectory(nullptr);
        // HACK: Similarly prevent calling gDirectory->Append
        PredEvts[istage][idet1][idet2]->AddDirectory(false);
      }
    }
  }
  h_comb = new TH1F();
};

PredSet::~PredSet()
{
  for (int istage = 0; istage < Nstage; ++istage) {
    for (int idet1 = 0; idet1 < Ndetectors; ++idet1) {
      for (int idet2 = 0; idet2 < Ndetectors; ++idet2) {
        delete PredEvts[istage][idet1][idet2];
      }
    }
  }
  delete h_comb;
};
