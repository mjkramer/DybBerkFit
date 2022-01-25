#include "Binning.h"
#include "Config.h"
#include "DataSet.h"
#include "Paths.h"
#include "Spectrum.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"

#include <iostream>

using namespace Config;

// note that mixing angles here are ZERO, not nominal
void genPredictedIBD(Double_t s2t13 = 0, Double_t dm2ee = -1,
                     const char* dataset_filename = Paths::nominal_toyconfig(),
                     Double_t s2t14 = 0, Double_t dm241 = -1)
{
  // XXX fine binning, why?
  const int n_evis_bins = Binning::n_evis_fine();
  const double* evis_bins = Binning::evis_fine();

  Char_t name[1024];
  TH1F* h_nominal_ad[Nstage][Ndetectors];
  TH1F* h_ad[Nstage][Ndetectors];
  for (int istage = 0; istage < Nstage; ++istage) {
    for (int idet = 0; idet < Ndetectors; ++idet) {
      cout << "Creating histogram for Stage# " << istage + 1 << ": AD# "
           << idet + 1 << endl;
      sprintf(name, "h_stage%i_ad%i", istage + 1, idet + 1);
      h_ad[istage][idet] = new TH1F(name, name, n_evis_bins, evis_bins);
      sprintf(name, "h_nominal_stage%i_ad%i", istage + 1, idet + 1);
      h_nominal_ad[istage][idet] = new TH1F(name, name, n_evis_bins, evis_bins);
    }
  }

  // Create Expectations
  DataSet* mydata = new DataSet();
  mydata->load(dataset_filename);

  // Create Nominal (i.e. no random variation) Expectations
  DataSet* mydata_nominal = new DataSet();
  mydata_nominal->load(Paths::nominal_toyconfig());

  // Change oscillation parameters for testing
  if (s2t13 >= 0) {
    cout << "Changing theta" << endl;
    mydata->setDouble("sinSq2Theta13", s2t13);
    mydata_nominal->setDouble("sinSq2Theta13", s2t13);
  }

  if (dm2ee > 0) {
    mydata->setDouble("deltaMSqee", dm2ee);
    mydata_nominal->setDouble("deltaMSqee", dm2ee);
  }

  if (s2t14 >= 0) {
    mydata->setDouble("sinSq2Theta14", s2t14);
    mydata_nominal->setDouble("sinSq2Theta14", s2t14);
  }

  if (dm241 > 0) {
    mydata->setDouble("deltaMSq41", dm241);
    mydata_nominal->setDouble("deltaMSq41", dm241);
  }


  // Create Predictor (needed for livetimes, efficiencies, target masses... etc)
  Predictor* myPred = new Predictor();

  const char* Theta13InputsLocation[3] = {Paths::input(0), Paths::input(1),
    Paths::input(2)};

  cout << "There are " << Nstage << " stages" << endl;

  for (int istage = 0; istage < Nstage; istage++) {
    myPred->LoadMainData(Theta13InputsLocation[istage]);
  }


  // Create nominal Spectrum
  Spectrum* spectrumNormNominal = new Spectrum();
  spectrumNormNominal->passPredictor(myPred);
  // fixme: should use distances from Predictor (from FluxCalculator) in order
  // to avoid duplication
  spectrumNormNominal->loadDistances(Paths::baselines());
  spectrumNormNominal->initialize(mydata_nominal);

  spectrumNormNominal->loadBgSpecForToy();

  // Prepare destination file
  TFile* outfile = new TFile(Paths::predicted_ibd(), "RECREATE");

  // Generate nominal spectrum
  // cout << "-------------------- Nominal Spectrum -------------------" <<
  // endl;
  spectrumNormNominal->updateAntinu();

  for (int istage = 0; istage < Nstage; ++istage) {
    for (int idet = 0; idet < Ndetectors; ++idet) {
      h_nominal_ad[istage][idet]->Reset();
      for (int ibin = 0; ibin < spectrumNormNominal->nSamples(); ++ibin) {
        // Don't fill below 0.7 MeV
        if (spectrumNormNominal->energyArray(idet)[ibin] < 0.7)
          continue;

        h_nominal_ad[istage][idet]->Fill(
            spectrumNormNominal->energyArray(idet)[ibin],
            spectrumNormNominal->positronDetectedArray(istage, idet)[ibin] *
                spectrumNormNominal->binWidth());
        // cout << "The added positron at energy " <<
        // spectrumNormNominal->energyArray(idet)[ibin] << " is " <<
        // spectrumNormNominal->positronDetectedArray(idet)[ibin] << endl;//tmp
      } // ibin loop

      // cout << "AD" << idet+1 << ": " << h_nominal_ad[idet]->Integral() <<
      // endl;//tmp
    }
  }

  outfile->cd();

  for (int istage = 0; istage < Nstage; ++istage) {
    for (int idet = 0; idet < Ndetectors; ++idet) {
      h_nominal_ad[istage][idet]->Write();
    }
  }

  outfile->Write();
  outfile->Close();

  gDirectory->Delete("*");

} // end of genPredictedIBD
