#include "genToySpectraTree.H"

#include "Binning.h"
#include "Config.h"
#include "DataSet.h"
#include "Paths.h"
#include "Spectrum.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom3.h"
#include "TTree.h"
#include <iostream>

#ifdef __CLING__
extern int omp_get_thread_num();
#else
#include <omp.h>
#endif

using namespace Config;

void genToySpectraTree(TString dataset_filename, TString output_filename,
                       double s2t13, double dm2ee, double s2t14, double dm241)
{
  int nToys = 1000;
  // string opt="_allerrors";

  // // define output histograms

  const int n_evis_bins = Binning::n_evis();
  double* evis_bins = Binning::evis();

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
  mydata->load(dataset_filename.Data());

  // Create Nominal (i.e. no random variation) Expectations
  DataSet* mydata_nominal = new DataSet();
  mydata_nominal->load(Paths::nominal_toyconfig());

  // Change oscillation parameters for testing
  if (s2t13 >= 0) {
    mydata->setDouble("sinSq2Theta13", s2t13);
    mydata_nominal->setDouble("sinSq2Theta13", s2t13);
  } else
    s2t13 = mydata_nominal->getDouble("sinSq2Theta13");

  if (dm2ee > 0) {
    mydata->setDouble("deltaMSqee", dm2ee);
    mydata_nominal->setDouble("deltaMSqee", dm2ee);
  } else
    dm2ee = mydata_nominal->getDouble("deltaMSqee");

  if (s2t14 >= 0) {
    mydata->setDouble("sinSq2Theta14", s2t14);
    mydata_nominal->setDouble("sinSq2Theta14", s2t14);
  } else
    s2t14 = mydata_nominal->getDouble("sinSq2Theta14");

  if (dm241 > 0) {
    mydata->setDouble("deltaMSq41", dm241);
    mydata_nominal->setDouble("deltaMSq41", dm241);
  } else
    dm241 = mydata_nominal->getDouble("deltaMSq41");


  const char* Theta13InputsLocation[3] = {Paths::input(0), Paths::input(1),
    Paths::input(2)};

  // Create Predictor (needed for livetimes, efficiencies, target masses... etc)
  static Predictor* myPred;
#pragma omp threadprivate(myPred)
#pragma omp parallel
  {
    myPred = new Predictor();
    for (int istage = 0; istage < Nstage; istage++) {
      myPred->LoadMainData(Theta13InputsLocation[istage]);
    }
  } // parallel

  cout << "Data loaded" << endl;

  // Create nominal Spectrum
  Spectrum* spectrumNormNominal = new Spectrum();

  cout << "Something with spectrum done" << endl;

  spectrumNormNominal->passPredictor(myPred);
  // fixme: should use distances from Predictor (from FluxCalculator) in order
  // to avoid duplication
  spectrumNormNominal->loadDistances(Paths::baselines());
  spectrumNormNominal->initialize(mydata_nominal);

  TString AccidentalSpectrumLocation[3] = {
      Paths::acc_spectra(0), Paths::acc_spectra(1), Paths::acc_spectra(2)};

  spectrumNormNominal->loadBgSpecForToy(AccidentalSpectrumLocation,
                                        Paths::li9(), Paths::amc(), Paths::fn(),
                                        Paths::aln());


  // Create Spectrum for toy
  static Spectrum* spectrumNorm;
#pragma omp threadprivate(spectrumNorm)
#pragma omp parallel
  {
    spectrumNorm = new Spectrum();
    spectrumNorm->passPredictor(myPred);
    // fixme: should use distances from Predictor (from FluxCalculator) in order
    // to avoid duplication
    spectrumNorm->loadDistances(Paths::baselines());
    spectrumNorm->initialize(mydata);
    spectrumNorm->loadBgSpecForToy(AccidentalSpectrumLocation,
                                   Paths::li9(), Paths::amc(), Paths::fn(),
                                   Paths::aln());
    // NB: Adding 1 is important because 0 means random (UUID) seed
    spectrumNorm->setRandomSeed(1 + omp_get_thread_num());
  } // parallel


  // Prepare destination file
  TFile* outfile = new TFile(output_filename.Data(), "RECREATE");

  TTree* tree_true_pars = new TTree("true_pars", "true oscillation parameters");
  tree_true_pars->Branch("sin22theta13", &s2t13, "sin22theta13/D");
  tree_true_pars->Branch("deltamsqee", &dm2ee, "deltamsqee/D");
  tree_true_pars->Branch("sin22theta14", &s2t14, "sin22theta14/D");
  tree_true_pars->Branch("deltamsq41", &dm241, "deltamsq41/D");
  tree_true_pars->Fill();

  TTree* tree = new TTree("tr", "toy spectra");

  for (int istage = 0; istage < Nstage; ++istage) {
    for (int idet = 0; idet < Ndetectors; ++idet) {
      sprintf(name, "h_stage%i_ad%i", istage + 1, idet + 1);
      tree->Branch(name, "TH1F", h_ad[istage][idet]);
    }
  }


  // Generate nominal spectrum
  // cout << "-------------------- Nominal Spectrum -------------------" <<
  // endl;
  spectrumNormNominal->updateAntinu();
  spectrumNormNominal->updateBgDetected();

  for (int istage = 0; istage < Nstage; ++istage) {
    for (int idet = 0; idet < Ndetectors; ++idet) {
      h_nominal_ad[istage][idet]->Reset();
      for (int ibin = 0; ibin < spectrumNormNominal->nSamples(); ++ibin) {
        h_nominal_ad[istage][idet]->Fill(
            spectrumNormNominal->energyArray(idet)[ibin],
            spectrumNormNominal->positronDetectedArray(istage, idet)[ibin] *
                spectrumNormNominal->binWidth());
        if (istage == 2)
          cout << "The added positron at energy "
               << spectrumNormNominal->energyArray(idet)[ibin] << " is "
               << spectrumNormNominal->positronDetectedArray(istage, idet)[ibin]
               << endl; // tmp
      }                 // ibin loop
      for (int ibin = 0; ibin < spectrumNormNominal->nSamplesBkg(); ++ibin) {
        h_nominal_ad[istage][idet]->Fill(
            spectrumNormNominal->energyArrayBkg(idet)[ibin],
            spectrumNormNominal->bgDetectedArray(istage, idet)[ibin]);
        if (istage == 2)
          cout << "The added bg at energy "
               << spectrumNormNominal->energyArray(idet)[ibin] << " is "
               << spectrumNormNominal->bgDetectedArray(istage, idet)[ibin]
               << endl; // tmp
      }


      // cout << "AD" << idet+1 << ": " << h_nominal_ad[idet]->Integral() <<
      // endl;//tmp
    }
  }

  TRandom3* generator = new TRandom3();

  // Generate toys
#pragma omp parallel for
  for (int itoy = 0; itoy < nToys; itoy++) {
    cout << "-------------------- Generating toy " << itoy
         << "-------------------" << endl;
    spectrumNorm->updateAntinu();
    spectrumNorm->updateBgDetected();

    // this part is added to inflate low energy bin uncertainty
    double rand_bin[8][3] = {}; // 8 detectors and 3 bins
    double uncertainty = Config::lowBinInflation;
    if (uncertainty) {
      for (int idet = 0; idet < Ndetectors; ++idet) {
        for (int ibin = 0; ibin < 3; ++ibin) {
          rand_bin[idet][ibin] = generator->Gaus(0., uncertainty);
          if (rand_bin[idet][ibin] < -1.)
            rand_bin[idet][ibin] = -1.;
          cout << "Rand is " << rand_bin[idet][ibin] << endl;
        }
      }
    }

#pragma omp critical
    {
      for (int istage = 0; istage < Nstage; ++istage) {
        for (int idet = 0; idet < Ndetectors; ++idet) {
          h_ad[istage][idet]->Reset();
          for (int ibin = 0; ibin < spectrumNorm->nSamples(); ++ibin) {
            h_ad[istage][idet]->Fill(
                spectrumNorm->energyArray(idet)[ibin],
                spectrumNorm->positronDetectedArray(istage, idet)[ibin] *
                    spectrumNorm->binWidth());
            // cout << "The added positron at energy " <<
            // spectrumNorm->energyArray(idet)[ibin] << " is " <<
            // spectrumNorm->positronDetectedArray(idet)[ibin] << endl;//tmp

          } // ibin loop
          for (int ibin = 0; ibin < spectrumNorm->nSamplesBkg(); ++ibin) {
            h_ad[istage][idet]->Fill(
                spectrumNorm->energyArrayBkg(idet)[ibin],
                spectrumNorm->bgDetectedArray(istage, idet)[ibin]);
            // cout << "The added bg at energy " <<
            // spectrumNorm->energyArray(idet)[ibin] << " is " <<
            // spectrumNorm->bgDetectedArray(idet)[ibin] << endl;//tmp

          } // ibin loop

          for (int ibin = 0; ibin < 3; ++ibin) {
            float bincontent = h_ad[istage][idet]->GetBinContent(ibin + 1);
            h_ad[istage][idet]->SetBinContent(
                ibin + 1, bincontent * (1. + rand_bin[idet][ibin]));
          }

          // cout << "AD" << idet+1 << ": " << h_ad[idet]->Integral() <<
          // endl;//tmp

        } // idet loop
      }   // istage loop
      tree->Fill();
    } // critical

  } // itoy loop

  outfile->cd();

  for (int istage = 0; istage < Nstage; ++istage) {
    for (int idet = 0; idet < Ndetectors; ++idet) {
      h_nominal_ad[istage][idet]->Write();
    }
  }

  outfile->Write();
  outfile->Close();

  gDirectory->Delete("*");

} // end of genToySpectraTree
