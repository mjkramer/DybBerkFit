#include "Binning.h"
#include "Config.h"
#include "DataSet.h"
#include "Paths.h"
#include "Predictor.h"

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TMinuit.h"
#include "TTree.h"

#include <fstream>
#include <iostream>
#include <map>
#include <stdio.h>
#include <stdlib.h>

using namespace Config;

// XXX replace 37 with Binning::n_evis() (which we must make constexpr?)
const int cov_matrix_dimension = 16 * Nstage * 37; // 37 evis bins

Double_t M[cov_matrix_dimension][cov_matrix_dimension];

void build_covmatrix(const Char_t* toymc_filename,
                     const Char_t* output_filename, Int_t bkg_flag)
{
  cout << "Dimension of cov matrix is " << cov_matrix_dimension << "x"
       << cov_matrix_dimension << endl;

  const Int_t n_evis_bins = Binning::n_evis();
  double* evis_bins = Binning::evis();

  const Int_t n_enu_bins = Binning::n_enu();
  double* enu_bins = Binning::enu();

  const char* Theta13InputsLocation[3] = {Paths::input(0), Paths::input(1),
    Paths::input(2)};

  Int_t ntoys = -1;

  static Predictor* myPred;
  static FluxCalculator* fluxcalc;
#pragma omp threadprivate(myPred, fluxcalc)
#pragma omp parallel
  {
    myPred = new Predictor();
    for (int istage = 0; istage < Nstage; istage++) {
      myPred->LoadMainData(Theta13InputsLocation[istage]);
    }

    // flux calculator in super-hist mode
    fluxcalc = new FluxCalculator(Paths::baselines(), Paths::histogram());
    myPred->EnterFluxCalculator(fluxcalc);

    myPred->LoadPredictedIBD(Paths::predicted_ibd());

    ntoys = myPred->LoadToyIBDSpec(toymc_filename);
#pragma omp master
    cout << "number of toy samples: " << ntoys << endl;

    myPred->LoadToyMCEntry(0, false);

    myPred->LoadBgSpec();

    myPred->SetEvisBins(n_evis_bins, evis_bins);
    myPred->SetEnuBins(n_enu_bins, enu_bins);

    myPred->LoadEvisToEnuMatrix(Paths::response());
  }


  DataSet* mydata_nominal = new DataSet();
  mydata_nominal->load(Paths::nominal_toyconfig());
  double sinSq2Theta13 = mydata_nominal->getDouble("sinSq2Theta13");
  double deltam2_ee_default = mydata_nominal->getDouble("deltaMSqee");

  // Double_t sinSq2Theta13 = 0.0;
  // Double_t hierarchy=1;//-1 for inverted
  // Double_t deltam2_32_default=2.41e-3;//eV2
  // Double_t deltam2_21_default=7.50e-5;//eV2
  // Double_t
  // deltam2_31_default=deltam2_32_default+hierarchy*deltam2_21_default;
  // Double_t deltam2_ee_default=deltam2_32_default + (hierarchy*5.17e-5);

  PredSet* predset_nominal = 0;

  //  predset =
  //  pred->MakeOneSuperPrediction(sin22t13[step_dm2][step],dm213[step_dm2][step],false);

  // Double_t M[16*Nstage*n_evis_bins][16*Nstage*n_evis_bins]; //Declaring M
  // inside here causes overstack problem!

  Double_t Nobs_nominal[Nstage][16 * n_evis_bins];
  Double_t Npred_nominal[Nstage][16 * n_evis_bins];

  for (Int_t i = 0; i < 16 * Nstage * n_evis_bins; i++) {
    for (Int_t j = 0; j < 16 * Nstage * n_evis_bins; j++) {
      M[i][j] = 0;
    }
  }

  // First, read nominal values
  myPred->LoadToyMCNominalSpec();
  predset_nominal = myPred->MakeOneSuperPrediction(
      sinSq2Theta13, deltam2_ee_default, 0, -1, false);

  for (int istage = 0; istage < Nstage; ++istage) {
    for (int idet = 4; idet < 8; ++idet) {
      for (int idet2 = 0; idet2 < 4; ++idet2) {
        for (Int_t i = 0; i < n_evis_bins; i++) {
          Int_t iii = ((idet - 4) * 4 + idet2) * n_evis_bins + i;

          // Set nonsensical values for inactive ADs
          if ((NdetectorsConfig[istage][idet] == 0) ||
              (NdetectorsConfig[istage][idet2] == 0))
            Nobs_nominal[istage][iii] = 0;
          else
            Nobs_nominal[istage][iii] =
                myPred->GetCorrEvtsSpec(istage, idet)->GetBinContent(i + 1);

          Npred_nominal[istage][iii] =
              predset_nominal->GetPred(istage, idet, idet2)
                             ->GetBinContent(i + 1);
          //}
          // if(istage==0 && (idet2==3 || idet==7)){
          //  cout << "Near: " << idet2 << " Far: " << idet << endl;
          //  cout << "Nominal: " << i << "\t" << Nobs_nominal[istage][iii] <<
          //  "\t" << Npred_nominal[istage][iii] << endl;
          //}
        }
      }
    }
  }

  if (ntoys > 1000)
    ntoys = 1000;

#pragma omp parallel for
  for (Int_t itoy = 0; itoy < ntoys; itoy++) {
    myPred->LoadToyMCEntry(itoy, true);
    PredSet* predset = myPred->MakeOneSuperPrediction(
        sinSq2Theta13, deltam2_ee_default, 0, -1, false);

    Double_t Nobs[Nstage][16 * n_evis_bins];
    Double_t Npred[Nstage][16 * n_evis_bins];

    for (int istage = 0; istage < Nstage; ++istage) {
      for (int idet = 4; idet < 8; ++idet) {
        for (int idet2 = 0; idet2 < 4; ++idet2) {
          for (Int_t i = 0; i < n_evis_bins; i++) {
            Int_t iii = ((idet - 4) * 4 + idet2) * n_evis_bins + i;

            // Set nonsensical values for inactive ADs
            if ((NdetectorsConfig[istage][idet] == 0) ||
                (NdetectorsConfig[istage][idet2] == 0))
              Nobs[istage][iii] = 0;
            else
              Nobs[istage][iii] =
                  myPred->GetCorrEvtsSpec(istage, idet)->GetBinContent(i + 1);

            Npred[istage][iii] =
                predset->GetPred(istage, idet, idet2)->GetBinContent(i + 1);

            // if (Npred[iii] < 0.01)
            // if(istage==0 && (idet2==3 || idet==7)){
            // cout << "Near: " << idet2 << " Far: " << idet << endl;
            // cout << "Nominal: " << i << "\t" << Nobs_nominal[istage][iii] <<
            // "\t" << Npred_nominal[istage][iii] << endl; cout << "Non-nominal:
            // " << i << "\t" << Nobs[istage][iii] << "\t" << Npred[istage][iii]
            // << endl;
            //}
          }
        }
      }
    }


    for (int istage = 0; istage < Nstage; ++istage) {
      for (Int_t i = 0; i < 16 * n_evis_bins; i++) {
        if (Npred_nominal[istage][i] == 0)
          continue;

        for (int jstage = 0; jstage < Nstage; ++jstage) {
          for (Int_t j = 0; j < 16 * n_evis_bins; j++) {
            if (Npred_nominal[jstage][j] == 0)
              continue;

            //        M[i][j] += (Nobs[i] - Npred[i])*(Nobs[j] -
            //        Npred[j])/Npred[i]/Npred[j];

            if (bkg_flag == 0) {
              if (Npred[istage][i] == 0)
                continue;
              if (Npred[jstage][j] == 0)
                continue;

#pragma omp atomic
              M[i + 16 * n_evis_bins * istage][j + 16 * n_evis_bins * jstage] +=
                  (Nobs[istage][i] - Npred[istage][i] -
                   Nobs_nominal[istage][i] + Npred_nominal[istage][i]) *
                  (Nobs[jstage][j] - Npred[jstage][j] -
                   Nobs_nominal[jstage][j] + Npred_nominal[jstage][j]) /
                  Npred[istage][i] / Npred[jstage][j];

            } else { // do not normalize by the predictions, because the
                     // predicted number of events becomes almost 0 or sometimes
                     // negative for some background-dominated bins. It screws
                     // up the covariance matrix calculation.
#pragma omp atomic
              M[i + 16 * n_evis_bins * istage][j + 16 * n_evis_bins * jstage] +=
                  (Nobs[istage][i] - Npred[istage][i] -
                   Nobs_nominal[istage][i] + Npred_nominal[istage][i]) *
                  (Nobs[jstage][j] - Npred[jstage][j] -
                   Nobs_nominal[jstage][j] + Npred_nominal[jstage][j]);
            }
          }
        }
      }
    }

  } // End of toy loop

  ofstream outf(output_filename);

  // cout << " ============ Covariance Matrics =================" << endl;
  for (Int_t i = 0; i < 16 * Nstage * n_evis_bins; i++) {
    for (Int_t j = 0; j < 16 * Nstage * n_evis_bins; j++) {
      M[i][j] /= (Double_t)ntoys;
      // cout << M[i][j] << endl;
      outf << M[i][j] << "\n";
    }
    // cout << endl;
    outf << "\n";
  }
  // cout << x << endl;
  outf.close();
}
