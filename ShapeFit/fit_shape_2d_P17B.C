#include "Binning.h"
#include "Config.h"
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
#include "TROOT.h"
#include "TStyle.h"
#include "TTree.h"

#include <fstream>
#include <iostream>
#include <map>
#include <stdio.h>
#include <stdlib.h>

// The values we get from Minuit on 2020_01_26@yolo3
// (just for diagnostic comparison)
const Double_t s22t13_nom = 0.08437;
const Double_t dm2ee_nom = 0.0024583;

using namespace Config;

// 4 nModes, 16 MaxPredictions, 3 Nstages, 37 MaxBins
Double_t final_covmatrix[4][16 * 3 * 37][16 * 3 * 37];
Double_t final_covmatrix_sum[4][16 * 37][16 * 37];

Predictor* pred;
#pragma omp threadprivate(pred)

void minuit_fcn(Int_t& npar, Double_t* gin, Double_t& f, Double_t* x,
                Int_t iflag)
{ // function for minuit minimization
  Double_t sin22t13 = x[0];
  Double_t dm2_ee = x[1];
  // cout << "sin22t13: " << sin22t13 << " dm2_ee: " << dm2_ee << endl;
  f = pred->CalculateChi2Cov(sin22t13, dm2_ee, 0, -1);
}

Double_t CalculateChi2(TH1F* h_s22t13, TH1F* h_data);


void fit_shape_2d_P17B(
    Int_t PeriodFlag = -1, //(0=6AD, 1=8AD, 2=7AD, -1=6+8+7AD, 7=6+8AD)
    bool isMC = false)
{
  const char* savefilename = Paths::fit_result();

  const int n_evis_bins = Binning::n_evis();
  double* evis_bins = Binning::evis();
  double* enu_bins = Binning::enu();

  TString sig_spectra_filename[3] = {
      Paths::sig_spectra(0), Paths::sig_spectra(1), Paths::sig_spectra(2)};
  Char_t Theta13InputsLocation[3][1024];

  strcpy(Theta13InputsLocation[0], Paths::input(0));
  strcpy(Theta13InputsLocation[1], Paths::input(1));
  strcpy(Theta13InputsLocation[2], Paths::input(2));

  Double_t stat_factor = 1;

  bool isNominalMC = true;

  int Nevts = 1;
  if (isMC) {
    if (!isNominalMC)
      Nevts = 500;
  }
  gStyle->SetOptStat(0);
  gStyle->SetLabelSize(0.08, "XY");
  gStyle->SetTitleSize(0.04, "XY");

  //  Char_t toymc_filename[256] =
  //  "./toyfiles/nominal/toySpectra_allsys_and_stat.root";

  static FluxCalculator* fluxcalc;
#pragma omp threadprivate(fluxcalc)

#pragma omp parallel
  {
    pred = new Predictor();

    pred->SetStage(PeriodFlag);

    // flux calculator in super-hist mode
    fluxcalc = new FluxCalculator(Paths::baselines(), Paths::histogram());

    pred->EnterFluxCalculator(fluxcalc);

    for (int istage = 0; istage < Nstage; istage++) {
      pred->LoadMainData(Theta13InputsLocation[istage]);
    }

    pred->LoadPredictedIBD(Paths::predicted_ibd());

    pred->LoadIBDSpec(sig_spectra_filename);

    // load bg afterwards since here is when correct events
    pred->LoadBgSpec();


    // pred->SetStatFactor(stat_factor);
    // pred->Set8ADStatFactor(extra_days);

    pred->SetEvisBins(n_evis_bins, evis_bins);
    pred->SetEnuBins(Binning::n_enu(), enu_bins);

    pred->LoadEvisToEnuMatrix(Paths::response());
    // pred->LoadEvisToEnuMatrix("matrix_evis_to_enu_fine");

    pred->LoadCovMatrix(Paths::sig_covmatrix(), Paths::bg_covmatrix());

    //-->Prediction at t13=0
    PredSet* _predset_0 = pred->MakeOneSuperPrediction(0, -1, 0, -1, true);
  }

  //*************************************

  //-->Best fits histograms
  TH2F* h_s22t13_dm2 =
      new TH2F("h_s22t13_dm2", "", 100, 0, 0.2, 100, 1.5e-3, 3.5e-3);
  h_s22t13_dm2->GetXaxis()->SetTitle("sin^{2}(2#theta_{13})");
  h_s22t13_dm2->GetYaxis()->SetTitle("#Delta m^{2}");
  h_s22t13_dm2->SetMarkerStyle(20);
  h_s22t13_dm2->SetMarkerSize(0.8);

  // TH1F* h_noosc = predset_0->GetCombinedPred();
  /*
  //data (combined all "FOUR" far ADs)
  TH1F *h_data[2];

  for(Int_t istage=0;istage<Nstage;istage++){
  h_data[istage] =
  (TH1F*)pred->GetCorrEvtsSpec(istage,4)->Clone(LeakStr("h_data_stage%i",istage));
  h_data[istage]->Add(pred->GetCorrEvtsSpec(istage,5));
  h_data[istage]->Add(pred->GetCorrEvtsSpec(istage,6));
  h_data[istage]->Add(pred->GetCorrEvtsSpec(istage,7));
  cout << "Stage: " << istage+1 << " h_data: " << h_data[istage]->Integral() <<
  endl;
  }
  */
  // Do the minimization manually
  //----------------------------

  //++best parameters
  Double_t minchi2 = 1e8;
  Double_t bests22t13 = 0;
  Double_t bestdm2ee = 0;
  //  PredSet bestpred;

  // const Int_t nsteps = 31;
  // Double_t s22t13start=0.07;
  // Double_t s22t13end=0.10;

  // const Int_t nsteps_dm2 = 31  ;
  // Double_t dm2eestart=2.1e-3;
  // Double_t dm2eeend=2.8e-3;

  // const Int_t nsteps = 101;
  // Double_t s22t13start=0.06;
  // Double_t s22t13end=0.11;

  // const Int_t nsteps_dm2 = 101;
  // Double_t dm2eestart=2.1e-3;
  // Double_t dm2eeend=2.9e-3;

  Double_t chi2result[nsteps_dm2][nsteps];
  Double_t dchi2result[nsteps_dm2][nsteps];
  Double_t sin22t13[nsteps_dm2][nsteps];
  Double_t dm2ee[nsteps_dm2][nsteps];

  Double_t chi2_min = 1e10;
  Double_t s2t_min = 1e10;
  Double_t dm2_min = 1e10;


  TH1F* h_minchi2 = new TH1F("h_minchi2", "", 100, 0, 300);
  h_minchi2->GetXaxis()->SetTitle("Minimum #chi^{2}");
  h_minchi2->GetYaxis()->SetTitle("fake experiments");

  const Double_t s22t13_binwidth = (s22t13end - s22t13start) / (nsteps - 1);
  const Double_t s22t13_histstart = s22t13start - s22t13_binwidth/2;
  const Double_t s22t13_histend = s22t13end + s22t13_binwidth/2;
  const Double_t dm2ee_binwidth = (dm2eeend - dm2eestart) / (nsteps_dm2 - 1);
  const Double_t dm2ee_histstart = dm2eestart - dm2ee_binwidth/2;
  const Double_t dm2ee_histend = dm2eeend + dm2ee_binwidth/2;

  TH2F* h_chi2_all = new TH2F("h_chi2_all", "", nsteps, s22t13_histstart, s22t13_histend,
                              nsteps_dm2, dm2ee_histstart, dm2ee_histend);
  h_chi2_all->GetXaxis()->SetTitle("sin^2(2#theta_{13})");
  h_chi2_all->GetYaxis()->SetTitle("#Delta m^2");
  TH2F* h_chi2_temp = (TH2F*)h_chi2_all->Clone("h_chi2_temp");

  TDirectory* dir = gDirectory;

  // savefile
  TFile* savefile = new TFile(savefilename, "RECREATE");
  TTree* tr = new TTree("tr_fit", "fit results");
  tr->Branch("chi2_map", &dchi2result[0][0],
             LeakStr("chi2_map[%d][%d]/D", nsteps_dm2, nsteps));
  tr->Branch("chi2_min", &chi2_min, "chi2_min/D");
  tr->Branch("dm2_min", &dm2_min, "dm2_min/D");
  tr->Branch("s2t_min", &s2t_min, "s2t_min/D");


  TH2F* h_chi2_map = new TH2F("h_chi2_map", "h_chi2_map", nsteps, s22t13_histstart,
                              s22t13_histend, nsteps_dm2, dm2ee_histstart, dm2ee_histend);


  dir->cd();

  const Int_t NPars = 2;
  TMinuit* minu = new TMinuit(NPars);
  // minu->SetPrintLevel(-1);
  minu->SetFCN(minuit_fcn);


  for (int ievt = 0; ievt < Nevts; ++ievt) {
    minchi2 = 1e10;
    bests22t13 = 0;
    bestdm2ee = 0;
    h_chi2_temp->Reset();

    if (isMC) {
      if (isNominalMC)
        //        pred->LoadToyMCEntry(ievt,true);
        pred->LoadToyMCNominalSpec();
      else
        pred->LoadToyMCEntry(ievt, true);
    }

    Int_t ierflag;
    Double_t arglist[2];
    arglist[0] = 1.0;

    // Fit by Minuit
    minu->mnexcm("SET ERR", arglist, 1, ierflag);

    minu->mnparm(0, "SinSq2Theta13",
                 0.082, // 34
                 0.001, 0, 0.2, ierflag);
    minu->mnparm(1, "DeltaMSqee",
                 0.0024, // 4
                 0.0001, 0.0015, 0.0035, ierflag);

    // minu->FixParameter(0);
    // minu->FixParameter(1);
    arglist[0] = 10000;
    arglist[1] = 1.0;
    //   minu->mnexcm("SIMPLEX", arglist ,2,ierflag);
    //   std::cout << "================ SIMPLEX finished with error flag = " <<
    //   ierflag  << " ============================" << std::endl; if (ierflag
    //   != 0) return;
    arglist[1] = 1.0;
    //    minu->mnexcm("MINIMIZE", arglist ,2,ierflag);

    minu->mnexcm("MIGRAD", arglist, 2, ierflag); // Default fitter
    // std::cout << "================ MIGRAD finished with error flag = " <<
    // ierflag  << " ============================" << std::endl;

    Double_t fpar, ferr;
    minu->GetParameter(0, fpar, ferr);
    bests22t13 = fpar;
    minu->GetParameter(1, fpar, ferr);
    bestdm2ee = fpar;
    Double_t best_pars[2] = {bests22t13, bestdm2ee};
    Double_t* grad = nullptr; // not used by minuit_fcn

    minu->Eval(2, grad, minchi2, best_pars, 0);

    cout << "======== fit results (for checking) :" << bests22t13 << " "
         << bestdm2ee << " " << minchi2 << endl;

#pragma omp parallel for
    for (int step_dm2 = 0; step_dm2 < nsteps_dm2; ++step_dm2) {
#pragma omp critical
      cout << "Running dm2_step=" << step_dm2 << " out of " << nsteps_dm2
           << endl;

      for (int step = 0; step < nsteps; ++step) {
        sin22t13[step_dm2][step] =
            (s22t13end - s22t13start) * step * 1. / (nsteps - 1) + s22t13start;
        dm2ee[step_dm2][step] =
            (dm2eeend - dm2eestart) * step_dm2 * 1. / (nsteps_dm2 - 1) +
            dm2eestart;

        chi2result[step_dm2][step] = pred->CalculateChi2Cov(
            sin22t13[step_dm2][step], dm2ee[step_dm2][step], 0, -1);
#pragma omp critical
        h_chi2_temp->SetBinContent(step + 1, step_dm2 + 1,
                                   chi2result[step_dm2][step]);
      }
    } // loop over t13

    cout << "Event: " << ievt << "; --> best fit: " << bests22t13 << ","
         << bestdm2ee << "; minchi2=" << minchi2 << endl; // tmp
    h_s22t13_dm2->Fill(bests22t13, bestdm2ee);
    h_minchi2->Fill(minchi2);

    // subtract minchi2 from h_chi2_temp
    for (int step_dm2 = 0; step_dm2 < nsteps_dm2; ++step_dm2) {
      for (int step = 0; step < nsteps; ++step) {
        h_chi2_temp->SetBinContent(
            step + 1, step_dm2 + 1,
            h_chi2_temp->GetBinContent(step + 1, step_dm2 + 1) - minchi2);
      }
    }

    h_chi2_all->Add(h_chi2_temp, 1);

    for (int step = 0; step < nsteps; ++step) {
      for (int step_dm2 = 0; step_dm2 < nsteps_dm2; ++step_dm2) {
        h_chi2_map->SetBinContent(step + 1, step_dm2 + 1,
                                  chi2result[step_dm2][step] - minchi2);
        dchi2result[step_dm2][step] = chi2result[step_dm2][step] - minchi2;
      }
    }

    chi2_min = minchi2;
    s2t_min = bests22t13;
    dm2_min = bestdm2ee;

    tr->Fill();
  } // loop over toys

  PredSet* tmp_set;
  PredSet* bestpred = new PredSet();

  cout << "Making best prediction" << endl;
  tmp_set = pred->MakeOneSuperPrediction(bests22t13, bestdm2ee, 0, -1, true);
  cout << "Best prediction done!" << endl;
  for (Int_t istage = 0; istage < Nstage; istage++)
    for (Int_t idet = 4; idet < 8; idet++)
      for (Int_t idet2 = 0; idet2 < 4; idet2++)
        bestpred->SetPred(
            istage, idet, idet2,
            (TH1F*)tmp_set->GetPred(istage, idet, idet2)->Clone());

  TCanvas* c0 = new TCanvas("c0", "c0", 600, 600);
  h_chi2_map->Draw("zcol");


  cout << "======================== Results ============================"
       << endl;
  cout << "Best sin2(2theta13): " << bests22t13 << endl;
  cout << "Best delta m2(ee): " << bestdm2ee << endl;
  cout << "Best chi2: " << minchi2 << endl; // tmp
  cout << "Best prediction: " << endl;

  for (Int_t istage = 0; istage < Nstage; istage++) {
    cout << "Stage#: " << istage + 1 << endl;
    bestpred->PrintToScreen(istage);
  }

  Double_t* tmp_vector;

  const Int_t nModes = 4;
  const Int_t nPredictions[nModes] = {
      16, 2, 4, 1}; // number of predictions to fit for each mode

  Double_t final_pred[nModes][16 * Nstage * n_evis_bins];
  Double_t final_pred_null[nModes][16 * Nstage * n_evis_bins];
  Double_t final_pred_nom[nModes][16 * Nstage * n_evis_bins];
  Double_t final_obs[nModes][16 * Nstage * n_evis_bins];
  Double_t final_obserror[nModes][16 * Nstage * n_evis_bins];

  Double_t final_pred_sum[nModes][16 * n_evis_bins];
  Double_t final_pred_null_sum[nModes][16 * n_evis_bins];
  Double_t final_pred_nom_sum[nModes][16 * n_evis_bins];
  Double_t final_obs_sum[nModes][16 * n_evis_bins];
  Double_t final_obserror_sum[nModes][16 * n_evis_bins];

  savefile->cd();
  TH1D::SetDefaultSumw2();

  TString hnames[nModes][16];

  hnames[0][0] = TString("AD5 / AD1");
  hnames[0][1] = TString("AD5 / AD2");
  hnames[0][2] = TString("AD5 / AD3");
  hnames[0][3] = TString("AD5 / AD4");
  hnames[0][4] = TString("AD6 / AD1");
  hnames[0][5] = TString("AD6 / AD2");
  hnames[0][6] = TString("AD6 / AD3");
  hnames[0][7] = TString("AD6 / AD4");
  hnames[0][8] = TString("AD7 / AD1");
  hnames[0][9] = TString("AD7 / AD2");
  hnames[0][10] = TString("AD7 / AD3");
  hnames[0][11] = TString("AD7 / AD4");
  hnames[0][12] = TString("AD8 / AD1");
  hnames[0][13] = TString("AD8 / AD2");
  hnames[0][14] = TString("AD8 / AD3");
  hnames[0][15] = TString("AD8 / AD4");

  hnames[1][0] = TString("AD5+AD6+AD7+AD8 / AD1+AD2");
  hnames[1][1] = TString("AD5+AD6+AD7+AD8 / AD3+AD4");

  hnames[2][0] = TString("AD5+AD6+AD7+AD8 / AD1");
  hnames[2][1] = TString("AD5+AD6+AD7+AD8 / AD2");
  hnames[2][2] = TString("AD5+AD6+AD7+AD8 / AD3");
  hnames[2][3] = TString("AD5+AD6+AD7+AD8 / AD4");

  hnames[3][0] = TString("AD5+AD6+AD7+AD8 / AD1+AD2+AD3+AD4");


  TH1D* h_final_obs[Nstage][nModes][16];
  TH1D* h_final_pred[Nstage][nModes][16];
  TH1D* h_final_pred_null[Nstage][nModes][16];
  TH1D* h_final_pred_nom[Nstage][nModes][16];

  TH1D* h_final_obs_sum[nModes][16];
  TH1D* h_final_pred_sum[nModes][16];
  TH1D* h_final_pred_null_sum[nModes][16];
  TH1D* h_final_pred_nom_sum[nModes][16];

  TH2D* h_final_covmatrix[nModes]; // actually the correlation matrix
  TH2D* h_final_covmatrix_really[nModes]; // this is really the covmatrix

  Int_t n_rebinned_evis_bins = pred->GetNumRebinnedEvisBins();
  Double_t* rebinned_evis_bins = pred->GetRebinnedEvisBins();


  for (Int_t iMode = 0; iMode < nModes; iMode++) {
    h_final_covmatrix[iMode] = new TH2D(
        LeakStr("h_final_covmatrix_mode%d", iMode), "Final covariance matrix",
        nPredictions[iMode] * Nstage * n_rebinned_evis_bins, 0,
        nPredictions[iMode] * Nstage * n_rebinned_evis_bins,
        nPredictions[iMode] * Nstage * n_rebinned_evis_bins, 0,
        nPredictions[iMode] * Nstage * n_rebinned_evis_bins);
    h_final_covmatrix_really[iMode] = new TH2D(
        LeakStr("h_final_covmatrix_really_mode%d", iMode), "Final covariance matrix (really)",
        nPredictions[iMode] * Nstage * n_rebinned_evis_bins, 0,
        nPredictions[iMode] * Nstage * n_rebinned_evis_bins,
        nPredictions[iMode] * Nstage * n_rebinned_evis_bins, 0,
        nPredictions[iMode] * Nstage * n_rebinned_evis_bins);

    for (Int_t i = 0; i < nPredictions[iMode]; i++) {
      h_final_obs_sum[iMode][i] = new TH1D(
          LeakStr("h_final_obs_sum_mode%d_%d", iMode, i), hnames[iMode][i].Data(),
          n_rebinned_evis_bins, rebinned_evis_bins);
      h_final_obs_sum[iMode][i]->GetXaxis()->SetTitle("Visible energy");
      h_final_obs_sum[iMode][i]->GetYaxis()->SetTitle("Events per day");

      h_final_pred_sum[iMode][i] = new TH1D(
          LeakStr("h_final_pred_sum_mode%d_%d", iMode, i), hnames[iMode][i].Data(),
          n_rebinned_evis_bins, rebinned_evis_bins);
      h_final_pred_null_sum[iMode][i] = new TH1D(
          LeakStr("h_final_pred_null_sum_mode%d_%d", iMode, i),
          hnames[iMode][i].Data(), n_rebinned_evis_bins, rebinned_evis_bins);
      h_final_pred_nom_sum[iMode][i] = new TH1D(
          LeakStr("h_final_pred_nom_sum_mode%d_%d", iMode, i),
          hnames[iMode][i].Data(), n_rebinned_evis_bins, rebinned_evis_bins);

      for (Int_t istage = 0; istage < Nstage; istage++) {
        h_final_obs[istage][iMode][i] = new TH1D(
            LeakStr("h_final_obs_stage%d_mode%d_%d", istage, iMode, i),
            hnames[iMode][i].Data(), n_rebinned_evis_bins, rebinned_evis_bins);
        h_final_obs[istage][iMode][i]->GetXaxis()->SetTitle("Visible energy");
        h_final_obs[istage][iMode][i]->GetYaxis()->SetTitle("Events per day");

        h_final_pred[istage][iMode][i] = new TH1D(
            LeakStr("h_final_pred_stage%d_mode%d_%d", istage, iMode, i),
            hnames[iMode][i].Data(), n_rebinned_evis_bins, rebinned_evis_bins);
        h_final_pred_null[istage][iMode][i] = new TH1D(
            LeakStr("h_final_pred_null_stage%d_mode%d_%d", istage, iMode, i),
            hnames[iMode][i].Data(), n_rebinned_evis_bins, rebinned_evis_bins);
        h_final_pred_nom[istage][iMode][i] = new TH1D(
            LeakStr("h_final_pred_nom_stage%d_mode%d_%d", istage, iMode, i),
            hnames[iMode][i].Data(), n_rebinned_evis_bins, rebinned_evis_bins);
      }
    }
  }

  TH1D* h_final_obs_ratio[Nstage][nModes][16];
  TH1D* h_final_pred_ratio[Nstage][nModes][16];
  TH1D* h_final_pred_null_ratio[Nstage][nModes][16];
  TH1D* h_final_pred_nom_ratio[Nstage][nModes][16];

  TH1D* h_final_obs_ratio_sum[nModes][16];
  TH1D* h_final_pred_ratio_sum[nModes][16];
  TH1D* h_final_pred_null_ratio_sum[nModes][16];
  TH1D* h_final_pred_nom_ratio_sum[nModes][16];
  dir->cd();

  PredSet* nullpred = new PredSet();

  cout << "Making zero prediction" << endl;
  tmp_set = pred->MakeOneSuperPrediction(0, -1, 0, -1, true);
  for (Int_t istage = 0; istage < Nstage; istage++)
    for (Int_t idet = 4; idet < 8; idet++)
      for (Int_t idet2 = 0; idet2 < 4; idet2++)
        nullpred->SetPred(
            istage, idet, idet2,
            (TH1F*)tmp_set->GetPred(istage, idet, idet2)->Clone());

  PredSet* nompred = new PredSet();

  cout << "Making nominal prediction" << endl;
  tmp_set = pred->MakeOneSuperPrediction(s22t13_nom, dm2ee_nom, 0, -1, true);
  for (Int_t istage = 0; istage < Nstage; istage++)
    for (Int_t idet = 4; idet < 8; idet++)
      for (Int_t idet2 = 0; idet2 < 4; idet2++)
        nompred->SetPred(
            istage, idet, idet2,
            (TH1F*)tmp_set->GetPred(istage, idet, idet2)->Clone());

  // try to dump final covariance matrix
  ofstream fout_cmat("final_covmatrix_mode1.txt");
  // -1 => both near and far site stat errors
  tmp_vector = pred->GetFinalCovMatrix(-1, 1);
  for (Int_t i = 0; i < nPredictions[1] * Nstage * n_rebinned_evis_bins; i++) {
    for (Int_t j = 0; j < nPredictions[1] * Nstage * n_rebinned_evis_bins;
         j++) {
      fout_cmat
          << " "
          << tmp_vector[i * nPredictions[1] * Nstage * n_rebinned_evis_bins +
                        j];
    }
    fout_cmat << "\n";
  }
  fout_cmat.close();


  /////

  for (Int_t iMode = 0; iMode < nModes; iMode++) {
    // cout << "iMode: " << iMode << endl;

    // 0 => no far site stat error
    tmp_vector = pred->GetFinalCovMatrix(0, iMode);

    for (Int_t i = 0; i < nPredictions[iMode] * n_rebinned_evis_bins; i++) {
      for (Int_t j = 0; j < nPredictions[iMode] * n_rebinned_evis_bins; j++) {
        final_covmatrix_sum[iMode][i][j] = 0;
      }
    }

    for (Int_t istage = 0; istage < Nstage; istage++) {
      for (Int_t jstage = 0; jstage < Nstage; jstage++) {
        for (Int_t i = 0; i < nPredictions[iMode] * n_rebinned_evis_bins; i++) {
          for (Int_t j = 0; j < nPredictions[iMode] * n_rebinned_evis_bins;
               j++) {
            final_covmatrix
                [iMode][istage * nPredictions[iMode] * n_rebinned_evis_bins + i]
                [jstage * nPredictions[iMode] * n_rebinned_evis_bins +
                 j] = tmp_vector
                    [(istage * nPredictions[iMode] * n_rebinned_evis_bins + i) *
                         nPredictions[iMode] * Nstage * n_rebinned_evis_bins +
                     jstage * nPredictions[iMode] * n_rebinned_evis_bins + j];

            //      cout << " " << final_covmatrix[i][j];
          }
          //    cout << endl;
        }
      }
    }

    tmp_vector = pred->GetFinalCovMatrixSum(0, iMode);
    for (Int_t i = 0; i < nPredictions[iMode] * n_rebinned_evis_bins; i++) {
      for (Int_t j = 0; j < nPredictions[iMode] * n_rebinned_evis_bins; j++) {
        final_covmatrix_sum[iMode][i][j] +=
            tmp_vector[i * nPredictions[iMode] * n_rebinned_evis_bins + j];
      }
    }


    tmp_vector = pred->GetFinalPred(bestpred, iMode);
    for (Int_t i = 0; i < nPredictions[iMode] * Nstage * n_rebinned_evis_bins;
         i++) {
      final_pred[iMode][i] = tmp_vector[i];
    }

    tmp_vector = pred->GetFinalObs(iMode);
    for (Int_t i = 0; i < nPredictions[iMode] * Nstage * n_rebinned_evis_bins;
         i++) {
      final_obs[iMode][i] = tmp_vector[i];
    }

    tmp_vector = pred->GetFinalObsError(iMode);
    for (Int_t i = 0; i < nPredictions[iMode] * Nstage * n_rebinned_evis_bins;
         i++) {
      final_obserror[iMode][i] = tmp_vector[i];
    }

    tmp_vector = pred->GetFinalPred(nullpred, iMode);
    for (Int_t i = 0; i < nPredictions[iMode] * Nstage * n_rebinned_evis_bins;
         i++) {
      final_pred_null[iMode][i] = tmp_vector[i];
    }

    tmp_vector = pred->GetFinalPred(nompred, iMode);
    for (Int_t i = 0; i < nPredictions[iMode] * Nstage * n_rebinned_evis_bins;
         i++) {
      final_pred_nom[iMode][i] = tmp_vector[i];
    }

    tmp_vector = pred->GetFinalPredSum(bestpred, iMode);
    for (Int_t i = 0; i < nPredictions[iMode] * n_rebinned_evis_bins; i++) {
      final_pred_sum[iMode][i] = tmp_vector[i];
    }

    tmp_vector = pred->GetFinalObsSum(iMode);
    for (Int_t i = 0; i < nPredictions[iMode] * n_rebinned_evis_bins; i++) {
      final_obs_sum[iMode][i] = tmp_vector[i];
    }

    tmp_vector = pred->GetFinalObsErrorSum(iMode);
    for (Int_t i = 0; i < nPredictions[iMode] * n_rebinned_evis_bins; i++) {
      final_obserror_sum[iMode][i] = tmp_vector[i];
    }

    tmp_vector = pred->GetFinalPredSum(nullpred, iMode);
    for (Int_t i = 0; i < nPredictions[iMode] * n_rebinned_evis_bins; i++) {
      final_pred_null_sum[iMode][i] = tmp_vector[i];
    }

    tmp_vector = pred->GetFinalPredSum(nompred, iMode);
    for (Int_t i = 0; i < nPredictions[iMode] * n_rebinned_evis_bins; i++) {
      final_pred_nom_sum[iMode][i] = tmp_vector[i];
    }

    for (Int_t iPred = 0; iPred < nPredictions[iMode]; iPred++) {
      for (Int_t istage = 0; istage < Nstage; istage++) {
        for (Int_t i = 0; i < n_rebinned_evis_bins; i++) {
          // cout << "evis: " << i << endl;

          h_final_obs[istage][iMode][iPred]->SetBinContent(
              i + 1, final_obs[iMode][(istage * nPredictions[iMode] + iPred) *
                                          n_rebinned_evis_bins +
                                      i]);
          h_final_obs[istage][iMode][iPred]->SetBinError(
              i + 1,
              final_obserror[iMode][(istage * nPredictions[iMode] + iPred) *
                                        n_rebinned_evis_bins +
                                    i]);

          // cout << "h_final_obs: " <<
          // final_obs[iMode][(istage*nPredictions[iMode]+iPred)*n_rebinned_evis_bins+i]
          // << ", " <<
          // final_obserror[iMode][(istage*nPredictions[iMode]+iPred)*n_rebinned_evis_bins+i]
          // << endl;

          h_final_pred[istage][iMode][iPred]->SetBinContent(
              i + 1, final_pred[iMode][(istage * nPredictions[iMode] + iPred) *
                                           n_rebinned_evis_bins +
                                       i]);
          h_final_pred[istage][iMode][iPred]->SetBinError(
              i + 1,
              sqrt(final_covmatrix[iMode]
                                  [(istage * nPredictions[iMode] + iPred) *
                                       n_rebinned_evis_bins +
                                   i][(istage * nPredictions[iMode] + iPred) *
                                          n_rebinned_evis_bins +
                                      i]));
          // cout << "h_final_pred: " <<
          // final_pred[iMode][(istage*nPredictions[iMode]+iPred)*n_rebinned_evis_bins+i]
          // << ", " <<
          // sqrt(final_covmatrix[iMode][(istage*nPredictions[iMode]+iPred)*n_rebinned_evis_bins+i][(istage*nPredictions[iMode]+iPred)*n_rebinned_evis_bins+i])
          // << endl;

          h_final_pred_null[istage][iMode][iPred]->SetBinContent(
              i + 1,
              final_pred_null[iMode][(istage * nPredictions[iMode] + iPred) *
                                         n_rebinned_evis_bins +
                                     i]);
          h_final_pred_null[istage][iMode][iPred]->SetBinError(i + 1, 0);

          // cout << "h_final_pred_null: " <<
          // final_pred_null[iMode][(istage*nPredictions[iMode]+iPred)*n_rebinned_evis_bins+i]
          // << endl;

          h_final_pred_nom[istage][iMode][iPred]->SetBinContent(
              i + 1,
              final_pred_nom[iMode][(istage * nPredictions[iMode] + iPred) *
                                         n_rebinned_evis_bins +
                                     i]);
          h_final_pred_nom[istage][iMode][iPred]->SetBinError(i + 1, 0);
        }


        h_final_pred[istage][iMode][iPred]->SetLineStyle(2);
        h_final_pred[istage][iMode][iPred]->SetLineColor(2);
        h_final_pred[istage][iMode][iPred]->SetFillStyle(1001);
        //  h_final_pred[iMode][iPred]->SetFillColor(kGray);
        h_final_pred[istage][iMode][iPred]->SetFillColor(2);
        h_final_pred_null[istage][iMode][iPred]->SetLineColor(4);
        h_final_pred_nom[istage][iMode][iPred]->SetLineColor(5);

        savefile->cd();

        h_final_obs_ratio[istage][iMode][iPred] =
            (TH1D*)h_final_obs[istage][iMode][iPred]->Clone(LeakStr(
                "h_final_obs_ratio_stage%d_mode%d_%d", istage, iMode, iPred));
        h_final_obs_ratio[istage][iMode][iPred]->GetYaxis()->SetTitle(
            "Ratio to null osc. prediction");
        h_final_pred_ratio[istage][iMode][iPred] =
            (TH1D*)h_final_pred[istage][iMode][iPred]->Clone(LeakStr(
                "h_final_pred_ratio_stage%d_mode%d_%d", istage, iMode, iPred));
        h_final_pred_null_ratio[istage][iMode][iPred] =
            (TH1D*)h_final_pred_null[istage][iMode][iPred]->Clone(
                LeakStr("h_final_pred_null_ratio_stage%d_mode%d_%d", istage, iMode,
                     iPred));
        h_final_pred_nom_ratio[istage][iMode][iPred] =
            (TH1D*)h_final_pred_nom[istage][iMode][iPred]->Clone(
                LeakStr("h_final_pred_nom_ratio_stage%d_mode%d_%d", istage, iMode,
                     iPred));

        if (h_final_pred_null[istage][iMode][iPred]->Integral() >
            0) { // Check to not do division by zero
          h_final_obs_ratio[istage][iMode][iPred]->Divide(
              h_final_pred_null[istage][iMode][iPred]);
          h_final_pred_ratio[istage][iMode][iPred]->Divide(
              h_final_pred_null[istage][iMode][iPred]);
          h_final_pred_null_ratio[istage][iMode][iPred]->Divide(
              h_final_pred_null[istage][iMode][iPred]);
          h_final_pred_nom_ratio[istage][iMode][iPred]->Divide(
              h_final_pred_null[istage][iMode][iPred]);
        } else {
          h_final_obs_ratio[istage][iMode][iPred]->Scale(0);
          h_final_pred_ratio[istage][iMode][iPred]->Scale(0);
          h_final_pred_null_ratio[istage][iMode][iPred]->Scale(0);
          h_final_pred_nom_ratio[istage][iMode][iPred]->Scale(0);
        }

        dir->cd();
      }
    }

    for (Int_t iPred = 0; iPred < nPredictions[iMode]; iPred++) {
      for (Int_t i = 0; i < n_rebinned_evis_bins; i++) {
        h_final_obs_sum[iMode][iPred]->SetBinContent(
            i + 1, final_obs_sum[iMode][iPred * n_rebinned_evis_bins + i]);
        h_final_obs_sum[iMode][iPred]->SetBinError(
            i + 1, final_obserror_sum[iMode][iPred * n_rebinned_evis_bins + i]);

        h_final_pred_sum[iMode][iPred]->SetBinContent(
            i + 1, final_pred_sum[iMode][iPred * n_rebinned_evis_bins + i]);
        h_final_pred_sum[iMode][iPred]->SetBinError(
            i + 1,
            sqrt(final_covmatrix_sum[iMode][iPred * n_rebinned_evis_bins + i]
                                    [iPred * n_rebinned_evis_bins + i]));

        h_final_pred_null_sum[iMode][iPred]->SetBinContent(
            i + 1,
            final_pred_null_sum[iMode][iPred * n_rebinned_evis_bins + i]);
        h_final_pred_null_sum[iMode][iPred]->SetBinError(i + 1, 0);

        h_final_pred_nom_sum[iMode][iPred]->SetBinContent(
            i + 1,
            final_pred_nom_sum[iMode][iPred * n_rebinned_evis_bins + i]);
        h_final_pred_nom_sum[iMode][iPred]->SetBinError(i + 1, 0);
      }


      h_final_pred_sum[iMode][iPred]->SetLineStyle(2);
      h_final_pred_sum[iMode][iPred]->SetLineColor(2);
      h_final_pred_sum[iMode][iPred]->SetFillStyle(1001);
      //  h_final_pred[iMode][iPred]->SetFillColor(kGray);
      h_final_pred_sum[iMode][iPred]->SetFillColor(2);
      h_final_pred_null_sum[iMode][iPred]->SetLineColor(4);
      h_final_pred_nom_sum[iMode][iPred]->SetLineColor(5);

      savefile->cd();

      h_final_obs_ratio_sum[iMode][iPred] =
          (TH1D*)h_final_obs_sum[iMode][iPred]->Clone(
              LeakStr("h_final_obs_ratio_sum_mode%d_%d", iMode, iPred));
      h_final_obs_ratio_sum[iMode][iPred]->GetYaxis()->SetTitle(
          "Ratio to null osc. prediction");
      h_final_pred_ratio_sum[iMode][iPred] =
          (TH1D*)h_final_pred_sum[iMode][iPred]->Clone(
              LeakStr("h_final_pred_ratio_sum_mode%d_%d", iMode, iPred));
      h_final_pred_null_ratio_sum[iMode][iPred] =
          (TH1D*)h_final_pred_null_sum[iMode][iPred]->Clone(
              LeakStr("h_final_pred_null_ratio_sum_mode%d_%d", iMode, iPred));
      h_final_pred_nom_ratio_sum[iMode][iPred] =
          (TH1D*)h_final_pred_nom_sum[iMode][iPred]->Clone(
              LeakStr("h_final_pred_nom_ratio_sum_mode%d_%d", iMode, iPred));

      if (h_final_pred_null_sum[iMode][iPred]->Integral() >
          0) { // Check to not do division by zero
        h_final_obs_ratio_sum[iMode][iPred]->Divide(
            h_final_pred_null_sum[iMode][iPred]);
        h_final_pred_ratio_sum[iMode][iPred]->Divide(
            h_final_pred_null_sum[iMode][iPred]);
        h_final_pred_null_ratio_sum[iMode][iPred]->Divide(
            h_final_pred_null_sum[iMode][iPred]);
        h_final_pred_nom_ratio_sum[iMode][iPred]->Divide(
            h_final_pred_null_sum[iMode][iPred]);
      } else {
        h_final_obs_ratio_sum[iMode][iPred]->Scale(0);
        h_final_pred_ratio_sum[iMode][iPred]->Scale(0);
        h_final_pred_null_ratio_sum[iMode][iPred]->Scale(0);
        h_final_pred_nom_ratio_sum[iMode][iPred]->Scale(0);
      }

      dir->cd();
    }


    for (Int_t i = 0; i < nPredictions[iMode] * Nstage * n_rebinned_evis_bins;
         i++) {
      for (Int_t j = 0; j < nPredictions[iMode] * Nstage * n_rebinned_evis_bins;
           j++) {
        if (final_covmatrix[iMode][i][i] == 0)
          continue;
        if (final_covmatrix[iMode][j][j] == 0)
          continue;

        h_final_covmatrix[iMode]->SetBinContent(
            i + 1, j + 1,
            final_covmatrix[iMode][i][j] / sqrt(final_covmatrix[iMode][i][i] *
                                                final_covmatrix[iMode][j][j]));
        h_final_covmatrix_really[iMode]->SetBinContent(
            i + 1, j + 1,
            final_covmatrix[iMode][i][j]);
      }
    }
  }
  /*
    TCanvas * c2 = new TCanvas("c2","c2",1200,1000);

    TCanvas * c3 = new TCanvas("c3","c3",600,600);

    c2->Print("./fit_result_files/fit_shape_2d.pdf["); //Open PDF file

    for (Int_t iMode = 0; iMode < nModes; iMode++){
    c3->Clear();
    for(Int_t istage=0;istage<Nstage;istage++){
    c2->Clear();

    if (iMode == 0){
    c2->Divide(4,4);
    }
    if (iMode == 1){
    c2->SetWindowSize(1200,600);
    c2->Divide(2,1);
    }
    if (iMode == 2){
    c2->SetWindowSize(1200,1000);
    c2->Divide(2,2);
    }
    if (iMode == 3){
    c2->SetWindowSize(1200,1000);
    c2->Divide(1,1);
    }

    for (Int_t iPred = 0; iPred <nPredictions[iMode]; iPred++){
    c2->cd(iPred+1);
    TVirtualPad * pad = c2->GetPad(iPred+1);
    pad->Divide(1,2,0,0);
    pad->cd(1);
    pad->GetPad(1)->SetRightMargin(0.05);;
    pad->GetPad(1)->SetLeftMargin(0.15);;
    h_final_obs[istage][iMode][iPred]->SetTitle(";;Events per day");
    h_final_obs[istage][iMode][iPred]->SetTitleOffset(1.5,"Y");
    h_final_obs[istage][iMode][iPred]->GetXaxis()->SetRangeUser(0.7,7.9);
    h_final_obs[istage][iMode][iPred]->Draw();
    h_final_pred_null[istage][iMode][iPred]->Draw("hist same");
    h_final_pred[istage][iMode][iPred]->Draw("e2 same");
    h_final_obs[istage][iMode][iPred]->Draw("same");

    pad->cd(2);
    pad->GetPad(2)->SetRightMargin(0.05);;
    pad->GetPad(2)->SetLeftMargin(0.15);;
    h_final_pred_ratio[istage][iMode][iPred]->SetTitle(";;Ratio");
    h_final_pred_ratio[istage][iMode][iPred]->GetXaxis()->SetRangeUser(0.7,7.9);
    h_final_pred_ratio[istage][iMode][iPred]->SetTitleOffset(1.5,"Y");
    h_final_pred_ratio[istage][iMode][iPred]->SetMinimum(0.8);
    h_final_pred_ratio[istage][iMode][iPred]->SetMaximum(1.2);
    h_final_pred_ratio[istage][iMode][iPred]->Draw("e2");
    h_final_pred_null_ratio[istage][iMode][iPred]->Draw("hist same");
    h_final_obs_ratio[istage][iMode][iPred]->Draw("e0 same");

    }

    c2->Print("./fit_result_files/fit_shape_2d.pdf");
    }

    c2->Clear();

    if (iMode == 0){
    c2->Divide(4,4);
    }
    if (iMode == 1){
    c2->SetWindowSize(1200,600);
    c2->Divide(2,1);
    }
    if (iMode == 2){
    c2->SetWindowSize(1200,1000);
    c2->Divide(2,2);
    }
    if (iMode == 3){
    c2->SetWindowSize(1200,1000);
    c2->Divide(1,1);
    }

    for (Int_t iPred = 0; iPred <nPredictions[iMode]; iPred++){
    c2->cd(iPred+1);
    TVirtualPad * pad = c2->GetPad(iPred+1);
    pad->Divide(1,2,0,0);
    pad->cd(1);
    pad->GetPad(1)->SetRightMargin(0.05);;
    pad->GetPad(1)->SetLeftMargin(0.15);;
    h_final_obs_sum[iMode][iPred]->SetTitle(";;Events per day");
    h_final_obs_sum[iMode][iPred]->SetTitleOffset(1.5,"Y");
    h_final_obs_sum[iMode][iPred]->GetXaxis()->SetRangeUser(0.7,7.9);
    h_final_obs_sum[iMode][iPred]->Draw();
    h_final_pred_null_sum[iMode][iPred]->Draw("hist same");
    h_final_pred_sum[iMode][iPred]->Draw("e2 same");
    h_final_obs_sum[iMode][iPred]->Draw("same");

    pad->cd(2);
    pad->GetPad(2)->SetRightMargin(0.05);;
    pad->GetPad(2)->SetLeftMargin(0.15);;
    h_final_pred_ratio_sum[iMode][iPred]->SetTitle(";;Ratio");
    h_final_pred_ratio_sum[iMode][iPred]->GetXaxis()->SetRangeUser(0.7,7.9);
    h_final_pred_ratio_sum[iMode][iPred]->SetTitleOffset(1.5,"Y");
    h_final_pred_ratio_sum[iMode][iPred]->SetMinimum(0.8);
    h_final_pred_ratio_sum[iMode][iPred]->SetMaximum(1.2);
    h_final_pred_ratio_sum[iMode][iPred]->Draw("e2");
    h_final_pred_null_ratio_sum[iMode][iPred]->Draw("hist same");
    h_final_obs_ratio_sum[iMode][iPred]->Draw("e0 same");

    }

    c2->Print("./fit_result_files/fit_shape_2d.pdf");

    //c3->cd();
    //h_final_covmatrix[iMode]->Draw("zcol");
    //c3->Print("./fit_result_files/fit_shape_2d.pdf");

    }

    c2->Print("./fit_result_files/fit_shape_2d.pdf]"); //Close PDF file
  */
  /////////////////////////////////////////////////////////////////////////////////////////


  TH2F* h_chi2_all_cl = (TH2F*)h_chi2_all->Clone("h_chi2_all_cl");
  if (Nevts > 1) {
    TCanvas* can_bestfits = new TCanvas("can_bestfits", "Best fits");
    can_bestfits->Divide(2, 2);
    can_bestfits->cd(1);
    h_s22t13_dm2->SetStats(0);
    h_s22t13_dm2->Draw();
    can_bestfits->cd(2);
    h_minchi2->Draw();
    can_bestfits->cd(3);
    h_chi2_all->SetStats(0);
    h_chi2_all->Scale(1. / Nevts);
    h_chi2_all->Draw("colz");
    can_bestfits->cd(4);
    double contours[2];
    contours[0] = 2.3;
    contours[1] = 4.61; // for 3sigma is 9.21
    // double colors[2]={4,2};
    h_chi2_all_cl->SetContour(2, contours);
    h_chi2_all_cl->Draw("cont1");
    h_s22t13_dm2->Draw("same");
  }

  savefile->cd();
  h_chi2_map->Write();
  // chi2graph->Write();
  // dchi2graph->Write();
  h_s22t13_dm2->Write();
  h_minchi2->Write();
  savefile->Write();
  h_chi2_all->Write();
  h_chi2_all_cl->Write();

  // Save corrected events spectra
  Char_t name[1024];
  for (Int_t istage = 0; istage < Nstage; istage++) {
    for (Int_t idet = 0; idet < Ndetectors; idet++) {
      sprintf(name, "CorrIBDEvts_stage%i_AD%i", istage, idet + 1);
      pred->GetCorrEvtsSpec(istage, idet)->Clone(name)->Write();

      sprintf(name, "CorrAccEvtsSpec_stage%i_AD%i", istage, idet + 1);
      pred->GetCorrAccEvtsSpec(istage, idet)->Clone(name)->Write();

      sprintf(name, "CorrLi9EvtsSpec_stage%i_AD%i", istage, idet + 1);
      pred->GetCorrLi9EvtsSpec(istage, idet)->Clone(name)->Write();

      sprintf(name, "CorrAmcEvtsSpec_stage%i_AD%i", istage, idet + 1);
      pred->GetCorrAmcEvtsSpec(istage, idet)->Clone(name)->Write();

      sprintf(name, "CorrFnEvtsSpec_stage%i_AD%i", istage, idet + 1);
      pred->GetCorrFnEvtsSpec(istage, idet)->Clone(name)->Write();

      sprintf(name, "CorrAlnEvtsSpec_stage%i_AD%i", istage, idet + 1);
      pred->GetCorrAlnEvtsSpec(istage, idet)->Clone(name)->Write();

      sprintf(name, "CorrMuonDecayEvtsSpec_stage%i_AD%i", istage, idet + 1);
      pred->GetCorrMuonDecayEvtsSpec(istage, idet)->Clone(name)->Write();
    }
  }

  for (Int_t iMode = 0; iMode < nModes; iMode++)
    h_final_covmatrix[iMode]->Write();

  savefile->Close();

  cout << "The end!" << endl;

  // Something (destructors?) making it take forever to exit, so...
  quick_exit(0);


  // dump final stat error and weighted mean coefficients for later use
  /*
    ofstream fout_wmean("wmean_coeff.txt");
    for (Int_t i = 0; i < n_rebinned_evis_bins; i++){
    for (Int_t j = 0; j < 4; j++){
    fout_wmean << "\t"<<pred->GetWmeanCoeff(i,j);
    }
    fout_wmean << endl;
    }
    fout_wmean.close();



    ofstream fout_stat("covmatrix_nearsite_stat.txt");
    tmp_vector = pred->GetFinalCovMatrix(0,bestpred,0);
    for (Int_t i = 0; i < nPredictions[0]*n_rebinned_evis_bins; i++){
    for (Int_t j = 0; j < nPredictions[0]*n_rebinned_evis_bins; j++){
    fout_stat << pred->GetFracStatError(i,j) << endl;

    }
    fout_stat << endl;
    }
    fout_stat.close();


    ofstream fout_stat_total("covmatrix_total_stat.txt");
    tmp_vector = pred->GetFinalCovMatrix(-1,bestpred,0);
    for (Int_t i = 0; i < nPredictions[0]*n_rebinned_evis_bins; i++){
    for (Int_t j = 0; j < nPredictions[0]*n_rebinned_evis_bins; j++){
    fout_stat_total << pred->GetFracStatError(i,j) << endl;

    }
    fout_stat_total << endl;
    }
    fout_stat_total.close();

    ofstream fout_bg_total("covmatrix_frac_bg.txt");
    for (Int_t i = 0; i < nPredictions[0]*n_rebinned_evis_bins; i++){
    for (Int_t j = 0; j < nPredictions[0]*n_rebinned_evis_bins; j++){
    fout_bg_total << pred->GetFracBgError(i,j) << endl;

    }
    fout_bg_total << endl;
    }
    fout_bg_total.close();
  */

} // end of Test

Double_t CalculateChi2(TH1F* h_s22t13, TH1F* h_data)
{
  Double_t chi2out = 0;
  int nbins = 0;
  for (int ibin = 1; ibin < h_s22t13->GetXaxis()->FindBin(7); ++ibin) {
    nbins++;

    double bdata = h_data->GetBinContent(ibin);
    double bs22t13 = h_s22t13->GetBinContent(ibin);
    double edata = h_data->GetBinError(ibin);
    if (edata != 0 && bdata != 0) {
      chi2out += pow(bdata - bs22t13, 2) * 1. / pow(edata, 2);
    }
  }

  return chi2out; //*1./(nbins-1);
}
