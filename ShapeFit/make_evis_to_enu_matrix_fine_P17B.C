#include "Binning.h"
#include "Config.h"

#include <TFile.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TString.h>
#include <fstream>
#include <iostream>

using namespace Config;
using namespace std;

void make_evis_to_enu_matrix_fine_P17B()
{
  TString input_file = response_root_filename;
  TString output_file = response_filename;
  TString output_file_rateonly = response_filename_rateonly;

  const Int_t n_evis_bins = Binning::n_evis();
  const double* evis_bins = Binning::evis();

  const Int_t n_enu_bins = Binning::n_enu();
  const double* enu_bins = Binning::enu();

  //  TH1F * h_enu_tmp = new TH1F("h_enu_tmp","h_enu_tmp",39,enu_bins);

  TFile* f = new TFile(input_file.Data());
  TH2F* h_orig = (TH2F*)f->Get("h_evis_vs_enu_ad1");
  ofstream fout(output_file.Data());

  for (Int_t iEvisBin = 0; iEvisBin < n_evis_bins; iEvisBin++) {
    //    h_enu_tmp->Reset();
    TH1D* htmp = h_orig->ProjectionX(Form("h_enu_%d", iEvisBin), iEvisBin + 1,
                                     iEvisBin + 1);

    TH1D* htmp_rebin = (TH1D*)htmp->Rebin(
        n_enu_bins, Form("h_enu_rebin_%d", iEvisBin), enu_bins);
    Double_t norm = htmp_rebin->Integral();
    htmp_rebin->Scale(1. / norm);


    for (Int_t iEnuBin = 0; iEnuBin < n_enu_bins; iEnuBin++) {
      cout << htmp_rebin->GetBinContent(iEnuBin + 1) << endl;
      fout << htmp_rebin->GetBinContent(iEnuBin + 1) << endl;
    }
    cout << endl;
    fout << endl;

    // // for debugging
    // TCanvas * c1 = new TCanvas();
    // htmp_rebin->Draw();
    // c1->Update();
    // Char_t dummy[20];
    // cin >> dummy;
  }
  fout.close();

  ofstream fout_rateonly(output_file_rateonly.Data());
  for (Int_t iEvisBin = 0; iEvisBin < n_evis_bins; iEvisBin++) {
    //    h_enu_tmp->Reset();
    TH1D* htmp = h_orig->ProjectionX(
        Form("h_enu_%d", iEvisBin)); // add everything together

    TH1D* htmp_rebin = (TH1D*)htmp->Rebin(
        n_enu_bins, Form("h_enu_rebin_%d", iEvisBin), enu_bins);
    Double_t norm = htmp_rebin->Integral();
    htmp_rebin->Scale(1. / norm);


    for (Int_t iEnuBin = 0; iEnuBin < n_enu_bins; iEnuBin++) {
      cout << htmp_rebin->GetBinContent(iEnuBin + 1) << endl;
      fout_rateonly << htmp_rebin->GetBinContent(iEnuBin + 1) << endl;
    }
    cout << endl;
    fout_rateonly << endl;

    // // for debugging
    // TCanvas * c1 = new TCanvas();
    // htmp_rebin->Draw();
    // c1->Update();
    // Char_t dummy[20];
    // cin >> dummy;
  }
  fout_rateonly.close();
}
