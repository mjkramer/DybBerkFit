#include "Config.h"
#include "Paths.h"

#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TH2F.h>
#include <TList.h>
#include <TObjArray.h>
#include <TStyle.h>
#include <TTree.h>

#include <iostream>

using namespace Config;
using namespace std;

const int nS2T = nsteps_s22t14;
const int nDM2 = nsteps_dm214_all;

void fill_thresh(TFile* f_thresh, const char* hname,
                 Double_t dchi2_thresh[nDM2][nS2T])
{
  auto h_thresh = f_thresh->Get<TH2F>(hname);
  for (int idm2 = 0; idm2 < nDM2; ++idm2) {
    for (int is2t = 0; is2t < nS2T; ++is2t) {
      dchi2_thresh[idm2][is2t] = h_thresh->GetBinContent(is2t+1, idm2+1);
      // if (dchi2_thresh[idm2][is2t] < 0) dchi2_thresh[idm2][is2t] = 0;
    }
  }
}

void fill_dchi2_rel(Double_t chi2_map[1][1][nDM2][nS2T],
                    Double_t dchi2_thresh[nDM2][nS2T],
                    TH2D* h_dchi2_rel)
{
  for (int idm2 = 0; idm2 < nDM2; ++idm2) {
    bool foundLT1 = false, foundGT1 = false;
    for (int is2t = 0; is2t < nS2T; ++is2t) {
      float dchi2_rel;
      const float thresh = dchi2_thresh[idm2][is2t];
      const float chi2 = chi2_map[0][0][idm2][is2t];
      if (thresh == -1) {         // did we skip this point?
        if (foundLT1 or foundGT1) { // right edge?
          dchi2_rel = 1001;
        } else {                // left edge
          dchi2_rel = 0;
        }
      } else {                  // thresh != -1
        dchi2_rel = chi2 / thresh;
        if (dchi2_rel < 1) foundLT1 = true;
        else foundGT1 = true;
      }
      h_dchi2_rel->SetBinContent(is2t+1, idm2+1, dchi2_rel);
    } // s2t loop
    if (not (foundLT1 and foundGT1)) {
      cout << "WARNING: Missing needed grid points for idm2 = " << idm2 << endl;
    }
  }   // dm2 loop
}

void make_data_contours_FC()
{
  gStyle->SetPalette(1);

  Double_t s22t14_bins[nS2T + 1];
  Double_t dm214_bins[nDM2 + 1];

  for (int step_dm214 = 0; step_dm214 < nsteps_dm214; ++step_dm214) {
    dm214_bins[step_dm214] =
        exp(log(dm214start) + log_dm214_step * ((Double_t)step_dm214 - 0.5));
  }

  for (int step_dm214 = 0; step_dm214 < nsteps_dm214_2 + 1; ++step_dm214) {
    dm214_bins[step_dm214 + nsteps_dm214] =
        dm214start_2 + lin_dm214_step * ((Double_t)step_dm214 - 0.5);
  }

  for (int step_s22t14 = 0; step_s22t14 < nsteps_s22t14 + 1; ++step_s22t14) {
    s22t14_bins[step_s22t14] =
        exp(log(s22t14start) + log_s22t14_step * ((Double_t)step_s22t14 - 0.5));
  }

  TH2D *h_dchi2_rel_90cl =
      new TH2D("h_dchi2_rel_90cl", "h_dchi2_rel_90cl", nS2T,
               s22t14_bins, nDM2, dm214_bins);

  TH2D *h_dchi2_rel_95cl =
      new TH2D("h_dchi2_rel_95cl", "h_dchi2_rel_95cl", nS2T,
               s22t14_bins, nDM2, dm214_bins);

  Double_t dchi2_thresh_90cl[nDM2][nS2T];
  Double_t dchi2_thresh_95cl[nDM2][nS2T];

  auto f_thresh = new TFile(Paths::outpath("thresholds_FC.root"));
  fill_thresh(f_thresh, "h90_thresh", dchi2_thresh_90cl);
  fill_thresh(f_thresh, "h95_thresh", dchi2_thresh_95cl);

  auto f_data = new TFile(Paths::outpath("fit_shape_3d.root"));
  auto tr = f_data->Get<TTree>("tr_fit");
  Double_t chi2_map[1][1][nDM2][nS2T];
  Double_t s2t_best, dm2_best, s2t14_best, dm214_best;
  tr->SetBranchAddress("chi2_map", chi2_map);
  tr->SetBranchAddress("s2t_min",&s2t_best);
  tr->SetBranchAddress("dm2_min",&dm2_best);
  tr->SetBranchAddress("dm241_min",&dm214_best);
  tr->SetBranchAddress("s2t14_min",&s2t14_best);
  tr->GetEntry(0);

  fill_dchi2_rel(chi2_map, dchi2_thresh_90cl, h_dchi2_rel_90cl);
  fill_dchi2_rel(chi2_map, dchi2_thresh_95cl, h_dchi2_rel_95cl);

  TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);
  c1->SetLeftMargin(0.15);
  c1->SetRightMargin(0.05);

  Double_t contours[3] = {-1e100, 1, 1000};
  h_dchi2_rel_90cl->SetTitleOffset(2.0, "Y");
  h_dchi2_rel_90cl->SetTitle(";sin^{2}(2#theta_{14});#Deltam^{2}_{41} (eV^{2})");
  h_dchi2_rel_90cl->SetContour(3, contours);
  h_dchi2_rel_95cl->SetContour(3, contours);
  //  hchi2_map->Draw("cont2");

  TGraph *g_best = new TGraph(1, &s2t14_best, &dm214_best);
  g_best->SetMarkerStyle(29);
  g_best->SetMarkerColor(2);
  g_best->SetMarkerSize(2);

  h_dchi2_rel_90cl->Draw("CONT LIST");
  gPad->Update();
  c1->GetPad(0)->Update();
  TObjArray *conts =
      (TObjArray *)gROOT->GetListOfSpecials()->FindObject("contours");

  TList *list_90cl = (TList *)conts->At(0)->Clone("list_90cl");

  h_dchi2_rel_95cl->Draw("CONT LIST");
  gPad->Update();
  c1->GetPad(0)->Update();
  conts = (TObjArray *)gROOT->GetListOfSpecials()->FindObject("contours");

  TList *list_95cl = (TList *)conts->At(0)->Clone("list_95cl");

  cout << list_90cl->GetSize() << " " << list_95cl->GetSize() << endl;

  //  TCanvas * c2 = new TCanvas("c2","c2",600,600);
  c1->Clear();
  c1->SetGrid();
  c1->SetLeftMargin(0.15);
  c1->SetRightMargin(0.05);
  TH1F *hh = c1->DrawFrame(s22t14start, dm214start, s22t14end, dm214end_2);
  hh->SetTitleOffset(2.0, "Y");
  hh->SetTitle(";sin^{2}(2#theta_{14});#Deltam^{2}_{41}");

  // for (Int_t i = 0; i < list_90cl->GetSize(); i++) {
  for (Int_t i = 0; i < 1; i++) {
    auto cont_90cl = (TGraph *)list_90cl->At(i);
    cont_90cl->SetLineColor(1);
    cont_90cl->SetLineWidth(2);
    cont_90cl->Draw("L");
  }
  // for (Int_t i = 0; i < list_95cl->GetSize(); i++) {
  for (Int_t i = 0; i < 1; i++) {
    auto cont_95cl = (TGraph *)list_95cl->At(i);
    cont_95cl->SetLineColor(2);
    cont_95cl->SetLineWidth(2);
    cont_95cl->Draw("L");
  }
  g_best->Draw("P");

  c1->SetLogx();
  c1->SetLogy();

  TFile *fout = new TFile(Paths::outpath("data_contours_FC.root"), "RECREATE");

  g_best->SetName("dm241_vs_sin22theta14_best");

  g_best->Write();
  // cont_90cl->Write();
  // cont_95cl->Write();
  for (Int_t i = 0; i < list_90cl->GetSize(); i++) {
    auto cont_90cl = (TGraph *)list_90cl->At(i);
    cont_90cl->SetName(Form("dm241_vs_sin22theta14_90cl_%d", i));
    cont_90cl->Write();
  }
  for (Int_t i = 0; i < list_95cl->GetSize(); i++) {
    auto cont_95cl = (TGraph *)list_95cl->At(i);
    cont_95cl->SetName(Form("dm241_vs_sin22theta14_95cl_%d", i));
    cont_95cl->Write();
  }
  // gnew_68->Write();
  // gnew_90->Write();
  h_dchi2_rel_90cl->Write();
  h_dchi2_rel_95cl->Write();
  fout->Close();
}
