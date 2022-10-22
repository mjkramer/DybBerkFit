#include "Config.h"
#include "Paths.h"

#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TList.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TStyle.h>
#include <TTree.h>

#include <fstream>

using namespace Config;
using namespace std;

const char* INPUTFILE = "/global/cfs/cdirs/dayabay/scratch/mkramer/sterile_fits/FC/full_v4v5v3v1_extraIWS_newE@extraIws_reallyOn_newNonUni_alphasNeutrons@bcw/fit_shape_3d_band.root";

void make_contours_CLs_band()
{
  gStyle->SetPalette(1);

  const double cls_limit = 0.05;

  const int nS2T = nsteps_s22t14;
  const int nDM2 = nsteps_dm214_all;

  Double_t chi2_map[1][1][nDM2][nS2T];
  Double_t s2t_best;
  Double_t dm2_best;

  Double_t s2t14_best;
  Double_t dm214_best;
  Double_t chi2_null;
  Double_t chi2_min;

  Double_t s22t14_bins[nS2T+1];
  Double_t dm214_bins[nDM2+1];

  for (int step_dm214=0; step_dm214<nsteps_dm214; ++step_dm214){
    dm214_bins[step_dm214] = exp(log(dm214start) + log_dm214_step*((Double_t)step_dm214-0.5));
  }

  for (int step_dm214=0; step_dm214<nsteps_dm214_2+1; ++step_dm214){
    dm214_bins[step_dm214+nsteps_dm214] = dm214start_2 + lin_dm214_step*((Double_t)step_dm214-0.5);
  }

  for (int step_s22t14=0; step_s22t14<nsteps_s22t14+1; ++step_s22t14){
    s22t14_bins[step_s22t14] = exp(log(s22t14start) + log_s22t14_step*((Double_t)step_s22t14-0.5));
  }

  Double_t dchi2_ave_3n[nDM2][nS2T];
  Double_t dchi2_ave_4n[nDM2][nS2T];

  Double_t cls_asimov[nDM2][nS2T];
  Double_t cls_map[nDM2][nS2T];

  ifstream fin_3n(Paths::outpath("asimov_dchi2_3n.txt"));
  for (int iline = 0; iline < nS2T*nDM2; ++iline) {
    int id, is;
    fin_3n >> id >> is >> dchi2_ave_3n[id][is];
  }

  ifstream fin_4n(Paths::outpath("asimov_dchi2_4n.txt"));
  for (int iline = 0; iline < nS2T*nDM2; ++iline) {
    int id, is;
    fin_4n >> id >> is >> dchi2_ave_4n[id][is];
  }

  TH2D* h_cls = new TH2D("h_cls", "h_cls", nS2T, s22t14_bins,
                          nDM2, dm214_bins);

  TFile* f_data = new TFile(INPUTFILE);

  TTree * tr = (TTree*)f_data->Get("tr_fit");
  tr->SetBranchAddress("chi2_map",chi2_map);
  tr->SetBranchAddress("s2t_min",&s2t_best);
  tr->SetBranchAddress("dm2_min",&dm2_best);
  tr->SetBranchAddress("dm241_min",&dm214_best);
  tr->SetBranchAddress("s2t14_min",&s2t14_best);
  tr->SetBranchAddress("chi2_null",&chi2_null);
  tr->SetBranchAddress("chi2_min",&chi2_min);
  // tr->GetEntry(0);

  TCanvas * c1 = new TCanvas("c1","c1",600,600);

  const char* outfilename = "band_contours.root";
  TFile * fout = new TFile(Paths::outpath(outfilename),
                          "RECREATE");


  const int nfits = tr->GetEntries();
  for (int entry = 0; entry < nfits; ++entry) {
    tr->GetEntry(entry);
    h_cls->Reset();

    for (Int_t iDM2 = 0; iDM2 < nDM2;  iDM2 ++) {
      for (Int_t iS2T = 0; iS2T < nS2T;  iS2T ++) {
        // Double_t dchi2_obs = chi2_map[0][0][iDM2][iS2T] + chi2_min;
        Double_t dchi2_obs = chi2_map[0][0][iDM2][iS2T] - chi2_null+chi2_min;
        //dchi2_obs = dchi2_ave_3n[iDM2][iS2T]; // in case of expected

        Double_t p0 = 0.5 + 0.5 * TMath::Erf((dchi2_obs - dchi2_ave_3n[iDM2][iS2T])
                                            /TMath::Sqrt(8*TMath::Abs(dchi2_ave_3n[iDM2][iS2T])));

        Double_t p1 = 0.5 + 0.5 * TMath::Erf((dchi2_obs + dchi2_ave_4n[iDM2][iS2T])
                                            /TMath::Sqrt(8*TMath::Abs(dchi2_ave_4n[iDM2][iS2T])));
        if (1-p0 > 0){
          cls_map[iDM2][iS2T] = (1-p1)/(1-p0);
        }else{
          cls_map[iDM2][iS2T] = 0;
        }

        h_cls->SetBinContent(iS2T+1, iDM2+1, cls_map[iDM2][iS2T]);
      }
    }

    c1->SetLeftMargin(0.15);
    c1->SetRightMargin(0.05);

    Double_t contours[3] = {-1000,cls_limit,1000};
    h_cls->SetTitleOffset(2.0,"Y");
    h_cls->SetTitle(";sin^{2}(2#theta_{14});#Deltam^{2}_{41} (eV^{2})");
    h_cls->SetContour(3,contours);

    TGraph * g_best = new TGraph(1, &s2t14_best, &dm214_best);
    g_best->SetMarkerStyle(29);
    g_best->SetMarkerColor(2);
    g_best->SetMarkerSize(2);

    h_cls->Draw("CONT LIST");
    gPad->Update();
    c1->GetPad(0)->Update();
    TObjArray *conts = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");

    TList *list_68 = (TList*)conts->At(0)->Clone("list_68");

    //  TCanvas * c2 = new TCanvas("c2","c2",600,600);
    c1->Clear();
    c1->SetGrid();
    c1->SetLeftMargin(0.15);
    c1->SetRightMargin(0.05);
    TH1F * hh = c1->DrawFrame(s22t14start,dm214start,s22t14end,dm214end_2);
    hh->SetTitleOffset(2.0,"Y");
    hh->SetTitle(";sin^{2}(2#theta_{14});#Deltam^{2}_{41}");

    // gPad->SetLogy();
    // gPad->SetLogx();

    for (Int_t i = 0; i < list_68->GetSize(); i++){
      TGraph* cont_68 = (TGraph*)list_68->At(i);
      cont_68->SetLineColor(1);
      cont_68->SetLineWidth(2);
      cont_68->Draw("L");
    }
    g_best->Draw("P");

    g_best->SetName("dm241_vs_sin22theta14_best");

    // g_best->Write();
    // cont_68->Write();
    // cont_90->Write();
    for (Int_t i = 0; i < list_68->GetSize(); i++){
      TGraph * cont_68 = (TGraph*)list_68->At(i);
      cont_68->SetName(Form("dm241_vs_sin22theta14_cls_95cl_%d_%d",i, entry));
      cont_68->Write();
    }
  }

  // h_cls->Write();

  double x[nsteps_dm214_all];
  double y[nsteps_dm214_all];
  double limits[nfits];
  for (int iy = 0; iy < nsteps_dm214_all; ++iy) {
    y[iy] = h_cls->GetYaxis()->GetBinCenter(iy + 1);
    // y[iy] = 0.5 * (dm214_bins[iy] + dm214_bins[iy+1]);
    for (int ifit = 0; ifit < nfits; ++ifit) {
      auto g0 = fout->Get<TGraph>(Form("dm241_vs_sin22theta14_cls_95cl_0_%d", ifit));
      // std::cout << g0 << std::endl;
      TGraph g_inv(g0->GetN(), g0->GetY(), g0->GetX());
      // XXX
      // limits[ifit] = g_inv.Eval(y[iy], nullptr, "S");
      limits[ifit] = g_inv.Eval(y[iy]);
    }
    // XXX
    // x[iy] = TMath::Median(nfits, limits);
    std::sort(limits, limits + nfits);
    if (nfits % 2)
      x[iy] = limits[nfits/2];
    else
      x[iy] = 0.5 * (limits[nfits/2 - 1] + limits[nfits/2]);
  }

  auto g = new TGraph(nsteps_dm214_all, x, y);
  g->SetName("dm241_vs_sin22theta14_cls_95cl_median");
  g->Write();

  fout->Close();
}
