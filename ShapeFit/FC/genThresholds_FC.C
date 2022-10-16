#include "Config.h"
#include "Paths.h"

#include <TH2F.h>
#include <TFile.h>

#include <iostream>
#include <fstream>

using namespace Config;

void genThresholds_FC()
{
  const Int_t nq=2;
  Double_t xq[nq];
  Double_t yq[nq];

  //Define percentile that we are interested in
  xq[0]=0.90; //90th percentile
  xq[1]=0.95; //95th percentile

  const int nS2T = nsteps_s22t14;
  const int nDM2 = nsteps_dm214_all;

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

  TH2F *h95_thresh = new TH2F("h95_thresh","",nsteps_s22t14,s22t14_bins,nsteps_dm214_all,dm214_bins);
  h95_thresh->GetXaxis()->SetTitle("sin^{2}2#theta_{14}");
  h95_thresh->GetYaxis()->SetTitle("#Delta m^{2}_{41}");
  h95_thresh->GetYaxis()->SetTitleOffset(1.5);
  h95_thresh->GetXaxis()->CenterTitle();
  h95_thresh->GetYaxis()->CenterTitle();

  TH2F *h90_thresh = (TH2F*)h95_thresh->Clone("h90_thresh");

  for (int is2t = 0; is2t < nS2T; ++is2t) {
    double s2t = exp(log(s22t14start) + log_s22t14_step * is2t);
    for (int idm2 = 0; idm2 < nDM2; ++idm2) {
      double dm2;
      if (idm2 < nsteps_dm214) {
        dm2 = exp(log(dm214start) + log_dm214_step * idm2);
      } else {
        dm2 = dm214start_2 + lin_dm214_step * (idm2 - nsteps_dm214);
      }

      const char *toyconfig = "allsys_w_dm2ee_and_stat";
      auto infilename = Paths::outpath(
          "FC_fits/"
          "FC_%s_s2t13_%4.4f_dm2ee_%5.5f_s2t14_%4.4f_dm214_%5.5f.root",
          toyconfig, S22T13, DM2EE, s2t, dm2);

      const bool file_exists = std::ifstream(infilename).good();
      if (not file_exists) {
        h90_thresh->Fill(s2t, dm2, -1);
        h95_thresh->Fill(s2t, dm2, -1);
        continue;
      }

      TFile infile(infilename);

      auto h_dchi2 = infile.Get<TH1F>("h_dchi2");
      if (h_dchi2 == nullptr) {
        std::cout << "WARNING: Missing h_dchi2 in " << infilename << std::endl;
        h90_thresh->Fill(s2t, dm2, -1);
        h95_thresh->Fill(s2t, dm2, -1);
        continue;
      }
      h_dchi2->GetQuantiles(nq, yq, xq);
      h90_thresh->Fill(s2t, dm2, yq[0]);
      h95_thresh->Fill(s2t, dm2, yq[1]);
    }
  }

  TFile outfile(Paths::outpath("thresholds_FC.root"), "RECREATE");
  h90_thresh->Write();
  h95_thresh->Write();
}
