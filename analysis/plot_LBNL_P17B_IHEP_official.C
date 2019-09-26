{
  gStyle->SetPalette(1);
  gStyle->SetLabelFont(62, "XY");
  gStyle->SetTitleFont(62, "XY");
  gStyle->SetLegendFont(62);
  gROOT->ForceStyle();
  const Int_t nS2T = 101;       // IHEP-official
  const Int_t nDM2 = 101;       // IHEP-official
  // const Int_t nS2T = 31;        // mine
  // const Int_t nDM2 = 31;        // mine

  Double_t chi2_map[nDM2][nS2T];
  Double_t chi2_s2t[nS2T];
  Double_t chi2_dm2[nDM2];
  Double_t s2t_best;
  Double_t dm2_best;
  Double_t chi2_best = 1e10;
  Double_t chi2_min_from_tree = 1e10;

  Double_t dm_corr = 0.; // 5.17e-5;

  Double_t s2t_min = 0.06;      // IHEP-official
  Double_t s2t_max = 0.11;      // IHEP-official
  Double_t dm2_min = 2.1e-3 - dm_corr; // IHEP-official
  Double_t dm2_max = 2.9e-3 - dm_corr; // IHEP-official
  // Double_t s2t_min = 0.07;      // mine
  // Double_t s2t_max = 0.10;      // mine
  // Double_t dm2_min = 2.1e-3 - dm_corr; // mine
  // Double_t dm2_max = 2.8e-3 - dm_corr; // mine

  Double_t s2t_step = (s2t_max - s2t_min) / (nS2T - 1);
  Double_t dm2_step = (dm2_max - dm2_min) / (nDM2 - 1);

  // TFile* bestfit_file = new TFile("../fit_shape_2d_2017Model_P17B_LBNL_BCWbin.root","OPEN");
  // TFile* bestfit_file = new TFile("../fit_shape_2d_2017Model_P17B_LBNL.root","OPEN");
  // TFile* bestfit_file = new TFile("../fit_shape_2d_SCNLModel_P17B_LBNL_BCWbin.root","OPEN");
  // TFile* bestfit_file = new TFile("../fit_shape_2d_SCNLModel_P17B_LBNL.root","OPEN");
  // TFile* bestfit_file = new TFile("../fit_shape_2d_2017Model_P17B_BCW_BCWbin.root","OPEN");
  // TFile* bestfit_file = new TFile("../fit_shape_2d_2017Model_postP15A_LBNL_BCWbin.root","OPEN");
  // TFile* bestfit_file = new TFile("../fit_shape_2d_2017Model_P17B_LBNL_0.2_inflated.root","OPEN");
  // TFile* bestfit_file = new TFile("../fit_shape_2d_2017Model_P17B_IHEP_BCWbin_inverted.root","OPEN");
  // TFile* bestfit_file = new TFile("../fit_shape_2d_2017Model_P17B_IHEP_BCWbin_statonly.root","OPEN");

  // TFile *bestfit_file = new TFile("../ShapeFit/fit_result_files/" "fit_shape_2d_2017Model_P17B_LBNL_0.1_inflated.root", "OPEN");
  // TFile* bestfit_file = new TFile("/global/project/projectdirs/dayabay/scratch/beda/Theta13Analysis2017_Based_on_2016/ShapeFit/fit_result_files/fit_shape_2d_2017Model_P17B_LBNL_0.1_inflated.root","OPEN");
  // TFile* bestfit_file = new TFile("/global/project/projectdirs/dayabay/scratch/beda/Theta13Analysis2017_Based_on_2016/ShapeFit/fit_result_files/fit_shape_2d_2017Model_P17B_LBNL.root","OPEN");
  TFile* bestfit_file = new TFile("/project/projectdirs/dayabay/scratch/beda/Theta13Analysis2017_Based_on_2016/ShapeFit/fit_result_files/fit_shape_2d_2017Model_P17B_IHEP_BCWbin.root","OPEN"); // official
  

  TFile *fout = new TFile("./LBNL_contour_P17B.root", "recreate");

  TH2D *hchi2_map =
      new TH2D("hchi2_map", "hchi2_map", nS2T, s2t_min - 0.5 * s2t_step,
               s2t_max + 0.5 * s2t_step, nDM2, dm2_min - 0.5 * dm2_step,
               dm2_max + 0.5 * dm2_step);

  TH1D *hchi2_s2t =
      new TH1D("hchi2_s2t", "hchi2_s2t", nS2T, s2t_min - 0.5 * s2t_step,
               s2t_max + 0.5 * s2t_step);

  TH1D *hchi2_dm2 =
      new TH1D("hchi2_dm2", "hchi2_dm2", nDM2, dm2_min - 0.5 * dm2_step,
               dm2_max + 0.5 * dm2_step);

  Double_t sin22t13[nS2T];
  Double_t dm2[nDM2];

  for (int step = 0; step < nS2T; ++step) {
    sin22t13[step] = (s2t_max - s2t_min) * step * 1. / (nS2T - 1) + s2t_min;
  }
  for (int step_dm2 = 0; step_dm2 < nDM2; ++step_dm2) {
    dm2[step_dm2] =
        (dm2_max - dm2_min) * step_dm2 * 1. / (nDM2 - 1) + dm2_min; //-dm_corr;
  }

  //**********READ decoherence_result.txt *****************//
  // Best Fit
  TTree *tr_full = (TTree *)bestfit_file->Get("tr_fit");
  tr_full->SetBranchAddress("s2t_min", &s2t_best);
  tr_full->SetBranchAddress("dm2_min", &dm2_best);
  tr_full->SetBranchAddress("chi2_min", &chi2_min_from_tree);
  // dm2_best -= dm_corr;
  tr_full->GetEntry(0);
  cout << s2t_best << endl;
  cout << "Real chi2 min from tr_fit is " << chi2_min_from_tree << endl;

  TH2F *h_chi2_map = (TH2F *)bestfit_file->Get("h_chi2_map");
  // ifstream input_file("./LBNL_6+8AD.txt");

  for (Int_t id = 0; id < nDM2; id++) {
    for (Int_t is = 0; is < nS2T; is++) {

      // input_file >> id >> is >> chi2_map[id][is];

      chi2_map[id][is] = h_chi2_map->GetBinContent(is + 1, id + 1);
      // Find smallest chi2 in map
      if (chi2_map[id][is] < chi2_best) {
        chi2_best = chi2_map[id][is];
        // s2t_best = s2t_min + is * s2t_step;
        // dm2_best = dm2_min + id * dm2_step;
      }
    }
  }

  cout << chi2_best << " but not using this" << endl;

  chi2_best = 0.; // since it is already corrected from file

  for (Int_t id = 0; id < nDM2; id++) {
    for (Int_t is = 0; is < nS2T; is++) {
      hchi2_map->SetBinContent(is + 1, id + 1, chi2_map[id][is] - chi2_best);
    }
  }
  hchi2_map->Write();

  // Now get best projection to remove theta13 dependency
  Double_t projected_chi2;
  Int_t is_best = 0;
  for (Int_t id = 0; id < nDM2; id++) {

    projected_chi2 = 1e10; // Reset to absurdly high value

    // Scan over theta13 value
    for (Int_t is = 0; is < nS2T; is++) {

      if (chi2_map[id][is] < projected_chi2) {
        projected_chi2 = chi2_map[id][is];
        is_best = is;
      }

    } // t13

    // Now save best projected value into TH1D
    chi2_dm2[id] = chi2_map[id][is_best] - chi2_best;
    hchi2_dm2->SetBinContent(id + 1, chi2_dm2[id]);

  } // dm2

  hchi2_dm2->Write();

  Int_t id_best = 0;
  for (Int_t is = 0; is < nS2T; is++) {

    projected_chi2 = 1e10; // Reset to absurdly high value

    // Scan over dm2 value
    for (Int_t id = 0; id < nDM2; id++) {

      if (chi2_map[id][is] < projected_chi2) {
        projected_chi2 = chi2_map[id][is];
        id_best = id;
      }

    } // dm2

    // Now save best projected value into TH1D
    chi2_s2t[is] = chi2_map[id_best][is] - chi2_best;
    hchi2_s2t->SetBinContent(is + 1, chi2_s2t[is]);

  } // t13

  hchi2_s2t->Write();

  /*


  Char_t output_filename[256];
  //  sprintf(output_filename, "contours_2d_%s_toyMC.root",opt);
  //  sprintf(output_filename, "contours_2d_toyMC.root");
  //  sprintf(output_filename, "contours_2d_datafit_p12b.root");
  sprintf(output_filename,
  "contours_2d_datafit_physics_NL2014_unblinded_6AD_8AD.root");


  TTree * tr_full = (TTree*)f_full->Get("tr_fit");
  tr_full->SetBranchAddress("chi2_map",chi2_map);
  tr_full->SetBranchAddress("s2t_min",&s2t_best);
  tr_full->SetBranchAddress("dm2_min",&dm2_best);
  tr_full->SetBranchAddress("chi2_min",&dm2_best);
  tr_full->GetEntry(0);
  chi2_best = 1e8;
  for (Int_t iDM2 = 0; iDM2 < nDM2;  iDM2 ++){
    for (Int_t iS2T = 0; iS2T < nS2T;  iS2T ++){
      hchi2_map->SetBinContent(iS2T+1, iDM2+1, chi2_map[iDM2][iS2T]);

      if(chi2_best > chi2_map[iDM2][iS2T])
        chi2_best = chi2_map[iDM2][iS2T];

    }
  }
  */

  TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);
  c1->SetLeftMargin(0.15);
  c1->SetRightMargin(0.05);
  // Best fit, 1sigma, 90%, 99%
  // Double_t contours[4] = {-1, 2.30, 4.61,9.21};

  // Best fit, 1sigma, 2sigma,3sigma
  Double_t contours[4] = {-1, 2.30, 6.18, 11.83};

  hchi2_map->SetTitleOffset(2.0, "Y");
  hchi2_map->SetTitle(";sin^{2}(2#theta_{13});#Deltam^{2}_{ee}");
  hchi2_map->SetContour(4, contours);
  //  hchi2_map->Draw("cont2");

  // dm2_best = 1000*(dm2_best);

  cout << "s2t_best: " << s2t_best << endl;
  cout << "dm2_best: " << dm2_best << endl;
  cout << "chi2_best: " << chi2_best << endl;

  TGraph *g_best = new TGraph(1, &s2t_best, &dm2_best);
  g_best->SetMarkerStyle(29);
  g_best->SetMarkerColor(3);
  g_best->SetMarkerSize(1);

  hchi2_map->Draw("CONT LIST");
  gPad->Update();
  c1->GetPad(0)->Update();
  TObjArray *conts =
      (TObjArray *)gROOT->GetListOfSpecials()->FindObject("contours");
  // TObjArray *conts = c1->GetPad(0)->FindObject("contours");
  //  cout << conts << endl;
  Int_t ncontours = conts->GetSize();
  cout << ncontours << endl;

  TList *list_68 = (TList *)conts->At(0);
  TList *list_90 = (TList *)conts->At(1);
  TList *list_99 = (TList *)conts->At(2);

  cout << list_68->GetSize() << " " << list_90->GetSize() << " "
       << list_99->GetSize() << endl;

  TGraph *cont_68 = (TGraph *)list_68->First();
  TGraph *cont_90 = (TGraph *)list_90->First();
  TGraph *cont_99 = (TGraph *)list_99->First();

  Int_t n_68 = cont_68->GetN();
  Double_t *s2t_68 = cont_68->GetX();
  Double_t *dm2_68 = cont_68->GetY();

  Int_t n_90 = cont_90->GetN();
  Double_t *s2t_90 = cont_90->GetX();
  Double_t *dm2_90 = cont_90->GetY();

  Int_t n_99 = cont_99->GetN();
  Double_t *s2t_99 = cont_99->GetX();
  Double_t *dm2_99 = cont_99->GetY();

  Double_t tan2_68[10000];
  Double_t tan2_90[10000];
  Double_t tan2_99[10000];

  for (Int_t i = 0; i < n_68; i++) {
    Double_t cos2t = TMath::Sqrt(1 - s2t_68[i]);
    tan2_68[i] = (1 - cos2t) / (1 + cos2t);
    //    cout << " " << i << " " << s2t_68[i] << " " << tan2_68[i] << endl;
  }
  for (Int_t i = 0; i < n_90; i++) {
    Double_t cos2t = TMath::Sqrt(1 - s2t_90[i]);
    tan2_90[i] = (1 - cos2t) / (1 + cos2t);
    //    cout << " " << i << " " << s2t_68[i] << " " << tan2_68[i] << endl;
  }
  for (Int_t i = 0; i < n_99; i++) {
    Double_t cos2t = TMath::Sqrt(1 - s2t_99[i]);
    tan2_99[i] = (1 - cos2t) / (1 + cos2t);
    //    cout << " " << i << " " << s2t_68[i] << " " << tan2_68[i] << endl;
  }

  TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
  c2->SetGrid();
  c2->SetLeftMargin(0.15);
  c2->SetRightMargin(0.05);
  // TH1F * hh = c2->DrawFrame(0.00,0.0015,0.20,0.0035);
  TH1F *hh = c2->DrawFrame(0.06, 0.0021, 0.11, 0.0029);
  hh->SetTitleOffset(2.0, "Y");
  // hh->SetTitle(";sin^{2}(2#theta_{13});#Deltam^{2}_{ee} (10^{-3} eV^{2})");
  hh->SetTitle(";sin^{2}(2#theta_{13});#Deltam^{2}_{ee} [eV^{2}])");
  // hh->SetTitle(";sin^{2}(2#theta_{13});#sigma_{rel}");

  for (Int_t i = 0; i < list_68->GetSize(); i++) {
    cont_68 = (TGraph *)list_68->At(i);
    cont_68->SetLineColor(4);
    cont_68->SetLineStyle(1);
    cont_68->SetLineWidth(2);
    cont_68->Draw("L");
  }
  for (Int_t i = 0; i < list_90->GetSize(); i++) {
    cont_90 = (TGraph *)list_90->At(i);
    cont_90->SetLineColor(1);
    cont_90->SetLineStyle(1);
    cont_90->SetLineWidth(2);
    cont_90->Draw("L");
  }
  for (Int_t i = 0; i < list_99->GetSize(); i++) {
    cont_99 = (TGraph *)list_99->At(i);
    cont_99->SetLineColor(2);
    cont_99->SetLineStyle(1);
    cont_99->SetLineWidth(2);
    cont_99->Draw("L");
  }
  g_best->Draw("P");
  /*
  TFile* file_p15a=new TFile("./LBNL_contour.root","READ");
  c2->cd();
  TGraph* gr_best_p15a=(TGraph*)file_p15a->Get("dm2_ee_vs_s2t13_best");
  gr_best_p15a->Draw("SAME P");
  TGraph* gr_68=(TGraph*)file_p15a->Get("dm2_ee_vs_s2t13_cont68");
  gr_68->Draw("SAME");
  TGraph* gr_95=(TGraph*)file_p15a->Get("dm2_ee_vs_s2t13_cont95");
  gr_95->Draw("SAME");
  TGraph* gr_99=(TGraph*)file_p15a->Get("dm2_ee_vs_s2t13_cont99");
  gr_99->Draw("SAME");


  c2->Print("LBNL_contour_s2t_dm2ee.pdf");*/

  TCanvas *c3 = new TCanvas("c3", "c3", 600, 600);
  c3->SetLeftMargin(0.15);
  c3->SetRightMargin(0.05);
  TH1F *hhh = c3->DrawFrame(0, 0, 0.09, 0.0055);
  hhh->SetTitleOffset(2.0, "Y");
  hhh->SetTitle(";tan^{2}(#theta_{13});#Deltam^{2}_{ee}");
  TGraph *gnew_68 = new TGraph(n_68, tan2_68, dm2_68);
  gnew_68->SetLineColor(2);
  gnew_68->SetLineWidth(2);
  gnew_68->Draw("L");
  TGraph *gnew_90 = new TGraph(n_90, tan2_90, dm2_90);
  gnew_90->SetLineColor(4);
  gnew_90->SetLineWidth(2);
  gnew_90->Draw("L");
  TGraph *gnew_99 = new TGraph(n_99, tan2_99, dm2_99);
  gnew_99->SetLineColor(4);
  gnew_99->SetLineWidth(2);
  gnew_99->Draw("L");

  // TFile * fout = new TFile(output_filename,"recreate");

  g_best->SetName("dm2_ee_vs_s2t13_best");
  cont_68->SetName("dm2_ee_vs_s2t13_cont68");
  cont_90->SetName("dm2_ee_vs_s2t13_cont95");
  cont_99->SetName("dm2_ee_vs_s2t13_cont99");
  // gnew_68->SetName("dm2_vs_tan2theta13_cont68");
  // gnew_90->SetName("dm2_vs_tan2theta13_cont90");
  // gnew_99->SetName("dm2_vs_tan2theta13_cont99");

  g_best->Write();
  cont_68->Write();
  cont_90->Write();
  cont_99->Write();
  // gnew_68->Write();
  // gnew_90->Write();
  // gnew_99->Write();

  c2->SaveAs("LBNL_contour_P17B.pdf");

  const Int_t maxthre = 3;

  Double_t dchi2_thre[maxthre] = {1, 4, 9};

  Double_t allowed_min_s[maxthre];
  Double_t allowed_max_s[maxthre];

  Double_t allowed_min_d[maxthre];
  Double_t allowed_max_d[maxthre];

  for (Int_t i = 0; i < nS2T - 1; i++) {
    for (Int_t ithre = 0; ithre < maxthre; ithre++) {
      if (chi2_s2t[i] > dchi2_thre[ithre] &&
          chi2_s2t[i + 1] < dchi2_thre[ithre]) {
        allowed_min_s[ithre] =
            ((dchi2_thre[ithre] - chi2_s2t[i + 1]) * sin22t13[i] +
             (chi2_s2t[i] - dchi2_thre[ithre]) * sin22t13[i + 1]) /
            (chi2_s2t[i] - chi2_s2t[i + 1]);
      }
      if (chi2_s2t[i] < dchi2_thre[ithre] &&
          chi2_s2t[i + 1] > dchi2_thre[ithre]) {
        allowed_max_s[ithre] =
            ((dchi2_thre[ithre] - chi2_s2t[i + 1]) * sin22t13[i] +
             (chi2_s2t[i] - dchi2_thre[ithre]) * sin22t13[i + 1]) /
            (chi2_s2t[i] - chi2_s2t[i + 1]);
      }
    }
  }

  for (Int_t i = 0; i < nDM2 - 1; i++) {
    for (Int_t ithre = 0; ithre < maxthre; ithre++) {
      if (chi2_dm2[i] > dchi2_thre[ithre] &&
          chi2_dm2[i + 1] < dchi2_thre[ithre]) {
        allowed_min_d[ithre] =
            ((dchi2_thre[ithre] - chi2_dm2[i + 1]) * dm2[i] +
             (chi2_dm2[i] - dchi2_thre[ithre]) * dm2[i + 1]) /
            (chi2_dm2[i] - chi2_dm2[i + 1]);
      }
      if (chi2_dm2[i] < dchi2_thre[ithre] &&
          chi2_dm2[i + 1] > dchi2_thre[ithre]) {
        allowed_max_d[ithre] =
            ((dchi2_thre[ithre] - chi2_dm2[i + 1]) * dm2[i] +
             (chi2_dm2[i] - dchi2_thre[ithre]) * dm2[i + 1]) /
            (chi2_dm2[i] - chi2_dm2[i + 1]);
      }
    }
  }

  cout << " sin^2(2theta_13) = " << s2t_best << " - "
       << s2t_best - allowed_min_s[0] << " + " << allowed_max_s[0] - s2t_best
       << endl;

  cout << " delta m^2_ee or delta m^2_32 check carefully! = "
       << dm2_best /*-5.17e-5*/ // - dm_corr
       << " - " << dm2_best - allowed_min_d[0] << " + "
       << allowed_max_d[0] - dm2_best << endl;

  cout << " chi2 = " << chi2_best << endl;

  fout->Close();

  TGraph *cont68_better = new TGraph();
  double x, y;
  for (int ipoint = 0; ipoint < cont_68->GetN(); ++ipoint) {
    cont_68->GetPoint(ipoint, x, y);
    cont68_better->SetPoint(ipoint, x, y * 1000);
  }
  TGraph *cont90_better = new TGraph();
  for (int ipoint = 0; ipoint < cont_90->GetN(); ++ipoint) {
    cont_90->GetPoint(ipoint, x, y);
    cont90_better->SetPoint(ipoint, x, y * 1000);
  }
  TGraph *cont99_better = new TGraph();
  for (int ipoint = 0; ipoint < cont_99->GetN(); ++ipoint) {
    cont_99->GetPoint(ipoint, x, y);
    cont99_better->SetPoint(ipoint, x, y * 1000);
  }

  TGraphAsymmErrors *bestfit_better = new TGraphAsymmErrors();
  bestfit_better->SetPoint(0, s2t_best, dm2_best * 1000);
  bestfit_better->SetPointError(0, s2t_best - allowed_min_s[0],
                                allowed_max_s[0] - s2t_best,
                                dm2_best * 1000 - allowed_min_d[0] * 1000,
                                allowed_max_d[0] * 1000 - dm2_best * 1000);

  TGraph *gr_dm2_profile = new TGraph();
  for (ipoint = 0; ipoint < nDM2; ++ipoint) {
    gr_dm2_profile->SetPoint(ipoint, chi2_dm2[ipoint], dm2[ipoint] * 1000);
    cout << chi2_dm2[ipoint] << " " << dm2[ipoint] * 1000 << endl;
  }
  gr_dm2_profile->SetTitle("");
  gr_dm2_profile->SetLineColor(kBlue);
  gr_dm2_profile->SetLineWidth(2);

  TGraph *gr_s2t_profile = new TGraph(nS2T, sin22t13, chi2_s2t);
  gr_s2t_profile->SetTitle("");
  gr_s2t_profile->SetLineColor(kBlue);
  gr_s2t_profile->SetLineWidth(2);

  TCanvas *c_limits = new TCanvas("c_limits", "c_limits", 900, 900);
  c_limits->Divide(2, 2);
  c_limits->cd(3); // contours
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.14);
  TH1F *small_frame =
      gPad->DrawFrame(0.065, 0.00215 * 1000, 0.105, 0.00285 * 1000);

  small_frame->SetTitle(
      ";sin^{2}(2#theta_{13});#Deltam^{2}_{ee} [#times10^{-3} eV^{2}]");
  small_frame->GetXaxis()->CenterTitle();
  small_frame->GetYaxis()->CenterTitle();
  small_frame->GetXaxis()->SetTitleOffset(1.);
  small_frame->GetYaxis()->SetTitleOffset(1.19);
  small_frame->GetXaxis()->SetTitleSize(0.06);
  small_frame->GetYaxis()->SetTitleSize(0.06);
  small_frame->GetXaxis()->SetNdivisions(505);
  small_frame->GetXaxis()->SetLabelSize(0.042);
  small_frame->GetYaxis()->SetLabelSize(0.042);
  gPad->SetPad(0., 0., 0.7, 0.7);

  gPad->SetTopMargin(0.);
  gPad->SetRightMargin(0.);
  // cont99_better->SetFillStyle(3001);
  // cont90_better->SetFillStyle(3001);
  // cont68_better->SetFillStyle(3001);
  cont68_better->SetFillColor(46);
  cont90_better->SetFillColor(8);
  cont99_better->Draw("SAME F");
  cont99_better->SetFillColor(9);
  cont90_better->Draw("SAME F");
  cont68_better->Draw("SAME F");
  gPad->SetGrid(1);

  small_frame->Draw("sameaxis");
  gPad->RedrawAxis("g");

  // gStyle->SetHatchesSpacing(100.);

  bestfit_better->Draw("SAME P");

  c_limits->cd(4); // blank
  gPad->SetPad(0.9, 0.9, 1., 1.);

  c_limits->cd(2); // s2t profile
  gPad->SetPad(0., 0.7, 1., 1.);
  gPad->SetBottomMargin(0.);
  gPad->SetRightMargin(0.3);
  gPad->SetLeftMargin(0.105);

  TH1F *h2 = new TH1F("h2", "h2", 100, 0.065, 0.105);
  h2->GetYaxis()->SetRangeUser(0., 10.);
  h2->Draw("Y+");
  h2->SetTitle("");
  h2->GetYaxis()->SetNdivisions(505);

  gr_s2t_profile->GetYaxis()->SetRangeUser(0., 10.);
  gr_s2t_profile->GetYaxis()->SetTitle("#Delta#chi^{2}");
  gr_s2t_profile->Draw("SAME L");
  gr_s2t_profile->GetYaxis()->SetNdivisions(505);

  TLine *l1_2 = new TLine(allowed_min_s[0], 0., allowed_min_s[0], 1.);
  TLine *l2_2 = new TLine(allowed_max_s[0], 0., allowed_max_s[0], 1.);
  TLine *l3_2 = new TLine(allowed_min_s[0], 1., allowed_max_s[0], 1.);
  l1_2->SetLineStyle(2);
  l2_2->SetLineStyle(2);
  l3_2->SetLineStyle(2);
  l1_2->SetLineColor(kRed);
  l2_2->SetLineColor(kRed);
  l3_2->SetLineColor(kRed);
  l1_2->Draw("SAME");
  l2_2->Draw("SAME");
  l3_2->Draw("SAME");

  TLine *l1_3 = new TLine(allowed_min_s[1], 0., allowed_min_s[1], 4.);
  TLine *l2_3 = new TLine(allowed_max_s[1], 0., allowed_max_s[1], 4.);
  TLine *l3_3 = new TLine(allowed_min_s[1], 4., allowed_max_s[1], 4.);
  l1_3->SetLineStyle(2);
  l2_3->SetLineStyle(2);
  l3_3->SetLineStyle(2);
  l1_3->SetLineColor(kRed);
  l2_3->SetLineColor(kRed);
  l3_3->SetLineColor(kRed);
  l1_3->Draw("SAME");
  l2_3->Draw("SAME");
  l3_3->Draw("SAME");

  TLine *l1 = new TLine(allowed_min_s[2], 0., allowed_min_s[2], 9.);
  TLine *l2 = new TLine(allowed_max_s[2], 0., allowed_max_s[2], 9.);
  TLine *l3 = new TLine(allowed_min_s[2], 9., allowed_max_s[2], 9.);
  l1->SetLineStyle(2);
  l2->SetLineStyle(2);
  l3->SetLineStyle(2);
  l1->SetLineColor(kRed);
  l2->SetLineColor(kRed);
  l3->SetLineColor(kRed);
  l1->Draw("SAME");
  l2->Draw("SAME");
  l3->Draw("SAME");

  gPad->SetGrid(1);

  c_limits->cd(4); // dm2 profile
  gPad->SetPad(0.7, 0., 1., 1.0);
  gPad->SetTopMargin(0.3);
  gPad->SetBottomMargin(0.098);
  gPad->SetGrid(1);
  gPad->SetLeftMargin(0);
  // gPad->SetTopMargin(0);
  TH1F *h = new TH1F("h", "h", 1, 0.001, 10.);
  h->GetYaxis()->SetRangeUser(2.15, 2.85);
  h->SetBinContent(1, 0.5);
  h->Draw("X+");

  h->GetXaxis()->SetNdivisions(505);
  h->SetTitle("");
  h->SetStats(0);
  h->GetXaxis()->SetLabelOffset(0.007);
  h->GetXaxis()->SetLabelSize(0.07);
  h->GetXaxis()->SetTitle("#Delta#chi^{2}");
  h->GetXaxis()->SetTitleSize(0.08);
  h->GetXaxis()->SetLabelFont(62);
  h->GetXaxis()->SetTitleOffset(0.4);
  gr_dm2_profile->Draw("SAME L");
  // gr_dm2_profile->GetXaxis()->SetRangeUser(0.015,10.);
  // gr_dm2_profile->GetXaxis()->SetTitle("#Delta#chi^{2}");

  TGaxis *ax = new TGaxis(0., 2.85, 0., 3.157, 0.001, 10., 505, "+L");
  ax->Draw("SAME");
  ax->SetLabelOffset(0.03);
  ax->SetTitle("#Delta#chi^{2}");
  ax->SetTitleOffset(1.5);
  ax->SetLabelSize(0.07);
  ax->SetTitleSize(0.08);
  // ax->SetTitleFont(62);

  TLine *ll1_2 =
      new TLine(0., allowed_min_d[0] * 1000, 1., allowed_min_d[0] * 1000);
  TLine *ll2_2 =
      new TLine(0., allowed_max_d[0] * 1000, 1., allowed_max_d[0] * 1000);
  TLine *ll3_2 =
      new TLine(1., allowed_min_d[0] * 1000, 1., allowed_max_d[0] * 1000);
  ll1_2->SetLineStyle(2);
  ll2_2->SetLineStyle(2);
  ll3_2->SetLineStyle(2);
  ll1_2->SetLineColor(kRed);
  ll2_2->SetLineColor(kRed);
  ll3_2->SetLineColor(kRed);
  ll1_2->Draw("SAME");
  ll2_2->Draw("SAME");
  ll3_2->Draw("SAME");

  TLine *ll1_1 =
      new TLine(0., allowed_min_d[1] * 1000, 4., allowed_min_d[1] * 1000);
  TLine *ll2_1 =
      new TLine(0., allowed_max_d[1] * 1000, 4., allowed_max_d[1] * 1000);
  TLine *ll3_1 =
      new TLine(4., allowed_min_d[1] * 1000, 4., allowed_max_d[1] * 1000);
  ll1_1->SetLineStyle(2);
  ll2_1->SetLineStyle(2);
  ll3_1->SetLineStyle(2);
  ll1_1->SetLineColor(kRed);
  ll2_1->SetLineColor(kRed);
  ll3_1->SetLineColor(kRed);
  ll1_1->Draw("SAME");
  ll2_1->Draw("SAME");
  ll3_1->Draw("SAME");

  TLine *ll1_3 =
      new TLine(0., allowed_min_d[2] * 1000, 9., allowed_min_d[2] * 1000);
  TLine *ll2_3 =
      new TLine(0., allowed_max_d[2] * 1000, 9., allowed_max_d[2] * 1000);
  TLine *ll3_3 =
      new TLine(9., allowed_min_d[2] * 1000, 9., allowed_max_d[2] * 1000);
  ll1_3->SetLineStyle(2);
  ll2_3->SetLineStyle(2);
  ll3_3->SetLineStyle(2);
  ll1_3->SetLineColor(kRed);
  ll2_3->SetLineColor(kRed);
  ll3_3->SetLineColor(kRed);
  ll1_3->Draw("SAME");
  ll2_3->Draw("SAME");
  ll3_3->Draw("SAME");

  c_limits->Print("./pics/limits.pdf");
}
