void make_evis_to_enu_matrix_fine_P17B(){
  TString input_file = "../outputs/evis_to_enu_fine_2017Model_p17b.root";
  TString output_file = "matrix_evis_to_enu_fine_2017Model_P17B.txt";
  TString output_file_rateonly = "matrix_evis_to_enu_rateonly_fine_2017_Model_P17B.txt";
  const Int_t n_evis_bins = 37;
   
  Double_t evis_bins[38]; // Single bins between 0.7 and 1.0 MeV. 0.2 MeV bins from 1.0 to 8.0 MeV. Single bin between 8.0 and 12 MeV. total 37 bins
  evis_bins[0] = 0.7;
  for (Int_t i = 0; i < 36; i++){
    evis_bins[i+1] = 0.2 *i + 1.0;
  }
  evis_bins[37] = 12.0;

  const Int_t n_enu_bins = 156;
  
  Double_t enu_bins[156+1]; // 156 bins between 1.8 and 9.6 MeV
  for (Int_t i = 0; i < 156+1; i++){
    enu_bins[i] = 0.05 * i + 1.8;
  }

  
  //  TH1F * h_enu_tmp = new TH1F("h_enu_tmp","h_enu_tmp",39,enu_bins);
  
  TFile * f = new TFile(input_file.Data());
  TH2F * h_orig = (TH2F*)f->Get("h_evis_vs_enu_ad1");
  ofstream fout(output_file.Data());
  
  for (Int_t iEvisBin = 0; iEvisBin < n_evis_bins; iEvisBin++){
    //    h_enu_tmp->Reset();
    TH1D * htmp = h_orig->ProjectionX(Form("h_enu_%d",iEvisBin),iEvisBin+1,iEvisBin+1);

    TH1D * htmp_rebin = (TH1D*)htmp->Rebin(n_enu_bins,Form("h_enu_rebin_%d",iEvisBin),enu_bins);
    Double_t norm = htmp_rebin->Integral();
    htmp_rebin->Scale(1./norm);


    for (Int_t iEnuBin = 0; iEnuBin < n_enu_bins; iEnuBin++){
      cout << htmp_rebin->GetBinContent(iEnuBin+1) << endl;
      fout << htmp_rebin->GetBinContent(iEnuBin+1) << endl;
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
  for (Int_t iEvisBin = 0; iEvisBin < n_evis_bins; iEvisBin++){
    //    h_enu_tmp->Reset();
    TH1D * htmp = h_orig->ProjectionX(Form("h_enu_%d",iEvisBin));// add everything together
    
    TH1D * htmp_rebin = (TH1D*)htmp->Rebin(n_enu_bins,Form("h_enu_rebin_%d",iEvisBin),enu_bins);
    Double_t norm = htmp_rebin->Integral();
    htmp_rebin->Scale(1./norm);


    for (Int_t iEnuBin = 0; iEnuBin < n_enu_bins; iEnuBin++){
      cout << htmp_rebin->GetBinContent(iEnuBin+1) << endl;
      fout_rateonly << htmp_rebin->GetBinContent(iEnuBin+1) << endl;
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
