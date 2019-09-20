{

  //Style
  gStyle->SetOptStat(0);
  gStyle->SetTitleOffset(1.6,"y");
  
  TH1D::SetDefaultSumw2();
  
  //Binning for fit histo
  const Int_t n_evis_bins = 37;
  Double_t evis_bins[n_evis_bins+1]; // Single bins between 0.7 and 1.0 MeV. 0.2 MeV bins from 1.0 to 8.0 MeV. Single bin between 8.0 and 12 MeV. total 37 bins
  evis_bins[0] = 0.7;
  for (Int_t i = 0; i < 36; i++){
    evis_bins[i+1] = 0.2 *i + 1.0;
  }
  evis_bins[37] = 12.0;

  //Binning for toy histo
  const Int_t n_evis_bins_toy = 240;
  Double_t evis_bins_toy[241]; // The toy spectra currently spans 0-12MeV in 240 bins
  for (Int_t i = 0; i <= n_evis_bins_toy; i++){
    evis_bins_toy[i] = 0.05 *i;
  }
  
  const Int_t n_evis_bins_fine = 113;
  Double_t evis_bins_fine[n_evis_bins_fine+1]; // Single bins between 0.7 and 1.0 MeV. 0.2 MeV bins from 1.0 to 8.0 MeV. Single bin between 8.0 and 12 MeV. total 37 bins
  for (Int_t i = 0; i < n_evis_bins_fine; i++){
    evis_bins_fine[i] = 0.1 *i + 0.7;
  }
  evis_bins_fine[n_evis_bins_fine] = 12.0;

 
  TH1F * h_sum;// = new TH1D("h_sum","h_sum",n_evis_bins_fine,evis_bins_fine);
  TH1D * h_sum_fit = new TH1D("h_sum_fit","h_sum_fit",n_evis_bins_fine,evis_bins_fine);

  TH1D * h_toy = new TH1D("h_toy","h_toy",n_evis_bins_toy,evis_bins_toy);

  //  TH1D * h_rebin = new TH1D("h_rebin","h_rebin",n_evis_bins,evis_bins);
  
  TDirectory * dir = gDirectory;
  
  Double_t livedays[3] = {229,187,236}; // number of simulated live days, accorging to Doc-6779
  TFile * f = new TFile("corrAmCSpec-DocDB-8617-v3.root");
 
  //style for h_sum
  h_sum = (TH1F*)f->Get("hCorrAmCPromptSpec");
  h_sum->SetTitle("");
  h_sum->GetXaxis()->SetTitle("E_{prompt} (MeV)");

  TCanvas * c3 = new TCanvas("c3","c3",600,600);
  h_sum->GetXaxis()->SetRangeUser(0,8);
  h_sum->Fit("expo");
  //  h_sum->Fit("pol2");
  //  h_sum->Draw();

  TF1 *ffit = (TF1*)h_sum->GetFunction("expo");

  for (Int_t i = 0; i < h_sum_fit->GetNbinsX(); i++){
    Double_t x = h_sum_fit->GetBinCenter(i+1);
    Double_t w = h_sum_fit->GetBinWidth(i+1);
    if (x < 5.5)
      h_sum_fit->SetBinContent(i+1,ffit->Eval(x) * w);
  }

  //fill toy histogram
  for (Int_t i = 0; i < h_toy->GetNbinsX(); i++){
    Double_t x = h_toy->GetBinCenter(i+1);
    Double_t w = h_toy->GetBinWidth(i+1);
    if (x>0.7 && x < 5.5)
      h_toy->SetBinContent(i+1,ffit->Eval(x) * w);
  }
  
  TCanvas * c4 = new TCanvas("c4","c4",600,600);

  //rebin for fit histo
  TH1D* h_rebin = (TH1D*)h_sum_fit->Rebin(n_evis_bins,"h_rebin",&evis_bins[0]);
  h_rebin->SetTitle("");
  h_rebin->GetXaxis()->SetTitle("E_{prompt} (MeV)");
  h_rebin->Draw();
  h_rebin->Scale(100);
  h_rebin->GetYaxis()->SetTitle("arbitrary units");

  TFile * fout = new TFile("amc_spectrum.root","recreate");
  h_sum->Write();
  h_sum_fit->Write();
  h_rebin->Write();
  h_toy->Write();
  fout->Close();
  
}
