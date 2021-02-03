void Nominal_Li9He8_Based_on_IHEP(){
  TFile* file_li9=new TFile("toyli9spec_IHEPmodel_v1.root","READ");
  TFile* file_he8=new TFile("toyhe8spec_IHEPmodel_v1.root","READ");

  TH1F* h_li9=(TH1F*)(file_li9->Get("h_eVisAllSmeared"))->Clone("h_li9");
  TH1F* h_he8=(TH1F*)(file_he8->Get("h_eVisAllSmeared"))->Clone("h_he8");

  float frac=0.055;

  float li9norm=h_li9->Integral();
  float he8norm=h_he8->Integral();

  cout<<li9norm<<" "<<he8norm<<endl;

  h_he8->Scale(frac/(1.-frac)*li9norm/he8norm);


  cout<<h_li9->Integral()/(h_li9->Integral()+h_he8->Integral())<<" "<<h_he8->Integral()/(h_li9->Integral()+h_he8->Integral())<<endl;

  TH1F* h_combined=(TH1F*)h_li9->Clone("h_nominal");
  h_combined->SetTitle("h_nominal");
  h_combined->Reset();
  for(int ibin=0; ibin<h_combined->GetNbinsX(); ++ibin){
    h_combined->SetBinContent(ibin+1,h_li9->GetBinContent(ibin+1)+h_he8->GetBinContent(ibin+1));
    // if(h_combined->GetBinCenter(ibin+1)<0.7) h_combined->SetBinContent(ibin+1,0.);
  }

  TCanvas* canv=new TCanvas("canv","canv",500,500);
  h_combined->Draw();
  h_li9->SetLineColor(kBlue);
  h_li9->Draw("SAME");
  h_he8->SetLineColor(kGreen);
  h_he8->Draw("SAME");

  cout<<h_combined->Integral()<<endl;

  //h_combined->Scale(1./h_combined->Integral());

  TFile* outputfile=new TFile("./8he9li_nominal_spectrum_IHEP.root","RECREATE");
  // XXX for comparison against 8he9li_nominal_spectrum.root
  h_combined->Scale(0.99549488 / h_combined->Integral());
  h_combined->Write();
  outputfile->Close();

  TFile* file=new TFile("8he9li_nominal_spectrum_P15A.root","READ");
  TH1F* h_old_spectra=(TH1F*)file->Get("h_nominal");
  cout<<h_old_spectra->Integral()<<endl;
}
