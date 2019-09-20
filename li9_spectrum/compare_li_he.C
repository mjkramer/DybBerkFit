void compare_li_he(){
    TFile* li_ihep=new TFile("./toyli9spec_IHEPmodel_v1.root","READ");
    TH1F* h_ihep_li=(TH1F*)li_ihep->Get("h_eVisAllSmeared");
    TFile* li_bcw=new TFile("./toyli9spec_BCWmodel_v1.root","READ");
    TH1F* h_bcw_li=(TH1F*)li_bcw->Get("h_eVisAllSmeared");
    TFile* li_lbnl=new TFile("./8he9li_nominal_spectrum_frac0.0.root","READ");
    TH1F* h_lbnl_li=(TH1F*)li_lbnl->Get("h_nominal");
    h_lbnl_li->Scale(500000);
    
    TFile* he_ihep=new TFile("./toyhe8spec_IHEPmodel_v1.root","READ");
    TH1F* h_ihep_he=(TH1F*)he_ihep->Get("h_eVisAllSmeared");
    TFile* he_bcw=new TFile("./toyhe8spec_BCWmodel_v1.root","READ");
    TH1F* h_bcw_he=(TH1F*)he_bcw->Get("h_eVisAllSmeared");
    TFile* he_lbnl=new TFile("./8he9li_nominal_spectrum_frac1.0.root","READ");
    TH1F* h_lbnl_he=(TH1F*)he_lbnl->Get("h_nominal");
    h_lbnl_he->Scale(500000);
    
    TCanvas* c_li=new TCanvas("c_li","c_li",500,500);
    h_ihep_li->Draw();
    h_bcw_li->Draw("SAME");
    h_bcw_li->SetLineColor(kGreen);
    h_lbnl_li->Draw("SAME");
    h_lbnl_li->SetLineColor(kBlue);
    c_li->Print("./pics/li_comparison.pdf");
    
    TCanvas* c_he=new TCanvas("c_he","c_he",500,500);
    h_ihep_he->Draw();
    h_bcw_he->Draw("SAME");
    h_bcw_he->SetLineColor(kGreen);
    h_lbnl_he->Draw("SAME");
    h_lbnl_he->SetLineColor(kBlue);
    c_he->Print("./pics/he_comparison.pdf");
}
