#include <iostream>

#include "TFile.h"
#include "TH1.h"


//Note: makes a dumb flat spectrum for fast-neutrons using the binning needed for the toy mc. When a better spectrum is ready please replace. 
void make_flat_spectrum(){
  
  //histogram for toy 
  const Int_t n_evis_toy = 240;
  const float emin=0;
  const float emax=12;
  TH1D *h_toy = new TH1D("h_toy","h_toy",n_evis_toy,emin,emax);
  for(int ibin=0;ibin<n_evis_toy;++ibin){
    double x=h_toy->GetBinCenter(ibin+1);
    if(x>0.7) h_toy->SetBinContent(ibin+1,1);
  }

  //rebinned histogram for fit 
  const Int_t n_evis_bins = 37;
  Double_t evis_bins[38]; // Single bins between 0.7 and 1.0 MeV. 0.2 MeV bins from 1.0 to 8.0 MeV. Single bin between 8.0 and 12 MeV. total 37 bins
  evis_bins[0] = 0.7;
  for (Int_t i = 0; i < 36; i++){
    evis_bins[i+1] = 0.2 *i + 1.0;
  }
  evis_bins[37] = 12.0;
  TH1D * h_fit = (TH1D*)h_toy->Rebin(n_evis_bins,"h_fit",&evis_bins[0]);
  
  //saving to file
  TFile *fout = new TFile("fn_spectrum.root","RECREATE");
  h_toy->Write();
  h_fit->Write();
  
  

}//end of make_flat_spectrum
