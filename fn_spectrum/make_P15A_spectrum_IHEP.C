#include <iostream>

#include "TFile.h"
#include "TH1.h"


//Note: makes a dumb flat spectrum for fast-neutrons using the binning needed for the toy mc. When a better spectrum is ready please replace. 
void make_P15A_spectrum_IHEP(){
  
  Double_t E0[8] = 
    {68.68,
     68.68,
     59.16,
     59.16,
     67.91,
     67.91,
     67.91,
     67.91
    };

  //histogram for toy 
  const Int_t n_evis_toy = 240;
  const float emin=0;
  const float emax=12;

  Char_t name_fine[1024];
  Char_t name[1024];

  TH1D *h_toy[8];

  for (Int_t i=0;i<8;i++){
    
    sprintf(name_fine,"h_%iAD_fn_fine",i+1);

    h_toy[i] = new TH1D(name_fine,name_fine,n_evis_toy,emin,emax);
    for(int ibin=0;ibin<n_evis_toy;++ibin){
      double x=h_toy[i]->GetBinCenter(ibin+1);
      if(x>0.7) h_toy[i]->SetBinContent(ibin+1,pow(x/E0[i],-(x/E0[i])));
    }

    h_toy[i]->Scale(1.0/h_toy[i]->Integral());

  }

  //rebinned histogram for fit 
  const Int_t n_evis_bins = 37;
  Double_t evis_bins[38]; // Single bins between 0.7 and 1.0 MeV. 0.2 MeV bins from 1.0 to 8.0 MeV. Single bin between 8.0 and 12 MeV. total 37 bins
  evis_bins[0] = 0.7;
  for (Int_t i = 0; i < 36; i++){
    evis_bins[i+1] = 0.2 *i + 1.0;
  }
  evis_bins[37] = 12.0;

  TH1D * h_fit[8];

  for (Int_t i=0;i<8;i++){
    
    sprintf(name,"h_%iAD_fn",i+1);

    h_fit[i] = (TH1D*)h_toy[i]->Rebin(n_evis_bins,name,&evis_bins[0]);
  }

  //saving to file
  TFile *fout = new TFile("P15A_fn_spectrum_IHEP.root","RECREATE");

  for (Int_t i=0;i<8;i++){
    h_toy[i]->Write();
    h_fit[i]->Write();
  }
  

}//end of make_flat_spectrum
