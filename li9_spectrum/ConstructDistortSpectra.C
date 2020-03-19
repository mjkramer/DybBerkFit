#include <iostream>
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TRandom3.h"
//#include "readtoyli9tree.C"
#include "readtoyhe8tree.C"

// varies neutron, alpha quenching by established amounts. Then sums Li9 with He8, also varying relative proportion between these two. Produces two files: nominal 8He/9Li, and another one with all the distortion functions

const int Nexp=250;

void ConstructDistortSpectra(){

  TRandom3 *gRandom3 = new TRandom3();
  gRandom3->SetSeed(0.1);

  //INPUT ******************
  double cutoff=0.7;//<--energy cutoff
  double he8li9_frac = 0.055;
  double he8li9_frac_var = 0.;
  double neutron_quenching=1;
  double neutron_quenching_var=1;
  double alpha_quenching=1;
  double alpha_quenching_var=1;

  TFile *nlfile = new TFile("Model1.root","READ");//<--to feed to readtoyli9/he8spectra.C macros
  string nominalname="8he9li_nominal_spectrum_frac_0.055.root";
  string distortname="8he9li_distort_neutron100_alpha100_frac0.055_N250.root";

  // ***********************

  //construct nominal spectrum ***************************
  TH1F *h9li_nominal = readtoyli9tree(nlfile,"",neutron_quenching,alpha_quenching);
  TH1F *h8he_nominal = readtoyhe8tree(nlfile,"",neutron_quenching);
  //-->scale to number of entries
  h9li_nominal->Scale(1./h9li_nominal->Integral());
  h8he_nominal->Scale(1./h8he_nominal->Integral());

  TH1F *h_nominal = (TH1F*)h9li_nominal->Clone("h_nominal");
  h_nominal->Scale(1-he8li9_frac);
  h_nominal->Add(h8he_nominal,he8li9_frac);

  //assert(h_nominal->Integral()>0.99 && h_nominal->Integral()<1.01);

  //apply cutoff
  for(int ibin=1;ibin<=h_nominal->GetXaxis()->GetNbins();++ibin){
    if(h_nominal->GetBinCenter(ibin)<cutoff) h_nominal->SetBinContent(ibin,0);
  }

  TFile *nominalfile = new TFile(nominalname.c_str(),"RECREATE");
  nominalfile->cd();
  h_nominal->Write();

  //make distortion functions ****************************
  TH1F *h_distort = (TH1F*)h_nominal->Clone("h_distort");
  TFile *distortfile = new TFile(distortname.c_str(),"RECREATE");
  TTree *tr_distort = new TTree("tr_distort","Tree containing distortion functions for 8He/9Li");
  tr_distort->Branch("h_distort","TH1F",h_distort);

  for(int iexp=0;iexp<Nexp;++iexp){

    cout << "In fake experiment " << iexp+1 << "/" << Nexp << endl;

    //do necessary fluctuations
    double exp_alpha_quenching=-1;
    double exp_neutron_quenching=-1;
    double exp_frac=-1;
    while(exp_neutron_quenching<0){
      exp_neutron_quenching = gRandom3->Gaus(neutron_quenching,neutron_quenching_var);
    }
    while(exp_alpha_quenching<0){
      exp_alpha_quenching = gRandom3->Gaus(alpha_quenching,alpha_quenching_var);
    }
    while(exp_frac<0){
      exp_frac = gRandom3->Gaus(he8li9_frac,he8li9_frac_var);
    }

    //regenerate spectra
    cout << "Parameters: " << exp_neutron_quenching << "," << exp_alpha_quenching << "," << exp_frac << endl;
    TH1F *h9li_temp = readtoyli9tree(nlfile,"",exp_neutron_quenching,exp_alpha_quenching);
    cout << "Generated li9" << endl;
    TH1F *h8he_temp = readtoyhe8tree(nlfile,"",exp_neutron_quenching);
    cout << "Generated he8" << endl;
    h9li_temp->Scale(1./h9li_temp->Integral());
    h8he_temp->Scale(1./h8he_temp->Integral());

    h_distort->Reset();
    h_distort->Add(h9li_temp,1-exp_frac);
    h_distort->Add(h8he_temp,exp_frac);
    h_distort->Divide(h_nominal);

    distortfile->cd();
    tr_distort->Fill();

    delete h9li_temp;
    delete h8he_temp;


  }//end of fake experiment generation

  distortfile->cd();
  tr_distort->Write();

  TCanvas *can = new TCanvas("can","can");
  can->cd();
  h9li_nominal->Draw();
  h8he_nominal->Draw("same");
  h_nominal->Draw("same");

}
