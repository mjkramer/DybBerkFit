#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TF1.h"
#include "TSystem.h"
#include "readtoyli9tree_2018model.C" //need Nneutron from there, as well as reso_func_bcw

//When savestr=="" it's intended to work as not showing plots; if savestr!="" then will show plots and save things to a file
TH1F* readtoyhe8tree(TFile *nlfile,string savestr="", float scale_neutron_quenching=1.){

  gStyle->SetOptStat(0);
  TRandom3 *gRandom3 = new TRandom3();

  // INPUT *******************************
  TFile *infile = new TFile("toyhe8spec_v1.root","READ");
  TTree *tr = (TTree*)infile->Get("tr");
  //TFile *nlfile = new TFile("nl_models/Model1.root","READ");
  nlfile->cd();

  //scale_neutron_quenching=1.;

  ifstream infile_alpha("./BCW-Model/alpha.dat");
  ifstream infile_neutron("./BCW-Model/neutron.dat");

  TF1 * reso_func = new TF1("reso_func",reso_func_bcw,0,20,3);
  reso_func->SetParameters(0.075,0.009); // based on Bryce's TN

  // *************************************

  TFile *savefile;
  if(savestr!="") savefile = new TFile(savestr.c_str(),"RECREATE");

  //electronics
  TGraph *elec_nl = (TGraph*)nlfile->Get("electronicsNL");

  //read and plot BCW model (quenching) ********************************
  double energy_in,model_in;
  double Alpha_energy[Nalpha],Alpha_model[Nalpha];
  double Neutron_energy[Nneutron],Neutron_model[Nneutron];
  int cont=0;
  while(infile_alpha >> energy_in >> model_in){
    //cout << "reading " << energy_in << model_in << endl;
    Alpha_energy[cont]=energy_in;
    Alpha_model[cont]=model_in;
    ++cont;
  }
  cont=0;
  while(infile_neutron >> energy_in >> model_in){
    //cout << "reading " << energy_in << model_in << endl;
    Neutron_energy[cont]=energy_in;
    Neutron_model[cont]=model_in;
    ++cont;
  }

  // infile_alpha.close();
  //infile_neutron.close();

  //right now getting alphas and neutrons from bcw, electron and gammas from ihep
  TGraph *graph_bcw_alpha = new TGraph(Nalpha,Alpha_energy,Alpha_model);
  graph_bcw_alpha->SetLineColor(11);
  graph_bcw_alpha->SetMarkerColor(11);
  graph_bcw_alpha->GetXaxis()->SetTitle("True energy (MeV)");
  graph_bcw_alpha->GetYaxis()->SetTitle("E_{vis}/E_{true}");

  TGraph *graph_ihep_electron = (TGraph*)nlfile->Get("electronScintNL");

  TGraph *graph_bcw_neutron = new TGraph(Nneutron,Neutron_energy,Neutron_model);
  graph_bcw_neutron->SetLineColor(4);
  graph_bcw_neutron->SetMarkerColor(4);
  graph_bcw_neutron->GetXaxis()->SetTitle("True energy (MeV)");
  graph_bcw_neutron->GetYaxis()->SetTitle("E_{vis}/E_{true}");

  TGraph *graph_ihep_gamma = (TGraph*)nlfile->Get("gammaScintNL");

  TCanvas *can_nl = new TCanvas("can_nl","Quenching non-lin");
  can_nl->cd();
  graph_ihep_electron->SetTitle("NL Model (IHEP & BCW)");
  graph_ihep_electron->Draw("APL");
  graph_bcw_alpha->Draw("PL same");
  graph_bcw_neutron->Draw("PL same");
  graph_ihep_gamma->Draw("PL same");
  graph_ihep_electron->GetYaxis()->SetRangeUser(0,1.2);
  // *****************************************

  TCanvas *can_elec = new TCanvas("can_elec","BCW electronics non-lin");
  can_elec->cd();
  elec_nl->Draw();

  //Set branches
  double Te,Tn,Tgamma;
  tr->SetBranchAddress("Te",&Te);
  tr->SetBranchAddress("Tn",&Tn);
  tr->SetBranchAddress("Tgamma",&Tgamma);
  int nentries = tr->GetEntries();

  //Define histograms
  TH1F* h_eVis = new TH1F("h_eVis","",240,0,12);
  h_eVis->GetXaxis()->SetTitle("E_{vis} (MeV)");
  h_eVis->GetYaxis()->SetTitle("# of entries");
  h_eVis->SetLineColor(4);
  TH1F *h_eVisAllSmeared = (TH1F*)h_eVis->Clone("h_eVisAllSmeared");
  h_eVisAllSmeared->SetLineColor(2);
  TH1F *h_eVisElectronSmeared = (TH1F*)h_eVis->Clone("h_eVisElectronSmeared");
  h_eVisElectronSmeared->SetLineColor(kGreen+1);
  TH1F *h_eVisNeutronSmeared = (TH1F*)h_eVis->Clone("h_eVisNeutronSmeared");
  h_eVisNeutronSmeared->SetLineColor(4);
  TH1F *h_eVisGammaSmeared = (TH1F*)h_eVis->Clone("h_eVisGammaSmeared");
  h_eVisGammaSmeared->SetLineColor(4);

  TH2F *h_eVis_eTrue = new TH2F("h_eVis_eTrue","",500,0,12,500,0,1);
  h_eVis_eTrue->GetXaxis()->SetTitle("Electron true energy (MeV)");
  h_eVis_eTrue->GetYaxis()->SetTitle("Total visible energy (MeV)");

  //loop over tree
  int srCtr=0;
  while(tr->GetEntry(srCtr)){
    ++srCtr;

    //display progress
    int progress=0;
    progress = srCtr*100/nentries;
    cout << "\r" << progress << "%" << flush;

    //apply energy model to particles
    double eVisElectron = Te*graph_ihep_electron->Eval(Te);
    double eVisNeutron = Tn*graph_bcw_neutron->Eval(Tn);
    double eVisGamma = Tgamma*graph_ihep_gamma->Eval(Tgamma);

    //sum all energies
    double eVisAll = eVisElectron+eVisNeutron+eVisGamma;

    //apply electronics nl
    eVisElectron *= elec_nl->Eval(eVisElectron);
    eVisNeutron *= elec_nl->Eval(eVisNeutron);
    eVisGamma *= elec_nl->Eval(eVisGamma);
    eVisAll = eVisElectron+eVisNeutron+eVisGamma;

    //apply resolution
    double sigma = reso_func->Eval(eVisElectron);
    double eVisElectronSmeared = gRandom3->Gaus(eVisElectron,sigma);
    sigma = reso_func->Eval(eVisNeutron);
    double eVisNeutronSmeared = gRandom3->Gaus(eVisNeutron,sigma);
    sigma = reso_func->Eval(eVisGamma);
    double eVisGammaSmeared = gRandom3->Gaus(eVisGamma,sigma);
    sigma =reso_func->Eval(eVisAll);
    double eVisAllSmeared = gRandom3->Gaus(eVisAll,sigma);

    //fill histograms
    h_eVisAllSmeared->Fill(eVisAllSmeared);
    h_eVisNeutronSmeared->Fill(eVisNeutronSmeared);
    h_eVisElectronSmeared->Fill(eVisElectronSmeared);
    h_eVisGammaSmeared->Fill(eVisGammaSmeared);

    //h_eVis_eTrue->Fill(Te,eVisAll);
    h_eVis_eTrue->Fill(eVisAll,eVisElectron*1/eVisAll);

  }//end of main loop

  //TCanvas *can = new TCanvas("can","");
  //h_eVisAll->Draw();
  //h_eVisNeutronFluct1->Draw("same");
  //h_eVisNeutronFluct2->Draw("same");

  TCanvas *c2 = new TCanvas("c2","");
  c2->cd();
  h_eVisElectronSmeared->Draw();
  h_eVisAllSmeared->Draw("same");

  if(savestr!=""){
    savefile->cd();
    h_eVisElectronSmeared->Write();
    h_eVisNeutronSmeared->Write();
    h_eVisGammaSmeared->Write();
    h_eVisAllSmeared->Write();
  }

  if(savestr==""){
    delete h_eVisElectronSmeared;
    delete h_eVisNeutronSmeared;
    delete h_eVisGammaSmeared;
    delete gRandom3;
    delete reso_func;
    delete can_nl;
    delete can_elec;
    delete c2;
    delete graph_bcw_neutron;
    delete h_eVis;
    delete h_eVis_eTrue;
    infile->Close();
  }

  if (savestr != "") savefile->Write();
  return h_eVisAllSmeared;

}//end of readtoyhe8tree

/*
//From Yasu
Double_t nl_func_ihep_electron(Double_t * x, Double_t * par){
Double_t e_true = x[0];
Double_t escale_par = par[5]; // Add flat energy scale parameter
Double_t escale_offset = par[6]; // Add fixed energy offset

Double_t nl_fac_511 = par[4]; // magic factor read from Soren's slides


// if e_true  < 0, assume e_true = 0 and calculate non-linearlity factor
// Then apply to the true positron + gamma energy
// this is needed since true energy sometimes become below 0 after IAV correction.

if (e_true < 0){
e_true = 0;
}

Double_t e_nonlinear_fac
= ((par[0] + par[3]*e_true)/(1+par[1]*exp(-par[2]*e_true))*e_true)
/ (e_true);

//  Double_t e_nonlinear = x[0] * e_nonlinear_fac * escale_par + escale_offset;

return e_nonlinear_fac;

}*/
