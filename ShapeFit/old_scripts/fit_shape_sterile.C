
#include <iostream>
#include <map>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TMinuit.h"
#include "TMath.h"
#include "TGraph.h"
#include "Predictor.h"
#include "TTree.h"


Predictor *pred;


void minuit_fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *x, Int_t iflag){ // function for minuit minimization
  Double_t sin22t13 = x[0];
  Double_t dm2_31 = x[1];
  Double_t sin22t14 = x[2];
  Double_t dm2_41 = x[3];
  //  f =  pred->CalculateChi2Cov(sin22t13,dm2_31+7.5e-5);
  f =  pred->CalculateChi2Cov(sin22t13,dm2_31+7.5e-5,sin22t14,dm2_41);
}


Double_t CalculateChi2(TH1F *h_s22t13,TH1F* h_data);

void fit_shape_sterile(bool isMC = false,
                       Char_t sig_matrix_filename[256] = "covariance_matrices_unified/matrix_sigsys_full.txt",
                       Char_t bg_matrix_filename[256] = "covariance_matrices_unified/matrix_bgsys_full.txt",
                       Char_t savefilename[256] = "fit_results_data.root",
                       Char_t sig_spectra_filename[256] = "./Spectra/ibd_eprompt_shapes.root",
                       Double_t stat_factor = 1){


  bool isNominalMC = false;

  
  int Nevts=1;
  if(isMC){
    if (!isNominalMC)
      Nevts=3;
  }

  gStyle->SetLabelSize(0.08,"XY");
  gStyle->SetTitleSize(0.04,"XY");

  //  Char_t toymc_filename[256] = "./toyfiles/nominal/toySpectra_allsys_and_stat.root";
  
  
  // Create Predictor
  pred = new Predictor();
  
  //  FluxCalculator* fluxcalc = new FluxCalculator("./Distances/unblinded_baseline.txt","./Flux/SuperHistograms_32week_unblinded_fine_huber-french.root");//<--flux calculator in super-hist mode 
  //  FluxCalculator* fluxcalc = new FluxCalculator("./Distances/unblinded_baseline.txt","./Flux/SuperHistograms_32week_unblinded_fine_ill-vogel.root");//<--flux calculator in super-hist mode 
  FluxCalculator* fluxcalc = new FluxCalculator("./Distances/unblinded_baseline.txt","./Flux/SuperHistograms_32week_unblinded_fine_huber-french.root");//<--flux calculator in super-hist mode 


  pred->EnterFluxCalculator(fluxcalc);

  pred->LoadMainData("./Inputs/Theta13-inputs_32week_inclusive.txt");
  //  pred->LoadMainData("./Inputs/Theta13-inputs_20week_inclusive.txt");
    
  //Int_t ntoys = pred->LoadToyIBDSpec(toymc_filename);
  //cout << "number of toy samples: " << ntoys << endl;

  if(!isMC) {
    pred->LoadIBDSpec(sig_spectra_filename);
    //pred->LoadToyIBDSpec(toymc_filename);
    //pred->LoadToyMCNominalSpec();
  }
  else {
    //    pred->LoadToyIBDSpec(toymc_filename);
    pred->LoadToyIBDSpec(sig_spectra_filename);
    pred->LoadToyMCEntry(0,false);//<--2nd argument is so that it does not attempt to correct before loading Bg Spec. 
  }
  
  pred->LoadBgSpec("./Spectra/accidental_eprompt_shapes.root",
                   "../li9_spectrum/8he9li_nominal_spectrum.root",
                   "../amc_spectrum/amc_spectrum.root",
                   "../fn_spectrum/fn_spectrum.root",
                   "../alpha-n-spectrum/result-DocDB9667.root");//<---load bg afterwards since here is when correct events
  
  
  pred->SetStatFactor(stat_factor);

  const Int_t n_evis_bins = 37;

  Double_t evis_bins[n_evis_bins+1]; // Single bins between 0.7 and 1.0 MeV. 0.2 MeV bins from 1.0 to 8.0 MeV. Single bin between 8.0 and 12 MeV. total 37 bins
  evis_bins[0] = 0.7;
  for (Int_t i = 0; i < n_evis_bins-1; i++){
    evis_bins[i+1] = 0.2 *i + 1.0;
  }
  evis_bins[n_evis_bins] = 12.0;

  
  const Int_t n_enu_bins = 156;
  Double_t enu_bins[n_enu_bins+1]; // testing fine bins
  for (Int_t i = 0; i < n_enu_bins+1; i++){
    enu_bins[i] = 0.05 * i + 1.8;
  }
  // const Int_t n_enu_bins = 39;
  // Double_t enu_bins[n_enu_bins+1]; // 39 bins between 1.8 and 9.6 MeV
  // for (Int_t i = 0; i < n_enu_bins+1; i++){
  //   enu_bins[i] = 0.2 * i + 1.8;
  // }

  pred->SetEvisBins(n_evis_bins,evis_bins);
  pred->SetEnuBins(n_enu_bins,enu_bins);

  
  //pred->LoadEvisToEnuMatrix("matrix_evis_to_enu_rateonly_unified_fine_p12e_unblinded.txt");
  pred->LoadEvisToEnuMatrix("matrix_evis_to_enu_unified_fine_p12e_unblinded.txt");
  //  pred->LoadEvisToEnuMatrix("matrix_evis_to_enu_unified_fine.txt");

  pred->LoadCovMatrix(sig_matrix_filename,bg_matrix_filename);

  //*************************************

  //-->Best fits histograms
  TH2F *h_s22t13_dm2 = new TH2F("h_s22t13_dm2","",100,0,0.2,100,1.5e-3,3.5e-3);
  h_s22t13_dm2->GetXaxis()->SetTitle("sin^{2}(2#theta_{13})");
  h_s22t13_dm2->GetYaxis()->SetTitle("#Delta m^{2}");
  h_s22t13_dm2->SetMarkerStyle(20);
  h_s22t13_dm2->SetMarkerSize(0.8);
  
  //-->Prediction at t13=0
  cout << "Prediction at t13=0" << endl;
  PredSet *predset_0 = pred->MakeOneSuperPrediction(0,-1,true);
  TH1F* h_noosc = predset_0->GetCombinedPred();

  //data (combined all three far ADs)
  TH1F *h_data = (TH1F*)pred->GetCombCorrEvtsSpec(3)->Clone("h_data");
  h_data->Add(pred->GetCombCorrEvtsSpec(4));
  h_data->Add(pred->GetCombCorrEvtsSpec(5));


  //Do the minimization manually
  //----------------------------
  
  //++best parameters
  Double_t minchi2=1e8;
  Double_t bests22t13=0;
  Double_t bestdm213=0;
  Double_t bests22t14=0;
  Double_t bestdm214=0;
  //  PredSet bestpred;
  
  const Int_t nsteps = 1;
  // //  const Int_t nsteps = 201;
  // Double_t s22t13start=0.00;
  // Double_t s22t13end=0.20;
  //  const Int_t nsteps = 101;
  Double_t s22t13start=0.05;
  Double_t s22t13end=0.15;

  //  const Int_t nsteps_dm2 = 101;
  const Int_t nsteps_dm2 = 1;
  Double_t dm213start=1.5e-3;
  Double_t dm213end=3.5e-3;


  const Int_t nsteps_s22t14 = 2;
  //  const Int_t nsteps_s22t14 = 3;
  Double_t s22t14start=0.0;
  Double_t s22t14end=1.0;

  const Int_t nsteps_dm214 = 2;
  //  const Int_t nsteps_dm214 = 3;
  Double_t dm214start=1.0e-4;
  Double_t dm214end=10.0;

  const Int_t nsteps_dm214_fit = 56;

  
  Double_t chi2result[nsteps_dm2][nsteps][nsteps_dm214][nsteps_s22t14];
  Double_t dchi2result[nsteps_dm2][nsteps][nsteps_dm214][nsteps_s22t14];
  Double_t sin22t13[nsteps_dm2][nsteps];
  Double_t dm213[nsteps_dm2][nsteps];

  Double_t chi2_min = 1e6;
  Double_t s2t_min = 1e6;
  Double_t dm2_min = 1e6;
  Double_t s2t14_min = 1e6;
  Double_t dm241_min = 1e6;

  Double_t dm241_ini[nsteps_dm214_fit];
  Double_t dm241_min_local[nsteps_dm214_fit];
  Double_t chi2_min_local[nsteps_dm214_fit];
  Int_t fit_status[nsteps_dm214_fit];

  
  Double_t s22t14_bins[nsteps_s22t14+1];
  Double_t dm214_bins[nsteps_dm214+1];
  
  for(int step_dm214=0;step_dm214<nsteps_dm214+1;++step_dm214){
    dm214_bins[step_dm214]=exp(log(dm214end/dm214start)*((Double_t)step_dm214-0.5)*1./(nsteps_dm214-1)+log(dm214start));
  }
  for(int step_s22t14=0;step_s22t14<nsteps_s22t14+1;++step_s22t14){
    s22t14_bins[step_s22t14]=(s22t14end-s22t14start)*((Double_t)step_s22t14-0.5)*1./(nsteps_s22t14-1)+s22t14start;      
  }

  
  TH1F *h_minchi2 = new TH1F("h_minchi2","",100,0,300);
  h_minchi2->GetXaxis()->SetTitle("Minimum #chi^{2}");
  h_minchi2->GetYaxis()->SetTitle("fake experiments");

  //  TH2F *h_chi2_all = new TH2F("h_chi2_all","",nsteps,s22t13start,s22t13end,nsteps_dm2,dm213start,dm213end);
  //  TH2F *h_chi2_all = new TH2F("h_chi2_all","",nsteps_s22t14,s22t14start,s22t14end,nsteps_dm214,dm214start,dm214end);
  TH2D *h_chi2_all = new TH2D("h_chi2_all","",nsteps_s22t14,s22t14_bins,nsteps_dm214,dm214_bins);
  h_chi2_all->GetXaxis()->SetTitle("sin^2(2#theta_{13})");
  h_chi2_all->GetYaxis()->SetTitle("#Delta m^2");
  TH2D *h_chi2_temp = (TH2D*)h_chi2_all->Clone("h_chi2_temp");

  TDirectory * dir = gDirectory;

  
  //savefile
  TFile *savefile = new TFile(savefilename,"RECREATE");
  TTree * tr = new TTree("tr_fit","fit results");
  tr->Branch("chi2_map",&dchi2result[0][0][0][0],Form("chi2_map[%d][%d][%d][%d]/D",nsteps_dm2,nsteps,nsteps_dm214,nsteps_s22t14)); 
  tr->Branch("chi2_min",&chi2_min,"chi2_min/D"); 
  tr->Branch("dm2_min",&dm2_min,"dm2_min/D"); 
  tr->Branch("s2t_min",&s2t_min,"s2t_min/D"); 
  // new branched for sterile mixing 
  tr->Branch("dm241_min",&dm241_min,"dm241_min/D"); 
  tr->Branch("s2t14_min",&s2t14_min,"s2t14_min/D"); 

  tr->Branch("dm241_ini",&dm241_ini[0],Form("dm241_ini[%d]/D",nsteps_dm214_fit)); 
  tr->Branch("dm241_min_local",&dm241_min_local[0],Form("dm241_min_local[%d]/D",nsteps_dm214_fit)); 
  tr->Branch("chi2_min_local",&chi2_min_local[0],Form("chi2_min_local[%d]/D",nsteps_dm214_fit)); 
  tr->Branch("fit_status",&fit_status[0],Form("fit_status[%d]/I",nsteps_dm214_fit)); 
  
  
  TH2D* h_chi2_map = new TH2D("h_chi2_map","h_chi2_map",
                              nsteps_s22t14,s22t14_bins,nsteps_dm214,dm214_bins);


  dir->cd();

  const Int_t NPars = 4;
  TMinuit * minu = new TMinuit(NPars);
  minu->SetPrintLevel(-1);
  minu->SetFCN(minuit_fcn);


  Int_t ierflag;
  Double_t arglist[2];
  arglist[0] = 1.0;
  
  // Fit by Minuit
  minu->mnexcm("SET ERR",arglist,1,ierflag);
  
  minu->mnparm(0,
               "SinSq2Theta13",
               0.089,
               0.001,
               0,0.2,ierflag);
  minu->mnparm(1,
               "DeltaMSq31",
               0.00241,
               0.0001,
               0.0015,0.0035,ierflag);
  minu->mnparm(2,
               "SinSq2Theta14",
               0.01,
               0.001,
               0,0.2,ierflag);
  minu->mnparm(3,
               "DeltaMSq41",
               0.1,
               0.0001,
               1e-4,10.0,ierflag);
  

  
  //  PredSet *predset=0;
  //  for(int ievt=0;ievt<Nevts;++ievt){
  for(int ievt=0;ievt<80;++ievt){
    
    minchi2=1e8;
    bests22t13=0;
    bestdm213=0;
    bests22t14=0;
    bestdm214=0;
    
    h_chi2_temp->Reset();
    
    if(isMC){
      if (isNominalMC)
        //        pred->LoadToyMCEntry(ievt,true);
        pred->LoadToyMCNominalSpec();
      else
        pred->LoadToyMCEntry(ievt,true);
    }


    Double_t sin22t14;
    Double_t dm214;

    Double_t * grad;
    Double_t fpar,ferr;
    // semi parameter scan:
    for(int step_dm214=0;step_dm214<nsteps_dm214_fit;++step_dm214){
      //      dm214=exp(log(dm214end/dm214start)*step_dm214*1./(nsteps_dm214_fit-1)+log(dm214start));  
      //      dm214=(dm214end-dm214start)*step_dm214*1./(nsteps_dm214-1)+dm214start;  

      //      dm214=(0.501-0.001)*step_dm214*1./(nsteps_dm214_fit-1)+0.001;  // linear scan with 5e-3 steps
      //      dm214=exp(log(1.0/0.01)*step_dm214*1./(nsteps_dm214_fit-1)+log(0.01));  
      // attempt to optimize the initial point distributions
      if (step_dm214 < 3){
        dm214= 0.003 * step_dm214+0.001;
      }else if (step_dm214 < 11){
        dm214= 0.005 * (step_dm214 - 3) +0.01;
      }else{
        dm214= 0.01 * (step_dm214 - 11) +0.05;
      }
      
      cout << " ========== " << step_dm214 << " / " << nsteps_dm214_fit <<" (initial dm2 = " << dm214 << " ) ========== " << endl;
      minu->mnparm(0,
                   "SinSq2Theta13",
                   0.089,
                   0.001,
                   0,0.2,ierflag);
      minu->mnparm(1,
                   "DeltaMSq31",
                   0.00241,
                   0.0001,
                   0.0015,0.0035,ierflag);
      minu->mnparm(2,
                   "SinSq2Theta14",
                   0.02,
                   0.01,
                   0,0.2,ierflag);
      minu->mnparm(3,
                   "DeltaMSq41",
                   dm214,
                   0.5*dm214,
                   1e-4,10.0,ierflag);
      minu->FixParameter(1);
      //      minu->FixParameter(3);
      
      arglist[0] = 10000;
      arglist[1] = 1.0;
      // let only sin22theta13 and sin22theta14 float
      minu->mnexcm("MIGRAD", arglist ,2,ierflag);
      //      Int_t  fit_status = minu->GetStatus();
      //      minu->mnexcm("IMPROVE", arglist ,1,ierflag);
      
      Double_t pars[4];
      minu->GetParameter(0,fpar,ferr);
      pars[0] = fpar;
      minu->GetParameter(1,fpar,ferr);
      pars[1]= fpar;
      minu->GetParameter(2,fpar,ferr);
      pars[2] = fpar;
      minu->GetParameter(3,fpar,ferr);
      pars[3] = fpar;
      Double_t minchi2_tmp = 1e6;
      minu->Eval(4,grad,minchi2_tmp,pars,0);

      cout << " Status = " << ierflag << ", chi2 = " << minchi2_tmp << ", "
           << pars[0] << " " << pars[1] << " " << pars[2] << " " << pars[3] << endl;
      dm241_ini[step_dm214] = dm214;
      dm241_min_local[step_dm214] = pars[3];
      chi2_min_local[step_dm214] = minchi2_tmp;
      fit_status[step_dm214] = ierflag;

      
      if (minchi2_tmp < minchi2){
        minchi2 = minchi2_tmp;
        bests22t13 = pars[0];
        bestdm213 =  pars[1];
        bests22t14 = pars[2];
        bestdm214 =  pars[3];
      }
    }
    
   Double_t best_pars[4] = {bests22t13,bestdm213,bests22t14,bestdm214};
    
    cout << "======== fit results (for checking) :"
         << " " <<  bests22t13 << " " << bestdm213 
         << " " <<  bests22t14 << " " << bestdm214
         << " " << minchi2 << endl;

    
    for(int step_dm2=0;step_dm2<nsteps_dm2;++step_dm2){
      for(int step=0;step<nsteps;++step){
        for(int step_dm214=0;step_dm214<nsteps_dm214;++step_dm214){
          cout << " ========== " << step_dm214 << " / " << nsteps_dm214 <<" ========== " << endl;
          for(int step_s22t14=0;step_s22t14<nsteps_s22t14;++step_s22t14){
            if (nsteps == 1)
              sin22t13[step_dm2][step]=0.089;
            else
              sin22t13[step_dm2][step]=(s22t13end-s22t13start)*step*1./(nsteps-1)+s22t13start;
              
            if (nsteps_dm2 == 1)
              dm213[step_dm2][step]=0.00241;
            else
              dm213[step_dm2][step]=(dm213end-dm213start)*step_dm2*1./(nsteps_dm2-1)+dm213start;  
        
            sin22t14=(s22t14end-s22t14start)*step_s22t14*1./(nsteps_s22t14-1)+s22t14start;      
            //            dm214=(dm214end-dm214start)*step_dm214*1./(nsteps_dm214-1)+dm214start;
            // log scale
            dm214=exp(log(dm214end/dm214start)*step_dm214*1./(nsteps_dm214-1)+log(dm214start));  

            minu->mnparm(0,
                         "SinSq2Theta13",
                         0.089,
                         0.001,
                         0,0.2,ierflag);
            minu->mnparm(1,
                         "DeltaMSq31",
                         0.00241,
                         0.0001,
                         0.0015,0.0035,ierflag);
            minu->mnparm(2,
                         "SinSq2Theta14",
                         sin22t14,
                         0.001,
                         0,0.2,ierflag);
            minu->mnparm(3,
                         "DeltaMSq41",
                         dm214,
                         0.0001,
                         1e-4,10.0,ierflag);
            minu->FixParameter(1);
            minu->FixParameter(2);
            minu->FixParameter(3);

            // let only sin22theta13 float
            arglist[0] = 10000;
            arglist[1] = 1.0;
            minu->mnexcm("MIGRAD", arglist ,2,ierflag);
            

            Double_t pars[4];
            minu->GetParameter(0,fpar,ferr);
            pars[0] = fpar;
            minu->GetParameter(1,fpar,ferr);
            pars[1]= fpar;
            minu->GetParameter(2,fpar,ferr);
            pars[2] = fpar;
            minu->GetParameter(3,fpar,ferr);
            pars[3] = fpar;
            Double_t minchi2_tmp = 1e6;
            minu->Eval(4,grad,minchi2_tmp,pars,0);
            
            chi2result[step_dm2][step][step_dm214][step_s22t14]=minchi2_tmp;
            
            //            chi2result[step_dm2][step][step_dm214][step_s22t14]=pred->CalculateChi2Cov(sin22t13[step_dm2][step],dm213[step_dm2][step]+7.5e-5,sin22t14,dm214);
            //            delete predset;
          }
        }
      }
    }//loop over t13
    minu->Release(2);
    minu->Release(3);
    
    cout << "Event: " << ievt << "; --> best fit: " << bests22t13 << "," << bestdm213 << "; minchi2=" << minchi2 << endl;//tmp
    h_s22t13_dm2->Fill(bests22t13,bestdm213);
    h_minchi2->Fill(minchi2);

    for(int step_dm214=0;step_dm214<nsteps_dm214;++step_dm214){
      for(int step_s22t14=0;step_s22t14<nsteps_s22t14;++step_s22t14){
        h_chi2_map->SetBinContent(step_s22t14+1,step_dm214+1,chi2result[0][0][step_dm214][step_s22t14]-minchi2);
        dchi2result[0][0][step_dm214][step_s22t14] = chi2result[0][0][step_dm214][step_s22t14]-minchi2;
      }
    }
    
    chi2_min = minchi2;
    s2t_min = bests22t13;
    dm2_min = bestdm213;
    s2t14_min = bests22t14;
    dm241_min = bestdm214;

    tr->Fill();

    
  }//loop over toys

  PredSet * bestpred = pred->MakeOneSuperPrediction(bests22t13,bestdm213+7.5e-5,bests22t14,bestdm214,false);


  TCanvas * c0 = new TCanvas ("c0","c0",600,600);
  h_chi2_map->Draw("zcol");


  cout << "======================== Results ============================" << endl;
  cout << "Best sin2(2theta13): " << bests22t13 << endl;
  cout << "Best delta m2(31): " << bestdm213 << endl;
  cout << "Best chi2: " << minchi2 <<  endl;//tmp
  cout << "Best prediction: " << endl;
  bestpred->PrintToScreen();

  Double_t * tmp_vector;

  const Int_t nModes = 4;
  const Int_t nPredictions[nModes] = {9,2,3,1}; // number of predictions to fit for each mode

  Double_t final_covmatrix[nModes][9*n_evis_bins][9*n_evis_bins];
  Double_t final_pred[nModes][9*n_evis_bins];
  Double_t final_pred_null[nModes][9*n_evis_bins];
  Double_t final_obs[nModes][9*n_evis_bins];
  Double_t final_obserror[nModes][9*n_evis_bins];

  
  savefile->cd();
  TH1D::SetDefaultSumw2();

  TString hnames[nModes][9];
  // hnames[0][0] = TString("AD4 compared w/ AD1");
  // hnames[0][1] = TString("AD4 compared w/ AD2");
  // hnames[0][2] = TString("AD4 compared w/ AD3");
  // hnames[0][3] = TString("AD5 compared w/ AD1");
  // hnames[0][4] = TString("AD5 compared w/ AD2");
  // hnames[0][5] = TString("AD5 compared w/ AD3");
  // hnames[0][6] = TString("AD6 compared w/ AD1");
  // hnames[0][7] = TString("AD6 compared w/ AD2");
  // hnames[0][8] = TString("AD6 compared w/ AD3");

  // hnames[1][0] = TString("AD4+AD5+AD6 compared w/ AD1+AD2");
  // hnames[1][1] = TString("AD4+AD5+AD6 compared w/ AD3");

  // hnames[2][0] = TString("AD4+AD5+AD6 compared w/ AD1");
  // hnames[2][1] = TString("AD4+AD5+AD6 compared w/ AD2");
  // hnames[2][2] = TString("AD4+AD5+AD6 compared w/ AD3");

  // hnames[3][0] = TString("AD4+AD5+AD6 compared w/ AD1+AD2+AD3");

  hnames[0][0] = TString("AD4 / AD1");
  hnames[0][1] = TString("AD4 / AD2");
  hnames[0][2] = TString("AD4 / AD3");
  hnames[0][3] = TString("AD5 / AD1");
  hnames[0][4] = TString("AD5 / AD2");
  hnames[0][5] = TString("AD5 / AD3");
  hnames[0][6] = TString("AD6 / AD1");
  hnames[0][7] = TString("AD6 / AD2");
  hnames[0][8] = TString("AD6 / AD3");

  hnames[1][0] = TString("AD4+AD5+AD6 / AD1+AD2");
  hnames[1][1] = TString("AD4+AD5+AD6 / AD3");

  hnames[2][0] = TString("AD4+AD5+AD6 / AD1");
  hnames[2][1] = TString("AD4+AD5+AD6 / AD2");
  hnames[2][2] = TString("AD4+AD5+AD6 / AD3");

  hnames[3][0] = TString("AD4+AD5+AD6 / AD1+AD2+AD3");

  
  TH1D * h_final_obs[nModes][9];
  TH1D * h_final_pred[nModes][9];
  TH1D * h_final_pred_null[nModes][9];
  TH2D * h_final_covmatrix[nModes];

  Int_t n_rebinned_evis_bins = pred->GetNumRebinnedEvisBins();
  Double_t * rebinned_evis_bins = pred->GetRebinnedEvisBins();
  for (Int_t iMode = 0; iMode < nModes; iMode++){
    for (Int_t i = 0; i < nPredictions[iMode]; i++){
      h_final_obs[iMode][i] = new TH1D(Form("h_final_obs_mode%d_%d",iMode,i),hnames[iMode][i].Data(),n_rebinned_evis_bins,rebinned_evis_bins);
      h_final_obs[iMode][i]->GetXaxis()->SetTitle("Visible energy");
      h_final_obs[iMode][i]->GetYaxis()->SetTitle("Events per day");
      h_final_pred[iMode][i] = new TH1D(Form("h_final_pred_mode%d_%d",iMode,i),hnames[iMode][i].Data(),n_rebinned_evis_bins,rebinned_evis_bins);
      h_final_pred_null[iMode][i] = new TH1D(Form("h_final_pred_null_mode%d_%d",iMode,i),hnames[iMode][i].Data(),n_rebinned_evis_bins,rebinned_evis_bins);
    }
    h_final_covmatrix[iMode]= new TH2D(Form("h_final_covmatrix_mode%d",iMode),"Final covariance matrix",
                                       nPredictions[iMode]*n_rebinned_evis_bins,0,nPredictions[iMode]*n_rebinned_evis_bins,
                                       nPredictions[iMode]*n_rebinned_evis_bins,0,nPredictions[iMode]*n_rebinned_evis_bins);
  }

  TH1D * h_final_obs_ratio[nModes][9];
  TH1D * h_final_pred_ratio[nModes][9];
  TH1D * h_final_pred_null_ratio[nModes][9];
  
  dir->cd();

  PredSet * nullpred = pred->MakeOneSuperPrediction(0,-1,false);

  // try to dump final covariance matrix
  
  ofstream fout_cmat("final_covmatrix_mode1.txt");
  tmp_vector = pred->GetFinalCovMatrix(-1,bestpred,1);
  for (Int_t i = 0; i < nPredictions[1]*n_rebinned_evis_bins; i++){
    for (Int_t j = 0; j < nPredictions[1]*n_rebinned_evis_bins; j++){
      fout_cmat << " " << tmp_vector[i*nPredictions[1]*n_rebinned_evis_bins+j];
    }
   fout_cmat << endl;
  }
  fout_cmat.close();


  /////
  
  for (Int_t iMode = 0; iMode < nModes; iMode++){
    tmp_vector = pred->GetFinalCovMatrix(0,bestpred,iMode);

    for (Int_t i = 0; i < nPredictions[iMode]*n_rebinned_evis_bins; i++){
      for (Int_t j = 0; j < nPredictions[iMode]*n_rebinned_evis_bins; j++){
        final_covmatrix[iMode][i][j] = tmp_vector[i*nPredictions[iMode]*n_rebinned_evis_bins+j];
          //      cout << " " << final_covmatrix[i][j];
          }
      //    cout << endl;
    }
    tmp_vector = pred->GetFinalPred(bestpred,iMode);
    for (Int_t i = 0; i < nPredictions[iMode]*n_rebinned_evis_bins; i++){
      final_pred[iMode][i] = tmp_vector[i];
    }
    
    
    tmp_vector = pred->GetFinalObs(iMode);
    for (Int_t i = 0; i < nPredictions[iMode]*n_rebinned_evis_bins; i++){
      final_obs[iMode][i] = tmp_vector[i];
    }
    
    tmp_vector = pred->GetFinalObsError(iMode);
    for (Int_t i = 0; i < nPredictions[iMode]*n_rebinned_evis_bins; i++){
      final_obserror[iMode][i] = tmp_vector[i];
    }
    
    
    tmp_vector = pred->GetFinalPred(nullpred,iMode);
    for (Int_t i = 0; i < nPredictions[iMode]*n_rebinned_evis_bins; i++){
      final_pred_null[iMode][i] = tmp_vector[i];
    }

  
    for (Int_t iPred = 0; iPred <nPredictions[iMode]; iPred++){
      for (Int_t i = 0; i <n_rebinned_evis_bins; i++){
        h_final_obs[iMode][iPred]->SetBinContent(i+1,final_obs[iMode][iPred*n_rebinned_evis_bins+i]);
        h_final_obs[iMode][iPred]->SetBinError(i+1,final_obserror[iMode][iPred*n_rebinned_evis_bins+i]);
        
        
        h_final_pred[iMode][iPred]->SetBinContent(i+1,final_pred[iMode][iPred*n_rebinned_evis_bins+i]);
        h_final_pred[iMode][iPred]->SetBinError(i+1,sqrt(final_covmatrix[iMode][iPred*n_rebinned_evis_bins+i][iPred*n_rebinned_evis_bins+i]));
        
        h_final_pred_null[iMode][iPred]->SetBinContent(i+1,final_pred_null[iMode][iPred*n_rebinned_evis_bins+i]);
        h_final_pred_null[iMode][iPred]->SetBinError(i+1,0);
      }

      h_final_pred[iMode][iPred]->SetLineStyle(2);
      h_final_pred[iMode][iPred]->SetLineColor(2);
      h_final_pred[iMode][iPred]->SetFillStyle(1001);
      //  h_final_pred[iMode][iPred]->SetFillColor(kGray);
      h_final_pred[iMode][iPred]->SetFillColor(2);
      h_final_pred_null[iMode][iPred]->SetLineColor(4);

      savefile->cd();
      
      h_final_obs_ratio[iMode][iPred] = (TH1D*)h_final_obs[iMode][iPred]->Clone(Form("h_final_obs_ratio_mode%d_%d",iMode,iPred));
      h_final_obs_ratio[iMode][iPred]->GetYaxis()->SetTitle("Ratio to null osc. prediction");
      h_final_pred_ratio[iMode][iPred] = (TH1D*)h_final_pred[iMode][iPred]->Clone(Form("h_final_pred_ratio_mode%d_%d",iMode,iPred));
      h_final_pred_null_ratio[iMode][iPred] = (TH1D*)h_final_pred_null[iMode][iPred]->Clone(Form("h_final_pred_null_ratio_mode%d_%d",iMode,iPred));
      h_final_obs_ratio[iMode][iPred]->Divide(h_final_pred_null[iMode][iPred]);
      h_final_pred_ratio[iMode][iPred]->Divide(h_final_pred_null[iMode][iPred]);
      h_final_pred_null_ratio[iMode][iPred]->Divide(h_final_pred_null[iMode][iPred]);
      
      dir->cd();

      
      
    }


    for (Int_t i = 0; i <nPredictions[iMode]*n_rebinned_evis_bins; i++){
      for (Int_t j = 0; j <nPredictions[iMode]*n_rebinned_evis_bins; j++){
        h_final_covmatrix[iMode]->SetBinContent(i+1,j+1,final_covmatrix[iMode][i][j]/sqrt(final_covmatrix[iMode][i][i]*final_covmatrix[iMode][j][j]));
      }
    }
    
    

  }  

  TCanvas * c2 = new TCanvas("c2","c2",1200,1000);
  TCanvas * c3 = new TCanvas("c3","c3",600,600);

  c2->Print("fit_results.pdf[");
  for (Int_t iMode = 0; iMode < nModes; iMode++){
    c2->Clear();
    c3->Clear();
    if (iMode == 0){
      c2->Divide(3,3);
    }
    if (iMode == 1){
      c2->SetWindowSize(1200,600);
      c2->Divide(2,1);
    }
    if (iMode == 2){
      c2->SetWindowSize(1200,600);
      c2->Divide(3,1);
    }
    if (iMode == 3){
      c2->SetWindowSize(1200,1000);
      c2->Divide(1,1);
    }

    for (Int_t iPred = 0; iPred <nPredictions[iMode]; iPred++){
      c2->cd(iPred+1);
      TVirtualPad * pad = c2->GetPad(iPred+1);
      pad->Divide(1,2,0,0);
      pad->cd(1);
      pad->GetPad(1)->SetRightMargin(0.05);;
      pad->GetPad(1)->SetLeftMargin(0.15);;
      h_final_obs[iMode][iPred]->SetTitle(";;Events per day");
      h_final_obs[iMode][iPred]->SetTitleOffset(1.5,"Y");
      h_final_obs[iMode][iPred]->GetXaxis()->SetRangeUser(0.7,7.9);
      h_final_obs[iMode][iPred]->Draw();
      h_final_pred_null[iMode][iPred]->Draw("hist same");
      h_final_pred[iMode][iPred]->Draw("e2 same");
      h_final_obs[iMode][iPred]->Draw("same");
  
      pad->cd(2);
      pad->GetPad(2)->SetRightMargin(0.05);;
      pad->GetPad(2)->SetLeftMargin(0.15);;
      h_final_pred_ratio[iMode][iPred]->SetTitle(";;Ratio");
      h_final_pred_ratio[iMode][iPred]->GetXaxis()->SetRangeUser(0.7,7.9);
      h_final_pred_ratio[iMode][iPred]->SetTitleOffset(1.5,"Y");
      h_final_pred_ratio[iMode][iPred]->SetMinimum(0.8);
      h_final_pred_ratio[iMode][iPred]->SetMaximum(1.2);
      h_final_pred_ratio[iMode][iPred]->Draw("e2");
      h_final_pred_null_ratio[iMode][iPred]->Draw("hist same");
      h_final_obs_ratio[iMode][iPred]->Draw("e same");
        
    }

    c2->Print("fit_results.pdf");
    c3->cd();
      
    h_final_covmatrix[iMode]->Draw("zcol");
    c3->Print("fit_results.pdf");
  
  }
  c2->Print("fit_results.pdf]");
  
  /////////////////////////////////////////////////////////////////////////////////////////
  

  TH2F *h_chi2_all_cl = (TH2F*)h_chi2_all->Clone("h_chi2_all_cl");
  if(Nevts>1){
    TCanvas *can_bestfits = new TCanvas("can_bestfits","Best fits");
    can_bestfits->Divide(2,2);
    can_bestfits->cd(1);
    h_s22t13_dm2->SetStats(0);
    h_s22t13_dm2->Draw();
    can_bestfits->cd(2);
    h_minchi2->Draw();
    can_bestfits->cd(3);
    h_chi2_all->SetStats(0);
    h_chi2_all->Scale(1./Nevts);
    h_chi2_all->Draw("colz");
    can_bestfits->cd(4);
    double contours[2];contours[0]=2.3;contours[1]=4.61;//for 3sigma is 9.21
    //double colors[2]={4,2};
    h_chi2_all_cl->SetContour(2,contours);
    h_chi2_all_cl->Draw("cont1");
    h_s22t13_dm2->Draw("same");
  }
  
  savefile->cd();
  // h_chi2_map->Write();
  // chi2graph->Write();
  // dchi2graph->Write();
  h_s22t13_dm2->Write();
  h_minchi2->Write();
  savefile->Write();
  h_chi2_all->Write();
  h_chi2_all_cl->Write();
  savefile->Close();


  // dump final stat error and weighted mean coefficients for later use

  ofstream fout_wmean("wmean_coeff.txt");
  for (Int_t i = 0; i < n_rebinned_evis_bins; i++){
    for (Int_t j = 0; j < 3; j++){
      fout_wmean << "\t"<<pred->GetWmeanCoeff(i,j);
    }
    fout_wmean << endl;
  }
  fout_wmean.close();
  


  ofstream fout_stat("covmatrix_nearsite_stat.txt");
  tmp_vector = pred->GetFinalCovMatrix(0,bestpred,0);
  for (Int_t i = 0; i < nPredictions[0]*n_rebinned_evis_bins; i++){
    for (Int_t j = 0; j < nPredictions[0]*n_rebinned_evis_bins; j++){
      fout_stat << pred->GetFracStatError(i,j) << endl;
      
    }
   fout_stat << endl;
  }
  fout_stat.close();
  
  
  ofstream fout_stat_total("covmatrix_total_stat.txt");
  tmp_vector = pred->GetFinalCovMatrix(-1,bestpred,0);
  for (Int_t i = 0; i < nPredictions[0]*n_rebinned_evis_bins; i++){
    for (Int_t j = 0; j < nPredictions[0]*n_rebinned_evis_bins; j++){
      fout_stat_total << pred->GetFracStatError(i,j) << endl;
      
    }
   fout_stat_total << endl;
  }
  fout_stat_total.close();
  
  ofstream fout_bg_total("covmatrix_frac_bg.txt");
  for (Int_t i = 0; i < nPredictions[0]*n_rebinned_evis_bins; i++){
    for (Int_t j = 0; j < nPredictions[0]*n_rebinned_evis_bins; j++){
      fout_bg_total << pred->GetFracBgError(i,j) << endl;
      
    }
   fout_bg_total << endl;
  }
  fout_bg_total.close();
  
  
}//end of Test

Double_t CalculateChi2(TH1F *h_s22t13,TH1F* h_data){

  Double_t chi2out=0; 
  int nbins=0;
  for(int ibin=1;ibin<h_s22t13->GetXaxis()->FindBin(7);++ibin){
    nbins++;
 
    double bdata=h_data->GetBinContent(ibin); 
    double bs22t13=h_s22t13->GetBinContent(ibin);
    double edata=h_data->GetBinError(ibin);
    if(edata!=0 && bdata!=0){
      chi2out+=pow(bdata-bs22t13,2)*1./pow(edata,2);
    }
  }

  return chi2out;//*1./(nbins-1);
}
