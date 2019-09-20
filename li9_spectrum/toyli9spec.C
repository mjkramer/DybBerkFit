//toy MC to derive true 9Li spectra (beta, neutron, and two alphas) 
//jpochoa@lbl.gov 

//notes to self: 
// Some of the br's are clearly wrong (some not add to 1)
// Best results so far obtained removing newdeltaQ<0 requirement, using fixed newdeltaQ for He5 only and not using Be9 states around 11 MeV (both of them) 

#include <iostream>
#include "TStyle.h"
#include "TRandom.h"
#include "TMath.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"

const int Nbins_fermi = 1000;//<--nbins for fermi beta spectra
const int Ndecays=1e5;//number of fake experiments (i.e. Li9 decays)
const int be9z=4;//atomic number of li9 daughter 

const double li9state=13.6067;//<--MeV

//Be9 LEVELS: from Dan's thesis (and/or Nyman et al. table 2)
const int Nbe9levels = 4;//number of Be9 states Li9 decays into which are not stable
const double be9levels[Nbe9levels]={11.81,11.283,2.78,2.4294};//<--energy levels of Be9 states Li9 decays into and which are not stable
const double be9br[Nbe9levels]={0.03,0.01,0.16,0.30};//<--branching ratios for Be9 states Li9 decays into and which are not stable; 
const double be9widths[Nbe9levels]={0.4,0.575,1.08,0.000077};//<-- widths of be9 states in MeV (from Nyman et al.)

// Be8 LEVELS: from TUNL data (and also Nyman et al. table 1 and p221. of Tilley et al.)
// Note: levels are measured with respect to g.s. of Be9
const int Nbe8levels = 2;//number of Be8 levels that Be9 from Li9 can go into 
const double be8levels[Nbe8levels]={1.665,3.03+1.665};
const double be8widths[Nbe8levels]={5.57e-6,1.51};//note: first level is ground state, second level is 2+

// He5 LEVELS: from TUNL data (and also Nyman et al. table 1 and p221. of Tilley et al.)
// Note: levels are measured with respect to g.s. of Be9
const int Nhe5levels = 2;//number of He5 levels that Be9 from Li9 can go into
const double he5levels[Nhe5levels]={2.467,2.467+1.27};
const double he5widths[Nhe5levels]={0.648,5.57};//note: first level is ground state, second level is 1/2-

// final LEVELS: from TUNL data
// Note: these levels are measured with respect to their respective ground states. Note also that they have to be consistent with other values (i.e. He5Final+he5level(g.s.) better be Be8level(g.s.)+Be8Final). 
const double Be8Final = -0.092;
const double He5Final = -0.894;

// Branching Ratios from Be9 into Be8/He5 (depend on Be9 state, of course). These come mostly from Brown et al. (Table 1). In some (most) cases they are guesstimated (especially for the 11.81 MeV Be9 state). For each Nbe9level, the sum of both be8 and he5 should be 1. 
double br_be8[Nbe9levels][Nbe8levels]={{0.0,0.48},{0.05,0.45},{0.33,0.34},{0.05,0.92}};
double br_he5[Nbe9levels][Nhe5levels]={{0.02,0.50},{0.50,0},{0.33,0.},{0.03,0}};

// Qval for Be8 -> alpha + alpha and for He5 -> n + alpha
//double QvalBe8 = 0.092;//from TUNL, in MeVs
//double QvalHe5 = 0.798;//from TUNL, in MeVs

//parameters that change even less (never?) 
const double m_e = 0.511;//<--mass of the electron in MeV
const double m_be9 = 9.01218 * 931.49;//<--mass of Be9 in MeV
const double m_n = 939.56;//<--mass of the neutron in MeV 
const double m_alpha = 4.0026032*931.49;//<--mass of alpha particle in MeV (from wikipedia)
const double m_be8 = 8.0053051*931.49;//<--mass of the be8 isotope in MeV (from wikipedia)
const double m_he5 = 5.01222*931.49;//<--mass of the he5 isotope in MeV (from wikipedia)

TF1 *func_be8[Nbe8levels];
TF1 *func_he5[Nhe5levels];
TF1 *h_fermispec[Nbe9levels];

TRandom *ran = new TRandom();

//Coulomb corrections
//x is kinetic energy;p[0] is Q; p[1] is m_e
double FermiFunc(double *x, double *p){

  double mom=sqrt(x[0]*x[0]+2*p[1]*x[0]);
  double Etot=x[0]+p[1];
  
  //Coulomb corrections (from BCW's tech note DocDB-8768-v4, p27)
  double firstterm=Etot*1./mom;
  double alpha=-0.0846+0.0248*be9z+0.000237*pow(be9z*1.0,2.0);
  if(x[0]<1.2*p[1]){
    alpha=-0.811+0.0446*be9z+0.000108*pow(be9z*1.0,2.0);
  }
  double beta=0.0115+0.000358*be9z-0.0000617*pow(be9z*1.0,2.0);
  if(x[0]<1.2*p[1]){
    beta = 0.673-0.0182*be9z+0.0000638*pow(be9z*1.0,2.0);
  }  
  double fermifunc = firstterm*exp(alpha+beta*sqrt(x[0]*1./p[1]));
  
  //Fermi's spectra: pE(Q-T)^2 (normalization is not important, as randomly sampling from these hists)
  // p = sqrt(T*T + 2*m_e*T);
  // E = T + m_e;
  double fermispec = mom*Etot*(p[0]-x[0])*(p[0]-x[0]);

  return fermispec*fermifunc;

};

TF1* FermiSpec(double Q){
   
  TF1 *fermi = new TF1("fermi",FermiFunc,0,Q,2);
  fermi->SetParameter(0,Q);
  fermi->SetParameter(1,m_e);
    
  return fermi;

};

double RelBeta(double mass,double Kin){
  
  //beta = p / E;
  return sqrt(Kin*Kin+2*mass*Kin)*1./(Kin+mass);
  
}

double RelGamma(double relbeta){

  //gamma = 1/sqrt(1-beta^2);
  return 1./sqrt(1-pow(relbeta,2));

}

//----------------------------------------------------------------
//Handle the decays Li9->Be8+n (and Be8->alpha+alpha ) 
//Output kinetic energies go into Ekin (0=n,1=alpha1,2=alpha2)
void Be8Decay(int be9state, int be8state, double *Ekin, TF1 *func_in){

  // Some comments: 
  // be8state=0 is easiest case, i.e. decay into be8+n, where be8 is in ground-state; this is always kinematically allowed, and gives practically mono-energetic m except for puny width of Be9 parent state and recoil of Be9
  // be8state=1 is decay into be8+n where be8 is in 2+ excited state @ 2.94MeV 
  // If decay is kinematically forbidden then do other one.
  double totalavailable=be9levels[be9state]-(be8levels[0]+Be8Final);//<--Be8Final is measured w.r.t. Be8 g.s. 
  double deltaQ = -1;
  double newdeltaQ = -1;
  double be9level_fluct = -1;
  double be8level_fluct = -1;
  double Tn_cm = - 1;
  
  //be9level_fluct cannot go above initial Li9 energy and below energy of Be8 
  while(deltaQ<0 || be9level_fluct < be8levels[0]+Be8Final || be9level_fluct > li9state || be8level_fluct > li9state || be8level_fluct < be8levels[0]+Be8Final){// || newdeltaQ < 0){
    be9level_fluct = ran->BreitWigner(be9levels[be9state],be9widths[be9state]);
    if(be9state==3){
      be8level_fluct = func_in->GetRandom();
    } else {
      be8level_fluct = ran->BreitWigner(be8levels[be8state],be8widths[be8state]);
    }
    deltaQ=be9level_fluct-be8level_fluct;
    newdeltaQ=totalavailable-deltaQ;//<--alternative newdeltaQ calculation; superseded below 
  }
  //cout << "expected vs. obtained " << be9levels[be9state]-be8levels[be8state] << "," << deltaQ << endl;//tmp
  
  //relativistic expression (overkill...)
  //Ekin[0] = (pow(deltaQ,2)+2*m_be8*deltaQ)*1./(2*m_n+2*deltaQ+2*m_be8);
  //non-relativistic expression
  Tn_cm = deltaQ*m_be8*1./(m_be8+m_n);
  //Apply boost from recoil of Be9
  double Tbe9_recoil = fabs(0.01*ran->Gaus(0,1));
  double beta_be9 = RelBeta(m_be9,Tbe9_recoil);
  double gamma_be9 = RelGamma(beta_be9);
  double beta_n_cm = RelBeta(m_n,Tn_cm);
  double theta_cm = ran->Uniform(0,2*TMath::Pi());
  Ekin[0]=gamma_be9*(beta_be9*beta_n_cm*TMath::Cos(theta_cm)+1.0)*(Tn_cm+m_n)-m_n;
  double Tbe8 = Tbe9_recoil+deltaQ-Ekin[0];
    
  //-->subsequent decay into alpha+alpha
  //cout << "new and old " << newdeltaQ << "," << -1*Be8Final << endl;//tmp
  newdeltaQ = be8level_fluct-be8levels[0]-Be8Final; 
  //double Talpha1_cm = (pow(newdeltaQ,2)+2*m_alpha*newdeltaQ)*1./(2*m_alpha+2*newdeltaQ+2*m_alpha);//<--relativistic (overkill)
  double Talpha1_cm = newdeltaQ*0.5;//two alphas, each takes half
  double beta_be8 = RelBeta(m_be8,Tbe8);
  double gamma_be8 = RelGamma(beta_be8);
  double beta_alpha1_cm = RelBeta(m_alpha,Talpha1_cm);
  // apply Lorentz Boost with Be8 energy  (m + Tn) = (m+Tn_cm) * gamma * (beta_he5 * beta*n_cm cos(theta_cm)+1)	 
  theta_cm = ran->Uniform(0,2*TMath::Pi());
  Ekin[1] = gamma_be8*(beta_be8*beta_alpha1_cm*TMath::Cos(theta_cm)+1.0)*(Talpha1_cm+m_alpha)-m_alpha;
  Ekin[2] = Tbe8+newdeltaQ-Ekin[1];

  //Get electron energy
  double endpoint=li9state-be9level_fluct;
  if(be9state==3){//2.43MeV state has negligible width, so no need to regenerate beta spectrum everytime
    Ekin[3] = h_fermispec[be9state]->GetRandom();
  } else{ 
    TF1 *betaspec = (TF1*)FermiSpec(li9state-be9level_fluct);
    Ekin[3] = betaspec->GetRandom();
    delete betaspec;
  }
  //--> quicker option (but not strictly energy conserving) 
  //Ekin[3] = h_fermispec[be9state]->GetRandom();
  Ekin[4]=endpoint;//<--residual; used for energy conservation check purposes
  
}//Li9->Be8+n

//--------------------------------------------------------------
//Handle the decays Li9->He5+alpha (and He5->n+alpha)
//Output kinetic energies go into Ekin (0=n,1=alpha1,2=alpha2)
void He5Decay(int be9state, int he5state, double *Ekin, TF1 *func_in){

  //Do He5
  double totalavailable=be9levels[be9state]-(he5levels[0]+He5Final);//<--He5Final measured w.r.t. g.s. of He5
  double deltaQ=-1;
  double be9level_fluct=-1;
  double he5level_fluct=-1;
  double newdeltaQ = -1;
  //be9 state cannot fluctuate above initial 9Li energy or below breakup level  
  while(deltaQ<0 || be9level_fluct > li9state || be9level_fluct < he5levels[0]+He5Final || he5level_fluct > li9state || he5level_fluct< he5levels[0]+He5Final){// || newdeltaQ < 0){
    be9level_fluct = ran->BreitWigner(be9levels[be9state],be9widths[be9state]);
    he5level_fluct = ran->BreitWigner(he5levels[he5state],he5widths[he5state]);
    //he5level_fluct = func_in->GetRandom();
    deltaQ=be9level_fluct-he5level_fluct;
    newdeltaQ=totalavailable-deltaQ;
  }

  //Apply boost from recoil of Be9
  double Talpha1_cm = deltaQ*m_he5*1./(m_he5+m_alpha);
  double Tbe9_recoil = fabs(0.01*ran->Gaus(0,1));
  double beta_be9 = RelBeta(m_be9,Tbe9_recoil);
  double gamma_be9 = RelGamma(beta_be9);
  double beta_alpha1_cm = RelBeta(m_alpha,Talpha1_cm);
  double theta_cm = ran->Uniform(0,2*TMath::Pi());
  Ekin[1]=gamma_be9*(beta_be9*beta_alpha1_cm*TMath::Cos(theta_cm)+1.0)*(Talpha1_cm+m_alpha)-m_alpha;
  double The5 = Tbe9_recoil+deltaQ-Ekin[1];

  //if no boost 
  //Ekin[1] = (pow(deltaQ,2)+2*m_he5*deltaQ)*1./(2*m_alpha+2*deltaQ+2*m_he5);
  //Ekin[1]=deltaQ*m_he5*1./(m_he5+m_alpha);
  //double The5 = deltaQ-Ekin[1];
  
  //-->subsequent decay of He5 into n+alpha

  newdeltaQ = he5level_fluct-he5levels[0]-He5Final;//+ he5level_fluct;
  //cout << "I calculated " << m_he5-m_alpha-m_n << endl;//tmp
  //kinetic energy of neutron in rest frame of He5
  double Tn_cm = newdeltaQ*m_alpha*1./(m_alpha+m_n); 
  double beta_he5 = RelBeta(m_he5,The5);
  double gamma_he5 = RelGamma(beta_he5);
  double beta_n_cm = RelBeta(m_n,Tn_cm);
  // apply Lorentz Boost with He5 energy  
  theta_cm = ran->Uniform(0,2*TMath::Pi());
  Ekin[0] = gamma_he5*(beta_he5*beta_n_cm*cos(theta_cm)+1)*(Tn_cm+m_n)-m_n;
  Ekin[2] = The5+newdeltaQ-Ekin[0];  

  //Get electron energy
  double endpoint=li9state-be9level_fluct;
  TF1 *betaspec = (TF1*)FermiSpec(li9state-be9level_fluct);
  Ekin[3] = betaspec->GetRandom();
  delete betaspec;
  //--> quicker option (but not strictly energy conserving) 
  //Ekin[3] = h_fermispec[be9state]->GetRandom();
  Ekin[4]=endpoint;//<--residual; used for energy conservation check purposes

}//Li9->He5+alpha

void toyli9spec(){

  //setup ******************
  gStyle->SetOptStat(0);
  TFile *savefile = new TFile("toyli9spec.root","RECREATE");
  TTree *tr = new TTree("tr","toy li9 decay MC");
  double be9energy, Te, Tn, Talpha1, Talpha2;
  tr->Branch("be9energy",&be9energy,"be9energy/D");
  tr->Branch("Te",&Te,"Te/D");
  tr->Branch("Tn",&Tn,"Tn/D");
  tr->Branch("Talpha1",&Talpha1,"Talpha1/D");
  tr->Branch("Talpha2",&Talpha2,"Talpha2/D");
  
  //************************

  //Make functions for sampling
  func_be8[0] = new TF1("func_be8_0","TMath::BreitWigner(x,[0],[1])",be8levels[0]-6*be8widths[0],be8levels[0]+6*be8widths[0]);
  func_be8[0]->SetParameter(0,be8levels[0]);
  func_be8[0]->SetParameter(1,be8widths[0]);
  func_be8[1] = new TF1("func_be8_1","TMath::BreitWigner(x,[0],[1])*TMath::Gaus(x,[2],[3])",Be8Final,be8levels[1]+6*be8widths[1]);
  //func_be8[1] = new TF1("func_be8_1","TMath::BreitWigner(x,[0],[1])",be8levels[1]-6*be8widths[1],be8levels[1]+6*be8widths[1]);
  func_be8[1]->SetParameter(0,be8levels[1]);
  func_be8[1]->SetParameter(1,be8widths[1]);
  func_be8[1]->SetParameter(2,2.2);
  func_be8[1]->SetParameter(3,0.15);
  //Note to self: ROOT is weird. When changing range of functions above functions the weird triangular shape of results change. That is why try to make it as narrow as reasonable. 

  /*
  //func_he5[0] = new TF1("func_he5_1","TMath::BreitWigner(x,[0],[1])*TMath::Gaus(x,[2],[3])",He5Final,he5levels[0]+6*he5widths[0]);
  func_he5[0] = new TF1("func_he5_1","TMath::BreitWigner(x,[0],[1])",He5Final,he5levels[0]+6*he5widths[0]);
  func_he5[0]->SetParameter(0,he5levels[0]);
  func_he5[0]->SetParameter(1,he5widths[0]);
  func_he5[0]->SetParameter(2,0);
  func_he5[0]->SetParameter(3,1.0);
  func_he5[1] = new TF1("func_he5_0","TMath::BreitWigner(x,[0],[1])",he5levels[0]-6*he5widths[0],he5levels[0]+6*he5widths[0]);
  func_he5[1]->SetParameter(0,he5levels[0]);
  func_he5[1]->SetParameter(1,he5widths[0]);
  */

  //Make Fermi beta decay spectra for all decay modes and plot them
  TH1F *h_bkgd = new TH1F("h_bkgd","",100,0,12);
  h_bkgd->GetXaxis()->SetTitle("E_{true} (MeV)");
  TCanvas *can_fermispec = new TCanvas("can_fermispec","");
  can_fermispec->Divide(3,2);
  for(int ist=0;ist<Nbe9levels;ist++){
    can_fermispec->cd(ist+1);
    h_fermispec[ist] = (TF1*)FermiSpec(li9state-be9levels[ist]);
    h_bkgd->Draw();
    h_fermispec[ist]->Draw("same");
  }
  can_fermispec->Draw();

  //Make probability histogram to draw decay mode 
  TH1F *h_be9decay = new TH1F("h_be9decay","",Nbe9levels,0,Nbe9levels);
  for(int ist=0;ist<Nbe9levels;ist++){
    h_be9decay->SetBinContent(ist+1,be9br[ist]);
  }

  //Destination histograms
  TH1F *h_be9state = new TH1F("h_be9state","",Nbe9levels,0,Nbe9levels);
  h_be9state->GetXaxis()->SetTitle("E_{true} (MeV)");
  h_be9state->GetXaxis()->SetTitleSize(0.05);
  TH1F *h_be9beta = new TH1F("h_be9beta","",Nbins_fermi,0,li9state*1.1);
  h_be9beta->GetXaxis()->SetTitle("E_{true} (MeV)");
  h_be9beta->GetXaxis()->SetTitleSize(0.05);
  TH1F *h_be9n = new TH1F("h_be9n","",Nbins_fermi,0,li9state*0.65);
  h_be9n->GetXaxis()->SetTitle("E_{true} (MeV)");
  h_be9n->GetXaxis()->SetTitleSize(0.05);
  TH1F *h_be9alpha = new TH1F("h_be9alpha","",Nbins_fermi,0,li9state*0.65);
  h_be9alpha->GetXaxis()->SetTitle("E_{true} (MeV)");
  h_be9alpha->GetXaxis()->SetTitleSize(0.05);
  h_be9alpha->SetLineColor(2);
  
  Char_t histname[1024];
  TH1F *h_Tn[Nbe9levels][Nbe8levels+Nhe5levels];
  TH1F *h_Talpha[Nbe9levels][Nbe8levels+Nhe5levels];
  TH1F *h_Tenergy[Nbe9levels][Nbe8levels+Nhe5levels];
  for(int ibe9=0;ibe9<Nbe9levels;++ibe9){
    for(int ibe8=0;ibe8<Nbe8levels+Nhe5levels;++ibe8){
      sprintf(histname,"h_Tn_be9_%i_dest_%i",ibe9,ibe8);
      h_Tn[ibe9][ibe8] = (TH1F*)h_be9n->Clone(histname);
      h_Tn[ibe9][ibe8]->GetXaxis()->SetTitle("True neutron energy (MeV)");
      h_Tn[ibe9][ibe8]->GetXaxis()->SetLabelSize(0.09);
      h_Tn[ibe9][ibe8]->GetYaxis()->SetLabelSize(0.09);
      h_Tn[ibe9][ibe8]->GetXaxis()->SetTitleSize(0.09);
      h_Tn[ibe9][ibe8]->GetYaxis()->SetTitleSize(0.09);
      h_Tn[ibe9][ibe8]->GetXaxis()->SetTitleOffset(0.95);
      sprintf(histname,"h_Talpha_be9_%i_dest_%i",ibe9,ibe8);
      h_Talpha[ibe9][ibe8] = (TH1F*)h_be9n->Clone(histname);
      h_Talpha[ibe9][ibe8]->SetLineColor(2);
      h_Talpha[ibe9][ibe8]->GetXaxis()->SetTitle("True #alpha energy (MeV)");
      h_Talpha[ibe9][ibe8]->GetXaxis()->SetTitleSize(0.09);
      h_Talpha[ibe9][ibe8]->GetYaxis()->SetTitleSize(0.09);
      h_Talpha[ibe9][ibe8]->GetXaxis()->SetLabelSize(0.09);
      h_Talpha[ibe9][ibe8]->GetYaxis()->SetLabelSize(0.09);
      h_Talpha[ibe9][ibe8]->GetXaxis()->SetTitleOffset(0.95);
      sprintf(histname,"h_Tenergy_be9_%i_dest_%i",ibe9,ibe8);
      h_Tenergy[ibe9][ibe8] = (TH1F*)h_be9beta->Clone(histname);
      h_Tenergy[ibe9][ibe8]->SetLineColor(kGreen+2);
      h_Tenergy[ibe9][ibe8]->GetXaxis()->SetTitle("Q+#alpha+n energy (MeV)");
      h_Tenergy[ibe9][ibe8]->GetXaxis()->SetTitleSize(0.09);
      h_Tenergy[ibe9][ibe8]->GetYaxis()->SetTitleSize(0.09);
      h_Tenergy[ibe9][ibe8]->GetXaxis()->SetLabelSize(0.09);
      h_Tenergy[ibe9][ibe8]->GetYaxis()->SetLabelSize(0.09);
      h_Tenergy[ibe9][ibe8]->GetXaxis()->SetTitleOffset(0.95);
    }
  }
   
  //Main Li9 decay loop
  int progress=0;
  for(int iexp=0;iexp<Ndecays;++iexp){

    //print progress
    progress = iexp*100/Ndecays;
    cout << "\r" << "Progress: " << progress << "\%" << flush;
    
    //Decay from Li9 into Be9 -----------------------------------
    //determine which decay mode (i.e. which be9 state fell into)
    int be9state = (int)h_be9decay->GetRandom();
    h_be9state->Fill(be9state);
              
    //------------------------------------------------------------
    //Handle decays of Be9

    //if(be9state==0 or be9state==1) continue;

    //Depending on be9 state make probability histogram of all possible decays and sample it
    TH1F *h_prob = new TH1F("h_prob","",Nbe8levels+Nhe5levels,0,Nbe8levels+Nhe5levels);
    for(int iprob=0;iprob<Nbe8levels;++iprob){
      h_prob->SetBinContent(iprob+1,br_be8[be9state][iprob]);
    }
    for(int iprob=Nbe8levels;iprob<Nbe8levels+Nhe5levels;++iprob){
      h_prob->SetBinContent(iprob+1,br_he5[be9state][iprob-Nbe8levels]);
    }
    int state = (int)h_prob->GetRandom();
    bool isbe8=true;
    if(state>=Nbe8levels) isbe8=false;
    delete h_prob;
    
    //Do the actual decays 
    double Ekin[5]={0};
    if(isbe8){
      Be8Decay(be9state,state,Ekin,func_be8[state]);
    } else {
      He5Decay(be9state,state-Nbe8levels,Ekin,func_he5[state-Nbe8levels]);    
    }
    
    //Extract variables
    be9energy = be9levels[be9state];
    Te = Ekin[3];
    Tn = Ekin[0];
    Talpha1 = Ekin[1];
    Talpha2 = Ekin[2];
    tr->Fill();

    //Fill histograms  
    h_be9beta->Fill(Ekin[3]);
    h_be9n->Fill(Ekin[0]);
    h_be9alpha->Fill(Ekin[1]);
    h_be9alpha->Fill(Ekin[2]);
    h_Tn[be9state][state]->Fill(Ekin[0]);
    h_Talpha[be9state][state]->Fill(Ekin[1]);
    h_Talpha[be9state][state]->Fill(Ekin[2]);
    h_Tenergy[be9state][state]->Fill(Ekin[0]+Ekin[1]+Ekin[2]+Ekin[4]);
    //cout << "I am summing " << Ekin[0] << "," << Ekin[1] << "," << Ekin[2] << "," << Ekin[4] << endl;//tmp
    //cout << "total: " << Ekin[0]+Ekin[1]+Ekin[2]+Ekin[4] << endl;
    
  }
  
  TCanvas *can_be9st = new TCanvas("can_be9st","Be9 states after Li9 decay");
  can_be9st->Divide(2,2);
  can_be9st->cd(1);
  h_be9state->Scale(1./h_be9state->Integral());
  h_be9state->Draw();
  //h_be9decay->Draw();
  can_be9st->cd(2);
  h_be9beta->Draw();
  can_be9st->cd(3)->SetLogy(1);
  h_be9n->Draw();
  can_be9st->cd(4)->SetLogy(1);
  h_be9alpha->Draw();

  TCanvas *can_Tn = new TCanvas("can_Tn","");
  can_Tn->Divide(Nbe8levels+Nhe5levels,Nbe9levels);
  int totalpads=(Nbe8levels+Nhe5levels)*Nbe9levels;
  for(int ibe9=0;ibe9<Nbe9levels;++ibe9){
    for(int ibe8=0;ibe8<Nbe8levels+Nhe5levels;++ibe8){
      int padnum = totalpads - (ibe9+1)*(Nbe8levels+Nhe5levels) + (ibe8+1);
      can_Tn->cd(padnum);
      h_Tn[ibe9][ibe8]->Draw();
      can_Tn->cd(padnum)->SetLogy(1);
    }
  }
  
  TCanvas *can_Talpha = new TCanvas("can_Talpha","");
  can_Talpha->Divide(Nbe8levels+Nhe5levels,Nbe9levels);
  for(int ibe9=0;ibe9<Nbe9levels;++ibe9){
    for(int ibe8=0;ibe8<Nbe8levels+Nhe5levels;++ibe8){
      int padnum = totalpads - (ibe9+1)*(Nbe8levels+Nhe5levels) + (ibe8+1);
      can_Talpha->cd(padnum);
      h_Talpha[ibe9][ibe8]->Draw();
      can_Talpha->cd(padnum)->SetLogy(1);
    }
  }

  TCanvas *can_Tenergy = new TCanvas("can_Tenergy","");
  can_Tenergy->Divide(Nbe8levels+Nhe5levels,Nbe9levels);
  for(int ibe9=0;ibe9<Nbe9levels;++ibe9){
    for(int ibe8=0;ibe8<Nbe8levels+Nhe5levels;++ibe8){
      int padnum = totalpads - (ibe9+1)*(Nbe8levels+Nhe5levels) + (ibe8+1);
      can_Tenergy->cd(padnum);
      h_Tenergy[ibe9][ibe8]->Draw();
      can_Tenergy->cd(padnum)->SetLogy(1);
    }
  }

  //Write file
  savefile->cd();
  tr->Write();
  

}//end of toyli9spec

// Junk ************************ 

  

