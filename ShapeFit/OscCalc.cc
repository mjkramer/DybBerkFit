#include "OscCalc.h"

using namespace std;

OscCalc::OscCalc(){

  //default parameters
  hierarchy=1;//-1 for inverted
  deltam2_32_default=2.44e-3;//eV2
  deltam2_21_default=7.53e-5;//eV2
  deltam2_41_default=1.0;
  deltam2_factor_default=5.17e-5;//eV2, conversion factor between deltam2_ee and deltam2_32.
  sin22t12_default=0.851;//from pdg-2018
  theta14_default=0;

  //set operating parameters to default parameters
  deltam2_21=deltam2_21_default;
  deltam2_32=deltam2_32_default;
  deltam2_31=deltam2_32+hierarchy*deltam2_21;
  deltam2_41=deltam2_41_default;
  deltam2_factor=deltam2_factor_default;
  //assuming 4th state (sterile) is heavier than the others
  deltam2_42=deltam2_41-deltam2_21;//<--2 is heavier than 1
  deltam2_43=deltam2_41-hierarchy*deltam2_31;//<--depends on hierarchy
  sin22t12=sin22t12_default;
  theta14=theta14_default;

  //Load default unoscillated spectrum
  //TFile *anspecfile = new TFile("./UnoscSpectrum/UnoscSpectrum.root","READ");
  //anspecdef = (TH1F*)anspecfile->Get("anspec");

}

OscCalc::~OscCalc(){}

void OscCalc::SetDeltaM2_21(double deltam2){

  deltam2_21 = deltam2;

}

void OscCalc::SetDeltaM2_32(double deltam2){

  deltam2_32=deltam2;
  deltam2_31=deltam2_32+hierarchy*deltam2_21;
  deltam2_ee=deltam2_32+hierarchy*deltam2_factor;

}

void OscCalc::SetDeltaM2_31(double deltam2){

  deltam2_31=deltam2;
  deltam2_32=deltam2_31-hierarchy*deltam2_21;
  deltam2_ee=deltam2_32+hierarchy*deltam2_factor;

}

void OscCalc::SetDeltaM2_ee(double deltam2){

  deltam2_ee=deltam2;
  deltam2_32=deltam2_ee-hierarchy*deltam2_factor;
  deltam2_31=deltam2_32+hierarchy*deltam2_21;

}

void OscCalc::SetDeltaM2_41(double deltam2){

  deltam2_41=deltam2;
  deltam2_42=deltam2_41-deltam2_21;
  deltam2_43=deltam2_41-hierarchy*deltam2_31;//<--depends on hierarchy

}

void OscCalc::SetTheta14(double theta_in){

  theta14=theta_in;

}

void OscCalc::SetTheta12(double theta_in){

  sin22t12=pow(TMath::Sin(2*theta_in),2);

}

void OscCalc::SetHierarchy(int hie){

  hierarchy=hie;
  //leave deltam2_ee the same and recalculate deltam2_32 and deltam2_31 from it
  deltam2_32=deltam2_ee-hierarchy*deltam2_factor;
  deltam2_31=deltam2_32+hierarchy*deltam2_21;

}

double OscCalc::GetDeltaM2_21(){

  return deltam2_21;

}

double OscCalc::GetDeltaM2_32(){

  return deltam2_32;

}

double OscCalc::GetDeltaM2_31(){

  return deltam2_31;

}

double OscCalc::GetDeltaM2_ee(){

  return deltam2_ee;

}

double OscCalc::GetDeltaM2_41(){

  return deltam2_41;

}

double OscCalc::GetTheta14(){

  return theta14;

}

double OscCalc::GetTheta12(){

  return TMath::ASin(sqrt(sin22t12))*0.5;

}

double OscCalc::DeltaTerm(double deltam2_in, double L, double E){

  return pow(sin(1.267*deltam2_in*L/E),2);

}
double OscCalc::DeltaTermInt(double deltam2_in, double L, double E){
  Double_t phase = 1.267*deltam2_in*L/E;

  return 1./1.267/deltam2_in *  (0.5*phase - 0.25 * sin(2*phase));

}

Double_t OscCalc::OscProb(double L, double E, double s22t13){

  double theta13=TMath::ASin(sqrt(s22t13))*0.5;
  double theta12=TMath::ASin(sqrt(sin22t12))*0.5;

  //3-neutrino flavor calculation using equation 9 in the TDR
  /*double term1=sin22t12*pow(cos(theta13),4)*pow(sin(1.267*deltam2_21*L/E),2);
    double term2=pow(cos(theta12)*sin(2*theta13)*sin(1.267*deltam2_31*L/E),2);
    double term3=pow(sin(theta12)*sin(2*theta13)*sin(1.267*deltam2_32*L/E),2);
    return 1-(term1+term2+term3); */

  //4-neutrino flavor calculation according to Eq. 6 of DocDB-9296-v1
  double Ue1=cos(theta14)*cos(theta13)*cos(theta12);
  double Ue2=cos(theta14)*cos(theta13)*sin(theta12);
  double Ue3=cos(theta14)*sin(theta13);
  double Ue4=sin(theta14);

  double term1=4*pow(Ue1,2)*pow(Ue2,2)*DeltaTerm(deltam2_21,L,E)+4*pow(Ue1,2)*pow(Ue3,2)*DeltaTerm(deltam2_31,L,E);
  double term2=4*pow(Ue2,2)*pow(Ue3,2)*DeltaTerm(deltam2_32,L,E);
  double term3=4*pow(Ue1,2)*pow(Ue4,2)*DeltaTerm(deltam2_41,L,E)+4*pow(Ue2,2)*pow(Ue4,2)*DeltaTerm(deltam2_42,L,E);
  double term4=4*pow(Ue3,2)*pow(Ue4,2)*DeltaTerm(deltam2_43,L,E);

  return 1-term1-term2-term3-term4;

}

Double_t OscCalc::OscProbInt(double L, double E, double s22t13, int term){

  double theta13=TMath::ASin(sqrt(s22t13))*0.5;
  double theta12=TMath::ASin(sqrt(sin22t12))*0.5;

  //4-neutrino flavor calculation according to Eq. 6 of DocDB-9296-v1
  double Ue1=cos(theta14)*cos(theta13)*cos(theta12);
  double Ue2=cos(theta14)*cos(theta13)*sin(theta12);
  double Ue3=cos(theta14)*sin(theta13);
  double Ue4=sin(theta14);

  double result = 0;

  double term1;
  double term2;
  double term3;
  double term4;
  double term5;
  double term6;

  if (term < 0){
    term1=4*pow(Ue1,2)*pow(Ue2,2)*DeltaTermInt(deltam2_21,L,E);
    term2=4*pow(Ue1,2)*pow(Ue3,2)*DeltaTermInt(deltam2_31,L,E);
    term3=4*pow(Ue2,2)*pow(Ue3,2)*DeltaTermInt(deltam2_32,L,E);
    term4=4*pow(Ue1,2)*pow(Ue4,2)*DeltaTermInt(deltam2_41,L,E);
    term5=4*pow(Ue2,2)*pow(Ue4,2)*DeltaTermInt(deltam2_42,L,E);
    term6=4*pow(Ue3,2)*pow(Ue4,2)*DeltaTermInt(deltam2_43,L,E);
  }else{
    term1=DeltaTermInt(deltam2_21,L,E);
    term2=DeltaTermInt(deltam2_31,L,E);
    term3=DeltaTermInt(deltam2_32,L,E);
    term4=DeltaTermInt(deltam2_41,L,E);
    term5=DeltaTermInt(deltam2_42,L,E);
    term6=DeltaTermInt(deltam2_43,L,E);
  }

  if (term < 0){ // all terms included
    result = L/E-term1-term2-term3-term4-term5-term6;
  }else if (term == 0){
    result = L/E;
  }else if (term == 1){
    result = -term1;
  }else if (term == 2){
    result = -term2;
  }else if (term == 3){
    result = -term3;
  }else if (term == 4){
    result = -term4;
  }else if (term == 5){
    result = -term5;
  }else if (term == 6){
    result = -term6;
  }

  return result;

}

Double_t OscCalc::OscFunc(double L, double s22t13, TH1F *anspec){

  TH1F *anspecosc = (TH1F*)anspec->Clone("anspecosc");
  anspecosc->SetLineColor(2);
  for(int pts=1;pts<=anspec->GetXaxis()->GetNbins();++pts){

    double E=anspec->GetXaxis()->GetBinCenter(pts);
    anspecosc->SetBinContent(pts,anspec->GetBinContent(pts)*OscProb(L,E,s22t13));

  }

  double intrat = anspecosc->Integral()*1./anspec->Integral();
  delete anspecosc;
  return intrat;

}

Double_t OscCalc::OscFunc(double L, double s22t13){

  TH1F *anspecosc = (TH1F*)anspecdef->Clone("anspecosc");
  anspecosc->SetLineColor(2);
  for(int pts=1;pts<=anspecdef->GetXaxis()->GetNbins();++pts){

    double E=anspecdef->GetXaxis()->GetBinCenter(pts);
    anspecosc->SetBinContent(pts,anspecdef->GetBinContent(pts)*OscProb(L,E,s22t13));

  }

  double intrat = anspecosc->Integral()*1./anspecdef->Integral();
  delete anspecosc;
  return intrat;

}

TH1F* OscCalc::OscSpec(double L, double s22t13, TH1F *anspec){

  TH1F *anspecosc = (TH1F*)anspec->Clone("anspecosc");
  anspecosc->SetLineColor(2);
  anspecosc->Reset();//<--should have no effect; just a precaution
  for(int pts=1;pts<=anspec->GetXaxis()->GetNbins();++pts){

    double E=anspec->GetXaxis()->GetBinCenter(pts);
    anspecosc->SetBinContent(pts,anspec->GetBinContent(pts)*OscProb(L,E,s22t13));

  }

  return anspecosc;

}

void OscCalc::OscSpec(double L, double s22t13, TH1F *anspec, TH1F* anspecosc){

  //  anspecosc = (TH1F*)anspec->Clone("anspecosc");
  anspecosc->SetLineColor(2);
  anspecosc->Reset();//<--should have no effect; just a precaution
  for(int pts=1;pts<=anspec->GetXaxis()->GetNbins();++pts){

    double E=anspec->GetXaxis()->GetBinCenter(pts);
    anspecosc->SetBinContent(pts,anspec->GetBinContent(pts)*OscProb(L,E,s22t13));

  }

  // return anspecosc;

}

void OscCalc::OscSpecBinInt(double L, double s22t13, TH1F *anspec, TH1F* anspecosc, int term){ // bin integrated oscillation prediction

  //  anspecosc = (TH1F*)anspec->Clone("anspecosc");
  anspecosc->SetLineColor(2);
  anspecosc->Reset();//<--should have no effect; just a precaution


  int ndiv = 10;// make finer divisions by this number in order to calculate oscillation probability more precisely for high-frequency oscillation search

  for(int pts=1;pts<=anspec->GetXaxis()->GetNbins();++pts){

    double E_binw = anspec->GetXaxis()->GetBinWidth(pts);

    double E_low = anspec->GetXaxis()->GetBinLowEdge(pts);
    double E_high = E_low + E_binw;
    double E_center = 0.5*(E_low+E_high);

    double N_center = anspec->GetBinContent(pts);
    double N_low = 0;
    if (pts > 1) {
      N_low = anspec->GetBinContent(pts-1);
    }
    double N_high = 0;
    if (pts < anspec->GetXaxis()->GetNbins()) {
      N_high = anspec->GetBinContent(pts+1);
    }

    double E_step = (E_high-E_low)/(Double_t)ndiv;
    double N_osc_sum = 0;
    for (int idiv = 0; idiv < ndiv; idiv++){
      double E_low_fine = E_low + E_step*idiv;
      double E_high_fine = E_low + E_step*(idiv+1);
      double E_center_fine = 0.5*(E_low_fine+E_high_fine);
      double N_center_fine = 0;
      double dE = E_center_fine - E_center;

      if (dE < 0){
        N_center_fine = ((E_binw+dE) * N_center - dE * N_low)/E_binw;
      }else{
        N_center_fine = ((E_binw-dE) * N_center + dE * N_high)/E_binw;
      }

      double OscProbAve = 1;


      if (E_low_fine > 0)
        OscProbAve = (OscProbInt(L,E_low_fine,s22t13,term) - OscProbInt(L,E_high_fine,s22t13,term))/(L/E_low_fine - L/E_high_fine);
      else
        OscProbAve = 1;

      N_osc_sum += N_center_fine * OscProbAve;
    }

    anspecosc->SetBinContent(pts,N_osc_sum/(double)ndiv);
  }

  // return anspecosc;

}


double OscCalc::OscSpecSample(double L, double s22t13, TH1F *anspec, double Eneu){

  int ibin = anspec->GetXaxis()->FindBin(Eneu);
  double val=anspec->GetBinContent(ibin);
  val*=OscProb(L,Eneu,s22t13);

  return val;

}

Double_t OscCalc::OscSpecInt(double L, double s22t13, TH1F *anspec){

  TH1F *anspecosc =  OscSpec(L, s22t13, anspec);
  double theint = anspecosc->Integral();
  delete anspecosc;
  return 1e30*theint;

}
