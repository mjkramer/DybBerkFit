#include <iostream>
#include "TObject.h"
#include "TH1.h"
//#include "FluxCalculator.h"
#include "PredSet.h"
//This object holds all the data for one time period (observed events, livetimes, muon efficiencies, backgrounds... etc). 

using namespace std;

class TimePeriodData : public TObject {

public:
  TimePeriodData();
  ~TimePeriodData();
  void CorrectEvts(bool print=false);//<---correct for bg, livetimes and muon efficiencies
  void CorrectSpec(bool print=false);//<---correct for bg, livetimes and muon efficiencies
  void PrintToScreen();
  
  Double_t Livetime[Ndetectors];//<--in days? (doesn't matter as long as use same units)
  Double_t MuonVetoEff[Ndetectors];//<--muon veto efficiency (percent of time that did not cut out?)
  Double_t DMCEff[Ndetectors];
  Double_t ObsEvts[Ndetectors];
  Double_t ErrEvts[Ndetectors];
  Double_t CorrEvts[Ndetectors];
  Double_t BgEvts[Ndetectors];
  Double_t AccEvts[Ndetectors];
  Double_t Li9Evts[Ndetectors];
  Double_t AmcEvts[Ndetectors];
  Double_t AlnEvts[Ndetectors];
  Double_t FnEvts[Ndetectors];
  Double_t AccErr[Ndetectors];
  Double_t Li9Err[Ndetectors];
  Double_t AmcErr[Ndetectors];
  Double_t AlnErr[Ndetectors];
  Double_t FnErr[Ndetectors];
  Double_t CorrBgEvts[Ndetectors];
  Double_t ErrBgEvts[Ndetectors]; // statistical uncertainties for the backgrounds 
  Double_t TargetMass[Ndetectors];
  Double_t BgEvtsLiv[Ndetectors];

  //Note to self: do not define as many histograms for bg as have equivalents for numbers above for simplicity and in order to not have to clone them a lot. 
  TH1F *ObsEvtsSpec[Ndetectors];
  TH1F *CorrBgEvtsSpec[Ndetectors];
  TH1F *CorrAccEvtsSpec[Ndetectors];
  TH1F *CorrLi9EvtsSpec[Ndetectors];
  TH1F *CorrAmcEvtsSpec[Ndetectors];
  TH1F *CorrFnEvtsSpec[Ndetectors];
  TH1F *CorrAlnEvtsSpec[Ndetectors];
  TH1F *CorrEvtsSpec[Ndetectors];
  
private:
  bool firstcorrection;
  int nwarnings;
  int warlimit;//<--limit of consistency warnings to print
  
public:
  ClassDef(TimePeriodData,1);
 
};
