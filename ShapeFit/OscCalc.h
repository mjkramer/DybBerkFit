#include <iostream>
#include "TObject.h"
#include "TFile.h"
#include "TH1F.h"
#include "TMath.h"

class OscCalc : public TObject {

public: 
  OscCalc();
  ~OscCalc();
  Double_t OscProb(double L, double E, double s22t13);
  Double_t OscProbInt(double L, double E, double s22t13, int term = -1);
  Double_t OscFunc(double L, double s22t13, TH1F *anspec);
  Double_t OscFunc(double L, double s22t13);
  TH1F *OscSpec(double L, double s22t13, TH1F *anspec);
  void OscSpec(double L, double s22t13, TH1F *anspec, TH1F *anspeosc);
  void OscSpecBinInt(double L, double s22t13, TH1F *anspec, TH1F *anspeosc, int term = -1);
  Double_t OscSpecInt(double L, double s22t13, TH1F *anspec);
  double OscSpecSample(double L, double s22t13, TH1F *anspec, double Eneu);
  void SetDeltaM2_32(double deltam2);
  void SetDeltaM2_31(double deltam2);
  void SetDeltaM2_41(double deltam2);
  void SetDeltaM2_ee(double deltam2);
  void SetDeltaM2_21(double deltam2);
  void SetTheta12(double theta_in);
  void SetTheta14(double theta_in);
  void SetHierarchy(int hie);//1 for NH, -1 for NH
  double GetDeltaM2_32();
  double GetDeltaM2_31();
  double GetDeltaM2_41();
  double GetDeltaM2_ee();
  double GetDeltaM2_21();
  double GetTheta12();
  double GetTheta14();
  //note: default parameters are only for those independent ones (for example, knowing hierarchy and deltam2_32 you can get deltam2_31)
  double deltam2_32_default;
  double deltam2_21_default;
  double deltam2_41_default;
  double deltam2_factor_default;
  double sin22t12_default;
  double theta14_default;
  int hierarchy;
     
private:
  double deltam2_32;
  double deltam2_21;
  double deltam2_31;
  double deltam2_ee;
  double deltam2_41;
  double deltam2_42;
  double deltam2_43;
  double deltam2_factor;
  double sin22t12;
  double theta14;
  TH1F *anspecdef;
  bool loadedtable;
  double DeltaTerm(double deltam2_in,double L, double E);
  double DeltaTermInt(double deltam2_in,double L, double E);

public:
  

public:
  ClassDef(OscCalc,1);

};

// Call the ClassImp macro to give the ABC class RTTI and
// full I/O capabilities.
#if !defined(__CINT__)
ClassImp(OscCalc);
#endif

