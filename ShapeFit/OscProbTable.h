#include <iostream>
#include <string.h>
#include "TObject.h"
#include "TH1.h"
#include "Predictor.h"
#include "Ranger.h"

const int max_energy_bins=240;
const int max_evis_bins=40;
const int max_steps_dmee=1;
const int max_steps_dm41=1001;
const int num_cores = 6;
const int num_terms = 7;


class OscProbTable : public TObject {

public:
  OscProbTable(int method = 1);
  OscProbTable(Predictor *pred_in, int method = 1);
  ~OscProbTable();
  void passPredictor(Predictor *pred_in);
  void SetDMeeRange(int nsteps, double min, double max);
  void SetDM41Range(int nsteps, double min, double max, bool logscale = false);
  void SetDM41Range2(int nsteps, double min, double max, bool logscale = false);
  //  TH1F *PlotExtrapolationHist(int idet_far, int idet_near, int is22t13, int is22t14, int idmee, int idm41);
  void MakeAllPredictionsQuick(double s22t13, double dmee, double s22t14=-1, double dm41=-1, bool print=false);//so it's not called the same as in Predictor, but it's definitely a copy except for the difference of making things using the extrap table; note order of arguments is same as there, to avoid confusion
  void MakeOneSuperPredictionQuick(double sin2theta13, double dm2, Double_t s22t14=-1, Double_t dm2_41=-1, bool print=false);//same comment as for function above
  double CalculateChi2CovQuick(Double_t sin22t13, Double_t dm2_ee,Double_t sin22t14=-1, Double_t dm2_41=-1, Int_t stage=-1);

  void GenerateTable();

  void WriteTable(const char * outfilename);
  void ReadTable(const char * infilename);

  void SetMethod(int method);
  
  
private:
  Predictor *pred;

  bool generatedTable;

  Ranger *ranger_dmee;
  Ranger *ranger_dm41;
  Ranger *ranger_dm41_2; // second dm2 ranger
    
  
  //humongous array with all the info 
  double OscProbTableLookup[Nstage][Ndetectors][max_steps_dmee][max_steps_dm41][max_energy_bins][num_terms];
  double M_evis_to_enu[max_evis_bins][max_energy_bins];
  
  PredSet *m_superpred;

  int m_method;


public:
  ClassDef(OscProbTable,1);

};

// Call the ClassImp macro to give the ABC class RTTI and
// full I/O capabilities.
#if !defined(__CINT__)
ClassImp(OscProbTable);
#endif
