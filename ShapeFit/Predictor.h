#ifndef PREDICTOR_H
#define PREDICTOR_H


#include "TObject.h"
#include "TGraph.h"
#include "TFile.h"
#include "TTree.h"
#include "TimePeriodData.h"

//const int Nstage=Nstage;//<--as long as this number is larger than the actual number of periods, you're OK.

const Int_t max_n_enu_bins = 156; // number of true enu bins
const Int_t max_n_evis_bins = 37; // number of evis bins

const Int_t nNearHalls = 2; //number of near halls

class Predictor : public TObject {

public:
  Predictor();
  ~Predictor();
  //general stuff
  void EnterFluxCalculator(FluxCalculator *fluxcalc_in);
  void LoadMainData(const Char_t *mainmatrixname);//<-- if call this have to call LoadBgSpec as well as that is where corrections are applied
  void LoadPredictedIBD(const Char_t *nibdname);
  void LoadIBDSpec(TString *ibdspecname);
  Int_t LoadToyIBDSpec(const Char_t *toyibdspecname);
  void LoadToyMCEntry(Int_t i, bool correct=true);
  void LoadToyMCNominalSpec();
  void LoadBgSpec(TString *accspecname,
                  const Char_t *li9specname,
                  const Char_t *amcspecname,
                  const Char_t *fnspecname,
                  const Char_t *alnspecname);
  FluxCalculator *  GetFluxCalculator(){return fluxcalc;}

  void MakePrediction(double sin2theta13, double dm2, Double_t sin22t14, Double_t dm2_41, int tperiod, TimePeriodData &tperdat); //, PredSet *predout);//<--use dm2=-1 to get default value
  void MakeAllPredictions(double sin2theta13, double dm2=-1, Double_t sin22t14 = 0, Double_t dm2_41 = -1, bool print=false);
  PredSet *MakeOneSuperPrediction(double sin2theta13, double dm2=-1,Double_t sin22t14 = 0,Double_t dm2_41 = -1, bool print=false);
  TimePeriodData *GetTimePeriodData(int istage);
  void CombineData();//this is a sub-function of MakeOneSuperPrediction
  void CombinePredictions();//this is a sub-function of MakeOneSuperPrediction

  //Change number of events
  void EnterObsEvts(double obsad1,double obsad2, double obsad3, double obsad4, double obsad5, double obsad6, double obsad7, double obsad8, int istage);

  //Return data
  TH1F *GetCombCorrEvtsSpec(int idet);
  TH1F *GetCorrEvtsSpec(int istage, int idet);

  TH1F *GetCorrAccEvtsSpec(int istage, int idet);
  TH1F *GetCorrLi9EvtsSpec(int istage, int idet);
  TH1F *GetCorrAmcEvtsSpec(int istage, int idet);
  TH1F *GetCorrFnEvtsSpec(int istage, int idet);
  TH1F *GetCorrAlnEvtsSpec(int istage, int idet);

  void LoadEvisToEnuMatrix(const Char_t *evis_to_enu_matrix_name);

  void SetStatFactor(Double_t fac); // set factor to artificailly increase statitstics, by reducing stat error by 1/sqrt(stat_factor)


  //Chi-squares
  Double_t CalculateChi2Cov(Double_t sin22t13, Double_t dm2_ee, Double_t sin22t14 = 0, Double_t dm2_41 = -1);
  Double_t CalculateChi2Cov();

  Double_t CalculateChi2CovRate(Double_t sin22t13, Double_t dm2_ee, Double_t sin22t14 = 0, Double_t dm2_41 = -1);
  Double_t CalculateChi2CovRate();


  //Covariance matrix stuff
  void LoadCovMatrix(const Char_t *covmatrixname_sig, const Char_t *covmatrixname_bg);
  void AddandScaleCovMatrix(Int_t type = -1);
  void InvertMatrix();
  void InvertRateMatrix();
  void CalculateStatError();
  void CalculateNearSiteStatError();

  void SetStage(Int_t i=-1);

  // Get final covariance matrix for given predictions
  Double_t * GetFinalCovMatrix(Int_t type, Int_t mode = 1);
  Double_t * GetFinalCovMatrixSum(Int_t type, Int_t mode = 1);
  Double_t * GetFinalObs(Int_t mode = 1);
  Double_t * GetFinalObsError(Int_t mode = 1);
  Double_t * GetFinalPred(PredSet *evtset, Int_t mode = 1);

  Double_t * GetFinalObsSum(Int_t mode = 1);
  Double_t * GetFinalObsErrorSum(Int_t mode = 1);
  Double_t * GetFinalPredSum(PredSet *evtset, Int_t mode = 1);

  void SetEvisBins(Int_t n, Double_t * bins, Int_t rebin_fac = 1);
  void SetEnuBins(Int_t n, Double_t * bins);

  Double_t * GetRebinnedEvisBins();
  Int_t GetNumRebinnedEvisBins();

  //Double_t GetWmeanCoeff(Int_t istage, Int_t i, Int_t idet);
  Double_t GetFracStatError(Int_t i, Int_t j);
  Double_t GetFracBgError(Int_t i, Int_t j);


  void FixCovMatrix(Double_t sin22t13, Double_t dm2_ee,Double_t sin22t14, Double_t dm2_41);

public:
  FluxCalculator *fluxcalc;

  Double_t * CombineMatrix(Int_t mode = 1, Bool_t MakeSum = false);

  //object to hold data sets at each time period
  TimePeriodData tdper[Nstage];

  TH1F *CorrEvtsCoreSpec[Nstage][Ndetectors][max_n_evis_bins][Ncores];

  Double_t EvisToEnuMatrix(Int_t ievis, Int_t ienu){return M_evis_to_enu[ievis][ienu];}

  Int_t GetNEvisBins(){return n_evis_bins;}
  Int_t GetNEnuBins(){return n_enu_bins;}

  Double_t* GetEnuBins(){return &enu_bins[0];}


private:
  Double_t LookupOscFunc[Ndetectors][Ndetectors][NpointsS];
  Double_t LookupS[NpointsS];

  //object to hold predictions for each time period
  PredSet * predper;

  //Actual events (after corrections, all combined into one)
  Double_t CombCorrEvts[Ndetectors];
  Double_t CombErrEvts[Ndetectors];
  Double_t CombLivetime[Ndetectors];
  Double_t CombCorrBgEvts[Ndetectors];

  TH1F *PredEvtsTrueSpec[Nstage][Ndetectors][Ndetectors][max_n_evis_bins];
  TH1F *CombCorrEvtsSpec[Ndetectors];
  TH1F *CombCorrBgEvtsSpec[Ndetectors];
  TH1F *CombCorrAccEvtsSpec[Ndetectors];
  TH1F *CombCorrLi9EvtsSpec[Ndetectors];
  TH1F *CombCorrFnEvtsSpec[Ndetectors];
  TH1F *CombCorrAmCEvtsSpec[Ndetectors];

  TH1F *Spec[Nstage][Ndetectors];

  TH1F *CorrEvtsTrueSpec[Nstage][Ndetectors][max_n_evis_bins];
  TH1F *PredEvtsSpec[Nstage][Ndetectors][Ndetectors];

  Double_t M_evis_to_enu[max_n_evis_bins][max_n_enu_bins];

  Double_t M[MaxPredictions*Nstage*max_n_evis_bins][MaxPredictions*Nstage*max_n_evis_bins]; // covariance matrix for the fit
  Double_t M_scaled[MaxPredictions*Nstage*max_n_evis_bins][MaxPredictions*Nstage*max_n_evis_bins]; // scaled covariance matrix for the fit
  Double_t M_inv[MaxPredictions*Nstage*max_n_evis_bins][MaxPredictions*Nstage*max_n_evis_bins]; // inverted covariance matrix for the fit
  Double_t M_rate_inv[MaxPredictions*Nstage][MaxPredictions*Nstage]; // inverted covariance matrix for the rate-only fit

  Double_t M_sig_sys[MaxPredictions*Nstage*max_n_evis_bins][MaxPredictions*Nstage*max_n_evis_bins];
  Double_t M_bg_sys[MaxPredictions*Nstage*max_n_evis_bins][MaxPredictions*Nstage*max_n_evis_bins];
  Double_t M_bg_sys_frac[MaxPredictions*Nstage*max_n_evis_bins][MaxPredictions*Nstage*max_n_evis_bins];
  Double_t M_stat[MaxPredictions*Nstage*max_n_evis_bins][MaxPredictions*Nstage*max_n_evis_bins];
  Double_t M_stat_frac[MaxPredictions*Nstage*max_n_evis_bins][MaxPredictions*Nstage*max_n_evis_bins]; // fractional size of the statistical error

  Double_t mode1_coeff[Nstage][Ndetectors]; // coefficients for taking weighted mean.
  Double_t NearEH_coeff[Nstage][Ndetectors]; // coefficients for taking weighted sum of EH1+EH2
  Double_t FarAD_stage_coeff[Nstage][Ndetectors]; // coefficients for taking weighted sum of Far AD over different periods
  Double_t FarEH_stage_coeff[Nstage][Ndetectors]; // coefficients for taking weighted sum of Far EH over different periods

  Double_t nIBD[Nstage][Ndetectors]; //Store predicted # IBD from ToyMC

  Double_t chi2_map[1000];

  Int_t stage; //0=6AD, 1=8AD, -1=6+8AD period

  Int_t nsteps;
  Double_t s2t_step;
  Bool_t loadedosctable;

  Bool_t BgMatrixScaled;
  Bool_t isMC;

  Int_t combine_mode;


  Char_t dummyname[1024];


  Int_t n_evis_bins ;
  Int_t n_enu_bins ;

  Int_t n_evis_bins_rebin ;
  Double_t evis_bins_rebin[max_n_evis_bins+1];

  Double_t evis_bins[max_n_evis_bins+1];
  Double_t enu_bins[max_n_enu_bins+1];

  Int_t evis_rebin_map[max_n_evis_bins];



  Bool_t FirstToySample;
  Bool_t FirstMakePrediction;
  Bool_t FirstMakeSuperPrediction;
  TH1F *hweightedsum;


  Double_t final_obs[MaxPredictions*Nstage*max_n_evis_bins];
  Double_t final_obserror[MaxPredictions*Nstage*max_n_evis_bins];
  Double_t final_pred[MaxPredictions*Nstage*max_n_evis_bins];
  Double_t final_covmatrix[MaxPredictions*Nstage*max_n_evis_bins*MaxPredictions*Nstage*max_n_evis_bins];

  Double_t final_obs_sum[MaxPredictions*max_n_evis_bins];
  Double_t final_obserror_sum[MaxPredictions*max_n_evis_bins];
  Double_t final_pred_sum[MaxPredictions*max_n_evis_bins];
  Double_t final_covmatrix_sum[MaxPredictions*max_n_evis_bins*MaxPredictions*max_n_evis_bins];

  Double_t stat_factor[Nstage]; // factor to artificailly increase statitstics, by reducing stat error by 1/sqrt(stat_factor)
  Double_t days[Nstage];

  Double_t extra_days_8ad; // extra days for 8AD data taking period to reduce stat error
  Double_t days_6ad; //
  Double_t days_8ad; //
  Double_t days_7ad; //

  TFile *m_toyinfilespec;
  TTree * m_tr_toy;

  PredSet * superpred;

  Bool_t CalculateMode1Coeff;
  Bool_t CalculateNearEHCoeff;
  Bool_t CalculateStageCoeff;
  Bool_t RecalculateCovMatrix;

public:
  ClassDef(Predictor,1);

private:
  friend class OscProbTable;    // Only used by fit_shape_sterile_hybrid_scan
};

#endif
