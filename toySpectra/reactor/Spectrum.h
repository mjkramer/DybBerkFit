#ifndef SPECTRUM_H
#define SPECTRUM_H

#include "Predictor.h"

#define MAX_SAMPLES 3000
//#define MAX_SAMPLES 240
//const int Ncores=6;
//const int Ndetectors=6;

class DataSet;
class OscCalc;
class IsotopeTable;
class CrossSectionTable;
class CoreSpectrum;
class TH1F;
class TF1;
class TRandom3;
class TGraph;

// Non-linearlity function developed  by Liangjiang and Soeren
Double_t nl_func_ihep(Double_t * x, Double_t * par);


double NonLinearity(double depositE,
                   double delta_a1,
                   double delta_a2,
                   double delta_a3,
                   double delta_a4,
                    double delta_a5);

// LBNL non linearity model
Double_t nl_func_lbnl(Double_t * x, Double_t * par); 


// Resolution function
Double_t reso_func_bcw(Double_t * x, Double_t * par);



const Int_t n_bcw_positron_nl = 1000;

double m_lbnl_positron_nl_e[300];
double m_lbnl_positron_nl_fac[300];
double m_lbnl_positron_nl_err[3][300];

const Int_t n_unified_nl_points = 500;

// conversion from neutrino energy to positron energy
  Double_t enu_to_epositron(Double_t * x, Double_t * par); 
  Double_t enu_to_epositron_0th(Double_t * x, Double_t * par); 


class Spectrum {
 public:
  Spectrum();
  virtual ~Spectrum();
  // Return the nuebar spectrum (events / MeV)
  double antiNu(int istage, int idet, double energy);
  // Return the positron true spectrum (events / MeV)
  double positronTrue(int istage, int idet, double energy);
  // Return the positron detected spectrum (events / MeV)
  double positronDetected(int istage, int idet, double energy);
  // Return the energy axis for spectra
  double* energyArray(int idet);
  // Return the nuebar spectrum (events / MeV)
  double* energyArrayBkg(int idet);
  // Return the nuebar spectrum (events / MeV)
  double* antiNuArray(int istage, int idet);
  // Return the nuebar spectrum, no oscillation (events / MeV)
  double* antiNuNoOscArray(int istage, int idet);
  // Return the positron true spectrum (events / MeV)
  double* positronTrueArray(int istage, int idet);
  // Return the positron true spectrum after IAV distortion (events / MeV)
  double* positronIavDistortedArray(int istage, int idet);
  // Return the positron non-linear spectrum (events / MeV)
  double* positronNLArray(int istage, int idet);
  // Return the positron detected spectrum (events / MeV)
  double* positronDetectedArray(int istage, int idet);
  // Return the non-linearity corrected positron detected spectrum (events / MeV)
  double* positronNLCorrectedArray(int istage, int idet);
  // Return the background detected spectrum
  double* bgDetectedArray(int istage, int idet);
  // Return the non-linearity corrected background detected spectrum
  double* bgNLCorrectedArray(int istage, int idet);
  // Initialize the spectrum
  void initialize(DataSet* data);
  // Load distances
  void loadDistances(const char *distancematrixname);
  // load bg spectra (they have different binning to the ones used by ShapeFit, which are loaded by Predictor)
  void loadBgSpecForToy(TString *accspecname, 
			const Char_t *li9specname, 
			const Char_t *amcspecname,
			const Char_t *fnspecname,
			const Char_t *alnspecname);
  void setBgRemoveFlag(bool acc_flag, bool li9_flag, bool fn_flag, 
                       bool amc_flag, bool aln_flag);
  
  // Update Oscillation
  //void setDeltaMSqee(double deltaMSqee);
  void setOscillation(double sinSq2Th12, double sinSq2Th13, 
                      double deltaMSq21, double deltaMSqee);
  // Update the energy scale
  void setEnergyScale(double alpha, double beta);
  // Number of samples
  int nSamples();
  // Bin width (used for normalizatoin)
  double binWidth();
  // Oscillation parameter
  int nSamplesBkg();
  // Bin width (used for normalizatoin)
  double binWidthBkg();
  // Oscillation parameter
  double deltaMSqee();
  // Rebin spectrum to lower resolution -> speed up calculation
  void rebin(int nCombine);
  // Generate MC data histogram with statistical uncertainties
  //  TH1F* monteCarlo();
  void setRandomSeed(unsigned int i);

  // Update internal data
  //void updateOscillation();
  void updateAntinu(int icore_select = -1);
  void updatePositronTrue(double eNu_min = -1, double eNu_max = -1);
  void updatePositronDetected();
  void correctNonLinearity();
  void updateBgDetected();
  
  // Pass Predictor object (containing input data, spectra and other goodies)
  void passPredictor(Predictor* pred_in);

  // Distorsion curves
  TF1* getDistortionCurve(double amount);
  TF1 *getFNDistortionCurve(double amount);

  TF1 *amcfunc;//<--function with AmC background shape


  
 private:
  
  //Bg spectra
  TH1F *CorrAccEvtsSpec[Nstage][Ndetectors];
  TH1F *CorrAmcEvtsSpec[Nstage][Ndetectors];
  TH1F *CorrAlnEvtsSpec[Nstage][Ndetectors];
  TH1F *CorrLi9EvtsSpec[Nstage][Ndetectors];
  TH1F *CorrFnEvtsSpec[Nstage][Ndetectors];
  
  void setRandomSolarOscPars(); 
  void setRandomDm2ee();

  void setRandomReactorPower(); 
  void setRandomIavDistortion();
  void setRandomScintiNonLinear(); 
  void setRandomIhepScintiNonLinear(); 
  void setRandomBcwScintiNonLinear(); 
  void setRandomLbnlScintiNonLinear(); 
  void setRandomUnifiedScintiNonLinear(); 
  void setRandomElecNonLinear(); 
  void setRandomResolution();  
  void setRandomRelativeEnergyScale(); 
  void setRandomAbsoluteEnergyScale(); 
  void setRandomRelativeEnergyOffset(); 
  void setRandomAbsoluteEnergyOffset(); 

  void setRandomDetectorEfficiency();
    
 private:
  double m_eMin; // Minimum energy in spectra
  double m_eMax; // Maximum energy in spectra
  int m_nSamples; // Samples in current data
  double m_binWidth; // Convenience variable for spectra sampling

  int m_nSamplesBkg; // Samples in current data
  double m_binWidthBkg; // Convenience variable for spectra sampling
  
  double m_detectorDistance[Ndetectors][Ncores]; // Distance from reactor to detector
  double m_alpha; // Detector energy non-linearity parameter
  double m_alpha_nominal; // Detector energy non-linearity parameter
  double m_alpha_error; // Detector energy non-linearity parameter

  double m_beta; // Detector energy linearity parameter 
  double m_beta_nominal; // Detector energy linearity parameter 

  double m_detectorResolution; // Detector energy resolution parameter
  double m_detectorResolution_nominal; // Detector energy resolution parameter
  double m_detectorResolution_error; // Detector energy resolution parameter
  double m_detectorResolution_error_uncorr; // Detector energy resolution parameter
  double m_detectorResolution_bias[Ndetectors]; // Random biases of the detector resolution.

  //double m_detectorEfficiency;// Detector efficiency
  //double m_runningTime;
  //double m_detectorSize;

  double m_detectorEfficiency[Nstage][Ndetectors];// Detector efficiency
  double m_detectorEfficiency_dmc[Nstage][Ndetectors];// Detector efficiency for DMC cut
  double m_detectorEfficiency_muveto[Nstage][Ndetectors];// Detector efficiency for muon veto

  // other detector efficiency, common for all detectors
  double m_detectorEfficiency_Dt;
  double m_detectorEfficiency_Ep;
  double m_detectorEfficiency_Ed_nominal;
  double m_detectorEfficiency_flash;
  double m_detectorEfficiency_nGd;
  double m_detectorEfficiency_spill;

  double m_detectorEfficiency_Ed[Ndetectors]; // can be different depending on relative energy scale shift

  double m_runningTime[Nstage][Ndetectors];
  double m_detectorSize[Ndetectors];

  
  double m_targetProtonsPerKton;
  double m_nominalReactorPower[Nstage][Ncores];
  double m_reactorPower[Nstage][Ncores];
  double m_reactorPowerError[Nstage][Ncores];
  double m_sinSq2Th12;
  double m_sinSq2Th14;  
  double m_sinSq2Th13;
  double m_deltaMSq21;
  double m_deltaMSq41;
  double m_deltaMSqee;

  double m_sinSq2Th12_nominal;
  double m_sinSq2Th12_error;
  double m_sinSq2Th14_nominal;
  double m_sinSq2Th14_error;

  double m_deltaMSqee_nominal;
  double m_deltaMSqee_error;
  double m_deltaMSq21_nominal;
  double m_deltaMSq21_error;
  double m_deltaMSq41_nominal;
  double m_deltaMSq41_error;
  
  double m_energy[Ndetectors][MAX_SAMPLES];
  double m_energy_bkg[Ndetectors][MAX_SAMPLES];
  double m_antiNuSpectrumNoOsc[Nstage][Ndetectors][MAX_SAMPLES];
  double m_survivalProbability[Nstage][Ndetectors][MAX_SAMPLES];
  double m_antiNuSpectrum[Nstage][Ndetectors][MAX_SAMPLES];
  double m_positronTrueSpectrum[Nstage][Ndetectors][MAX_SAMPLES];
  double m_positronIavDistortedSpectrum[Nstage][Ndetectors][MAX_SAMPLES];
  double m_positronNLSpectrum[Nstage][Ndetectors][MAX_SAMPLES];
  double m_positronDetectedSpectrum[Nstage][Ndetectors][MAX_SAMPLES];
  double m_positronNLCorrectedSpectrum[Nstage][Ndetectors][MAX_SAMPLES];
  double m_bgDetectedSpectrum[Nstage][Ndetectors][MAX_SAMPLES];
  double m_bgNLCorrectedSpectrum[Nstage][Ndetectors][MAX_SAMPLES];

  double m_iav_corr_orig[240][240]; // original IAV correction
  double m_iav_frac_orig[240];

  double m_nominal_iav_corr[MAX_SAMPLES][MAX_SAMPLES];
  double m_nominal_iav_frac[MAX_SAMPLES];

  double m_iav_error; // relative uncertainty of the IAV thickness
  double m_iav_corr[Ndetectors][MAX_SAMPLES][MAX_SAMPLES];
  double m_iav_frac[Ndetectors][MAX_SAMPLES];


  double m_abs_escale;
  double m_abs_escale_nominal;
  double m_abs_escale_error;
  double m_rel_escale[Ndetectors];
  double m_rel_escale_nominal[Ndetectors];
  double m_rel_escale_error[Ndetectors];
  

  double m_useIhepNonLinearModel;
  double m_useBcwNonLinearModel;
  double m_useLbnlNonLinearModel;
  double m_useUnifiedNonLinearModel;
  double m_use2015NonLinearModel;
  double m_correlateUnifiedNonlinearPars;


  double m_ihep_nl_par[5];
  double m_ihep_nl_par_nominal[5];
  double m_ihep_nl_par_error[5];
  double m_ihep_nl_par_covmatrix[5][5];
  double m_ihep_nl_par_covmatrix_l[5][5];
  
  double m_bcw_nl_par[5];
  double m_bcw_nl_par_nominal[5];
  double m_bcw_nl_par_error[5];
  double m_bcw_nl_par_covmatrix[5][5];
  double m_bcw_nl_par_covmatrix_l[5][5];
  

  double m_bcw_elec_nl_par[5];
  double m_bcw_elec_nl_par_nominal[5];
  double m_bcw_elec_nl_par_error[5];
  
  
  double m_lbnl_nl_par[3];
  double m_lbnl_nl_par_nominal[3];
  double m_lbnl_nl_par_error[3];


  double m_unified_nl_par[10];
  double m_unified_nl_par_nominal[10];
  double m_unified_nl_par_error[10];

  
  double m_detectorEfficiency_rel_error;

  
  double m_rel_eoffset[Nstage][Ndetectors];
  double m_abs_eoffset;
  double m_rel_eoffset_error;
  double m_abs_eoffset_error;
  
  
  OscCalc* m_osccalc;
  //  IsotopeTable*  m_isotopes;
  IsotopeTable*  m_isotopes[Nstage][Ncores]; // Isotope tables for each different cores
  CrossSectionTable*  m_xsec;


  CoreSpectrum* m_corespectrum[Nstage];
  
  double m_randomizeSolarOscPars;
  double m_randomizeDm2ee;

  double m_randomizeIsotopeFraction;
  double m_randomizeReactorPower;
  double m_randomizeIavDistortion;
  double m_randomizeScintiNonLinear;
  double m_randomizeElecNonLinear;
  double m_randomizeResolution;
  double m_randomizeRelativeEnergyScale;
  double m_randomizeAbsoluteEnergyScale;
  double m_randomizeRelativeEnergyOffset;
  double m_randomizeAbsoluteEnergyOffset;

  double m_randomizeDetectorEfficiency;

  double m_randomizeCoreSpectra;

  double m_varyAccBg;
  double m_varyAmcBg;
  double m_varyFnBg;
  double m_varyLi9Bg;
  double m_varyAlnBg;

  double m_distortAccBg;
  double m_distortAmcBg;
  double m_distortFnBg;
  string m_distortLi9Bg;
  TFile *m_file_distortLi9Bg;
  TTree *m_tree_distortLi9Bg;
  TH1F *m_func_distortLi9Bg;
  int m_entries_distortLi9Bg;
  double m_distortAlnBg;

  double m_statisticalFluctuation;


  bool m_removeAccBg;
  bool m_removeLi9Bg;
  bool m_removeFnBg;
  bool m_removeAmcBg;
  bool m_removeAlnBg;


  void extractPredictorData();
  void loadIavCorrection(const char *iavcorrectionname);
  
  TRandom3 * ran;
  Predictor *pred;

  
  // Detector response non-linearlity function
  TF1 * nl_func;


  // Detector resolution function
  TF1 * reso_func;

  TF1 * enu_to_epositron_func;

  // Non-linearlity function developed  by Xin
  Double_t nl_func_bcw(Double_t * x, Double_t * par);
  double m_bcw_positron_nl_e[n_bcw_positron_nl];
  double m_bcw_positron_nl_fac[n_bcw_positron_nl];
  
  TGraph* g_bcw_elec_nl_error[2];

  // Unified non linearity model
  Double_t nl_func_unified(Double_t * x, Double_t * par); 

  int m_num_unified_nl_pars;
  
  double m_unified_positron_nl_e[n_unified_nl_points];
  double m_unified_positron_nl_fac[n_unified_nl_points];
  double m_unified_positron_nl_err[10][n_unified_nl_points];
  
  double m_unified_nl_par_covmatrix[10][10];
  double m_unified_nl_par_covmatrix_l[10][10];

  
  bool UseChristineModel;

  double  m_useBcwFluxUncertainty;
  double  m_useAbInitioSpectra;

  
};

#endif // SPECTRUM_H


