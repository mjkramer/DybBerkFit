#ifndef CORESPECTRUM_H
#define CORESPECTRUM_H

#include "Predictor.h"
#include "CrossSectionTable.h"
//#include "IsotopeTable.h"
const int MAX_CORESPECTRA_SAMPLES = 240;
const int MAX_ABINITIO_CORESPECTRA_SAMPLES = 12000;
const  int m_nBcwBins = 26;

class TRandom3;

class CoreSpectrum {
 public:
  CoreSpectrum();
  virtual ~CoreSpectrum();
  // Return the antineutrino energy spectrum for this fission isotope
  //  Units: antineutrinos MeV^-1 fission^-1
  double antiNuSpectrum(unsigned int isotopeId, double e_nu);
  double IBDSpectrum(unsigned int isotopeId, double e_nu);
  Bool_t loadSpectra(const char* filename_nom,const char* filename_mcov);
  Bool_t loadAbInitioSpectra(const char* filename_nom);
  Bool_t setFlatSpectra();
  void setRandomSeed(unsigned int seed);
  void setRandomAntiNuSpectra();

  Bool_t loadCovMatrixBCW(const char* filename_mcov);
  void setRandomIBDSpectraBCW();
  double IBDSpectrumBCW(unsigned int coreId, double e_nu);
  
  double eMin();
  double eMax();
 private:
  double m_eMin; // Minimum Enu energy in spectrum
  double m_eMax; // Maximum Enu energy in spectrum
  int m_nSamples; // Samples in current data
  double m_binWidth; // Convenience variable for resolution
  /* double m_dNdE[Ncores][MAX_CORESPECTRA_SAMPLES]; // Anti-nu spectrum data */
  /* double m_dNdE_nom[Ncores][MAX_CORESPECTRA_SAMPLES]; // Anti-nu spectrum data */
  double m_dNdE[Ncores][MAX_ABINITIO_CORESPECTRA_SAMPLES]; // Anti-nu spectrum data
  double m_dNdE_nom[Ncores][MAX_ABINITIO_CORESPECTRA_SAMPLES]; // Anti-nu spectrum data
  double m_dNdE_mcov[Ncores*MAX_CORESPECTRA_SAMPLES][Ncores*MAX_CORESPECTRA_SAMPLES]; // Covariance matrix for Anti-nu spectrum data
  double L[Ncores*MAX_CORESPECTRA_SAMPLES][Ncores*MAX_CORESPECTRA_SAMPLES]; // lower triangle of the covariance matrix

  double m_dIBDdE_nom[Ncores][m_nBcwBins]; // IBD spectra data with bcw bins
  double m_dIBDdE[Ncores][m_nBcwBins]; // IBD spectra data with bcw bins

  double L_BCW[Ncores*m_nBcwBins][Ncores*m_nBcwBins]; // lower triangle of the covariance matrix
  double m_bcw_bincenter[m_nBcwBins];
  double m_bcw_bins[m_nBcwBins+1];

  //  double m_dNdExXsec[Ncores][MAX_CORESPECTRA_SAMPLES]; // Anti-nu spectrum x xsec

  bool isCoreSpectrumLoaded;
  bool isBcwCovMatrixLoaded;

  TRandom3 * ran ;

  CrossSectionTable * m_xsec;
  bool loadCrossSectionTable(const char* filename_xsec);
  
  bool isAbInitioSpectraUsed;
  
};

#endif // CORESPECTRUM_H
