#pragma once

const int MAX_ISOTOPE_ID = 5;
const int MAX_ISOTOPE_SAMPLES = 1000;

class TRandom3;

class IsotopeTable {
public:
  IsotopeTable();
  virtual ~IsotopeTable();
  // Return the antineutrino energy spectrum for this fission isotope
  //  Units: antineutrinos MeV^-1 fission^-1
  bool isActive(unsigned int isotopeId);
  double meanEnergyPerFission(unsigned int isotopeId);
  double nominalFissionFraction(unsigned int isotopeId);
  double FissionFraction(unsigned int isotopeId);
  double antiNuSpectrum(unsigned int isotopeId, double e_nu);
  void loadIsotopes(const char* filename);
  void loadSpectra(const char* filename);
  void loadIsotopeCovMatrix(const char* filename);
  void setRandomSeed(unsigned int seed);
  void setRandomFissionFractionCorr();
  void setRandomFissionFraction();
  double eMin();
  double eMax();
private:
  bool m_isActive[MAX_ISOTOPE_ID]; // Is isotope loaded and active?
  double m_meanEnergyPerFission[MAX_ISOTOPE_ID]; // Mean energy per fission
  double m_nominalFissionFraction[MAX_ISOTOPE_ID]; // Nominal fission fraction
  double m_FissionFraction[MAX_ISOTOPE_ID]; // Place for random fission fraction
  double m_eMin; // Minimum Enu energy in spectrum
  double m_eMax; // Maximum Enu energy in spectrum
  double m_nSamples; // Samples in current data
  double m_binWidth; // Convenience variable for resolution
  double m_dNdE[MAX_ISOTOPE_ID][MAX_ISOTOPE_SAMPLES]; // Anti-nu spectrum data

  double m_FissionFractionErr1D[MAX_ISOTOPE_ID]; // 1 dimensional relative error of the fission fraction
  double m_FissionFractionErr[MAX_ISOTOPE_ID][MAX_ISOTOPE_ID]; // Covariance matrix for the fission fraction
  int m_nDimCovMatrix; // dimension of the covariance matrix

  double L[MAX_ISOTOPE_ID][MAX_ISOTOPE_ID];

  bool isIsotopeTableLoaded;
  bool isIsotopeCovMatrixLoaded;
  TRandom3 * ran ;
};
