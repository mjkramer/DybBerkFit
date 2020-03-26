#pragma once

#define MAX_XSEC_SAMPLES 1050

class CrossSectionTable {
public:
  CrossSectionTable();
  virtual ~CrossSectionTable();
  // Return the inverse beta decay cross section
  //  Units: interactions cm^2 antineutrino^-1 target-proton^-1
  double inverseBetaDecay(double e_nu);
  int load(const char* filename);
  int loadOneFile(const char* filename, double* xData, double* yData);
  double eMin();
  double eMax();
private:
  double m_eMin; // Minimum Enu energy in spectrum
  double m_eMax; // Maximum Enu energy in spectrum
  double m_nSamples; // Samples in current data
  double m_binWidth; // Convenience variable for resolution
  double m_xsec[MAX_XSEC_SAMPLES]; // Anti-nu spectrum data
};
