#include "CrossSectionTable.h"

#include <string>
#include <iostream>
#include <fstream>

CrossSectionTable::CrossSectionTable()
: m_eMin(0),
  m_eMax(0),
  m_nSamples(0),
  m_binWidth(1.0)
{
  for(int idx=0; idx<MAX_XSEC_SAMPLES; idx++){
    m_xsec[idx] = 0;
  }
}

CrossSectionTable::~CrossSectionTable()
{
  ;
}

double CrossSectionTable::inverseBetaDecay(double e_nu)
{
  // Return the inverse beta decay cross section
  //  Units: interactions cm^2 antineutrino^-1 target-proton^-1
  if(e_nu<m_eMin || e_nu>=m_eMax) return 0.0;
  // By default, use nearest sample point
  // FIXME: could consider alternate interpolations or splines,
  //        but must measure impact on speed
  //int binIdx = (int)((e_nu-m_eMin)/m_binWidth);
  // By default, use linear interpolation between nearest sample points
  int binIdxLow = (int)((e_nu-m_eMin)/m_binWidth);
  double value = m_xsec[binIdxLow];
  static bool linearInterp = true;
  if(linearInterp){
    int binIdxHigh = binIdxLow+1;
    double binLowE = m_eMin + binIdxLow*m_binWidth;
    double dE = e_nu - binLowE;
    double slope = (m_xsec[binIdxHigh] 
		    - m_xsec[binIdxLow])/m_binWidth;
    value += dE*slope;
  }
  return value;
}

int CrossSectionTable::load(const char* filename)
{
  double xData[MAX_XSEC_SAMPLES];
  double yData[MAX_XSEC_SAMPLES];
  int nSamples = this->loadOneFile(filename,xData,yData);
  if(nSamples<=0) return -1;

  // Check uniform binning of data
  double eMin = xData[0];
  double eMax = xData[nSamples-1];
  double binWidth = xData[1]-xData[0];
  double delta = binWidth - (eMax-eMin)/nSamples;
  double epsilon = 0.0000001;
  if(delta*delta > epsilon){
    // Doesn't appear to be uniform binning
    std::cout << "Error: binning not uniform in " << filename << std::endl;
    return -1;
  }
  // Copy data into cross-section array
  m_nSamples = nSamples;
  m_eMin = eMin;
  m_eMax = eMax;
  m_binWidth = (eMax-eMin)/(nSamples-1);
  for(int idx=0; idx<nSamples; idx++){
    m_xsec[idx] = yData[idx];
  }
  return 0;
}

int CrossSectionTable::loadOneFile(const char* filename, double* xData, 
				   double* yData)
{
  // Load data from one file
  ifstream fileData(filename);
  if(!fileData.is_open() || !fileData.good()){
    std::cout << "Error: Failed to open data file " << filename << std::endl;
    return -1; // FAILURE
  }
  std::string line;
  double xTmp=0;
  double yTmp=0;
  unsigned int curSample = 0;
  while(MAX_XSEC_SAMPLES > curSample){
    if(fileData.peek()=='#'){
      // Skip lines starting with '#'
      getline(fileData,line);
      continue;
    }else{
      fileData >> xTmp >> yTmp;
    }
    if(!fileData.good()) break;
    xData[curSample] = xTmp;
    yData[curSample] = yTmp;
    curSample++;
  }
  return curSample;
}

double CrossSectionTable::eMin()
{
  return m_eMin;
}

double CrossSectionTable::eMax()
{
  return m_eMax;
}

