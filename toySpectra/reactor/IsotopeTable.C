#include "IsotopeTable.h"

#include <string>
#include <iostream>
#include <fstream>

#include <TRandom.h>
#include <TRandom3.h>
#include <TMatrixD.h>
#include <TDecompChol.h>


IsotopeTable::IsotopeTable()
  : m_eMin(0),
    m_eMax(0),
    m_nSamples(0),
    m_binWidth(1.0)
{
  for(int isoId=0; isoId<MAX_ISOTOPE_ID; isoId++){
    m_isActive[isoId] = false;
    m_meanEnergyPerFission[isoId] = 0;
    m_nominalFissionFraction[isoId] = 0;
    m_FissionFractionErr1D[isoId] = 0;

    for(int idx=0; idx<MAX_ISOTOPE_SAMPLES; idx++){
      m_dNdE[isoId][idx] = 0;
    }
  }
  isIsotopeTableLoaded = false;
  isIsotopeCovMatrixLoaded = false;
  ran = new TRandom3();


}

IsotopeTable::~IsotopeTable()
{
  ;
}

bool IsotopeTable::isActive(unsigned int isotopeId)
{
  return m_isActive[isotopeId];
}

double IsotopeTable::meanEnergyPerFission(unsigned int isotopeId)
{
  return m_meanEnergyPerFission[isotopeId];
}

double IsotopeTable::nominalFissionFraction(unsigned int isotopeId)
{
  return m_nominalFissionFraction[isotopeId];
}

double IsotopeTable::FissionFraction(unsigned int isotopeId)
{
  return m_FissionFraction[isotopeId];
}


double IsotopeTable::antiNuSpectrum(unsigned int isotopeId, double e_nu)
{
  // Return the anti-neutrino production
  // Units: [neutrinos MeV^-1 fission^-1]
  if(!m_isActive[isotopeId]) return 0.0;
  if(e_nu<m_eMin || e_nu>=m_eMax) return 0.0;
  // By default, use linear interpolation between nearest sample points
  int binIdxLow = (int)((e_nu-m_eMin)/m_binWidth);
  double value = m_dNdE[isotopeId][binIdxLow];
  static bool linearInterp = true;
  if(linearInterp){
    int binIdxHigh = binIdxLow+1;
    double binLowE = m_eMin + binIdxLow*m_binWidth;
    double dE = e_nu - binLowE;
    double slope = (m_dNdE[isotopeId][binIdxHigh]
                    - m_dNdE[isotopeId][binIdxLow])/m_binWidth;
    value += dE*slope;
  }
  return value;
}

void IsotopeTable::loadIsotopes(const char* filename)
{
  // Load isotope data file
  ifstream fileData(filename);
  if(!fileData.is_open() || !fileData.good()){
    std::cout << "IsotopeTable::loadIsotopes: "
              << "Error: Failed to open data file " << filename << std::endl;
    return;
  }
  std::string line;
  unsigned int isotopeId=0;
  double meanEnergyPerFission=0;
  double nominalFissionFraction=0;
  while(true){
    if(fileData.peek()=='#'){
      // Skip lines starting with '#'
      getline(fileData,line);
      continue;
    }else{
      fileData >> isotopeId >> meanEnergyPerFission >> nominalFissionFraction;
    }
    if(!fileData.good()) break;
    m_isActive[isotopeId]=true;
    m_meanEnergyPerFission[isotopeId]=meanEnergyPerFission;
    m_nominalFissionFraction[isotopeId]=nominalFissionFraction;
    m_FissionFraction[isotopeId]=nominalFissionFraction;
  }
  isIsotopeTableLoaded = true;
  return;
}

void IsotopeTable::loadIsotopeCovMatrix(const char* filename)
{

  if (!isIsotopeTableLoaded){
    cout << "you must load isotope table first!!! " << endl;
    return;
  }

  // Load isotope data file
  ifstream fileData(filename);
  if(!fileData.is_open() || !fileData.good()){
    std::cout << "IsotopeTable::loadIsotopeCovMatrix: "
              << "Error: Failed to open data file " << filename << std::endl;
    return;
  }
  std::string line;
  double relErrors[MAX_ISOTOPE_ID];

  // First line should be number of isotopes
  while(true){
    if(fileData.peek()=='#'){
      // Skip lines starting with '#'
      getline(fileData,line);
      continue;
    }else{
      fileData >> m_nDimCovMatrix;
      getline(fileData,line);
      cout << "\tDimension of the covariance matrix: " << m_nDimCovMatrix << endl;
      break;
    }
  }

  // Second line should be an array of the relative errors
  while(true){
    if(fileData.peek()=='#'){
      // Skip lines starting with '#'
      getline(fileData,line);
      continue;
    }else{
      for (int i = 0; i < m_nDimCovMatrix; i++)
        fileData >> relErrors[i];
      getline(fileData,line);
      cout << "\tRelative errors: ";
      for (int i = 0; i < m_nDimCovMatrix; i++)
        cout << " " <<  relErrors[i];
      cout << endl;

      break;
    }
  }
  for (int i = 0; i < m_nDimCovMatrix; i++){
    m_FissionFractionErr1D[i+1] = relErrors[i];
  }

  // Finally, we read the covariance matrix
  while(true){
    if(fileData.peek()=='#'){
      // Skip lines starting with '#'
      getline(fileData,line);
      continue;
    }else{
      for (int i = 0; i < m_nDimCovMatrix; i++)
        for (int j = 0; j < m_nDimCovMatrix; j++)
          fileData >> m_FissionFractionErr[i][j];
      cout << "\tCovariance matrix:" << endl;
      for (int i = 0; i < m_nDimCovMatrix; i++){
        cout << "\t";
        for (int j = 0; j < m_nDimCovMatrix; j++){
          cout << "\t" <<  m_FissionFractionErr[i][j];
        }
        cout << endl;
      }
      break;

    }
  }

  // since the fission fraction table start from 1 though 4, a special care is taken
  for (int i = 0; i < m_nDimCovMatrix; i++){
    for (int j = 0; j < m_nDimCovMatrix; j++){
      m_FissionFractionErr[i][j]
        = m_FissionFractionErr[i][j] * relErrors[i]*relErrors[j]
        * m_nominalFissionFraction[i+1] * m_nominalFissionFraction[j+1];
    }
  }


  TMatrixD covmatrix(MAX_ISOTOPE_ID,MAX_ISOTOPE_ID,&m_FissionFractionErr[0][0]);

  covmatrix.ResizeTo(m_nDimCovMatrix,m_nDimCovMatrix);
  //  covmatrix.Print();
  // //  ematrix_fakedata.Print();
  // Double_t tmpvec[NBins];
  // TMatrixD ranvec(NBins,1);
  // TMatrixD parvec(NBins,1);

  TDecompChol chol(covmatrix);
  chol.Decompose();

  TMatrixD cmat(chol.GetU());
  //  cmat.Print();
  TMatrixD tcmat(cmat.Transpose(cmat));
  // std::cout << "Cholesky matrix---------" << std::endl;
  // tcmat.Print();

  // Finally, add zeros to column and ring 0, to match the fission fraction data

  double * tmp_matrix = tcmat.GetMatrixArray();

  for (int i = 0; i < MAX_ISOTOPE_ID; i++){
    for (int j = 0; j < MAX_ISOTOPE_ID; j++){
      if (i == 0 || j == 0)
        L[i][j] = 0;
      else
        L[i][j] = tmp_matrix[(i-1)*m_nDimCovMatrix + j-1];
      //      cout << "\t" << L[i][j] ;
    }
    //    cout << endl;
  }



  isIsotopeCovMatrixLoaded = true;

  return;
}

void IsotopeTable::setRandomSeed(unsigned int seed)
{
  ran->SetSeed(seed);
}

void IsotopeTable::setRandomFissionFractionCorr()
{
  // set correlated fission fraction using the covariance matrix
  double ranvec[MAX_ISOTOPE_ID];
  for (int i = 0; i < MAX_ISOTOPE_ID; i++){
    ranvec[i] = ran->Gaus(0,1);
    //    cout << "\t" << ranvec[i];
  }
  //  cout << endl;

  double sum_nominal = 0;
  double sum_new = 0;

  for (int i = 0; i < MAX_ISOTOPE_ID; i++){
    m_FissionFraction[i] = m_nominalFissionFraction[i];
    for (int j = 0; j < MAX_ISOTOPE_ID; j++){
      m_FissionFraction[i] += L[i][j] * ranvec[j];
      //      cout << i << " " << j << " " << L[i][j] << " " << ranvec[j] << endl;
    }
    sum_nominal += m_nominalFissionFraction[i];
    sum_new += m_FissionFraction[i];
  }
  for (int i = 0; i < MAX_ISOTOPE_ID; i++){ // renormalize
    m_FissionFraction[i] = m_FissionFraction[i]*sum_nominal/sum_new;
  }

  // for (int i = 0; i < MAX_ISOTOPE_ID; i++){
  //   cout << i << "\t" << m_nominalFissionFraction[i] << "\t" << m_FissionFraction[i] << endl;
  // }
  // cout << "Sum: \t" << sum_nominal << "\t" << sum_new << endl;
  // cout << endl;


}

void IsotopeTable::setRandomFissionFraction()
{
  // set random fission fraction assuming they are uncorrelated (i.e. not using covariance matrix.

  double sum_nominal = 0;
  double sum_new = 0;

  for (int i = 0; i < MAX_ISOTOPE_ID; i++){
    m_FissionFraction[i]
      = (1 + ran->Gaus(0,1) * m_FissionFractionErr1D[i])
      * m_nominalFissionFraction[i];
  }

}



void IsotopeTable::loadSpectra(const char* filename)
{
  // Load isotope spectra data file
  ifstream fileData(filename);
  if(!fileData.is_open() || !fileData.good()){
    std::cout << "IsotopeTable::loadSpectra: "
              << "Error: Failed to open data file " << filename << std::endl;
    return;
  }
  std::string line;
  double eMin=0;
  double eMax=0;
  double eNu=0;
  double dNdE_U235=0;
  double dNdE_U238=0;
  double dNdE_Pu239=0;
  double dNdE_Pu241=0;
  unsigned int curSample = 0;
  double binWidth = 0;
  while(MAX_ISOTOPE_SAMPLES > curSample){
    if(fileData.peek()=='#'){
      // Skip lines starting with '#'
      getline(fileData,line);
      continue;
    }else{
      fileData >> eNu >> dNdE_U235 >> dNdE_U238 >> dNdE_Pu239 >> dNdE_Pu241;
    }
    if(!fileData.good()) break;
    if(curSample==0){
      eMin=eNu;
      eMax=eNu;
    }
    if(eNu>eMax) eMax=eNu;
    if(curSample==1) binWidth = eMax-eMin;
    m_dNdE[1][curSample] = dNdE_U235;
    m_dNdE[2][curSample] = dNdE_U238;
    m_dNdE[3][curSample] = dNdE_Pu239;
    m_dNdE[4][curSample] = dNdE_Pu241;
    curSample++;
  }
  // Check uniform binning
  double epsilon = 0.0000001;
  double delta = binWidth - (eMax-eMin)/curSample;
  if(delta*delta > epsilon){
    // Doesn't appear to be uniform binning
    std::cout << "IsotopeTable: Error: binning not uniform in " << filename
              << std::endl;
    return;
  }
  m_eMin = eMin;
  m_eMax = eMax;
  m_nSamples = curSample;
  m_binWidth = binWidth;
  return;
}

double IsotopeTable::eMin()
{
  return m_eMin;
}

double IsotopeTable::eMax()
{
  return m_eMax;
}
