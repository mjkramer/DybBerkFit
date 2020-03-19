#include "CoreSpectrum.h"

#include <string>
#include <iostream>
#include <fstream>

#include <TRandom.h>
#include <TRandom3.h>
#include <TMatrixD.h>
#include <TDecompChol.h>


CoreSpectrum::CoreSpectrum()
  : m_eMin(0),
    m_eMax(0),
    m_nSamples(0),
    m_binWidth(1.0)
{

  isCoreSpectrumLoaded = false;
  isAbInitioSpectraUsed = false;
  ran = new TRandom3();

  m_xsec = new CrossSectionTable();



}

CoreSpectrum::~CoreSpectrum()
{
  ;
}


double CoreSpectrum::antiNuSpectrum(unsigned int coreId, double e_nu)
{
  // Return the anti-neutrino production
  // Units: [neutrinos MeV^-1 s^-1 1E-18]
  if(e_nu>=m_eMax) return 0.0;
  // By default, use linear interpolation between nearest sample points


  // Christine's spectra starts from 1.85 MeV, while we want spectra from 1.80 MeV. So changed it so that it can do linear extrapolation below 1.85 MeV.

  int binIdxLow = 0;
  if (e_nu > m_eMin)
    binIdxLow = (int)((e_nu-m_eMin)/m_binWidth);

  double value = m_dNdE[coreId][binIdxLow];
  static bool linearInterp = true;
  if(linearInterp){
    int binIdxHigh = binIdxLow+1;
    double binLowE = m_eMin + binIdxLow*m_binWidth;
    double dE = e_nu - binLowE;
    double slope = (m_dNdE[coreId][binIdxHigh]
                    - m_dNdE[coreId][binIdxLow])/m_binWidth;
    value += dE*slope;
  }
  return value * 1.0e18;

}

double CoreSpectrum::IBDSpectrum(unsigned int coreId, double e_nu)
{
  // Return the anti-neutrino production
  // Units: [neutrinos MeV^-1 fission^-1]
  if(e_nu>=m_eMax) return 0.0;
  // By default, use linear interpolation between nearest sample points


  // Christine's spectra starts from 1.85 MeV, while we want spectra from 1.80 MeV. So changed it so that it can do linear extrapolation below 1.85 MeV.

  int binIdxLow = 0;
  if (e_nu > m_eMin)
    binIdxLow = (int)((e_nu-m_eMin)/m_binWidth);
  int binIdxHigh = binIdxLow+1;


  double value_low = m_dNdE[coreId][binIdxLow] * m_xsec->inverseBetaDecay(m_eMin + binIdxLow*m_binWidth);
  double value_high = m_dNdE[coreId][binIdxHigh] * m_xsec->inverseBetaDecay(m_eMin + binIdxHigh*m_binWidth);
  double value = value_low;
  //  cout << "=============================" << coreId << " " << e_nu << " " << binIdxLow  << " " << value_low << " " << m_dNdE[coreId][binIdxLow] << endl;

  static bool linearInterp = true;
  if(linearInterp){
    double binLowE = m_eMin + binIdxLow*m_binWidth;
    double dE = e_nu - binLowE;
    double slope = (value_high - value_low) /m_binWidth;
    value += dE*slope;
  }
  if (value < 0) value = 0;

  return value * 1.0e18;

}

double CoreSpectrum::IBDSpectrumBCW(unsigned int coreId, double e_nu){

  double value = 0;

  if(e_nu>=m_eMax) value = 0.0;
  // By default, use linear interpolation between nearest sample points


  int binIdx = 0;
  if (e_nu < m_bcw_bins[0] || e_nu > m_bcw_bins[m_nBcwBins])
    value = 0;
  else{
    for (int ibin = 0; ibin < m_nBcwBins; ibin++){
      if (e_nu > m_bcw_bins[ibin] && e_nu <=m_bcw_bins[ibin+1]){
        binIdx = ibin;
        value = m_dIBDdE[coreId][binIdx]/m_dIBDdE_nom[coreId][binIdx]*IBDSpectrum(coreId,e_nu);
        break;
      }
    }

  }

  return value;

}




void CoreSpectrum::setRandomSeed(unsigned int seed)
{
  ran->SetSeed(seed);
}


void CoreSpectrum::setRandomAntiNuSpectra(){

  if (isAbInitioSpectraUsed){
    return;
  }

  double ranvec[Ncores*MAX_CORESPECTRA_SAMPLES];

  for (int i = 0; i < Ncores*m_nSamples; i++){
    ranvec[i] = ran->Gaus(0,1);
  }

  for (int iCore = 0; iCore < Ncores; iCore++){
    for (int iSample = 0; iSample < m_nSamples; iSample++){
      m_dNdE[iCore][iSample] = m_dNdE_nom[iCore][iSample];
      for (int jCore = 0; jCore < Ncores; jCore++){
        for (int jSample = 0; jSample < m_nSamples; jSample++){
          m_dNdE[iCore][iSample]
            += L[iCore*m_nSamples+iSample][jCore*m_nSamples+jSample]*ranvec[jCore*m_nSamples+jSample];
        }
      }
      //      cout << m_dNdE_nom[iCore][iSample] << "\t" << m_dNdE[iCore][iSample] << endl;

    }
  }

}

void CoreSpectrum::setRandomIBDSpectraBCW(){
  // fluctuate using reactor flux covariance matrix by BCW.
  // Note that BCW uses different binnings.
  // Therefore, first calculate variation in BCW bins, and then make linear interpolation

  double ranvec[Ncores*m_nBcwBins];
  for (int i = 0; i < Ncores* m_nBcwBins; i++){
    ranvec[i] = ran->Gaus(0,1);
  }

  for (int iCore = 0; iCore < Ncores; iCore++){
    for (int iSample = 0; iSample < m_nBcwBins; iSample++){
      m_dIBDdE[iCore][iSample] = m_dIBDdE_nom[iCore][iSample];
      for (int jCore = 0; jCore < Ncores; jCore++){
        for (int jSample = 0; jSample < m_nBcwBins; jSample++){
          m_dIBDdE[iCore][iSample]
            += L_BCW[iCore*m_nBcwBins+iSample][jCore*m_nBcwBins+jSample]*ranvec[jCore*m_nBcwBins+jSample];
        }
      }
      //      cout << m_dIBDdE_nom[iCore][iSample] << "\t" << m_dIBDdE[iCore][iSample] << endl;

    }
  }

}


Bool_t CoreSpectrum::loadSpectra(const char* filename_nom,const char* filename_mcov)
{
  // Load isotope spectra data file
  ifstream fileData(filename_nom);
  if(!fileData.is_open() || !fileData.good()){
    std::cout << "CoreSpectrum::loadSpectra: "
              << "Error: Failed to open data file " << filename_nom << std::endl;
    return false;
  }
  std::cout << "CoreSpectrum::loadSpectra: "
            << "Reading nominal neutrino spectra from " << filename_nom << std::endl;

  std::string line;
  double eMin=0;
  double eMax=0;
  double eNu=0;
  double dNdE[Ncores]={0};
  unsigned int curSample = 0;
  double binWidth = 0;
  while(MAX_CORESPECTRA_SAMPLES > curSample){
    if(fileData.peek()=='#'){
      // Skip lines starting with '#'
      getline(fileData,line);
      continue;
    }else{
      fileData >> eNu;
      for (int i = 0; i < Ncores; i++){
        fileData >> dNdE[i];
      }
    }
    if(!fileData.good()) break;
    if(curSample==0){
      eMin=eNu;
      eMax=eNu;
    }
    if(eNu>eMax) eMax=eNu;
    if(curSample==1) binWidth = eMax-eMin;
    for (int i = 0; i < Ncores; i++){
      m_dNdE_nom[i][curSample] = dNdE[i];
      m_dNdE[i][curSample] = dNdE[i];
    }
    curSample++;
  }
  // Check uniform binning
  double epsilon = 0.0000001;
  double delta = binWidth - (eMax-eMin)/curSample;
  if(delta*delta > epsilon){
    // Doesn't appear to be uniform binning
    std::cout << "CoreSpectrum: Error: binning not uniform in " << filename_nom
              << std::endl;
    return false;
  }
  m_eMin = eMin;
  m_eMax = eMax;
  m_nSamples = curSample;
  m_binWidth = binWidth;

  // Read covarianvematrix

  ifstream fileData_mcov(filename_mcov);
  if(!fileData_mcov.is_open() || !fileData_mcov.good()){
    std::cout << "CoreSpectrum::loadSpectra: "
              << "Error: Failed to open data file " << filename_mcov << std::endl;
    return false;
  }
  std::cout << "CoreSpectrum::loadSpectra: "
            << "Reading neutrino spectra covariance matrix from " << filename_mcov << std::endl;

  for (int i = 0; i < m_nSamples*Ncores; i++){
    for (int j = 0; j < m_nSamples*Ncores; j++){
      fileData_mcov >> m_dNdE_mcov[i][j];
      //      cout << "\t" << m_dNdE_mcov[i][j] << endl;

    }
    //    cout << endl;
  }

  // Increase diagonal part by a tiny bit (0.1%) to make decomposition stable
  // for (int i = 0; i < Ncores; i++){
  //   for (int j = 0; j < m_nSamples; j++){
  //     m_dNdE_mcov[i*m_nSamples+j][i*m_nSamples+j] *=1.000001;
  //     //      m_dNdE_mcov[i*m_nSamples+j][i*m_nSamples+j] *=1.001;
  //     //      m_dNdE_mcov[i*m_nSamples+j][i*m_nSamples+j] += m_dNdE_nom[i][curSample]*0.000001;
  //   }
  // }



  TMatrixD covmatrix(Ncores*MAX_CORESPECTRA_SAMPLES,Ncores*MAX_CORESPECTRA_SAMPLES,&m_dNdE_mcov[0][0]);
  covmatrix.ResizeTo(m_nSamples*Ncores,m_nSamples*Ncores);

  //  covmatrix.Print();
  // //  ematrix_fakedata.Print();
  // Double_t tmpvec[NBins];
  // TMatrixD ranvec(NBins,1);
  // TMatrixD parvec(NBins,1);

  TDecompChol chol(covmatrix);
  Bool_t decomp_success = chol.Decompose();

  TMatrixD cmat(chol.GetU());
  //  cmat.Print();
  TMatrixD tcmat(cmat.Transpose(cmat));
  //  std::cout << "Cholesky matrix---------" << std::endl;
  //  tcmat.Print();

  double * tmp_matrix = tcmat.GetMatrixArray();

  for (int i = 0; i < m_nSamples*Ncores; i++){
    for (int j = 0; j <  m_nSamples*Ncores; j++){
      L[i][j] = tmp_matrix[i*m_nSamples*Ncores + j];
      //      cout << "\t" << L[i][j] << endl;
    }
    //    cout << endl;
  }


  isCoreSpectrumLoaded = true;

  return decomp_success;
}
//////////
Bool_t CoreSpectrum::loadAbInitioSpectra(const char* filename_nom)
{
  // Load isotope spectra data file
  ifstream fileData(filename_nom);
  if(!fileData.is_open() || !fileData.good()){
    std::cout << "CoreSpectrum::loadSpectra: "
              << "Error: Failed to open data file " << filename_nom << std::endl;
    return false;
  }
  std::cout << "CoreSpectrum::loadSpectra: "
            << "Reading AbInitio nominal neutrino spectra from " << filename_nom << std::endl;

  std::string line;
  double eMin=0;
  double eMax=0;
  double eNu=0;
  double dNdE[Ncores]={0};
  unsigned int curSample = 0;
  double binWidth = 0;
  while(MAX_ABINITIO_CORESPECTRA_SAMPLES > curSample){
    if(fileData.peek()=='#'){
      // Skip lines starting with '#'
      getline(fileData,line);
      continue;
    }else{
      fileData >> eNu;
      for (int i = 0; i < Ncores; i++){
        fileData >> dNdE[i];
      }
    }
    if(!fileData.good()) break;
    if(curSample==0){
      eMin=eNu;
      eMax=eNu;
    }
    if(eNu>eMax) eMax=eNu;
    if(curSample==1) binWidth = eMax-eMin;
    //    cout << eNu ;
    for (int i = 0; i < Ncores; i++){
      m_dNdE_nom[i][curSample] = dNdE[i];
      m_dNdE[i][curSample] = dNdE[i];
      //      cout << "\t" << dNdE[i];
    }
    //    cout << endl;
    curSample++;
  }
  // Check uniform binning
  double epsilon = 0.0000001;
  double delta = binWidth - (eMax-eMin)/curSample;
  if(delta*delta > epsilon){
    // Doesn't appear to be uniform binning
    std::cout << "CoreSpectrum: Error: binning not uniform in " << filename_nom
              << std::endl;
    return false;
  }
  m_eMin = eMin;
  m_eMax = eMax;
  m_nSamples = curSample;
  m_binWidth = binWidth;

  cout << "Finished reading ab initio spectra: eMin = " << m_eMin << " eMax = " << m_eMax << " nSamples = " << m_nSamples << " binWidth = "<< m_binWidth << endl;

  isAbInitioSpectraUsed = true;
  return true;
}
//////////
Bool_t CoreSpectrum::setFlatSpectra() // set flat reactor spectra for debugging
{
  for (int iSample = 0; iSample <  m_nSamples; iSample++){
    for (int i = 0; i < Ncores; i++){
      m_dNdE_nom[i][iSample] = 10.0;
      m_dNdE[i][iSample] = 10.0;
    }
  }
  return true;
}

////////////////////////////////////////
Bool_t CoreSpectrum::loadCovMatrixBCW(const char* filename_mcov)
{

  m_bcw_bins[0] = 1.807;
  for (Int_t i = 1; i < m_nBcwBins; i++){
    m_bcw_bins[i] = 2.125 + ((double)i-1)*0.25;
  }
  m_bcw_bins[m_nBcwBins] = 12.775;

  for (Int_t i = 0; i < m_nBcwBins; i++){
    m_bcw_bincenter[i] = 0.5*(m_bcw_bins[i] + m_bcw_bins[i+1]);
  }

  for (int iCore = 0; iCore < Ncores; iCore++){
    for (int ibin = 0; ibin < m_nBcwBins; ibin++){
      m_dIBDdE_nom[iCore][ibin] = 0;
      for (int iSample = 0; iSample < m_nSamples; iSample++){
        Double_t e = ((double)iSample + 0.5)*m_binWidth + m_eMin;
        if(e >= m_bcw_bins[ibin] && e < m_bcw_bins[ibin+1]){
          m_dIBDdE_nom[iCore][ibin] += IBDSpectrum(iCore,e)*1.0e-18*m_binWidth/(m_bcw_bins[ibin+1]-m_bcw_bins[ibin]);
        }
      }
      m_dIBDdE[iCore][ibin] = m_dIBDdE_nom[iCore][ibin];
    }
  }

  // Read covariancematrix

  TFile *  fileData_mcov = new TFile(filename_mcov);
  if(!fileData_mcov->IsOpen()){
    std::cout << "CoreSpectrum::loadSpectra: "
              << "Error: Failed to open data file " << filename_mcov << std::endl;
    return false;
  }
  std::cout << "CoreSpectrum::loadSpectra: "
            << "Reading neutrino spectra covariance matrix from " << filename_mcov << std::endl;


  TMatrixDSym * covmatrix  = (TMatrixDSym*)fileData_mcov->Get("covarEbins");
  if (covmatrix->GetNcols() != m_nBcwBins*Ncores){
    std::cout << "Inconsistent number bins in BCW reactor flux covariance matrix!" << std::endl;
    return false;
  }

  double covmatrix_abs[Ncores*m_nBcwBins][Ncores*m_nBcwBins];
  double * tmp_covmatrix_rel = covmatrix->GetMatrixArray();

  for (int iCore = 0; iCore < Ncores; iCore++){
    for (int ibin = 0; ibin < m_nBcwBins; ibin++){
      int iii = iCore*m_nBcwBins+ibin;
      for (int jCore = 0; jCore < Ncores; jCore++){
        for (int jbin = 0; jbin < m_nBcwBins; jbin++){
          int jjj = jCore*m_nBcwBins+jbin;
          covmatrix_abs[iii][jjj]=tmp_covmatrix_rel[iii*Ncores*m_nBcwBins+jjj]*m_dIBDdE_nom[iCore][ibin]*m_dIBDdE_nom[jCore][jbin];
        }
      }
    }
  }

  TMatrixD cmat_abs(Ncores*m_nBcwBins,Ncores*m_nBcwBins);
  cmat_abs.SetMatrixArray(&covmatrix_abs[0][0]);

  //  cmat_abs.Print();

  TDecompChol chol(cmat_abs);
  Bool_t decomp_success = chol.Decompose();

  TMatrixD cmat(chol.GetU());
  //  cmat.Print();
  TMatrixD tcmat(cmat.Transpose(cmat));
  //  std::cout << "Cholesky matrix---------" << std::endl;
  //  tcmat.Print();

  double * tmp_matrix = tcmat.GetMatrixArray();

  for (int i = 0; i < m_nBcwBins*Ncores; i++){
    for (int j = 0; j <  m_nBcwBins*Ncores; j++){
      L_BCW[i][j] = tmp_matrix[i*m_nBcwBins*Ncores + j];
      //      cout << "\t" << L[i][j] << endl;
    }
    //    cout << endl;
  }

  isBcwCovMatrixLoaded = true;


  return decomp_success;
}

double CoreSpectrum::eMin()
{
  return m_eMin;
}

double CoreSpectrum::eMax()
{
  return m_eMax;
}


bool CoreSpectrum::loadCrossSectionTable(const char* filename_xsec){
  m_xsec->load(filename_xsec);
  return true;
}
