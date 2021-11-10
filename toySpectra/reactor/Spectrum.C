#include "Spectrum.h"

#include "Config.h"
#include "CoreSpectrum.h"
#include "CrossSectionTable.h"
#include "DataSet.h"
#include "IsotopeTable.h"
#include "Paths.h"
#include "Utils.h"

#include <TCanvas.h>
#include <TDecompChol.h>
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include <TMatrixD.h>
#include "TRandom3.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <cassert>

using namespace std;
using namespace Config;

Spectrum::Spectrum() :
    m_eMin(0), m_eMax(0), m_nSamples(0), m_binWidth(0), m_nSamplesBkg(0),
    m_binWidthBkg(0), m_alpha(0), m_beta(0), m_detectorResolution(0)
{
  m_osccalc = new OscCalc();


  m_randomizeIsotopeFraction = -1;
  m_randomizeReactorPower = -1;
  m_randomizeIavDistortion = -1;


  // latest value from doc-9401-v7
  m_detectorEfficiency_Dt = 0.9870;
  m_detectorEfficiency_Ep = 0.9981;
  // m_detectorEfficiency_Ed_nominal = 0.927;
  m_detectorEfficiency_flash = 0.9998;
  m_detectorEfficiency_nGd = 0.842;
  m_detectorEfficiency_spill = 1.049;

  // The (AD-dependent) nominal delayed cut efficiency is now taken from the
  // Theta13-inputs file. See extractPredictorData.

  // for (Int_t idet = 0; idet < Ndetectors; idet++) {
  //   m_detectorEfficiency_Ed[idet] = m_detectorEfficiency_Ed_nominal;
  // }

  ran = new TRandom3();

  // FIXME: hard-coded detector efficiency uncertainty
  m_detectorEfficiency_rel_error = 0.0011; // DocDB-10956


  UseChristineModel = true;

  // Caution: These reactor power and their errors are NOT used if we chose to
  // use Chistine's covariance matrix

  // FIXME: hard-coded reactor power.
  // It's an average of first 140 days of data
  for (Int_t istage = 0; istage < Nstage; istage++) {
    if (istage == 0) {
      m_nominalReactorPower[istage][0] = 2.178;
      m_nominalReactorPower[istage][1] = 2.876;
      m_nominalReactorPower[istage][2] = 2.389;
      m_nominalReactorPower[istage][3] = 2.444;
      m_nominalReactorPower[istage][4] = 2.820;
      m_nominalReactorPower[istage][5] = 2.787;
    } else {
      for (int i = 0; i < 6; i++) {
        m_nominalReactorPower[istage][i] = 2.9;
      }
    }
  }

  for (Int_t istage = 0; istage < Nstage; istage++) {
    for (int i = 0; i < Ncores; i++) {
      m_reactorPower[istage][i] = m_nominalReactorPower[istage][i];
    }
  }

  // FIXME: hard-coded reactor power uncertainty (relative).
  for (Int_t istage = 0; istage < Nstage; istage++) {
    m_reactorPowerError[istage][0] = 0.005;
    m_reactorPowerError[istage][1] = 0.005;
    m_reactorPowerError[istage][2] = 0.005;
    m_reactorPowerError[istage][3] = 0.005;
    m_reactorPowerError[istage][4] = 0.005;
    m_reactorPowerError[istage][5] = 0.005;
  }

  m_iav_error = 0.04; // nominal value: 4%, based on Bryce's document


  ///////////////////////////////////////////////////////////////////

  reso_func = new TF1("reso_func", reso_func_bcw, 0, 20, 3);
  // reso_func->SetParameters(0.0148,0.0869,0.0271); // based on Doc-8982 //old
  reso_func->SetParameters(0.016, 0.081, 0.026); // based on 1230 days paper
  for (int ii = 0; ii < Ndetectors; ii++) {
    m_detectorResolution_bias[ii] = 0;
  }


  ///////////////////////////////////////////////////////////////////

  //  enu_to_epositron_func = new
  //  TF1("enu_to_epositron_func",this,&Spectrum::enu_to_epositron,0,20,1,"Spectrum","enu_to_epositron");
  enu_to_epositron_func =
      new TF1("enu_to_epositron_func", enu_to_epositron, 0, 20, 1);
  enu_to_epositron_func->SetParameter(0, -2);

  //  enu_to_epositron_func = new
  //  TF1("enu_to_epositron_func",&Spectrum::enu_to_epositron_0th,0,20,0,"Spectrum","enu_to_epositron_0th");

  ///////////////////////////////////////////////////////////////////


  // dy default, all backgrounds are included into the spectra.
  m_removeAccBg = false;
  m_removeLi9Bg = false;
  m_removeFnBg = false;
  m_removeAmcBg = false;
  m_removeAlnBg = false;
}

Spectrum::~Spectrum()
{
  if (m_osccalc) {
    delete m_osccalc;
    m_osccalc = 0;
  }
}

int Spectrum::nSamples()
{
  return m_nSamples;
}

double Spectrum::binWidth()
{
  return m_binWidth;
}

int Spectrum::nSamplesBkg()
{
  return m_nSamplesBkg;
}

double Spectrum::binWidthBkg()
{
  return m_binWidthBkg;
}

double Spectrum::antiNu(int istage, int idet, double energy_nu)
{
  // Return the anti-neutrino spectrum
  // Units: [antineutrinos MeV^-1]
  if (energy_nu < m_eMin || energy_nu >= m_eMax)
    return 0.0;
  // By default, use nearest sample point
  // FIXME: could consider alternate interpolations or splines,
  //        but must measure impact on speed
  int binIdx = (int)((energy_nu - m_eMin) / m_binWidth);
  return m_antiNuSpectrum[istage][idet][binIdx];
}

double Spectrum::positronTrue(int istage, int idet, double energy_e)
{
  // Return the true positron spectrum (including annihilation energy)
  // Units: [positrons MeV^-1]
  if (energy_e < m_eMin || energy_e >= m_eMax)
    return 0.0;
  // By default, use nearest sample point
  // FIXME: could consider alternate interpolations or splines,
  //        but must measure impact on speed
  int binIdx = (int)((energy_e - m_eMin) / m_binWidth);
  return m_positronTrueSpectrum[istage][idet][binIdx];
}

double Spectrum::positronDetected(int istage, int idet, double energy_e)
{
  // Return the detected positron spectrum (including detector effects)
  // Units: [positrons MeV^-1]
  if (energy_e < m_eMin || energy_e >= m_eMax)
    return 0.0;
  // By default, use nearest sample point
  // FIXME: could consider alternate interpolations or splines,
  //        but must measure impact on speed
  int binIdx = (int)((energy_e - m_eMin) / m_binWidth);
  return m_positronDetectedSpectrum[istage][idet][binIdx];
}

double* Spectrum::energyArray(int idet)
{
  // Return an array representing the energy axis for spectra
  return m_energy[idet];
}
double* Spectrum::energyArrayBkg(int idet)
{
  // Return an array representing the energy axis for spectra
  return m_energy_bkg[idet];
}

double* Spectrum::antiNuArray(int istage, int idet)
{
  // Return the raw antineutrino spectrum data
  return m_antiNuSpectrum[istage][idet];
}

double* Spectrum::antiNuNoOscArray(int istage, int idet)
{
  // Return the raw antineutrino spectrum data
  return m_antiNuSpectrumNoOsc[istage][idet];
}

double* Spectrum::positronTrueArray(int istage, int idet)
{
  // Return the true positron spectrum data
  return m_positronTrueSpectrum[istage][idet];
}
double* Spectrum::positronIavDistortedArray(int istage, int idet)
{
  // Return the true positron spectrum data afte Iav distortion
  return m_positronIavDistortedSpectrum[istage][idet];
}

double* Spectrum::positronNLArray(int istage, int idet)
{
  // Return the detected positron spectrum data
  return m_positronNLSpectrum[istage][idet];
}

double* Spectrum::positronDetectedArray(int istage, int idet)
{
  // Return the detected positron spectrum data
  return m_positronDetectedSpectrum[istage][idet];
}

double* Spectrum::positronNLCorrectedArray(int istage, int idet)
{
  // Return the non-linearity corrected detected positron spectrum data
  return m_positronNLCorrectedSpectrum[istage][idet];
}

double* Spectrum::bgDetectedArray(int istage, int idet)
{
  // Return the detected background spectrum data
  return m_bgDetectedSpectrum[istage][idet];
}
double* Spectrum::bgNLCorrectedArray(int istage, int idet)
{
  // Return the non-linearity corrected detected background spectrum data
  return m_bgNLCorrectedSpectrum[istage][idet];
}

// Get the oscillation parameter
double Spectrum::deltaMSqee()
{
  return m_osccalc->GetDeltaM2_ee();
}

// Update the oscillation parameters
/*void Spectrum::setDeltaMSqee(double deltaMSqee)
  {
  if(deltaMSqee == m_oscillation->deltaMSqee()) return;
  m_oscillation->setDeltaMSqee(deltaMSqee);
  this->updateOscillation();
  }*/

void Spectrum::setOscillation(double sinSq2Th12, double sinSq2Th13,
                              double deltaMSq21, double deltaMSqee)
{
  m_sinSq2Th12 = sinSq2Th12;
  m_sinSq2Th13 = sinSq2Th13;
  m_deltaMSq21 = deltaMSq21;
  m_deltaMSqee = deltaMSqee;

  m_osccalc->SetDeltaM2_ee(deltaMSqee);

  m_osccalc->SetDeltaM2_21(deltaMSq21);
  m_osccalc->SetTheta12(TMath::ASin(sqrt(sinSq2Th12)) * 0.5);
}

void Spectrum::setRandomSolarOscPars()
{
  m_sinSq2Th12 = m_sinSq2Th12_nominal + m_sinSq2Th12_error * ran->Gaus(0, 1);
  m_deltaMSq21 = m_deltaMSq21_nominal + m_deltaMSq21_error * ran->Gaus(0, 1);

  m_osccalc->SetDeltaM2_ee(m_deltaMSqee); //<--just in case

  m_osccalc->SetDeltaM2_21(m_deltaMSq21);
  m_osccalc->SetTheta12(TMath::ASin(sqrt(m_sinSq2Th12)) * 0.5);
}

void Spectrum::setRandomDm2ee()
{
  m_deltaMSqee = m_deltaMSqee_nominal + m_deltaMSqee_error * ran->Gaus(0, 1);

  m_osccalc->SetDeltaM2_ee(m_deltaMSqee);
  m_osccalc->SetDeltaM2_21(m_deltaMSq21); //<--just in case
  m_osccalc->SetTheta12(TMath::ASin(sqrt(m_sinSq2Th12)) *
                        0.5); //<--just in case
}

void Spectrum::setEnergyScale(double alpha, double beta)
{
  if (alpha == m_alpha && beta == m_beta)
    return;
  m_alpha = alpha;
  m_beta = beta;

  cout << "WARNING: do not use setEenrgyScale() function anymore... " << endl;
  this->updatePositronDetected();
}

/*void Spectrum::updateOscillation()
  {
  // Updated internal data based on oscillation parameter change
  double conv_km_to_m = 1000.;
  double totalCounts = 0;
  for(int idx=0; idx<m_nSamples; idx++){
  if(m_energy[idx]==0) continue;
  m_survivalProbability[idx] =
  m_oscillation->survivalProbability(m_energy[idx],
  m_detectorDistance[][]*conv_km_to_m);
  m_antiNuSpectrum[idx] = (m_survivalProbability[idx]
  * m_antiNuSpectrumNoOsc[idx]);
  totalCounts += m_antiNuSpectrum[idx] * m_binWidth;
  }
  //cout << "Updated total antinu Osc counts: " << totalCounts << endl;
  this->updatePositronTrue();
  }*/


void Spectrum::updateAntinu(int icore_select)
{
  double avgEnergyPerFission[Nstage][Ncores] = {{0}};


  if (UseChristineModel || m_useAbInitioSpectra > 0) {
    if (m_randomizeCoreSpectra > 0) {
      for (int istage = 0; istage < Nstage; istage++) {
        // YN: Add option to BCW flux uncertainty
        if (m_useBcwFluxUncertainty > 0) {
          m_corespectrum[istage]->setRandomIBDSpectraBCW();
        } else {
          m_corespectrum[istage]->setRandomAntiNuSpectra();
        }
      }
    }
  } /*else{
      if (m_randomizeIsotopeFraction > 0)
      for(int istage=0;istage<Nstage;++istage)
      for(int icore=0;icore<Ncores;++icore)
      m_isotopes[istage][icore]->setRandomFissionFraction();

      for(int istage=0;istage<Nstage;++istage){
      for(int icore=0;icore<Ncores;++icore){
      avgEnergyPerFission[istage][icore] = 0;

      for(int isoIdx=1; isoIdx<5; isoIdx++){
      // avgEnergyPerFission += (m_isotopes->nominalFissionFraction(isoIdx)
      //         * m_isotopes->meanEnergyPerFission(isoIdx));
      avgEnergyPerFission[istage][icore] +=
      (m_isotopes[istage][icore]->FissionFraction(isoIdx)
      * m_isotopes[istage][icore]->meanEnergyPerFission(isoIdx));
      }
      }
      }
      }*/


  if (m_randomizeReactorPower > 0)
    this->setRandomReactorPower();

  if (m_randomizeSolarOscPars > 0)
    this->setRandomSolarOscPars();

  if (m_randomizeDm2ee > 0)
    this->setRandomDm2ee();

  double conv_mev_per_s_gw = 6.2415e21;
  double conv_km2_per_cm2 = 1.0000e-10;
  double conv_s_per_year = 24 * 60 * 60 * 365.;
  double conv_km_to_m = 1000.;
  for (int istage = 0; istage < Nstage; ++istage) {
    for (int idet = 0; idet < Ndetectors; idet++) {
      double totalCounts = 0;
      for (int idx = 0; idx < m_nSamples; idx++) {
        m_energy[idet][idx] = (0.5 + idx) * m_binWidth + m_eMin;
        m_antiNuSpectrum[istage][idet][idx] = 0;
        m_antiNuSpectrumNoOsc[istage][idet][idx] = 0;

        for (int icore = 0; icore < Ncores; ++icore) {
          if (icore_select >= 0 && icore != icore_select)
            continue; // select only one core if specified.

          float oscproba = 1;
          if (m_energy[idx] != 0) {
            // oscproba=m_osccalc->OscProb(m_detectorDistance[idet][icore]*conv_km_to_m,
            // m_energy[idet][idx], m_sinSq2Th13);
            double Elow = m_energy[idet][idx] - m_binWidth * 0.5;
            double Ehigh = m_energy[idet][idx] + m_binWidth * 0.5;
            // super ad-hock introduction of baseline smearing.
            double Lwidth = 0.0; // in km
            double Llow = m_detectorDistance[idet][icore] - Lwidth;
            double Lhigh = m_detectorDistance[idet][icore] + Lwidth;
            // double prob_bef =
            // m_osccalc->OscProbInt(m_detectorDistance[idet][icore]*conv_km_to_m,Elow,m_sinSq2Th13);
            // double prob_aft =
            // m_osccalc->OscProbInt(m_detectorDistance[idet][icore]*conv_km_to_m,Ehigh,m_sinSq2Th13);
            double prob_bef =
                m_osccalc->OscProbInt(Lhigh * conv_km_to_m, Elow, m_sinSq2Th13);
            double prob_aft =
                m_osccalc->OscProbInt(Llow * conv_km_to_m, Ehigh, m_sinSq2Th13);
            //    oscproba = (prob_bef -
            //    prob_aft)/(m_detectorDistance[idet][icore]*conv_km_to_m/Elow -
            //    m_detectorDistance[idet][icore]*conv_km_to_m/Ehigh);
            oscproba = (prob_bef - prob_aft) / (Lhigh * conv_km_to_m / Elow -
                                                Llow * conv_km_to_m / Ehigh);
            // this integration assumes flat distribution in L/E space
          }

          if (UseChristineModel || m_useAbInitioSpectra > 0) {
            // YN: Adding BCW flux uncertainty
            double flux = 0;
            if (m_useBcwFluxUncertainty > 0) {
              flux = m_corespectrum[istage]->IBDSpectrumBCW(
                  icore, m_energy[idet][idx]);
            } else {
              flux = m_corespectrum[istage]->IBDSpectrum(icore,
                                                         m_energy[idet][idx]);
            }


            m_antiNuSpectrum[istage][idet][idx] +=
                (m_reactorPower[istage][icore] /
                 m_nominalReactorPower[istage][icore] *
                 m_runningTime[istage][idet] * m_detectorSize[idet] *
                 m_targetProtonsPerKton
                 // * m_corespectrum->antiNuSpectrum(icore,m_energy[idet][idx])
                 // * m_xsec->inverseBetaDecay(m_energy[idet][idx])
                 * flux * m_detectorEfficiency[istage][idet] *
                 conv_km2_per_cm2 * conv_s_per_year /
                 (4 * TMath::Pi() * m_detectorDistance[idet][icore] *
                  m_detectorDistance[idet][icore]) *
                 oscproba);


            // for producing the superhistogram, we don't need to apply detector
            // efficiencies, etc
            m_antiNuSpectrumNoOsc[istage][idet][idx] +=
                (m_reactorPower[istage][icore] /
                 m_nominalReactorPower[istage][icore] *
                 m_runningTime[istage][1] // changed to 1 since AD1 was taken
                                          // from analysis
                 * m_detectorSize[0] *
                 m_targetProtonsPerKton
                 // * m_corespectrum->antiNuSpectrum(icore,m_energy[idet][idx])
                 // * m_xsec->inverseBetaDecay(m_energy[idet][idx])
                 * flux * conv_s_per_year);


          } /*else{

 for(int isoIdx=1; isoIdx<5; isoIdx++){

 m_antiNuSpectrum[idet][idx] += (m_runningTime[idet]
 * m_detectorSize[idet]
 * m_targetProtonsPerKton
 * m_reactorPower[icore]
 * m_isotopes[icore]->FissionFraction(isoIdx)
 * m_isotopes[icore]->antiNuSpectrum(isoIdx, m_energy[idet][idx])
 * m_xsec->inverseBetaDecay(m_energy[idet][idx])
 * m_detectorEfficiency[idet]
 * conv_mev_per_s_gw
 * conv_km2_per_cm2
 * conv_s_per_year
 / (4 * TMath::Pi()
 * m_detectorDistance[idet][icore]
 *m_detectorDistance[idet][icore])
 / avgEnergyPerFission[icore]
 * oscproba
 );


 // for producing the superhistogram, we don't need to apply detector
 efficiencies, etc m_antiNuSpectrumNoOsc[idet][idx] += (m_runningTime[0]
 * m_detectorSize[0]
 * m_targetProtonsPerKton
 * m_reactorPower[icore]
 * m_isotopes[icore]->FissionFraction(isoIdx)
 * m_isotopes[icore]->antiNuSpectrum(isoIdx,
 m_energy[idet][idx])
 * m_xsec->inverseBetaDecay(m_energy[idet][idx])
 * conv_mev_per_s_gw
 * conv_s_per_year
 / avgEnergyPerFission[icore]
 );
 }//isotope loop
 }*/
        }   // core loop
        totalCounts += m_antiNuSpectrum[istage][idet][idx] * m_binWidth;
      } // samples loop
      // cout << "Total counts [NoOsc]: " << totalCounts << endl;
    } // detector loop
  }   // Stage loop
  // this->updateOscillation();

  this->updatePositronTrue();
}


//  void Spectrum::updatePositronTrue(int eNuIdx_select)
void Spectrum::updatePositronTrue(double eNu_min, double eNu_max)
{
  // Updated true positron spectrum
  //   Convert neutrino energy to positron energy


  if (m_randomizeIavDistortion > 0) {
    this->setRandomIavDistortion();
  }


  ///////////////
  for (int istage = 0; istage < Nstage; istage++) {
    for (int idet = 0; idet < Ndetectors; idet++) {
      double totalCounts = 0;

      // // OLD code that simply shift energy by 0.8 MeV
      // int offset = (int)(0.8 / m_binWidth);
      // for(int idx=0; idx<m_nSamples; idx++){
      //   int eNuIdx = (int)(idx + offset);
      //   if(eNuIdx<0 || eNuIdx>=m_nSamples){
      //     m_positronTrueSpectrum[idet][idx] = 0;
      //     continue;
      //   }
      //   if (eNuIdx_select >= 0 && eNuIdx != eNuIdx_select){ // select single
      //   true Enu bin if specified.
      //     m_positronTrueSpectrum[idet][idx] = 0;
      //     continue;
      //   }

      //   m_positronTrueSpectrum[idet][idx] = m_antiNuSpectrum[idet][eNuIdx];
      //   totalCounts += m_positronTrueSpectrum[idet][idx] * m_binWidth;
      // }


      for (int idx = 0; idx < m_nSamples; idx++) {
        double e_positron = m_energy[idet][idx];
        if (e_positron < 1.022) {
          m_positronTrueSpectrum[istage][idet][idx] = 0;
          continue;
        }
        double e_nu = enu_to_epositron_func->GetX(e_positron);
        int e_nu_Idx = (int)((e_nu - m_eMin) / m_binWidth);
        if (e_nu_Idx < 0 || e_nu_Idx >= m_nSamples) {
          m_positronTrueSpectrum[istage][idet][idx] = 0;
          continue;
        }
        if (eNu_min >= 0 &&
            (e_nu < eNu_min ||
             e_nu > eNu_max)) { // select single true Enu bin if specified.
          m_positronTrueSpectrum[istage][idet][idx] = 0;
          continue;
        }
        double binScaling = enu_to_epositron_func->Derivative(e_nu, 0, 0.0001);
        double dNdE = 0;
        int sign = 1;
        double dE = e_nu - (0.5 + e_nu_Idx) * m_binWidth;
        if (dE < 0)
          sign = -1;

        if (e_nu_Idx == m_nSamples - 1) {
          dNdE = m_antiNuSpectrum[istage][idet][e_nu_Idx];
        } else {
          dNdE = (TMath::Abs(dE) / m_binWidth) *
                     m_antiNuSpectrum[istage][idet][e_nu_Idx + sign] +
                 (1 - TMath::Abs(dE) / m_binWidth) *
                     m_antiNuSpectrum[istage][idet][e_nu_Idx];
        }

        m_positronTrueSpectrum[istage][idet][idx] = dNdE / binScaling;
      }

      // ///////////////

      // cout << "Updated total positron true counts: " << totalCounts << endl;

      // Apply IAV correction
      for (int idx = 0; idx < m_nSamples; idx++) {
        m_positronIavDistortedSpectrum[istage][idet][idx] = 0;
      }
      for (int idx = 0; idx < m_nSamples; idx++) {
        for (int jdx = 0; jdx < idx + 1; jdx++) {
          m_positronIavDistortedSpectrum[istage][idet][jdx] +=
              m_iav_corr[idet][idx][jdx] *
              m_positronTrueSpectrum[istage][idet][idx];
        }
        //      cout << m_positronIAVDistortedSpectrum[idet][idx] << endl;
        // m_positronIavDistortedSpectrum[idet][idx]
        //   = m_positronTrueSpectrum[idet][idx];
      }
    }
  }

  this->updatePositronDetected();
}

void Spectrum::updatePositronDetected()
{
  // Updated detected positron spectrum
  //   Convert positron true energy to detected energy using
  //   simple detector response model

  // randomize various parameters if specified
  if (m_randomizeScintiNonLinear > 0) {
    if (m_useIhepNonLinearModel > 0)
      this->setRandomIhepScintiNonLinear();
    else if (m_useBcwNonLinearModel > 0)
      this->setRandomBcwScintiNonLinear();
    else if (m_useLbnlNonLinearModel > 0)
      this->setRandomLbnlScintiNonLinear();
    else if (m_useUnifiedNonLinearModel > 0 ||
             m_use2015NonLinearModel > 0) // original unified model and 2015
                                          // updates has the same structure
      this->setRandomUnifiedScintiNonLinear();
    else
      this->setRandomScintiNonLinear();
  }
  if (m_randomizeElecNonLinear > 0)
    this->setRandomElecNonLinear();

  if (m_randomizeResolution > 0)
    this->setRandomResolution();

  if (m_randomizeRelativeEnergyScale > 0)
    this->setRandomRelativeEnergyScale();

  if (m_randomizeAbsoluteEnergyScale > 0)
    this->setRandomAbsoluteEnergyScale();


  if (m_randomizeRelativeEnergyOffset > 0)
    this->setRandomRelativeEnergyOffset();

  if (m_randomizeAbsoluteEnergyOffset > 0)
    this->setRandomAbsoluteEnergyOffset();

  if (m_randomizeDetectorEfficiency > 0)
    this->setRandomDetectorEfficiency(); // this must be called after relative
                                         // energy scale uncertainty


  for (int istage = 0; istage < Nstage; istage++) {
    for (int idet = 0; idet < Ndetectors; idet++) {
      // First: Apply detector non-linear energy scale
      double totalCounts = 0;
      if (m_useIhepNonLinearModel || m_useBcwNonLinearModel ||
          m_useLbnlNonLinearModel || m_useUnifiedNonLinearModel ||
          m_use2015NonLinearModel) {
        if (m_useIhepNonLinearModel) {
          // reset non-linear parameters
          nl_func->SetParameters(m_ihep_nl_par[0], m_ihep_nl_par[1],
                                 m_ihep_nl_par[2], m_ihep_nl_par[3],
                                 m_ihep_nl_par[4],
                                 m_abs_escale * m_rel_escale[idet],
                                 m_abs_eoffset + m_rel_eoffset[istage][idet]);
        } else if (m_useBcwNonLinearModel) {
          nl_func->SetParameters(m_bcw_elec_nl_par[0], m_bcw_elec_nl_par[1],
                                 m_bcw_elec_nl_par[2], m_bcw_elec_nl_par[3],
                                 m_bcw_elec_nl_par[4],
                                 m_abs_escale * m_rel_escale[idet],
                                 m_abs_eoffset + m_rel_eoffset[istage][idet]);


        } else if (m_useLbnlNonLinearModel) {
          nl_func->SetParameters(m_lbnl_nl_par[0], m_lbnl_nl_par[1],
                                 m_lbnl_nl_par[2],
                                 m_abs_escale * m_rel_escale[idet],
                                 m_abs_eoffset + m_rel_eoffset[istage][idet]);

        } else if (m_useUnifiedNonLinearModel > 0 ||
                   m_use2015NonLinearModel > 0) {
          for (Int_t i = 0; i < m_num_unified_nl_pars; i++) {
            nl_func->SetParameter(i, m_unified_nl_par[i]);
          }

          nl_func->SetParameter(m_num_unified_nl_pars,
                                m_abs_escale * m_rel_escale[idet]);
          nl_func->SetParameter(m_num_unified_nl_pars + 1,
                                m_abs_eoffset + m_rel_eoffset[istage][idet]);
        }


        for (int idx = 0; idx < m_nSamples; idx++) {
          double energy_nl = m_energy[idet][idx];
          double energy_true = nl_func->GetX(energy_nl);
          int eTrueIdx = (int)((energy_true - m_eMin) / m_binWidth);
          if (eTrueIdx < 0 || eTrueIdx >= m_nSamples) {
            m_positronNLSpectrum[istage][idet][idx] = 0;
            continue;
          }

          double binScaling = nl_func->Derivative(energy_true);
          double dNdE = 0;
          int sign = 1;
          double dEtrue = energy_true - (0.5 + eTrueIdx) * m_binWidth;
          if (dEtrue < 0)
            sign = -1;

          if (eTrueIdx == m_nSamples - 1) {
            dNdE = m_positronIavDistortedSpectrum[istage][idet][eTrueIdx];
            // }if(TMath::Abs(m_positronIavDistortedSpectrum[idet][eTrueIdx+sign]
            // -m_positronIavDistortedSpectrum[idet][eTrueIdx])/m_positronIavDistortedSpectrum[idet][eTrueIdx]
            // > 0.5){
            //   dNdE = m_positronIavDistortedSpectrum[idet][eTrueIdx];
          } else {
            // dNdE = (m_positronIavDistortedSpectrum[idet][eTrueIdx]
            //         +((dEtrue/m_binWidth)
            //           *(m_positronIavDistortedSpectrum[idet][eTrueIdx+1]
            //             -m_positronIavDistortedSpectrum[idet][eTrueIdx])));

            // change the interpolation method better
            dNdE = (TMath::Abs(dEtrue) / m_binWidth) *
                       m_positronIavDistortedSpectrum[istage][idet]
                                                     [eTrueIdx + sign] +
                   (1 - TMath::Abs(dEtrue) / m_binWidth) *
                       m_positronIavDistortedSpectrum[istage][idet][eTrueIdx];
          }
          //        cout << energy_nl << "\t" << energy_true << "\t" <<
          //        dNdE/binScaling << endl;

          m_positronNLSpectrum[istage][idet][idx] = dNdE / binScaling;
          //        m_positronNLSpectrum[idet][idx] =
          //        m_positronIavDistortedSpectrum[idet][idx];
          totalCounts += m_positronNLSpectrum[istage][idet][idx] * m_binWidth;

          // cout << energy_true << "\t" << energy_nl << "\t"
          //      << dNdE << "\t" << binScaling << endl;
        }


      } else {
        // Original formula:
        // E_NL = (1 - alpha/(E_true*E_true)) * beta * E_true
        // Alpha: non-liner effect at low energy
        // Beta: overall energy scale

        double alpha = m_alpha;
        double beta = m_abs_escale * m_rel_escale[idet];

        if (beta == 1.0 && alpha == 0.0) {
          for (int idx = 0; idx < m_nSamples; idx++) {
            m_positronNLSpectrum[istage][idet][idx] =
                m_positronIavDistortedSpectrum[istage][idet][idx];
          }
        } else {
          for (int idx = 0; idx < m_nSamples; idx++) {
            double energy_nl = m_energy[idet][idx];
            double energy_true =
                0.5 * ((energy_nl / beta) +
                       TMath::Sqrt(4 * alpha +
                                   ((energy_nl * energy_nl) / (beta * beta))));
            int eTrueIdx = (int)((energy_true - m_eMin) / m_binWidth);
            if (eTrueIdx < 0 || eTrueIdx >= m_nSamples) {
              m_positronNLSpectrum[istage][idet][idx] = 0;
              continue;
            }

            double binScaling =
                (1 + alpha / (energy_true * energy_true)) * beta;
            double dNdE = 0;
            double dEtrue = energy_true - eTrueIdx * m_binWidth;
            if (eTrueIdx == m_nSamples - 1) {
              dNdE = m_positronIavDistortedSpectrum[istage][idet][eTrueIdx];
            } else {
              // Linear interpolate between points
              dNdE =
                  (m_positronIavDistortedSpectrum[istage][idet][eTrueIdx] +
                   ((dEtrue / m_binWidth) *
                    (m_positronIavDistortedSpectrum[istage][idet]
                                                   [eTrueIdx + 1] -
                     m_positronIavDistortedSpectrum[istage][idet][eTrueIdx])));
            }
            m_positronNLSpectrum[istage][idet][idx] = dNdE / binScaling;
            totalCounts += m_positronNLSpectrum[istage][idet][idx] * m_binWidth;
            // m_positronDetectedSpectrum[idx] = m_positronNLSpectrum[idx];
          }
          // cout << "Updated total positron NL counts: " << totalCounts <<
          // endl;
        }
      }


      // Second: Apply detector resolution
      totalCounts = 0;
      double resolutionRange = 8; // sigma
      for (int idx = 0; idx < m_nSamples; idx++) {
        m_positronDetectedSpectrum[istage][idet][idx] = 0;
      }

      for (int idx = 0; idx < m_nSamples; idx++) {
        if (m_detectorResolution == 0) {
          m_positronDetectedSpectrum[istage][idet][idx] =
              m_positronNLSpectrum[istage][idet][idx];
          continue;
        }
        if (m_positronNLSpectrum[istage][idet][idx] == 0)
          continue;
        double energy_nl = m_energy[idet][idx];
        //      double sigma = m_detectorResolution * TMath::Sqrt(energy_nl);
        double sigma =
            (reso_func->Eval(energy_nl) + m_detectorResolution_bias[idet]) *
            energy_nl;
        double minDetE = energy_nl - resolutionRange * sigma;
        double maxDetE = energy_nl + resolutionRange * sigma;
        int minDetEIdx = (int)((minDetE - m_eMin) / m_binWidth);
        int maxDetEIdx = (int)((maxDetE - m_eMin) / m_binWidth);
        if (minDetEIdx < 0)
          minDetEIdx = 0;
        if (maxDetEIdx >= m_nSamples)
          maxDetEIdx = m_nSamples - 1;
        for (int detIdx = minDetEIdx; detIdx <= maxDetEIdx; detIdx++) {
          if (detIdx == 0)
            continue;
          double gausFactor =
              TMath::Gaus((energy_nl - m_energy[idet][detIdx]), 0, sigma, true);
          m_positronDetectedSpectrum[istage][idet][detIdx] +=
              (m_positronNLSpectrum[istage][idet][idx] * gausFactor *
               m_binWidth);
        }
      }

      for (int idx = 0; idx < m_nSamples; idx++) {
        totalCounts +=
            m_positronDetectedSpectrum[istage][idet][idx] * m_binWidth;
      }
      // cout << "Updated total positron Det counts: " << totalCounts << endl;
    }
  }

  // Finally add statistical fluctuation

  if (m_statisticalFluctuation) {
    for (int istage = 0; istage < Nstage; ++istage) {
      for (int idet = 0; idet < Ndetectors; ++idet) {
        for (int ibin = 0; ibin < m_nSamples; ++ibin) {
          m_positronDetectedSpectrum[istage][idet][ibin] =
              ((double)(ran->Poisson(
                  m_positronDetectedSpectrum[istage][idet][ibin] *
                  m_binWidth))) /
              m_binWidth;
        }
      }
    }
  }
}

void Spectrum::correctNonLinearity()
{
  for (int istage = 0; istage < Nstage; istage++) {
    for (int idet = 0; idet < Ndetectors; idet++) {
      // First: Apply detector non-linear energy scale
      double totalCounts = 0;
      if (m_useIhepNonLinearModel || m_useBcwNonLinearModel ||
          m_useLbnlNonLinearModel || m_useUnifiedNonLinearModel ||
          m_use2015NonLinearModel) {
        if (m_useIhepNonLinearModel) {
          // reset non-linear parameters
          nl_func->SetParameters(
              m_ihep_nl_par_nominal[0], m_ihep_nl_par_nominal[1],
              m_ihep_nl_par_nominal[2], m_ihep_nl_par_nominal[3],
              m_ihep_nl_par_nominal[4],
              m_abs_escale_nominal * m_rel_escale_nominal[idet], 0);
        } else if (m_useBcwNonLinearModel) {
          nl_func->SetParameters(
              m_bcw_elec_nl_par_nominal[0], m_bcw_elec_nl_par_nominal[1],
              m_bcw_elec_nl_par_nominal[2], m_bcw_elec_nl_par_nominal[3],
              m_bcw_elec_nl_par_nominal[4],
              m_abs_escale_nominal * m_rel_escale_nominal[idet], 0);


        } else if (m_useLbnlNonLinearModel) {
          nl_func->SetParameters(
              m_lbnl_nl_par_nominal[0], m_lbnl_nl_par_nominal[1],
              m_lbnl_nl_par_nominal[2],
              m_abs_escale_nominal * m_rel_escale_nominal[idet], 0);

        } else if (m_useUnifiedNonLinearModel > 0 ||
                   m_use2015NonLinearModel > 0) {
          for (Int_t i = 0; i < m_num_unified_nl_pars; i++) {
            nl_func->SetParameter(i, m_unified_nl_par_nominal[i]);
          }

          nl_func->SetParameter(m_num_unified_nl_pars,
                                m_abs_escale_nominal *
                                    m_rel_escale_nominal[idet]);
          nl_func->SetParameter(m_num_unified_nl_pars + 1, 0);
        }


        for (int idx = 0; idx < m_nSamples; idx++) {
          double energy_true = m_energy[idet][idx];
          double energy_det = nl_func->Eval(energy_true);

          int eDetIdx = (int)((energy_det - m_eMin) / m_binWidth);
          if (eDetIdx < 0 || eDetIdx >= m_nSamples) {
            m_positronNLCorrectedSpectrum[istage][idet][idx] = 0;
            m_bgNLCorrectedSpectrum[istage][idet][idx] = 0;
            continue;
          }

          double binScaling = nl_func->Derivative(energy_true);

          double dNdE = 0;
          int sign = 1;
          double dEdet = energy_det - (0.5 + eDetIdx) * m_binWidth;
          if (dEdet < 0)
            sign = -1;

          if (eDetIdx == m_nSamples - 1) {
            dNdE = m_positronDetectedSpectrum[istage][idet][eDetIdx];
          } else {
            dNdE =
                (TMath::Abs(dEdet) / m_binWidth) *
                    m_positronDetectedSpectrum[istage][idet][eDetIdx + sign] +
                (1 - TMath::Abs(dEdet) / m_binWidth) *
                    m_positronDetectedSpectrum[istage][idet][eDetIdx];
          }

          m_positronNLCorrectedSpectrum[istage][idet][idx] = dNdE * binScaling;
          totalCounts +=
              m_positronNLCorrectedSpectrum[istage][idet][idx] * m_binWidth;

          // correct background spectrum as well
          if (eDetIdx == m_nSamples - 1) {
            dNdE = m_bgDetectedSpectrum[istage][idet][eDetIdx];
          } else {
            dNdE = (TMath::Abs(dEdet) / m_binWidth) *
                       m_bgDetectedSpectrum[istage][idet][eDetIdx + sign] +
                   (1 - TMath::Abs(dEdet) / m_binWidth) *
                       m_bgDetectedSpectrum[istage][idet][eDetIdx];
          }

          m_bgNLCorrectedSpectrum[istage][idet][idx] = dNdE * binScaling;
        }


      } else {
        // Not supported at this moment!
      }
    }
  }
}


void Spectrum::loadDistances(const char* distancematrixname)
{
  string dummyLine;
  string thead;
  float d2, d1, l2, l1, l4, l3;
  //-->Distances
  cout << " Distances ++++++++++++++++++++++++++++++++++++++" << endl;
  ifstream disfile(distancematrixname);
  getline(disfile, dummyLine);
  while (disfile >> thead >> d1 >> d2 >> l1 >> l2 >> l3 >> l4) {
    cout << thead << "\t" << d1 << "\t" << d2 << "\t" << l1 << "\t" << l2
         << "\t" << l3 << "\t" << l4 << endl; // tmp
    int adnum = atoi(thead.substr(2, 1).c_str());
    // note: must convert into kms
    m_detectorDistance[adnum - 1][0] = d1 * 0.001;
    m_detectorDistance[adnum - 1][1] = d2 * 0.001;
    m_detectorDistance[adnum - 1][2] = l1 * 0.001;
    m_detectorDistance[adnum - 1][3] = l2 * 0.001;
    m_detectorDistance[adnum - 1][4] = l3 * 0.001;
    m_detectorDistance[adnum - 1][5] = l4 * 0.001;
  }
}

void Spectrum::loadBgSpecForToy(TString* accspecname, const Char_t* li9specname,
                                const Char_t* amcspecname,
                                const Char_t* fnspecname,
                                const Char_t* alnspecname)
{
  cout << "Loading bg spectra..." << endl;
  Char_t name[1024];
  TDirectory* dir = gDirectory;

  //(accidentals)
  for (int istage = 0; istage < Nstage; ++istage) {
    m_accspec[istage] = new TFile(accspecname[istage].Data(), "READ");
    dir->cd();
    for (int idet = 0; idet < Ndetectors; ++idet) {
      sprintf(name, "h_accidental_eprompt_fine_inclusive_eh%d_ad%d",
              detConfigEH[idet], detConfigAD[idet]);

      Char_t name2[1024];
      sprintf(name2, "CorrAccEvtsSpec_stage%i_eh%i_ad%d", istage + 1,
              detConfigEH[idet], detConfigAD[idet]);

      CorrAccEvtsSpec[istage][idet] =
          (TH1F*)m_accspec[istage]->Get(name)->Clone(name2);

      // cout << "tdper[iweek].ObsEvtsSpec[idet]->Integral(): " <<
      // tdper[0].ObsEvtsSpec[idet]->Integral() << endl;
    }
  }

  //-->I stopped here
  for (int istage = 0; istage < Nstage; ++istage) {
    for (int idet = 0; idet < Ndetectors; ++idet) {
      for (Int_t ibin = 0; ibin < CorrAccEvtsSpec[istage][idet]->GetNbinsX();
           ibin++) {
        CorrAccEvtsSpec[istage][idet]->SetBinError(ibin + 1, 0);
      }

      if (CorrAccEvtsSpec[istage][idet]->Integral() > 0)
        CorrAccEvtsSpec[istage][idet]->Scale(
            pred->tdper[istage].AccEvts[idet] /
            CorrAccEvtsSpec[istage][idet]->Integral());
    }
    // m_accspec[istage]->Close();
  }

  cout << "--> loaded accidental spectra" << endl;

  //(li9/he8)
  m_li9spec = new TFile(li9specname, "READ");
  dir->cd();
  for (int istage = 0; istage < Nstage; ++istage) {
    for (int idet = 0; idet < Ndetectors; ++idet) {
      sprintf(name, "CorrLi9EvtsSpec_ad%d", idet);
      CorrLi9EvtsSpec[istage][idet] =
          (TH1F*)m_li9spec->Get("h_nominal")->Clone(name);

      for (Int_t ibin = 0; ibin < CorrLi9EvtsSpec[istage][idet]->GetNbinsX();
           ibin++) {
        CorrLi9EvtsSpec[istage][idet]->SetBinError(ibin + 1, 0);
      }

      CorrLi9EvtsSpec[istage][idet]->Scale(
          pred->tdper[istage].Li9Evts[idet] /
          CorrLi9EvtsSpec[istage][idet]->Integral());
    }
  }
  // m_li9spec->Close();
  cout << "--> loaded Li9 spectra" << endl;

  //(amc)
  m_amcspec = new TFile(amcspecname, "READ");

  dir->cd();
  for (int istage = 0; istage < Nstage; ++istage) {
    for (int idet = 0; idet < Ndetectors; ++idet) {
      sprintf(name, "CorrAmcEvtsSpec_ad%d", idet);
      amcfunc = (TF1*)m_amcspec->Get("expo")->Clone();
      amcfunc->SetRange(0.7, 9.0);
      CorrAmcEvtsSpec[istage][idet] =
          (TH1F*)CorrLi9EvtsSpec[istage][idet]->Clone(name);
      CorrAmcEvtsSpec[istage][idet]->Reset();
      CorrAmcEvtsSpec[istage][idet]->Add(amcfunc);

      // old
      // CorrAmcEvtsSpec[idet] = (TH1F*)m_amcspec->Get("h_toy")->Clone(name);

      for (Int_t ibin = 0; ibin < CorrAmcEvtsSpec[istage][idet]->GetNbinsX();
           ibin++) {
        CorrAmcEvtsSpec[istage][idet]->SetBinError(ibin + 1, 0);
      }

      CorrAmcEvtsSpec[istage][idet]->Scale(
          pred->tdper[istage].AmcEvts[idet] /
          CorrAmcEvtsSpec[istage][idet]->Integral(
              1, CorrAmcEvtsSpec[istage][idet]->GetXaxis()->GetNbins()));
    }
  }
  // m_amcspec->Close();
  cout << "--> loaded AmC spectra" << endl;

  //(fn)
  m_fnspec = new TFile(fnspecname, "READ");
  dir->cd();
  for (int istage = 0; istage < Nstage; ++istage) {
    for (int idet = 0; idet < Ndetectors; ++idet) {
      sprintf(name, "CorrFnEvtsSpec_ad%d", idet);
      Char_t nameFn[1024];
      sprintf(nameFn, "h_%dAD_fn_fine", idet + 1);

      CorrFnEvtsSpec[istage][idet] = (TH1F*)m_fnspec->Get(nameFn)->Clone(name);

      for (Int_t ibin = 0; ibin < CorrFnEvtsSpec[istage][idet]->GetNbinsX();
           ibin++) {
        CorrFnEvtsSpec[istage][idet]->SetBinError(ibin + 1, 0);
      }
      CorrFnEvtsSpec[istage][idet]->Scale(
          pred->tdper[istage].FnEvts[idet] /
          CorrFnEvtsSpec[istage][idet]->Integral());
    }
  }
  // m_fnspec->Close();
  cout << "--> loaded fn spectra" << endl;

  //(aln)
  int AlphaAD[8] = {1, 2, 3, 8, 4, 5, 6, 7}; // Hack to get 8AD# aligned
  m_alnspec = new TFile(alnspecname, "READ");
  dir->cd();

  for (int istage = 0; istage < Nstage; ++istage) {
    for (int idet = 0; idet < Ndetectors; ++idet) {
      sprintf(name, "CorrAlnEvtsSpec_ad%d", idet);
      CorrAlnEvtsSpec[istage][idet] =
          (TH1F*)CorrLi9EvtsSpec[istage][idet]->Clone(name);
      CorrAlnEvtsSpec[istage][idet]->Reset();

      Char_t alnhistname[1024];
      sprintf(alnhistname, "AD%i;1",
              AlphaAD[idet]); //<-- the file has also a TCanvas with the same
                              // name (Ugly fix to use proper alpha-n)

      TH1F* htemp;
      htemp = (TH1F*)m_alnspec->Get(alnhistname)->Clone(name);

      for (Int_t ibin = 0; ibin < htemp->GetNbinsX(); ibin++) {
        // Only want spectrum if between 0.7 to 12 MeV
        if (htemp->GetXaxis()->GetBinCenter(ibin + 1) > 0.7 &&
            htemp->GetXaxis()->GetBinCenter(ibin + 1) < 12.0)
          CorrAlnEvtsSpec[istage][idet]->Fill(htemp->GetBinCenter(ibin + 1),
                                              htemp->GetBinContent(ibin + 1));
      }

      delete htemp;

      // assign zero errors
      for (Int_t ibin = 0; ibin < CorrAlnEvtsSpec[istage][idet]->GetNbinsX();
           ibin++) {
        CorrAlnEvtsSpec[istage][idet]->SetBinError(ibin + 1, 0);
      }

      CorrAlnEvtsSpec[istage][idet]->Scale(
          pred->tdper[istage].AlnEvts[idet] /
          CorrAlnEvtsSpec[istage][idet]->Integral());
    }
  }
  // m_alnspec->Close();
  cout << "--> loaded alpha-n spectra" << endl;

  cout << "done loading bg spectra!" << endl;

  // un-correcting the spectrum normalization....
  for (int istage = 0; istage < Nstage; ++istage) {
    for (int idet = 0; idet < Ndetectors; ++idet) {
      float factor = pred->tdper[istage].MuonVetoEff[idet] *
                     pred->tdper[istage].DMCEff[idet] *
                     pred->tdper[istage].Livetime[idet] *
                     // pred->tdper[istage].DelayedEff[idet] *
                     pred->tdper[istage].TargetMass[idet] /
                     pred->tdper[istage].TargetMass[0];

      CorrAccEvtsSpec[istage][idet]->Scale(factor);
      CorrAmcEvtsSpec[istage][idet]->Scale(factor);
      CorrLi9EvtsSpec[istage][idet]->Scale(factor);
      CorrFnEvtsSpec[istage][idet]->Scale(factor);
      CorrAlnEvtsSpec[istage][idet]->Scale(factor);
    }
  }
  /*
  //Make array
  for(int istage=0;istage<Nstage;++istage){
  for(int idet=0;idet<Ndetectors;++idet){
  assert(CorrAccEvtsSpec[istage][idet]->GetXaxis()->GetNbins()==m_nSamplesBkg);
  assert(CorrAmcEvtsSpec[istage][idet]->GetXaxis()->GetNbins()==m_nSamplesBkg);
  assert(CorrLi9EvtsSpec[istage][idet]->GetXaxis()->GetNbins()==m_nSamplesBkg);
  assert(CorrFnEvtsSpec[istage][idet]->GetXaxis()->GetNbins()==m_nSamplesBkg);
  assert(CorrAlnEvtsSpec[istage][idet]->GetXaxis()->GetNbins()==m_nSamplesBkg);
  }
  }*/

  for (int idet = 0; idet < Ndetectors; idet++) {
    for (int idx = 0; idx < m_nSamplesBkg; idx++) {
      m_energy_bkg[idet][idx] = (0.5 + idx) * m_binWidthBkg + m_eMin;
    }
  }

  cout << "updateBgDetected" << endl; // XXX
  this->updateBgDetected();
  cout << "END updateBgDetected" << endl; // XXX

} // end of LoadBgSpec

void Spectrum::setBgRemoveFlag(bool acc_flag, bool li9_flag, bool fn_flag,
                               bool amc_flag, bool aln_flag)
{
  m_removeAccBg = acc_flag;
  m_removeLi9Bg = li9_flag;
  m_removeFnBg = fn_flag;
  m_removeAmcBg = amc_flag;
  m_removeAlnBg = aln_flag;

  this->updateBgDetected();
}

void Spectrum::updateBgDetected()
{
  // Make a temporary copy of bg histograms which will have the fluctuations (do
  // not forget to delete or will have memory leak!)
  TH1F* CopyAccEvtsSpec[Nstage][Ndetectors];
  TH1F* CopyAmcEvtsSpec[Nstage][Ndetectors];
  TH1F* CopyAlnEvtsSpec[Nstage][Ndetectors];
  TH1F* CopyLi9EvtsSpec[Nstage][Ndetectors];
  TH1F* CopyFnEvtsSpec[Nstage][Ndetectors];
  Char_t name[1024];
  for (int istage = 0; istage < Nstage; ++istage) {
    for (int idet = 0; idet < Ndetectors; ++idet) {
      sprintf(name, "h_acc_copy_%i", idet);
      CopyAccEvtsSpec[istage][idet] =
          (TH1F*)CorrAccEvtsSpec[istage][idet]->Clone(name);
      sprintf(name, "h_amc_copy_%i", idet);
      CopyAmcEvtsSpec[istage][idet] =
          (TH1F*)CorrAmcEvtsSpec[istage][idet]->Clone(name);
      sprintf(name, "h_li9_copy_%i", idet);
      CopyLi9EvtsSpec[istage][idet] =
          (TH1F*)CorrLi9EvtsSpec[istage][idet]->Clone(name);
      sprintf(name, "h_fn_copy_%i", idet);
      CopyFnEvtsSpec[istage][idet] =
          (TH1F*)CorrFnEvtsSpec[istage][idet]->Clone(name);
      sprintf(name, "h_aln_copy_%i", idet);
      CopyAlnEvtsSpec[istage][idet] =
          (TH1F*)CorrAlnEvtsSpec[istage][idet]->Clone(name);
    }
  }
  // ****************************************************
  // Variation (i.e. scale whole bg spectra up/down by systematic error)
  // ****************************************************
  // define and initialize
  double scale_factor_acc[Nstage][Ndetectors];
  double scale_factor_amc[Nstage][Ndetectors];
  double scale_factor_li9[Nstage][Ndetectors];
  double scale_factor_fn[Nstage][Ndetectors];
  double scale_factor_aln[Nstage][Ndetectors];
  double corr_random[Nhalls] = {0};

  for (int istage = 0; istage < Nstage; ++istage) {
    for (int idet = 0; idet < Ndetectors; ++idet) {
      scale_factor_acc[istage][idet] = 1;
      scale_factor_amc[istage][idet] = 1;
      scale_factor_li9[istage][idet] = 1;
      scale_factor_fn[istage][idet] = 1;
      scale_factor_aln[istage][idet] = 1;
    }
  }
  // Fluctuate each background depending on how it makes the most sense
  // 1) For accidentals, fluctuate each AD independenlty
  // 2) For li9 and fn, fluctuate ADs in same site with same factor
  // 3) For amc, fluctuate all ADs in a correlated way
  // 4) For alpha-n, for now fluctuating all ADs independently

  // Accidental and alpha-n are entirely uncorrelated
  for (int istage = 0; istage < Nstage; ++istage) {
    for (int idet = 0; idet < Ndetectors; ++idet) {
      if (m_varyAccBg > 0) {
        if (pred->tdper[istage].AccEvts[idet] > 0)
          scale_factor_acc[istage][idet] =
              (1 + pred->tdper[istage].AccErr[idet] * 1. /
                       pred->tdper[istage].AccEvts[idet] * ran->Gaus(0, 1));
        else
          scale_factor_acc[istage][idet] = 0;
      }

      if (m_varyAlnBg > 0) {
        if (pred->tdper[istage].AlnEvts[idet] > 0)
          scale_factor_aln[istage][idet] =
              (1 + pred->tdper[istage].AlnErr[idet] * 1. /
                       pred->tdper[istage].AlnEvts[idet] * ran->Gaus(0, 1));
        else
          scale_factor_aln[istage][idet] = 0;
      }
    }
  }

  // amc are treated completely correlated
  corr_random[0] = ran->Gaus(0, 1);
  for (int istage = 0; istage < Nstage; ++istage) {
    for (int idet = 0; idet < Ndetectors; ++idet) {
      if (m_varyAmcBg > 0) {
        if (pred->tdper[istage].AmcEvts[idet] > 0)
          scale_factor_amc[istage][idet] =
              (1 + pred->tdper[istage].AmcErr[idet] * 1. /
                       pred->tdper[istage].AmcEvts[idet] * corr_random[0]);
        else
          scale_factor_amc[istage][idet] = 0;
      }
    }
  }

  // li are correlated between detectors in the same hall
  for (int ihall = 0; ihall < Nhalls; ++ihall) {
    corr_random[ihall] = ran->Gaus(0, 1);
  }

  for (int istage = 0; istage < Nstage; ++istage) {
    for (int idet = 0; idet < Ndetectors; ++idet) {
      if (m_varyLi9Bg > 0) {
        if (pred->tdper[istage].Li9Evts[idet] > 0)
          scale_factor_li9[istage][idet] =
              (1 + pred->tdper[istage].Li9Err[idet] * 1. /
                       pred->tdper[istage].Li9Evts[idet] *
                       corr_random[detConfigEH[idet] - 1]);
        else
          scale_factor_li9[istage][idet] = 0;
      }
    }
  }

  // fn are correlated between detectors in the same hall
  for (int ihall = 0; ihall < Nhalls; ++ihall) {
    corr_random[ihall] = ran->Gaus(0, 1);
  }

  for (int istage = 0; istage < Nstage; ++istage) {
    for (int idet = 0; idet < Ndetectors; ++idet) {
      if (m_varyFnBg > 0) {
        if (pred->tdper[istage].FnEvts[idet] > 0)
          scale_factor_fn[istage][idet] =
              (1 + pred->tdper[istage].FnErr[idet] * 1. /
                       pred->tdper[istage].FnEvts[idet] *
                       corr_random[detConfigEH[idet] - 1]);
        else
          scale_factor_fn[istage][idet] = 0;
      }
    }
  }

  // force set background contribution to 0 if the flag is set
  for (int istage = 0; istage < Nstage; ++istage) {
    for (int idet = 0; idet < Ndetectors; ++idet) {
      if (m_removeAccBg)
        scale_factor_acc[istage][idet] = 0;
      if (m_removeAlnBg)
        scale_factor_aln[istage][idet] = 0;
      if (m_removeAmcBg)
        scale_factor_amc[istage][idet] = 0;
      if (m_removeLi9Bg)
        scale_factor_li9[istage][idet] = 0;
      if (m_removeFnBg)
        scale_factor_fn[istage][idet] = 0;
    }
  }

  // *******************************************************
  // Distortion (i.e. distorting shape without affecting rate)
  // concerning correlations:
  // - accidentals: fluctuate shape of all ADs together
  // - amc: fluctuate shape of all ADs together
  // - li9: fluctuate shape of all ADs together
  // - fn: fluctuate shape of ADs in same site together
  // - aln: fluctuate shape of all ADs together
  // *******************************************************
  if (m_distortAccBg > 0) {
    TF1* func_acc = getDistortionCurve(m_distortAccBg);
    for (int istage = 0; istage < Nstage; ++istage) {
      for (int idet = 0; idet < Ndetectors; ++idet) {
        CopyAccEvtsSpec[istage][idet]->Multiply(func_acc);
        if (CopyAccEvtsSpec[istage][idet]->Integral() > 0)
          CopyAccEvtsSpec[istage][idet]->Scale(
              CorrAccEvtsSpec[istage][idet]->Integral() * 1. /
              CopyAccEvtsSpec[istage][idet]->Integral());
      }
    }
    delete func_acc;
  }

  if (m_distortAmcBg > 0) {
    TF1* amccopy = (TF1*)amcfunc->Clone("amccopy");
    double invpar1 =
        1. / amcfunc->GetParameter(1) * ran->Gaus(1, m_distortAmcBg);
    amccopy->SetParameter(1, 1. / invpar1);
    for (int istage = 0; istage < Nstage; ++istage) {
      for (int idet = 0; idet < Ndetectors; ++idet) {
        CopyAmcEvtsSpec[istage][idet]->Reset();
        CopyAmcEvtsSpec[istage][idet]->Add(amccopy);
        if (CopyAmcEvtsSpec[istage][idet]->Integral() > 0)
          CopyAmcEvtsSpec[istage][idet]->Scale(
              CorrAmcEvtsSpec[istage][idet]->Integral() * 1. /
              CopyAmcEvtsSpec[istage][idet]->Integral());
      }
    }
    delete amccopy;
  }

  if (m_distortLi9Bg != "null") {
    int ientry = int(ran->Uniform(0, m_entries_distortLi9Bg));
    TH1F* func_li;
    m_tree_distortLi9Bg->GetEntry(ientry);
    func_li = (TH1F*)m_func_distortLi9Bg->Clone();

    for (int istage = 0; istage < Nstage; ++istage) {
      for (int idet = 0; idet < Ndetectors; ++idet) {
        CopyLi9EvtsSpec[istage][idet]->Multiply(func_li);
        if (CopyLi9EvtsSpec[istage][idet]->Integral() > 0)
          CopyLi9EvtsSpec[istage][idet]->Scale(
              CorrLi9EvtsSpec[istage][idet]->Integral() * 1. /
              CopyLi9EvtsSpec[istage][idet]->Integral());
      }
    }
    delete func_li;
  }

  if (m_distortFnBg > 0) {
    // eh1
    TF1* func_fn = getFNDistortionCurve(m_distortFnBg);
    for (int istage = 0; istage < Nstage; ++istage) {
      for (int idet = 0; idet < 2; ++idet) {
        CopyFnEvtsSpec[istage][idet]->Multiply(func_fn);
        if (CopyFnEvtsSpec[istage][idet]->Integral() > 0)
          CopyFnEvtsSpec[istage][idet]->Scale(
              CorrFnEvtsSpec[istage][idet]->Integral() * 1. /
              CopyFnEvtsSpec[istage][idet]->Integral());
      }
    }
    delete func_fn;
    // eh2
    func_fn = getFNDistortionCurve(m_distortFnBg);
    for (int istage = 0; istage < Nstage; ++istage) {
      for (int idet = 2; idet < 4; ++idet) {
        CopyFnEvtsSpec[istage][idet]->Multiply(func_fn);
        if (CopyFnEvtsSpec[istage][idet]->Integral() > 0)
          CopyFnEvtsSpec[istage][idet]->Scale(
              CorrFnEvtsSpec[istage][idet]->Integral() * 1. /
              CopyFnEvtsSpec[istage][idet]->Integral());
      }
    }
    delete func_fn;
    // eh3
    func_fn = getFNDistortionCurve(m_distortFnBg);
    for (int istage = 0; istage < Nstage; ++istage) {
      for (int idet = 4; idet < 8; ++idet) {
        CopyFnEvtsSpec[istage][idet]->Multiply(func_fn);
        if (CopyFnEvtsSpec[istage][idet]->Integral() > 0)
          CopyFnEvtsSpec[istage][idet]->Scale(
              CorrFnEvtsSpec[istage][idet]->Integral() * 1. /
              CopyFnEvtsSpec[istage][idet]->Integral());
      }
    }
    delete func_fn;
  }

  if (m_distortAlnBg > 0) {
    TF1* func_aln = getDistortionCurve(m_distortAlnBg);
    for (int istage = 0; istage < Nstage; ++istage) {
      for (int idet = 0; idet < Ndetectors; ++idet) {
        CopyAlnEvtsSpec[istage][idet]->Multiply(func_aln);
        if (CopyAlnEvtsSpec[istage][idet]->Integral() > 0)
          CopyAlnEvtsSpec[istage][idet]->Scale(
              CorrAlnEvtsSpec[istage][idet]->Integral() * 1. /
              CopyAlnEvtsSpec[istage][idet]->Integral());
      }
    }
    delete func_aln;
  }


  // ********************************************************
  // Returning the bg spectra
  // ********************************************************

  for (int istage = 0; istage < Nstage; ++istage) {
    for (int idet = 0; idet < Ndetectors; ++idet) {
      for (int ibin = 0; ibin < m_nSamplesBkg; ++ibin) {
        m_bgDetectedSpectrum[istage][idet][ibin] = 0;
        m_bgDetectedSpectrum[istage][idet][ibin] +=
            scale_factor_acc[istage][idet] *
            CopyAccEvtsSpec[istage][idet]->GetBinContent(ibin + 1);
        m_bgDetectedSpectrum[istage][idet][ibin] +=
            scale_factor_amc[istage][idet] *
            CopyAmcEvtsSpec[istage][idet]->GetBinContent(ibin + 1);
        m_bgDetectedSpectrum[istage][idet][ibin] +=
            scale_factor_li9[istage][idet] *
            CopyLi9EvtsSpec[istage][idet]->GetBinContent(ibin + 1);
        m_bgDetectedSpectrum[istage][idet][ibin] +=
            scale_factor_fn[istage][idet] *
            CopyFnEvtsSpec[istage][idet]->GetBinContent(ibin + 1);
        m_bgDetectedSpectrum[istage][idet][ibin] +=
            scale_factor_aln[istage][idet] *
            CopyAlnEvtsSpec[istage][idet]->GetBinContent(ibin + 1);
      }

      delete CopyAccEvtsSpec[istage][idet];
      delete CopyAmcEvtsSpec[istage][idet];
      delete CopyLi9EvtsSpec[istage][idet];
      delete CopyFnEvtsSpec[istage][idet];
      delete CopyAlnEvtsSpec[istage][idet];
    }
  }

  if (m_statisticalFluctuation) {
    for (int istage = 0; istage < Nstage; ++istage) {
      for (int idet = 0; idet < Ndetectors; ++idet) {
        for (int ibin = 0; ibin < m_nSamplesBkg; ++ibin) {
          m_bgDetectedSpectrum[istage][idet][ibin] =
              (double)(ran->Poisson(m_bgDetectedSpectrum[istage][idet][ibin]));
          //          =
          //          ((double)(ran->Poisson(m_bgDetectedSpectrum[idet][ibin]*m_binWidth)))/m_binWidth;;
        }
      }
    }
  }
}


void Spectrum::loadIavCorrection(const char* iavcorrectionname)
{
  m_iavCorrFile = new TFile(iavcorrectionname);
  TH2F* Correction;
  Correction = (TH2F*)m_iavCorrFile->Get("Correction_LS")->Clone();

  cout << "reading Iav correction file from " << iavcorrectionname << endl;
  for (int i = 0; i < 240;
       i++) { // i: true positron energy bin; j: distorted energy bin
    // Total events in input spectrum bin i
    m_iav_frac_orig[i] = 0;
    if (Correction->Integral(i + 1, i + 1, 0, -1) > 0) {
      for (int j = 0; j < i + 1; j++) {
        m_iav_corr_orig[i][j] = Correction->GetBinContent(i + 1, j + 1) /
                                Correction->Integral(i + 1, i + 1, 0, -1);
        if (i != j) {
          m_iav_frac_orig[i] += m_iav_corr_orig[i][j];
        }
      }
    } else {
      for (int j = 0; j < i + 1; j++) {
        if (i == j) {
          m_iav_corr_orig[i][j] = 1;
        } else {
          m_iav_corr_orig[i][j] = 0;
        }
      }
      m_iav_frac_orig[i] = 0;
    }
  }

  delete Correction;

  /*Remove IAV here
    for(int i=0;i<240;i++){ // i: true positron energy bin; j: distorted energy
    bin for(int j=0;j<i+1-11;j++){

    m_iav_corr_orig[i][j] = 0;

    }

    }*/

  Double_t rescale_fac = 240.0 / (Double_t)m_nSamples;

  // rescale IAV correction into finer binning
  for (int i = 0; i < m_nSamples;
       i++) { // i: true positron energy bin; j: distorted energy bin
    int i_orig = (int)(rescale_fac * i);
    m_nominal_iav_frac[i] = m_iav_frac_orig[i_orig];
    for (int j = 0; j < m_nSamples; j++) {
      int j_orig = (int)(rescale_fac * j);
      if (i_orig == j_orig) {
        if (i == j) {
          m_nominal_iav_corr[i][j] = m_iav_corr_orig[i_orig][i_orig];
        }
      } else {
        m_nominal_iav_corr[i][j] =
            m_iav_corr_orig[i_orig][j_orig] * rescale_fac;
      }
    }
  }

  for (int idet = 0; idet < Ndetectors; idet++) {
    for (int i = 0; i < m_nSamples;
         i++) { // i: true positron energy bin; j: distorted energy bin
      m_iav_frac[idet][i] = m_nominal_iav_frac[i];
      for (int j = 0; j < m_nSamples; j++) {
        m_iav_corr[idet][i][j] = m_nominal_iav_corr[i][j];
      }
    }
  }

  // TCanvas * ccc = new TCanvas ("ccc","ccc");
  // TH2D * h_iav_corr = new
  // TH2D("h_iav_corr","h_iav_corr",m_nSamples,0,12,m_nSamples,0,12); for(int
  // i=0;i<m_nSamples;i++){ // i: true positron energy bin; j: distorted energy
  // bin
  //   for(int j=0;j<m_nSamples;j++){
  //     h_iav_corr->SetBinContent(i+1,j+1,m_iav_corr[0][i][j]);
  //   }
  // }
  // h_iav_corr->Draw("zcol");
  // ccc->Update();
}


void Spectrum::initialize(DataSet* data)
{
  //  Setup isotope antinu spectra
  m_xsec = new CrossSectionTable();
  m_xsec->load(data->getString("crossSectionFilename"));


  m_useBcwFluxUncertainty = data->getDouble("useBcwFluxUncertainty");
  m_useAbInitioSpectra = data->getDouble("useAbInitioSpectra");
  if (m_useAbInitioSpectra > 0)
    UseChristineModel = false;


  if (UseChristineModel) {
    Bool_t decomp_result;
    for (Int_t istage = 0; istage < Nstage; ++istage) {
      m_corespectrum[istage] = new CoreSpectrum();

      decomp_result = m_corespectrum[istage]->loadSpectra(
          Paths::reactor_spectrum(istage),
          Paths::reactor_covmatrix());


      if (!decomp_result) {
        cout << "Failed to load flux matrix! Exiting" << endl;
        exit(0);
      }


      m_corespectrum[istage]->loadCrossSectionTable(
          data->getString("crossSectionFilename"));

      // YN: Reading BCW flux uncertainty

      if (m_useBcwFluxUncertainty > 0) {
        cout << "Loading BCW flux covariance matrix at " << Paths::bcw_flux()
             << endl;
        Bool_t decomp_result_bcw =
          m_corespectrum[istage]->loadCovMatrixBCW(Paths::bcw_flux());
        if (!decomp_result_bcw) {
          cout << "Failed to load BCW flux matrix! Exiting" << endl;
          exit(0);
        }
      }
    }

  } else if (m_useAbInitioSpectra) {
    for (Int_t istage = 0; istage < Nstage; ++istage) {
      m_corespectrum[istage] = new CoreSpectrum();
      m_corespectrum[istage]->loadAbInitioSpectra(Paths::abinitio_spectra());
      m_corespectrum[istage]->loadCrossSectionTable(
          data->getString("crossSectionFilename"));
    }
  }


  /*else{
    for (Int_t istage=0;istage<Nstage;++istage){
    for (int icore = 0; icore < Ncores; icore++){
    m_isotopes[istage][icore] = new IsotopeTable();
    m_isotopes[istage][icore]->loadIsotopes(data->getString(Form("fissionIsotopeFilename%d",icore+1)));
    m_isotopes[istage][icore]->loadIsotopeCovMatrix(data->getString("fissionIsotopeCovMatrixFilename"));
    m_isotopes[istage][icore]->loadSpectra(data->getString("fissionSpectraFilename"));
    }
    }
    }*/

  m_eMin = data->getDouble("minimumEnergy");
  m_eMax = data->getDouble("maximumEnergy");
  m_nSamples = int(data->getDouble("nSpectrumSamples"));
  m_binWidth = (m_eMax - m_eMin) / m_nSamples;
  m_nSamplesBkg = int(data->getDouble("nSpectrumSamplesBkg"));
  m_binWidthBkg = (m_eMax - m_eMin) / m_nSamplesBkg;
  m_alpha = data->getDouble("detectorNonlinearAlpha");
  m_alpha_nominal = m_alpha;
  m_alpha_error = data->getDouble("detectorNonlinearAlphaError");
  //  m_beta = data->getDouble("detectorLinearBeta"); // not used for now.....


  m_abs_escale = data->getDouble("detectorAbsoluteEnergyScale");
  m_abs_escale_nominal = m_abs_escale;
  m_abs_escale_error = data->getDouble("detectorAbsoluteEnergyScaleError");


  loadIavCorrection(data->getString("iavCorrectionFilename"));

  // FIXME: hard_coded relative energy scale uncertainty
  for (Int_t idet = 0; idet < Ndetectors; idet++) {
    m_rel_escale[idet] = 1.0;
    m_rel_escale_error[idet] = 0.0020; // 0.20% (May 2014)
    m_rel_escale_nominal[idet] = m_rel_escale[idet];
  }

  // FIXME: hard_coded energy offset uncertainty (currently, set to 0.08 MeV)
  m_abs_eoffset = 0.0;
  m_abs_eoffset_error = 0.08; //  (MeV)

  for (Int_t istage = 0; istage < Nstage; ++istage) {
    for (Int_t idet = 0; idet < Ndetectors; idet++) {
      m_rel_eoffset[istage][idet] = 0.0;
    }
  }
  m_rel_eoffset_error = 0.013; //  (MeV)


  m_iav_error = data->getDouble("detectorIavThicknessError");

  m_detectorResolution = data->getDouble("detectorResolution") / 100.; // Pct
  m_detectorResolution_error =
      data->getDouble("detectorResolutionError") / 100.; // Pct
  m_detectorResolution_error_uncorr =
      data->getDouble("detectorResolutionErrorUncorr") / 100.; // Pct

  m_detectorResolution_nominal = m_detectorResolution;

  // m_detectorEfficiency = data->getDouble("detectorEfficiency");
  // m_runningTime = data->getDouble("runningTime");
  // m_detectorSize = data->getDouble("detectorSize");
  m_targetProtonsPerKton = data->getDouble("nominalTargetProtons");
  // m_reactorPower = data->getDouble("reactorPower");
  m_sinSq2Th12_nominal = data->getDouble("sinSq2Theta12");
  m_sinSq2Th12_error = data->getDouble("sinSq2Theta12Err");
  m_sinSq2Th14_nominal = data->getDouble("sinSq2Theta14");
  m_sinSq2Th14_error = data->getDouble("sinSq2Theta14Err");
  m_sinSq2Th13 = data->getDouble("sinSq2Theta13");
  m_deltaMSq21_nominal = data->getDouble("deltaMSq21");
  m_deltaMSq21_error = data->getDouble("deltaMSq21Err");
  m_deltaMSqee_nominal = data->getDouble("deltaMSqee");
  m_deltaMSqee_error = data->getDouble("deltaMSqeeErr");
  m_deltaMSq41_nominal = data->getDouble("deltaMSq41");
  m_deltaMSq41_error = data->getDouble("deltaMSq41Err");

  m_deltaMSqee = m_deltaMSqee_nominal;
  m_deltaMSq41 = m_deltaMSq41_nominal;
  m_deltaMSq21 = m_deltaMSq21_nominal;
  m_sinSq2Th12 = m_sinSq2Th12_nominal;
  m_sinSq2Th14 = m_sinSq2Th14_nominal;

  m_randomizeSolarOscPars = data->getDouble("randomizeSolarOscPars");
  m_randomizeDm2ee = data->getDouble("randomizeDm2ee");

  m_randomizeIsotopeFraction = data->getDouble("randomizeIsotopeFraction");
  m_randomizeReactorPower = data->getDouble("randomizeReactorPower");
  m_randomizeIavDistortion = data->getDouble("randomizeIavDistortion");
  m_randomizeScintiNonLinear = data->getDouble("randomizeScintiNonLinear");
  m_randomizeElecNonLinear = data->getDouble("randomizeElecNonLinear");
  m_randomizeResolution = data->getDouble("randomizeResolution");
  m_randomizeRelativeEnergyScale =
      data->getDouble("randomizeRelativeEnergyScale");
  m_randomizeAbsoluteEnergyScale =
      data->getDouble("randomizeAbsoluteEnergyScale");

  m_randomizeRelativeEnergyOffset =
      data->getDouble("randomizeRelativeEnergyOffset");
  m_randomizeAbsoluteEnergyOffset =
      data->getDouble("randomizeAbsoluteEnergyOffset");

  m_randomizeDetectorEfficiency =
      data->getDouble("randomizeDetectorEfficiency");

  // randomize antinu spectra from each core based on Christine's covariance
  // matrix
  m_randomizeCoreSpectra = data->getDouble("randomizeCoreSpectra");

  m_varyAccBg = data->getDouble("varyAccBg");
  m_varyAmcBg = data->getDouble("varyAmcBg");
  m_varyFnBg = data->getDouble("varyFnBg");
  m_varyLi9Bg = data->getDouble("varyLi9Bg");
  m_varyAlnBg = data->getDouble("varyAlnBg");

  m_distortAccBg = data->getDouble("distortAccBg");
  m_distortAmcBg = data->getDouble("distortAmcBg");
  m_distortLi9Bg = data->getString("distortLi9Bg");
  if (m_distortLi9Bg != "null") {
    m_file_distortLi9Bg = new TFile(m_distortLi9Bg.c_str(), "READ");
    m_tree_distortLi9Bg = (TTree*)m_file_distortLi9Bg->Get("tr_distort");
    m_entries_distortLi9Bg = m_tree_distortLi9Bg->GetEntries();
    cout << "The distortion tree has " << m_entries_distortLi9Bg << " entries"
         << endl;
    m_tree_distortLi9Bg->SetBranchAddress("h_distort", &m_func_distortLi9Bg);
  }
  m_distortFnBg = data->getDouble("distortFnBg");
  m_distortAlnBg = data->getDouble("distortAlnBg");

  m_statisticalFluctuation = data->getDouble("statisticalFluctuation");

  m_useIhepNonLinearModel = data->getDouble("useIhepNonLinearModel");

  m_useBcwNonLinearModel = data->getDouble("useBcwNonLinearModel");

  m_useLbnlNonLinearModel = data->getDouble("useLbnlNonLinearModel");

  m_useUnifiedNonLinearModel = data->getDouble("useUnifiedNonLinearModel");

  m_use2015NonLinearModel = data->getDouble("use2015NonLinearModel");

  m_correlateUnifiedNonlinearPars =
      data->getDouble("correlateUnifiedNonlinearPars");


  m_osccalc->SetDeltaM2_ee(m_deltaMSqee);
  m_osccalc->SetDeltaM2_21(m_deltaMSq21);
  //  m_osccalc->SetTheta13(TMath::ASin(sqrt(m_sinSq2Th13))*0.5);
  m_osccalc->SetTheta12(TMath::ASin(sqrt(m_sinSq2Th12)) * 0.5);

  m_osccalc->SetDeltaM2_41(m_deltaMSq41);
  m_osccalc->SetTheta14(TMath::ASin(sqrt(m_sinSq2Th14)) * 0.5);


  if (m_useIhepNonLinearModel > 0) {
    cout << "Using IHEP non-linearity model." << endl;


    // FIXME: hard-coded IHEP non-linear function
    // References: Doc-8942 (LiangJiang) and Doc-8493 (Soeren), which include
    // electronics non-linearlity Parameters are taken from Soeren's note and
    // assume using the peak charge Currently, the correlateion between the
    // parameters are neglected, but can be implemented easily.

    // parameter for the default charge
    m_ihep_nl_par_nominal[0] = 1.061;
    m_ihep_nl_par_nominal[1] = 0.17;
    m_ihep_nl_par_nominal[2] = 1.81;
    m_ihep_nl_par_nominal[3] = 0.0000;
    m_ihep_nl_par_nominal[4] =
        0.933; // magic factor for 0.511 MeV gamma non-linearlity factor, read
               // from Soren's slides

    // parameter for the peak charge
    // m_ihep_nl_par_nominal[0] = 1.078;
    // m_ihep_nl_par_nominal[1] = 0.24;
    // m_ihep_nl_par_nominal[2] = 1.94;
    // m_ihep_nl_par_nominal[3] = -0.0005;
    // m_ihep_nl_par_nominal[4] = 0.90;// magic factor for 0.511 MeV gamma
    // non-linearlity factor, read from Soren's slides

    m_ihep_nl_par_error[0] = 0;
    m_ihep_nl_par_error[1] = 0.08; // super preriminary symetric uncertainty
    m_ihep_nl_par_error[2] = 0.8;  // super preriminary symetric uncertainty
    m_ihep_nl_par_error[3] = 0.0;
    m_ihep_nl_par_error[4] = 0.01; // from 68Ge energy scale uncertainty

    // hard-coded covariance matrix for IHEP non-linear model.

    for (Int_t i = 0; i < 5; i++) {
      for (Int_t j = 0; j < 5; j++) {
        if (i == j)
          m_ihep_nl_par_covmatrix[i][j] = 1;
        else
          m_ihep_nl_par_covmatrix[i][j] = 0;
      }
    }

    m_ihep_nl_par_covmatrix[1][2] = 1 / sqrt(2);
    m_ihep_nl_par_covmatrix[2][1] =
        1 / sqrt(2); // manually put the covariance between p1 and p2


    TMatrixD covmatrix(5, 5, &m_ihep_nl_par_covmatrix[0][0]);

    // cout << "Covariance matrix for IHEP non-liner parameters:" << endl;
    // covmatrix.Print();

    TDecompChol chol(covmatrix);
    chol.Decompose();

    TMatrixD cmat(chol.GetU());
    //  cmat.Print();
    TMatrixD tcmat(cmat.Transpose(cmat));
    // std::cout << "Cholesky matrix---------" << std::endl;
    // tcmat.Print();

    double* tmp_matrix = tcmat.GetMatrixArray();

    for (int i = 0; i < 5; i++) {
      for (int j = 0; j < 5; j++) {
        m_ihep_nl_par_covmatrix_l[i][j] = tmp_matrix[i * 5 + j];
        //      cout << "\t" << m_ihep_nl_par_covmatrix_l[i][j];
      }
      //    cout << endl;
    }


    for (Int_t i = 0; i < 5; i++) {
      m_ihep_nl_par[i] = m_ihep_nl_par_nominal[i];
    }

    ///////////// End of IHEP non-linear parameter setting //////////////////

    nl_func = new TF1("nl_func", nl_func_ihep, 0, 12, 7);
    nl_func->SetParameters(m_ihep_nl_par[0], m_ihep_nl_par[1], m_ihep_nl_par[2],
                           m_ihep_nl_par[3], m_ihep_nl_par[4], 1.0, 0.0);


  } else if (m_useBcwNonLinearModel > 0) {
    cout << "Using BCW non-linearity model." << endl;

    ///////////////////////////////////////////////////////////////////

    /// Settting for the new BCW non-linearity function. Feb. 20, 2013
    /// Updated to the new BCW model (Mar.26, 2013)
    /// Further updated to the "final" BCW model (Apr. 1, 2013)


    ifstream bcw_positron_data(Paths::bcw_positron_data());

    for (Int_t i = 0; i < n_bcw_positron_nl; i++) {
      bcw_positron_data >> m_bcw_positron_nl_e[i] >> m_bcw_positron_nl_fac[i];
    }
    bcw_positron_data.close();

    ifstream bcw_elec_data(Paths::bcw_elec_data());

    for (Int_t i = 0; i < 3; i++) {
      bcw_elec_data >> m_bcw_elec_nl_par_nominal[i] >>
          m_bcw_elec_nl_par_error[i];
      m_bcw_elec_nl_par[i] = m_bcw_elec_nl_par_nominal[i];
    }
    bcw_elec_data.close();

    TFile* bcw_ele_err_file = new TFile(Paths::bcw_ele_err());
    g_bcw_elec_nl_error[0] = (TGraph*)bcw_ele_err_file->Get("g_up")->Clone();
    g_bcw_elec_nl_error[1] =
        (TGraph*)bcw_ele_err_file->Get("g_down")->Clone();
    bcw_ele_err_file->Close();

    for (Int_t i = 3; i < 5; i++) { // those describe additional uncertainty for
                                    // bcw model that are described in Doc-XXXX
      m_bcw_elec_nl_par[i] = 0;
      m_bcw_elec_nl_par[i] = m_bcw_elec_nl_par_nominal[i];
      m_bcw_elec_nl_par_error[i] = 1;
    }


    ///////////// End of BCW non-linear parameter setting //////////////////


    nl_func = new TF1("nl_func", this, &Spectrum::nl_func_bcw, 0, 12, 7,
                      "Spectrum", "nl_func_bcw");
    nl_func->SetParameters(m_bcw_elec_nl_par[0], m_bcw_elec_nl_par[1],
                           m_bcw_elec_nl_par[2], m_bcw_elec_nl_par[3],
                           m_bcw_elec_nl_par[4], 1.0, 0.0);


  } else if (m_useLbnlNonLinearModel > 0) {
    cout << "Using LBNL non-linearity model." << endl;


    /// Settting for the LBNL non-linearity function. Mar. 24, 2013

    ifstream lbnl_positron_data(Paths::lbnl_positron_data());
    if (!lbnl_positron_data.is_open()) {
      cout << "Error: cannot find LBNL non-linearity curve!!!" << endl;
      exit(0);
    }

    Double_t total_err;

    for (Int_t i = 0; i < 300; i++) {
      lbnl_positron_data >> m_lbnl_positron_nl_e[i] >>
          m_lbnl_positron_nl_fac[i] >> total_err >>
          m_lbnl_positron_nl_err[0][i] >> m_lbnl_positron_nl_err[1][i] >>
          m_lbnl_positron_nl_err[2][i];
      // cout << "\t" << m_lbnl_positron_nl_e[i] << "\t"
      // <<m_lbnl_positron_nl_fac[i] << "\t" <<total_err
      //      << "\t" <<m_lbnl_positron_nl_err[0][i]
      //      << "\t" <<m_lbnl_positron_nl_err[1][i]
      //      << "\t" <<m_lbnl_positron_nl_err[2][i] << endl;
    }
    for (Int_t i = 0; i < 3; i++) {
      m_lbnl_nl_par[i] = 0;
      m_lbnl_nl_par_nominal[i] = 0;
      m_lbnl_nl_par_error[i] = 1.0;
    }

    ///////////// End of LBNL non-linear parameter setting //////////////////


    nl_func = new TF1("nl_func", nl_func_lbnl, 0, 12, 5);
    nl_func->SetParameters(m_lbnl_nl_par[0], m_lbnl_nl_par[1], m_lbnl_nl_par[2],
                           1.0, 0.0);


  } else if (m_useUnifiedNonLinearModel > 0 || m_use2015NonLinearModel > 0) {
    if (m_useUnifiedNonLinearModel > 0)
      cout << "Using UNIFIED non-linearity model." << endl;
    else
      cout << "Using UNIFIED non-linearity model (2015 update)." << endl;


    /// Settting for the unified Non-linearity model summarized by Soeren.


    TString unified_nl_filename;

    m_num_unified_nl_pars = 4; // 5 marginal curves in the final inflated model

    TString nominal_graph_name;
    //    TString pull_graph_name[m_num_unified_nl_pars];
    TString pull_graph_name[10];


    if (m_use2015NonLinearModel > 0) {
      unified_nl_filename = Paths::unified_nl();
      nominal_graph_name = "nominal";
      for (Int_t i = 0; i < m_num_unified_nl_pars; i++) {
        pull_graph_name[i] = LeakStr("pull%d", i);
      }
    } else {
      unified_nl_filename = Paths::unified_nl_final();
      nominal_graph_name = "positron_0";
      for (Int_t i = 0; i < m_num_unified_nl_pars; i++) {
        pull_graph_name[i] = LeakStr("positron_%d", i + 1);
      }
    }


    for (Int_t i = 0; i < m_num_unified_nl_pars; i++) {
      m_unified_nl_par[i] = 0.0;
      m_unified_nl_par_nominal[i] = 0.0;
      m_unified_nl_par_error[i] = 1.0;
    }

    TGraph* g_unified_positron_nl;
    TGraph* g_unified_positron_nl_pulls[10];

    TFile* unified_nl_file = new TFile(unified_nl_filename.Data());

    if (!unified_nl_file->IsOpen()) {
      cout << "Error: cannot find the unified non-linearity curve!!!" << endl;
      exit(0);
    }

    cout << "Reading unified non-linearity model from " << unified_nl_filename
          << endl;

    g_unified_positron_nl =
        (TGraph*)unified_nl_file->Get(nominal_graph_name.Data())->Clone();
    for (Int_t i = 0; i < m_num_unified_nl_pars; i++) {
      g_unified_positron_nl_pulls[i] =
          (TGraph*)unified_nl_file->Get(pull_graph_name[i].Data())->Clone();
    }

    unified_nl_file->Close();

    // Copy into arrays to speed up
    for (Int_t ie = 0; ie < n_unified_nl_points; ie++) {
      Double_t e = 1.022 + 0.02 * ie;
      m_unified_positron_nl_e[ie] = e;
      m_unified_positron_nl_fac[ie] = g_unified_positron_nl->Eval(e);
      for (Int_t i = 0; i < m_num_unified_nl_pars; i++) {
        m_unified_positron_nl_err[i][ie] =
            g_unified_positron_nl_pulls[i]->Eval(e) -
            m_unified_positron_nl_fac[ie];
      }
    }


    nl_func = new TF1("nl_func", this, &Spectrum::nl_func_unified, 0, 12,
                      m_num_unified_nl_pars + 2, "Spectrum", "nl_func_unified");

    for (Int_t i = 0; i < m_num_unified_nl_pars; i++) {
      m_unified_nl_par_nominal[i] =
          data->getDouble(Form("defaultNonlinearPar%d", i));
    }

    for (Int_t i = 0; i < m_num_unified_nl_pars; i++) {
      nl_func->SetParameter(i, m_unified_nl_par[i]);
    }

    nl_func->SetParameter(m_num_unified_nl_pars, 1.0);
    nl_func->SetParameter(m_num_unified_nl_pars + 1, 0.0);


    TMatrixD covmatrix_unified(10, 10);

    // Read covariance matrix
    if (m_correlateUnifiedNonlinearPars) {
      double* mat_tmp;

      TFile* f_covmatrix =
          new TFile(data->getString("NonlinearCovmatrixFilename"));
      f_covmatrix->ls();
      mat_tmp =
          ((TMatrixD*)f_covmatrix->Get("coeffmatrix"))->GetMatrixArray();
      f_covmatrix->Close();

      for (Int_t i = 0; i < m_num_unified_nl_pars; i++) {
        for (Int_t j = 0; j < m_num_unified_nl_pars; j++) {
          m_unified_nl_par_covmatrix[i][j] =
              mat_tmp[i * m_num_unified_nl_pars + j];
        }
      }
    } else {
      for (Int_t i = 0; i < m_num_unified_nl_pars; i++) {
        for (Int_t j = 0; j < m_num_unified_nl_pars; j++) {
          if (i == j)
            m_unified_nl_par_covmatrix[i][j] = 1;
          else
            m_unified_nl_par_covmatrix[i][j] = 0;
        }
      }
    }

    covmatrix_unified.SetMatrixArray(&m_unified_nl_par_covmatrix[0][0]);
    covmatrix_unified.ResizeTo(m_num_unified_nl_pars, m_num_unified_nl_pars);
    cout << "Covariance matrix for UNIFIED non-liner parameters:" << endl;
    covmatrix_unified.Print();

    TDecompChol chol_unified(covmatrix_unified);
    chol_unified.Decompose();

    TMatrixD cmat_unified(chol_unified.GetU());
    //  cmat.Print();
    TMatrixD tcmat_unified(cmat_unified.Transpose(cmat_unified));
    // std::cout << "Cholesky matrix---------" << std::endl;
    // tcmat.Print();

    double* tmp_matrix_unified = tcmat_unified.GetMatrixArray();

    for (int i = 0; i < m_num_unified_nl_pars; i++) {
      for (int j = 0; j < m_num_unified_nl_pars; j++) {
        m_unified_nl_par_covmatrix_l[i][j] =
            tmp_matrix_unified[i * m_num_unified_nl_pars + j];
        //      cout << "\t" << m_unified_nl_par_covmatrix_l[i][j];
      }
      //    cout << endl;
    }

    // for (Int_t i = 0; i < m_num_unified_nl_pars; i++){
    //   m_unified_nl_par[i] = m_unified_nl_par_nominal[i];
    // }
  }


  this->updateAntinu();
}

void Spectrum::rebin(int nCombine)
{
  for (int istage = 0; istage < Nstage; istage++) {
    for (int idet = 0; idet < Ndetectors; idet++) {
      for (int idx = 0; idx < m_nSamples; idx++) {
        int newBinIdx = idx / nCombine;
        if (idx % nCombine == 0) {
          m_energy[idet][newBinIdx] = m_energy[idet][idx];
          m_antiNuSpectrum[istage][idet][newBinIdx] =
              m_antiNuSpectrum[istage][idet][idx] / nCombine;
          m_positronIavDistortedSpectrum[istage][idet][newBinIdx] =
              m_positronIavDistortedSpectrum[istage][idet][idx] / nCombine;
          m_positronNLSpectrum[istage][idet][newBinIdx] =
              m_positronNLSpectrum[istage][idet][idx] / nCombine;
          m_positronDetectedSpectrum[istage][idet][newBinIdx] =
              m_positronDetectedSpectrum[istage][idet][idx] / nCombine;
        } else {
          m_antiNuSpectrum[istage][idet][newBinIdx] +=
              m_antiNuSpectrum[istage][idet][idx] / nCombine;
          m_positronTrueSpectrum[istage][idet][newBinIdx] +=
              m_positronTrueSpectrum[istage][idet][idx] / nCombine;
          m_positronIavDistortedSpectrum[istage][idet][newBinIdx] +=
              m_positronIavDistortedSpectrum[istage][idet][idx] / nCombine;
          m_positronNLSpectrum[istage][idet][newBinIdx] +=
              m_positronNLSpectrum[istage][idet][idx] / nCombine;
          m_positronDetectedSpectrum[istage][idet][newBinIdx] +=
              m_positronDetectedSpectrum[istage][idet][idx] / nCombine;
        }
      }
    }
  }
  m_nSamples /= nCombine;
  m_binWidth *= nCombine;
}

// TH1F* Spectrum::monteCarlo()
// {

//   // FIXME: need to make histograms for all detectors
//   TH1F* mcResult = new TH1F("mcResult","",m_nSamples,m_eMin,m_eMax);
//   for(int idet=0; idet<Ndetectors; idet++){
//     for(int idx=0; idx<m_nSamples; idx++){
//       if(m_positronDetectedSpectrum[idx]>0){
//         mcResult->Fill(m_energy[idet][idx]+m_binWidth/2.,
//                        gRandom->Poisson(m_positronDetectedSpectrum[idx]
//                                         * m_binWidth));
//       }
//     }
//   }
//   return mcResult;
// }


void Spectrum::passPredictor(Predictor* pred_in)
{
  pred = pred_in;
  this->extractPredictorData();
}

void Spectrum::extractPredictorData()
{
  for (int istage = 0; istage < Nstage; istage++) {
    for (int idet = 0; idet < Ndetectors; idet++) {
      m_runningTime[istage][idet] = pred->tdper[istage].Livetime[idet] * 1. /
                                    365; // convert from days to years
      m_detectorSize[idet] = pred->tdper[istage].TargetMass[idet] *
                             1e-6; // convert from kg to ktons

      m_detectorEfficiency_Ed_nominal[istage][idet] =
        pred->tdper[istage].DelayedEff[idet];
      m_detectorEfficiency_Ed[istage][idet] =
        m_detectorEfficiency_Ed_nominal[istage][idet];

      m_detectorEfficiency[istage][idet] =
          pred->tdper[istage].DMCEff[idet] *
          pred->tdper[istage].MuonVetoEff[idet];
      m_detectorEfficiency[istage][idet] *=
          m_detectorEfficiency_Dt * m_detectorEfficiency_Ep *
          m_detectorEfficiency_Ed[istage][idet] * m_detectorEfficiency_flash *
          m_detectorEfficiency_nGd * m_detectorEfficiency_spill;
    }
  }
}


void Spectrum::setRandomSeed(unsigned int i)
{
  if (UseChristineModel || m_useAbInitioSpectra) {
    for (int istage = 0; istage < Nstage; istage++)
      m_corespectrum[istage]->setRandomSeed(i + 100000);
  } /*else{
      for (int icore = 0; icore < Ncores; icore++){
      m_isotopes[icore]->setRandomSeed(i+100000); // set same random seeds for
      all cores since variation of the fission fraction are correlated between
      ADs for the first order
      }
      }*/
  ran->SetSeed(i);
}


void Spectrum::setRandomReactorPower()
{
  for (int istage = 0; istage < Nstage; istage++) {
    for (int i = 0; i < Ncores; i++) {
      m_reactorPower[istage][i] =
          (1 + m_reactorPowerError[istage][i] * ran->Gaus(0, 1)) *
          m_nominalReactorPower[istage][i];
    }
  }
}


void Spectrum::setRandomIavDistortion()
{
  double factor[Ndetectors];
  for (int idet = 0; idet < Ndetectors; idet++) {
    factor[idet] = 1. + m_iav_error * ran->Gaus(0, 1);
  }

  for (int idet = 0; idet < Ndetectors; idet++) {
    for (int idx = 0; idx < m_nSamples; idx++) {
      m_iav_frac[idet][idx] = factor[idet] * m_nominal_iav_frac[idx];
      m_iav_corr[idet][idx][idx] = 1 - m_iav_frac[idet][idx];
      for (int jdx = 0; jdx < m_nSamples; jdx++) {
        if (idx != jdx) {
          m_iav_corr[idet][idx][jdx] =
              factor[idet] * m_nominal_iav_corr[idx][jdx];
        }
      }
    }

    // // for debugging. Checking the sum....
    // for(int idx=0; idx<m_nSamples; idx++){
    //   double test = 0;
    //   for(int jdx=0; jdx<m_nSamples; jdx++){
    //     test += m_iav_corr[idet][idx][jdx];
    //   }
    //   cout << "---- " << test << endl;
    // }
  }
}


void Spectrum::setRandomScintiNonLinear()
{ // randomize alpha
  m_alpha = m_alpha_nominal + m_alpha_error * ran->Gaus(0, 1);
}

void Spectrum::setRandomIhepScintiNonLinear()
{ // randomize ihep nonlinear model parameters and load to the function

  while (true) {
    bool ParSetOK = true;
    Double_t ranvec_uniform[5];
    Double_t ranvec[5];

    for (Int_t i = 0; i < 5; i++) {
      ranvec_uniform[i] = ran->Gaus(0, 1);
    }
    for (Int_t i = 0; i < 5; i++) {
      ranvec[i] = 0;
      for (Int_t j = 0; j < 5; j++) {
        ranvec[i] += m_ihep_nl_par_covmatrix_l[i][j] * ranvec_uniform[j];
      }
    }
    for (Int_t i = 0; i < 5; i++) {
      m_ihep_nl_par[i] =
          m_ihep_nl_par_nominal[i] + m_ihep_nl_par_error[i] * ranvec[i];
    }
    // m_ihep_nl_par[2]
    //   = m_ihep_nl_par_nominal[2]
    //   + (m_ihep_nl_par[1] - m_ihep_nl_par_nominal[1])
    //   * m_ihep_nl_par_error[2]/m_ihep_nl_par_error[1]; // assume p1 and p2
    //   are 100% correlated

    if (m_ihep_nl_par[1] < 0 || m_ihep_nl_par[2] < 0) {
      ParSetOK = false;
      cout << "Found a negative parameter for IHEP non-linearlity function. "
              "Regenerating a random parameter set...."
           << endl;
    }
    if (ParSetOK)
      break;
  }


  // cout << "IHEP non-linaerlity parameters: ";
  // for (Int_t i = 0; i < 5; i++)
  //   cout << "\t" << m_ihep_nl_par[i] ;
  // cout << endl;
}

void Spectrum::setRandomBcwScintiNonLinear()
{ // randomize ihep nonlinear model parameters and load to the function


  while (true) {
    bool ParSetOK = true;
    // Double_t ranvec_uniform[5];
    // Double_t ranvec[5];

    // for (Int_t i = 0; i < 5; i++){
    //   ranvec_uniform[i] = ran->Gaus(0,1);
    // }
    // for (Int_t i = 0; i < 5; i++){
    //   ranvec[i] = 0;
    //   for (Int_t j = 0; j < 5; j++){
    //     ranvec[i] += m_bcw_nl_par_covmatrix_l[i][j]*ranvec_uniform[j];
    //   }
    // }
    // for (Int_t i = 0; i < 5; i++){
    //   m_bcw_nl_par[i] = m_bcw_nl_par_nominal[i] + m_bcw_nl_par_error[i] *
    //   ranvec[i];

    // }

    // // if (m_bcw_nl_par[1] < 0 || m_bcw_nl_par[2] < 0){
    // //   ParSetOK = false;
    // //   cout << "Found a negative parameter for BCW non-linearlity function.
    // Regenerating a random parameter set...." << endl;
    // // }

    for (Int_t i = 0; i < 5; i++) {
      m_bcw_elec_nl_par[i] = m_bcw_elec_nl_par_nominal[i] +
                             m_bcw_elec_nl_par_error[i] * ran->Gaus(0, 1);
    }
    if (ParSetOK)
      break;
  }
}

void Spectrum::setRandomLbnlScintiNonLinear()
{ // randomize LBNL nonlinear model parameters and load to the function
  for (Int_t i = 0; i < 3; i++) {
    m_lbnl_nl_par[i] =
        m_lbnl_nl_par_nominal[i] + m_lbnl_nl_par_error[i] * ran->Gaus(0, 1);
  }
}

void Spectrum::setRandomUnifiedScintiNonLinear()
{
  Double_t ranvec_uniform[m_num_unified_nl_pars];
  Double_t ranvec[m_num_unified_nl_pars];

  for (Int_t i = 0; i < m_num_unified_nl_pars; i++) {
    ranvec_uniform[i] = ran->Gaus(0, 1);
  }
  for (Int_t i = 0; i < m_num_unified_nl_pars; i++) {
    ranvec[i] = 0;
    for (Int_t j = 0; j < m_num_unified_nl_pars; j++) {
      ranvec[i] += m_unified_nl_par_covmatrix_l[i][j] * ranvec_uniform[j];
    }
  }
  for (Int_t i = 0; i < m_num_unified_nl_pars; i++) {
    m_unified_nl_par[i] =
        m_unified_nl_par_nominal[i] + m_unified_nl_par_error[i] * ranvec[i];
  }

  // for (Int_t i = 0; i < m_num_unified_nl_pars; i++){
  //   m_unified_nl_par[i] = m_unified_nl_par_nominal[i] +
  //   m_unified_nl_par_error[i] * ran->Gaus(0,1);
  // }
}


void Spectrum::setRandomElecNonLinear()
{
  // currently, assumes that it is covered by setRandomScintiNonLinear().
}

void Spectrum::setRandomResolution()
{
  // m_detectorResolution
  //   = m_detectorResolution_nominal + m_detectorResolution_error *
  //   ran->Gaus(0,1);
  double corr_bias = m_detectorResolution_error * ran->Gaus(0, 1);
  for (int idet = 0; idet < Ndetectors; idet++) {
    m_detectorResolution_bias[idet] =
        corr_bias + m_detectorResolution_error_uncorr * ran->Gaus(0, 1);
  }
}

void Spectrum::setRandomDetectorEfficiency()
{
  for (int idet = 0; idet < Ndetectors; idet++) {
    // Ensures the same detector gets the same efficiency over different periods
    Double_t ran_det_eff =
        (1 + m_detectorEfficiency_rel_error * ran->Gaus(0, 1));

    for (int istage = 0; istage < Nstage; istage++) {
      m_detectorEfficiency[istage][idet] =
          ran_det_eff * pred->tdper[istage].DMCEff[idet] *
          pred->tdper[istage].MuonVetoEff[idet] * m_detectorEfficiency_Dt *
          m_detectorEfficiency_Ep * m_detectorEfficiency_Ed[istage][idet] *
          m_detectorEfficiency_flash * m_detectorEfficiency_nGd *
          m_detectorEfficiency_spill;
    }
  }
}

void Spectrum::setRandomRelativeEnergyScale()
{
    for (int idet = 0; idet < Ndetectors; idet++) {
      double rel_escale_shift = m_rel_escale_error[idet] * ran->Gaus(0, 1);
      m_rel_escale[idet] = m_rel_escale_nominal[idet] + rel_escale_shift;

      for (int istage = 0; istage < Nstage; istage++) {
        // recalculate detection efficiency
        m_detectorEfficiency_Ed[istage][idet] =
            m_detectorEfficiency_Ed_nominal[istage][idet] *
            (1.0 +
            rel_escale_shift *
                0.36); // 0.36 is a magic factor that convert energy scale shift to
                        // delayed energy cut efficiency (0.072/0.2) in DocDB-10956

          m_detectorEfficiency[istage][idet] =
              pred->tdper[istage].DMCEff[idet] *
              pred->tdper[istage].MuonVetoEff[idet] * m_detectorEfficiency_Dt *
              m_detectorEfficiency_Ep * m_detectorEfficiency_Ed[istage][idet] *
              m_detectorEfficiency_flash * m_detectorEfficiency_nGd *
              m_detectorEfficiency_spill;
    }
  }
}

void Spectrum::setRandomAbsoluteEnergyScale()
{
  m_abs_escale = m_abs_escale_nominal + m_abs_escale_error * ran->Gaus(0, 1);
}


void Spectrum::setRandomRelativeEnergyOffset()
{
  for (int istage = 0; istage < Nstage; istage++) {
    for (int idet = 0; idet < Ndetectors; idet++) {
      m_rel_eoffset[istage][idet] = m_rel_eoffset_error * ran->Gaus(0, 1);
    }
  }
}

void Spectrum::setRandomAbsoluteEnergyOffset()
{
  m_abs_eoffset = m_abs_eoffset_error * ran->Gaus(0, 1);
}


// For now using a straight line to distort spectra; basically fluctuating slope
TF1* Spectrum::getDistortionCurve(double amount)
{
  TF1* func = new TF1("func", "TMath::Abs([0]+[1]*x)", m_eMin, m_eMax);
  // ran->SetSeed(0);
  double slope = amount * ran->Gaus(0, 1);
  double anchor_point = 3.5;
  // want offset to be set by requiring func @ anchor point to be 1
  double offset = (1 - slope * anchor_point);
  func->SetParameter(0, offset);
  func->SetParameter(1, slope);
  return func;
}

// Use a different function for fast-neutrons
TF1* Spectrum::getFNDistortionCurve(double amount)
{
  TF1* func = new TF1("func", "[0]/(0.2*pow(x,0.1))+[1]", m_eMin, m_eMax);
  // ran->SetSeed(0);
  double par = amount * ran->Gaus(0, 1);
  func->SetParameter(0, par);
  func->SetParameter(1, 0);
  // set offset so that ira func(m_eMax)=1;
  func->SetParameter(1, 1 - 1 * func->Eval(m_eMax));

  return func;
}

Double_t nl_func_ihep(Double_t* x, Double_t* par)
{
  Double_t e_true = x[0] - 1.022;  // subtract annihiration gamma energy
  Double_t escale_par = par[5];    // Add flat energy scale parameter
  Double_t escale_offset = par[6]; // Add fixed energy offset

  Double_t nl_fac_511 = par[4]; // magic factor read from Soren's slides


  // if e_true  < 0, assume e_true = 0 and calculate non-linearlity factor
  // Then apply to the true positron + gamma energy
  // this is needed since true energy sometimes become below 0 after Iav
  // correction.

  if (e_true < 0) {
    e_true = 0;
  }

  Double_t e_nonlinear_fac =
      ((par[0] + par[3] * e_true) / (1 + par[1] * exp(-par[2] * e_true)) *
           e_true +
       1.022 * nl_fac_511) /
      (e_true + 1.022);

  Double_t e_nonlinear = x[0] * e_nonlinear_fac * escale_par + escale_offset;

  return e_nonlinear;
}


Double_t reso_func_bcw(Double_t* x, Double_t* par)
{
  Double_t e_orig = x[0];
  Double_t e_sigma = 1.0;

  if (e_orig > 0) {
    e_sigma = TMath::Sqrt(par[0] * par[0] + par[1] * par[1] / e_orig +
                          par[2] * par[2] / e_orig / e_orig);
  }

  return e_sigma;
}


// Double_t nl_func_bcw(Double_t * x, Double_t * par){
//   Double_t edep = x[0];
//   Double_t escale_par = par[5]; // Add flat energy scale parameter
//   Double_t escale_offset = par[6]; // Add fixed energy offset

//   Double_t e_nonlinear = edep *
//   NonLinearity(edep,par[0],par[1],par[2],par[3],par[4]) * escale_par +
//   escale_offset;

//   return e_nonlinear;

// }


// OLD BCW non-linearityr function
double NonLinearity(double depositE, double delta_a1, double delta_a2,
                    double delta_a3, double delta_a4, double delta_a5)
{ // non-linearity function by Xin

  Double_t energy = 0;
  if (depositE < 1.022) {
    depositE = 1.022;
  }

  if (depositE > 0) {
    energy = 1. / depositE;
  }

  Double_t a1 = 1.10476 * (1 + delta_a1);
  Double_t a2 = -0.303366 * (1 + delta_a2);
  Double_t a3 = 0.118501 * (1 + delta_a3);
  Double_t a4 = 1.01746 * (1 + delta_a4);
  Double_t a5 = 0.158503 * (1 + delta_a5);
  Double_t result;

  Double_t a6 =
      (a1 + a2 * 0.5 + a3 * pow(0.5, 2) - a4 - a5 * 0.5) / pow(0.5, 2);

  if (energy >= 0.5) {
    result = a1 + a2 * energy + a3 * pow(energy, 2);
  } else {
    result = a4 + a5 * energy + a6 * pow(energy, 2);
  }

  return result;
}


// NEW BCW non-linearityr function. Feb. 20, 2013
// NEW updated BCW non-linearityr function. Mar. 26, 2013
Double_t Spectrum::nl_func_bcw(Double_t* x, Double_t* par)
{
  double e_positron_true = x[0];
  double escale_par = par[5];    // Add flat energy scale parameter
  double escale_offset = par[6]; // Add fixed energy offset

  double scinti_nl_fac = 1;
  if (e_positron_true < m_bcw_positron_nl_e[0])
    scinti_nl_fac = m_bcw_positron_nl_fac[0];
  else if (e_positron_true > m_bcw_positron_nl_e[n_bcw_positron_nl - 1])
    scinti_nl_fac = m_bcw_positron_nl_fac[n_bcw_positron_nl - 1];
  else {
    for (Int_t i = 0; i < n_bcw_positron_nl - 1; i++) {
      if (e_positron_true >= m_bcw_positron_nl_e[i] &&
          e_positron_true < m_bcw_positron_nl_e[i + 1]) {
        scinti_nl_fac = ((m_bcw_positron_nl_e[i + 1] - e_positron_true) *
                             m_bcw_positron_nl_fac[i] +
                         (e_positron_true - m_bcw_positron_nl_e[i]) *
                             m_bcw_positron_nl_fac[i + 1]) /
                        (m_bcw_positron_nl_e[i + 1] - m_bcw_positron_nl_e[i]);
        break;
      }
    }
  }

  double visibleE = scinti_nl_fac * e_positron_true;
  //  double electronicsCorrection = exp(par[0] + par[1] * visibleE) +
  //  exp(par[2] + par[3] * visibleE);

  double err_offset = 0;
  double err_shift = 0;

  double elec_err_min_x = 3.45; // MeV
  double par3_up = g_bcw_elec_nl_error[0]->Eval(elec_err_min_x) - 1;
  double par3_down = g_bcw_elec_nl_error[1]->Eval(elec_err_min_x) - 1;
  if (par[3] > 0)
    err_offset = par[3] * par3_up;
  else
    err_offset = par[3] * par3_down;

  double par4_up = 0;
  double par4_down = 0;
  if (visibleE > elec_err_min_x) {
    par4_up = g_bcw_elec_nl_error[0]->Eval(visibleE) -
              g_bcw_elec_nl_error[0]->Eval(elec_err_min_x);
    par4_down = g_bcw_elec_nl_error[1]->Eval(visibleE) -
                g_bcw_elec_nl_error[1]->Eval(elec_err_min_x);
  } else {
    par4_up = g_bcw_elec_nl_error[1]->Eval(visibleE) -
              g_bcw_elec_nl_error[1]->Eval(elec_err_min_x);
    par4_down = g_bcw_elec_nl_error[0]->Eval(visibleE) -
                g_bcw_elec_nl_error[0]->Eval(elec_err_min_x);
  }
  if (par[4] > 0)
    err_shift = par[4] * par4_up;
  else
    err_shift = par[4] * par4_down;

  double electronicsCorrection =
      exp(par[0] + par[1] * visibleE) + par[2] + err_offset + err_shift;

  double final_energy =
      visibleE * electronicsCorrection * escale_par + escale_offset;

  return final_energy;
}


// LBNL non-linearityr function, based on beta-gamma spectra
Double_t nl_func_lbnl(Double_t* x, Double_t* par)
{
  // par[0]: flat energy scale shift for electron
  // par[1]: size of energy scale shift propotional to exp(-1.5*eVis) for
  // electron par[2]: size of energy scale shift due to Ge68 calibration point

  double e_positron_true = x[0];
  double escale_par = par[3];    // Add flat energy scale parameter
  double escale_offset = par[4]; // Add fixed energy offset


  double scinti_nl_fac = 1;
  double err[3];

  if (e_positron_true < m_lbnl_positron_nl_e[0]) {
    scinti_nl_fac = m_lbnl_positron_nl_fac[0];
    for (Int_t ierr = 0; ierr < 3; ierr++) {
      err[ierr] = m_lbnl_positron_nl_err[ierr][0];
    }
  } else if (e_positron_true > m_lbnl_positron_nl_e[299]) {
    scinti_nl_fac = m_lbnl_positron_nl_fac[299];
    for (Int_t ierr = 0; ierr < 3; ierr++) {
      err[ierr] = m_lbnl_positron_nl_err[ierr][299];
    }
  } else {
    for (Int_t i = 0; i < 299; i++) {
      if (e_positron_true >= m_lbnl_positron_nl_e[i] &&
          e_positron_true < m_lbnl_positron_nl_e[i + 1]) {
        scinti_nl_fac = ((m_lbnl_positron_nl_e[i + 1] - e_positron_true) *
                             m_lbnl_positron_nl_fac[i] +
                         (e_positron_true - m_lbnl_positron_nl_e[i]) *
                             m_lbnl_positron_nl_fac[i + 1]) /
                        (m_lbnl_positron_nl_e[i + 1] - m_lbnl_positron_nl_e[i]);
        for (Int_t ierr = 0; ierr < 3; ierr++) {
          err[ierr] = ((m_lbnl_positron_nl_e[i + 1] - e_positron_true) *
                           m_lbnl_positron_nl_err[ierr][i] +
                       (e_positron_true - m_lbnl_positron_nl_e[i]) *
                           m_lbnl_positron_nl_err[ierr][i + 1]) /
                      (m_lbnl_positron_nl_e[i + 1] - m_lbnl_positron_nl_e[i]);
        }

        break;
      }
    }
  }
  double random_nl_fac = scinti_nl_fac;


  for (Int_t ierr = 0; ierr < 3; ierr++) {
    random_nl_fac += par[ierr] * err[ierr];
  }
  double visibleE = random_nl_fac * e_positron_true;
  //  cout << "NL " << e_positron_true << " " << visibleE/e_positron_true <<
  //  endl;

  return visibleE * escale_par + escale_offset;
}

// Unified non-linearityr function
Double_t Spectrum::nl_func_unified(Double_t* x, Double_t* par)
{
  double e_positron_true = x[0];
  double escale_par =
      par[m_num_unified_nl_pars]; // Add flat energy scale parameter
  double escale_offset =
      par[m_num_unified_nl_pars + 1]; // Add fixed energy offset

  double scinti_nl_fac = 1;
  double err[10];

  if (e_positron_true < m_unified_positron_nl_e[0]) {
    scinti_nl_fac = m_unified_positron_nl_fac[0];
    for (Int_t ierr = 0; ierr < m_num_unified_nl_pars; ierr++) {
      err[ierr] = m_unified_positron_nl_err[ierr][0];
    }
  } else if (e_positron_true >
             m_unified_positron_nl_e[n_unified_nl_points - 1]) {
    scinti_nl_fac = m_unified_positron_nl_fac[n_unified_nl_points - 1];
    for (Int_t ierr = 0; ierr < m_num_unified_nl_pars; ierr++) {
      err[ierr] = m_unified_positron_nl_err[ierr][n_unified_nl_points - 1];
    }
  } else {
    for (Int_t i = 0; i < n_unified_nl_points - 1; i++) {
      if (e_positron_true >= m_unified_positron_nl_e[i] &&
          e_positron_true < m_unified_positron_nl_e[i + 1]) {
        scinti_nl_fac =
            ((m_unified_positron_nl_e[i + 1] - e_positron_true) *
                 m_unified_positron_nl_fac[i] +
             (e_positron_true - m_unified_positron_nl_e[i]) *
                 m_unified_positron_nl_fac[i + 1]) /
            (m_unified_positron_nl_e[i + 1] - m_unified_positron_nl_e[i]);
        for (Int_t ierr = 0; ierr < m_num_unified_nl_pars; ierr++) {
          err[ierr] =
              ((m_unified_positron_nl_e[i + 1] - e_positron_true) *
                   m_unified_positron_nl_err[ierr][i] +
               (e_positron_true - m_unified_positron_nl_e[i]) *
                   m_unified_positron_nl_err[ierr][i + 1]) /
              (m_unified_positron_nl_e[i + 1] - m_unified_positron_nl_e[i]);
        }

        break;
      }
    }
  }
  double random_nl_fac = scinti_nl_fac;

  for (Int_t ierr = 0; ierr < m_num_unified_nl_pars; ierr++) {
    //    cout << "\t" << par[ierr];
    random_nl_fac += par[ierr] * err[ierr];
  }
  //  cout << endl;


  // double random_nl_fac = g_unified_positron_nl->Eval(e_positron_true);

  // for (Int_t i = 0; i < m_num_unified_nl_pars; i++){
  //   random_nl_fac += m_unified_nl_par[i] *
  //   g_unified_positron_nl_pulls[i]->Eval(e_positron_true);
  // }


  double visibleE = random_nl_fac * e_positron_true;

  return visibleE * escale_par + escale_offset;
}


//////////////////////////////////

Double_t enu_to_epositron(Double_t* x, Double_t* par)
{
  Double_t Enu = x[0];
  Double_t costheta = par[0];

  Double_t Me = 0.510999; // MeV
  Double_t Mn = 939.565;  // MeV
  Double_t Mp = 938.272;  // MeV
  Double_t Delta = Mn - Mp;

  Double_t Ee0 = Enu - Delta;
  // if (Ee0 < Me){    // not allowed!
  //   return -1;
  // }

  Double_t gamma0 = Ee0 / Me;
  Double_t beta0 = sqrt(1 - 1 / gamma0 / gamma0);
  Double_t y2 = (Delta * Delta - Me * Me) / 2.;

  if (costheta < -1 || costheta > 1) { // then use average cos theta
    costheta = -0.034 * beta0 + 2.4 * Enu / Mp;
  }

  Double_t Ee1 = Ee0 * (1 - Enu / Mp * (1.0 - costheta * beta0)) - y2 / Mp;

  //   if (Ee1 < Me){    // not allowed!
  //     return -1;
  //   }

  return Ee0 + Me; // Add Me to include annihiration gamma energy
}

//////////////////////////////////


Double_t enu_to_epositron_0th(Double_t* x, Double_t* par)
{
  Double_t Enu = x[0];

  Double_t Me = 0.510999; // MeV
  Double_t Mn = 939.565;  // MeV
  Double_t Mp = 938.272;  // MeV
  Double_t Delta = Mn - Mp;


  Double_t Ee0 = Enu - Delta;
  // if (Ee0 < Me){    // not allowed!
  //   return -1;
  // }

  return Ee0 + Me; // Add Me to include annihiration gamma energy
}
