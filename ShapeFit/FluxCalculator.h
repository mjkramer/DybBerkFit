#pragma once

#include "Config.h"
#include "OscCalc.h"
#include "Utils.h"

#include "TFile.h"
#include "TH1.h"
#include "TObject.h"

#include <iostream>
#include <map>

using namespace Config;

// This class calculates flux matrix elements and extrapolation factors
// jpochoa, May 2012

class FluxCalculator : public TObject {
public:
  // Parameters for OscFuncQuick function
  static constexpr int NpointsS = 501;     // sinSq2t13
  static constexpr double Smin = 0.0;
  static constexpr double Smax = 0.2;
  static constexpr int NpointsD = 41;      // dm2ee
  // NB: if NpointsD <= 1 it will automatically use the default value
  // in osccalc; if NpointsD is not zero but query OscSpecIntQuick with
  // dm2=-1 it will find the closest point to the default value in osccalc
  static constexpr double Dmin = 1.5e-3;
  static constexpr double Dmax = 3.5e-3;

public:
  FluxCalculator();
  ~FluxCalculator();
  FluxCalculator(const Char_t* distancesfile, const Char_t* weeklyfluxdataname,
                 int nweeks_in);
  FluxCalculator(
      const Char_t* distancesfile,
      const Char_t* superhistname); //<--this one reads superhists and thus can
                                    // only work in inclusive mode
  double GetDistance(int idet, int icore);
  void LoadDistances(const Char_t* distancesmatrixname);
  void LoadWeeklyFlux(const Char_t* weeklyfluxdataname, int nweeks);

  std::map<int, TH1F*> CalculateFluxHistRow(int idet, double s22t13, int iweek,
                                            double dm2_ee = -1,
                                            double s22t14 = 0,
                                            double dm2_41 = -1, int term = -1);
  std::map<int, TH1F*> ExtrapolationFactorRow(int idet1, int idet2,
                                              double s22t13, int iweek,
                                              double dm2_ee = -1,
                                              double s22t14 = 0,
                                              double dm2_41 = -1);

  TH1F* h_flux[Nstage][Ncores]; // no longer used
  TH1F* h_super[Nstage][Ndetectors][Ncores];
  OscCalc* osccalc;

  // Super-matrix
  //  double *SuperMatrixRow(int idet, double s22t13, int iweek, double dm2=-1);
  // void LoadSuperMatrix(Char_t *fluxmatrixname);  //<--- now discontinued
  // bool useSuperMatrix; //<-- now discontinued
  void LoadSuperHistograms(const Char_t* superhistsname);
  bool useSuperHists;
  double PercContr[Ndetectors][Ncores];


  //-->speeding up the code (to go to slow but precise replace OscSpecIntQuick
  // by OscSpecInt in .cc file)
  TH1F* OscSpecQuick(int idet, int icore, double s22t13, int iweek, double dm2);
  void OscSpecQuick(int idet, int icore, double s22t13, int iweek, double dm2,
                    TH1F* hout);

  int nweeks;

private:
  double LookupOscSpecInt[Ndetectors][Ncores][Nstage][NpointsS][NpointsD];
  double LookupS[NpointsS];
  double LookupD[NpointsD];
  TH1F* LookupOscSpec[Ndetectors][Ncores][Nstage][NpointsS][NpointsD];

  Double_t Distance[Ndetectors][Ncores]; //<---indices are ADs,Cores (D1, D2,
                                         // L1, L2, L3, L4)
  TFile* WeeklyFluxData;
  TFile* SuperHistData;
  Char_t filename[1024];
  TH1F* FluxSpec;
  // std::map<int, TH1F*> fluxmap;


  Char_t dummyname[1024];

  std::map<int, TH1F*> mapout;
  TH1F* htotal;
  TH1F* hout[Ndetectors];
  TH1F* hout_near;

  Bool_t FirstTime;

public:
  ClassDef(FluxCalculator, 1);
};
