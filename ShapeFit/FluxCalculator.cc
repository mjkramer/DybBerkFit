#include "FluxCalculator.h"
#include <fstream>
#include <iostream>

using namespace std;

FluxCalculator::FluxCalculator()
{
  osccalc = new OscCalc();
  useSuperHists = false;
};

// Mode with weekly periods
FluxCalculator::FluxCalculator(const Char_t* distancesfile,
                               const Char_t* weeklyfluxdataname, int nweeks_in)
{
  nweeks = nweeks_in;
  useSuperHists = false;

  cout << "----------- Creating flux calculator ---------------" << endl;
  LoadDistances(distancesfile);
  cout << " --> Distances loaded successfully" << endl;
  osccalc = new OscCalc();
  cout << " --> OscCalc created successfully" << endl;
  LoadWeeklyFlux(weeklyfluxdataname, nweeks);
  cout << " --> Loaded flux file successfully" << endl;

  FirstTime = true;
};

// Mode with super-histograms
FluxCalculator::FluxCalculator(const Char_t* distancesfile,
                               const Char_t* superhistname)
{
  nweeks = 1; //<--have to use inclusive mode if use super-histograms
  useSuperHists =
      true; //<-- this mode is by definition the one with super-histograms

  cout << "----------- Creating flux calculator ---------------" << endl;
  cout << "Note: using super-histogram mode of flux-calculator, so need "
          "inclusive mode"
       << endl;
  LoadDistances(distancesfile);
  cout << " --> Distances loaded successfully" << endl;
  osccalc = new OscCalc();
  cout << " --> OscCalc created successfully" << endl;
  LoadSuperHistograms(superhistname);
  cout << " --> Super-histograms loaded successfully" << endl;

  FirstTime = true;
}

FluxCalculator::~FluxCalculator()
{
  delete osccalc;
  delete WeeklyFluxData;
  delete FluxSpec;
  for (int istage = 0; istage < Nstage; ++istage) {
    for (int icore = 0; icore < Ncores; ++icore) {
      delete h_flux[istage][icore];
      for (int idetector = 0; idetector < Ndetectors; ++idetector)
        delete h_super[istage][idetector][icore];
    }
  }
};

void FluxCalculator::LoadDistances(const Char_t* distancematrixname)
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
    Distance[adnum - 1][0] = d1;
    Distance[adnum - 1][1] = d2;
    Distance[adnum - 1][2] = l1;
    Distance[adnum - 1][3] = l2;
    Distance[adnum - 1][4] = l3;
    Distance[adnum - 1][5] = l4;
  }
}

void FluxCalculator::LoadWeeklyFlux(const Char_t* weeklyfluxdataname,
                                    int nweeks_in)
{
#pragma omp single copyprivate(WeeklyFluxData)
  WeeklyFluxData = new TFile(weeklyfluxdataname, "READ");

  Char_t filenametemp[1024];
  for (int iweek = 0; iweek < nweeks_in; ++iweek) {
    for (int icore = 0; icore < Ncores; ++icore) {
      sprintf(filenametemp, "Week%i/%i", iweek, icore);
      h_flux[icore][iweek] = (TH1F*)WeeklyFluxData->Get(filenametemp);
    }
  }

  cout << "----------------- WARNING: using the same flux shape for all cores "
          "for testing -------------------------"
       << endl;

} // end of LoadWeeklyFlux

std::map<int, TH1F*>
FluxCalculator::CalculateFluxHistRow(int idet, double s22t13, int istage,
                                     double dm2, double s22t14, double dm2_41,
                                     int term)
{
  // Sanity check
  // if(useSuperHists && iweek>0){
  //  cout << "Ay caramba!!! You are using super-hist mode but looping over an
  //  iweek larger than 0!! Cannot be. I recommend you abort!" << endl;
  //}

  // Set parameters for osccalc
  if (dm2 != -1) {
    osccalc->SetDeltaM2_ee(dm2);
  } else {
    dm2 = osccalc->GetDeltaM2_ee();
  }
  if (dm2_41 != -1) {
    osccalc->SetDeltaM2_41(dm2_41);
  } else {
    dm2_41 = osccalc->GetDeltaM2_41();
  }
  osccalc->SetTheta14(TMath::ASin(sqrt(s22t14)) * 0.5);

  if (FirstTime) {
    for (int icore = 0; icore < Ncores; ++icore) {
      if (!useSuperHists) {
        hout[icore] =
            (TH1F*)h_flux[istage][icore]->Clone(Form("hout_%d", icore));
      } else {
        hout[icore] =
            (TH1F*)h_super[istage][idet][icore]->Clone(Form("hout_%d", icore));
      }
    }
    htotal = (TH1F*)hout[0]->Clone("htotal");
    hout_near = (TH1F*)hout[0]->Clone("hout_near");
  }
  FirstTime = false;

  htotal->Reset();
  for (int icore = 0; icore < Ncores; ++icore) {
    //    sprintf(dummyname,"FluxHist_%i_%i",idet,icore);
    //    hout[icore] =
    //    (TH1F*)OscSpecQuick(idet,icore,s22t13,iweek)->Clone("clone");
    //    OscSpecQuick(idet,icore,s22t13,iweek,dm2,hout[icore]); // not using
    //    for now, since it dosent speed up for one fit

    if (!useSuperHists) {
      osccalc->OscSpec(Distance[idet][icore], s22t13, h_flux[istage][icore],
                       hout[icore]);
    } else {
      osccalc->OscSpecBinInt(Distance[idet][icore], s22t13,
                             h_super[istage][idet][icore], hout[icore], term);
    }

    // SuperMatrix now obsolete and replaced by super-histograms, which are a
    // more correct treatment if (useSuperMatrix){
    // hout[icore]->Divide(h_flux[istage][icore]);
    // hout[icore]->Scale(PercContr[idet][icore]);
    //}else{

    hout[icore]->Scale(1. / pow(Distance[idet][icore], 2));
    htotal->Add(hout[icore], 1);
  }

  // Spit out reactor contribution
  for (int icore = 0; icore < Ncores; ++icore) {
    if (term < 0)
      hout[icore]->Divide(htotal);
    // if term is specified, do not divide by the total, so that one can compare
    // with oter terms and detecgtors
    mapout[icore] = hout[icore];
  }

  return mapout;

} // end of CalculateFluxRow

double FluxCalculator::GetDistance(int idet, int icore)
{
  if (idet >= Ndetectors) {
    cout << "WARNING!!! idet=" << idet
         << " larger than current number of detectors" << endl;
    return -1;
  }
  if (icore >= Ncores) {
    cout << "WARNING!!! icore=" << icore
         << " larger than current number of detectors" << endl;
    return -1;
  }

  return Distance[idet][icore];
}

std::map<int, TH1F*>
FluxCalculator::ExtrapolationFactorRow(int idet_far, int idet_near,
                                       double s22t13, int istage, double dm2,
                                       double s22t14, double dm2_41)
{
  // Set parameters for osccalc
  if (dm2 != -1) {
    osccalc->SetDeltaM2_ee(dm2);
  } else {
    dm2 = osccalc->GetDeltaM2_ee();
  }
  if (dm2_41 != -1) {
    osccalc->SetDeltaM2_41(dm2_41);
  } else {
    dm2_41 = osccalc->GetDeltaM2_41();
  }
  osccalc->SetTheta14(TMath::ASin(sqrt(s22t14)) * 0.5);

  if (FirstTime) {
    for (int icore = 0; icore < Ncores; ++icore) {
      if (!useSuperHists) {
        hout[icore] =
            (TH1F*)h_flux[istage][icore]->Clone(Form("hout_%d", icore));
      } else {
        hout[icore] =
            (TH1F*)h_super[istage][0][icore]->Clone(Form("hout_%d", icore));
      }
    }

    htotal = (TH1F*)hout[0]->Clone("htotal");
    hout_near = (TH1F*)hout[0]->Clone("hout_near");
  }
  FirstTime = false;

  for (int icore = 0; icore < Ncores; ++icore) {
    if (!useSuperHists) {
      osccalc->OscSpec(Distance[idet_far][icore], s22t13, h_flux[istage][icore],
                       hout[icore]);
      osccalc->OscSpec(Distance[idet_near][icore], s22t13,
                       h_flux[istage][icore], hout_near);
    } else {
      if (idet_near < 0) { // then, calculate oscillation probabilities at a
                           // particular detector
        osccalc->OscSpecBinInt(Distance[idet_far][icore], s22t13,
                               h_super[istage][idet_far][icore], hout[icore]);
        osccalc->OscSpecBinInt(Distance[idet_far][icore], s22t13,
                               h_super[istage][idet_far][icore], hout_near, 0);
      } else {
        osccalc->OscSpecBinInt(Distance[idet_far][icore], s22t13,
                               h_super[istage][idet_far][icore], hout[icore]);
        osccalc->OscSpecBinInt(Distance[idet_near][icore], s22t13,
                               h_super[istage][idet_near][icore], hout_near);
      }
      // osccalc->OscSpec(Distance[idet_far][icore],s22t13,
      // h_super[istage][idet_far][icore],hout[icore]);
      // osccalc->OscSpec(Distance[idet_near][icore],s22t13,
      // h_super[istage][idet_near][icore],hout_near);
    }

    hout[icore]->Divide(hout_near);


    if (idet_near < 0) {
      // then, do not need to scale, since the distances are the same
    } else {
      hout[icore]->Scale(pow(Distance[idet_near][icore], 2) * 1. /
                         pow(Distance[idet_far][icore], 2));
    }
    mapout[icore] = hout[icore];
  }
  return mapout;
}


TH1F* FluxCalculator::OscSpecQuick(int idet, int icore, double s22t13,
                                   int iweek, double dm2)
{
  //-->find best point for sin22t13
  double bestS = (s22t13 - Smin) * (NpointsS - 1) * 1. / (Smax - Smin);
  int bestSi = (int)floor(bestS + 0.5);
  //-->find best point for dm2
  int bestDi = 0;
  if (dm2 == -1) {
    dm2 = osccalc->GetDeltaM2_ee();
  }
  if (NpointsD > 1) {
    double bestD = (dm2 - Dmin) * (NpointsD - 1) * 1. / (Dmax - Dmin);
    bestDi = (int)floor(bestD + 0.5);
  }

  //-->evaluate function
  return LookupOscSpec[idet][icore][iweek][bestSi][bestDi];
}

void FluxCalculator::OscSpecQuick(int idet, int icore, double s22t13, int iweek,
                                  double dm2, TH1F* h_osc)
{
  TH1F* h = (TH1F*)OscSpecQuick(idet, icore, s22t13, iweek, dm2);
  Int_t nbins = h->GetNbinsX();
  for (Int_t ibin = 0; ibin < nbins; ibin++) {
    h_osc->SetBinContent(ibin + 1, h->GetBinContent(ibin + 1));
  }
}

// Note: super-histograms are created with GenerateSuperHistograms.C.
// Super-histograms are the neutrino flux (binned in true neutrino energy) seen
// in one AD in an entire period (.e.g. 20 weeks), using the weekly efficiencies
// and livetimes to weigh each element in the sum.
void FluxCalculator::LoadSuperHistograms(const Char_t* superhistsname)
{
  Char_t name1[1024];
  Char_t name2[1024];
#pragma omp single copyprivate(infile)
  SuperHistData = new TFile(superhistsname, "READ");
  for (int istage = 0; istage < Nstage; ++istage) {
    for (int idet = 0; idet < Ndetectors; ++idet) {
      for (int icore = 0; icore < Ncores; ++icore) {
        sprintf(name1, "h_super_istage%i_idet%i_icore%i", istage, idet, icore);
        sprintf(name2, "super_istage%i_idet%i_icore%i", istage, idet, icore);
        h_super[istage][idet][icore] =
            (TH1F*)SuperHistData->Get(name1)->Clone(name2);
      }
    }
  }
}
