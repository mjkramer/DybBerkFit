#define LBNL_FIT_STERILE

#include "Binning.h"
#include "Config.h"
#include "FluxCalculator.h"
#include "Paths.h"
#include "Predictor.h"

#include <TMinuit.h>

#include <fstream>

Predictor *pred;

int PeriodFlag = -1;//(0=6AD, 1=8AD, 2=7AD, -1=6+8+7AD)

void minuit_fcn(int &npar, double *gin, double &f, double *x, int iflag){ // function for minuit minimization
  double sin22t13 = x[0];
  double dm2_ee = x[1];
  double sin22t14 = x[2];
  double dm2_41 = x[3];
  f =  pred->CalculateChi2Cov(sin22t13,dm2_ee,sin22t14,dm2_41);
  //f = oscprobtab->CalculateChi2CovQuick(sin22t13,dm2_ee,sin22t14,dm2_41,PeriodFlag);
}

void DoMinuitFit(TMinuit *minu, double dm214, double s2tt14=-1,
                 bool fixedSterileParams=false);

void fit_shape_3d()
{
  pred = new Predictor;
  pred->SetStage(PeriodFlag);
  auto fluxcalc = new FluxCalculator(Paths::baselines(), Paths::histogram());
  pred->EnterFluxCalculator(fluxcalc);
  for (int istage = 0; istage < Nstage; ++istage)
    pred->LoadMainData(Paths::input(istage));
  pred->LoadPredictedIBD(Paths::predicted_ibd());
  pred->LoadIBDSpec(Paths::all_sig_spectra().data());
  pred->LoadBgSpec();
  pred->SetEvisBins(Binning::n_evis(), Binning::evis());
  pred->SetEnuBins(Binning::n_enu(), Binning::enu());
  pred->LoadEvisToEnuMatrix(Paths::response());
  pred->LoadCovMatrix(Paths::sig_covmatrix(), Paths::bg_covmatrix(),
                      Paths::dm2ee_covmatrix());

  // override values from Config.h
  const int nsteps = 1;
  const double s22t13start = S22T13;
  const double s22t13end = S22T13;

  const int nsteps_dm2 = 1;
  const double dm2eestart = DM2EE;
  const double dm2eeend = DM2EE;

  const int nS2T = nsteps_s22t14;
  const int nDM2 = nsteps_dm214_all;

  const Int_t nsteps_dm214_fit = 56;

  Double_t chi2result[nsteps_dm2][nsteps][nDM2][nS2T];
  Double_t dchi2result[nsteps_dm2][nsteps][nDM2][nS2T];
  Double_t sin22t13[nsteps_dm2][nsteps];
  Double_t dm213[nsteps_dm2][nsteps];

  Double_t chi2_min = 1e6;
  Double_t s2t_min = 1e6;
  Double_t dm2_min = 1e6;
  Double_t s2t14_min = 1e6;
  Double_t dm241_min = 1e6;

  Double_t s22t14_bins[nS2T+1];
  Double_t dm214_bins[nDM2+1];

  for (int step_dm214=0; step_dm214<nsteps_dm214; ++step_dm214){
    dm214_bins[step_dm214] = exp(log(dm214start) + log_dm214_step*((Double_t)step_dm214-0.5));
  }

  for (int step_dm214=0; step_dm214<nsteps_dm214_2+1; ++step_dm214){
    dm214_bins[step_dm214+nsteps_dm214] = dm214start_2 + lin_dm214_step*((Double_t)step_dm214-0.5);
  }

  for (int step_s22t14=0; step_s22t14<nsteps_s22t14+1; ++step_s22t14){
    s22t14_bins[step_s22t14] = exp(log(s22t14start) + log_s22t14_step*((Double_t)step_s22t14-0.5));
  }

  TDirectory * dir = gDirectory;

  TFile *savefile = new TFile(Paths::outpath("fit_shape_3d.root"), "RECREATE");
  TTree * tr = new TTree("tr_fit","fit results");
  tr->Branch("chi2_map",&dchi2result[0][0][0][0],Form("chi2_map[%d][%d][%d][%d]/D",nsteps_dm2,nsteps,nsteps_dm214,nsteps_s22t14));
  tr->Branch("chi2_min",&chi2_min,"chi2_min/D");
  tr->Branch("dm2_min",&dm2_min,"dm2_min/D");
  tr->Branch("s2t_min",&s2t_min,"s2t_min/D");
  // new branched for sterile mixing
  tr->Branch("dm241_min",&dm241_min,"dm241_min/D");
  tr->Branch("s2t14_min",&s2t14_min,"s2t14_min/D");

  dir->cd();

  for (Int_t iDM2 = 0; iDM2 < nDM2;  iDM2 ++) {
    double dm214 = dm214_bins[iDM2];
    for (Int_t iS2T = 0; iS2T < nS2T;  iS2T ++) {
      double s22t14 = s22t14_bins[iS2T];
    }
  }

  auto minu = new TMinuit(4);
  minu->SetPrintLevel(-1);
  minu->SetFCN(minuit_fcn);

  for(int ievt=0;ievt<1;++ievt){
    double minchi2=1e8;
    double bests22t13=0;
    double bestdm213=0;
    double bests22t14=0;
    double bestdm214=0;

    Double_t sin22t14;
    Double_t dm214;

    Double_t * grad;
    Double_t fpar,ferr;
    // semi parameter scan:
    for(int step_dm214=0;step_dm214<nsteps_dm214_fit;++step_dm214){
      //      dm214=exp(log(dm214end/dm214start)*step_dm214*1./(nsteps_dm214_fit-1)+log(dm214start));
      //      dm214=(dm214end-dm214start)*step_dm214*1./(nsteps_dm214-1)+dm214start;

      //      dm214=(0.501-0.001)*step_dm214*1./(nsteps_dm214_fit-1)+0.001;  // linear scan with 5e-3 steps
      //      dm214=exp(log(1.0/0.01)*step_dm214*1./(nsteps_dm214_fit-1)+log(0.01));
      // attempt to optimize the initial point distributions
      if (step_dm214 < 3){
        dm214= 0.003 * step_dm214+0.001;
      }else if (step_dm214 < 11){
        dm214= 0.005 * (step_dm214 - 3) +0.01;
      }else{
        dm214= 0.01 * (step_dm214 - 11) +0.05;
      }

      cout << " ========== " << step_dm214 << " / " << nsteps_dm214_fit <<" (initial dm2 = " << dm214 << " ) ========== " << endl;
      DoMinuitFit(minu, dm214, 0.02);

      Double_t pars[4];
      minu->GetParameter(0,fpar,ferr);
      pars[0] = fpar;
      minu->GetParameter(1,fpar,ferr);
      pars[1]= fpar;
      minu->GetParameter(2,fpar,ferr);
      pars[2] = fpar;
      minu->GetParameter(3,fpar,ferr);
      pars[3] = fpar;
      Double_t minchi2_tmp = 1e6;
      minu->Eval(4,grad,minchi2_tmp,pars,0);

      if (minchi2_tmp < minchi2){
        minchi2 = minchi2_tmp;
        bests22t13 = pars[0];
        bestdm213 =  pars[1];
        bests22t14 = pars[2];
        bestdm214 =  pars[3];
      }
    }

    Double_t best_pars[4] = {bests22t13,bestdm213,bests22t14,bestdm214};

    cout << "======== fit results (for checking) :"
         << " " <<  bests22t13 << " " << bestdm213
         << " " <<  bests22t14 << " " << bestdm214
         << " " << minchi2 << endl;


    for(int step_dm2=0;step_dm2<nsteps_dm2;++step_dm2){
      for(int step=0;step<nsteps;++step){
        for(int step_dm214=0;step_dm214<nsteps_dm214;++step_dm214){
          cout << " ========== " << step_dm214 << " / " << nsteps_dm214 <<" ========== " << endl;
          for(int step_s22t14=0;step_s22t14<nsteps_s22t14;++step_s22t14){
            if (nsteps == 1)
              sin22t13[step_dm2][step]=S22T13;
            else
              sin22t13[step_dm2][step]=(s22t13end-s22t13start)*step*1./(nsteps-1)+s22t13start;

            if (nsteps_dm2 == 1)
              dm213[step_dm2][step]=DM2EE;
            // else
            //   dm213[step_dm2][step]=(dm213end-dm213start)*step_dm2*1./(nsteps_dm2-1)+dm213start;

            sin22t14=(s22t14end-s22t14start)*step_s22t14*1./(nsteps_s22t14-1)+s22t14start;
            //            dm214=(dm214end-dm214start)*step_dm214*1./(nsteps_dm214-1)+dm214start;
            // log scale
            dm214=exp(log(dm214end/dm214start)*step_dm214*1./(nsteps_dm214-1)+log(dm214start));

            DoMinuitFit(minu, dm214, sin22t14, true);

            Double_t pars[4];
            minu->GetParameter(0,fpar,ferr);
            pars[0] = fpar;
            minu->GetParameter(1,fpar,ferr);
            pars[1]= fpar;
            minu->GetParameter(2,fpar,ferr);
            pars[2] = fpar;
            minu->GetParameter(3,fpar,ferr);
            pars[3] = fpar;
            Double_t minchi2_tmp = 1e6;
            minu->Eval(4,grad,minchi2_tmp,pars,0);

            chi2result[step_dm2][step][step_dm214][step_s22t14]=minchi2_tmp;

            //            chi2result[step_dm2][step][step_dm214][step_s22t14]=pred->CalculateChi2Cov(sin22t13[step_dm2][step],dm213[step_dm2][step]+7.5e-5,sin22t14,dm214);
            //            delete predset;
          }
        }
      }
    }//loop over t13

    cout << "Event: " << ievt << "; --> best fit: " << bests22t13 << "," << bestdm213 << "; minchi2=" << minchi2 << endl;//tmp

    for(int step_dm214=0;step_dm214<nsteps_dm214;++step_dm214){
      for(int step_s22t14=0;step_s22t14<nsteps_s22t14;++step_s22t14){
        dchi2result[0][0][step_dm214][step_s22t14] = chi2result[0][0][step_dm214][step_s22t14]-minchi2;
      }
    }

    chi2_min = minchi2;
    s2t_min = bests22t13;
    dm2_min = bestdm213;
    s2t14_min = bests22t14;
    dm241_min = bestdm214;

    tr->Fill();


  }//loop over toys

}

void DoMinuitFit(TMinuit* minu, double dm214, double s22t14,
                 bool fixedSterileParams)
{
  int ierflag;
  minu->mnparm(0, "SinSq2Theta13", S22T13, 0.01, 0, 0.2, ierflag);
  minu->mnparm(1, "DeltaMSqee", DM2EE, 0.0001, 0.0015, 0.0035, ierflag);
  minu->mnparm(2, "SinSq2Theta14",
               s22t14 == -1 ? 0.02 : s22t14,
               0.01, 0, 1.0, ierflag);
  minu->mnparm(3, "DeltaMSq41", dm214, 0.5*dm214, 1e-4, 10.0, ierflag);

  minu->FixParameter(1);

  if (fixedSterileParams) {
    minu->FixParameter(2);
    minu->FixParameter(3);
  } else {
    minu->Release(2);
    minu->Release(3);
  }

  double arglist[2];
  arglist[0] = 10000;
  arglist[1] = 1.0;
  minu->mnexcm("MIGRAD", arglist, 2, ierflag);
}