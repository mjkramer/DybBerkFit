#include "Binning.h"
#include "Config.h"
#include "FluxCalculator.h"
#include "OscProbTable.h"
#include "Paths.h"
#include "Predictor.h"

#include <TMinuit.h>

#include <cstdio>
#include <fstream>
#include <iostream>

static Predictor *pred;
#pragma omp threadprivate(pred)
static OscProbTable *oscprobtab;
#pragma omp threadprivate(oscprobtab)

int PeriodFlag = -1; //(0=6AD, 1=8AD, 2=7AD, -1=6+8+7AD)

void minuit_fcn(int &npar, double *gin, double &f, double *x,
                int iflag) { // function for minuit minimization
  double sin22t13 = x[0];
  double dm2_ee = x[1];
  double sin22t14 = x[2];
  double dm2_41 = x[3];
  // f =  pred->CalculateChi2Cov(sin22t13,dm2_ee,sin22t14,dm2_41,PeriodFlag);
  f = oscprobtab->CalculateChi2CovQuick(sin22t13, dm2_ee, sin22t14, dm2_41,
                                        PeriodFlag);
}

void DoMinuitFit(TMinuit *minu, double dm214, double s2tt14 = -1);

void fit_shape_3d_FC(const char *toyfilename, const char *outfilename, int ntoys = -1) {
  std::setbuf(stdout, nullptr);

#pragma omp parallel
  {
    pred = new Predictor;
    pred->SetStage(PeriodFlag);
    auto fluxcalc = new FluxCalculator(Paths::baselines(), Paths::histogram());
    pred->EnterFluxCalculator(fluxcalc);
    for (int istage = 0; istage < Nstage; ++istage)
      pred->LoadMainData(Paths::input(istage));
    pred->LoadPredictedIBD(Paths::predicted_ibd()); // dummy, not used
    // pred->LoadBgSpec();
    pred->SetEvisBins(Binning::n_evis(), Binning::evis());
    pred->SetEnuBins(Binning::n_enu(), Binning::enu());
    pred->LoadEvisToEnuMatrix(Paths::response());
    pred->LoadCovMatrix(Paths::sig_covmatrix(), Paths::bg_covmatrix(),
                        Paths::dm2ee_covmatrix());
    const int ntoys_all = pred->LoadToyIBDSpec(toyfilename);
#pragma omp master
    if (ntoys == -1)
      ntoys = ntoys_all;
    pred->LoadBgSpec();
  }

  // override values from Config.h
  const int nsteps = 1;
  const double s22t13start = S22T13;
  const double s22t13end = S22T13;

  const int nsteps_dm2 = 1;
  const double dm2eestart = DM2EE;
  const double dm2eeend = DM2EE;

  Ranger *ranger_dm41 = new Ranger();
  Ranger *ranger_dm41_2 = new Ranger();
  ranger_dm41->nsteps = nsteps_dm214;
  ranger_dm41->min = dm214start;
  ranger_dm41->max = dm214end;
  ranger_dm41->setLogScale();
  ranger_dm41_2->nsteps = nsteps_dm214_2;
  ranger_dm41_2->min = dm214start_2;
  ranger_dm41_2->max = dm214end_2;

  TFile *infile = new TFile(toyfilename, "READ");
  TTree *true_pars = (TTree *)infile->Get("true_pars");
  double trues22t13 = -1;
  double truedm2ee = -1;
  double trues22t14 = -1;
  double truedm241 = -1;
  true_pars->SetBranchAddress("sin22theta13", &trues22t13);
  true_pars->SetBranchAddress("deltamsqee", &truedm2ee);
  true_pars->SetBranchAddress("sin22theta14", &trues22t14);
  true_pars->SetBranchAddress("deltamsq41", &truedm241);
  true_pars->GetEntry(0);
  infile->Close();

#pragma omp parallel
  {
    oscprobtab = new OscProbTable(pred);
    oscprobtab->SetDMeeRange(nsteps_dm2, dm2eestart, dm2eeend);
    oscprobtab->SetDM41Range(nsteps_dm214, dm214start, dm214end, true);
    oscprobtab->SetDM41Range2(nsteps_dm214_2, dm214start_2, dm214end_2);
    oscprobtab->ReadTable(Paths::outpath("OscProbTable.txt"));
  }

  static TMinuit *minu;
#pragma omp threadprivate(minu)
#pragma omp parallel
  {
    minu = new TMinuit(4);
    minu->SetPrintLevel(-1);
    minu->SetFCN(minuit_fcn);
  }

  TFile *outfile = new TFile(outfilename, "recreate");
  TH1F *h_dchi2 = new TH1F("h_dchi2", "", 1150, -15, 100);

  //   static bool FixedTheMatrix = false;
  // #pragma omp threadprivate(FixedTheMatrix)

  // Since ntoys is (possibly) set by the master thread above, insert a barrier
#pragma omp barrier

#pragma omp parallel for
  for (int itoy = 0; itoy < ntoys; ++itoy) {
    pred->LoadToyMCEntry(itoy);

    pred->SetSin22t13Step(20, 0.00, 0.20); // Set here!
    pred->FixCovMatrix(S22T13, DM2EE, 0., 0.1e-3);

    double minchi2 = 1e8;

    // semi parameter scan:
    for (int step_dm214 = 0; step_dm214 < nsteps_dm214_all; ++step_dm214) {
      cout << " ========== " << step_dm214 << " / " << nsteps_dm214_all
           << " ========== " << endl;

      double dm214;
      if (step_dm214 < nsteps_dm214)
        dm214 = ranger_dm41->returnVal(step_dm214);
      else
        dm214 = ranger_dm41_2->returnVal(step_dm214 - nsteps_dm214);

      DoMinuitFit(minu, dm214);

      Double_t *grad;
      Double_t fpar, ferr;
      Double_t pars[4];
      Double_t err[4];

      minu->GetParameter(0, fpar, ferr);
      pars[0] = fpar;
      err[0] = ferr;
      minu->GetParameter(1, fpar, ferr);
      pars[1] = fpar;
      err[1] = ferr;
      minu->GetParameter(2, fpar, ferr);
      pars[2] = fpar;
      err[2] = ferr;
      minu->GetParameter(3, fpar, ferr);
      pars[3] = fpar;
      err[3] = ferr;
      Double_t minchi2_tmp = 1e6;
      minu->Eval(4, grad, minchi2_tmp, pars, 0);

      if (minchi2_tmp < minchi2) {
        minchi2 = minchi2_tmp;
      }
    }

    DoMinuitFit(minu, truedm241, trues22t14);
    double bests22t13_true14, bests22t13_true14_err;
    minu->GetParameter(0, bests22t13_true14, bests22t13_true14_err);
    // double
    // chi2_true=pred->CalculateChi2Cov(bests22t13_true14,truedm2ee,trues22t14,truedm241);
    double chi2_true = oscprobtab->CalculateChi2CovQuick(
        bests22t13_true14, truedm2ee, trues22t14, truedm241, PeriodFlag);

#pragma omp critical
    {
      h_dchi2->Fill(chi2_true - minchi2);
    }
  }

  outfile->cd();
  h_dchi2->Write();
  outfile->Close();
}

// If s22t14 is -1, float it; otherwise, fix it
void DoMinuitFit(TMinuit* minu, double dm214, double s22t14)
{
  int ierflag;
  minu->mnparm(0, "SinSq2Theta13", S22T13, 0.01, 0, 0.2, ierflag);
  minu->mnparm(1, "DeltaMSqee", DM2EE, 0.0001, 0.0015, 0.0035, ierflag);
  minu->mnparm(2, "SinSq2Theta14",
               s22t14 == -1 ? 0.02 : s22t14,
               0.01, 0, 1.0, ierflag);
  minu->mnparm(3, "DeltaMSq41", dm214, 0.5*dm214, 1e-4, 10.0, ierflag);

  minu->FixParameter(1);
  minu->FixParameter(3);

  if (s22t14 != -1)
    minu->FixParameter(2);
  else
    minu->Release(2);

  double arglist[2];
  arglist[0] = 10000;
  arglist[1] = 1.0;
  minu->mnexcm("MIGRAD", arglist, 2, ierflag);
}
