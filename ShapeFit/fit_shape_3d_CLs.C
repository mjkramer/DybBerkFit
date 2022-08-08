#define LBNL_FIT_STERILE

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

static int PeriodFlag = -1;//(0=6AD, 1=8AD, 2=7AD, -1=6+8+7AD)

static void minuit_fcn(int &npar, double *gin, double &f, double *x, int iflag){ // function for minuit minimization
  double sin22t13 = x[0];
  double dm2_ee = x[1];
  double sin22t14 = x[2];
  double dm2_41 = x[3];
  //f =  pred->CalculateChi2Cov(sin22t13,dm2_ee,sin22t14,dm2_41,PeriodFlag);
  f = oscprobtab->CalculateChi2CovQuick(sin22t13,dm2_ee,sin22t14,dm2_41,PeriodFlag);
}

static void DoMinuitFit(TMinuit *minu, double dm214, double s2tt14=-1);

void fit_shape_3d_CLs(bool fit4nuSamples=false, int igrid=-1)
{
  std::setbuf(stdout, nullptr);

  const int step_s22t14_start = igrid != -1 ? igrid : 0;
  const int step_s22t14_end = igrid != -1 ? (igrid+1) : nsteps_s22t14;

#pragma omp parallel
  {
    pred = new Predictor;
    pred->SetStage(PeriodFlag);
    auto fluxcalc = new FluxCalculator(Paths::baselines(), Paths::histogram());
    pred->EnterFluxCalculator(fluxcalc);
    for (int istage = 0; istage < Nstage; ++istage)
      pred->LoadMainData(Paths::input(istage));
    // NB: For the 4nu case, these are not the correct PredictedIBD values, but
    // we don't care since we're not looking at the summed matrix
    pred->LoadPredictedIBD(Paths::predicted_ibd());
    pred->LoadIBDSpec(Paths::all_sig_spectra().data()); // XXX
    pred->LoadBgSpec();
    pred->SetEvisBins(Binning::n_evis(), Binning::evis());
    pred->SetEnuBins(Binning::n_enu(), Binning::enu());
    pred->LoadEvisToEnuMatrix(Paths::response());
    pred->LoadCovMatrix(Paths::sig_covmatrix(), Paths::bg_covmatrix(),
                        Paths::dm2ee_covmatrix());

    if (fit4nuSamples) {
      auto path =
          Paths::outpath("toys_parscans/toySpectra_parscans_nominal.root");
      pred->LoadToyIBDSpec(path);
    } else {
      pred->LoadToyIBDSpec(Paths::toytree("sigsys"));
      pred->LoadToyMCNominalSpec();
    }
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

#pragma omp parallel
  {
    oscprobtab = new OscProbTable(pred);
    oscprobtab->SetDMeeRange(nsteps_dm2, dm2eestart, dm2eeend);
    oscprobtab->SetDM41Range(nsteps_dm214, dm214start, dm214end, true);
    oscprobtab->SetDM41Range2(nsteps_dm214_2, dm214start_2, dm214end_2);
    oscprobtab->ReadTable(Paths::outpath("OscProbTable.txt"));
  }

  static TMinuit* minu;
#pragma omp threadprivate(minu)
#pragma omp parallel
  {
    minu = new TMinuit(4);
    minu->SetPrintLevel(-1);
    minu->SetFCN(minuit_fcn);
  }

  static bool FixedTheMatrix = false;
#pragma omp threadprivate(FixedTheMatrix)

  ofstream fout(Paths::outpath("asimov_dchi2_%s.txt",
                               fit4nuSamples ? "4n" : "3n"));

#pragma omp parallel for
  for (int step_dm214 = 0; step_dm214 < nsteps_dm214_all; ++step_dm214) {
#pragma omp critical
    cout << "========" << step_dm214 << "/" << nsteps_dm214_all << "========"
         << endl;

    double dm214;
    if (step_dm214 < nsteps_dm214)
      dm214=ranger_dm41->returnVal(step_dm214);
    else
      dm214=ranger_dm41_2->returnVal(step_dm214-nsteps_dm214);

    for (int step_s22t14 = step_s22t14_start; step_s22t14 < step_s22t14_end; ++step_s22t14) {
      double s22t14 = exp(log(s22t14start) + log_s22t14_step*step_s22t14);
      double bests22t13, bests22t13_err;
      if (fit4nuSamples) {
        pred->LoadToyMCNominalSpec(Form("_s2t13_%4.4f_dm2ee_%5.5f_s2t14_%4.4f_dm214_%5.5f",
                                        S22T13, DM2EE, s22t14, dm214));
      }
      if (not FixedTheMatrix) {
        pred->SetSin22t13Step(20, 0.00, 0.20); //Set here!
        pred->FixCovMatrix(S22T13, DM2EE, 0., 0.1e-3);
        FixedTheMatrix = true;
      }
      // We are fitting toy data to the _opposite_ model
      const double model_s22t14 = fit4nuSamples ? 0 : s22t14;
      DoMinuitFit(minu, dm214, model_s22t14);
      minu->GetParameter(0, bests22t13, bests22t13_err);
      double chi2 = oscprobtab->CalculateChi2CovQuick(bests22t13, DM2EE,
                                                      model_s22t14, dm214,
                                                      PeriodFlag);
#pragma omp critical
      {
        fout << step_dm214 << " " << step_s22t14 << " " << chi2;
        if (fit4nuSamples)
          fout << " " << bests22t13;
        fout << endl;
      }

    }
  }

  fout.close();
  quick_exit(0);
}

// If s22t14 is -1, float it; otherwise, fix it
static void DoMinuitFit(TMinuit* minu, double dm214, double s22t14)
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
