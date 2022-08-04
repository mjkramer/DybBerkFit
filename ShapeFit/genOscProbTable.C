#include "Binning.h"
#include "FluxCalculator.h""
#include "OscProbTable.h"
#include "Paths.h"
#include "Predictor.h"

void genOscProbTable()
{
  auto pred = new Predictor();
  auto fluxcalc = new FluxCalculator(Paths::baselines(), Paths::histogram());
  pred->EnterFluxCalculator(fluxcalc);

  const Int_t n_evis_bins = Binning::n_evis();
  double* evis_bins = Binning::evis();

  const Int_t n_enu_bins = Binning::n_enu();
  double* enu_bins = Binning::enu();

  pred->SetEvisBins(n_evis_bins, evis_bins);
  pred->SetEnuBins(n_enu_bins, enu_bins);

  pred->LoadEvisToEnuMatrix(Paths::response());

  // const Int_t nsteps = 1;
  // Double_t s22t13start=0.084;
  // Double_t s22t13end=0.084;

  const Int_t nsteps_dm2 = 1;
  // Double_t dm2eestart=2.53e-3;
  // Double_t dm2eeend=2.53e-3;
  // consistent with data_file: (XXX update to FD result?)
  Double_t dm2eestart=2.501e-3;
  Double_t dm2eeend=2.501e-3;

  // const Int_t nsteps_s22t14 = 1;
  // Double_t s22t14start=0.0;
  // Double_t s22t14end=0.0;

  //grid in dm241 has to be same as that of toy MC generation -----------
  // (presumably "toy MC generation" => the asimov dchi2)
  const Int_t nsteps_dm214 = 120;
  const Int_t nsteps_dm214_2 = 100;
  const Int_t nsteps_dm214_all = nsteps_dm214 + nsteps_dm214_2;

  Double_t old_log_dm214_step = log(0.05/0.0005) / 99.0; //Set log step first
  Double_t dm214start=exp(log(0.0005)+old_log_dm214_step*(-20));
  Double_t dm214end=0.05;

  Double_t dm214start_2=0.05+ 0.25/100;
  Double_t dm214end_2=0.3;

  auto oscprobtab = new OscProbTable(pred);
  oscprobtab->SetDMeeRange(nsteps_dm2,dm2eestart,dm2eeend);
  oscprobtab->SetDM41Range(nsteps_dm214,dm214start,dm214end,true);
  oscprobtab->SetDM41Range2(nsteps_dm214_2,dm214start_2,dm214end_2);

  oscprobtab->GenerateTable();
  oscprobtab->WriteTable(Paths::outpath("OscProbTable.txt"));
}
