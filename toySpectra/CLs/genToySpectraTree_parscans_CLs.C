#include "Binning.h"
#include "Config.h"
#include "DataSet.h"
#include "Paths.h"
#include "Spectrum.h"

#include "TFile.h"
#include "TH1.h"

#include <cassert>
#include <iostream>

#ifdef __CLING__
extern int omp_get_thread_num();
#else
#include <omp.h>
#endif

using namespace Config;

void genToySpectraTree_parscans_CLs(const char *dataset_name,
                                    int itask,
                                    int ntasks)
{
  const int n_evis_bins = Binning::n_evis();
  double* evis_bins = Binning::evis();

  Char_t name[1024];
  TH1F* h_nominal_ad[Nstage][Ndetectors];
  for (int istage = 0; istage < Nstage; ++istage) {
    for (int idet = 0; idet < Ndetectors; ++idet) {
      cout << "Creating histogram for Stage# " << istage + 1 << ": AD# "
           << idet + 1 << endl;
      sprintf(name, "h_nominal_stage%i_ad%i", istage + 1, idet + 1);
      h_nominal_ad[istage][idet] = new TH1F(name, name, n_evis_bins, evis_bins);
    }
  }

  DataSet* mydata = new DataSet();
  mydata->load(Paths::toyconfig(dataset_name));

  const char *Theta13InputsLocation[3] = {Paths::input(0), Paths::input(1),
                                          Paths::input(2)};
  // Create Predictor (needed for livetimes, efficiencies, target masses... etc)
  static Predictor *myPred;
#pragma omp threadprivate(myPred)
#pragma omp parallel
  {
    myPred = new Predictor();
    for (int istage = 0; istage < Nstage; istage++) {
      myPred->LoadMainData(Theta13InputsLocation[istage]);
    }
  } // parallel

  static Spectrum *spectrumNormNominal;
#pragma omp threadprivate(spectrumNormNominal)
#pragma omp parallel
  {
    spectrumNormNominal = new Spectrum();
    spectrumNormNominal->passPredictor(myPred);
    // fixme: should use distances from Predictor (from FluxCalculator) in order
    // to avoid duplication
    spectrumNormNominal->loadDistances(Paths::baselines());
    spectrumNormNominal->initialize(mydata);
    spectrumNormNominal->loadBgSpecForToy();
    // NB: Adding 1 is important because 0 means random (UUID) seed
    // spectrumNormNominal->setRandomSeed(1 + omp_get_thread_num());
  } // parallel

  // override values from Config.h
  const int nsteps = 1;
  const double s22t13start = S22T13;
  const double s22t13end = S22T13;

  const int nsteps_dm2 = 1;
  const double dm2eestart = DM2EE;
  const double dm2eeend = DM2EE;

  // Uncomment the commented lines in the block below if nsteps(_dm2) > 1
  // And comment the following assert:
  assert(nsteps == 1 && nsteps_dm2 == 1);

  const int npoints = nsteps * nsteps_dm2 * nsteps_s22t14 * nsteps_dm214_all;

  Double_t s2t13_array[npoints + 1];
  Double_t dm2ee_array[npoints + 1];
  Double_t s2t14_array[npoints + 1];
  Double_t dm214_array[npoints + 1];
  Int_t dm214_id[npoints + 1];
  Int_t ipoint = 0;
  for (Int_t it14 = 0; it14 < nsteps_s22t14; it14++) {
    for (Int_t im41 = 0; im41 < nsteps_dm214_all; im41++) {
      for (Int_t it13 = 0; it13 < nsteps; it13++) {
        for (Int_t imee = 0; imee < nsteps_dm2; imee++) {
          s2t13_array[ipoint] = s22t13start;
          dm2ee_array[ipoint] = dm2eestart;

          // Uncomment lines lines if nsteps(_dm2) > 1
          // They cause a warning when nsteps is 1, so we commented them

          // if (nsteps > 1)
          //   s2t13_array[ipoint] =
          //       s22t13start + it13 * (s22t13end - s22t13start) / (nsteps - 1);
          // if (nsteps_dm2 > 1)
          //   dm2ee_array[ipoint] =
          //       dm2eestart + imee * (dm2eeend - dm2eestart) / (nsteps_dm2 - 1);

          s2t14_array[ipoint] = exp(log(s22t14start) + log_s22t14_step * it14);

          if (im41 < nsteps_dm214) {
            dm214_array[ipoint] = exp(log(dm214start) + log_dm214_step * im41);
          } else {
            dm214_array[ipoint] =
                dm214start_2 + lin_dm214_step * (im41 - nsteps_dm214);
          }
          dm214_id[ipoint] = im41 + 1;

          ipoint++;
        }
      }
    }
  }

  s2t13_array[npoints] = S22T13;
  dm2ee_array[npoints] = DM2EE;
  s2t14_array[npoints] = 0.0;
  dm214_array[npoints] = 0.0;
  dm214_id[npoints] = 0.0;

  // Note that the REAL number of points is npoints+1

  int points_per_task = (npoints+1) / ntasks;
  int ipar_start = itask * points_per_task;
  int ipar_end = ipar_start + points_per_task;
  // in case npoints+1 is not divisible by ntasks
  if ( itask == ntasks - 1 ) ipar_end = npoints+1;

  auto output_filename = Paths::outpath("toys_parscans/toySpectra_parscans_%s_%02d.root", dataset_name, itask);
  TFile *outfile = new TFile(output_filename, "RECREATE");

#pragma omp parallel for
  for (int ipar = ipar_start; ipar < ipar_end; ++ipar) {
    spectrumNormNominal->setOscillationSterile(
        s2t13_array[ipar], dm2ee_array[ipar], s2t14_array[ipar],
        dm214_array[ipar]);

    spectrumNormNominal->updateAntinu();
    spectrumNormNominal->updateBgDetected();

#pragma omp critical
    {
      for (int istage = 0; istage < Nstage; ++istage) {
        for (int idet = 0; idet < Ndetectors; ++idet) {
          auto name =
              Form("h_nominal_stage%i_ad%i_s2t13_%4.4f_dm2ee_%5.5f_s2t14_%4.4f_"
                   "dm214_%5.5f",
                   istage + 1, idet + 1, s2t13_array[ipar], dm2ee_array[ipar],
                   s2t14_array[ipar], dm214_array[ipar]);
          h_nominal_ad[istage][idet]->SetName(name);
          h_nominal_ad[istage][idet]->SetTitle(name);
          h_nominal_ad[istage][idet]->Reset();

          for (int ibin = 0; ibin < spectrumNormNominal->nSamples(); ++ibin) {
            h_nominal_ad[istage][idet]->Fill(
                spectrumNormNominal->energyArray(idet)[ibin],
                spectrumNormNominal->positronDetectedArray(istage, idet)[ibin] *
                    spectrumNormNominal->binWidth());
          }

          for (int ibin = 0; ibin < spectrumNormNominal->nSamplesBkg();
               ++ibin) {
            h_nominal_ad[istage][idet]->Fill(
                spectrumNormNominal->energyArrayBkg(idet)[ibin],
                spectrumNormNominal->bgDetectedArray(istage, idet)[ibin]);
          }

          outfile->cd();
          h_nominal_ad[istage][idet]->Write();
          cout << "WROTE " << name << endl;
        }
      }
      cout << endl;
    }
  }

  outfile->Write();
  outfile->Close();
}
