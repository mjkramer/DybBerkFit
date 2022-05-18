#include "Binning.h"
#include "Config.h"
#include "DataSet.h"
#include "Paths.h"
#include "Spectrum.h"

#include "TFile.h"
#include "TH1.h"

void genAsimovToy(double s2t13, double dm2ee)
{
  const int n_evis_bins = Binning::n_evis();
  double* evis_bins = Binning::evis_fine();

  TH1F* h_eprompt[Nstage][Ndetectors];

  for (int istage = 0; istage < Nstage; ++istage) {
    for (int idet = 0; idet < Ndetectors; ++idet) {
      const char* name = Form("h_ibd_eprompt_inclusive_eh%d_ad%d",
                              detConfigEH[idet], detConfigAD[idet]);
      h_eprompt[istage][idet] = new TH1F(name, name, n_evis_bins, evis_bins);
    }
  }

  DataSet* mydata_nominal = new DataSet;
  mydata_nominal->load(Paths::nominal_toyconfig());
  mydata_nominal->setDouble("sinSq2Theta13", s2t13);
  mydata_nominal->setDouble("deltaMSqee", dm2ee);

  const char* Theta13InputsLocation[3] =
    {Paths::input(0), Paths::input(1), Paths::input(2)};

  Predictor* myPred = new Predictor;
  for (int istage = 0; istage < Nstage; ++istage)
    myPred->LoadMainData(Theta13InputsLocation[istage]);

  Spectrum* spectrumNormNominal = new Spectrum;
  spectrumNormNominal->passPredictor(myPred);
  spectrumNormNominal->loadDistances(Paths::baselines());
  spectrumNormNominal->initialize(mydata_nominal);
  spectrumNormNominal->loadBgSpecForToy();

  spectrumNormNominal->updateAntinu();
  spectrumNormNominal->updateBgDetected();

  for (int istage = 0; istage < Nstage; ++istage) {
    for (int idet = 0; idet < Ndetectors; ++idet) {
      for (int ibin = 0; ibin < spectrumNormNominal->nSamples(); ++ibin) {
        h_eprompt[istage][idet]
          ->Fill(spectrumNormNominal->energyArray(idet)[ibin],
                 spectrumNormNominal->positronDetectedArray(istage, idet)[ibin]
                 * spectrumNormNominal->binWidth());
      }

      for (int ibin = 0; ibin < spectrumNormNominal->nSamplesBkg(); ++ibin) {
        // do not multiply by binWidth for backgrounds
        h_eprompt[istage][idet]
          ->Fill(spectrumNormNominal->energyArrayBkg(idet)[ibin],
                 spectrumNormNominal->bgDetectedArray(istage, idet)[ibin]);
      }
    }
  }

  for (int istage = 0; istage < Nstage; ++istage) {
    TFile* outfile = new TFile(Paths::sig_spectra(istage), "RECREATE");
    outfile->cd();
    for (int idet = 0; idet < Ndetectors; ++idet) {
      h_eprompt[istage][idet]->Write();
    }
    outfile->Write();
    outfile->Close();
  }
}
