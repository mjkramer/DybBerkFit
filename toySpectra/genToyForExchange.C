#include "Binning.h"
#include "Config.h"
#include "DataSet.h"
#include "Paths.h"
#include "Spectrum.h"

#include "TFile.h"
#include "TH1.h"

// output_filename_format should be, e.g., "LBNL_%s_Oscillation1.root".
// We will replace %s with 6AD, etc.
void genToyForExchange(const char* output_filename_format,
                       double s2t13, double dm2ee)
{
  const int n_evis_bins = Binning::n_evis_fine();
  double* evis_bins = Binning::evis_fine();

  TH1F* h_nominal_ad[Nstage][Ndetectors];

  for (int istage = 0; istage < Nstage; ++istage) {
    for (int idet = 0; idet < Ndetectors; ++idet) {
      const char* name = Form("h_nominal_stage%i_ad%i", istage + 1, idet + 1);
      h_nominal_ad[istage][idet] = new TH1F(name, name, n_evis_bins, evis_bins);
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
        h_nominal_ad[istage][idet]
          ->Fill(spectrumNormNominal->energyArray(idet)[ibin],
                 spectrumNormNominal->positronDetectedArray(istage, idet)[ibin]
                 * spectrumNormNominal->binWidth());
      }
      for (int ibin = 0; ibin < spectrumNormNominal->nSamplesBkg(); ++ibin) {
        h_nominal_ad[istage][idet]
          ->Fill(spectrumNormNominal->energyArrayBkg(idet)[ibin],
                 spectrumNormNominal->bgDetectedArray(istage, idet)[ibin]);
        // do not multiply by binWidth for backgrounds
      }
    }
  }

  const char* istage2name[3] = {"6AD", "8AD", "7AD"};

  for (int istage = 0; istage < Nstage; ++istage) {
    const char* outname = Form(output_filename_format, istage2name[istage]);
    TFile *outfile = new TFile(outname, "RECREATE");
    outfile->cd();
    for (int idet = 0; idet < Ndetectors; ++idet)
      h_nominal_ad[istage][idet]->Write();
    outfile->Write();
    outfile->Close();
  }
}
