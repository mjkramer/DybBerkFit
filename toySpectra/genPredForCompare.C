#include "Config.h"
#include "DataSet.h"
#include "Paths.h"
#include "Spectrum.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"

#include <fstream>
#include <iostream>


using namespace Config;


void genPredForCompare(bool noOffEq = false, bool noSNF = false)
{
  const char* outfilename;
  if (noOffEq && noSNF) {
    Paths::gReactorSuffix = "";
    outfilename = "predForCompare_noOffEq_noSNF.root";
  } else if (noOffEq) {
    Paths::gReactorSuffix = "_SNF";
    outfilename = "predForCompare_noOffEq.root";
  } else if (noSNF) {
    Paths::gReactorSuffix = "_nonEq";
    outfilename = "predForCompare_noSNF.root";
  } else {
    Paths::gReactorSuffix = "_SNF_nonEq";
    outfilename = "predForCompare.root";
  }

  TH1F* h_nu[Nstage][Ndetectors];
  TH1F* h_true[Nstage][Ndetectors];
  TH1F* h_dep[Nstage][Ndetectors];
  TH1F* h_rec[Nstage][Ndetectors];

  const int nbins = 240;

  for (int istage = 0; istage < Nstage; ++istage) {
    for (int idet = 0; idet < Ndetectors; ++idet) {
      const char* name;
      name = Form("h_nu_stage%d_ad%d", istage+1, idet+1);
      h_nu[istage][idet] = new TH1F(name, name, nbins, 0, 12);
      name = Form("h_true_stage%d_ad%d", istage+1, idet+1);
      h_true[istage][idet] = new TH1F(name, name, nbins, 0, 12);
      name = Form("h_dep_stage%d_ad%d", istage+1, idet+1);
      h_dep[istage][idet] = new TH1F(name, name, nbins, 0, 12);
      name = Form("h_rec_stage%d_ad%d", istage+1, idet+1);
      h_rec[istage][idet] = new TH1F(name, name, nbins, 0, 12);
    }
  }

  DataSet* mydata = new DataSet();
  mydata->load(Paths::nominal_toyconfig());
  mydata->setDouble("sinSq2Theta13", 0);
  mydata->setDouble("sinSq2Theta12", 0);
  mydata->setDouble("sinSq2Theta14", 0);

  Predictor* myPred = new Predictor();

  const char* Theta13InputsLocation[3] = {Paths::input(0), Paths::input(1),
    Paths::input(2)};

  for (int istage = 0; istage < Nstage; istage++)
    myPred->LoadMainData(Theta13InputsLocation[istage]);

  Spectrum* spectrum = new Spectrum();
  spectrum->passPredictor(myPred);
  spectrum->loadDistances(Paths::baselines());
  spectrum->initialize(mydata);
  spectrum->loadBgSpecForToy();

  TFile* outfile = new TFile(outfilename, "RECREATE");

  spectrum->updateAntinu();

  for (int istage = 0; istage < Nstage; ++istage) {
    for (int idet = 0; idet < Ndetectors; ++idet) {
      for (int ibin = 0; ibin < spectrum->nSamples(); ++ibin) {
        const double E = spectrum->energyArray(idet)[ibin];
        const double W = spectrum->binWidth();
        // NB: In genSuperHistograms we don't multiply by W, but it doesn't
        // matter here since we're interested in shape comparisons only
        // (absolute normalization is irrelevant).
        h_nu[istage][idet]->Fill(E,
          W * spectrum->antiNuNoOscArray(istage, idet)[ibin]);
        // In the other "gens" we do multiply by W.
        h_true[istage][idet]->Fill(E,
          W * spectrum->positronTrueArray(istage, idet)[ibin]);
        h_dep[istage][idet]->Fill(E,
          W * spectrum->positronIavDistortedArray(istage, idet)[ibin]);
        h_rec[istage][idet]->Fill(E,
          W * spectrum->positronDetectedArray(istage, idet)[ibin]);
      }
    }
  }

  outfile->cd();

  for (int istage = 0; istage < Nstage; ++istage) {
    for (int idet = 0; idet < Ndetectors; ++idet) {
      h_nu[istage][idet]->Write();
      h_true[istage][idet]->Write();
      h_dep[istage][idet]->Write();
      h_rec[istage][idet]->Write();
    }
  }
}
