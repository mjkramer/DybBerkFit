#include <iostream>

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TTree.h"

#include "DataSet.h"
#include "Spectrum.h"
#include "Config.h"
#include "Binning.h"

using namespace Config;
using namespace std;

void genEvisToEnuMatrix(Double_t  s2t13 = -1, Double_t  dm2ee = -1, Double_t  s2t14 = -1, Double_t  dm241 = -1){

  // string opt="_allerrors";

  // // define output histograms

  //  const Int_t n_evis_bins = 52;

  const Int_t n_evis_bins = Binning::n_evis();
  const double* evis_bins = Binning::evis();

  Char_t name[1024];
  TH1F *h_nominal_ad[Ndetectors];
  TH2F *h_evis_vs_enu_ad[Ndetectors];
  //for(int idet=0;idet<Ndetectors;++idet){
  for(int idet=0;idet<1;++idet){
    cout << "Creating histogram for AD# " << idet+1 << endl;
    sprintf(name,"h_nominal_ad%i",idet+1);
    h_nominal_ad[idet] = new TH1F (name,name,n_evis_bins,evis_bins);
    sprintf(name,"h_evis_vs_enu_ad%i",idet+1);
    h_evis_vs_enu_ad[idet] = new TH2F (name,name,240,0,12,n_evis_bins,evis_bins);
  }

  // Create Nominal (i.e. no random variation) Expectations
  DataSet *mydata_nominal = new DataSet();
  mydata_nominal->load(nominal_noosc_dataset_filename);

  // Change oscillation parameters for testing
  if (s2t13 >= 0){
    mydata_nominal->setDouble("sinSq2Theta13",s2t13);
  } else
    s2t13 = mydata_nominal->getDouble("sinSq2Theta13");

  if (dm2ee > 0){
    mydata_nominal->setDouble("deltaMSqee",dm2ee);
  } else
    dm2ee = mydata_nominal->getDouble("deltaMSqee");

  if (s2t14 >= 0){
    mydata_nominal->setDouble("sinSq2Theta14",s2t14);
  } else
    s2t14 = mydata_nominal->getDouble("sinSq2Theta14");

  if (dm241 > 0){
    mydata_nominal->setDouble("deltaMSq41",dm241);
  } else
    dm241 = mydata_nominal->getDouble("deltaMSq41");


  // Create Predictor (needed for livetimes, efficiencies, target masses... etc)
  Predictor *myPred = new Predictor();

  const char* Theta13InputsLocation[3] = {input_filename0, input_filename1, input_filename2};
  for(int istage=0;istage<Nstage;istage++){
    myPred->LoadMainData(Theta13InputsLocation[istage]);
  }

  // Create nominal Spectrum
  Spectrum *spectrumNormNominal = new Spectrum();
  spectrumNormNominal->passPredictor(myPred);
  //fixme: should use distances from Predictor (from FluxCalculator) in order to avoid duplication
  spectrumNormNominal->loadDistances(baselines_filename);
  spectrumNormNominal->initialize(mydata_nominal);

  TString AccidentalSpectrumLocation[3] = {acc_spectra_filename0,acc_spectra_filename1,acc_spectra_filename2};

  spectrumNormNominal->loadBgSpecForToy(AccidentalSpectrumLocation,
                                        li9_filename,
                                        amc_filename,
                                        fn_filename,
                                        aln_filename);

  //Prepare destination file
  TFile *outfile = new TFile(response_root_filename, "RECREATE");


  //Generate nominal spectrum

  spectrumNormNominal->updateAntinu();
  //spectrumNormNominal->updateBgDetected();
  /*
  //for(int idet=0;idet<Ndetectors;++idet){
  for(int idet=0;idet<1;++idet){
  h_nominal_ad[idet]->Reset();
  for(int ibin=0;ibin<spectrumNormNominal->nSamples();++ibin){

  h_nominal_ad[idet]->Fill(spectrumNormNominal->energyArray(idet)[ibin],
  spectrumNormNominal->positronDetectedArray(0,idet)[ibin]*spectrumNormNominal->binWidth()
  + spectrumNormNominal->bgDetectedArray(0,idet)[ibin]
  );
  }//ibin loop
  for(int ibin=0;ibin<spectrumNormNominal->nSamplesBkg();++ibin){
  h_nominal_ad[idet]->Fill(spectrumNormNominal->energyArrayBkg(idet)[ibin],
  spectrumNormNominal->bgDetectedArray(0,idet)[ibin]
  );
  }


  cout << "AD" << idet+1 << ": " << h_nominal_ad[idet]->Integral() << endl;//tmp
  }
  */

  cout << "Generating 2D histograms. This takes a while..." << endl;
  //  create 2D histograms
  //for(int idet=0;idet<Ndetectors;++idet){
  for(int idet=0;idet<1;++idet){
    //    for(int ibin_enu=0;ibin_enu<spectrumNormNominal->nSamples();++ibin_enu){
    cout << "AD" << idet+1 << endl;
    TAxis * xa = h_evis_vs_enu_ad[idet]->GetXaxis();
    for(int ibin_enu=0;ibin_enu< xa->GetNbins();++ibin_enu){
      cout << ibin_enu << " / " << xa->GetNbins() << endl;
      spectrumNormNominal->updatePositronTrue(xa->GetBinLowEdge(ibin_enu+1),xa->GetBinUpEdge(ibin_enu+1));

      for(int ibin=0;ibin<spectrumNormNominal->nSamples();++ibin){
        //        h_evis_vs_enu_ad[idet]->Fill(spectrumNormNominal->energyArray(idet)[ibin_enu],
        h_evis_vs_enu_ad[idet]->Fill(xa->GetBinCenter(ibin_enu+1),
                                     spectrumNormNominal->energyArray(idet)[ibin],
                                     spectrumNormNominal->positronDetectedArray(0,idet)[ibin]*spectrumNormNominal->binWidth()
                                     ); // don't include background spectra for this true enu_vs_evis convertion.

      }
    }//idet loop

  } //ibin_enu loop


  outfile->cd();
  //for(int idet=0;idet<Ndetectors;++idet){
  for(int idet=0;idet<1;++idet){
    h_nominal_ad[idet]->Write();
    h_evis_vs_enu_ad[idet]->Write();
  }

  outfile->Write();
  outfile->Close();

}//end of genToySpectraTree
