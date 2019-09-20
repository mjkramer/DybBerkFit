#include <iostream>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TTree.h"


#include "DataSet.h"
#include "Spectrum.h"


void genEvisToEnuMatrix(TString nominal_dataset_filename = "./data_file/dyb_data_v1_nominal_noosc.txt",
                        TString output_filename="../outputs/evis_to_enu_fine_2017Model_p17b_IHEP_BCWbin.root",
                        //CHANGE BINNING IF NAME IS WITH BCWbin
                        //do not forget change binning scheme accordingly to the name
			Double_t  s2t13 = -1,
			Double_t  dm2ee = -1,
			Double_t  s2t14 = -1,
			Double_t  dm241 = -1){

  //TString nominal_dataset_filename = "./data_file_unified_nl_p12e_unblinded/dyb_data_v1_nominal_noosc.txt",

  int nToys = 0;
  // string opt="_allerrors";

  // // define output histograms

  //  const Int_t n_evis_bins = 52;
    
    /*
  const Int_t n_evis_bins = 37;
  Double_t evis_bins[n_evis_bins+1]; // Single bins between 0.7 and 1.0 MeV. 0.2 MeV bins from 1.0 to 8.0 MeV. Single bin between 8.0 and 12 MeV. total 37 bins
  evis_bins[0] = 0.7;
  for (Int_t i = 0; i < n_evis_bins-1; i++){
    evis_bins[i+1] = 0.2 *i + 1.0;
  }
  evis_bins[n_evis_bins] = 12.0;*/
    
    //in case of BCW binning
   const Int_t n_evis_bins = 26;
    Double_t evis_bins[n_evis_bins+1]; // Single bins between 0.7 and 1.0 MeV. 0.2 MeV bins from 1.0 to 8.0 MeV. Single bin between 8.0 and 12 MeV. total 37 bins
    evis_bins[0] = 0.7;
    for (Int_t i = 0; i < n_evis_bins-1; i++){
        evis_bins[i+1] = 0.25 *i + 1.3;
    }
    evis_bins[n_evis_bins] = 12.0;


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
  mydata_nominal->load(nominal_dataset_filename);

  // Change oscillation parameters for testing
  if (s2t13 >= 0){
    mydata_nominal->setDouble("sinSq2Theta13",s2t13);
  }

  if (dm2ee > 0){
    mydata_nominal->setDouble("deltaMSqee",dm2ee);
  }

  if (s2t14 >= 0){
    mydata_nominal->setDouble("sinSq2Theta14",s2t14);
  }

  if (dm241 > 0){
    mydata_nominal->setDouble("deltaMSq41",dm241);
  }

  
  // Create Predictor (needed for livetimes, efficiencies, target masses... etc)
  Predictor *myPred = new Predictor();
  //myPred->LoadMainData("../ShapeFit/Inputs/Theta13-inputs_32week_inclusive.txt"); 

  Char_t Theta13InputsLocation[3][1024] = {"../ShapeFit/Inputs/Theta13-inputs_P17B_inclusive_6ad_IHEP.txt","../ShapeFit/Inputs/Theta13-inputs_P17B_inclusive_8ad_IHEP.txt","../ShapeFit/Inputs/Theta13-inputs_P17B_inclusive_7ad_IHEP.txt"};
  for(int istage=0;istage<Nstage;istage++){
    myPred->LoadMainData(Theta13InputsLocation[istage]); 
  }
  
  // Create nominal Spectrum 
  Spectrum *spectrumNormNominal = new Spectrum();
  spectrumNormNominal->passPredictor(myPred);
  //fixme: should use distances from Predictor (from FluxCalculator) in order to avoid duplication
  spectrumNormNominal->loadDistances("./unblinded_baseline.txt");                    
  spectrumNormNominal->initialize(mydata_nominal);

  TString AccidentalSpectrumLocation[3] = {"../ShapeFit/Spectra/accidental_eprompt_shapes_6ad_IHEP_BCWbin.root","../ShapeFit/Spectra/accidental_eprompt_shapes_8ad_IHEP_BCWbin.root","../ShapeFit/Spectra/accidental_eprompt_shapes_7ad_IHEP_BCWbin.root"};
  
    spectrumNormNominal->loadBgSpecForToy(AccidentalSpectrumLocation,
                                          "../li9_spectrum/8he9li_nominal_spectrum.root",
                                          "../amc_spectrum/amc_spectrum.root",
                                          "../fn_spectrum/P15A_fn_spectrum_IHEP.root",
                                          "../alpha-n-spectrum/result-DocDB9667.root");


  
  //Prepare destination file
  TFile *outfile = new TFile(output_filename.Data(),"RECREATE");

  
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
