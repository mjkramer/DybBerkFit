#include <iostream>
#include <fstream>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TTree.h"


#include "DataSet.h"
#include "Spectrum.h"


void genSuperHistograms(TString dataset_filename = "./data_file/dyb_data_v1_nominal.txt",
                        //TString output_filename="../ShapeFit/Flux/SuperHistograms_P17B_p15a_fine_huber-french.root"){
                        TString output_filename="../ShapeFit/Flux/SuperHistograms_P17B_2017Model_fine_huber-french_IHEP_official.root"){
  
  TString output_supermatrix_filename = "SuperMatrix_P15A_P15Aperiod_toy.txt";
  

  TH1F * h[Nstage][Ndetectors][Ncores];
  
  
  // Create Expectations
  DataSet *mydata = new DataSet();
  mydata->load(dataset_filename.Data());

  // Create Predictor (needed for livetimes, efficiencies, target masses... etc)
  Predictor *myPred = new Predictor();
 
  //Char_t Theta13InputsLocation[2][1024] = {"../ShapeFit/Inputs/Theta13-inputs_P15A_inclusive_6ad.txt","../ShapeFit/Inputs/Theta13-inputs_P15A_inclusive_8ad_p14a.txt"};
  Char_t Theta13InputsLocation[3][1024] = {"../ShapeFit/Inputs/Theta13-inputs_P17B_inclusive_6ad_IHEP.txt","../ShapeFit/Inputs/Theta13-inputs_P17B_inclusive_8ad_IHEP.txt","../ShapeFit/Inputs/Theta13-inputs_P17B_inclusive_7ad_IHEP.txt"};

  for(int istage=0;istage<Nstage;istage++){
    myPred->LoadMainData(Theta13InputsLocation[istage]); 
  }
  

  // Create nominal Spectrum 
  Spectrum *spectrumNorm = new Spectrum();
  spectrumNorm->passPredictor(myPred);
  //fixme: should use distances from Predictor (from FluxCalculator) in order to avoid duplication
  spectrumNorm->loadDistances("./unblinded_baseline.txt");                    
  spectrumNorm->initialize(mydata);
  //TString AccidentalSpectrumLocation[2] = {"../ShapeFit/Spectra/accidental_eprompt_shapes_6ad.root","../ShapeFit/Spectra/accidental_eprompt_shapes_8ad_p14a.root"};
    TString AccidentalSpectrumLocation[3] = {"../ShapeFit/Spectra/accidental_eprompt_shapes_6ad_IHEP_BCWbin.root","../ShapeFit/Spectra/accidental_eprompt_shapes_8ad_IHEP_BCWbin.root","../ShapeFit/Spectra/accidental_eprompt_shapes_7ad_IHEP_BCWbin.root"};

  spectrumNorm->loadBgSpecForToy(AccidentalSpectrumLocation,
				 "../li9_spectrum/8he9li_nominal_spectrum.root",
                                 "../amc_spectrum/amc_spectrum.root",
                                 "../fn_spectrum/P15A_fn_spectrum_IHEP.root",
                                 "../alpha-n-spectrum/result-DocDB9667.root");
  

  // load distances hare as well to construct traditional supermatrix
  
  Double_t m_detectorDistance[Ndetectors][Ncores];

  string dummyLine;
  string thead;
  float d2,d1,l2,l1,l4,l3;
  //-->Distances  
  cout << " Distances ++++++++++++++++++++++++++++++++++++++" << endl;
  ifstream disfile("./unblinded_baseline.txt");
  getline(disfile, dummyLine);
  while(disfile >> thead >> d1 >> d2 >> l1 >> l2 >> l3 >> l4){

    cout << thead << "\t" << d1 << "\t" << d2 << "\t" << l1 << "\t" << l2 << "\t" << l3 << "\t" << l4 << endl;//tmp
    int adnum=atoi(thead.substr(2,1).c_str());
    //note: must convert into kms
    m_detectorDistance[adnum-1][0]=d1*0.001;
    m_detectorDistance[adnum-1][1]=d2*0.001;
    m_detectorDistance[adnum-1][2]=l1*0.001;
    m_detectorDistance[adnum-1][3]=l2*0.001;
    m_detectorDistance[adnum-1][4]=l3*0.001;
    m_detectorDistance[adnum-1][5]=l4*0.001;
  }
  

  //Double_t M[Nstage][Ndetectors][Ncores];

  
  
  //Prepare destination file
  TFile *outfile = new TFile(output_filename.Data(),"RECREATE");

  Char_t name[1024];
  for(int istage=0;istage<Nstage;++istage){
  for(int idet=0;idet<Ndetectors;++idet){
    for(int icore=0;icore<Ncores;++icore){
      cout << "Creating histogram for Stage# " << istage+1 << " AD# " << idet+1 << " from Core " << icore << endl;
      sprintf(name,"h_super_istage%d_idet%d_icore%d",istage,idet,icore);
      //h[istage][idet][icore] = new TH1F (name,name,41,1.8,10);
      h[istage][idet][icore] = new TH1F (name,name,164,1.8,10);
    }
  }
  }
  
  //Generate nominal spectrum
  
  for(int istage=0;istage<Nstage;++istage){
  for(int icore=0;icore<Ncores;++icore){
    spectrumNorm->updateAntinu(icore);
    for(int idet=0;idet<Ndetectors;++idet){
      for(int ibin=0;ibin<spectrumNorm->nSamples();++ibin){
        
        h[istage][idet][icore]->Fill(spectrumNorm->energyArray(idet)[ibin],
				     spectrumNorm->antiNuNoOscArray(istage,idet)[ibin]);
      }//ibin loop


      // construct supermatrix element
      //M[istage][idet][icore] = h[istage][idet][icore]->Integral()/pow(m_detectorDistance[idet][icore],2);
      
    }
    //    cout << "AD" << idet+1 << ": " << h_nominal_ad[idet]->Integral() << endl;//tmp
  }
  }
  
  outfile->Write();
  outfile->Close();



  /*
  // normalize and dump supermatrix
  ofstream fout_supermatrix(output_supermatrix_filename.Data());
  fout_supermatrix << "#\tD1\tD2\tL1\tL2\tL3\tL4" << endl;
  
  for(int idet=0;idet<Ndetectors;++idet){
    Double_t norm = 0;
    for(int icore=0;icore<Ncores;++icore)
      norm+= M[idet][icore];
    fout_supermatrix << Form("AD%d",idet+1);
    for(int icore=0;icore<Ncores;++icore){
      M[idet][icore] /= norm;
      fout_supermatrix << "\t" << M[idet][icore];
    }
    fout_supermatrix << endl;
  }
  */
  

}//end of genToySpectraTree