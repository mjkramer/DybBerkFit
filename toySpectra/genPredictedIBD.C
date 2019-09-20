#include <iostream>
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TTree.h"
#include "DataSet.h"
#include "Spectrum.h"

void genPredictedIBD(Double_t  s2t13 = 0,
		     Double_t  dm2ee = -1,
             //TString output_filename="../ShapeFit/PredictedIBD/PredictedIBD_bestfit_P17B.root",
		     TString output_filename="../ShapeFit/PredictedIBD/PredictedIBD_P17B_2017Model_fine_huber-french_IHEP_official.root",
             //TString output_filename="../ShapeFit/PredictedIBD/PredictedIBD_P17B_fine_huber-french.root",
		     TString nominal_dataset_filename = "./data_file/dyb_data_v1_nominal.txt",
		     TString dataset_filename = "./data_file/dyb_data_v1_nominal.txt",
		     Double_t  s2t14 = 0,
		     Double_t  dm241 = -1){
 
 
    //  const Int_t n_evis_bins = 52;
    
    
    const Int_t n_evis_bins = 240;
    Double_t evis_bins[n_evis_bins+1]; // Single bins between 0.7 and 1.0 MeV. 0.2 MeV bins from 1.0 to 8.0 MeV. Single bin between 8.0 and 12 MeV. total 37 bins
    evis_bins[0] = 0.0;
    for (Int_t i = 0; i < n_evis_bins-1; i++){
        evis_bins[i+1] = 0.05 *i + 0.05;
    }
    evis_bins[n_evis_bins] = 12.0;
    /*
    const Int_t n_evis_bins = 26;
    Double_t evis_bins[n_evis_bins+1]; // Single bins between 0.7 and 1.0 MeV. 0.2 MeV bins from 1.0 to 8.0 MeV. Single bin between 8.0 and 12 MeV. total 37 bins
    evis_bins[0] = 0.7;
    for (Int_t i = 0; i < n_evis_bins-1; i++){
        evis_bins[i+1] = 0.25*i + 1.3;
    }
    evis_bins[n_evis_bins] = 12.0;*/
    
    
    Char_t name[1024];
    TH1F *h_nominal_ad[Nstage][Ndetectors];
    TH1F *h_ad[Nstage][Ndetectors];
    for(int istage=0;istage<Nstage;++istage){
        for(int idet=0;idet<Ndetectors;++idet){
            cout << "Creating histogram for Stage# " << istage+1 << ": AD# " << idet+1 << endl;
            sprintf(name,"h_stage%i_ad%i",istage+1,idet+1);
            h_ad[istage][idet] = new TH1F (name,name,n_evis_bins,evis_bins);
            sprintf(name,"h_nominal_stage%i_ad%i",istage+1,idet+1);
            h_nominal_ad[istage][idet] = new TH1F (name,name,n_evis_bins,evis_bins);
        }
    }
    
    // Create Expectations
    DataSet *mydata = new DataSet();
    mydata->load(dataset_filename.Data());
    
    // Create Nominal (i.e. no random variation) Expectations
    DataSet *mydata_nominal = new DataSet();
    mydata_nominal->load(nominal_dataset_filename);
    
    // Change oscillation parameters for testing
    if (s2t13 >= 0){
        cout<<"Changing theta"<<endl;
        mydata->setDouble("sinSq2Theta13",s2t13);
        mydata_nominal->setDouble("sinSq2Theta13",s2t13);
    }
    
    if (dm2ee > 0){
        mydata->setDouble("deltaMSqee",dm2ee);
        mydata_nominal->setDouble("deltaMSqee",dm2ee);
    }
    
    if (s2t14 >= 0){
        mydata->setDouble("sinSq2Theta14",s2t14);
        mydata_nominal->setDouble("sinSq2Theta14",s2t14);
    }
    
    if (dm241 > 0){
        mydata->setDouble("deltaMSq41",dm241);
        mydata_nominal->setDouble("deltaMSq41",dm241);
    }
    
    
    // Create Predictor (needed for livetimes, efficiencies, target masses... etc)
    Predictor *myPred = new Predictor();
    
    //Char_t Theta13InputsLocation[2][1024] = {"../ShapeFit/Inputs/Theta13-inputs_P15A_inclusive_6ad.txt","../ShapeFit/Inputs/Theta13-inputs_P15A_inclusive_8ad_p14a.txt"};
    Char_t Theta13InputsLocation[3][1024] = {"../ShapeFit/Inputs/Theta13-inputs_P17B_inclusive_6ad_IHEP.txt","../ShapeFit/Inputs/Theta13-inputs_P17B_inclusive_8ad_IHEP.txt","../ShapeFit/Inputs/Theta13-inputs_P17B_inclusive_7ad_IHEP.txt"};
    
    cout<<"There are "<<Nstage<<" stages"<<endl;
    
    for(int istage=0;istage<Nstage;istage++){
        myPred->LoadMainData(Theta13InputsLocation[istage]);
    }
    
    
    
    // Create nominal Spectrum
    Spectrum *spectrumNormNominal = new Spectrum();
    spectrumNormNominal->passPredictor(myPred);
    //fixme: should use distances from Predictor (from FluxCalculator) in order to avoid duplication
    spectrumNormNominal->loadDistances("./unblinded_baseline.txt");
    spectrumNormNominal->initialize(mydata_nominal);
    
    //TString AccidentalSpectrumLocation[2] = {"../ShapeFit/Spectra/accidental_eprompt_shapes_6ad.root","../ShapeFit/Spectra/accidental_eprompt_shapes_8ad_p14a.root"};
    TString AccidentalSpectrumLocation[3] = {"../ShapeFit/Spectra/accidental_eprompt_shapes_6ad_IHEP_BCWbin.root","../ShapeFit/Spectra/accidental_eprompt_shapes_8ad_IHEP_BCWbin.root","../ShapeFit/Spectra/accidental_eprompt_shapes_7ad_IHEP_BCWbin.root"};
    
    spectrumNormNominal->loadBgSpecForToy(AccidentalSpectrumLocation,
                                          "../li9_spectrum/8he9li_nominal_spectrum.root",
                                          "../amc_spectrum/amc_spectrum.root",
                                          "../fn_spectrum/P15A_fn_spectrum_IHEP.root",
                                          "../alpha-n-spectrum/result-DocDB9667.root");
    
    
    //Prepare destination file
    TFile *outfile = new TFile(output_filename.Data(),"RECREATE");
    
    //Generate nominal spectrum
    // cout << "-------------------- Nominal Spectrum -------------------" << endl;
    spectrumNormNominal->updateAntinu();
    
    for(int istage=0;istage<Nstage;++istage){
        for(int idet=0;idet<Ndetectors;++idet){
            h_nominal_ad[istage][idet]->Reset();
            for(int ibin=0;ibin<spectrumNormNominal->nSamples();++ibin){
                
                //Don't fill below 0.7 MeV
                if (spectrumNormNominal->energyArray(idet)[ibin] < 0.7) continue;
                
                h_nominal_ad[istage][idet]->Fill(spectrumNormNominal->energyArray(idet)[ibin],
                                                 spectrumNormNominal->positronDetectedArray(istage,idet)[ibin]*spectrumNormNominal->binWidth()
                                                 );
                // cout << "The added positron at energy " << spectrumNormNominal->energyArray(idet)[ibin] << " is " << spectrumNormNominal->positronDetectedArray(idet)[ibin] << endl;//tmp
            }//ibin loop
            
            //cout << "AD" << idet+1 << ": " << h_nominal_ad[idet]->Integral() << endl;//tmp
        }
    }
    
    outfile->cd();
    
    for(int istage=0;istage<Nstage;++istage){
        for(int idet=0;idet<Ndetectors;++idet){
            h_nominal_ad[istage][idet]->Write();
        }
    }
    
    outfile->Write();
    outfile->Close();
    
    gDirectory->Delete("*");
    
}//end of genPredictedIBD

