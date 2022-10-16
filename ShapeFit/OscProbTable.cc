#include <iostream>
#include <iomanip>
#include <assert.h>
#include <fstream>
#include <sstream>

#include "OscProbTable.h"

using namespace std;

OscProbTable::OscProbTable(int method){

  // NOTE: there are two methods defined.
  // Method 0: Store oscillation probability for each reconstructed energy bin. This method is much faster but the result is not correct in this version.
  // Method 1: Store oscillation probability for each true energy bin. This is slower than method 0, but following the original method in Predictor.cc, and the result is correct.
  
  generatedTable=false;
  ranger_dmee = new Ranger();
  ranger_dm41 = new Ranger();
  ranger_dm41_2 = new Ranger();
  //m_superpred = new PredSet();
  m_method = method;
  cout << "Warning: you used the default constructor of OscProbTable; You need to pass it a Predictor object or it won't work properly." << endl;
  if (m_method != 0 || m_method != 1){
    cout << "Warning: Not supported method for OscProbTable: method " << m_method << endl;
    cout << "         Force setting to method 1...." << endl;
    m_method = 1;
  }


};

OscProbTable::OscProbTable(Predictor *pred_in, int method){

  generatedTable=false;
  ranger_dmee = new Ranger();
  ranger_dm41 = new Ranger();
  ranger_dm41_2 = new Ranger();
  //m_superpred = new PredSet();
  pred = pred_in;
  m_method = method;
  if (m_method != 0 && m_method != 1){
    cout << "Warning: Not supported method for OscProbTable: method " << m_method << endl;
    cout << "         Force setting to method 1...." << endl;
    m_method = 1;
  }

};

void OscProbTable::passPredictor(Predictor *pred_in){

  pred = pred_in;
  generatedTable=false;

};

OscProbTable::~OscProbTable(){
  delete pred;

};

void OscProbTable::SetMethod(int method){
  m_method = method;
  if (m_method != 0 && m_method != 1){
    cout << "Warning: Not supported method for OscProbTable: method " << m_method << endl;
    cout << "         Force setting to method 1...." << endl;
    m_method = 1;
  }
}


void OscProbTable::SetDMeeRange(int nsteps_in, double min_in, double max_in){
 
  ranger_dmee->nsteps=nsteps_in;
  ranger_dmee->min=min_in;
  ranger_dmee->max=max_in;  
  generatedTable=false;//<--if change range, force regeneration of table

};

void OscProbTable::SetDM41Range(int nsteps_in, double min_in, double max_in, bool logscale){

  ranger_dm41->nsteps=nsteps_in;
  ranger_dm41->min=min_in;
  ranger_dm41->max=max_in;
  if (logscale) ranger_dm41->setLogScale();
  
  generatedTable=false;//<--if change range, force regeneration of table

};

void OscProbTable::SetDM41Range2(int nsteps_in, double min_in, double max_in, bool logscale){

  ranger_dm41_2->nsteps=nsteps_in;
  ranger_dm41_2->min=min_in;
  ranger_dm41_2->max=max_in;
  if (logscale) ranger_dm41_2->setLogScale();
  
  generatedTable=false;//<--if change range, force regeneration of table

};



void OscProbTable::GenerateTable(){

  
  if(ranger_dmee->nsteps==0 || ranger_dm41->nsteps==0){
    cout << "ERROR! ranges have not been properly set! Use SetXXXXRange() for all variables.... Aborting" << endl;
  }
  //Make sure spectra are corrected (i.e. bg removed... etc) before calculating extrapolation table


  // loading evis to enu matrix locally for faster access
  for(int ibin_evis=0; ibin_evis < pred->GetNEvisBins(); ++ibin_evis){
    for(int ibin_enu=0; ibin_enu < pred->GetNEnuBins(); ++ibin_enu){
      M_evis_to_enu[ibin_evis][ibin_enu] = pred->EvisToEnuMatrix(ibin_evis,ibin_enu);
    }
  }
        
  
  //Extrapolation table generation
  cout << "Starting extrapolation table calculation... this may take a while" << endl;
  
  Int_t istep = 0;
  Int_t total_dm41_steps = ranger_dm41->nsteps + ranger_dm41_2->nsteps;
  Int_t total_steps = ranger_dmee->nsteps * total_dm41_steps;
  
  for(int idmee=0;idmee<ranger_dmee->nsteps;++idmee){
    double dmee=ranger_dmee->returnVal(idmee);
    for(int idm41=0;idm41<total_dm41_steps;++idm41){
      double dm41 = 0;
      if (idm41 < ranger_dm41->nsteps){
        dm41=ranger_dm41->returnVal(idm41);
      }else{
        dm41=ranger_dm41_2->returnVal(idm41 - ranger_dm41->nsteps);
      }
      cout << "\r" << "Extrapolation table progress: " << (int)(100*(Double_t)istep/(Double_t)total_steps) << "%" << flush;
      istep++;

      //make predictions
      for(int istage = 0;istage<Nstage;istage++){
        for(int idet=0;idet<Ndetectors;++idet){

          for (int iterm = 0; iterm < num_terms; iterm++){
            std::map<int,TH1F*> mapflux = pred->GetFluxCalculator()->CalculateFluxHistRow(idet,1,istage,dmee,1,dm41,iterm);

            double * enu_bins = pred->GetEnuBins();
            if (m_method == 0){
              for(int ibin=1;ibin<=pred->GetNEvisBins();++ibin){
                OscProbTableLookup[istage][idet][idmee][idm41][ibin-1][iterm] = 0;
                for(int icore=0;icore<Ncores;++icore){
                  for(int ibin_enu=1;ibin_enu<=pred->GetNEnuBins();++ibin_enu){
                    double Etrue = 0.5*(enu_bins[ibin_enu-1]+enu_bins[ibin_enu]);
                    //                  if (Etrue > 8.7) Etrue = 8.7; // this is because the current flux histogram ends at 8.8 MeV
                    int thebin=mapflux[icore]->GetXaxis()->FindBin(Etrue);
                    double factor=mapflux[icore]->GetBinContent(thebin);
                    OscProbTableLookup[istage][idet][idmee][idm41][ibin-1][iterm] += M_evis_to_enu[ibin-1][ibin_enu-1] * factor;

                  }
                }
              }
            }else if (m_method == 1){
              for(int ibin_enu=1;ibin_enu<=pred->GetNEnuBins();++ibin_enu){
                OscProbTableLookup[istage][idet][idmee][idm41][ibin_enu-1][iterm] = 0;
                for(int icore=0;icore<Ncores;++icore){
                  double Etrue = 0.5*(enu_bins[ibin_enu-1]+enu_bins[ibin_enu]);
                  //                if (Etrue > 8.7) Etrue = 8.7; // this is because the current flux histogram ends at 8.8 MeV
                  int thebin=mapflux[icore]->GetXaxis()->FindBin(Etrue);
                  double factor=mapflux[icore]->GetBinContent(thebin);
                  OscProbTableLookup[istage][idet][idmee][idm41][ibin_enu-1][iterm] += factor;

                }
              }
            }
          }
        }//idet
      }//istage
    }//idm41
  }//idmee
  
  cout << " " << endl;
  
  generatedTable=true;

}//generateTable

void OscProbTable::WriteTable(const char * outfilename){
 
  ofstream outfile;
  outfile.open(outfilename);

  Int_t total_dm41_steps = ranger_dm41->nsteps + ranger_dm41_2->nsteps;
  
  for(int istage = 0;istage<Nstage;istage++){
    for(int idet=0;idet<Ndetectors;++idet){
      
      for(int idmee=0;idmee<ranger_dmee->nsteps;++idmee){
        for(int idm41=0;idm41<total_dm41_steps;++idm41){

          for(int ibin_enu=1;ibin_enu<=pred->GetNEnuBins();++ibin_enu){
            for (int iterm = 0; iterm < num_terms; iterm++){
              outfile << istage << "\t"
                      << idet << "\t"
                      << idmee << "\t"
                      << idm41 << "\t"
                      << ibin_enu << "\t"
                      << iterm << "\t"
                      << std::setprecision(16) << OscProbTableLookup[istage][idet][idmee][idm41][ibin_enu-1][iterm] <<"\n";
            }//iterm
          }//ibin_enu

        }//idm41
      }//idmee

    }//idet
  }//istage

  outfile.close();

}//WriteTable

void OscProbTable::ReadTable(const char * infilename){
 
  // loading evis to enu matrix locally for faster access
  for(int ibin_evis=0; ibin_evis < pred->GetNEvisBins(); ++ibin_evis){
    for(int ibin_enu=0; ibin_enu < pred->GetNEnuBins(); ++ibin_enu){
      M_evis_to_enu[ibin_evis][ibin_enu] = pred->EvisToEnuMatrix(ibin_evis,ibin_enu);
    }
  }

  cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "+++++++++++++++++++ Reading OscProbTable File +++++++++++++++++" << endl;

  cout << "Reading from: " << infilename << endl;

  string line;
  Int_t istage = 0, idet = 0, idmee = 0, idm41 = 0, ibin_enu = 0, iterm = 0;
  Double_t OscProbTableValue = 0;

  ifstream infile(infilename);
  while(!infile.eof()){

    getline(infile,line);
    string firstchar=line.substr(0,1);
    if(firstchar=="#") continue;//<-- ignore lines with comments
    
    istringstream iss(line);
    int column=0;
  
    while(iss){
      string sub; iss >> sub;
      if(column==0) istage=atof(sub.c_str());
      if(column==1) idet=atof(sub.c_str());
      if(column==2) idmee=atof(sub.c_str());
      if(column==3) idm41=atof(sub.c_str());
      if(column==4) ibin_enu=atof(sub.c_str());
      if(column==5) iterm=atof(sub.c_str());
      if(column==6) OscProbTableValue=atof(sub.c_str());
      column+=1;
    }//looping over columns

    OscProbTableLookup[istage][idet][idmee][idm41][ibin_enu-1][iterm] = OscProbTableValue;
   
  }//end of file

  generatedTable=true;

}//ReadTable

// TH1F *OscProbTable::PlotExtrapolationHist(int idet_far, int idet_near, int is22t13, int is22t14, int idmee, int idm41){

//   if(!generatedTable) GenerateTable();
  
//   TH1F *hout = (TH1F*)pred->tdper[0].CorrEvtsSpec[0]->Clone("hout");
//   hout->Reset();
//   for(int ibin=1;ibin<=hout->GetXaxis()->GetNbins();++ibin){
//     hout->SetBinContent(ibin,OscProbTableLookup[istage][idet_far][idet_near][is22t13][is22t14][idmee][idm41][ibin-1]);
//   }//ibin
  
//   return hout;

// }//PlotExtrapolationFactor

void OscProbTable::MakeAllPredictionsQuick(double s22t13, double dmee, double s22t14, double dm41, bool print){
 
  if(!generatedTable) GenerateTable();
  
  //crash if query is outside range of grid
  if(ranger_dmee->nsteps>1 && (dmee<ranger_dmee->min || dmee>ranger_dmee->max)){
    cout << "ERROR!! requested dmee in OscProbTable is outside of grid's range... aborting!" << endl;
    exit(0);
  }

  
  double dm41_test = dm41;
  double dm41_test_2 = dm41;
  if (ranger_dm41->getLogScale()){
    dm41_test = log10(dm41);
  }
  if (ranger_dm41_2->getLogScale()){
    dm41_test_2 = log10(dm41);
  }

  //-->find the best indices (note to self: could consider adding warning to Ranger::findIndex when you do not land on a grid point). 
  int idmee = ranger_dmee->findIndex(dmee);
  int idm41 = -1;
  if(ranger_dm41->nsteps>1 && (dm41_test>=ranger_dm41->min - 1.0e-10 &&  dm41_test<=ranger_dm41->max + 1.0e-10)){
    idm41 = ranger_dm41->findIndex(dm41);
  }else if (ranger_dm41_2->nsteps>1 && (dm41_test_2>=ranger_dm41_2->min - 1.0e-10 &&  dm41_test_2<=ranger_dm41_2->max + 1.0e-10)){
    idm41 = ranger_dm41->nsteps + ranger_dm41_2->findIndex(dm41);
  }else{
    cout << "ERROR!! requested dm41=" << dm41 << " in OscProbTable is outside of grid's range ["
         << ranger_dm41->min << "," << ranger_dm41->max << "] or ["
         << ranger_dm41_2->min << "," << ranger_dm41_2->max << "]... aborting!" << endl;
    cout << dm41_test << " " << dm41_test_2 << endl;
    exit(0);
  }
  
  double sin22t12 = 0.857;//eV2

  if (s22t14 != s22t14){
    cout << "in predictor, s22t14 is " << s22t14 << endl;
  }

  
  // protection for negative sin22theta values
  if (s22t14 < 0)s22t14 = 0;
  if (s22t13 < 0)s22t13 = 0;
  
  double theta14=TMath::ASin(sqrt(s22t14))*0.5;
  double theta13=TMath::ASin(sqrt(s22t13))*0.5;
  double theta12=TMath::ASin(sqrt(sin22t12))*0.5;


  //4-neutrino flavor calculation according to Eq. 6 of DocDB-9296-v1
  double Ue1=cos(theta14)*cos(theta13)*cos(theta12);
  double Ue2=cos(theta14)*cos(theta13)*sin(theta12);
  double Ue3=cos(theta14)*sin(theta13);
  double Ue4=sin(theta14);
  
  double term_fac[num_terms] = {
    1,
    4*pow(Ue1,2)*pow(Ue2,2),
    4*pow(Ue1,2)*pow(Ue3,2),
    4*pow(Ue2,2)*pow(Ue3,2),
    4*pow(Ue1,2)*pow(Ue4,2),
    4*pow(Ue2,2)*pow(Ue4,2),
    4*pow(Ue3,2)*pow(Ue4,2)
  };
  
  
  //-->do extrapolation using stored factors
  for(int istage = 0; istage<Nstage;istage++){
    for(int idet_far=4;idet_far<Ndetectors;++idet_far){
      for(int idet_near=0;idet_near<4;++idet_near){
        pred->PredEvtsSpec[istage][idet_far][idet_near]->Reset();
        for(int ibin=1;ibin<=pred->PredEvtsSpec[istage][idet_far][idet_near]->GetXaxis()->GetNbins();++ibin){
          if (m_method == 0){
            Double_t fac_far = 0;
            Double_t fac_near = 0;
            for (Int_t icore = 0; icore < num_cores; icore++){
              for (Int_t iterm = 0; iterm < num_terms; iterm++){
                fac_far += term_fac[iterm]*OscProbTableLookup[istage][idet_far][idmee][idm41][ibin-1][iterm];
                fac_near += term_fac[iterm]*OscProbTableLookup[istage][idet_near][idmee][idm41][ibin-1][iterm];
              }
            }
            double nevts = pred->tdper[istage].CorrEvtsSpec[idet_near]->GetBinContent(ibin);
            double tofill = nevts*fac_far/fac_near;;
            pred->PredEvtsSpec[istage][idet_far][idet_near]->SetBinContent(ibin,tofill);
          }else if (m_method == 1){
            Double_t fac_extrap = 0;
            Double_t norm = 0;
            for(int ibin_enu=1;ibin_enu<=pred->GetNEnuBins();++ibin_enu){
              Double_t fac_far = 0;
              Double_t fac_near = 0;
              Double_t fac_near_null =  term_fac[0]*OscProbTableLookup[istage][idet_near][idmee][idm41][ibin_enu-1][0];
              for (Int_t iterm = 0; iterm < num_terms; iterm++){
                fac_far += term_fac[iterm]*OscProbTableLookup[istage][idet_far][idmee][idm41][ibin_enu-1][iterm];
                fac_near += term_fac[iterm]*OscProbTableLookup[istage][idet_near][idmee][idm41][ibin_enu-1][iterm];
              }
              if (fac_near > 0){
                // Note by Y. Nakajima on Dec. 3, 2013
                // This part attempt to modify evis to enu conversion matrix badsed on oscillation probability
                // In order to make this work, M_evis_to_enu has to be generated assuming no oscillation at all.
                // Then, modified matrix, M'_evis_to_enu is calculated by:
                //              M' = M*fac_near/fac_near_null / norm,
                // where,
                //            norm = Sum_enubin (M * fac_near/fac_near_null)
                // This normalization factor 'norm' is introduced in ordeer to keep Sum_enubin (M') to unity.
                // Then, the extrapolation factor can be calculated as:
                //      fac_extrap = Sum_enubin (M' * fac_far/fac_near)
                //                 = Sum_enubin (M * fac_far/fac_near_null)

                fac_extrap += M_evis_to_enu[ibin-1][ibin_enu-1]*fac_far/fac_near_null;
                norm += M_evis_to_enu[ibin-1][ibin_enu-1]*fac_near/fac_near_null;
              }else{
                cout << "negative survival probability!!!!" << endl;
                cout << ibin_enu << " " << fac_near << " " << fac_far
                     << " " << s22t14 << " " << s22t13 <<  " " << dm41<< " " << dmee << endl;

              }

            }
            double nevts = pred->tdper[istage].CorrEvtsSpec[idet_near]->GetBinContent(ibin);
            if (nevts < 0){ // avoid negative number of events in prediction
              nevts = 0;
            }
            pred->PredEvtsSpec[istage][idet_far][idet_near]->SetBinContent(ibin,nevts*fac_extrap/norm);
          }
        }//ibin

        pred->predper->SetPred(istage,idet_far,idet_near,(TH1F*)pred->PredEvtsSpec[istage][idet_far][idet_near]);

      }//idet_near
    }//idet_far
  }//istage
  
  if(print){
    cout << "-----------------------------------------------------------" << endl;
    cout << "------------------- Predictions ---------------------------" << endl;
    
    for(int ii=0;ii<Nstage;++ii){
      //    predset[ii] = *(MakePrediction(sin2theta13,dm2,ii,tdper[ii]));
      cout << "Period #" << ii << endl;
      pred->predper->PrintToScreen(ii);
      
    }
  }

};//MakeAllPredictionsQuick

void OscProbTable::MakeOneSuperPredictionQuick(double sin2theta13, double dm2, Double_t s22t14, Double_t dm2_41, bool print){

  this->MakeAllPredictionsQuick(sin2theta13,dm2,s22t14,dm2_41, print);

  //PredSet *superpred = new PredSet();  
  
  if(print){
    pred->CombineData();
    pred->CombinePredictions();
    cout << "=================================================================" << endl;
    cout << "================= Super-prediction ==============================" << endl;
    pred->superpred->PrintToScreen(0);
    cout << "================= Combined Data =================================" << endl;
    for(int idet=0;idet<Ndetectors;++idet){
      cout << "AD" << idet+1 << ": " << pred->CombCorrEvtsSpec[idet]->Integral() << endl;//tmp
    }
  }

};

double OscProbTable::CalculateChi2CovQuick(Double_t sin22t13, Double_t dm2_ee,Double_t sin22t14, Double_t dm2_41, Int_t stage){

  this->MakeOneSuperPredictionQuick(sin22t13,dm2_ee,sin22t14,dm2_41,false);
  return pred->CalculateChi2Cov(sin22t13);

}//CalculateChi2Cov
