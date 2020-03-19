#include "TimePeriodData.h"
TimePeriodData::TimePeriodData(){
  for(int idet1=0;idet1<Ndetectors;++idet1){
    ObsEvts[idet1]=0;
    ErrEvts[idet1]=0;
    CorrEvts[idet1]=0;
    BgEvts[idet1]=0;
    AccEvts[idet1]=0;
    Li9Evts[idet1]=0;
    AmcEvts[idet1]=0;
    AlnEvts[idet1]=0;
    FnEvts[idet1]=0;
    AccErr[idet1]=0;
    Li9Err[idet1]=0;
    AmcErr[idet1]=0;
    AlnErr[idet1]=0;
    FnErr[idet1]=0;
    CorrBgEvts[idet1]=0;
    ErrBgEvts[idet1]=0;
    MuonVetoEff[idet1]=0;
    Livetime[idet1]=0;
    DMCEff[idet1]=0;
    TargetMass[idet1]=0;
    BgEvtsLiv[idet1]=0;
  }
  firstcorrection=true;
  nwarnings=0;
  warlimit=5;
};
TimePeriodData::~TimePeriodData(){};

void TimePeriodData::CorrectEvts(bool print){

  //Calculate statistical errors
  for(int idet=0;idet<Ndetectors;++idet){
    ErrEvts[idet]=sqrt(ObsEvts[idet]);
  }

  //Correct backgrounds for livetime and substract

  //cout << "Backgrounds after livetime " << endl;
  for(int idet=0;idet<Ndetectors;++idet){
    BgEvtsLiv[idet]=BgEvts[idet]*Livetime[idet];
    ObsEvts[idet]-=BgEvtsLiv[idet];

    //-->note: scaling to live-days because below re-scale down to 1-liveday
    AccEvts[idet]*=Livetime[idet];
    AmcEvts[idet]*=Livetime[idet];
    AlnEvts[idet]*=Livetime[idet];
    Li9Evts[idet]*=Livetime[idet];
    FnEvts[idet]*=Livetime[idet];
    AccErr[idet]*=Livetime[idet];
    AmcErr[idet]*=Livetime[idet];
    AlnErr[idet]*=Livetime[idet];
    Li9Err[idet]*=Livetime[idet];
    FnErr[idet]*=Livetime[idet];
  }

  //Scale back to 1 day of livetime and correct for muon efficiencies
  //cout << "Observed/Corrected events +++++++++++++++++++++++++++++++++" << endl;//tmp
  for(int idet=0;idet<Ndetectors;++idet){
    //scale to 1 day and mass of AD1

    float factor=0;
    if (MuonVetoEff[idet]*DMCEff[idet]*Livetime[idet]*TargetMass[idet] > 0)
      factor = 1./MuonVetoEff[idet]
        *1./DMCEff[idet]
        *1./Livetime[idet]
        *TargetMass[0]*1./TargetMass[idet];

    CorrEvts[idet]=(ObsEvts[idet])*factor;
    AccEvts[idet]*=factor;
    AmcEvts[idet]*=factor;
    Li9Evts[idet]*=factor;
    FnEvts[idet]*=factor;
    AlnEvts[idet]*=factor;

    AccErr[idet]*=factor;
    AmcErr[idet]*=factor;
    Li9Err[idet]*=factor;
    FnErr[idet]*=factor;
    AlnErr[idet]*=factor;

    ErrEvts[idet]=ErrEvts[idet]*factor;
    CorrBgEvts[idet]=BgEvtsLiv[idet]*factor;
    ErrBgEvts[idet] = sqrt(BgEvtsLiv[idet])*factor;
    if(print==true){
      cout << "AD" << idet+1 << ": " << CorrEvts[idet] << " +/- " << ErrEvts[idet] << "(stat)" << endl;//tmp
    }
  }
}//end of CorrectEvts

void TimePeriodData::CorrectSpec(bool print){

  //Set stat errors
  for(int idet=0;idet<Ndetectors;++idet){
    //cout << "ObsEvtsSpec[idet]: " << ObsEvtsSpec[idet]->Integral() << endl;
    for(int ibin=1;ibin<=ObsEvtsSpec[idet]->GetXaxis()->GetNbins();++ibin){
      ObsEvtsSpec[idet]->SetBinError(ibin,sqrt(ObsEvtsSpec[idet]->GetBinContent(ibin)));
    }//bins
  }//dets

  //Correct backgrounds for livetime and subtract from obs events
  Char_t name[1024];
  for(int idet=0;idet<Ndetectors;++idet){
    float factor=0;

    if (MuonVetoEff[idet]*DMCEff[idet]*Livetime[idet]*TargetMass[idet] > 0)
      factor = 1./MuonVetoEff[idet]
        *1./DMCEff[idet]
        *1./Livetime[idet]
        *TargetMass[0]*1./TargetMass[idet];

    //copy ObsEvtsSpec into CorrEvtsSpec without leaking memory (although this is a small memory leak, since it's only once per fake experiment)
    sprintf(name,"CorrEvtsSpec_AD%i",idet);
    if(firstcorrection){
      CorrEvtsSpec[idet]= (TH1F*)ObsEvtsSpec[idet]->Clone(name);
    } else {
      CorrEvtsSpec[idet]->Reset();//<-- just a precaution
      for(int ibin=1;ibin<=ObsEvtsSpec[idet]->GetXaxis()->GetNbins();++ibin){
        CorrEvtsSpec[idet]->SetBinContent(ibin,ObsEvtsSpec[idet]->GetBinContent(ibin));
        CorrEvtsSpec[idet]->SetBinError(ibin,ObsEvtsSpec[idet]->GetBinError(ibin));
      }
    }

    if (factor>0)
      CorrBgEvtsSpec[idet]->Scale(1/factor);
    else
      CorrBgEvtsSpec[idet]->Scale(0);

    CorrEvtsSpec[idet]->Add(CorrBgEvtsSpec[idet],-1);
    CorrEvtsSpec[idet]->Scale(factor);
    CorrBgEvtsSpec[idet]->Scale(factor);

  }//idet

  //Ensuring consistency
  if(nwarnings<warlimit){
    ++nwarnings;
    cout << "Ensuring consistency between numbers from spectra and text file" << endl;
    for(int idet=0;idet<Ndetectors;++idet){
      double intspec = CorrEvtsSpec[idet]->Integral();
      double intnum = CorrEvts[idet];
      if(print==true){
        cout << "AD" << idet+1 << ": " << intspec << " (spec); " << intnum << " (text)" << endl;//tmp
      }
      if(TMath::Abs(intspec - intnum) > 0.001 * intspec){
        cout << "Ay caramba!!! Discrepancy in number of IBD corrected events between spectra and text files; should abort..." << endl;
        //cout << "Enter one number and press enter to continue, although you should really hit Ctrl+C" << endl;
        //Note: abilitate lines below when have 9Li and other backgrounds
        //int catastrophe;
        //cin >> catastrophe;

      }
    }//idet

  }else if(nwarnings==warlimit){
    cout << "WILL NOT PERFORM ANY MORE CONSISTENCY CHECKS AND WILL NOT PRINT ANY MORE WARNINGS" << endl;
    ++nwarnings;
  }

  firstcorrection=false;
}

void TimePeriodData::PrintToScreen(){

  //for(int istage=0;istage<Nstage;++istage){
  // cout << "---Stage: " << istage+1 << " ---" << endl;

  cout << "Observed events: ";
  for(int idet=0;idet<Ndetectors;++idet){
    cout << ObsEvts[idet]+BgEvtsLiv[idet];

    if(idet == Ndetectors-1) cout << endl;
    else cout << ",";
  }

  cout << "Muon efficiencies: ";
  for(int idet=0;idet<Ndetectors;++idet){
    cout << MuonVetoEff[idet];

    if(idet == Ndetectors-1) cout << endl;
    else cout << ",";
  }

  cout << "DMC efficiencies: ";
  for(int idet=0;idet<Ndetectors;++idet){
    cout << DMCEff[idet];

    if(idet == Ndetectors-1) cout << endl;
    else cout << ",";
  }

  cout << "Livetime: ";
  for(int idet=0;idet<Ndetectors;++idet){
    cout << Livetime[idet];

    if(idet == Ndetectors-1) cout << endl;
    else cout << ",";
  }

  cout << "Backgrounds: ";
  for(int idet=0;idet<Ndetectors;++idet){
    cout << BgEvts[idet];

    if(idet == Ndetectors-1) cout << endl;
    else cout << ",";
  }

  cout << "Target Masses: ";
  for(int idet=0;idet<Ndetectors;++idet){
    cout << TargetMass[idet];

    if(idet == Ndetectors-1) cout << endl;
    else cout << ",";
  }

  cout << "Corrected Evts: ";
  for(int idet=0;idet<Ndetectors;++idet){
    cout << CorrEvts[idet];

    if(idet == Ndetectors-1) cout << endl;
    else cout << ",";
  }

  cout << "--------------------" << endl;
}
