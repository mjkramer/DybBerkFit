#include <iostream>
#include <fstream>
#include <sstream>
#include "Predictor.h"
#include "TMatrixD.h"
#include "TF1.h"

/*
  June 19: Yasu
  Fixed bug in the livetime weighting
  October 2012: Pedro
  Changed SuperMatrix by SuperHists approach
*/

Predictor::Predictor(){

  superpred = new PredSet();
  predper = new PredSet();

  //TimePeriodData's
  for(int ii=0;ii<Nstage;++ii){
    tdper[ii]=TimePeriodData();
  }

  //Variables
  BgMatrixScaled = false;
  isMC=false;

  FirstToySample = true;
  FirstMakePrediction = true;
  FirstMakeSuperPrediction = true;

  CalculateMode1Coeff = true;
  CalculateNearEHCoeff = true;
  CalculateStageCoeff = true;
  RecalculateCovMatrix = true;

  for (Int_t istage=0;istage<Nstage;istage++)
    stat_factor[istage] = 1;

  days_6ad = 189.88;
  days_8ad =  922.92;
  days_7ad =  200.;  //I do not know this number exactly

  days[0] = days_6ad;
  days[1] = days_8ad;
  days[2] = days_7ad;

  extra_days_8ad = 0;

  combine_mode = 1;
  // combine_mode = 0;

  // combine_mode = 0: No scaling, 16*Nstage*n_evis_bins dimensions
  // combine_mode = 1: Default scaling with  2*Nstage*n_evis_bins dimensions
  // combine_mode = 2: Scaling with  4*Nstage*n_evis_bins dimensions
  // combine_mode = 3: Ultimate scaling with  Nstage*n_evis_bins dimensions, which fit the sum of three far detectors against the weigthed mean of 3 near site predictions

}

Predictor::~Predictor(){};

void Predictor::SetStatFactor(Double_t fac){

  stat_factor[Nstage-1] = fac;

}

void Predictor::SetStage(Int_t i){

  stage = i;

}

void Predictor::EnterFluxCalculator(FluxCalculator *fluxcalc_in){

  fluxcalc = fluxcalc_in;

}

void Predictor::SetEvisBins(Int_t n, Double_t* bins, Int_t rebin_fac){
  n_evis_bins = n;
  for (Int_t i = 0; i < n+1; i++){
    evis_bins[i] = bins[i];
    evis_bins_rebin[i] = 0;
  }
  for (Int_t istage = 0; istage < Nstage; istage ++){
    for (Int_t idet = 0; idet < Ndetectors; idet ++){
      for (Int_t jdet = 0; jdet < Ndetectors; jdet ++){
        PredEvtsSpec[istage][idet][jdet] = new TH1F(Form("h_evis_Stage%i_AD%i_AD%i",istage,idet,jdet),
                                                    Form("h_evis_Stage%i_AD%i_AD%i",istage,idet,jdet),
                                                    n_evis_bins,evis_bins);
      }
    }
  }

  evis_bins_rebin[0] = evis_bins[0];

  n_evis_bins_rebin = (Int_t)(n_evis_bins-1)/rebin_fac + 1;
  cout << "Rebinning from " << n_evis_bins << " to " << n_evis_bins_rebin << endl;
  for (Int_t i = 0; i < n_evis_bins; i++){
    Int_t i_rebin =(Int_t)i/rebin_fac;
    evis_rebin_map[i] = i_rebin;
    cout << "rebin " << i << " " << evis_rebin_map[i] << endl;
    if (evis_bins_rebin[i_rebin+1] < evis_bins[i+1]){
      evis_bins_rebin[i_rebin+1] = evis_bins[i+1];
    }
  }
  cout << "Bin boundaries";
  for (Int_t i = 0; i < n_evis_bins_rebin+1; i++){
    cout << "\t" << evis_bins_rebin[i];
  }
  cout << endl;

}

Double_t * Predictor::GetRebinnedEvisBins(){ // return bins after rebinning
  return evis_bins_rebin;
}

Int_t Predictor::GetNumRebinnedEvisBins(){ // return bins after rebinning
  return n_evis_bins_rebin;
}


void Predictor::SetEnuBins(Int_t n, Double_t* bins){
  n_enu_bins = n;
  cout << "----------- enu bins-------------- " << endl;
  for (Int_t i = 0; i < n+1; i++){
    enu_bins[i] = bins[i];
    cout << "\t" << enu_bins[i];
  }
  cout << endl;

  for (Int_t istage = 0; istage < Nstage; istage++ ){
    for (Int_t idet = 0; idet < Ndetectors; idet++ ){
      for (Int_t ie = 0; ie < max_n_evis_bins; ie++ ){
        CorrEvtsTrueSpec[istage][idet][ie] = new TH1F(Form("h_true_enu_stage%d_det%d_e%d",istage,idet,ie),
                                                      Form("h_true_enu_stage%d_det%d_e%d",istage,idet,ie),
                                                      n_enu_bins,enu_bins);
      }
    }
  }

}

void Predictor::LoadMainData(const Char_t *mainmatrixname){

  cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "+++++++++++++++++++ Reading Main File +++++++++++++++++" << endl;

  string line;
  int linenum=0;//<---caution: only increments for lines that do not begin with #
  int istage=0;
  ifstream mainfile(mainmatrixname);
  while(!mainfile.eof()){

    getline(mainfile,line);
    string firstchar=line.substr(0,1);
    //cout << "The first char " << firstchar << endl;//tmp
    if(firstchar=="#") continue;//<-- ignore lines with comments

    //Special numbers
    if(linenum==0) {
      istage=atoi(line.c_str());
      cout << "Stage: " << istage << endl;
    }
    if(linenum==1){
      if(atoi(firstchar.c_str())==0) isMC=true;
      if(atoi(firstchar.c_str())==1) isMC=false;
      cout << "Simflag: " << isMC << " (1-->MC, 0-->Data)" << endl;
    }
    // The meat of the file <----------------
    if(linenum>2) {;//<--the "meat" begins at line=3
      cout << "reading " << line << endl;
      istringstream iss(line);
      int row=0;
      int column=0;
      double readvals[Ndetectors]={0};
      while(iss){
        string sub; iss >> sub;
        //if(column==0) istage=atoi(sub.c_str());
        if(column==1) row=atoi(sub.c_str());
        if(column>1 && sub!="") readvals[column-2]=atof(sub.c_str());
        column+=1;
      }//looping over columns

      //Fill stuff that was found ++++++++++++++++++
      //-->dates
      if(row==0) continue;//<--we don't care
      //-->obs events
      if(row==1){
        for(int ii=0;ii<Ndetectors;ii++){
          tdper[istage-1].ObsEvts[ii]=readvals[ii];
        }
      }
      //-->livetimes
      if(row==2){
        for(int ii=0;ii<Ndetectors;ii++){
          tdper[istage-1].Livetime[ii]=readvals[ii];
        }
      }
      //-->muon efficiencies
      if(row==3){
        for(int ii=0;ii<Ndetectors;ii++){
          tdper[istage-1].MuonVetoEff[ii]=readvals[ii];
        }
      }
      //-->dmc efficiencies
      if(row==4){
        for(int ii=0;ii<Ndetectors;ii++){
          tdper[istage-1].DMCEff[ii]=readvals[ii];
        }
      }
      //-->target masses
      if(row==8){
        for(int ii=0;ii<Ndetectors;ii++){
          tdper[istage-1].TargetMass[ii]=readvals[ii];
        }
      }
      //-->bg events
      if(row==9)continue;//<-- we don't care about total background

      //-->bg systematics
      if(row==10) continue;//<-- we don't care about background systematics
      //-->acc bg
      if(row==11){
        for(int ii=0;ii<Ndetectors;ii++){
          tdper[istage-1].AccEvts[ii]=readvals[ii];
          //in new input format have to correct for efficiencies here
          tdper[istage-1].AccEvts[ii]*=tdper[istage-1].DMCEff[ii]
            *tdper[istage-1].MuonVetoEff[ii];

        }
      }
      if(row==12){
        for(int ii=0;ii<Ndetectors;ii++){
          tdper[istage-1].AccErr[ii]=readvals[ii];
          //in new input format have to correct for efficiencies here
          tdper[istage-1].AccErr[ii]*=tdper[istage-1].DMCEff[ii]
            *tdper[istage-1].MuonVetoEff[ii];
        }
      }

      //-->li9 bg
      if(row==13){
        for(int ii=0;ii<Ndetectors;ii++){
          tdper[istage-1].Li9Evts[ii]=readvals[ii];
          //in new input format have to correct for efficiencies here
          tdper[istage-1].Li9Evts[ii]*=tdper[istage-1].DMCEff[ii]
            *tdper[istage-1].MuonVetoEff[ii];
        }
      }
      if(row==14){
        for(int ii=0;ii<Ndetectors;ii++){
          tdper[istage-1].Li9Err[ii]=readvals[ii];
          //in new input format have to correct for efficiencies here
          tdper[istage-1].Li9Err[ii]*=tdper[istage-1].DMCEff[ii]
            *tdper[istage-1].MuonVetoEff[ii];
        }
      }
      //-->fast-n bg
      if(row==15){
        for(int ii=0;ii<Ndetectors;ii++){
          tdper[istage-1].FnEvts[ii]=readvals[ii];
          //in new input format have to correct for efficiencies here
          tdper[istage-1].FnEvts[ii]*=tdper[istage-1].DMCEff[ii]
            *tdper[istage-1].MuonVetoEff[ii];
        }
      }
      if(row==16){
        for(int ii=0;ii<Ndetectors;ii++){
          tdper[istage-1].FnErr[ii]=readvals[ii];
          //in new input format have to correct for efficiencies here
          tdper[istage-1].FnErr[ii]*=tdper[istage-1].DMCEff[ii]
            *tdper[istage-1].MuonVetoEff[ii];
        }
      }
      //-->amc bg
      if(row==17){
        for(int ii=0;ii<Ndetectors;ii++){
          tdper[istage-1].AmcEvts[ii]=readvals[ii];
          //in new input format have to correct for efficiencies here
          tdper[istage-1].AmcEvts[ii]*=tdper[istage-1].DMCEff[ii]
            *tdper[istage-1].MuonVetoEff[ii];
        }
      }
      if(row==18){
        for(int ii=0;ii<Ndetectors;ii++){
          tdper[istage-1].AmcErr[ii]=readvals[ii];
          //in new input format have to correct for efficiencies here
          tdper[istage-1].AmcErr[ii]*=tdper[istage-1].DMCEff[ii]
            *tdper[istage-1].MuonVetoEff[ii];
        }
      }
      //-->aln bg
      if(row==19){
        for(int ii=0;ii<Ndetectors;ii++){
          tdper[istage-1].AlnEvts[ii]=readvals[ii];
          //in new input format have to correct for efficiencies here
          tdper[istage-1].AlnEvts[ii]*=tdper[istage-1].DMCEff[ii]
            *tdper[istage-1].MuonVetoEff[ii];
        }
      }
      if(row==20){
        for(int ii=0;ii<Ndetectors;ii++){
          tdper[istage-1].AlnErr[ii]=readvals[ii];
          //in new input format have to correct for efficiencies here
          tdper[istage-1].AlnErr[ii]*=tdper[istage-1].DMCEff[ii]
            *tdper[istage-1].MuonVetoEff[ii];
        }
      }

    }//only lines >1
    ++linenum;
  }//Reading main file

  //Add total background here
  for(int ii=0;ii<Ndetectors;ii++){
    tdper[istage-1].BgEvts[ii]=0;
    tdper[istage-1].BgEvts[ii]+=tdper[istage-1].AccEvts[ii];
    tdper[istage-1].BgEvts[ii]+=tdper[istage-1].AmcEvts[ii];
    tdper[istage-1].BgEvts[ii]+=tdper[istage-1].Li9Evts[ii];
    tdper[istage-1].BgEvts[ii]+=tdper[istage-1].FnEvts[ii];
    tdper[istage-1].BgEvts[ii]+=tdper[istage-1].AlnEvts[ii];
  }

  // for(int ii=0;ii<Nstage;++ii){
  cout << "Stage #" << istage << endl;
  tdper[istage-1].CorrectEvts(true);
  tdper[istage-1].PrintToScreen();//tmp
  //}
}//end of LoadMainMatrix

void Predictor::LoadPredictedIBD(const Char_t *nibdname){
  cout << "Loading # ibd ..." << endl;

  TFile *infile = new TFile(nibdname,"READ");

  Char_t name[1024];

  for(int istage=0;istage<Nstage;++istage){
    cout << "istage: " << istage << endl;
    for(int idet=0;idet<Ndetectors;++idet){

      sprintf(name,"h_nominal_stage%i_ad%i",istage+1,idet+1);
      TH1F *htemp = (TH1F*)infile->Get(name);
      nIBD[istage][idet] = htemp->Integral();
      cout << "idet: " << idet << ", nIBD: " << nIBD[istage][idet] << endl;
      htemp->Delete();
    }
  }
  cout << "done loading # ibd!" << endl;
}

void Predictor::LoadIBDSpec(TString *ibdspecname){

  cout << "Loading ibd spectra..." << endl;
  Char_t name[1024];
  Char_t name2[1024];
  //++ IBDs
  TFile *infilespec[Nstage];
  for(int istage=0;istage<Nstage;++istage){
    infilespec[istage] = new TFile(ibdspecname[istage].Data(),"READ");
    cout<<"File is "<<ibdspecname[istage].Data()<<endl;

    for(int idet=0;idet<Ndetectors;++idet){

      sprintf(name,"h_ibd_eprompt_inclusive_eh%i_ad%d",detConfigEH[idet],detConfigAD[idet]);
      sprintf(name2,"h_ibd_eprompt_inclusive_stage%i_eh%i_ad%d",istage+1,detConfigEH[idet],detConfigAD[idet]);

      //cout << "tdper[iweek].ObsEvtsSpec[idet]->Integral(): " << tdper[0].ObsEvtsSpec[idet]->Integral() << endl;
      tdper[istage].ObsEvtsSpec[idet] = (TH1F*)infilespec[istage]->Get(name)->Clone(name2);

    }
  }
  cout << "done loading ibd spectra!" << endl;
}

Int_t Predictor::LoadToyIBDSpec(const Char_t* filename){

  cout << "Loading toy ibd spectra..." << endl;

  m_toyinfilespec = new TFile(filename);
  m_tr_toy  = (TTree*)m_toyinfilespec->Get("tr");
  m_tr_toy->Print();

  //++ IBDs
  for(int istage=0;istage<Nstage;++istage){
    for(int idet=0;idet<Ndetectors;++idet){
      tdper[istage].ObsEvtsSpec[idet] = 0;
      m_tr_toy->SetBranchAddress(Form("h_stage%i_ad%i",istage+1,idet+1),&(tdper[istage].ObsEvtsSpec[idet]));
    }
  }
  cout << "done loading ibd spectra file!" << endl;

  isMC = true;
  return m_tr_toy->GetEntries();
}


void Predictor::LoadToyMCNominalSpec(){

  m_toyinfilespec->cd();
  for(int ii=0;ii<Nstage;++ii){
    cout << "Period #" << ii+1 << endl;
    for(int idet=0;idet<Ndetectors;++idet){
      tdper[ii].ObsEvtsSpec[idet] = (TH1F*)m_toyinfilespec->Get(Form("h_nominal_stage%i_ad%d",ii+1,idet+1));
    }
    tdper[ii].CorrectSpec(true);
    //tdper[ii].PrintToScreen();//tmp
  }

}

void Predictor::LoadToyMCEntry(Int_t i, bool correct){

  m_tr_toy->GetEntry(i);
  cout << "Loading toy MC entry : " << i << endl;
  for(int ii=0;ii<Nstage;++ii){
    if(correct==true) tdper[ii].CorrectSpec(true);
    //tdper[ii].PrintToScreen();//tmp
  }
  RecalculateCovMatrix = true; // Force to recalculate

}

void Predictor::LoadBgSpec(TString *accspecname,
                           const Char_t *li9specname,
                           const Char_t *amcspecname,
                           const Char_t *fnspecname,
                           const Char_t *alnspecname){

  cout << "Loading bg spectra..." << endl;
  Char_t name[1024];

  //(accidentals)
  TFile *accspec[Nstage];

  for(int istage=0;istage<Nstage;++istage){
    accspec[istage]= new TFile(accspecname[istage].Data(),"READ");

    for(int idet=0;idet<Ndetectors;++idet){

      sprintf(name,"h_accidental_eprompt_inclusive_eh%i_ad%i",detConfigEH[idet],detConfigAD[idet]);

      Char_t name2[1024];
      sprintf(name2,"CorrAccEvtsSpec_stage%i_eh%i_ad%d",istage+1,detConfigEH[idet],detConfigAD[idet]);



      tdper[istage].CorrAccEvtsSpec[idet] = (TH1F*)accspec[istage]->Get(name)->Clone(name2);

      for (Int_t ibin = 0; ibin < tdper[istage].CorrAccEvtsSpec[idet]->GetNbinsX();ibin++){
        tdper[istage].CorrAccEvtsSpec[idet]->SetBinError(ibin+1,0);
      }

    }
    //cout << "tdper[iweek].ObsEvtsSpec[idet]->Integral(): " << tdper[0].ObsEvtsSpec[idet]->Integral() << endl;
  }

  for(int istage=0;istage<Nstage;++istage){
    for(int idet=0;idet<Ndetectors;++idet){
      if (tdper[istage].CorrAccEvtsSpec[idet]->Integral() == 0)
        tdper[istage].CorrAccEvtsSpec[idet]->Scale(0);
      else
        tdper[istage].CorrAccEvtsSpec[idet]->Scale(tdper[istage].AccEvts[idet]*1./tdper[istage].CorrAccEvtsSpec[idet]->Integral());//<---Scale this histogram to expected number of events based on input file
    }
  }
  //accspec->Close();
  cout << "--> loaded accidental spectra" << endl;

  //(li9/he8)
  TFile *li9spec = new TFile(li9specname,"READ");
  for(int istage=0;istage<Nstage;++istage){
    for(int idet=0;idet<Ndetectors;++idet){
      sprintf(name,"h_nominal");
      TH1F *htemp = (TH1F*)li9spec->Get(name);
      sprintf(name,"CorrLi9EvtsSpec_stage%d_ad%d",istage,idet);
      tdper[istage].CorrLi9EvtsSpec[idet] = (TH1F*)tdper[istage].ObsEvtsSpec[idet]->Clone(name);
      tdper[istage].CorrLi9EvtsSpec[idet]->Reset();
      //fill into destination histogram
      for(Int_t ibin = 0; ibin < htemp->GetNbinsX();ibin++){
        tdper[istage].CorrLi9EvtsSpec[idet]->Fill(htemp->GetBinCenter(ibin+1),htemp->GetBinContent(ibin+1));
      }
      delete htemp;

      //set bin errors to zero
      for (Int_t ibin = 0; ibin < tdper[istage].CorrLi9EvtsSpec[idet]->GetNbinsX();ibin++){
        tdper[istage].CorrLi9EvtsSpec[idet]->SetBinError(ibin+1,0);
      }
      //cout << "The scaling factor is " << tdper[istage].Li9Evts[idet] << "," << tdper[istage].CorrLi9EvtsSpec[idet]->Integral()  << endl;//tmp

      //scale
      if (tdper[istage].CorrLi9EvtsSpec[idet]->Integral() == 0)
        tdper[istage].CorrLi9EvtsSpec[idet]->Scale(0);
      else
        tdper[istage].CorrLi9EvtsSpec[idet]->Scale(tdper[istage].Li9Evts[idet]*1./tdper[istage].CorrLi9EvtsSpec[idet]->Integral());//<---Scale this histogram to expected number of events based on input file
    }
  }
  //li9spec->Close();
  cout << "--> loaded li9/he8 spectra" << endl;

  //(fast-n)
  TFile *fnspec = new TFile(fnspecname,"READ");
  for(int istage=0;istage<Nstage;++istage){
    for(int idet=0;idet<Ndetectors;++idet){

      sprintf(name,"CorrFnEvtsSpec_stagek%d_ad%d",istage,idet);
      Char_t nameFn[1024];
      sprintf(nameFn,"h_%dAD_fn_fine",idet+1);
      tdper[istage].CorrFnEvtsSpec[idet] = (TH1F*)tdper[istage].ObsEvtsSpec[idet]->Clone(name);
      tdper[istage].CorrFnEvtsSpec[idet]->Reset();

      TH1F *htemp = (TH1F*)fnspec->Get(nameFn)->Clone(name);
      //fill into destination histogram
      for(Int_t ibin = 0; ibin < htemp->GetNbinsX();ibin++){
        tdper[istage].CorrFnEvtsSpec[idet]->Fill(htemp->GetBinCenter(ibin+1),htemp->GetBinContent(ibin+1));
      }
      delete htemp;

      for (Int_t ibin = 0; ibin < tdper[istage].CorrFnEvtsSpec[idet]->GetNbinsX();ibin++){
        tdper[istage].CorrFnEvtsSpec[idet]->SetBinError(ibin+1,0);
      }
      if (tdper[istage].CorrFnEvtsSpec[idet]->Integral() == 0)
        tdper[istage].CorrFnEvtsSpec[idet]->Scale(0);
      else
        tdper[istage].CorrFnEvtsSpec[idet]->Scale(tdper[istage].FnEvts[idet]*1./tdper[istage].CorrFnEvtsSpec[idet]->Integral());//<---Scale this histogram to expected number of events based on input file
    }
  }
  //fnspec->Close();
  cout << "--> loaded fast-n spectra" << endl;

  //(amc)
  TFile *amcspec = new TFile(amcspecname,"READ");
  for(int istage=0;istage<Nstage;++istage){
    for(int idet=0;idet<Ndetectors;++idet){

      TF1 *amcfunc = (TF1*)amcspec->Get("expo");
      sprintf(name,"CorrAmcEvtsSpec_stage%d_ad%d",istage,idet);
      tdper[istage].CorrAmcEvtsSpec[idet] = (TH1F*)tdper[istage].CorrFnEvtsSpec[idet]->Clone(name);
      tdper[istage].CorrAmcEvtsSpec[idet]->Reset();

      //Add function
      //Note: cannot just do Add(amcfunc) as bin sizes are different
      for(int ibin=1;ibin<=tdper[istage].CorrAmcEvtsSpec[idet]->GetXaxis()->GetNbins();++ibin){

        double lowedge=tdper[istage].CorrAmcEvtsSpec[idet]->GetXaxis()->GetBinLowEdge(ibin);
        double upedge=tdper[istage].CorrAmcEvtsSpec[idet]->GetXaxis()->GetBinUpEdge(ibin);
        double intval=amcfunc->Integral(lowedge,upedge);
        tdper[istage].CorrAmcEvtsSpec[idet]->SetBinContent(ibin,intval);

      }

      //old
      //tdper[istage].CorrAmcEvtsSpec[idet] = (TH1F*)amcspec->Get("h_rebin")->Clone(name);

      for (Int_t ibin = 0; ibin < tdper[istage].CorrAmcEvtsSpec[idet]->GetNbinsX();ibin++){
        tdper[istage].CorrAmcEvtsSpec[idet]->SetBinError(ibin+1,0);
      }

      if (tdper[istage].CorrAmcEvtsSpec[idet]->Integral() == 0)
        tdper[istage].CorrAmcEvtsSpec[idet]->Scale(0);
      else
        tdper[istage].CorrAmcEvtsSpec[idet]->Scale(tdper[istage].AmcEvts[idet]*1./tdper[istage].CorrAmcEvtsSpec[idet]->Integral());//<---Scale this histogram to expected number of events per liveday based on input file
    }
  }
  //amcspec->Close();
  cout << "--> loaded AmC spectra" << endl;

  //(alpha-n)
  int AlphaAD[8] = {1,2,3,8,4,5,6,7}; //IHEP AD convention
  TFile *alnspec = new TFile(alnspecname,"READ");
  for(int istage=0;istage<Nstage;++istage){
    for(int idet=0;idet<Ndetectors;++idet){

      sprintf(name,"AD%i;1",AlphaAD[idet]);
      TH1F *htemp = (TH1F*)alnspec->Get(name);
      sprintf(name,"CorrAlnEvtsSpec_stage%d_ad%d",istage,idet);
      tdper[istage].CorrAlnEvtsSpec[idet] = (TH1F*)tdper[istage].ObsEvtsSpec[idet]->Clone(name);
      tdper[istage].CorrAlnEvtsSpec[idet]->Reset();

      //fill into destination histogram
      for(Int_t ibin = 0; ibin < htemp->GetNbinsX();ibin++){
        tdper[istage].CorrAlnEvtsSpec[idet]->Fill(htemp->GetBinCenter(ibin+1),htemp->GetBinContent(ibin+1));
      }

      delete htemp;

      //assign zero errors
      for (Int_t ibin = 0; ibin < tdper[istage].CorrAlnEvtsSpec[idet]->GetNbinsX();ibin++){
        tdper[istage].CorrAlnEvtsSpec[idet]->SetBinError(ibin+1,0);
      }

      //scale
      if (tdper[istage].CorrAlnEvtsSpec[idet]->Integral() == 0)
        tdper[istage].CorrAlnEvtsSpec[idet]->Scale(0);
      else
        tdper[istage].CorrAlnEvtsSpec[idet]->Scale(tdper[istage].AlnEvts[idet]*1./tdper[istage].CorrAlnEvtsSpec[idet]->Integral());//<---Scale this histogram to expected number of events based on input file
    }
  }
  //alnspec->Close();
  cout << "--> loaded alpha-n spectra" << endl;

  cout << "done loading bg spectra!" << endl;

  //create sum of bg's spectra
  for(int istage=0;istage<Nstage;++istage){
    for(int idet=0;idet<Ndetectors;++idet){

      sprintf(name,"CorrBgEvtsSpec_stage%d_ad%d",istage,idet);
      tdper[istage].CorrBgEvtsSpec[idet] = (TH1F*)tdper[istage].CorrAccEvtsSpec[idet]->Clone(name);
      tdper[istage].CorrBgEvtsSpec[idet]->Add(tdper[istage].CorrLi9EvtsSpec[idet],1);
      tdper[istage].CorrBgEvtsSpec[idet]->Add(tdper[istage].CorrAmcEvtsSpec[idet],1);
      tdper[istage].CorrBgEvtsSpec[idet]->Add(tdper[istage].CorrFnEvtsSpec[idet],1);
      tdper[istage].CorrBgEvtsSpec[idet]->Add(tdper[istage].CorrAlnEvtsSpec[idet],1);

      tdper[istage].CorrBgEvtsSpec[idet]->SetDirectory(0);
      tdper[istage].CorrAccEvtsSpec[idet]->SetDirectory(0);
      tdper[istage].CorrLi9EvtsSpec[idet]->SetDirectory(0);
      tdper[istage].CorrAmcEvtsSpec[idet]->SetDirectory(0);
      tdper[istage].CorrFnEvtsSpec[idet]->SetDirectory(0);
      tdper[istage].CorrAlnEvtsSpec[idet]->SetDirectory(0);

    }
  }

  accspec[0]->Close();
  accspec[1]->Close();
  li9spec->Close();
  fnspec->Close();
  amcspec->Close();
  alnspec->Close();
  //cout << "FIXME: need correction code here for data" << endl;
  // apply corrections
  for(int istage=0;istage<Nstage;++istage){
    cout << "Period #" << istage+1 << endl;
    tdper[istage].CorrectSpec(true);
  }
  cout << "done correcting both ibd and bg spectra!" << endl;

}//end of LoadBgSpec

//Extract Bkg Spectrum here
TH1F* Predictor::GetCorrAccEvtsSpec(Int_t istage, Int_t idet){
  TH1F * dummy = (TH1F*)tdper[istage].CorrAccEvtsSpec[idet]->Clone();
  dummy->Scale(tdper[istage].Livetime[idet] * tdper[istage].MuonVetoEff[idet] * tdper[istage].DMCEff[idet]);
  return dummy;
  // return tdper[istage].CorrAccEvtsSpec[idet]->Scale(tdper[istage].Livetime[idet]);
}
TH1F* Predictor::GetCorrLi9EvtsSpec(Int_t istage, Int_t idet){
  TH1F * dummy = (TH1F*)tdper[istage].CorrLi9EvtsSpec[idet]->Clone();
  dummy->Scale(tdper[istage].Livetime[idet] * tdper[istage].MuonVetoEff[idet] * tdper[istage].DMCEff[idet]);
  return dummy;
  //return tdper[istage].CorrLi9EvtsSpec[idet]->Scale(tdper[istage].Livetime[idet]);
}
TH1F* Predictor::GetCorrAmcEvtsSpec(Int_t istage, Int_t idet){
  TH1F * dummy = (TH1F*)tdper[istage].CorrAmcEvtsSpec[idet]->Clone();
  dummy->Scale(tdper[istage].Livetime[idet] * tdper[istage].MuonVetoEff[idet] * tdper[istage].DMCEff[idet]);
  return dummy;
  //return tdper[istage].CorrAmcEvtsSpec[idet]->Scale(tdper[istage].Livetime[idet]);
}
TH1F* Predictor::GetCorrFnEvtsSpec(Int_t istage, Int_t idet){
  TH1F * dummy = (TH1F*)tdper[istage].CorrFnEvtsSpec[idet]->Clone();
  dummy->Scale(tdper[istage].Livetime[idet] * tdper[istage].MuonVetoEff[idet] * tdper[istage].DMCEff[idet]);
  return dummy;
  //return tdper[istage].CorrFnEvtsSpec[idet]->Scale(tdper[istage].Livetime[idet]);
}
TH1F* Predictor::GetCorrAlnEvtsSpec(Int_t istage, Int_t idet){
  TH1F * dummy = (TH1F*)tdper[istage].CorrAlnEvtsSpec[idet]->Clone();
  dummy->Scale(tdper[istage].Livetime[idet] * tdper[istage].MuonVetoEff[idet] * tdper[istage].DMCEff[idet]);
  return dummy;
  //return tdper[istage].CorrAlnEvtsSpec[idet]->Scale(tdper[istage].Livetime[idet]);
}
//Done!

void Predictor::EnterObsEvts(double obsad1,double obsad2, double obsad3, double obsad4, double obsad5, double obsad6, double obsad7, double obsad8,int istage){

  tdper[istage].ObsEvts[0]=obsad1;
  tdper[istage].ObsEvts[1]=obsad2;
  tdper[istage].ObsEvts[2]=obsad3;
  tdper[istage].ObsEvts[3]=obsad4;
  tdper[istage].ObsEvts[4]=obsad5;
  tdper[istage].ObsEvts[5]=obsad6;
  tdper[istage].ObsEvts[6]=obsad7;
  tdper[istage].ObsEvts[7]=obsad8;
  tdper[istage].CorrectEvts(false);

}

TH1F* Predictor::GetCombCorrEvtsSpec(int idet){

  return CombCorrEvtsSpec[idet];

}

TH1F* Predictor::GetCorrEvtsSpec(int istage, int idet){

  //TH1F * dummy = (TH1F*)tdper[istage].CorrEvtsSpec[idet]->Clone();
  //dummy->Scale(tdper[istage].Livetime[idet]);// * tdper[istage].MuonVetoEff[idet] * tdper[istage].DMCEff[idet]);
  //return dummy;
  //return tdper[istage].CorrEvtsSpec[idet]->Scale(tdper[istage].Livetime[idet]);
  return tdper[istage].CorrEvtsSpec[idet];
}

TimePeriodData *Predictor::GetTimePeriodData(int istage){

  return &tdper[istage];

}

void Predictor::MakePrediction(double sin2theta13, double dm2, double s22t14, double dm2_41, int istage, TimePeriodData &tperdat){ //, PredSet *predout){
  for(int idet=0;idet<4;++idet){
    // Convert from Evis to Eneu_true
    // The procedure is rather complex now in oreder to rescale Evis-to-Enu conversion matrix using the oscillation probability.
    std::map<int,TH1F*> mapflux_noosc = fluxcalc->CalculateFluxHistRow(idet,sin2theta13,istage,dm2,s22t14,dm2_41,0);
    // Get the inputs to calculate reactor flux contribution from each core.
    // This information is needed in order to properly rescale the matrix using oscillation probability.
    // For this, we should not assume any oscillation, because evis-to-enu matrix is generated without any oscillation.
    // The final argument of CalculateFluxHistRow, "0",  removes any oscillation effects.
    if (tperdat.CorrEvtsSpec[idet]->GetXaxis()->GetNbins() != n_evis_bins){
      cout << "Error Error Error! Inconsistent number of bins for the spectrum inputs " << tperdat.CorrEvtsSpec[idet]->GetXaxis()->GetNbins() << " vs. " <<  n_evis_bins << endl;
    }
    double factor_noosc[max_n_enu_bins][Ncores];
    for(int ibin_enu=0;ibin_enu< n_enu_bins; ++ibin_enu){
      double factor_norm=0;

      double Etrue=CorrEvtsTrueSpec[istage][idet][0]->GetXaxis()->GetBinCenter(ibin_enu+1);

      for(int icore=0;icore<Ncores;++icore){

        int thebin=mapflux_noosc[icore]->GetXaxis()->FindBin(Etrue);
        factor_noosc[ibin_enu][icore] = mapflux_noosc[icore]->GetBinContent(thebin);
        factor_norm += factor_noosc[ibin_enu][icore];

      }
      for(int icore=0;icore<Ncores;++icore){
        factor_noosc[ibin_enu][icore] /= factor_norm;
      }

    }
    std::map<int,TH1F*> mapextrap = fluxcalc->ExtrapolationFactorRow(idet,-1,sin2theta13,istage,dm2,s22t14,dm2_41);

    for(int ibin_evis=0;ibin_evis < n_evis_bins; ++ibin_evis){
      Double_t norm = 0;
      double oldval = tperdat.CorrEvtsSpec[idet]->GetBinContent(ibin_evis+1);
      if (oldval < 0){
        cout << "Predictor: Negative bin content after background subtraction for detector " << idet
             << ", bin " << ibin_evis << endl;
        cout << "Setting prediction for this bin to be 0" << endl;
        oldval = 0;
      }
      CorrEvtsTrueSpec[istage][idet][ibin_evis]->Reset();
      for(int ibin_enu=0;ibin_enu< n_enu_bins; ++ibin_enu){
        for(int icore=0;icore<Ncores;++icore){
          double Etrue=CorrEvtsTrueSpec[istage][idet][ibin_evis]->GetXaxis()->GetBinCenter(ibin_enu+1);
          int thebin=mapflux_noosc[icore]->GetXaxis()->FindBin(Etrue);
          double oscprob=mapextrap[icore]->GetBinContent(thebin);
          CorrEvtsTrueSpec[istage][idet][ibin_evis]
            ->AddBinContent(ibin_enu+1,oldval*M_evis_to_enu[ibin_evis][ibin_enu]*factor_noosc[ibin_enu][icore]*oscprob);

          norm += M_evis_to_enu[ibin_evis][ibin_enu]*factor_noosc[ibin_enu][icore]*oscprob;
        }
      }
      if (norm != 0)
        CorrEvtsTrueSpec[istage][idet][ibin_evis]->Scale(1/norm);

    }
    //Get events at each detector from each core (in terms of true E)
    std::map<int,TH1F*> mapflux = fluxcalc->CalculateFluxHistRow(idet,sin2theta13,istage,dm2,s22t14,dm2_41);
    for(int ibin_evis=0;ibin_evis < n_evis_bins; ++ibin_evis){
      for(int icore=0;icore<Ncores;++icore){
        if (FirstMakePrediction){
          //Need to ensure other istages are created
          for(int ii=0;ii<Nstage;++ii){
            sprintf(dummyname,"CorrEvtsCoreSpec_%i_%i_%i_%i",ii,idet,ibin_evis,icore);
            CorrEvtsCoreSpec[ii][idet][ibin_evis][icore]  = (TH1F*)CorrEvtsTrueSpec[ii][idet][ibin_evis]->Clone(dummyname);
          }
        }
        CorrEvtsCoreSpec[istage][idet][ibin_evis][icore]->Reset();
        //Note: multiplication has to be done manually in order to allow different bins between E hists and flux hists
        for(int ibin=1;ibin<=CorrEvtsTrueSpec[istage][idet][ibin_evis]->GetXaxis()->GetNbins();++ibin){
          double Etrue=CorrEvtsTrueSpec[istage][idet][ibin_evis]->GetXaxis()->GetBinCenter(ibin);
          //          if (Etrue > 8.7) Etrue = 8.7; // this is because the current flux histogram ends at 8.8 MeV
          int thebin=mapflux[icore]->GetXaxis()->FindBin(Etrue);
          double factor=mapflux[icore]->GetBinContent(thebin);
          CorrEvtsCoreSpec[istage][idet][ibin_evis][icore]->SetBinContent(ibin,CorrEvtsTrueSpec[istage][idet][ibin_evis]->GetBinContent(ibin)*factor);
        }
      }
    }
  }

  //Extrapolate
  for(int idet=4;idet<8;++idet){
    for(int idet2=0;idet2<4;++idet2){
      std::map<int,TH1F*> mapextrap = fluxcalc->ExtrapolationFactorRow(idet,idet2,sin2theta13,istage,dm2,s22t14,dm2_41);
      for(int ibin_evis=0;ibin_evis < n_evis_bins; ++ibin_evis){
        if (FirstMakePrediction){
          //Need to ensure other istages are created
          for(int ii=0;ii<Nstage;++ii){
            sprintf(dummyname,"PredEvtsTrueSpec_%i_%i_%i_%i",ii,idet,idet2,ibin_evis);
            PredEvtsTrueSpec[ii][idet][idet2][ibin_evis] = (TH1F*)CorrEvtsCoreSpec[ii][idet2][ibin_evis][0]->Clone(dummyname);
          }
        }
        PredEvtsTrueSpec[istage][idet][idet2][ibin_evis]->Reset();

        for(int icore=0;icore<Ncores;++icore){
          //same as above: have to do multiplication manually to allow for different bins
          for(int ibin=1;ibin<=CorrEvtsCoreSpec[istage][idet2][ibin_evis][icore]->GetXaxis()->GetNbins();++ibin){
            double Etrue=CorrEvtsCoreSpec[istage][idet2][ibin_evis][icore]->GetXaxis()->GetBinCenter(ibin);
            //            if (Etrue > 8.7) Etrue = 8.7; // this is because the current flux histogram ends at 8.8 MeV
            int thebin=mapextrap[icore]->GetXaxis()->FindBin(Etrue);
            double factor=mapextrap[icore]->GetBinContent(thebin);
            PredEvtsTrueSpec[istage][idet][idet2][ibin_evis]->AddBinContent(ibin,CorrEvtsCoreSpec[istage][idet2][ibin_evis][icore]->GetBinContent(ibin)*factor);
          }
        }
      }
    }
  }

  //Go back to Evis
  //PredSet *predout = new PredSet();
  for(int idet=4;idet<8;++idet){
    for(int idet2=0;idet2<4;++idet2){

      PredEvtsSpec[istage][idet][idet2]->Reset();

      //No prediction should be made if det or det2 were inactive
      if(NdetectorsConfig[istage][idet] != 0)
        if(NdetectorsConfig[istage][idet2] != 0)
          for(int ibin_evis=0;ibin_evis < n_evis_bins; ibin_evis++)
            PredEvtsSpec[istage][idet][idet2]->SetBinContent(ibin_evis+1,PredEvtsTrueSpec[istage][idet][idet2][ibin_evis]->Integral());

      //fill output object
      predper->SetPred(istage,idet,idet2,PredEvtsSpec[istage][idet][idet2]);
    }
  }

  FirstMakePrediction = false;

  //return predout;

}

void Predictor::CombineData(){

  for(int idet=0;idet<Ndetectors;++idet){
    double livsum=0;
    double weightedsum=0;
    double weightederr=0;
    double weightedbg=0;

    if (FirstMakeSuperPrediction){
      sprintf(dummyname,"CombCorrEvtsSpec_AD%i",idet);
      CombCorrEvtsSpec[idet] = (TH1F*)tdper[0].CorrEvtsSpec[idet]->Clone(dummyname);
      sprintf(dummyname,"CombCorrBgEvtsSpec_AD%i",idet);
      CombCorrBgEvtsSpec[idet] = (TH1F*)tdper[0].CorrBgEvtsSpec[idet]->Clone(dummyname);
      //note: do not set FirstMakeSuperPrediction to false as that is done below
    }
    CombCorrEvtsSpec[idet]->Reset();
    CombCorrBgEvtsSpec[idet]->Reset();

    for(int istage=0;istage<Nstage;++istage){

      float factor=tdper[istage].MuonVetoEff[idet]
        *tdper[istage].DMCEff[idet]
        *tdper[istage].Livetime[idet]
        *tdper[istage].TargetMass[idet]/tdper[istage].TargetMass[0];

      weightedsum+=tdper[istage].CorrEvts[idet]*factor;
      livsum+=factor;
      weightederr+=pow(tdper[istage].ErrEvts[idet]*factor,2);
      weightedbg+=tdper[istage].CorrBgEvts[idet]*factor;

      CombCorrEvtsSpec[idet]->Add(tdper[istage].CorrEvtsSpec[idet],factor);
      CombCorrBgEvtsSpec[idet]->Add(tdper[istage].CorrBgEvtsSpec[idet],factor);
    }

    if (livsum == 0){
      CombLivetime[idet]=0;
      CombCorrEvts[idet]=0;
      CombErrEvts[idet]=0;
      CombCorrBgEvts[idet]=0;
      CombCorrEvtsSpec[idet]->Scale(0);
      CombCorrBgEvtsSpec[idet]->Scale(0);
    }
    else{
      CombLivetime[idet]=livsum;
      CombCorrEvts[idet]=weightedsum*1./livsum;
      CombErrEvts[idet]=sqrt(weightederr)*1./livsum;
      CombCorrBgEvts[idet]=weightedbg*1./livsum;
      CombCorrEvtsSpec[idet]->Scale(1./livsum);
      CombCorrBgEvtsSpec[idet]->Scale(1./livsum);
    }

  }//for data

}//end of CombineData

void Predictor::CombinePredictions(){

  for(int idet=4;idet<8;++idet){
    for(int idet2=0;idet2<4;++idet2){

      double livsum=0;
      sprintf(dummyname,"PredEvtsSpec_%i_%i",idet,idet2);
      if (FirstMakeSuperPrediction){
        hweightedsum = (TH1F*)predper->GetPred(0,idet,idet2)->Clone(dummyname);
        FirstMakeSuperPrediction = false;
      }
      hweightedsum->Reset();

      for(int istage=0;istage<Nstage;++istage){

        //Skip over inactive AD
        if(NdetectorsConfig[istage][idet] == 0) continue;
        if(NdetectorsConfig[istage][idet2] == 0) continue;

        float factor=tdper[istage].MuonVetoEff[idet]
          *tdper[istage].DMCEff[idet]
          *tdper[istage].Livetime[idet]
          *tdper[istage].TargetMass[idet]/tdper[istage].TargetMass[0];
        hweightedsum->Add(predper->GetPred(istage,idet,idet2),factor);
        livsum+=factor;
      }

      if (livsum == 0)
        hweightedsum->Scale(0);
      else
        hweightedsum->Scale(1./livsum);

      superpred->SetPred(0,idet,idet2,hweightedsum);
    }
  }

}//CombinePredictions


PredSet *Predictor::MakeOneSuperPrediction(double sin2theta13, double dm2, double s22t14, double dm2_41, bool print){

  this->MakeAllPredictions(sin2theta13,dm2, s22t14, dm2_41,print);

  //PredSet *superpred = new PredSet();
  //print = true;
  if(print){
    CombineData();
    CombinePredictions();
    cout << "=================================================================" << endl;
    cout << "================= Super-prediction ==============================" << endl;
    superpred->PrintToScreen(0);
    cout << "================= Combined Data =================================" << endl;
    for(int idet=0;idet<Ndetectors;++idet){
      cout << "AD" << idet+1 << ": " << CombCorrEvtsSpec[idet]->Integral() << endl;//tmp
    }
  }

  return predper;

}//end of MakeOneSuperPrediction


void Predictor::MakeAllPredictions(double sin2theta13, double dm2, double s22t14, double dm2_41, bool print){
  //print = true;
  if(print){
    cout << "-----------------------------------------------------------" << endl;
    cout << "------------------- Predictions ---------------------------" << endl;
  }
  for(int istage=0;istage<Nstage;++istage){
    MakePrediction(sin2theta13,dm2,s22t14, dm2_41, istage,tdper[istage]); //,predper[istage]);

    if(print){
      cout << "Stage #" << istage << endl;
      predper->PrintToScreen(istage);
    }
  }
}


void Predictor::LoadEvisToEnuMatrix(const Char_t *evis_to_enu_matrix_name){
  ifstream matrix_file(evis_to_enu_matrix_name);
  if (!matrix_file.is_open()){
    cout << "Error: cannot find evis to enu convertion matrix " << evis_to_enu_matrix_name << endl;
    return;
  }
  for (Int_t ibinx = 0; ibinx < n_evis_bins; ibinx++ ){
    for (Int_t ibiny = 0; ibiny < n_enu_bins; ibiny++ ){
      matrix_file >>  M_evis_to_enu[ibinx][ibiny] ;
      //      cout <<  M_evis_to_enu[ibinx][ibiny] << endl;
    }
    //    cout <<  endl;
  }
}


void Predictor::LoadCovMatrix(const Char_t *covmatrixname_sig, const Char_t *covmatrixname_bg){

  string dummyLine;
  string thead;
  //-->Distances
  cout << " Signal coavariance matrix ++++++++++++++++++++++++++++++++++++++" << endl;
  ifstream covfile_sig(covmatrixname_sig);
  for (Int_t i = 0; i < MaxPredictions*Nstage*n_evis_bins; i++){
    for (Int_t j = 0; j < MaxPredictions*Nstage*n_evis_bins; j++){
      covfile_sig >> M_sig_sys[i][j];
      if (i%n_evis_bins == 0 && j%n_evis_bins == 0)
        cout <<  "\t" << M_sig_sys[i][j] ;
    }
    if (i%(Nstage*n_evis_bins) == 0)
      cout << endl;
  }

  cout << " Background coavariance matrix ++++++++++++++++++++++++++++++++++++++" << endl;
  ifstream covfile_bg(covmatrixname_bg);
  for (Int_t i = 0; i < MaxPredictions*Nstage*n_evis_bins; i++){
    for (Int_t j = 0; j < MaxPredictions*Nstage*n_evis_bins; j++){
      covfile_bg >> M_bg_sys[i][j];
      if (i%n_evis_bins == 0 && j%n_evis_bins == 0)
        cout << "\t" <<  M_bg_sys[i][j];
    }
    if (i%(Nstage*n_evis_bins) == 0)
      cout << endl;
  }
}

void Predictor::AddandScaleCovMatrix(Int_t type){

  for (Int_t i = 0; i < MaxPredictions*Nstage*n_evis_bins; i++){
    for (Int_t j = 0; j < MaxPredictions*Nstage*n_evis_bins; j++){
      M[i][j] = 0;
    }
  }

  if (type == -1)
    CalculateStatError();

  else if (type == 0) // remove far site stat error
    CalculateNearSiteStatError();

  for(int istage=0;istage<Nstage;istage++){
    for(int jstage=0;jstage<Nstage;jstage++){
      for(int idet=4;idet<8;++idet){
        for(int idet2=0;idet2<4;++idet2){
          for(int jdet=4;jdet<8;++jdet){
            for(int jdet2=0;jdet2<4;++jdet2){

              //Skip over inactive AD
              if(NdetectorsConfig[istage][idet] == 0) continue;
              if(NdetectorsConfig[istage][idet2] == 0) continue;

              if(NdetectorsConfig[jstage][jdet] == 0) continue;
              if(NdetectorsConfig[jstage][jdet2] == 0) continue;

              Int_t iii = (idet-4)*4 +idet2;
              Int_t jjj = (jdet-4)*4 +jdet2;
              for (Int_t i = 0; i < n_evis_bins; i++){
                for (Int_t j = 0; j < n_evis_bins; j++){

                  Int_t i2 = (istage*MaxPredictions+iii)*n_evis_bins+i;
                  Int_t j2 = (jstage*MaxPredictions+jjj)*n_evis_bins+j;

                  M_bg_sys_frac[i2][j2] = M_bg_sys[i2][j2]
                    /predper->GetPred(istage,idet,idet2)->GetBinContent(i+1)
                    /predper->GetPred(jstage,jdet,jdet2)->GetBinContent(j+1);

                  if (type == -1 || type == 0){ // all errors
                    M[i2][j2]
                      += M_sig_sys[i2][j2]
                      * predper->GetPred(istage,idet,idet2)->GetBinContent(i+1)
                      * predper->GetPred(jstage,jdet,jdet2)->GetBinContent(j+1)
                      + M_bg_sys[i2][j2] // already scaled
                      + M_stat[i2][j2];

                    /*
                      if(n_evis_bins==37){ //ensure to do this only for LBNL binning
                      if(i==j){ //only for diagonal elements
                      if(i<3){ //only for first three bins
                      //cout<<"Adding extra uncertainty"<<endl;
                      M[i2][j2]+=0.4*0.4* predper->GetPred(istage,idet,idet2)->GetBinContent(i+1)
                      * predper->GetPred(jstage,jdet,jdet2)->GetBinContent(j+1);  //adding 10%, 20%, 40% additional uncertainty
                      }
                      }
                      }*/


                  }else if (type ==1) {// signal systematics only
                    M[i2][j2]
                      += M_sig_sys[i2][j2] * predper->GetPred(istage,idet,idet2)->GetBinContent(i+1)
                      * predper->GetPred(jstage,jdet,jdet2)->GetBinContent(j+1);
                  }else if  (type ==2) {// backgrounds systematics only
                    M[i2][j2] += M_bg_sys[i2][j2]; // already scaed
                  }else if  (type ==3) {// near-site stat error only
                    M[i2][j2] +=  M_stat[i2][j2];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}


void Predictor::InvertMatrix(){

  Int_t npredictions = 2;

  if (combine_mode == 0)
    npredictions = 16;
  if (combine_mode == 1)
    npredictions = 2;
  if (combine_mode == 2)
    npredictions = 4;
  if (combine_mode == 3)
    npredictions = 1;


  Double_t * M_scaled_tmp = CombineMatrix(combine_mode);

  if (stage == -1){ //Use full covariance matrix
    TMatrixD * mat = new TMatrixD(npredictions*Nstage*n_evis_bins_rebin,npredictions*Nstage*n_evis_bins_rebin);
    mat->SetMatrixArray(M_scaled_tmp);
    //cout << "Matrix Determinant: " << mat->Determinant() << endl;
    //cout << "Matrix" << endl;
    //mat->Print();
    mat->Invert();
    //cout << "Inverted Matrix" << endl;
    //mat->Print();
    Double_t * m_tmp = mat->GetMatrixArray();
    for (Int_t i = 0; i < npredictions*Nstage*n_evis_bins_rebin; i++){
      for (Int_t j = 0; j < npredictions*Nstage*n_evis_bins_rebin; j++){
        M_inv[i][j] = m_tmp[i*npredictions*Nstage*n_evis_bins_rebin+j];
        //cout << " " << M_inv[i][j];
      }
      //cout << endl;
    }
    delete mat;
  }

  else if(stage == 7){ //added by Beda - 6+8 AD period
    int Stages_in_68AD=2;

    TMatrixD * mat = new TMatrixD(npredictions*Stages_in_68AD*n_evis_bins_rebin,npredictions*Stages_in_68AD*n_evis_bins_rebin);

    Double_t M_scaled_tmp_stage[npredictions*Stages_in_68AD*n_evis_bins_rebin*npredictions*Stages_in_68AD*n_evis_bins_rebin];


    for(int i=0;i<npredictions*Stages_in_68AD*n_evis_bins_rebin;++i){
      for(int j=0;j<npredictions*Stages_in_68AD*n_evis_bins_rebin;++j){
        M_scaled_tmp_stage[i*Stages_in_68AD*npredictions*n_evis_bins_rebin+j] = M_scaled_tmp[i*npredictions*Nstage*n_evis_bins_rebin+j];
        //check it this is correct stage->Stages_in_68AD
      }
    }


    mat->SetMatrixArray(M_scaled_tmp_stage);
    //cout << "Matrix Determinant: " << mat->Determinant() << endl;
    //cout << "Matrix" << endl;
    //mat->Print();

    //cout<<"Inverting normal matrix"<<endl;

    mat->Invert();
    //cout << "Inverted Matrix" << endl;
    //mat->Print();
    Double_t * m_tmp = mat->GetMatrixArray();

    for(int istage=0;istage<Nstage;istage++){
      for(int jstage=0;jstage<Nstage;jstage++){
        for (Int_t i = 0; i < npredictions*n_evis_bins_rebin; i++){
          for (Int_t j = 0; j < npredictions*n_evis_bins_rebin; j++){

            if(istage!=2 && jstage!=2)  //skip 7AD stage
              M_inv[istage*npredictions*n_evis_bins_rebin+i][jstage*npredictions*n_evis_bins_rebin+j] = m_tmp[(i+istage*npredictions*n_evis_bins_rebin)*npredictions*Stages_in_68AD*n_evis_bins_rebin+(jstage*npredictions*n_evis_bins_rebin+j)];
            else
              M_inv[istage*npredictions*n_evis_bins_rebin+i][jstage*npredictions*n_evis_bins_rebin+j] =0;
            //cout << " " << M_inv[i][j];
          }
          //cout << endl;
        }
      }
    }
    delete mat;
  }
  else{ //Need to be very careful on how to get the correct portion of the covmatrix
    TMatrixD * mat = new TMatrixD(npredictions*n_evis_bins_rebin,npredictions*n_evis_bins_rebin);

    Double_t M_scaled_tmp_stage[npredictions*n_evis_bins_rebin*npredictions*n_evis_bins_rebin];

    for(int i=0;i<npredictions*n_evis_bins_rebin;++i){
      for(int j=0;j<npredictions*n_evis_bins_rebin;++j){
        M_scaled_tmp_stage[i*npredictions*n_evis_bins_rebin+j] = M_scaled_tmp[(stage*npredictions*n_evis_bins_rebin+i)*npredictions*Nstage*n_evis_bins_rebin+(stage*npredictions*n_evis_bins_rebin+j)];
        //mat(i,j) = M_scaled_tmp[(stage*npredictions*n_evis_bins_rebin+i)*npredictions*Nstage*n_evis_bins_rebin+(stage*npredictions*n_evis_bins_rebin+j)];
      }
    }

    mat->SetMatrixArray(M_scaled_tmp_stage);
    //cout << "Matrix Determinant: " << mat->Determinant() << endl;
    //cout << "Matrix" << endl;
    //mat->Print();

    mat->Invert();
    //cout << "Inverted Matrix" << endl;
    //mat->Print();
    Double_t * m_tmp = mat->GetMatrixArray();

    for(int istage=0;istage<Nstage;istage++){
      for(int jstage=0;jstage<Nstage;jstage++){
        for (Int_t i = 0; i < npredictions*n_evis_bins_rebin; i++){
          for (Int_t j = 0; j < npredictions*n_evis_bins_rebin; j++){

            if(istage==stage && jstage==stage)
              M_inv[istage*npredictions*n_evis_bins_rebin+i][jstage*npredictions*n_evis_bins_rebin+j] = m_tmp[i*npredictions*n_evis_bins_rebin+j];
            else
              M_inv[istage*npredictions*n_evis_bins_rebin+i][jstage*npredictions*n_evis_bins_rebin+j] =0;
            //cout << " " << M_inv[i][j];
          }
          //cout << endl;
        }
      }
    }
    delete mat;
  }

}

void Predictor::InvertRateMatrix(){

  cout<<"Calculating rate matrix for some reason"<<endl;
  Int_t StartStage, EndStage, NumStage;

  if (stage == -1){ //Use full covariance matrix
    cout<<"Using full cov matrix"<<endl;
    StartStage = 0;
    EndStage = Nstage;
    NumStage = Nstage;
  }
  else if(stage==7){
    StartStage=0;
    EndStage=2;
    NumStage = 2;
  }
  else{
    StartStage = stage;
    EndStage = stage+1;
    NumStage = 1;
  }


  Int_t npredictions = 2;

  if (combine_mode == 0)
    npredictions = 16;
  if (combine_mode == 1)
    npredictions = 2;
  if (combine_mode == 2)
    npredictions = 4;
  if (combine_mode == 3)
    npredictions = 1;

  if(combine_mode==1){

    TMatrixD * mat = new TMatrixD(nNearHalls*NumStage,nNearHalls*NumStage);

    Double_t M_scaled_tmp[nNearHalls*NumStage*nNearHalls*NumStage];
    for (Int_t i = 0; i < nNearHalls*NumStage; i++){
      for (Int_t j = 0; j < nNearHalls*NumStage; j++){
        M_scaled[i][j] = 0;
        M_scaled_tmp[i*nNearHalls*NumStage+j] = 0;
      }
    }


    for(int istage=StartStage;istage<EndStage;istage++){
      for(int jstage=StartStage;jstage<EndStage;jstage++){
        for(int idet=4;idet<8;++idet){
          if(NdetectorsConfig[istage][idet] == 0) continue;
          for(int idet2=0;idet2<4;++idet2){
            if(NdetectorsConfig[istage][idet2] == 0) continue;
            for(int jdet=4;jdet<8;++jdet){
              if(NdetectorsConfig[jstage][jdet] == 0) continue;
              for(int jdet2=0;jdet2<4;++jdet2){
                if(NdetectorsConfig[jstage][jdet2] == 0) continue;

                Int_t iii = (idet-4)*4 +idet2;
                Int_t jjj = (jdet-4)*4 +jdet2;

                for(int ie=0;ie<n_evis_bins;++ie){
                  for(int je=0;je<n_evis_bins;++je){

                    if (stage==-1 || stage==7)  //I guess this is proper extension
                      M_scaled[nNearHalls*istage+(detConfigEH[idet2]-1)][nNearHalls*jstage+(detConfigEH[jdet2]-1)]
                        += M[(istage*MaxPredictions+iii)*n_evis_bins+ie][(jstage*MaxPredictions+jjj)*n_evis_bins+je];
                    else
                      M_scaled[(detConfigEH[idet2]-1)][(detConfigEH[jdet2]-1)]
                        += M[(istage*MaxPredictions+iii)*n_evis_bins+ie][(jstage*MaxPredictions+jjj)*n_evis_bins+je];
                  }
                }
              }
            }
          }
        }
      }
    }


    for (Int_t i = 0; i < nNearHalls*NumStage; i++){
      for (Int_t j = 0; j < nNearHalls*NumStage; j++){
        M_scaled_tmp[i*nNearHalls*NumStage+j] = M_scaled[i][j];
      }
    }
    mat->SetMatrixArray(&M_scaled_tmp[0]);
    //mat->Print();
    mat->Invert();
    //mat->Print();
    Double_t * m_tmp = mat->GetMatrixArray();

    for (Int_t istage = 0; istage < Nstage; istage++){
      for (Int_t jstage = 0; jstage < Nstage; jstage++){
        for (Int_t i = 0; i < nNearHalls; i++){
          for (Int_t j = 0; j < nNearHalls; j++){

            if (stage==-1)
              M_rate_inv[istage*nNearHalls+i][jstage*nNearHalls+j] = m_tmp[((istage*nNearHalls)+i)*nNearHalls*Nstage+(jstage*nNearHalls)+j];
            else if(stage==7){
              if(istage<2 && jstage<2){
                M_rate_inv[istage*nNearHalls+i][jstage*nNearHalls+j] = m_tmp[((istage*nNearHalls)+i)*nNearHalls*Nstage+(jstage*nNearHalls)+j];
              }
              else
                M_rate_inv[istage*nNearHalls+i][jstage*nNearHalls+j] = 0;
            }
            else{
              if (istage==stage && jstage==stage)
                M_rate_inv[istage*nNearHalls+i][jstage*nNearHalls+j] = m_tmp[i*nNearHalls+j];
              else
                M_rate_inv[istage*nNearHalls+i][jstage*nNearHalls+j] = 0;
            }

            //      cout << " " << M_inv[i][j];
          }
          //    cout << endl;
        }
      }
    }

    /*cout<<"Inverted matrix:"<<endl;
      for (Int_t i = 0; i < Nstage*npredictions; i++){
      for (Int_t j = 0; j < Nstage*npredictions; j++){
      cout<<M_rate_inv[i][j]<<" ";
      }
      cout<<endl;
      }*/

    delete mat;
  }
  if(combine_mode==0){
    if(n_evis_bins_rebin!=n_evis_bins) cout<<"Bin mismatch for some reason"<<endl;
    TMatrixD * mat = new TMatrixD(npredictions*Nstage,npredictions*Nstage);

    Double_t M_scaled_tmp[npredictions*Nstage*npredictions*Nstage];
    for (Int_t i = 0; i < npredictions*Nstage; i++){
      for (Int_t j = 0; j < npredictions*Nstage; j++){
        M_scaled[i][j] = 0;
        M_scaled_tmp[i*npredictions*Nstage+j] = 0;
      }
    }

    for(int istage=0;istage<Nstage;istage++){
      for(int jstage=0;jstage<Nstage;jstage++){
        for(int idet=4;idet<8;++idet){
          if(NdetectorsConfig[istage][idet] == 0) continue;
          for(int idet2=0;idet2<4;++idet2){
            if(NdetectorsConfig[istage][idet2] == 0) continue;
            for(int jdet=4;jdet<8;++jdet){
              if(NdetectorsConfig[jstage][jdet] == 0) continue;
              for(int jdet2=0;jdet2<4;++jdet2){
                if(NdetectorsConfig[jstage][jdet2] == 0) continue;

                Int_t iii = (idet-4)*4 +idet2;
                Int_t jjj = (jdet-4)*4 +jdet2;

                for(int ie=0;ie<n_evis_bins;++ie){
                  for(int je=0;je<n_evis_bins;++je){


                    M_scaled[npredictions*istage+iii][npredictions*jstage+jjj]
                      += M[(istage*MaxPredictions+iii)*n_evis_bins+ie][(jstage*MaxPredictions+jjj)*n_evis_bins+je];
                  }
                }
              }
            }
          }
        }
      }
    }

    for (Int_t istage = 0; istage < Nstage; istage++){
      for (Int_t jstage = 0; jstage < Nstage; jstage++){
        for (Int_t i = 0; i < npredictions; i++){
          for (Int_t j = 0; j < npredictions; j++){
            M_scaled_tmp[(istage*npredictions+i)*npredictions*Nstage+jstage*npredictions+j] = M_scaled[istage*npredictions+i][jstage*npredictions+j];
          }
        }
      }
    }
    mat->SetMatrixArray(&M_scaled_tmp[0]);
    //mat->Print();
    mat->Invert();
    //mat->Print();
    Double_t * m_tmp = mat->GetMatrixArray();

    for (Int_t istage = 0; istage < Nstage; istage++){
      for (Int_t jstage = 0; jstage < Nstage; jstage++){
        for (Int_t i = 0; i < npredictions; i++){
          for (Int_t j = 0; j < npredictions; j++){


            M_rate_inv[istage*npredictions+i][jstage*npredictions+j] = m_tmp[((istage*npredictions)+i)*npredictions*Nstage+(jstage*npredictions)+j];


            //      cout << " " << M_inv[i][j];
          }
          //    cout << endl;
        }
      }
    }

    delete mat;


    cout<<"Inverted matrix:"<<endl;
    for (Int_t i = 0; i < Nstage*npredictions; i++){
      for (Int_t j = 0; j < Nstage*npredictions; j++){
        cout<<M_rate_inv[i][j]<<" ";
      }
      cout<<endl;
    }


  }
  cout<<"Inversion done"<<endl;
}

Double_t Predictor::CalculateChi2Cov(Double_t sin22t13, Double_t dm2_ee,Double_t sin22t14, Double_t dm2_41){
  this->MakeAllPredictions(sin22t13,dm2_ee, sin22t14, dm2_41,false);
  return this->CalculateChi2Cov();
}

Double_t Predictor::CalculateChi2Cov(){

  Int_t npredictions = 2;

  if (combine_mode == 0)
    npredictions = 16;
  if (combine_mode == 1)
    npredictions = 2;
  if (combine_mode == 2)
    npredictions = 4;
  if (combine_mode == 3)
    npredictions = 1;

  if (RecalculateCovMatrix){
    AddandScaleCovMatrix();
    InvertMatrix();
  }

  Double_t N_pred[npredictions*Nstage*n_evis_bins];
  Double_t N_obs[npredictions*Nstage*n_evis_bins];

  for (Int_t i = 0; i < npredictions*Nstage*n_evis_bins; i++){
    N_pred[i] = 0;
    N_obs[i] = 0;
  }

  // Fill observation and prediction
  for(int istage=0;istage<Nstage;++istage){
    for(int idet=4;idet<8;++idet){
      for(int idet2=0;idet2<4;++idet2){

        //Skip over inactive AD
        if(NdetectorsConfig[istage][idet] == 0) continue;
        if(NdetectorsConfig[istage][idet2] == 0) continue;

        for (Int_t i = 0; i < n_evis_bins; i++){

          if (combine_mode == 0){
            Int_t iii = (idet-4)*4 +idet2;
            N_obs[(npredictions*istage+iii)*n_evis_bins_rebin+evis_rebin_map[i]] += tdper[istage].CorrEvtsSpec[idet]->GetBinContent(i+1);
            N_pred[(npredictions*istage+iii)*n_evis_bins_rebin+evis_rebin_map[i]] += predper->GetPred(istage,idet,idet2)->GetBinContent(i+1);
          }

          if (combine_mode == 1){
            N_obs[(npredictions*istage+(detConfigEH[idet2]-1))*n_evis_bins_rebin+evis_rebin_map[i]] +=  mode1_coeff[istage][idet2] * tdper[istage].CorrEvtsSpec[idet]->GetBinContent(i+1);
            N_pred[(npredictions*istage+(detConfigEH[idet2]-1))*n_evis_bins_rebin+evis_rebin_map[i]] += mode1_coeff[istage][idet2] * predper->GetPred(istage,idet,idet2)->GetBinContent(i+1);
          }

          if (combine_mode == 2){
            N_obs[(npredictions*istage+idet2)*n_evis_bins_rebin+evis_rebin_map[i]] += NearEH_coeff[istage][idet2] * tdper[istage].CorrEvtsSpec[idet]->GetBinContent(i+1);
            N_pred[(npredictions*istage+idet2)*n_evis_bins_rebin+evis_rebin_map[i]] += NearEH_coeff[istage][idet2] * predper->GetPred(istage,idet,idet2)->GetBinContent(i+1);
          }

          if (combine_mode == 3){
            N_obs[npredictions*istage*n_evis_bins_rebin+evis_rebin_map[i]] += NearEH_coeff[istage][idet2] * tdper[istage].CorrEvtsSpec[idet]->GetBinContent(i+1);
            N_pred[npredictions*istage*n_evis_bins_rebin+evis_rebin_map[i]] += NearEH_coeff[istage][idet2] * predper->GetPred(istage,idet,idet2)->GetBinContent(i+1);
          }

        }
      }
    }
  }

  //Choose which stage to minimize chi2
  //stage==-1 corresponds to minimizing over all stages

  Int_t startStage;
  Int_t endStage;

  if (stage==-1){
    startStage = 0;
    endStage = Nstage;
  }
  else if(stage==7){
    startStage=0;
    endStage=2;
  }
  else{
    startStage = stage;
    endStage = stage+1;
  }

  Double_t chi2out=0;
  for (Int_t i = npredictions*startStage*n_evis_bins_rebin; i < npredictions*endStage*n_evis_bins_rebin; i++){
    //cout << "i: " << i << " N_obs: " << N_obs[i]<< " N_pred: " << N_pred[i] << endl;
    for (Int_t j = npredictions*startStage*n_evis_bins_rebin; j < npredictions*endStage*n_evis_bins_rebin; j++){
      //cout << "j: " << j << " N_obs: " << N_obs[j]<< " N_pred: " << N_pred[j] << endl;

      //cout << "M_inv: " << M_inv[i][j] << " chi2: " << M_inv[i][j]*(N_obs[i]-N_pred[i])*(N_obs[j]-N_pred[j]) << endl;
      //for (Int_t i = 0; i < npredictions*Nstage*n_evis_bins_rebin; i++){
      //for (Int_t j = 0; j < npredictions*Nstage*n_evis_bins_rebin; j++){
      chi2out+= M_inv[i][j]
        * (N_obs[i]-N_pred[i])
        * (N_obs[j]-N_pred[j]);
    }
  }
  //cout << "Chi2: " << chi2out << endl;
  return chi2out;

}

Double_t Predictor::CalculateChi2CovRate(Double_t sin22t13, Double_t dm2_ee,Double_t sin22t14, Double_t dm2_41){
  this->MakeAllPredictions(sin22t13,dm2_ee, sin22t14, dm2_41,false);
  return this->CalculateChi2CovRate();
}

Double_t Predictor::CalculateChi2CovRate(){

  Int_t npredictions = 2;

  if (combine_mode == 0)
    npredictions = 16;
  if (combine_mode == 1)
    npredictions = 2;
  if (combine_mode == 2)
    npredictions = 4;
  if (combine_mode == 3)
    npredictions = 1;

  if (RecalculateCovMatrix){
    AddandScaleCovMatrix();
    InvertRateMatrix();
  }

  Double_t N_pred[npredictions*Nstage];
  Double_t N_obs[npredictions*Nstage];

  for (Int_t i = 0; i < npredictions*Nstage; i++){
    N_pred[i] = 0;
    N_obs[i] = 0;
  }

  cout<<"Here"<<endl;

  // Fill observation and predictions
  for(int istage=0;istage<Nstage;istage++){
    for(int idet=4;idet<8;++idet){
      for(int idet2=0;idet2<4;++idet2){

        //Skip over inactive AD
        if(NdetectorsConfig[istage][idet] == 0) continue;
        if(NdetectorsConfig[istage][idet2] == 0) continue;

        if(combine_mode==0){
          Int_t iii = (idet-4)*4 +idet2;
          for (Int_t i = 0; i < n_evis_bins; i++){
            N_obs[istage*npredictions+iii] += tdper[istage].CorrEvtsSpec[idet]->GetBinContent(i+1);
            N_pred[istage*npredictions+iii] += predper->GetPred(istage,idet,idet2)->GetBinContent(i+1);
          }
        }
        if(combine_mode==1){
          for (Int_t i = 0; i < n_evis_bins; i++){
            N_obs[istage*npredictions+(detConfigEH[idet2]-1)] += tdper[istage].CorrEvtsSpec[idet]->GetBinContent(i+1);
            N_pred[istage*nNearHalls+(detConfigEH[idet2]-1)] += predper->GetPred(istage,idet,idet2)->GetBinContent(i+1);
          }
        }
      }
    }
  }

  cout<<"Here2"<<endl;

  Double_t chi2out=0;
  for (Int_t istage = 0; istage < 3; istage++){
    for (Int_t jstage = 0; jstage < 3; jstage++){
      for (Int_t i = 0; i < npredictions; i++){
        for (Int_t j = 0; j < npredictions; j++){

          chi2out+= M_rate_inv[istage*npredictions+i][jstage*npredictions+j]
            * (N_obs[istage*npredictions+i]-N_pred[istage*npredictions+i])
            * (N_obs[jstage*npredictions+j]-N_pred[jstage*npredictions+j]);

        }
      }
    }
  }
  //cout << chi2out << " " << N_obs[0]  << " " << N_obs[1]  << " " << N_pred[0] << " " << N_pred[1] << endl;

  return chi2out;

}

void Predictor::FixCovMatrix(Double_t sin22t13, Double_t dm2_ee,Double_t sin22t14, Double_t dm2_41){
  this->MakeAllPredictions(sin22t13,dm2_ee, sin22t14, dm2_41,false);

  if (RecalculateCovMatrix){
    AddandScaleCovMatrix();
    InvertMatrix();
  }

  cout << "Calculated and fixed the covariance matrix using parameters: " << sin22t13 << " " <<  dm2_ee<< " " << sin22t14 << " " << dm2_41 << " Stage: " << stage << endl;
  RecalculateCovMatrix = false;

  return;

}

//////////////////////////////////////
Double_t * Predictor::GetFinalPred(PredSet *evtset, Int_t mode){

  Int_t npredictions = 2;

  if (mode == 0)
    npredictions = 16;
  if (mode == 1)
    npredictions = 2;
  if (mode == 2)
    npredictions = 4;
  if (mode == 3)
    npredictions = 1;

  for (Int_t i = 0; i < MaxPredictions*Nstage*n_evis_bins; i++){
    final_pred[i] = 0;
  }

  for(int istage=0;istage<Nstage;++istage){
    for(int idet=4;idet<8;++idet){
      for(int idet2=0;idet2<4;++idet2){
        Int_t iii = (idet-4)*4 +idet2;
        for (Int_t i = 0; i < n_evis_bins; i++){
          if (mode == 0)
            final_pred[(npredictions*istage+iii)*n_evis_bins_rebin+evis_rebin_map[i]] += evtset->GetPred(istage,idet,idet2)->GetBinContent(i+1);
          if (mode == 1)
            final_pred[(npredictions*istage+(detConfigEH[idet2]-1))*n_evis_bins_rebin+evis_rebin_map[i]] += mode1_coeff[istage][idet2] * evtset->GetPred(istage,idet,idet2)->GetBinContent(i+1);
          if (mode == 2)
            final_pred[(npredictions*istage+idet2)*n_evis_bins_rebin+evis_rebin_map[i]] += evtset->GetPred(istage,idet,idet2)->GetBinContent(i+1);
          if (mode == 3)
            final_pred[npredictions*istage*n_evis_bins_rebin+evis_rebin_map[i]] += NearEH_coeff[istage][idet2] * evtset->GetPred(istage,idet,idet2)->GetBinContent(i+1);

        }
      }
    }
  }
  return &final_pred[0];

}
//////////////////////////////////////

Double_t * Predictor::GetFinalObs(Int_t mode){

  Int_t npredictions = 2;

  if (mode == 0)
    npredictions = 16;
  if (mode == 1)
    npredictions = 2;
  if (mode == 2)
    npredictions = 4;
  if (mode == 3)
    npredictions = 1;

  for (Int_t i = 0; i < MaxPredictions*Nstage*n_evis_bins; i++){
    final_obs[i] = 0;
  }

  for(int istage=0;istage<Nstage;istage++){
    for(int idet=4;idet<8;++idet){
      for(int idet2=0;idet2<4;++idet2){
        Int_t iii = (idet-4)*4 +idet2;
        for (Int_t i = 0; i < n_evis_bins; i++){
          if (mode == 0)
            final_obs[(istage*npredictions+iii)*n_evis_bins_rebin+evis_rebin_map[i]] += tdper[istage].CorrEvtsSpec[idet]->GetBinContent(i+1);
          if (mode == 1)
            final_obs[(istage*npredictions+detConfigEH[idet2]-1)*n_evis_bins_rebin+evis_rebin_map[i]] += mode1_coeff[istage][idet2] * tdper[istage].CorrEvtsSpec[idet]->GetBinContent(i+1);
          if (mode == 2)
            final_obs[(istage*npredictions+idet2)*n_evis_bins_rebin+evis_rebin_map[i]] += tdper[istage].CorrEvtsSpec[idet]->GetBinContent(i+1);
          if (mode == 3)
            final_obs[istage*npredictions*n_evis_bins_rebin+evis_rebin_map[i]] += NearEH_coeff[istage][idet2] * tdper[istage].CorrEvtsSpec[idet]->GetBinContent(i+1);
        }
      }
    }
  }
  return &final_obs[0];

}

//////////////////////////////////////

Double_t * Predictor::GetFinalObsError(Int_t mode){
  Int_t npredictions = 2;

  if (mode == 0)
    npredictions = 16;
  if (mode == 1)
    npredictions = 2;
  if (mode == 2)
    npredictions = 4;
  if (mode == 3)
    npredictions = 1;

  for (Int_t i = 0; i < MaxPredictions*Nstage*n_evis_bins; i++){
    final_obserror[i] = 0;
  }

  for(int istage=0;istage<Nstage;istage++){
    for(int idet=4;idet<8;++idet){
      for(int idet2=0;idet2<4;++idet2){
        Int_t iii = (idet-4)*4 +idet2;
        for (Int_t i = 0; i < n_evis_bins; i++){
          if (mode == 0)
            final_obserror[(istage*npredictions+iii)*n_evis_bins_rebin+evis_rebin_map[i]] += pow(tdper[istage].CorrEvtsSpec[idet]->GetBinError(i+1),2);
          if (mode == 1)
            final_obserror[(istage*npredictions+detConfigEH[idet2]-1)*n_evis_bins_rebin+evis_rebin_map[i]] += mode1_coeff[istage][idet2] * pow(tdper[istage].CorrEvtsSpec[idet]->GetBinError(i+1),2);
          if (mode == 2)
            final_obserror[(istage*npredictions+idet2)*n_evis_bins_rebin+evis_rebin_map[i]] += pow(tdper[istage].CorrEvtsSpec[idet]->GetBinError(i+1),2);
          if (mode == 3)
            final_obserror[istage*npredictions*n_evis_bins_rebin+evis_rebin_map[i]] += NearEH_coeff[istage][idet2] * pow(tdper[istage].CorrEvtsSpec[idet]->GetBinError(i+1),2);
        }
      }
    }
  }

  for (Int_t i = 0; i < Nstage*npredictions*n_evis_bins_rebin; i++){
    final_obserror[i] = sqrt(final_obserror[i]);
  }
  return &final_obserror[0];

}

//////////////////////////////////////
//////////////////////////////////////
Double_t * Predictor::GetFinalPredSum(PredSet *evtset, Int_t mode){

  Int_t npredictions = 2;

  if (mode == 0)
    npredictions = 16;
  if (mode == 1)
    npredictions = 2;
  if (mode == 2)
    npredictions = 4;
  if (mode == 3)
    npredictions = 1;

  for (Int_t i = 0; i < MaxPredictions*n_evis_bins; i++){
    final_pred_sum[i] = 0;
  }

  for(int istage=0;istage<Nstage;++istage){
    for(int idet=4;idet<8;++idet){
      for(int idet2=0;idet2<4;++idet2){
        Int_t iii = (idet-4)*4 +idet2;
        for (Int_t i = 0; i < n_evis_bins; i++){
          if (mode == 0){
            final_pred_sum[iii*n_evis_bins_rebin+evis_rebin_map[i]] += FarAD_stage_coeff[istage][idet] * evtset->GetPred(istage,idet,idet2)->GetBinContent(i+1);
          }
          if (mode == 1)
            final_pred_sum[(detConfigEH[idet2]-1)*n_evis_bins_rebin+evis_rebin_map[i]] += FarEH_stage_coeff[istage][idet] * mode1_coeff[istage][idet2] * evtset->GetPred(istage,idet,idet2)->GetBinContent(i+1);
          if (mode == 2)
            final_pred_sum[idet2*n_evis_bins_rebin+evis_rebin_map[i]] += FarAD_stage_coeff[istage][idet] * evtset->GetPred(istage,idet,idet2)->GetBinContent(i+1);
          if (mode == 3)
            final_pred_sum[evis_rebin_map[i]] += FarEH_stage_coeff[istage][idet]* NearEH_coeff[istage][idet2] * evtset->GetPred(istage,idet,idet2)->GetBinContent(i+1);

        }
      }
    }
  }
  return &final_pred_sum[0];

}
//////////////////////////////////////

Double_t * Predictor::GetFinalObsSum(Int_t mode){

  Int_t npredictions = 2;

  if (mode == 0)
    npredictions = 16;
  if (mode == 1)
    npredictions = 2;
  if (mode == 2)
    npredictions = 4;
  if (mode == 3)
    npredictions = 1;

  for (Int_t i = 0; i < MaxPredictions*n_evis_bins; i++){
    final_obs_sum[i] = 0;
  }

  for(int istage=0;istage<Nstage;istage++){
    for(int idet=4;idet<8;++idet){
      for(int idet2=0;idet2<4;++idet2){
        Int_t iii = (idet-4)*4 +idet2;
        for (Int_t i = 0; i < n_evis_bins; i++){
          if (mode == 0)
            final_obs_sum[iii*n_evis_bins_rebin+evis_rebin_map[i]] += FarAD_stage_coeff[istage][idet] * tdper[istage].CorrEvtsSpec[idet]->GetBinContent(i+1);
          if (mode == 1)
            final_obs_sum[(detConfigEH[idet2]-1)*n_evis_bins_rebin+evis_rebin_map[i]] += FarEH_stage_coeff[istage][idet] * mode1_coeff[istage][idet2] * tdper[istage].CorrEvtsSpec[idet]->GetBinContent(i+1);
          if (mode == 2)
            final_obs_sum[idet2*n_evis_bins_rebin+evis_rebin_map[i]] += FarAD_stage_coeff[istage][idet] * tdper[istage].CorrEvtsSpec[idet]->GetBinContent(i+1);
          if (mode == 3)
            final_obs_sum[evis_rebin_map[i]] += FarEH_stage_coeff[istage][idet]* NearEH_coeff[istage][idet2] * tdper[istage].CorrEvtsSpec[idet]->GetBinContent(i+1);
        }
      }
    }
  }
  return &final_obs_sum[0];

}

//////////////////////////////////////

Double_t * Predictor::GetFinalObsErrorSum(Int_t mode){
  Int_t npredictions = 2;

  if (mode == 0)
    npredictions = 16;
  if (mode == 1)
    npredictions = 2;
  if (mode == 2)
    npredictions = 4;
  if (mode == 3)
    npredictions = 1;

  for (Int_t i = 0; i < MaxPredictions*n_evis_bins; i++){
    final_obserror_sum[i] = 0;
  }

  for(int istage=0;istage<Nstage;istage++){
    for(int idet=4;idet<8;++idet){
      for(int idet2=0;idet2<4;++idet2){
        Int_t iii = (idet-4)*4 +idet2;
        for (Int_t i = 0; i < n_evis_bins; i++){
          if (mode == 0)
            final_obserror_sum[iii*n_evis_bins_rebin+evis_rebin_map[i]] += pow(FarAD_stage_coeff[istage][idet] * tdper[istage].CorrEvtsSpec[idet]->GetBinError(i+1),2);
          if (mode == 1)
            final_obserror_sum[(detConfigEH[idet2]-1)*n_evis_bins_rebin+evis_rebin_map[i]] += mode1_coeff[istage][idet2] * pow(FarEH_stage_coeff[istage][idet] * tdper[istage].CorrEvtsSpec[idet]->GetBinError(i+1),2);
          if (mode == 2)
            final_obserror_sum[idet2*n_evis_bins_rebin+evis_rebin_map[i]] += pow(FarAD_stage_coeff[istage][idet] * tdper[istage].CorrEvtsSpec[idet]->GetBinError(i+1),2);
          if (mode == 3)
            final_obserror_sum[evis_rebin_map[i]] += NearEH_coeff[istage][idet2] * pow(FarEH_stage_coeff[istage][idet]* tdper[istage].CorrEvtsSpec[idet]->GetBinError(i+1),2);
        }
      }
    }
  }

  for (Int_t i = 0; i < npredictions*n_evis_bins_rebin; i++){
    final_obserror_sum[i] = sqrt(final_obserror_sum[i]);
  }
  return &final_obserror_sum[0];

}

//////////////////////////////////////

void Predictor::CalculateStatError(){
  for (int i = 0; i < MaxPredictions*Nstage*n_evis_bins; i++){
    for (int j = 0; j < MaxPredictions*Nstage*n_evis_bins; j++){
      M_stat[i][j] = 0;
      M_stat_frac[i][j] = 0;
    }
  }

  for(int istage=0;istage<Nstage;++istage){
    //for(int jstage=Nstage;jstage<Nstage;++jstage){
    int jstage = istage; //Only istage==jstage is important
    for(int idet=4;idet<8;++idet){
      for(int idet2=0;idet2<4;++idet2){
        for(int jdet=4;jdet<8;++jdet){
          for(int jdet2=0;jdet2<4;++jdet2){
            Int_t iii = (idet-4)*4 +idet2;
            Int_t jjj = (jdet-4)*4 +jdet2;

            if (idet2 == jdet2){ // sharing the same near detector during the same stage
              //cout << "sharing near detector " << iii << " " << jjj << endl;
              for (Int_t i = 0; i < n_evis_bins; i++){
                Int_t j = i; //Only i==j matters
                //Scale near site statistical error ( = sqrt(N_tot^2 + N_bg^2)) to the far site
                //"evtset-> GetEvts(idet,idet2) / CorrEvts[idet2]" gives the far/near ratio at given theta 13 value

                Double_t err_i = 0;
                if (tdper[istage].CorrEvtsSpec[idet2]->GetBinContent(i+1) > 0)
                  err_i = tdper[istage].CorrEvtsSpec[idet2]->GetBinError(i+1) * predper->GetPred(istage,idet,idet2)->GetBinContent(i+1) / tdper[istage].CorrEvtsSpec[idet2]->GetBinContent(i+1);
                Double_t err_j = 0;
                if (tdper[jstage].CorrEvtsSpec[jdet2]->GetBinContent(j+1) > 0)
                  err_j = tdper[jstage].CorrEvtsSpec[jdet2]->GetBinError(j+1) * predper->GetPred(jstage,jdet,jdet2)->GetBinContent(j+1) / tdper[jstage].CorrEvtsSpec[jdet2]->GetBinContent(j+1);


                M_stat[(MaxPredictions*istage+iii)*n_evis_bins+i][(MaxPredictions*jstage+jjj)*n_evis_bins+j] += err_i * err_j / stat_factor[istage];

                M_stat_frac[(MaxPredictions*istage+iii)*n_evis_bins+i][(MaxPredictions*jstage+jjj)*n_evis_bins+j] += err_i * err_j /predper->GetPred(istage,idet,idet2)->GetBinContent(i+1)/predper->GetPred(jstage,jdet,jdet2)->GetBinContent(j+1);

                //cout << "M_stat: " <<  M_stat[(MaxPredictions*istage+iii)*n_evis_bins+i][(MaxPredictions*jstage+jjj)*n_evis_bins+j] << endl;
              }

            }


            if (idet == jdet){ // sharing the same far detector during the same stage

              //cout << "sharing far detector " << iii << " " << jjj << endl;
              for (Int_t i = 0; i < n_evis_bins; i++){
                Int_t j = i; //Only j==i matters

                float factor=tdper[istage].MuonVetoEff[idet]
                  *tdper[istage].DMCEff[idet]
                  *tdper[istage].Livetime[idet]
                  *tdper[istage].TargetMass[idet]/tdper[istage].TargetMass[0];


                // Pseudo Pearson chi2
                Double_t err_i = 0;
                Double_t err_j = 0;

                if (factor > 0){
                  err_i = sqrt((predper->GetPred(istage,idet,idet2)->GetBinContent(i+1) + tdper[istage].CorrBgEvtsSpec[idet]->GetBinContent(i+1))*factor)/factor;
                  err_j = sqrt((predper->GetPred(jstage,jdet,jdet2)->GetBinContent(j+1) + tdper[jstage].CorrBgEvtsSpec[jdet]->GetBinContent(j+1))*factor)/factor;

                  M_stat[(MaxPredictions*istage+iii)*n_evis_bins+i][(MaxPredictions*jstage+jjj)*n_evis_bins+j] += err_i * err_j / stat_factor[istage];

                  if (predper->GetPred(istage,idet,idet2)->GetBinContent(i+1) != 0 && predper->GetPred(jstage,jdet,jdet2)->GetBinContent(j+1) != 0)
                    M_stat_frac[(MaxPredictions*istage+iii)*n_evis_bins+i][(MaxPredictions*jstage+jjj)*n_evis_bins+j] += err_i * err_j /predper->GetPred(istage,idet,idet2)->GetBinContent(i+1)/predper->GetPred(jstage,jdet,jdet2)->GetBinContent(j+1);

                }
                // Neyman chi2
                // to cross-check
                // Double_t err_i = tperdata->CorrEvtsSpec[idet]->GetBinError(i+1);
                // Double_t err_j = tperdata->CorrEvtsSpec[jdet]->GetBinError(j+1);



                //cout << "M_stat: " <<  M_stat[(MaxPredictions*istage+iii)*n_evis_bins+i][(MaxPredictions*jstage+jjj)*n_evis_bins+j] << endl;

              }
            }
          }
        }
      }
    }
  }
}


void Predictor::CalculateNearSiteStatError(){

  cout << "Calculating near site statistical error only...." << endl;
  for (int i = 0; i < MaxPredictions*Nstage*n_evis_bins; i++){
    for (int j = 0; j < MaxPredictions*Nstage*n_evis_bins; j++){
      M_stat[i][j] = 0;
      M_stat_frac[i][j] = 0;
    }
  }

  for(int istage=0;istage<Nstage;++istage){
    //for(int jstage=Nstage;jstage<Nstage;++jstage){
    int jstage = istage; //Only istage==jstage is important
    for(int idet=4;idet<8;++idet){
      for(int idet2=0;idet2<4;++idet2){
        for(int jdet=4;jdet<8;++jdet){
          for(int jdet2=0;jdet2<4;++jdet2){
            Int_t iii = (idet-4)*4 +idet2;
            Int_t jjj = (jdet-4)*4 +jdet2;

            if (idet2 == jdet2){ // sharing the same near detector during the same stage
              for (Int_t i = 0; i < n_evis_bins; i++){
                Int_t j = i; //Only i==j matters
                //Scale near site statistical error ( = sqrt(N_tot^2 + N_bg^2)) to the far site
                //"evtset-> GetEvts(idet,idet2) / CorrEvts[idet2]" gives the far/near ratio at given theta 13 value

                Double_t err_i = 0;
                if (tdper[istage].CorrEvtsSpec[idet2]->GetBinContent(i+1) > 0)
                  err_i = tdper[istage].CorrEvtsSpec[idet2]->GetBinError(i+1) * predper->GetPred(istage,idet,idet2)->GetBinContent(i+1) / tdper[istage].CorrEvtsSpec[idet2]->GetBinContent(i+1);
                Double_t err_j = 0;
                if (tdper[jstage].CorrEvtsSpec[jdet2]->GetBinContent(j+1) > 0)
                  err_j = tdper[jstage].CorrEvtsSpec[jdet2]->GetBinError(j+1) * predper->GetPred(jstage,jdet,jdet2)->GetBinContent(j+1) / tdper[jstage].CorrEvtsSpec[jdet2]->GetBinContent(j+1);


                M_stat[(MaxPredictions*istage+iii)*n_evis_bins+i][(MaxPredictions*jstage+jjj)*n_evis_bins+j] += err_i * err_j / stat_factor[istage];

                M_stat_frac[(MaxPredictions*istage+iii)*n_evis_bins+i][(MaxPredictions*jstage+jjj)*n_evis_bins+j] += err_i * err_j /predper->GetPred(istage,idet,idet2)->GetBinContent(i+1)/predper->GetPred(jstage,jdet,jdet2)->GetBinContent(j+1);
              }
            }
          }
        }
      }
    }
  }
}

Double_t * Predictor::GetFinalCovMatrix(Int_t type, Int_t mode){

  if (RecalculateCovMatrix){
    AddandScaleCovMatrix(type); // make matrix with a subset of error for display purpose
  }

  return CombineMatrix(mode,false);

}


Double_t * Predictor::GetFinalCovMatrixSum(Int_t type, Int_t mode){

  if (RecalculateCovMatrix){
    AddandScaleCovMatrix(type); // make matrix with a subset of error for display purpose
  }

  return CombineMatrix(mode,true);

}


Double_t * Predictor::CombineMatrix(Int_t mode, Bool_t MakeSum){

  Int_t npredictions = 2;

  if (CalculateMode1Coeff == true){

    for (Int_t istage=0;istage<Nstage;istage++){
      mode1_coeff[istage][0] = (double)NdetectorsConfig[istage][0] / (double)(NdetectorsConfig[istage][0] + NdetectorsConfig[istage][1]) ;
      mode1_coeff[istage][1] = (double)NdetectorsConfig[istage][1] / (double)(NdetectorsConfig[istage][0] + NdetectorsConfig[istage][1]) ;
      mode1_coeff[istage][2] = (double)NdetectorsConfig[istage][2] / (double)(NdetectorsConfig[istage][2] + NdetectorsConfig[istage][3]) ;
      mode1_coeff[istage][3] = (double)NdetectorsConfig[istage][3] / (double)(NdetectorsConfig[istage][2] + NdetectorsConfig[istage][3]) ;
    }
    CalculateMode1Coeff = false;
  }

  if (CalculateNearEHCoeff == true){

    for (Int_t istage=0;istage<Nstage;istage++){

      //Count total predicted IBD in Near EH during a single period
      Double_t TotalNearEH_IBD = 0;

      for(Int_t idet=0;idet<4;++idet)//Sum over 4 near ADs
        TotalNearEH_IBD += nIBD[istage][idet];

      //Now calculate weighting
      for(Int_t idet=0;idet<4;++idet)
        NearEH_coeff[istage][idet] = nIBD[istage][idet] / TotalNearEH_IBD;

    }//End of istage

    CalculateNearEHCoeff = false;
  }

  if (CalculateStageCoeff == true){

    Double_t TotalFarEH_IBD = 0; //Keep track of total EH3 IBD for all periods
    Double_t FarEH_IBD[Nstage]; //Keep track of total EH3 IBD for a single period
    for (Int_t istage=0;istage<Nstage;istage++)
      FarEH_IBD[istage] = 0;

    //First let's calculate weighting for individual Far ADs
    for(Int_t idet=4;idet<Ndetectors;++idet){

      Double_t TotalFarAD_IBD = 0;
      for (Int_t istage=0;istage<Nstage;istage++){
        TotalFarAD_IBD += nIBD[istage][idet];
        TotalFarEH_IBD += nIBD[istage][idet];
        FarEH_IBD[istage] += nIBD[istage][idet];
      }

      //Now do weighting
      for (Int_t istage=0;istage<Nstage;istage++){
        FarAD_stage_coeff[istage][idet] = nIBD[istage][idet] / TotalFarAD_IBD;
        cout << "idet: " << idet << ", istage: " << istage << ", nIBD: " << nIBD[istage][idet] << ", FarAD: " << FarAD_stage_coeff[istage][idet] << endl;
      }

    }

    //Now that we have total EH3 IBD for all period, can calculate EH3 period weighting
    for (Int_t istage=0;istage<Nstage;istage++){
      for(Int_t idet=4;idet<Ndetectors;++idet){
        FarEH_stage_coeff[istage][idet] = FarEH_IBD[istage] / TotalFarEH_IBD; //No dependence on det

        cout << "istage: " << istage << ", idet: " << idet << ", nIBD: " << nIBD[istage][idet] << ", FarEH: " << FarEH_stage_coeff[istage][idet] << endl;
      }
    }

    CalculateStageCoeff = false;
  }

  if (mode == 0)
    npredictions = 16;
  if (mode == 1)
    npredictions = 2;
  if (mode == 2)
    npredictions = 4;
  if (mode == 3)
    npredictions = 1;

  for (Int_t i = 0; i < MaxPredictions*Nstage*max_n_evis_bins; i++){
    for (Int_t j = 0; j < MaxPredictions*Nstage*max_n_evis_bins; j++){
      M_scaled[i][j] = 0;
    }
  }

  for(int istage=0;istage<Nstage;++istage){
    for(int jstage=0;jstage<Nstage;++jstage){
      for(int idet=4;idet<8;++idet){
        for(int idet2=0;idet2<4;++idet2){
          for(int jdet=4;jdet<8;++jdet){
            for(int jdet2=0;jdet2<4;++jdet2){
              Int_t iii = (idet-4)*4 +idet2;
              Int_t jjj = (jdet-4)*4 +jdet2;

              for(int ie=0;ie<n_evis_bins;++ie){
                for(int je=0;je<n_evis_bins;++je){

                  if (MakeSum){

                    if (mode == 0){// no combining
                      M_scaled[(iii)*n_evis_bins_rebin+evis_rebin_map[ie]][(jjj)*n_evis_bins_rebin+evis_rebin_map[je]]
                        += FarAD_stage_coeff[istage][idet]* FarAD_stage_coeff[jstage][jdet]* M[(MaxPredictions*istage+iii)*n_evis_bins+ie][(MaxPredictions*jstage+jjj)*n_evis_bins+je];
                    }
                    if (mode == 1){ //Combine to 2x2 matrix
                      M_scaled[(detConfigEH[idet2]-1)*n_evis_bins_rebin+evis_rebin_map[ie]][(detConfigEH[jdet2]-1)*n_evis_bins_rebin+evis_rebin_map[je]]
                        +=FarEH_stage_coeff[istage][idet] * FarEH_stage_coeff[jstage][jdet] * mode1_coeff[istage][idet2] * mode1_coeff[jstage][jdet2] * M[(MaxPredictions*istage+iii)*n_evis_bins+ie][(MaxPredictions*jstage+jjj)*n_evis_bins+je];
                    }
                    if (mode == 2){ // Combine to 4x4 matrix
                      M_scaled[(idet2)*n_evis_bins_rebin+evis_rebin_map[ie]][(jdet2)*n_evis_bins_rebin+evis_rebin_map[je]]
                        +=FarAD_stage_coeff[istage][idet] * FarAD_stage_coeff[jstage][jdet] * M[(MaxPredictions*istage+iii)*n_evis_bins+ie][(MaxPredictions*jstage+jjj)*n_evis_bins+je];
                    }
                    if (mode == 3){
                      M_scaled[evis_rebin_map[ie]][evis_rebin_map[je]] +=FarEH_stage_coeff[istage][idet] * FarEH_stage_coeff[jstage][jdet] * NearEH_coeff[istage][idet2] * NearEH_coeff[jstage][jdet2] * M[(MaxPredictions*istage+iii)*n_evis_bins+ie][(MaxPredictions*jstage+jjj)*n_evis_bins+je];
                    }
                  }

                  else{
                    if (mode == 0){// no combining
                      M_scaled[(npredictions*istage+iii)*n_evis_bins_rebin+evis_rebin_map[ie]][(npredictions*jstage+jjj)*n_evis_bins_rebin+evis_rebin_map[je]]
                        += M[(MaxPredictions*istage+iii)*n_evis_bins+ie][(MaxPredictions*jstage+jjj)*n_evis_bins+je];
                    }
                    if (mode == 1){ //Combine to 2x2 matrix
                      M_scaled[(npredictions*istage+(detConfigEH[idet2]-1))*n_evis_bins_rebin+evis_rebin_map[ie]][(npredictions*jstage+(detConfigEH[jdet2]-1))*n_evis_bins_rebin+evis_rebin_map[je]]
                        += mode1_coeff[istage][idet2] * mode1_coeff[jstage][jdet2] * M[(MaxPredictions*istage+iii)*n_evis_bins+ie][(MaxPredictions*jstage+jjj)*n_evis_bins+je];


                      if(n_evis_bins_rebin==37 && n_evis_bins==37){
                        //cout<<"Adding exra uncertainty"<<endl;
                        if(ie==je){
                          if(ie<3){
                            M_scaled[(npredictions*istage+(detConfigEH[idet2]-1))*n_evis_bins_rebin+evis_rebin_map[ie]][(npredictions*jstage+(detConfigEH[jdet2]-1))*n_evis_bins_rebin+evis_rebin_map[je]]=1.4*1.4*M_scaled[(npredictions*istage+(detConfigEH[idet2]-1))*n_evis_bins_rebin+evis_rebin_map[ie]][(npredictions*jstage+(detConfigEH[jdet2]-1))*n_evis_bins_rebin+evis_rebin_map[je]];
                          }
                        }
                      }

                    }
                    if (mode == 2){ // Combine to 4x4 matrix
                      M_scaled[(npredictions*istage+idet2)*n_evis_bins_rebin+evis_rebin_map[ie]][(npredictions*jstage+jdet2)*n_evis_bins_rebin+evis_rebin_map[je]]
                        += M[(MaxPredictions*istage+iii)*n_evis_bins+ie][(MaxPredictions*jstage+jjj)*n_evis_bins+je];
                    }
                    if (mode == 3){ // Combine to 1x1 matrix
                      M_scaled[(npredictions*istage)*n_evis_bins_rebin+evis_rebin_map[ie]][(npredictions*jstage)*n_evis_bins_rebin+evis_rebin_map[je]]
                        += NearEH_coeff[istage][idet2] * NearEH_coeff[jstage][jdet2] * M[(MaxPredictions*istage+iii)*n_evis_bins+ie][(MaxPredictions*jstage+jjj)*n_evis_bins+je];
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }


  if (MakeSum){
    for (Int_t i = 0; i < npredictions*n_evis_bins_rebin; i++){
      for (Int_t j = 0; j < npredictions*n_evis_bins_rebin; j++){
        final_covmatrix[i*npredictions*n_evis_bins_rebin+j] = M_scaled[i][j];
      }
    }
  }else{
    for (Int_t i = 0; i < npredictions*Nstage*n_evis_bins_rebin; i++){
      for (Int_t j = 0; j < npredictions*Nstage*n_evis_bins_rebin; j++){
        final_covmatrix[i*npredictions*Nstage*n_evis_bins_rebin+j] = M_scaled[i][j];
      }
    }

  }

  return &final_covmatrix[0];

}

//Double_t  Predictor::GetWmeanCoeff(Int_t istage, Int_t i, Int_t idet){
//  return wmean_coeff[istage][i][idet];
//}
Double_t Predictor::GetFracStatError(Int_t i, Int_t j){
  return M_stat_frac[i][j];
}
Double_t Predictor::GetFracBgError(Int_t i, Int_t j){
  return M_bg_sys_frac[i][j];
}
