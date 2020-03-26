#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include <TROOT.h>

using namespace std;

Bool_t AddSNF = true;
Bool_t AddNonEq = true;

Double_t snf_fraction = 0.003;

const Int_t nMaxIbdPoints = 1200;
Double_t ibdxsec[nMaxIbdPoints];
Double_t ibdenu[nMaxIbdPoints];
Int_t nIbdPoints;

// used for manually add off-equiliblium correction
Double_t f_u235 = 0.586;
Double_t f_p239 = 0.288;
Double_t f_p241 = 0.050;

Double_t e_noneq[5] = {2,2.5,3.0,3.5, 4.0};
// Double_t noneq_u235[5] = {5.7, 4.4, 1.5, 0.7, 0.1};
// Double_t noneq_p239[5] = {2.1, 1.7, 0.5, 0.0, 0.0};
// Double_t noneq_p241[5] = {1.9, 1.5, 0.5, 0.0, 0.0};

Double_t noneq_corr[4][5] = {{5.7, 4.4, 1.5, 0.7, 0.1},
                             {0.0, 0.0, 0.0, 0.0, 0.0},
                             {2.1, 1.7, 0.5, 0.0, 0.0},
                             {1.9, 1.5, 0.5, 0.0, 0.0}};
Double_t noneq_ave[5];

TString isotope[4] = {"U235","U238","Pu239","Pu241"};

Int_t load_ibd_xsec(TString filename){
  std::string line;
  Double_t enu_tmp;
  Double_t xsec_tmp;

  ifstream fin(filename.Data());

  nIbdPoints = 0;
  while(nMaxIbdPoints > nIbdPoints){
    if(fin.peek()=='#'){
      // Skip lines starting with '#'
      getline(fin,line);
      continue;
    }else{
      fin >> enu_tmp >> xsec_tmp;
    }
    if(!fin.good()) break;
    ibdenu[nIbdPoints] = enu_tmp;
    ibdxsec[nIbdPoints] = xsec_tmp;
    nIbdPoints++;
  }
  fin.close();
  return nIbdPoints;
}

Double_t get_ibd_xsec(Double_t enu){

  Double_t xsec = 0;
  if (enu >= ibdenu[nIbdPoints-1])
    xsec =  ibdxsec[nIbdPoints-1];
  else{
    for (Int_t i = 0; i < nIbdPoints-1; i++){
      if (enu >= ibdenu[i] && enu < ibdenu[i+1]){
        xsec = ((ibdenu[i+1] - enu) * ibdxsec[i] + (enu - ibdenu[i]) * ibdxsec[i+1])
          / (ibdenu[i+1] - ibdenu[i]);
        break;
      }
    }
  }
  return xsec;
}

Double_t get_noneq_corr(Int_t isotope_id, Double_t e){
  Double_t corr = 0;

  if (e > 4.0) return corr;

  Int_t id = (Int_t)((e - 2.0)/0.5);
  if (id < 0) id = 0;

  corr = ((e_noneq[id+1] - e)*noneq_corr[isotope_id][id] + (e-e_noneq[id])*noneq_corr[isotope_id][id+1])
    /(e_noneq[id+1]-e_noneq[id]);

  return 0.01*corr;
}

void make_combined_spectra() {


  load_ibd_xsec("../../toySpectra/reactor/Xsec1_2011.dat");

  Double_t flux_w_snf[6][220];
  Double_t flux_wo_snf[6][220];
  Double_t flux_snf[6][220];

  Double_t snf_percentage[6][220];

  Double_t flux_wo_snf_total[6];
  Double_t flux_wo_snf_total_ave = 0;
  Double_t flux_snf_total[6];
  Double_t enu[220];
  Double_t enu_tmp;
  Double_t enu_tmp2;

  Double_t flux_tmp;
  Double_t flux_tmp_corr;

  Int_t Nstage=3;

  //TString stage_name[5] = {"P15A_6AD","P15A_P14Aperiod_8AD","P15A_full_blinded_8AD"};
  TString stage_name[3] = {"P17B_6AD","P17B_8AD","P17B_7AD"};

  for (Int_t istage = 0; istage < Nstage; istage++){
    for (Int_t j = 0; j < 6; j++){
      flux_wo_snf_total[j] = 0;
      flux_snf_total[j] = 0;
      for (Int_t i = 0; i < 220; i ++){
        flux_wo_snf[j][i] = 0;
      }
    }


    for (Int_t ii = 0; ii < 4; ii++){
      ifstream fin(Form("reactor_%s_%s_unblinded.txt",stage_name[istage].Data(),isotope[ii].Data()));
      for (Int_t i = 0; i < 220; i ++){
        fin >> enu_tmp;
        enu[i] = enu_tmp;
        //	cout << enu_tmp << " " << get_noneq_corr(ii,enu_tmp) << endl;
        for (Int_t j = 0; j < 6; j++){
          fin >> flux_tmp;
          if (AddNonEq)
            flux_tmp_corr = flux_tmp * (1 + get_noneq_corr(ii,enu_tmp));
          else
            flux_tmp_corr = flux_tmp;
          flux_wo_snf[j][i] += flux_tmp_corr;
          flux_wo_snf_total[j] += flux_tmp_corr * get_ibd_xsec(enu_tmp);
          flux_wo_snf_total_ave += flux_tmp_corr * get_ibd_xsec(enu_tmp);
        }
      }
      fin.close();
    }
    flux_wo_snf_total_ave /= 6.0;

    // Read SNF spectra
    //additioanl percentage of SNF is taken from DocDB-11687
    double snf_coarse_energy[11];
    for(int ien=0; ien<11;++ien) snf_coarse_energy[ien]=1.75+0.25*ien;

    double snf_coarse_percentage[6][11]={{1.169,1.317,1.329,0.9747,0.8996,0.7766,0.3185,0.2475,0.03444,0.0001644,0.},{1.186,1.334,1.345,0.9816,0.9050,0.7812,0.3200,0.2487,0.03460,0.0001657,0.},{1.196,1.367,1.412,1.118,1.043,0.9007,0.3691,0.2870,0.03996,0.0002051,0.},{1.158,1.321,1.362,1.077,1.003,0.8666,0.3560,0.2769,0.003857,0.0002021,0.},{1.087,1.265,1.342,1.175,1.107,0.9545,0.3854,0.2998,0.04178,0.0002208,0.},{0.9552,1.115,1.184,1.009,0.9525,0.8192,0.3258,0.2532,0.03527,0.0002063}}; //in percent


    TGraph* gr_snf_percentage[6];
    TCanvas* c_snf=new TCanvas("c_snf","c_snf",1000,500);
    c_snf->Divide(3,2);

    for(int icore=0; icore<6; ++icore){
      gr_snf_percentage[icore]=new TGraph();
      for(int ienfine=0; ienfine<220; ++ienfine){
        double fine_energy=1.85+0.05*ienfine;
        int coarse_energy_bin=-1;

        for(int ien=0; ien<11;++ien){
          if(fine_energy<snf_coarse_energy[ien]){
            coarse_energy_bin=ien-1;
            break;
          }
        }

        cout<<"coarse_energy_bin is "<<coarse_energy_bin<<endl;
        if(coarse_energy_bin==-1) snf_percentage[icore][ienfine]=0.;
        else{
          snf_percentage[icore][ienfine]=snf_coarse_percentage[icore][coarse_energy_bin]+(fine_energy-snf_coarse_energy[coarse_energy_bin])*(snf_coarse_percentage[icore][coarse_energy_bin+1]-snf_coarse_percentage[icore][coarse_energy_bin])/(snf_coarse_energy[coarse_energy_bin+1]-snf_coarse_energy[coarse_energy_bin]);
        }
        gr_snf_percentage[icore]->SetPoint(ienfine,fine_energy,snf_percentage[icore][ienfine]);
      }
      c_snf->cd(icore+1);
      gr_snf_percentage[icore]->Draw("AL");
    }




    ofstream fout(Form("reactor_%s_SNF_nonEq.txt",stage_name[istage].Data()));


    //  cout << 1./flux_snf_total[2]*flux_wo_snf_total_ave << endl;
    for (Int_t i = 0; i < 220; i ++){
      fout <<  std::setprecision(3)  << enu[i]  << std::setprecision(12);
      cout << enu[i];
      for (Int_t j = 0; j < 6; j++){
        if (AddSNF)
          flux_w_snf[j][i] = flux_wo_snf[j][i]*(1.+snf_percentage[j][i]/100.);  // /100. due to values are in percentage
        else
          flux_w_snf[j][i] = flux_wo_snf[j][i];


        //cout <<  "\t" << snf_fraction * flux_snf[0][i] /flux_snf_total[0]*flux_wo_snf_total_ave;
        cout << "\t" << flux_w_snf[j][i];
        fout << "\t" << flux_w_snf[j][i];

      }
      cout << endl;
      fout << endl;
    }
  }

}
