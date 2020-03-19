#include <iostream>
#include "TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TPad.h"

using namespace std;


const int Nstages=3;
const int Ncores=6;

//6AD 24.12.2011-28.7.2012
//8AD 19.10.2012-28.7.2015
//7AD after

const float MeV_per_s_per_GW=6.2415E21;

const int start_week[Nstages]={0,43,265};
const int end_week[Nstages]={31,261,297};  //including this week
const int stage_id[3]={6,8,7}; //6AD 8AD later add 7AD

float nominal_power[Nstages][Ncores];
float average_power[Nstages][Ncores];

float average_fission_frac[Nstages][Ncores][4]; //for isotopes 0-U235, 1-U238, 2-Pu239, 3-Pu241

float average_power_times_fission_frac[Nstages][Ncores][4]; //for isotopes 0-U235, 1-U238, 2-Pu239, 3-Pu241

float MeV_per_fission[4];

float weekly_livetime[300];

int week_counter=0;

void LoadLiveTime(){
  FILE* spectrum_file=fopen("dbd_livetime_P17B.txt","r");
  int day_counter=0;
    
  int iday=-1;
  int ieh=-1;
  float temp_livetime=-1.;
  float weekly_average_livetime=0.;
  int real_dat_counter=0;
  bool do_read=true;
  while(1){
    float average_day_livetime=0.;
        
    if(do_read){
      fscanf(spectrum_file,"%d %d %f",&iday,&ieh,&temp_livetime);
      average_day_livetime+=temp_livetime;
      if(ieh!=1) cout<<"Problem with hall"<<endl;
      fscanf(spectrum_file,"%d %d %f",&iday,&ieh,&temp_livetime);
      average_day_livetime+=temp_livetime;
      if(ieh!=2) cout<<"Problem with hall"<<endl;
      fscanf(spectrum_file,"%d %d %f",&iday,&ieh,&temp_livetime);
      if(ieh!=3) cout<<"Problem with hall"<<endl;
            

            
      //cout<<iday<<" "<<day_counter<<" "<<week_counter<<" "<<week_counter*7-iday<<endl;
            
      average_day_livetime+=temp_livetime;
      average_day_livetime=average_day_livetime/3.;
    }
    if(iday==2076){
      weekly_livetime[week_counter]=weekly_average_livetime/(7.*86400.);
      cout<<week_counter<<" "<<weekly_livetime[week_counter]<<endl;
      week_counter++;
      weekly_average_livetime=0;
      break;
    }
        
    if(real_dat_counter!=iday){
      //cout<<"   "<<weekly_average_livetime<<endl;
      weekly_average_livetime+=0.;
      cout<<"   "<<weekly_average_livetime<<endl;
      day_counter++;
      do_read=false;
    }
    else{
      weekly_average_livetime+=average_day_livetime;
      average_day_livetime=0.;
      day_counter++;
      do_read=true;
    }
    //cout<<day_counter<<" "<<iday+1<<" "<<real_dat_counter<<endl;
        
        
    if(day_counter%7==0){
      cout<<"   "<<weekly_average_livetime<<endl;
      weekly_livetime[week_counter]=weekly_average_livetime/(7.*86400.);
      cout<<week_counter<<" "<<weekly_livetime[week_counter]<<endl;
      weekly_average_livetime=0;
      week_counter++;
      day_counter=0;
    }
        
    real_dat_counter++;
  }
    
  for(int iweek=0; iweek<week_counter; ++iweek) cout<<"Weekly livetime for week "<<iweek<<" is "<<weekly_livetime[iweek]<<endl;
}


//TGraph* gr_henoch_spectra[Nstages][Ncores][4];
float array_henoch_spectra[Nstages][Ncores][4][220];
float energy_steps[220];


void ReadHenochSpectra(){
  for(int istage=0; istage<1;++istage){
    for(int iisotope=0; iisotope<4; ++iisotope){
      //for(int icore=0; icore<Ncores;++icore) gr_henoch_spectra[istage][icore][iisotope]=new TGraph();
      int ipoint=0;
      char filepath[64];
      if(iisotope==0) sprintf(filepath,"./LBNL_reac_cov/isotope_spectra/reactor_P15A_%dAD_U235.txt",stage_id[istage]);
      if(iisotope==1) sprintf(filepath,"./LBNL_reac_cov/isotope_spectra/reactor_P15A_%dAD_U238.txt",stage_id[istage]);
      if(iisotope==2) sprintf(filepath,"./LBNL_reac_cov/isotope_spectra/reactor_P15A_%dAD_Pu239.txt",stage_id[istage]);
      if(iisotope==3) sprintf(filepath,"./LBNL_reac_cov/isotope_spectra/reactor_P15A_%dAD_Pu241.txt",stage_id[istage]);
      cout<<filepath<<endl;
      FILE* spectrum_file=fopen(filepath,"r");
      if(spectrum_file == NULL) cout<<"file's not open"<<endl;
      float temp_energy=-1;
      float temp_spect[Ncores];
      while(1){
        fscanf(spectrum_file,"%f %f %f %f %f %f %f",&temp_energy,&temp_spect[0],&temp_spect[1],&temp_spect[2],&temp_spect[3],&temp_spect[4],&temp_spect[5]);
        if(feof(spectrum_file)) break;
        for(int icore=0; icore<Ncores;++icore){
          //gr_henoch_spectra[istage][icore][iisotope]->SetPoint(ipoint,temp_energy,temp_spect[icore]);
          array_henoch_spectra[istage][icore][iisotope][ipoint]=temp_spect[icore];
          energy_steps[ipoint]=temp_energy;
        }
        //cout<<temp_energy<<" "<<temp_spect[0]<<endl;
        ipoint++;
      }
      cout<<"There are "<<ipoint<<" points in Henoch's data"<<endl;
    }
  }
}

float spectra[220][4];

int Npoints_in_spectra=0;

void ReadCoreSpectra(int option=1){
  char filepath[64];
  if(option==0) sprintf(filepath,"./BCW/fissionIsotopeSpectra_Huber_Linear.txt");
  if(option==1) sprintf(filepath,"./BCW/fissionIsotopeSpectra_Huber_Fit.txt");
  if(option==2) sprintf(filepath,"./BCW/fissionIsotopeSpectra_ILL_Muller_Linear.txt");
  if(option==3) sprintf(filepath,"./BCW/fissionIsotopeSpectra_ILL_Muller_Fit.txt");
  if(option==4) sprintf(filepath,"./LBNL_Spectra/fissionIsotopeSpectra_Huber_v0.txt");
  cout<<filepath<<endl;
    
  FILE* input_spectrum_file=fopen(filepath,"r");
  if(input_spectrum_file == NULL) cout<<"file's not open"<<endl;
  float temp_energy=-1;
  float temp_spect[4];
    
  char dummy_str4[128];
  /*if(option==0)*/ for(int i=0; i<8;++i) fgets(dummy_str4,128, input_spectrum_file);
  //cout<<dummy_str4<<endl;
  int counter=0;
  int ispec=0;
  while(1){
    fscanf(input_spectrum_file,"%f %f %f %f %f",&temp_energy,&temp_spect[0],&temp_spect[1],&temp_spect[2],&temp_spect[3]);
    if(feof(input_spectrum_file)) break;
    if(temp_energy>1.8 && counter%5==0 && temp_energy<=12.8){
      //if(ispec==0) cout<<temp_energy<<endl;
      for(int iisotope=0; iisotope<4; ++iisotope) spectra[ispec][iisotope]=temp_spect[iisotope];
      Npoints_in_spectra++;
      ispec++;
    }
    counter++;
  }
}

float frac_power_array[Ncores][300];
float f_frac_array[Ncores][300][4]; //for isotopes 0-U235, 1-U238, 2-Pu239, 3-Pu241
float f_average_energy_per_fission[Ncores][300]; //for isotopes 0-U235, 1-U238, 2-Pu239, 3-Pu241

int nweeks_unblinded=0;

void ReadWeeklyAvg(){
  LoadLiveTime();
  for(int istage=0; istage<Nstages; ++istage){
    for(int icore=0; icore<Ncores;++icore){
      nominal_power[istage][icore]=2.9; //GW
      average_power[istage][icore]=0.;
      for(int iisotope=0; iisotope<4; ++iisotope){
        average_fission_frac[istage][icore][iisotope]=0.;
        average_power_times_fission_frac[istage][icore][iisotope]=0.;
      }
    }
  }
  FILE* energy_per_fission=fopen("./fissionIsotopeTable_v1.txt","r");
  char dummy_str3[128];
  for(int i=0; i<6;++i) fgets(dummy_str3,128, energy_per_fission);
  int isotope_num=-1;
  float temp_MeV_per_fission;
  while(1){
    fscanf(energy_per_fission,"%d %f %*f",&isotope_num,&temp_MeV_per_fission);
    if(feof(energy_per_fission)) break;
    MeV_per_fission[isotope_num-1]=temp_MeV_per_fission;
    cout<<isotope_num-1<<" "<<MeV_per_fission[isotope_num-1]<<endl;
  }
  fclose(energy_per_fission);
    
    
  FILE* power_fission_frac_file=fopen("./WeeklyAvg/WeeklyAvg_P15A_v1_from_Henoch.txt","r");
  int week=-1;
  int core=-1;
  char dummy_str[64];
  char dummy_str2[64];
  int start_sec;
  int end_sec;
  float frac_power=-1.;
  float f_frac_u235=-1.;
  float f_frac_u238=-1.;
  float f_frac_pu239=-1.;
  float f_frac_pu241=-1.;

  while(1){
    fscanf(power_fission_frac_file,"%d %d %s %s %d %d",&week,&core,dummy_str,dummy_str2,&start_sec,&end_sec);

    fscanf(power_fission_frac_file,"%f %f %f %f %f",&frac_power,&f_frac_u235,&f_frac_u238,&f_frac_pu239,&f_frac_pu241);
    //cout<<dummy_str<<" "<<dummy_str2<<" "<<start_sec<<endl;
        
        
    if(feof(power_fission_frac_file)) break;
    //cout<<week<<" "<<core<<" "<<frac_power<<" "<<f_frac_u235<<" "<<f_frac_u238<<" "<<f_frac_pu239<<" "<<f_frac_pu241<<endl;
    frac_power_array[core-1][week]=frac_power;
        
    f_frac_array[core-1][week][0]=f_frac_u235;
    f_frac_array[core-1][week][1]=f_frac_u238;
    f_frac_array[core-1][week][2]=f_frac_pu239;
    f_frac_array[core-1][week][3]=f_frac_pu241;
        
    f_average_energy_per_fission[core-1][week]=0.;
    for(int iisotope=0; iisotope<4; ++iisotope) f_average_energy_per_fission[core-1][week]+=MeV_per_fission[iisotope]*f_frac_array[core-1][week][iisotope];
    cout<<"Average energy per fission "<<f_average_energy_per_fission[core-1][week]<<endl;
    /*f_frac_array[core-1][week][0]=0.563;   //fixing fission fractions
      f_frac_array[core-1][week][1]=0.076;
      f_frac_array[core-1][week][2]=0.295;
      f_frac_array[core-1][week][3]=0.058;*/
    if(core==1) nweeks_unblinded++; //count only once per list of cores
  }

    
  fclose(power_fission_frac_file);
    
  for(int istage=0; istage<Nstages; ++istage){
    float tot_stage_livetime=0.;
    for(int icore=0; icore<Ncores; ++icore){
      for(int iweek=0; iweek<nweeks_unblinded; ++iweek){
        if(iweek>=start_week[istage] && iweek<=end_week[istage]){
          //weekly_livetime[iweek]=1.;
          if(icore==0) tot_stage_livetime+=weekly_livetime[iweek];
          average_power[istage][icore]+=frac_power_array[icore][iweek]*weekly_livetime[iweek];
          for(int iisotope=0; iisotope<4; ++iisotope){
            average_fission_frac[istage][icore][iisotope]+=f_frac_array[icore][iweek][iisotope]*weekly_livetime[iweek];
            average_power_times_fission_frac[istage][icore][iisotope]+=f_frac_array[icore][iweek][iisotope]*frac_power_array[icore][iweek]*weekly_livetime[iweek];
          }
        }
      }
      average_power[istage][icore]=average_power[istage][icore]/tot_stage_livetime; //both start and end weeks are included
      for(int iisotope=0; iisotope<4; ++iisotope){
        average_fission_frac[istage][icore][iisotope]=average_fission_frac[istage][icore][iisotope]/tot_stage_livetime;
        average_power_times_fission_frac[istage][icore][iisotope]=average_power_times_fission_frac[istage][icore][iisotope]/tot_stage_livetime;
      }
    }
        
  }
    
  for(int istage=0; istage<Nstages; ++istage){
    cout<<"Stage "<<istage+1<<":"<<endl;
    cout<<"Average power (fraction of nominal): ";
    for(int icore=0; icore<Ncores; ++icore){
      cout<<average_power[istage][icore]<<" ";
    }
    cout<<endl;
        
    for(int iisotope=0; iisotope<4; ++iisotope){
      if(iisotope==0) cout<<"Average U235 fission fraction: ";
      if(iisotope==1) cout<<"Average U238 fission fraction: ";
      if(iisotope==2) cout<<"Average Pu239 fission fraction: ";
      if(iisotope==3) cout<<"Average Pu241 fission fraction: ";
      for(int icore=0; icore<Ncores; ++icore){
        cout<<average_fission_frac[istage][icore][iisotope]<<" ";
      }
      cout<<endl;
    }
    cout<<endl;
    for(int icore=0; icore<Ncores; ++icore){
      float control_sum=0.;
      for(int iisotope=0; iisotope<4; ++iisotope){
        control_sum+=average_fission_frac[istage][icore][iisotope];
      }
      cout<<control_sum<<" ";
    }
    cout<<endl<<endl;
        
  }
}


float array_beda_spectra[Nstages][Ncores][4][220];


void ProduceSpectraDwyer(int option=0){
  ReadCoreSpectra(option);
  ReadWeeklyAvg();
    
  for(int istage=0; istage<Nstages; ++istage){
    for(int icore=0; icore<Ncores; ++icore){
      for(int iisotope=0; iisotope<4; ++iisotope){
        for(int ipoint=0; ipoint<Npoints_in_spectra; ++ipoint){
          array_beda_spectra[istage][icore][iisotope][ipoint]=0;
        }
      }
    }
  }
    
  float weekly_livetime_sum[Nstages];
  for(int istage=0; istage<Nstages; ++istage) weekly_livetime_sum[istage]=0.;
    
  cout<<"There are "<<week_counter<<" weeks"<<endl;
    
  for(int istage=0; istage<Nstages; ++istage){
    for(int iweek=0; iweek<week_counter; ++iweek){
      if(iweek>=start_week[istage] && iweek<=end_week[istage]) weekly_livetime_sum[istage]+=weekly_livetime[iweek];
    }
  }
    
  for(int istage=0; istage<Nstages; ++istage){
    for(int icore=0; icore<Ncores; ++icore){
      for(int iweek=0; iweek<nweeks_unblinded; ++iweek){
        for(int ipoint=0; ipoint<Npoints_in_spectra; ++ipoint){ //from 1.85 to 12.8 MeV with 0.05 MeV step
          for(int iisotope=0; iisotope<4; ++iisotope){
            if(iweek>=start_week[istage] && iweek<=end_week[istage]){
              array_beda_spectra[istage][icore][iisotope][ipoint]+=weekly_livetime[iweek]*nominal_power[istage][icore]*frac_power_array[icore][iweek]*MeV_per_s_per_GW/f_average_energy_per_fission[icore][iweek]*f_frac_array[icore][iweek][iisotope]*spectra[ipoint][iisotope]/1.E18;
            }
          }
        }
      }
    }
  }
    
  float blinded_f_frac_array[4]={0.64,0.08,0.25,0.03};
  for(int istage=0; istage<Nstages; ++istage){ //for blinded flux
    for(int icore=0; icore<Ncores; ++icore){
      for(int iweek=nweeks_unblinded; iweek<week_counter; ++iweek){
        for(int ipoint=0; ipoint<Npoints_in_spectra; ++ipoint){ //from 1.85 to 12.8 MeV with 0.05 MeV step
          for(int iisotope=0; iisotope<4; ++iisotope){
            if(iweek>=start_week[istage] && iweek<=end_week[istage]){
              array_beda_spectra[istage][icore][iisotope][ipoint]+=weekly_livetime[iweek]*nominal_power[istage][icore]*1.*MeV_per_s_per_GW/204.5759*blinded_f_frac_array[iisotope]*spectra[ipoint][iisotope]/1.E18;
            }
          }
        }
      }
    }
  }
    
    
  for(int istage=0; istage<Nstages; ++istage) cout<<"Weekly livetime is "<<weekly_livetime_sum[istage]<<endl;
    
  for(int istage=0; istage<Nstages; ++istage){
    for(int icore=0; icore<Ncores; ++icore){
      for(int iisotope=0; iisotope<4; ++iisotope){
        for(int ipoint=0; ipoint<Npoints_in_spectra; ++ipoint){
          array_beda_spectra[istage][icore][iisotope][ipoint]=array_beda_spectra[istage][icore][iisotope][ipoint]/weekly_livetime_sum[istage];
        }
      }
    }
  }
}


TGraph* gr_ratio[6][4];

TGraph* gr_beda_spectra[6][4];
TGraph* gr_henoch_spectra[6][4];

void PrintSpectra(){
  ReadHenochSpectra(); //for energy binning
  ProduceSpectraDwyer();
    
  string period[Nstages]={"6AD","8AD","7AD"};
  string isotope[4]={"_U235","_U238","_Pu239","_Pu241"};
  string base="./isotope_spectra_by_Beda/reactor_P17B_";
  string txt=".txt";
  for(int istage=0; istage<Nstages;++istage){
    for(int iisotope=0; iisotope<4; ++iisotope){
      string filename=base+period[istage]+isotope[iisotope]+txt;
      cout<<filename<<endl;
      FILE* outputfile=fopen(filename.c_str(),"w");
      for(int ipoint=0; ipoint<Npoints_in_spectra; ++ipoint){
        fprintf(outputfile,"%.2f\t%f\t%f\t%f\t%f\t%f\t%f\n",energy_steps[ipoint],array_beda_spectra[istage][0][iisotope][ipoint],array_beda_spectra[istage][1][iisotope][ipoint],array_beda_spectra[istage][2][iisotope][ipoint],array_beda_spectra[istage][3][iisotope][ipoint],array_beda_spectra[istage][4][iisotope][ipoint],array_beda_spectra[istage][5][iisotope][ipoint]);
      }
      fclose(outputfile);
    }
  }
}


