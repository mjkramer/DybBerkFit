//#include <iostream>
//#include "TROOT.h"

void run_fit_shape_2d_P15A_toyMC(Int_t input = 0){

  gROOT->ProcessLine(".x LoadClasses.C");
  gROOT->ProcessLine(".L fit_shape_2d_P17B_toyMC.C+");

  //TString bg_options[nopts] = ;
  cout << "Start" << endl;

 int itoy=input%3;
    cout<<itoy<<endl;
 int iAnalysis=(int)(input-itoy)/3;


 TString Analysis[3] = {"LBNL","BCW","IHEP"};

 Char_t sig_spectra_filename0[1024];
 Char_t sig_spectra_filename1[1024];
 Char_t sig_spectra_filename2[1024];
 Char_t acc_spectra_filename0[1024];
 Char_t acc_spectra_filename1[1024];
 Char_t acc_spectra_filename2[1024];
 Char_t input_filename0[1024];
 Char_t input_filename1[1024];
 Char_t input_filename2[1024];
 Char_t histogram_filename[1024];
 Char_t save_filename[1024];
 Char_t fn_filename[1024];
 
 sprintf(sig_spectra_filename0,"../outputs/ToyMCsComparisonP17B/%s_6AD_Oscillation%d_LBNLbin.root",Analysis[iAnalysis].Data(),itoy+1);
 sprintf(sig_spectra_filename1,"../outputs/ToyMCsComparisonP17B/%s_8AD_Oscillation%d_LBNLbin.root",Analysis[iAnalysis].Data(),itoy+1);
 sprintf(sig_spectra_filename2,"../outputs/ToyMCsComparisonP17B/%s_7AD_Oscillation%d_LBNLbin.root",Analysis[iAnalysis].Data(),itoy+1);
 
 sprintf(acc_spectra_filename0,"./Spectra/accidental_eprompt_shapes_6ad_LBNL.root");
 sprintf(acc_spectra_filename1,"./Spectra/accidental_eprompt_shapes_8ad_LBNL.root");
 sprintf(acc_spectra_filename2,"./Spectra/accidental_eprompt_shapes_7ad_LBNL.root");
 
 sprintf(input_filename0,"./Inputs/Theta13-inputs_P17B_inclusive_6ad_%s.txt",Analysis[iAnalysis].Data());
 sprintf(input_filename1,"./Inputs/Theta13-inputs_P17B_inclusive_8ad_%s.txt",Analysis[iAnalysis].Data());
 sprintf(input_filename2,"./Inputs/Theta13-inputs_P17B_inclusive_7ad_%s.txt",Analysis[iAnalysis].Data());

 sprintf(histogram_filename,"./Flux/SuperHistograms_P17B_2017Model_fine_huber-french.root");
	
 sprintf(save_filename,"./fit_result_files/toyMC_P17B_%s_Oscillation%d.root",Analysis[iAnalysis].Data(),itoy+1);

 sprintf(fn_filename,"../fn_spectrum/P15A_fn_spectrum.root");

 cout << "sig0: " << sig_spectra_filename0 << endl;
 cout << "sig1: " << sig_spectra_filename1 << endl;
 cout << "sig2: " << sig_spectra_filename2 << endl;
 cout << "acc0: " << acc_spectra_filename0 << endl;
 cout << "acc1: " << acc_spectra_filename1 << endl;
 cout << "acc2: " << acc_spectra_filename2 << endl;
 cout << "Input0: " << input_filename0 << endl;
 cout << "Input1: " << input_filename1 << endl;
 cout << "Input2: " << input_filename2 << endl;
 cout << "Histogram: " << histogram_filename << endl;
 cout << "savefile: " << save_filename << endl;
 cout << "fn: " << fn_filename << endl;
 

 fit_shape_2d_P17B_toyMC(sig_spectra_filename0,
			 sig_spectra_filename1,
             sig_spectra_filename2,
			 acc_spectra_filename0,
			 acc_spectra_filename1,
             acc_spectra_filename2,
			 input_filename0,
			 input_filename1,
             input_filename2,
			 histogram_filename,
			 save_filename,
			 fn_filename);
 
}
