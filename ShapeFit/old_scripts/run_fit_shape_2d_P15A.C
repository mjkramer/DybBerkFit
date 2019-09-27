//#include <iostream>
//#include "TROOT.h"

void run_fit_shape_2d_P15A(Int_t input = 1){

  gROOT->ProcessLine(".x LoadClasses.C");
  gROOT->ProcessLine(".L fit_shape_2d_P15A.C+");

  //TString bg_options[nopts] = ;
  cout << "Start" << endl;
  Int_t i = input-1; //convert to index starting with 0

  Int_t flag = int(i % 3)-1;
  
  Int_t iAnalysis = int(i/3);
  TString Analysis[6] = {"LBNL","BCW","IHEP","Dubna_AdSimple","Dubna_AdScaled","Tsinghua"};
  
 TString Period[3] = {"6+8AD","6AD","8AD"};
   
 Char_t sig_spectra_filename0[1024];
 Char_t sig_spectra_filename1[1024];
 Char_t acc_spectra_filename0[1024];
 Char_t acc_spectra_filename1[1024];
 Char_t input_filename0[1024];
 Char_t input_filename1[1024];
 Char_t fn_filename[1024];
 Char_t save_filename[1024];
 
 sprintf(sig_spectra_filename0,"./Spectra/ibd_eprompt_shapes_6ad_%s.root",Analysis[iAnalysis].Data());
 sprintf(sig_spectra_filename1,"./Spectra/ibd_eprompt_shapes_8ad_%s.root",Analysis[iAnalysis].Data());
 
 sprintf(acc_spectra_filename0,"./Spectra/accidental_eprompt_shapes_6ad_%s.root",Analysis[iAnalysis].Data());
 sprintf(acc_spectra_filename1,"./Spectra/accidental_eprompt_shapes_8ad_%s.root",Analysis[iAnalysis].Data());
 
 sprintf(input_filename0,"./Inputs/Theta13-inputs_P15A_inclusive_6ad_%s.txt",Analysis[iAnalysis].Data());
 sprintf(input_filename1,"./Inputs/Theta13-inputs_P15A_inclusive_8ad_%s.txt",Analysis[iAnalysis].Data());


 sprintf(fn_filename,"../fn_spectrum/P15A_fn_spectrum_%s.root",Analysis[iAnalysis].Data());

 //sprintf(save_filename,"./fit_result_files/fit_shape_2d_P15A_%s_%s.root",Analysis[iAnalysis].Data(),Period[flag+1].Data());
 sprintf(save_filename,"./fit_result_files/fit_shape_2d_P15A_stat_only_%s_%s.root",Analysis[iAnalysis].Data(),Period[flag+1].Data());

 cout << "Flag: " << flag << endl;
 cout << "sig0: " << sig_spectra_filename0 << endl;
 cout << "sig1: " << sig_spectra_filename1 << endl;
 cout << "acc0: " << acc_spectra_filename0 << endl;
 cout << "acc1: " << acc_spectra_filename1 << endl;
 cout << "Input0: " << input_filename0 << endl;
 cout << "Input1: " << input_filename1 << endl;
 cout << "Fn: " << fn_filename << endl;
 cout << "save: " << save_filename << endl;


 fit_shape_2d_P15A(flag,
		   sig_spectra_filename0,
		   sig_spectra_filename1,
		   acc_spectra_filename0,
		   acc_spectra_filename1,
		   input_filename0,
		   input_filename1,
		   fn_filename,
		   save_filename);
 
}
