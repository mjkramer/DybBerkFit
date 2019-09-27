void rungenPredictedIBD(int x){

  gROOT->ProcessLine(".x LoadClasses.C");
  gROOT->ProcessLine(".L genPredictedIBD.C+");
  
  Int_t i = x - 1;
  
  Double_t s22t13 = 0;
  Double_t s22t13_max = 0.10;
  Double_t s22t13_min = 0.00;
  Int_t nsteps = 100;

  Double_t s22t13_bin = (s22t13_max - s22t13_min)/nsteps;

  s22t13 = s22t13_min + (i * s22t13_bin); 

  Char_t output_filename[256];
  sprintf(output_filename,"../ShapeFit/PredictedIBD/sin22t13/PredictedIBD_%i.root",i);
   
    cout << i << "\t" << s22t13 << "\t" << output_filename << endl;
    genPredictedIBD(s22t13, -1, output_filename);
    //}
}
