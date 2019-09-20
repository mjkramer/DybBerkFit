void run_build_covmatrix(int x = 1){
  
  gROOT->ProcessLine(".x LoadClasses.C");
  gROOT->ProcessLine(".L build_covmatrix.C+");
    //gROOT->ProcessLine(".L build_covmatrix_postP15A.C+");

  Int_t i = x-1;

  const Int_t nopts_sig = 9;
  TString options_sig[nopts_sig] = {
    "rel_escale", //0
    "scinti_nl",//1
    "iav",//2
    "resolutoin",//3
    "det_eff",//4
    "core_spectra",//5
    "reac_power",//6
    "solar_oscpars",//7
    "sigsys"//8
    // "scinti_nl_corr_positron",//9
    // "scinti_nl_corr_gamma",//10
  };


  const Int_t nopts_bkg = 9;
  TString options_bkg[nopts_bkg] = {
    //"distort_aln", // 0
    "distort_amc",
    "distort_fn",
    "distort_li9",
    "vary_acc", // 4
    "vary_aln",
    "vary_amc",
    "vary_fn",
    "vary_li9",
    "bgsys" // 9
  };
  
  //for (Int_t i = 0; i < nopts_sig; i++){
  if (i < nopts_sig){
    //TString toymc_filename = "../outputs/toySpectra_" + options_sig[i] + "_2017Model_p17b.root";
    //TString covmatrix_filename = "covariance_matrices/matrix_" + options_sig[i] + "_2017Model_P17B.txt";
      
      TString toymc_filename = "../outputs/toySpectra_" + options_sig[i] + "_2017Model_p17b_0.2_inflated.root";
      TString covmatrix_filename = "covariance_matrices/matrix_" + options_sig[i] + "_2017Model_P17B_0.2_inflated.txt";
      
    cout << toymc_filename << "\t" << covmatrix_filename << endl;
    build_covmatrix(toymc_filename.Data(),covmatrix_filename.Data(),0);
  }else{
  
    //for (Int_t i = 0; i < nopts_bkg; i++){
    Int_t iii = i - nopts_sig;
    //TString toymc_filename = "../outputs/toySpectra_" + options_bkg[iii] + "_2017Model_p17b.root";
    //TString covmatrix_filename = "covariance_matrices/matrix_" + options_bkg[iii] + "_2017Model_P17B.txt";
      
      TString toymc_filename = "../outputs/toySpectra_" + options_bkg[iii] + "_2017Model_p17b_0.2_inflated.root";
      TString covmatrix_filename = "covariance_matrices/matrix_" + options_bkg[iii] + "_2017Model_P17B_0.2_inflated.txt";
      
      
    //TString toymc_filename = "../outputs/toySpectra/nominal/toySpectra_" + options_bkg[iii] + ".root";
    // TString covmatrix_filename = "covariance_matrices/matrix_" + options_bkg[iii] + ".txt";
    cout << toymc_filename << "\t" << covmatrix_filename << endl;
    build_covmatrix(toymc_filename.Data(),covmatrix_filename.Data(),1);
  }
  
  
}
