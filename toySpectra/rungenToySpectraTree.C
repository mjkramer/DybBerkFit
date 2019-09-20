void rungenToySpectraTree(int x){

  gROOT->ProcessLine(".x LoadClasses.C");
  gROOT->ProcessLine(".L genToySpectraTree.C+");
   // gROOT->ProcessLine(".L genToySpectraTree_postP15A.C+");
  
  const Int_t nopts = 18;
  TString options[nopts] = {
    //"allsys",//1
    "sigsys",//2
    "bgsys",//3
    "allsys_and_stat",//5
    "det_eff",//6
    "iav",//7
    "reac_power",//8
    "rel_escale",//9
    "resolutoin",//10
    "scinti_nl",//11
    "distort_aln",//12
    "distort_amc",//13
    "distort_fn",//14
    "distort_li9",//15
    "vary_acc",//16
    "vary_aln",//17
    "vary_amc",//18
    "vary_fn",//19
    "vary_li9"//20
  };
  
  Int_t i = x - 1;
  
  // Set of variation only affected by the reactor flux covariance matrix
  /*
  const Int_t nopts = 2;
  TString options[nopts] = {
    //"allsys_and_stat", //0
    //"allsys",//1
    "sigsys",//2
    //"core_spectra",//5
    "bgsys"//3
    //"nominal"
  };
  */
  TString dataset_name_base = "dyb_data_v1_";
  
  //for (Int_t i = 0; i < nopts; i++){
  //for (Int_t i = 0; i < nopts; i++){ // remove BCW flux variation
    //  for (Int_t i = 22; i < 23; i++){
    TString nominal_dataset_filename = "data_file/dyb_data_v1_nominal.txt";
    TString dataset_filename = "data_file/" + dataset_name_base + options[i] + ".txt";
    TString output_filename = "../outputs/toySpectra_" + options[i] + ".root";
    cout << dataset_filename << "\t" << output_filename << endl;
    genToySpectraTree(nominal_dataset_filename,dataset_filename,output_filename);
    //}
}
