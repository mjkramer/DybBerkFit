#include "genToySpectraTree.H"
#include "Paths.h"

#include <TROOT.h>
#include <TString.h>

#include <iostream>

using namespace std;

void rungenToySpectraTree(int x)
{
  const Int_t nopts = 19;
  TString options[nopts] = {
      "allsys",          // 1
      "sigsys",          // 2
      "bgsys",           // 3
      "allsys_and_stat", // 5
      "det_eff",         // 6
      "iav",             // 7
      "reac_power",      // 8
      "rel_escale",      // 9
      "resolutoin",      // 10
      "scinti_nl",       // 11
      "distort_aln",     // 12
      "distort_amc",     // 13
      "distort_fn",      // 14
      "distort_li9",     // 15
      "vary_acc",        // 16
      "vary_aln",        // 17
      "vary_amc",        // 18
      "vary_fn",         // 19
      "vary_li9"         // 20
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

  TString dataset_filename = Paths::toyconfig(options[i]);
  TString output_filename = Paths::toytree(options[i]);
  cout << dataset_filename << "\t" << output_filename << endl;
  genToySpectraTree(dataset_filename, output_filename);
  //}
}
