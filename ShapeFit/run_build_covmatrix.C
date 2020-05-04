#include "Paths.h"

#include <TString.h>

#include <iostream>

using namespace std;

void build_covmatrix(const Char_t* toymc_filename,
                     const Char_t* output_filename, Int_t bkg_flag);

void run_build_covmatrix(int x = 1)
{
  Int_t i = x - 1;

  const Int_t nopts_sig = 11;
  TString options_sig[nopts_sig] = {
      "rel_escale",              // 0
      "scinti_nl",               // 1
      "iav",                     // 2
      "resolutoin",              // 3
      "det_eff",                 // 4
      "core_spectra",            // 5
      "reac_power",              // 6
      "solar_oscpars",           // 7
      "sigsys",                  // 8
      "scinti_nl_corr_positron", // 9
      "scinti_nl_corr_gamma",    // 10
  };


  const Int_t nopts_bkg = 10;
  TString options_bkg[nopts_bkg] = {
      "distort_aln", // 0
      "distort_amc", // 1
      "distort_fn",  // 2
      "distort_li9", // 3
      "vary_acc",    // 4
      "vary_aln",    // 5
      "vary_amc",    // 6
      "vary_fn",     // 7
      "vary_li9",    // 8
      "bgsys",       // 9
  };

  if (i < nopts_sig) {
    auto toymc_filename = Paths::toytree(options_sig[i]);
    auto covmatrix_filename = Paths::covmatrix(options_sig[i]);

    cout << toymc_filename << "\t" << covmatrix_filename << endl;
    build_covmatrix(toymc_filename, covmatrix_filename, 0);
  } else {
    Int_t iii = i - nopts_sig;

    auto toymc_filename = Paths::toytree(options_bkg[iii]);
    auto covmatrix_filename = Paths::covmatrix(options_bkg[iii]);

    cout << toymc_filename << "\t" << covmatrix_filename << endl;
    build_covmatrix(toymc_filename, covmatrix_filename, 1);
  }

  cout << "Bye bye!" << endl;
}
