#include "Paths.h"

#include <TString.h>

#include <iostream>
#include <cstdlib>

using namespace std;

void build_covmatrix(const Char_t* toymc_filename,
                     const Char_t* output_filename, Int_t bkg_flag);

void run_build_covmatrix(const char* option, int bkg_flag)
{
  auto toymc_filename = Paths::toytree(option);
  auto covmatrix_filename = Paths::covmatrix(option);

  cout << toymc_filename << "\t" << covmatrix_filename << endl;
  build_covmatrix(toymc_filename, covmatrix_filename, bkg_flag);

  cout << "Bye bye!" << endl;
  quick_exit(0);
}
