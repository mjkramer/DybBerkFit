#include "genToySpectraTree.H"
#include "Paths.h"

#include <TROOT.h>
#include <TString.h>

#include <iostream>

using namespace std;

void rungenToySpectraTree(const char* option)
{
  TString dataset_filename = Paths::toyconfig(option);
  TString output_filename = Paths::toytree(option);
  cout << dataset_filename << "\t" << output_filename << endl;
  genToySpectraTree(dataset_filename, output_filename);
}
