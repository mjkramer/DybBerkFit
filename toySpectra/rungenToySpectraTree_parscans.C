#include "Paths.h"

#include <TString.h>

void genToySpectraTree_parscans(TString dataset_filename, TString output_base, int igrid = -1);

void rungenToySpectraTree_parscans()
{
  genToySpectraTree_parscans(Paths::toyconfig("nominal"),
                             Paths::outpath("toys_parscan/toySpectra_allsys_w_dm2ee_and_stat"));
}
