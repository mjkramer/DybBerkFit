#pragma once

#include <TString.h>

void genToySpectraTree(TString dataset_filename, TString output_filename,
                       double s2t13 = -1, double dm2ee = -1, double s2t14 = -1,
                       double dm241 = -1, int nToys = 1000);
