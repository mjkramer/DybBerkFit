#pragma once

namespace Binning {

extern bool useBcwBinning();
extern bool useIhepBinning();
extern bool useLbnlBinning();
extern bool useFineBinning();

extern double min_energy();

extern double* evis();
extern double* evis_fine();
extern double* enu();

extern int n_evis();
extern int n_evis_fine();
extern int n_enu();

} // namespace Binning
