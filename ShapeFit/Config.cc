#include "Config.h"

#include "Binning.h"

using namespace Config;

static void init_nominal_osc_pars() __attribute__((constructor));
static void init_nominal_osc_pars()
{
  if (Binning::useBcwBinning()) {
    S22T13 = S22T13_NOMINAL_BCWBINS;
    DM2EE = DM2EE_NOMINAL_BCWBINS;
  } else {
    // For IHEP binning, just use the LBNL values, since IHEP binning is more
    // similar to LBNL in the lowest two bins.
    S22T13 = S22T13_NOMINAL_LBNLBINS;
    DM2EE = DM2EE_NOMINAL_LBNLBINS;
  }
}
