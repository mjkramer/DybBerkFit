#include "Binning.h"
#include "Config.h"

#include <cstdlib>              // getenv
#include <cstring>              // strcmp

namespace Binning {

static const int _n_evis_lbnl = 37;
static const int _n_evis_bcw = 26;
static const int _n_evis_fine = 240;
static const int _n_enu = 156;

static double _evis_lbnl[_n_evis_lbnl + 1];
static double _evis_bcw[_n_evis_bcw + 1];
static double _evis_fine[_n_evis_fine + 1];
static double _enu[_n_enu + 1];

static bool useBcwBinning()
{
  const char* val = getenv("LBNL_FIT_BINNING");
  return val && strcmp("BCW", val) == 0;
}

static double min_energy()
{
  const char* val = getenv("LBNL_FIT_EMIN");
  return val ? atof(val) : 0.7;
}

static void init_evis_lbnl() __attribute__((constructor));
static void init_evis_lbnl()
{
  // Single bins between 0.7 and 1.0 MeV. 0.2 MeV bins from 1.0 to 8.0 MeV.
  // Single bin between 8.0 and 12 MeV. total 37 bins
  _evis_lbnl[0] = min_energy();
  for (Int_t i = 0; i < _n_evis_lbnl - 1; i++) {
    _evis_lbnl[i + 1] = 0.2 * i + 1.0;
  }
  _evis_lbnl[_n_evis_lbnl] = 12.0;
}

static void init_evis_bcw() __attribute__((constructor));
static void init_evis_bcw()
{
  _evis_bcw[0] = min_energy();
  for (Int_t i = 0; i < _n_evis_bcw - 1; i++) {
    _evis_bcw[i + 1] = 0.25 * i + 1.3;
  }
  _evis_bcw[_n_evis_bcw] = 12.0;
}

// Previously used by genPredictedIBD
static void init_evis_fine() __attribute__((constructor));
static void init_evis_fine()
{
  _evis_fine[0] = 0.0;
  for (Int_t i = 0; i < _n_evis_fine - 1; i++) {
    _evis_fine[i + 1] = 0.05 * i + 0.05;
  }
  _evis_fine[_n_evis_fine] = 12.0;
}

static void init_enu() __attribute__((constructor));
static void init_enu()
{
  for (Int_t i = 0; i < _n_enu + 1; i++) {
    _enu[i] = 0.05 * i + 1.8;
  }
}

double* evis()
{
  return useBcwBinning() ? _evis_bcw : _evis_lbnl;
}

double* evis_fine()
{
  return _evis_fine;
}

double* enu()
{
  return _enu;
}

int n_evis()
{
  return useBcwBinning() ? _n_evis_bcw : _n_evis_lbnl;
}

int n_evis_fine()
{
  return _n_evis_fine;
}

int n_enu()
{
  return _n_enu;
}

} // namespace Binning
