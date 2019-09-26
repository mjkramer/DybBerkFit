#include <string>
#include <stdexcept>

#include "Binning.h"

namespace Binning {

static const int n_evis_lbnl = 37;
static const int n_evis_bcw = 26;
static const int n_evis_fine_ = 240;

static const int n_enu_ = 156;

static const char* binningUndef =
  "Please set environment variable LBNL_FTR_BINNING to either LBNL or BCW";

static double* evis_lbnl()
{
  double* evis_bins = new double[n_evis_lbnl+1]; // Single bins between 0.7 and 1.0 MeV. 0.2 MeV bins from 1.0 to 8.0 MeV. Single bin between 8.0 and 12 MeV. total 37 bins
  evis_bins[0] = 0.7;
  for (Int_t i = 0; i < n_evis_lbnl-1; i++){
    evis_bins[i+1] = 0.2 *i + 1.0;
  }
  evis_bins[n_evis_lbnl] = 12.0;
    
  return evis_bins;
}

static double* evis_bcw()
{
  double* evis_bins = new double[n_evis_bcw+1]; 
  evis_bins[0] = 0.7;
  for (Int_t i = 0; i < n_evis_bcw-1; i++){
    evis_bins[i+1] = 0.25 *i + 1.3;
  }
  evis_bins[n_evis_bcw] = 12.0;

  return evis_bins;
}

double* evis()
{
  const char* binning = getenv("LBNL_FTR_BINNING");
  if (binning) {
    if (std::string("LBNL") == binning)
      return evis_lbnl();
    if (std::string("BCW") == binning)
      return evis_bcw();
  }

  throw std::runtime_error(binningUndef);
}

int n_evis()
{
  const char* binning = getenv("LBNL_FTR_BINNING");
  if (binning) {
    if (std::string("LBNL") == binning)
      return n_evis_lbnl;
    if (std::string("BCW") == binning)
      return n_evis_bcw;
  }

  throw std::runtime_error(binningUndef);
}

// Previously used by genPredictedIBD
double* evis_fine()
{
  double* evis_bins = new double[n_evis_fine_+1];
  evis_bins[0] = 0.0;
  for (Int_t i = 0; i < n_evis_fine_-1; i++){
    evis_bins[i+1] = 0.05 *i + 0.05;
  }
  evis_bins[n_evis_fine_] = 12.0;
  
  return evis_bins;
}

int n_evis_fine()
{
  return n_evis_fine_;
}

double* enu()
{
  double* enu_bins = new double[n_enu_+1];
  for (Int_t i = 0; i < n_enu_+1; i++){
    enu_bins[i] = 0.05 * i + 1.8;
  }

  return enu_bins;
}

int n_enu()
{
  return n_enu_;
}

}
