#include <fstream>

#include "Config.h"

using namespace Config;

void dump_parscans_FC(const char* outpath)
{
    // override values from Config.h
  const int nsteps = 1;
  const double s22t13start = S22T13;
  const double s22t13end = S22T13;

  const int nsteps_dm2 = 1;
  const double dm2eestart = DM2EE;
  const double dm2eeend = DM2EE;

  // Uncomment the commented lines in the block below if nsteps(_dm2) > 1
  // And comment the following assert:
  assert(nsteps == 1 && nsteps_dm2 == 1);

  const int npoints = nsteps * nsteps_dm2 * nsteps_s22t14 * nsteps_dm214_all;

  Double_t s2t13_array[npoints + 1];
  Double_t dm2ee_array[npoints + 1];
  Double_t s2t14_array[npoints + 1];
  Double_t dm214_array[npoints + 1];
  Int_t ipoint = 0;
  for (Int_t it14 = 0; it14 < nsteps_s22t14; it14++) {
    for (Int_t im41 = 0; im41 < nsteps_dm214_all; im41++) {
      for (Int_t it13 = 0; it13 < nsteps; it13++) {
        for (Int_t imee = 0; imee < nsteps_dm2; imee++) {
          s2t13_array[ipoint] = s22t13start;
          dm2ee_array[ipoint] = dm2eestart;

          // Uncomment lines lines if nsteps(_dm2) > 1
          // They cause a warning when nsteps is 1, so we commented them

          // if (nsteps > 1)
          //   s2t13_array[ipoint] =
          //       s22t13start + it13 * (s22t13end - s22t13start) / (nsteps - 1);
          // if (nsteps_dm2 > 1)
          //   dm2ee_array[ipoint] =
          //       dm2eestart + imee * (dm2eeend - dm2eestart) / (nsteps_dm2 - 1);

          s2t14_array[ipoint] = exp(log(s22t14start) + log_s22t14_step * it14);

          if (im41 < nsteps_dm214) {
            dm214_array[ipoint] = exp(log(dm214start) + log_dm214_step * im41);
          } else {
            dm214_array[ipoint] =
                dm214start_2 + lin_dm214_step * (im41 - nsteps_dm214);
          }

          ipoint++;
        }
      }
    }
  }

  s2t13_array[npoints] = S22T13;
  dm2ee_array[npoints] = DM2EE;
  s2t14_array[npoints] = 0.0;
  dm214_array[npoints] = 0.0;

  std::ofstream outfile(outpath);

  for (int ipar = 0; ipar < npoints+1; ++ipar) {
    // outfile << s2t13_array[ipar] << " " << dm2ee_array[ipar] << " "
    //         << s2t14_array[ipar] << " " << dm214_array[ipar] << std::endl;
    outfile << Form("%4.4f %5.5f %4.4f %5.5f",
                    s2t13_array[ipar], dm2ee_array[ipar],
                    s2t14_array[ipar], dm214_array[ipar]) << std::endl;
  }
}
