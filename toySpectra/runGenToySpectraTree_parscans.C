#define LBNL_FIT_STERILE

#include "Config.h"
#include "Paths.h"
#include "genToySpectraTree.H"

#include <TROOT.h>

using namespace Config;

void rungenToySpectraTree_parscan(int igrid = 0)
{
  const int npoints = nsteps * nsteps_dm2 * nsteps_s22t14 * nsteps_dm214;

  Double_t s2t13_array[npoints+1];
  Double_t dm2ee_array[npoints+1];
  Double_t s2t14_array[npoints+1];
  Double_t dm241_array[npoints+1];
  Int_t dm241_id[npoints+1];
  Int_t ipoint = 0;
  for (Int_t it14 = 0; it14 < npoints_s2t14; it14++){
    for (Int_t im41 = 0; im41 < npoints_dm241; im41++){
      for (Int_t it13 = 0; it13 < npoints_s2t13; it13++){
        for (Int_t imee = 0; imee < npoints_dm2ee; imee++){
          s2t13_array[ipoint] = s2t13 + 0.02*it13;
          dm2ee_array[ipoint] = dm2ee + 0.0001*imee;
          s2t14_array[ipoint] = exp(log(s2t14_start) + log_s2t14_step*it14);
          if (im41 < npoints_dm241_1){
            dm241_array[ipoint] = exp(log(dm214_start) + log_dm214_step*im41);
          }
          else{
            dm241_array[ipoint] = dm214_start_2 +lin_dm214_step*(im41-npoints_dm241_1);
          }
          dm241_id[ipoint] = im41+1;
          ipoint++;
        }
      }
    }
  }

  s2t13_array[npoints] = s2t13;
  dm2ee_array[npoints] = dm2ee;
  s2t14_array[npoints] = 0.0;
  dm241_array[npoints] = 0.0;
  dm241_id[npoints] = 0.0;

  Int_t ipar_start = 0;
  Int_t ipar_end = npoints;
  if (igrid > 0 && igrid <= npoints){
    ipar_start = igrid-1;
    ipar_end = igrid;
  }else if (igrid == 0){
    ipar_start = npoints;
    ipar_end = npoints+1;
  }

  for (int ipar = ipar_start; ipar < ipar_end; ++ipar) {
    auto outname = Paths::outpath("toys_parscan/toySpectra_allsys_w_dm2ee_and_stat_s2t14_%4.4f_dm241_%5.5f.root",
                                  s22t14, dm214);
    genToySpectraTree(Paths::toyconfig("nominal"), outname, s2t13_array[ipar]; dm2ee_array[ipar],
                      s2t14_array[ipar], dm241_array[ipar], 0)
  }
}
