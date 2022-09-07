#define LBNL_FIT_STERILE

#include "Config.h"

extern void fit_shape_3d_CLs(bool fit4nuSamples=false, int igrid=-1);

void run_fit_shape_3d_CLs_4n()
{
  fit_shape_3d_CLs(true);
}
