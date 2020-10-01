#include "readtoyhe8tree.C"     // includes readtoyli9tree.C

#include <TFile.h>

// The output should be identical to toyXXspec_IHEPmodel_v1.root
void gen_IHEPmodel_v1()
{
  auto nlfile = new TFile("Model1.root");
  readtoyli9tree(nlfile, "toyli9spec_IHEPmodel_v1.root");
  readtoyhe8tree(nlfile, "toyhe8spec_IHEPmodel_v1.root");
}
