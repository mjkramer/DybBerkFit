#include "readtoyhe8tree_2018model.C" // includes readtoyli9tree_2018model.C

#include <TFile.h>

void gen_2018model()
{
  auto nlfile = new TFile("energymodel_Apr2018_oldE.root");
  readtoyli9tree(nlfile, "toyli9spec_2018model.root");
  readtoyhe8tree(nlfile, "toyhe8spec_2018model.root");
}
