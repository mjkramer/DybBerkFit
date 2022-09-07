#include <cstdlib>

using namespace std;

void genToySpectraTree_parscans_CLs(const char *dataset_filename, int itask = 0,
                                    int ntasks = 1);

void rungenToySpectraTree_parscans_CLs() {
  genToySpectraTree_parscans_CLs("nominal", atoi(getenv("SLURM_LOCALID")),
                                 atoi(getenv("SLURM_NTASKS")));
}
