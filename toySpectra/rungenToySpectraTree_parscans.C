#include <cstdlib>

using namespace std;

void genToySpectraTree_parscans(const char *dataset_filename, int itask = 0,
                                int ntasks = 1);

void rungenToySpectraTree_parscans() {
  genToySpectraTree_parscans("nominal", atoi(getenv("SLURM_LOCALID")),
                             atoi(getenv("SLURM_NTASKS")));
}
