// Relative paths assume that CWD is a direct child dir of the fitter root,
// e.g. that we are in ShapeFit or toySpectra

#include "Paths.h"

#include "PathUtils.cxx"

namespace Paths {

const char* gReactorSuffix = "_SNF_nonEq";

const char* input(int istage)
{
  return inpath("Theta13-inputs_P17B_inclusive_%s.txt",
                stage_lwc(istage));
}

const char* sig_spectra(int istage)
{
  return inpath("ibd_eprompt_shapes_%s.root",
                stage_lwc(istage));
}

const char* acc_spectra(int istage)
{
  return inpath("accidental_eprompt_shapes_%s.root",
                stage_lwc(istage));
}

const char* fit_result()
{
  return outpath("fit_shape_2d.root");
}

const char* histogram()
{
  return outpath("SuperHistograms.root");
}

const char* predicted_ibd()
{
  return outpath("PredictedIBD.root");
}

const char* response()
{
  return outpath("response/matrix_evis_to_enu_fine.txt");
}

const char* response_rateonly()
{
  return outpath("response/matrix_evis_to_enu_rateonly_fine.txt");
}

const char* response_root()
{
  return outpath("response/evis_to_enu_fine.root");
}

const char* covmatrix(const char* option)
{
  return outpath("covmatrices/matrix_%s.txt", option);
}

const char* sig_covmatrix()
{
  return covmatrix("sigsys");
}

const char* bg_covmatrix()
{
  return covmatrix("bgsys");
}

const char* li9()
{
  return "../li9_spectrum/8he9li_nominal_spectrum.root";
}

const char* amc()
{
  return "../amc_spectrum/amc_spectrum.root";
}

const char* fn(bool ihep)
{
  return ihep ? "../fn_spectrum/P15A_fn_spectrum_IHEP.root"
    : "../fn_spectrum/P15A_fn_spectrum.root";
}

const char* aln()
{
  return "../alpha-n-spectrum/result-DocDB9667.root";
}

const char* muon_decay()
{
  return "../muon_decay_spectrum/MuonDecaySpec.root";
}

const char* baselines()
{
  return "../toySpectra/unblinded_baseline.txt";
}

const char* toyconfig(const char* option)
{
  const char* dir = normalized_or("LBNL_TOY_CONFIG_DIR",
                                  "toySpectra/data_file");
  return LeakStr("%s/dyb_data_v1_%s.txt", dir, option);
}

const char* nominal_toyconfig()
{
  return toyconfig("nominal");
}

const char* nominal_fine_toyconfig()
{
  return toyconfig("nominal_fine");
}

const char* toytree(const char* option)
{
  return outpath("toys/toySpectra_%s.root", option);
}

const char* reactor_spectrum(int istage)
{
  return LeakStr("../ReactorPowerCalculator/isotope_spectra_v4v5v3v1_blinded/reactor_%s%s.txt",
                  stage_upc(istage), gReactorSuffix);
}

const char* reactor_covmatrix()
{
  return "../reactor_covmatrix/p15a/nNu_Mcov_combined_huber-french_u238cor.txt";
}

const char* unified_nl()
{
  const char* tag = checkenv("IBDSEL_USE_SCNL") ? "new" : "old";
  // doc-11611
  // return LeakStr("../toySpectra/unified_nl_data/energymodel_%s_v1.root",
  //                 tag);
  // doc-11646
  // return LeakStr("../toySpectra/unified_nl_data/energymodel_Apr2018_%sE.root",
  //                tag);
  // doc-12585
  return LeakStr("../toySpectra/unified_nl_data/energymodel_Jan2022_%sE_v1.root",
                 tag);
}

// ------------------------------ UNUSED ------------------------------

// NB: "final" is actually old
const char* unified_nl_final()
{
  return "../toySpectra/unified_nl_data/nl_models_final.root";
}

const char* lbnl_positron_data()
{
  return "../toySpectra/lbnl_nl_data/lbnl_positron_nl.txt";
}

const char* abinitio_spectra()
{
  return "../abinitio_spectra/v2-v4/nuSpec_Reactor.txt";
}

const char* bcw_flux()
{
  return "../reactor_covmatrix/bcw/covarMatrix_rawibd.root";
}

const char* bcw_positron_data()
{
  return "../toySpectra/bcw_nl_data/positron.dat";
}

const char* bcw_elec_data()
{
  return "../toySpectra/bcw_nl_data/par.dat";
}

const char* bcw_ele_err()
{
  return "../toySpectra/bcw_nl_data/ele_err.root";
}

} // namespace Paths
