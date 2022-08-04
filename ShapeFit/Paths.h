#pragma once

#include <TString.h>

#include <vector>

// These functions are allowed to leak memory, so don't call them too
// repeatedly.

namespace Paths {

extern const char* gReactorSuffix;

const char* outpath(const char* fmt, ...);
const char* inpath(const char* fmt, ...);

const char* fit_result();

const char* input(int istage);
const char* sig_spectra(int istage);
std::vector<TString> all_sig_spectra();
const char* acc_spectra(int istage);

const char* histogram();
const char* predicted_ibd();

const char* response();
const char* response_rateonly();
const char* response_root();

const char* covmatrix(const char* option);
const char* sig_covmatrix();
const char* bg_covmatrix();
const char* dm2ee_covmatrix();

const char* li9();
const char* amc();
// Counterintuitively, ihep=true corresponds to Bei-Zhen's numbers in doc-10948
// (make_P15A_spectrum_IHEP.C), which are for an LBNL-like multiplicity cut.
// Meanwhile, ihep=false corresponds to Xiangpan's numbers
// (make_P15A_spectrum.C), which are for an IHEP-like multiplicity cut. So we
// actually want ihep=true. If we are using an IHEP-like selection, we want
// ihep=false. Go figure.
const char* fn(bool ihep=true);
const char* aln();
const char* muon_decay();

const char* baselines();

const char* toyconfig(const char* option);
const char* nominal_toyconfig();
const char* nominal_fine_toyconfig();

const char* toytree(const char* option);

const char* reactor_spectrum(int istage);
const char* reactor_covmatrix();

const char* unified_nl();

// ------------------------------ UNUSED ------------------------------

// NB: "final" is actually old
const char* unified_nl_final();

const char* lbnl_positron_data();

const char* abinitio_spectra();

const char* bcw_flux();

const char* bcw_positron_data();
const char* bcw_elec_data();
const char* bcw_ele_err();

} // namespace Paths
