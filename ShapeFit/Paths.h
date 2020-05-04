#pragma once

// These functions return either static strings or volatile strings living
// inside the circular buffer of ROOT's ::Form.

namespace Paths {

const char* fit_result();

const char* input(int istage);
const char* sig_spectra(int istage);
const char* acc_spectra(int istage);

const char* histogram();
const char* predicted_ibd();

const char* response();
const char* response_rateonly();
const char* response_root();

const char* covmatrix(const char* option);
const char* sig_covmatrix();
const char* bg_covmatrix();

const char* li9();
const char* amc();
const char* fn(bool ihep=false);
const char* aln();

const char* baselines();

const char* toyconfig(const char* option);
const char* nominal_toyconfig();
const char* nominal_fine_toyconfig();

const char* toytree(const char* option);

const char* reactor_spectrum(int istage);
const char* reactor_covmatrix();

const char* unified_nl(bool scnl=false);

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
