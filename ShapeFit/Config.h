namespace Config {

const char* input_filename0 = "./Inputs/Theta13-inputs_P17B_inclusive_6ad.txt";
const char* input_filename1 = "./Inputs/Theta13-inputs_P17B_inclusive_8ad.txt";
const char* input_filename2 = "./Inputs/Theta13-inputs_P17B_inclusive_7ad.txt";

const char* sig_spectra_filename0 = "./Spectra/ibd_eprompt_shapes_6ad.root";
const char* sig_spectra_filename1 = "./Spectra/ibd_eprompt_shapes_8ad.root";
const char* sig_spectra_filename2 = "./Spectra/ibd_eprompt_shapes_7ad.root";

const char* acc_spectra_filename0 = "./Spectra/accidental_eprompt_shapes_6ad.root";
const char* acc_spectra_filename1 = "./Spectra/accidental_eprompt_shapes_8ad.root";
const char* acc_spectra_filename2 = "./Spectra/accidental_eprompt_shapes_7ad.root";

const char* histogram_filename = "./Flux/SuperHistograms_P17B_2017Model_fine_huber-french.root";
const char* predicted_ibd_filename = "./PredictedIBD/PredictedIBD_P17B_2017Model_fine_huber-french.root";

const char* response_filename = "./matrix_evis_to_enu_fine_2017Model_P17B.txt";
const char* response_root_filename =  "../outputs/evis_to_enu_fine_2017Model_p17b.root";
const char* response_filename_rateonly = "matrix_evis_to_enu_rateonly_fine_2017_Model_P17B.txt";

const char* sig_matrix_filename = "covariance_matrices/matrix_sigsys.txt";
const char* bg_matrix_filename = "covariance_matrices/matrix_bgsys.txt";

const char* li9_filename = "../li9_spectrum/8he9li_nominal_spectrum.root";
const char* amc_filename = "../amc_spectrum/amc_spectrum.root";
// const char* fn_filename = "../fn_spectrum/P15A_fn_spectrum.root";
const char* fn_filename = "../fn_spectrum/P15A_fn_spectrum_IHEP.root"; // XXX
const char* aln_filename = "../alpha-n-spectrum/result-DocDB9667.root";

const char* baselines_filename = "Distances/unblinded_baseline.txt";

// relative to toySpectra
const char* nominal_dataset_filename = "../toySpectra/data_file/dyb_data_v1_nominal.txt";
const char* nominal_noosc_dataset_filename = "../toySpectra/data_file/dyb_data_v1_nominal_noosc.txt";

const double lowBinInflation = 0; // was 0.1 originally

}
