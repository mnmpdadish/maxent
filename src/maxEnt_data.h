#pragma once


//#include "includeDef.h"
#include <string>
#include <list>
#include <vector>
#include <map>
#include <stdio.h>
#include <cmath>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <ctime>
#include "armadillo"
#include "utilities.h"
#include <cstring>
#include <cmath>

#include "gnuplot_pipe.h"


using namespace arma;

//MAIN INPUT FILE PARAMETERS

static string data_file_param("data file:");

enum Data_params_name {BOSON, TAU_GF, TEM, NORM_A, M_1, ERR_M1, M_2, ERR_M2, M_3, ERR_M3, CONST_PART_G};

static map<Data_params_name, string> Data_params( {	{BOSON,"bosonic data (yes/[no]):"},
	{TAU_GF,"imaginary time data (yes/[no]):"},
	{TEM,"temperature (in energy units, k_B=1):"},
	{NORM_A,"norm of spectral function:"},
	{M_1,"1st moment:"},
	{ERR_M1,"1st moment error:"},
	{M_2,"2nd moment:"},
	{ERR_M2,"2nd moment error:"},
	{M_3,"3rd moment:"},
	{ERR_M3,"3rd moment error:"},
	{CONST_PART_G,"static part of data:"}} );

enum Intput_files_params_name {INPUT_DIR, COL_GR, COL_GI, ERROR_FILE, COL_ERR_GR, COL_ERR_GI, COVAR_RE_RE_FILE, COVAR_IM_IM_FILE, COVAR_RE_IM_FILE, COL_G_TAU, COL_ERR_G_TAU, COVAR_TAU_FILE};

static map<Intput_files_params_name, string> Input_files_params( { {INPUT_DIR, "input directory:"},
	{COL_GR, "Re(G) column in data file (default: 2):"},
	{COL_GI, "Im(G) column in data file (default: 3):"},
	{ERROR_FILE, "error file:"},
	{COL_ERR_GR, "Re(G) column in error file (default: 2):"},
	{COL_ERR_GI, "Im(G) column in error file (default: 3):"},
	{COVAR_RE_RE_FILE, "re-re covariance file:"},
	{COVAR_IM_IM_FILE, "im-im covariance file:"},
	{COVAR_RE_IM_FILE, "re-im covariance file:"},
	{COL_G_TAU, "column of G(tau) in data file (default: 2):"},
	{COL_ERR_G_TAU, "column of G(tau) error in error file (default: 2):"},
	{COVAR_TAU_FILE, "imaginary time covariance file:"} } );


enum Grig_params_name {CUTOFF_WN, SPECTR_FUNC_WIDTH, SPECTR_FUNC_CENTER, W_ORIGIN, STEP_W, GRID_W_FILE, NON_UNIFORM_GRID, USE_GRID_PARAMS, PARAM_GRID_PARAMS};

static map<Grig_params_name, string> Grid_params( { {CUTOFF_WN, "Matsubara frequency cutoff (in energy units, k_B=1):"},
	{SPECTR_FUNC_WIDTH, "spectral function width:"},
	{SPECTR_FUNC_CENTER, "spectral function center:"},
	{W_ORIGIN, "real frequency grid origin:"},
	{STEP_W, "real frequency step:"},
	{GRID_W_FILE, "real frequency grid file:"},
	{NON_UNIFORM_GRID,"use non uniform grid in main spectral range (yes/[no]):"},
	{USE_GRID_PARAMS, "use parameterized real frequency grid (yes/[no]):"},
	{PARAM_GRID_PARAMS, "grid parameters (w_0 dw_0 w_1 dw_1 ... w_{N-1} dw_{N-1} w_N):"} } );


enum Preproc_comp_params_name {EVAL_MOMENTS, MAX_M, DEFAULT_MODEL_CENTER, DEFAULT_MODEL_WIDTH, DEFAULT_MODEL_SHAPE, DEFAULT_MODEL_FILE, INIT_SPECTR_FUNC_FILE};
//, INTERP_TYPE
//PEAKED_DEFAULT_MODEL,
static map<Preproc_comp_params_name, string> Preproc_comp_params( { {EVAL_MOMENTS, "evaluate moments (yes/[no]):"},
	{MAX_M, "maximum moment:"},
	{DEFAULT_MODEL_CENTER, "default model center (default: 1st moment):"},
	{DEFAULT_MODEL_WIDTH, "default model half width (default: standard deviation):"},
	{DEFAULT_MODEL_SHAPE, "default model shape parameter (default: 2):"},
	{DEFAULT_MODEL_FILE, "default model file:"},
	{INIT_SPECTR_FUNC_FILE, "initial spectral function file:"} } );
//{INTERP_TYPE, "interpolation type (spline (default), quad, lin):"}
//	{PEAKED_DEFAULT_MODEL, "use peaked default model (yes/no):"},

enum Preproc_exec_params_name {PREPROSSESS_ONLY, DISPL_PREP_FIGS, DISPL_ADV_PREP_FIGS, PRINT_OTHER_PARAMS};

static map<Preproc_exec_params_name, string> Preproc_exec_params( { {PREPROSSESS_ONLY, "preprocess only (yes/[no]):"},
	{DISPL_PREP_FIGS, "display preprocessing figures (yes/[no]):"},
	{DISPL_ADV_PREP_FIGS, "display advanced preprocessing figures (yes/[no]):"},
	{PRINT_OTHER_PARAMS, "print other parameters (yes/[no]):"} } );


enum Output_files_params_name {OUTPUT_DIR, OUTPUT_NAME_SUFFIX, ALPHA_SAVE_MAX, ALPHA_SAVE_MIN, W_SAMPLE};
//, SAVE_A_OPT_ONLY

static map<Output_files_params_name, string> Output_files_params( { {OUTPUT_DIR, "output directory:"},
	{OUTPUT_NAME_SUFFIX, "output file names suffix:"},
	{ALPHA_SAVE_MAX, "maximum alpha for which results are saved:"},
	{ALPHA_SAVE_MIN, "minimum alpha for which results are saved:"},
	{W_SAMPLE, "spectral function sample frequencies (w_1 w_2 ... w_N):"} } );
//{ALPHA_SAVE, "alpha value below which all spectral functions are saved:"},
//{SAVE_A_OPT_ONLY, "save only optimal spectrum (yes/[no]):"}


enum Optim_comp_params_name {ALPHA_INIT, ALPHA_MIN, ALPHA_OPT_MAX, ALPHA_OPT_MIN};
//CHI2_ALPHA_SMOOTH_RANGE

static map<Optim_comp_params_name, string> Optim_comp_params( { {ALPHA_INIT, "initial value of alpha:"},
	{ALPHA_MIN, "minimum value of alpha:"},
	{ALPHA_OPT_MAX, "maximum optimal alpha:"},
	{ALPHA_OPT_MIN, "minimum optimal alpha:"} } );
//{CHI2_ALPHA_SMOOTH_RANGE, "chi2 vs alpha smoothing range in log10 scale:"}

enum Optim_exec_params_name {N_ALPHA, INITIALIZE_MAXENT, INITIALIZE_PREPROC, INTERACTIVE_MODE};

static map<Optim_exec_params_name, string> Optim_exec_params( { {N_ALPHA, "number of values of alpha computed in one execution:"},
	{INITIALIZE_MAXENT, "initialize maxent (yes/[no]):"},
	{INITIALIZE_PREPROC, "initialize preprocessing (yes/[no]):"},
	{INTERACTIVE_MODE, "interactive mode ([yes]/no):"} } );


enum Optim_displ_params_name {PRINT_ALPHA, SHOW_OPTIMAL_ALPHA_FIGS, SHOW_LOWEST_ALPHA_FIGS, SHOW_ALPHA_CURVES, REF_SPECTR_FILE};
//DISPL_OPTIM_FIGS,

static map<Optim_displ_params_name, string> Optim_displ_params( { {PRINT_ALPHA, "print results at each value of alpha (yes/[no]):"},
	{SHOW_OPTIMAL_ALPHA_FIGS,"show optimal alpha figures ([yes]/no):"},
	{SHOW_LOWEST_ALPHA_FIGS,"show lowest alpha figures ([yes]/no):"},
	{SHOW_ALPHA_CURVES,"show alpha dependant curves ([yes]/no):"},
	{REF_SPECTR_FILE, "reference spectral function file:"} } );
//{DISPL_OPTIM_FIGS, "display minimization figures (yes/[no]):"},

// OTHER PARAMETERS

extern "C++"
{
	class OmegaMaxEnt_data: public utilities
    {
    public:
        OmegaMaxEnt_data(int arg_N, char *args[]);
        ~OmegaMaxEnt_data();
        
        void loop_run();
        
    private:
        string fileName; //mc
        
        
        
		//bool load_input_params();
		//bool load_other_params();
		bool load_data_file(mat &data_array, string file_name);
		bool preproc();
		bool set_G_omega_n_fermions();
		bool set_G_omega_n_bosons();
		bool set_covar_G_omega_n();
		bool set_covar_chi_omega_n();
		bool set_covar_Gtau();
		bool set_moments_fermions();
		bool set_moments_bosons();
		bool compute_moments_omega_n();
		bool compute_moments_chi_omega_n();
		bool compute_moments_tau_fermions();
		bool compute_moments_tau_bosons();
		bool test_low_energy_peak_fermions();
		bool test_low_energy_peak_bosons();
		bool test_low_energy_peak_chi();
		bool set_initial_spectrum();
		bool set_initial_spectrum_chi();
		bool set_grid_omega_from_file();
		bool set_grid_omega_from_file_chi();
		bool set_grid_from_params();
		bool set_grid_from_params_chi();
		bool set_wc();
		bool set_wc_chi();
		bool set_omega_grid();
		bool set_omega_grid_chi();
		bool set_smooth_step_grid();
		bool set_quad_step_grid();
		bool set_exp_step_grid();
		bool set_A_ref();
		bool truncate_G_omega_n();
		bool truncate_chi_omega_n();
		bool set_default_model();
		bool set_default_model_chi();
		bool default_model_val_G(vec x, vec x0, vec coeffs, vec gaussians_params, vec &dm);
		bool default_model_val_chi(vec x, vec x0, vec coeffs, vec gaussians_params, vec &dm);
		bool Kernel_G_fermions();
		bool Kernel_G_fermions_grid_transf();
		bool Kernel_G_fermions_grid_transf_1();
		bool Kernel_G_fermions_grid_transf_omega();
		bool Kernel_G_bosons();
		bool Kernel_chi();
		bool diagonalize_covariance();
		bool diagonalize_covariance_chi();
		bool Fourier_transform_G_tau();
		
		bool non_uniform_frequency_grid(rowvec w_steps, rowvec wlims, double w0, vec R, vec &grid);
		void spline_coeffs_rel(double *x0, double *V, int N0, double *coeffs);
		bool spline_val(vec x, vec x0, vec coeffs, vec &s);
		bool sum_gaussians(vec x, rowvec x0, rowvec s0, rowvec wgt, vec &F);
		bool sum_gaussians_chi(vec x, rowvec x0, rowvec s0, rowvec wgt, vec &F);
		bool general_normal(vec x, double x0, double s0, double p, vec &F);
		double general_normal_val(double x, void *par[]);
		bool spline_G_part(vec x, uvec ind_xlims, vec xs, vec F, vec &coeffs);
		bool spline_G_omega_u(vec x, uvec ind_xlims, vec xs, vec F, vec &coeffs);
		bool spline_matrix_G_part(vec x, uvec ind_xlims, vec xs, mat &M);
		bool spline_matrix_grid_transf(vec w0, mat &M);
		bool spline_matrix_grid_transf_G_part(vec x, uvec ind_xlims, vec xs, mat &M);
		bool spline_matrix_grid_transf_G_part_1(vec x, uvec ind_xlims, vec xs, mat &M);
		bool spline_val_grid_transf(vec x, vec x0, vec coeffs, vec &s);
		bool spline_chi_part(vec x, uvec ind_xlims, vec xs, vec F, vec &coeffs);
		bool spline_matrix_chi(vec x, uvec ind_xlims, vec xs, mat &M);
		bool spline_val_G_part(vec x, vec x0, uvec ind_xlims, vec xs, vec coeffs, vec &sv);
		bool spline_val_G_part_grid_transf(vec x, vec x0, uvec ind_xlims, vec xs, vec coeffs, vec &sv);
		bool spline_val_G_part_grid_transf_1(vec x, vec x0, uvec ind_xlims, vec xs, vec coeffs, vec &sv);
		bool spline_val_chi_part(vec x, vec x0, uvec ind_xlims, vec xs, vec coeffs, vec &sv);
		double spline_val_G_part(double x, vec &x0, uvec &ind_xlims, vec &xs, vec &coeffs);
		bool fit_circle_arc(vec x, vec y, vec &arc_params);
		
		void minimize();
		//void minimize_increase_alpha();
		
		double spline_val_G_part_int(double x, void *par[]);
		
		//void plot(graph_2D &g, vec x, vec y, char *xl, char *yl, char *attr=NULL);
        
        bool create_default_input_params_file();
        bool create_default_other_params_file();
		
        void init_variables();
        void free_variables();
        
        string input_params_file_name;
        string other_params_file_name;
        
        //! internal computation parameters
        int Nn_max, Nw_min, Nw_max, Nn_fit_max, Nn_fit_fin, Niter_dA_max, Nalpha_max_figs, Nwsamp;
        double f_w_range, f_SW_std_omega, f_width_grid_dens, tol_tem, tol_norm, tol_R_G0_Gbeta, tol_M1, tol_M2, tol_M3, default_error_G, err_norm, default_error_M, tol_mean_C1, tol_std_C1, tol_rdw, Rmin_Dw_dw, Rdw_max, RW_grid, RWD_grid, minDefM, f_alpha_init, R_width_ASmin, f_Smin, diff_chi2_max, tol_int_dA, rc2H, fc2, pow_alpha_step_init, pow_alpha_step_min, chi2_alpha_smooth_range, f_scale_lalpha_lchi2, FNfitTauW, std_norm_peak_max, varM2_peak_max, peak_weight_min, RMAX_dlchi2_lalpha, f_alpha_min, save_alpha_range,  Rmin_SW_dw, R_peak_width_dw, R_wncutoff_wr, R_Dw_dw, R_SW_wr, R_wmax_wr_min;
        //tol_quad, f_chi2save, f_chi2min, Nwn_test_metal, R_d2G_chi_peak,
        
        //! input parameters
        string input_dir_in, *input_dir_prec, input_dir, data_file_name_in, *data_file_name_prec, data_file_name, boson_in, *boson_prec, tau_GF_in, *tau_GF_prec, tem_in, *tem_prec, M0_in, *M0_prec, M1_in, *M1_prec, errM1_in, *errM1_prec, M2_in, *M2_prec, errM2_in, *errM2_prec, M3_in, *M3_prec, errM3_in, *errM3_prec, static_part_G_in, *static_part_G_prec, col_Gr_in, *col_Gr_prec, col_Gi_in, *col_Gi_prec, error_file_in, *error_file_prec, error_file, col_errGr_in, *col_errGr_prec, col_errGi_in, *col_errGi_prec, covar_re_re_file_in, *covar_re_re_file_prec, covar_re_re_file, covar_im_im_file_in, *covar_im_im_file_prec, covar_im_im_file, covar_re_im_file_in, *covar_re_im_file_prec, covar_re_im_file, col_Gtau_in, *col_Gtau_prec, col_errGtau_in, *col_errGtau_prec, covar_tau_file_in, *covar_tau_file_prec, covar_tau_file, cutoff_wn_in, *cutoff_wn_prec, SW_in, *SW_prec, SC_in, *SC_prec, w_origin_in, *w_origin_prec, step_omega_in, *step_omega_prec, grid_omega_file_in, *grid_omega_file_prec, grid_omega_file, use_grid_params_in, *use_grid_params_prec, omega_grid_params_in, *omega_grid_params_prec, eval_moments_in, *eval_moments_prec, maxM_in, *maxM_prec, def_model_file_in, *def_model_file_prec, def_model_file, init_spectr_func_file_in, *init_spectr_func_file_prec, init_spectr_func_file, interp_type_in, *interp_type_prec, default_model_center_in, *default_model_center_prec, default_model_width_in, *default_model_width_prec, default_model_shape_in, *default_model_shape_prec, non_uniform_grid_in, *non_uniform_grid_prec;
		//peaked_default_model_in, *peaked_default_model_prec,
		
        bool use_grid_params, use_const_dw, use_exp_step, displ_prep_figs, displ_adv_prep_figs, print_other_params, boson, tau_GF, initialize, initialize_maxent, execute_maxent, save_spec_func, print_alpha, displ_optim_figs, cov_diag, moments_provided, eval_moments, covm_diag, wc_exists, w_exists, SW_set, SC_set, peak_exists, read_params, read_other_params, params_loaded, other_params_loaded, M1_set, M2_set, main_spectral_region_set, A_ref_change, show_optimal_alpha_figs, show_lowest_alpha_figs, show_alpha_curves, preproc_complete, Du_constant, non_uniform_grid, w_origin_set, interactive_mode;
		//peaked_default_model, save_Aopt_only,
        
        string interp_type, output_dir_in, output_dir, output_dir_fin, *output_dir_prec, output_name_suffix, output_name_format, w_sample_in, *w_sample_prec, Nalpha_in, alpha_min_in, alpha_init_in, alpha_opt_max_in, alpha_opt_min_in, *alpha_init_prec, alpha_save_max_in, alpha_save_min_in, A_ref_file, A_ref_file_in, *A_ref_file_prec, def_model_output_file_name, A_opt_name_format, A_opt_err_name_format, output_G_format, output_error_format, auto_corr_error_G_format, output_G_opt_format, error_G_opt_format, auto_corr_error_G_opt_format, output_moments_format, output_moments_opt_format, chi2_vs_alpha_format, Asamp_vs_alpha_format, A_opt_name_rm, A_opt_err_name_rm, output_G_opt_rm, error_G_opt_rm, auto_corr_error_G_opt_rm, output_moments_opt_rm;
		//alpha_save_in,
		
        double tem, cutoff_wn, SW, SC, w_origin, step_omega, signG, alpha0, alpha0_default, alpha, pow_alpha_step, alpha_min_default, alpha_min, alpha_opt_max, alpha_opt_min, M0, errM0, M1, errM1, M2, errM2, M3, errM3, std_omega, wl, wr, w0l, w0r, dwl, dwr, dw_peak, M0t, M1n, default_model_width, default_model_center, default_model_shape, dlchi2_lalpha_min, dlchi2_lalpha_max, alpha_save_max, alpha_save_min, lchi2_lalpha_lgth, static_part_G;
		//chi2save, chi2min, alpha_save,
        
        uint col_Gr, col_errGr, col_errGi, col_Gtau, col_errGtau, Nalpha, Nn, Nn_all, indG_0, indG_f, NM, NMinput, NM_odd, NM_even, Nw, NwA, Nwc, jfit, ind_cutoff_wn, NGM, Nalpha_max, NAprec, ind0, Ntau, Nn_as_min;
		uvec n, n_all, Nw_lims;
		int maxM, maxM_default, col_Gi, ind_alpha_vec, NnC, ind_curv, ind_curv0;
        
//        time_t *input_params_file_time, *other_params_file_time;
		
		mat K, KGM, KG_V, KM, KM_V, COV, CRR, CII, CRI, COVM, COVMfit, Ctau, Ctau_all, green_data, error_data, grid_w_data, def_data, Aw_data, Aref_data, Aprec, Aw_samp;
		//HK,
		rowvec omega_grid_params, w_sample;
		uvec w_sample_ind;
		vec Gr, Gi, Gchi2, G_V, GM, wn, wn_all, errGr, errGi, errG, errGtau, M, M_V, errM, M_even, M_odd, Mfit, ws, A, A0, Amin, wc, w, wA, dwS, default_model, w_ref, A_ref, chi2_vec, alpha_vec, S_vec, M_ord, Gtau, tau, dlchi2_lalpha_1, curv_lchi2_lalpha_1, grid_dens;
		cx_vec G, G_all;
		cx_mat Kcx;
		
		bool wn_sign_change, wn_inversed;
		
//		QDateTime *time_params_file, *time_other_params_file;		
		time_t *time_params_file, *time_other_params_file;
		
		//variables defined by Maxime Charlebois:
		vector<vec> vectors_A;  // A
      vector<vec> vectors_w;  // w
      
      FILE * handle_gnuplot;  //handle to talk to gnuplot

    };
    
}

