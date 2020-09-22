#include "maxEnt_data.h"

//Contains constructor and destructor.


OmegaMaxEnt_data::OmegaMaxEnt_data(int arg_N, char *args[])
{
    //input_params_file_name.assign("para.dat");
    //other_params_file_name.assign("OmegaMaxEnt_other_params.dat");
	
	if (arg_N>1)
		fileName = args[1];
	else
		exit(1);

	
	//params_loaded=false;
	//other_params_loaded=false;
	preproc_complete=false;
	initialize=true;
	initialize_maxent=true;
	//print_other_params=false;
	interactive_mode=false;
	init_variables();
	ind_alpha_vec=0;
	print_alpha=false;



	col_Gr=2;
	col_Gi=3;
	col_errGr=2;
	col_errGi=3;
	col_Gtau=2;
	col_errGtau=3;
	
	tau_GF=false;
	non_uniform_grid=false;
	use_grid_params=false;
	eval_moments=false;
	default_model_center=0;
	default_model_width=1.0;
	default_model_shape=2.0;
	

	//additional stuff
	Nalpha_max=500;
					
	
	// 
	SC_set=false;
	SW_set=false;
	M1_set=false;
	M2_set=false;
	execute_maxent=true;
	A_ref_change=false;
	w_origin_set=false;
	
	
	// so called "other params_int" initialization;
	Nw_min=50;
   Nw_max=800;
   Nn_fit_max=200;
   Nn_fit_fin=300;
	Nn_as_min=10;
   Niter_dA_max=20;
   Nwsamp=11;
	
	// so called "other params_fl" initialization;
	f_SW_std_omega=3;
	f_w_range=20.0;
	Rmin_SW_dw=25;
	tol_tem=1.0e-8;
	tol_norm=0.1;
	tol_M1=0.05;
	tol_M2=0.05;
	tol_M3=0.1;
	default_error_G=1.0e-4;
	err_norm=1.0e-6;
	default_error_M=1.0e-4;
	tol_mean_C1=0.001;
	tol_std_C1=0.01;
	tol_rdw=1.0e-10;
	Rmin_Dw_dw=4.0;
	Rdw_max=5.0;
	RW_grid=5.0;
	RWD_grid=10.0;
	minDefM=1.0e-20;
	f_alpha_init=1.0e3;
	R_width_ASmin=0.05;
	f_Smin=1.0;
	diff_chi2_max=0.1;
	//tol_int_dA=1.0e-12;    // this is too intense for nothing
	//tol_int_dA=5e-2;
	tol_int_dA=1.0e-4;
	rc2H=1.0e12;
	fc2=1.0e20;
	pow_alpha_step_init=0.2;
	pow_alpha_step_min=0.001;
	chi2_alpha_smooth_range=0.2;
	f_scale_lalpha_lchi2=0.2;
	FNfitTauW=4.0;
	std_norm_peak_max=0.01;
	varM2_peak_max=0.01;
	peak_weight_min=1.0e-4;
	RMAX_dlchi2_lalpha=0.01;
	f_alpha_min=100;
	save_alpha_range=1;
	R_peak_width_dw=10;
	R_wncutoff_wr=10;
	R_Dw_dw=20;
	R_SW_wr=1;
	R_wmax_wr_min=3;
	
	handle_gnuplot = gpc_init_image();
}

OmegaMaxEnt_data::~OmegaMaxEnt_data()
{
      //free_variables();
}

