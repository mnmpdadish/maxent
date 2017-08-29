
#include "maxEnt_data.h"

bool OmegaMaxEnt_data::preproc()
{
	initialize_maxent=true;
	initialize=false;
	bool init_spectrum_exists=false, file_grid_set=false, param_grid_set=false;
	wc_exists=false;
	w_exists=false;
	peak_exists=false;
	main_spectral_region_set=false;
	
	if (!boson) maxM_default=3;
	else maxM_default=1;
	
	if (!maxM_in.size())
		maxM=maxM_default;
	
	int j;
	
	cout<<"\nPREPROCESSING\n\n";
	
	vec data_col1=abs(green_data.col(0));
	if (data_col1.min()==0 && !boson && !tau_GF)
	{
		cout<<"warning: there is a value equal to zero in the first column of your data file. If your data is bosonic, set parameter \"bosonic data\" to \"yes\" in subsection DATA PARAMETERS.\n";
		cout<<"If your data is a function of imaginary time, set \"imaginary time data\" to \"yes\".\n";
		//return false;
	}
	data_col1.reset();
	
	if (!boson)
	{
		ind0=1;
		if (!tau_GF)
		{
			cout<<"fermionic Green function given in Matsubara frequency\n";
			if (!set_G_omega_n_fermions())
			{
				cout<<"Green function definition failed\n";
				return false;
			}
			if (!set_covar_G_omega_n())
			{
				cout<<"covariance definition failed\n";
				return false;
			}
			if (!set_moments_fermions())
			{
				cout<<"moments definition failed\n";
				return false;
			}
			jfit=0;
			if ( (!moments_provided && maxM>0) || eval_moments )
			{
				if (!compute_moments_omega_n())
				{
					cout<<"computation of moments failed\n";
					return false;
				}
			}
			if (std_omega) cout<<"standard deviation of spectrum: "<<std_omega<<endl;
			//cout << "salut1";
		}
		else
		{
			cout<<"fermionic Green function given in imaginary time\n";
			if (!M1_in.size() || !M2_in.size())
			{
				cout<<"For a Green function given in imaginary time, the first two moments are required to obtain the Matsubara frequency Green function.\n";
				cout<<"The program will try to extract those moments from the behavior of the Green function around tau=0 and tau=beta.\n";
			}
			
			tau=green_data.col(0);
			Ntau=tau.n_rows-1;
			Gtau=green_data.col(col_Gtau-1);
			if (Gtau.max()*Gtau.min()<0)
			{
				cout<<"error: G(tau) must not change sign.\n";
				return false;
			}
			if (Gtau(0)>0) Gtau=-Gtau;
			Gchi2=Gtau;
			
			cout<<"Number of imaginary time slices: "<<Ntau<<endl;
			
			M0t=-Gtau(Ntau)-Gtau(0);
			if (M0_in.size())
			{
				M0=stod(M0_in);
				if (abs(M0t-M0)/M0>tol_norm)
					cout<<"warning: -G(0)-G(beta) not equal to provided norm.\n";
			}
			cout<<"-G(0)-G(beta): "<<M0t<<endl;
			
			if (!set_moments_fermions())
			{
				cout<<"moments definition failed\n";
				return false;
			}
			
			if (!tem_in.size())
			{
				tem=1.0/tau(Ntau);
			}
			else
			{
				if (abs(tem-1.0/tau(Ntau))>tol_tem)
					cout<<"warning: provided temperature does not match imaginary time in data file.\n";
				tau=linspace<vec>(0,Ntau,Ntau+1)/(Ntau*tem);
			}
			
			if (!set_covar_Gtau())
			{
				cout<<"imaginary-time covariance matrix definition failed\n";
				return false;
			}
			
			if ( (!moments_provided && maxM>0) || eval_moments )
			{
				if (!compute_moments_tau_fermions())
				{
					cout<<"computation of moments failed\n";
					return false;
				}
			}
			if (std_omega) cout<<"standard deviation of spectrum: "<<std_omega<<endl;
			
			if (!Fourier_transform_G_tau())
			{
				cout<<"error: unable to Fourier transform G(tau)\n";
				return false;
			}
			
		}
		
		if (maxM!=maxM_default)
		{
			if (maxM<M_ord(NM-1) && maxM>=0)
			{
				int j=NM-1;
				while (j>0 && M_ord(j)>maxM) j--;
				NM=j+1;
				M=M.rows(0,j);
				COVM=COVM.submat(0,0,NM-1,NM-1);
				if (covm_diag)	errM=errM.rows(0,NM-1);
				M_ord=M_ord.rows(0,NM-1);
				if (NM<3)
				{
					errM=sqrt(COVM.diag());
					covm_diag=true;
				}
			}
			else if (maxM<0)
			{
				M.reset();
				M_ord.reset();
				errM.reset();
				COVM.reset();
				NM=0;
				covm_diag=true;
			}
		}
		cout<<"number of moments imposed to the spectral function (including normalization): "<<NM<<endl;

		if (SW_set && SC_set)
		{
			wl=SC-SW/2;
			wr=SC+SW/2;
			main_spectral_region_set=true;
		}
		else if (SW_set && M1_set)
		{
			wl=M1-SW/2;
			wr=M1+SW/2;
			main_spectral_region_set=true;
		}
		else if (SC_set && std_omega)
		{
			wl=SC-std_omega;
			wr=SC+std_omega;
			main_spectral_region_set=true;
		}
		else if (std_omega && M1_set)
		{
			wl=M1-std_omega;
			wr=M1+std_omega;
			main_spectral_region_set=true;
		}
		else if (std_omega)
		{
			wl=-std_omega;
			wr=std_omega;
			main_spectral_region_set=true;
		}
		
		test_low_energy_peak_fermions();
		
		Du_constant=true;
		
		if (init_spectr_func_file.size())
		{
			init_spectrum_exists=set_initial_spectrum();
		}
		else if (grid_omega_file.size())
		{
			file_grid_set=set_grid_omega_from_file();
		}
		else if (use_grid_params && omega_grid_params.n_cols>2)
		{
			param_grid_set=set_grid_from_params();
		}
		
		if (!w_exists)
		{
			if (!wc_exists)
			{
				if (!set_wc())
				{
					cout<<"Real frequency grid definition failed.\n";
					return false;
				}
			}
			
			if (!Du_constant)
			{
				if (!set_smooth_step_grid())
				{
					cout<<"Real frequency grid definition failed.\n";
					return false;
				}
			}
			else
			{
				if (!set_omega_grid())
				{
					cout<<"Real frequency grid definition failed.\n";
					return false;
				}
			}
				
		}
		
		
		j=0;
		while (abs(w_origin-w(j+1))<abs(w_origin-w(j)) && j<Nw-2) j++;
		int jc=j;
		
		if (Nw<Nw_max)
		{
			cout<<"number of real frequencies in the grid: "<<Nw<<endl;
			cout<<"minimum frequency: "<<w(0)<<endl;
			cout<<"maximum frequency: "<<w(Nw-1)<<endl;
			cout<<"boundaries of main spectral range: "<<wl<<", "<<wr<<endl;
			cout<<"grid origin: "<<w_origin<<endl;
			cout<<"frequency step at the grid origin: "<<w(jc+1)-w(jc)<<endl;
		}
		else
		{
			cout<<"number of frequencies Nw="<<Nw<<" larger than Nw_max="<<Nw_max<<endl;
			cout<<"Use section \"FREQUENCY GRID PARAMETERS\" in parameter file to modify the number of real frequencies or increase parameter \"Nw_max\" in file \"OmegaMaxEnt_other_params.dat\".\n";
			return false;
		}
		
		if (cutoff_wn_in.size()) // && cutoff_wn<wn(Nn-1))
		{
			//cout<<"cutoff_wn: "<<cutoff_wn<<endl;
			j=Nn-1;
			while (wn(j)>cutoff_wn && j>1) j--;
			ind_cutoff_wn=j;
		}
		else if (jfit)
		{
			ind_cutoff_wn=jfit;
		}
		else
			ind_cutoff_wn=Nn-1;
		
		G_all=G;
		n_all=n;
		wn_all=wn;
		Nn_all=Nn;
		if (ind_cutoff_wn<Nn-1 || ind_cutoff_wn+1>Nw)
		{
			if (!truncate_G_omega_n())
			{
				cout<<"warning: truncation of G failed.\n";
				return false;
			}
		}
		cout<<"Number of Matsubara frequencies used in chi2: "<<Nn<<endl;
		
		if (!set_default_model())
		{
			cout<<"definition of default model failed\n";
			return false;
		}
		
		Kernel_G_fermions_grid_transf_1();
		
		mat norm_DM_tmp=KM.row(0)*default_model;
		default_model=M0*default_model/norm_DM_tmp(0);
		
		diagonalize_covariance();
		
		dwS=(w.rows(2,Nw-1)-w.rows(0,Nw-3))/2;
		
		default_model=exp(1)*(default_model.rows(ind0,Nw-2)+minDefM);
		
		vec A0_default=default_model/exp(1);
		if (!init_spectrum_exists) A0=A0_default;

		vec Pd=sqrt(4*PI*A0_default/dwS);
		mat P=diagmat(Pd);
		mat KGj=KG_V*P;
		mat U, V;
		vec sK, sK2;
		svd(U,sK,V,KGj);
		sK2=pow(sK,2);
		
		alpha0_default=f_alpha_init*sK2.max();
		
		rowvec w0ASmin(2);
		rowvec s0ASmin;
		rowvec wgtASmin=ones<rowvec>(2);
		if (std_omega)
		{
			w0ASmin(0)=SC-std_omega;
			w0ASmin(1)=SC+std_omega;
			s0ASmin=R_width_ASmin*std_omega*ones<rowvec>(2);
		}
		else
		{
			w0ASmin(0)=SC-SW/2;
			w0ASmin(1)=SC+SW/2;
			s0ASmin=R_width_ASmin*(SW/2)*ones<rowvec>(2);
		}
		vec ASmin;
		sum_gaussians(w, w0ASmin, s0ASmin, wgtASmin, ASmin);
		ASmin=ASmin.rows(1,Nw-2)+minDefM;
		
		double Smin=-sum(ASmin % dwS % log(ASmin/default_model))/(2*PI);
		
		if (Smin<0)
			alpha_min_default=2*Nn/(f_Smin*abs(Smin));
		else
			alpha_min_default=2*Nn/f_Smin;
		
		wA=w.rows(1,Nw-2);
		NwA=Nw-2;
		
		if (cov_diag) NnC=Nn/2;
		else NnC=Nn;
		
	}
	else if (col_Gi>0)
	{
		ind0=1;
		if (!tau_GF)
		{
			cout<<"bosonic Green function given in Matsubara frequency\n";
			if (!set_G_omega_n_bosons())
			{
				cout<<"Green function definition failed\n";
				return false;
			}
			if (!set_covar_G_omega_n())
			{
				cout<<"covariance definition failed\n";
				return false;
			}
			if (!set_moments_bosons())
			{
				cout<<"moments definition failed\n";
				return false;
			}
			jfit=0;
			if ( (!moments_provided && maxM>0) || eval_moments )
			{
				if (!compute_moments_omega_n())
				{
					cout<<"computation of moments failed\n";
					return false;
				}
			}
			if (std_omega) cout<<"standard deviation of spectrum: "<<std_omega<<endl;
		}
		else
		{
			cout<<"bosonic Green function given in imaginary time\n";
			if (!M1_in.size() || !M2_in.size())
			{
				cout<<"For a Green function given in imaginary time, the first two moments are required to obtain the Matsubara frequency Green function.\n";
				cout<<"The program will try to extract those moments from the behavior of the Green function around tau=0 and tau=beta.\n";
			}
			
			tau=green_data.col(0);
			Ntau=tau.n_rows-1;
			Gtau=green_data.col(col_Gtau-1);
			if (Gtau.max()*Gtau.min()<0)
			{
				cout<<"error: G(tau) must not change sign.\n";
				return false;
			}
			if (Gtau(0)>0) Gtau=-Gtau;
			Gchi2=Gtau;
			
			cout<<"Number of imaginary time slices: "<<Ntau<<endl;
			
			M0t=Gtau(Ntau)-Gtau(0);
			if (M0_in.size())
			{
				M0=stod(M0_in);
				if (M0t)
				{
					if (abs(M0t-M0)/M0t>tol_norm)
						cout<<"warning: norm of spectral function is different from provided one.\n";
				}
			}
			cout<<"G(0)-G(beta): "<<M0t<<endl;
			
			if (!set_moments_bosons())
			{
				cout<<"moments definition failed\n";
				return false;
			}
			
			if (!tem_in.size())
			{
				tem=1.0/tau(Ntau);
			}
			else
			{
				if (abs(tem-1.0/tau(Ntau))>tol_tem)
					cout<<"warning: provided temperature does not match imaginary time in data file.\n";
				tau=linspace<vec>(0,Ntau,Ntau+1)/(Ntau*tem);
			}
			
			if (!set_covar_Gtau())
			{
				cout<<"imaginary-time covariance matrix definition failed\n";
				return false;
			}
			
			M1n=abs(Gtau(0)/2+sum(Gtau.rows(1,Ntau-1))+Gtau(Ntau)/2)/(Ntau*tem);

			if ( (!moments_provided && maxM>0) || eval_moments )
			{
				if (!compute_moments_tau_bosons())
				{
					cout<<"computation of moments failed\n";
					return false;
				}
			}

			if (!Fourier_transform_G_tau())
			{
				cout<<"error: unable to Fourier transform G(tau)\n";
				return false;
			}

 			M1n=abs(Gr(0));
			if (!std_omega)
			{
				double var_omega=M1/M1n-pow(M0/M1n,2);
				std_omega=sqrt(var_omega);
			}
			
			if (!SC_set)
			{
				SC=M0/M1n;
				SC_set=true;
			}
			if (!SW_set)
			{
				SW=f_SW_std_omega*std_omega;
				SW_set=true;
			}
 			if (std_omega) cout<<"standard deviation of spectrum: "<<std_omega<<endl;

		}
		
		if (maxM!=maxM_default)
		{
			if (maxM<M_ord(NM-1) && maxM>=0)
			{
				int j=NM-1;
				while (j>0 && M_ord(j)>maxM) j--;
				NM=j+1;
				M=M.rows(0,j);
				COVM=COVM.submat(0,0,NM-1,NM-1);
				if (covm_diag)	errM=errM.rows(0,NM-1);
				M_ord=M_ord.rows(0,NM-1);
				if (NM<3)
				{
					errM=sqrt(COVM.diag());
					covm_diag=true;
				}
			}
			else if (maxM<0)
			{
				M.reset();
				M_ord.reset();
				errM.reset();
				COVM.reset();
				NM=0;
			}
		}
		cout<<"number of moments imposed to the spectral function: "<<NM<<endl;
		
		if (SW_set && SC_set)
		{
			wl=SC-SW/2;
			wr=SC+SW/2;
			main_spectral_region_set=true;
		}
		else if (SC_set && std_omega)
		{
			wl=SC-std_omega;
			wr=SC+std_omega;
			main_spectral_region_set=true;
		}
		else if (std_omega)
		{
			wl=-std_omega;
			wr=std_omega;
			main_spectral_region_set=true;
		}
		
		test_low_energy_peak_bosons();
		
		Du_constant=false;
		
		if (init_spectr_func_file.size())
		{
			init_spectrum_exists=set_initial_spectrum();
		}
		else if (grid_omega_file.size())
		{
			file_grid_set=set_grid_omega_from_file();
		}
		else if (use_grid_params && omega_grid_params.n_cols>2)
		{
			param_grid_set=set_grid_from_params();
		}
		
		if (!wc_exists)
		{
			if (!set_wc())
			{
				cout<<"Real frequency grid definition failed.\n";
				return false;
			}
		}
		
		if (!w_exists)
		{
			if (!set_omega_grid())
			{
				cout<<"Real frequency grid definition failed.\n";
				return false;
			}
		}
		
		j=0;
		while (abs(w_origin-w(j+1))<abs(w_origin-w(j)) && j<Nw-2) j=j+1;
		int jc=j;
		
		if (Nw<Nw_max)
		{
			cout<<"number of real frequencies in the grid: "<<Nw<<endl;
			cout<<"minimum frequency: "<<w(0)<<endl;
			cout<<"maximum frequency: "<<w(Nw-1)<<endl;
			cout<<"boundaries of main spectral range: "<<wl<<", "<<wr<<endl;
			cout<<"frequency step at the grid origin: "<<w(jc+1)-w(jc)<<endl;
		}
		else
		{
			cout<<"number of frequencies Nw="<<Nw<<" larger than Nw_max="<<Nw_max<<endl;
			cout<<"Use section \"FREQUENCY GRID PARAMETERS\" in parameter file to modify the number of real frequencies or increase parameter \"Nw_max\" in file \"OmegaMaxEnt_other_params.dat\".\n";
			return false;
		}
		
		if (cutoff_wn_in.size()) // && cutoff_wn<wn(Nn-1))
		{
			j=Nn-1;
			while (wn(j)>cutoff_wn && j>1) j--;
			ind_cutoff_wn=j;
		}
		else if (jfit)
		{
			ind_cutoff_wn=jfit;
		}
		else
			ind_cutoff_wn=Nn-1;
		
		G_all=G;
		n_all=n;
		wn_all=wn;
		Nn_all=Nn;
		if (ind_cutoff_wn<Nn-1 || ind_cutoff_wn+1>Nw)
		{
			if (!truncate_G_omega_n())
			{
				cout<<"warning: truncation of G failed.\n";
				return false;
			}
		}
		cout<<"Number of Matsubara frequencies used in chi2: "<<Nn<<endl;
		
		if (!set_default_model())
		{
			cout<<"definition of default model failed\n";
			return false;
		}
		
		Kernel_G_bosons();
		
		mat norm_DM_tmp=K.row(0)*default_model;
		default_model=M1n*default_model/abs(norm_DM_tmp(0));
		
		diagonalize_covariance();
		
		//HK=2*KGM.t()*KGM;
		
		default_model=exp(1)*(default_model.rows(ind0,Nw-2)+minDefM);
		
		dwS=(w.rows(2,Nw-1)-w.rows(0,Nw-3))/2;
		
		vec A0_default=default_model/exp(1);
		if (!init_spectrum_exists) A0=A0_default;
		
		vec Pd=sqrt(4*PI*A0_default/dwS);
		mat P=diagmat(Pd);
		mat KGj=KG_V*P;
		mat U, V;
		vec sK, sK2;
		svd(U,sK,V,KGj);
		sK2=pow(sK,2);
		
		alpha0_default=f_alpha_init*sK2.max();
		
		//vec Rchi2_S=(HK*default_model)/(2*dwS);
		//alpha0_default=f_alpha_init*max(abs(Rchi2_S));
		
		rowvec w0ASmin(2);
		rowvec s0ASmin;
		rowvec wgtASmin=ones<rowvec>(2);
		if (std_omega)
		{
			w0ASmin(0)=SC-std_omega;
			w0ASmin(1)=SC+std_omega;
			s0ASmin=R_width_ASmin*std_omega*ones<rowvec>(2);
		}
		else
		{
			w0ASmin(0)=SC-SW/2;
			w0ASmin(1)=SC+SW/2;
			s0ASmin=R_width_ASmin*(SW/2)*ones<rowvec>(2);
		}
		vec ASmin;
		sum_gaussians(w, w0ASmin, s0ASmin, wgtASmin, ASmin);
		ASmin=ASmin.rows(ind0,Nw-2)+minDefM;
		
		double Smin=-sum(ASmin % dwS % log(ASmin/default_model))/(2*PI);
		
		if (Smin<0)
			alpha_min_default=2*Nn/(f_Smin*abs(Smin));
		else
			alpha_min_default=2*Nn/f_Smin;
		
		wA=w.rows(1,Nw-2);
		NwA=Nw-2;
		
		if (cov_diag) NnC=Nn/2;
		else NnC=Nn;
	}
	else
	{
		w_origin=0;
		w_origin_set=true;
		ind0=0;
		if (!tau_GF)
		{
			cout<<"bosonic Green function given in Matsubara frequency\n";
			if (!set_G_omega_n_bosons())
			{
				cout<<"Green function definition failed\n";
				return false;
			}
			if (!set_covar_chi_omega_n())
			{
				cout<<"covariance definition failed\n";
				return false;
			}
			if (!set_moments_bosons())
			{
				cout<<"moments definition failed\n";
				return false;
			}
			
			jfit=0;
			if ( (!moments_provided && maxM>0) || eval_moments )
			{
				if (!compute_moments_chi_omega_n())
				{
					cout<<"computation of moments failed\n";
					return false;
				}
			}
			if (std_omega) cout<<"standard deviation of spectrum: "<<std_omega<<endl;
		}
		else
		{
			cout<<"bosonic Green function given in imaginary time\n";
			if (!M1_in.size() || !M2_in.size())
			{
				cout<<"For a Green function given in imaginary time, the first two moments are required to obtain the Matsubara frequency Green function.\n";
				cout<<"The program will try to extract those moments from the behavior of the Green function around tau=0 and tau=beta.\n";
			}
			
			tau=green_data.col(0);
			Ntau=tau.n_rows-1;
			Gtau=green_data.col(col_Gtau-1);
			if (Gtau.max()*Gtau.min()<0)
			{
				cout<<"error: G(tau) must not change sign.\n";
				return false;
			}
			if (Gtau(0)>0) Gtau=-Gtau;
			Gchi2=Gtau;
			
			cout<<"Number of imaginary time slices: "<<Ntau<<endl;
			
			M0t=Gtau(Ntau)-Gtau(0);
			if (M0_in.size())
			{
				M0=stod(M0_in);
				if (M0t)
				{
					if (abs(M0t-M0)/M0t>tol_norm)
						cout<<"warning: norm of spectral function is different from provided one.\n";
				}
			}
			cout<<"G(0)-G(beta): "<<M0t<<endl;
			
			if (!set_moments_bosons())
			{
				cout<<"moments definition failed\n";
				return false;
			}
			
			if (!tem_in.size())
			{
				tem=1.0/tau(Ntau);
			}
			else
			{
				if (abs(tem-1.0/tau(Ntau))>tol_tem)
					cout<<"warning: provided temperature does not match imaginary time in data file.\n";
				tau=linspace<vec>(0,Ntau,Ntau+1)/(Ntau*tem);
			}
			
			if (!set_covar_Gtau())
			{
				cout<<"imaginary-time covariance matrix definition failed\n";
				return false;
			}
			
			M1n=abs(Gtau(0)/2+sum(Gtau.rows(1,Ntau-1))+Gtau(Ntau)/2)/(Ntau*tem);
		
			if ( (!moments_provided && maxM>0) || eval_moments )
			{
				if (!compute_moments_tau_bosons())
				{
					cout<<"computation of moments failed\n";
					return false;
				}
			}
		
			if (!Fourier_transform_G_tau())
			{
				cout<<"error: unable to Fourier transform G(tau)\n";
				return false;
			}
			
			M1n=abs(Gr(0));
			if (!std_omega)
			{
				double var_omega=M1/M1n-pow(M0/M1n,2);
				std_omega=sqrt(var_omega);
			}
			
			if (!SC_set)
			{
				SC=0;
				SC_set=true;
			}
			if (!SW_set)
			{
				SW=f_SW_std_omega*std_omega;
				SW_set=true;
			}
			if (std_omega) cout<<"standard deviation of spectrum: "<<std_omega<<endl;
		}
		
		if (maxM!=maxM_default)
		{
			if (maxM<M_ord(NM-1) && maxM>0)
			{
				int j=NM-1;
				while (j>0 && M_ord(j)>maxM) j--;
				NM=j+1;
				M=M.rows(0,j);
				COVM=COVM.submat(0,0,NM-1,NM-1);
				if (covm_diag)	errM=errM.rows(0,NM-1);
				M_ord=M_ord.rows(0,NM-1);
				if (NM<3)
				{
					errM=sqrt(COVM.diag());
					covm_diag=true;
				}
			}
			else if (maxM<=0)
			{
				M.reset();
				M_ord.reset();
				errM.reset();
				COVM.reset();
				NM=0;
			}
		}
		cout<<"number of moments imposed to the spectral function: "<<NM<<endl;
		
		if (SW_set)
		{
			wr=SW/2;
			main_spectral_region_set=true;
		}
		else if (std_omega)
		{
			wr=std_omega/2;
			main_spectral_region_set=true;
		}
		
		test_low_energy_peak_chi();
		
		Du_constant=false;
		
		if (init_spectr_func_file.size())
		{
			init_spectrum_exists=set_initial_spectrum_chi();
		}
		else if (grid_omega_file.size())
		{
			file_grid_set=set_grid_omega_from_file_chi();
		}
		else if (use_grid_params && omega_grid_params.n_cols>2)
		{
			param_grid_set=set_grid_from_params_chi();
		}
		
		if (!wc_exists)
		{
			if (!set_wc_chi())
			{
				cout<<"Real frequency grid definition failed.\n";
				return false;
			}
		}
		
		if (!w_exists)
		{
			if (!set_omega_grid_chi())
			{
				cout<<"Real frequency grid definition failed.\n";
				return false;
			}
		}
		
		if (Nw<Nw_max)
		{
			cout<<"number of real frequencies in the grid: "<<Nw<<endl;
			cout<<"minimum frequency: "<<w(0)<<endl;
			cout<<"maximum frequency: "<<w(Nw-1)<<endl;
			cout<<"boundaries of main spectral range:  0, "<<wr<<endl;
			cout<<"frequency step at the grid origin: "<<w(1)-w(0)<<endl;
		}
		else
		{
			cout<<"number of frequencies Nw="<<Nw<<" larger than Nw_max="<<Nw_max<<endl;
			cout<<"Use section \"FREQUENCY GRID PARAMETERS\" in parameter file to modify the number of real frequencies or increase parameter \"Nw_max\" in file \"OmegaMaxEnt_other_params.dat\".\n";
			return false;
		}
		
		if (cutoff_wn_in.size()) // && cutoff_wn<wn(Nn-1))
		{
			j=Nn-1;
			while (wn(j)>cutoff_wn && j>1) j--;
			ind_cutoff_wn=j;
		}
		else if (jfit)
		{
			ind_cutoff_wn=jfit;
		}
		else
			ind_cutoff_wn=Nn-1;
		
		G_all=G;
		n_all=n;
		wn_all=wn;
		Nn_all=Nn;
		if (ind_cutoff_wn<Nn-1 || ind_cutoff_wn+1>Nw)
		{
			if (!truncate_chi_omega_n())
			{
				cout<<"warning: truncation of G failed.\n";
				return false;
			}
		}
		cout<<"Number of Matsubara frequencies used in chi2: "<<Nn<<endl;
		
		if (!set_default_model_chi())
		{
			cout<<"definition of default model failed\n";
			return false;
		}
		
		Kernel_chi();
		
		mat norm_DM_tmp=K.row(0)*default_model;
		default_model=M1n*default_model/abs(norm_DM_tmp(0));
		
		diagonalize_covariance_chi();
		
		//HK=2*KGM.t()*KGM;
		
		default_model=exp(1)*(default_model.rows(ind0,Nw-2)+minDefM);
		
		dwS.zeros(Nw-1);
		dwS.rows(1,Nw-2)=w.rows(2,Nw-1)-w.rows(0,Nw-3);
		dwS(0)=w(1)-w(0);
		
		vec A0_default=default_model/exp(1);
		if (!init_spectrum_exists) A0=A0_default;
		
		vec Pd=sqrt(4*PI*A0_default/dwS);
		mat P=diagmat(Pd);
		mat KGj=KG_V*P;
		mat U, V;
		vec sK, sK2;
		svd(U,sK,V,KGj);
		sK2=pow(sK,2);
		
		alpha0_default=f_alpha_init*sK2.max();
		
		//vec Rchi2_S=(HK*default_model)/(2*dwS);
		//alpha0_default=f_alpha_init*max(abs(Rchi2_S));

		rowvec w0ASmin(1);
		rowvec s0ASmin(1);
		rowvec wgtASmin(1);
		wgtASmin(0)=1;
		if (std_omega)
		{
			w0ASmin(0)=std_omega;
			s0ASmin(0)=R_width_ASmin*std_omega;
		}
		else
		{
			w0ASmin(0)=SW/2;
			s0ASmin(0)=R_width_ASmin*(SW/2);
		}
		vec ASmin;
		sum_gaussians_chi(w, w0ASmin, s0ASmin, wgtASmin, ASmin);
		ASmin=ASmin.rows(ind0,Nw-2)+minDefM;
		
		double Smin=-sum(ASmin % dwS % log(ASmin/default_model))/(2*PI);
		
		if (Smin<0)
			alpha_min_default=Nn/(f_Smin*abs(Smin));
		else
			alpha_min_default=Nn/f_Smin;
		
		wA=w.rows(ind0,Nw-2);
		NwA=Nw-1;
		
		if (cov_diag) NnC=Nn/2;
		else NnC=Nn;
	}
	
	return true;
}











void OmegaMaxEnt_data::minimize()
{
	char alpha_output[100], file_name[200];
	char alpha_output_format[]="%d \t alpha: % 1.4e,  Q: % 1.4e,  S: % 1.4e,  chi2: % 1.4e\n";
	double mean_int_dA_prec, mean_int_dA_prec2, A1min, chi2prec, Q, S;
	mat chi2, P, KGMj, U, V, mean_int_dA, M_save;
	vec A1, AS, c1, c2, grS, B, B2, Pd, sK, sK2(A.n_rows), D1, dA1, VPdA, dA_rel, DG;
	vec G_out, G_V_out, errIm, errRe, M_out, M_V_out, eigv_ind;
	uvec ind_c2_sat, ind_An, ind_Anul;
	int ind_alpha, iter_dA;
	
	int Nalpha_min=200;
	
	bool alpha_too_small=false;
	
	double pow_alpha0=log10(alpha0);
	double pow_alpha=log10(alpha);
	double alpha_prec=pow(10,pow_alpha+pow_alpha_step);
	
	double realmin=1e-100;
	
	double rADmin=realmin;
	vec Amin=rADmin*default_model;
	
	double rADchange=realmin;
	vec Achange=rADchange*default_model;
	
	svd(U,sK,V,KG_V);
	sK2.zeros();
	sK2.rows(0,sK.n_rows-1)=pow(sK,2);
	double alpha_c2_max=sK2.max()*rc2H;
	
//	mat absMK=abs(KGM.t()*KGM);
//	double alpha_c2_max=absMK.max()*rc2H;
	
//	cout<<"sK2.max(): "<<sK2.max()<<endl;
//	cout<<"alpha_c2_max: "<<alpha_c2_max<<endl;
		
	cout<<"Computing spectrum as a function of alpha...\n";
	
	vec lalpha, lchi2, arc_params;
	int j, jmin, jmax;
	double fs=f_scale_lalpha_lchi2;
	double smooth_length=chi2_alpha_smooth_range;
	double la, lc, ls;
	double xc, yc, x1, y1;
	
	c1=(log(rADchange)+1)*dwS;
	//vec c2_0=-fc2*c1/(2*exp(-1)*default_model);
	vec c2_alpha=2*PI*alpha_c2_max*ones<vec>(NwA);
	
	uword ind_max_dchi2_alpha;
	
	ind_alpha=1;
	DG=GM-KGM*A;
	chi2=DG.t()*DG;
	while (ind_alpha<=Nalpha && alpha>=alpha_min)
	{
		A1=A;
		
		c2=c2_alpha/alpha;
		//c2=c2_0;
		//ind_c2_sat=find((alpha*c2)/(2*PI)>alpha_c2_max);
		//if (ind_c2_sat.n_rows)
		//	c2.rows(ind_c2_sat)=(2*PI*alpha_c2_max/alpha)*ones<vec>(ind_c2_sat.n_rows);
		
		grS=dwS % log(A1/default_model) + dwS;
		ind_An=find(A1/default_model<rADchange);
		if (ind_An.n_rows)
		{
			grS.rows(ind_An)=2*c2.rows(ind_An) % (A1.rows(ind_An)-Achange.rows(ind_An))+c1.rows(ind_An);
		}
		
		B=KGM.t()*(GM-KGM*A1)-(alpha/(4*PI))*grS;
		
		Pd=sqrt(4*PI*A1/(alpha*dwS));
		if (ind_An.n_rows)
			Pd.rows(ind_An)=sqrt(2*PI/(alpha*c2(ind_An)));
		P=diagmat(Pd);
		
		KGMj=KGM*P;
		svd(U,sK,V,KGMj);
		
		B2=V.t()*(P*B);
		sK2.zeros();
		sK2.rows(0,sK.n_rows-1)=pow(sK,2);
		D1=sK2+1;
		VPdA=B2/D1;
		dA1=P*(V*VPdA);
		
		mean_int_dA=abs(dA1.t())*dwS;
		mean_int_dA_prec=2*mean_int_dA(0);
		mean_int_dA_prec2=mean_int_dA_prec;
		
		iter_dA=1;
		
		A1=A1+dA1;
		
		ind_Anul=find(A1==0);
		if (ind_Anul.n_rows)
		{
			A1.rows(ind_Anul)=Amin.rows(ind_Anul);
		//	A1.rows(ind_Anul)=arma::max(Amin.rows(ind_Anul),DBL_MIN);
		}
		
		A1min=min(A1-Amin);
		
		while ( (mean_int_dA(0)>tol_int_dA || A1min<0) && iter_dA<Niter_dA_max && (mean_int_dA(0)<mean_int_dA_prec || mean_int_dA(0)<mean_int_dA_prec2))
		{
			grS=dwS % log(A1/default_model) + dwS;
			ind_An=find(A1/default_model<rADchange);
			if (ind_An.n_rows)
			{
				grS.rows(ind_An)=2*c2.rows(ind_An) % (A1.rows(ind_An)-Achange.rows(ind_An))+c1.rows(ind_An);
			}
			
			B=KGM.t()*(GM-KGM*A1)-(alpha/(4*PI))*grS;
			
			Pd=sqrt(4*PI*A1/(alpha*dwS));
			if (ind_An.n_rows)
				Pd.rows(ind_An)=sqrt(2*PI/(alpha*c2(ind_An)));
			P=diagmat(Pd);
			
			KGMj=KGM*P;
			svd(U,sK,V,KGMj);
			
			B2=V.t()*(P*B);
			sK2.zeros();
			sK2.rows(0,sK.n_rows-1)=pow(sK,2);
			D1=sK2+1;
			VPdA=B2/D1;
			dA1=P*(V*VPdA);
			
			mean_int_dA_prec2=mean_int_dA_prec;
			mean_int_dA_prec=mean_int_dA(0);
			mean_int_dA=abs(dA1.t())*dwS;
			
			if (mean_int_dA(0)<mean_int_dA_prec)
			{
				A1=A1+dA1;
				A1min=min(A1-Amin);
			}
			
			ind_Anul=find(A1==0);
			if (ind_Anul.n_rows)
			{
				A1.rows(ind_Anul)=Amin.rows(ind_Anul);
			//	A1.rows(ind_Anul)=arma::max(Amin.rows(ind_Anul),DBL_MIN);
			}
			
			iter_dA++;
		}

/*
		if (mean_int_dA(0)>mean_int_dA_prec || iter_dA==Niter_dA_max)
		{
			cout<<"alpha= "<<alpha<<"  iteration: "<<iter_dA<<"  integrated absolute variation in A: "<<mean_int_dA_prec<<endl;
		}
*/		
		chi2prec=chi2(0);
		DG=GM-KGM*A1;
		chi2=DG.t()*DG;
		AS=A1;
		ind_An=find(AS/default_model<rADchange);
		if (ind_An.n_rows)
			AS.rows(ind_An)=Amin.rows(ind_An);
		S=-sum(AS % dwS % log(AS/default_model))/(2*PI);
		
		if (abs((chi2(0)-chi2prec)/chi2prec)>diff_chi2_max && alpha<alpha_prec && pow_alpha_step>2*pow_alpha_step_min)
			pow_alpha_step=pow_alpha_step/2;
		
		if ((chi2(0)<chi2prec && alpha<=alpha_prec) || (chi2(0)>chi2prec && alpha>alpha_prec) || pow_alpha==pow_alpha0)
		{
			Aprec.cols(0,NAprec-2)=Aprec.cols(1,NAprec-1);
			Aprec.col(NAprec-1)=A;
			A=A1;
			alpha_vec(ind_alpha_vec)=alpha;
			chi2_vec(ind_alpha_vec)=chi2(0);
			S_vec(ind_alpha_vec)=S;
			vec logA=log(A/default_model);
			ind_An=find(A/default_model<rADchange);
			if (ind_An.n_rows)
			{
				logA.rows(ind_An)=log(rADmin)*ones<vec>(ind_An.n_rows);
			}
			S_vec(ind_alpha_vec)=-sum(A % dwS % logA)/(2*PI);
			Aw_samp.row(ind_alpha_vec)=trans(A(w_sample_ind));
			
			if (print_alpha)
			{
				Q=chi2(0)-alpha*S_vec(ind_alpha_vec);
				sprintf(alpha_output,alpha_output_format,ind_alpha,alpha,Q,S_vec(ind_alpha_vec),chi2(0));
				cout<<alpha_output;
			}
			
			G_out=K*A;
			G_V_out=KG_V*A;
			
			if (!boson || col_Gi>0)
			{
				if (cov_diag)
				{
					uvec even_ind=linspace<uvec>(0,2*Nn-2,Nn);
					errRe=G_V.rows(even_ind)-G_V_out.rows(even_ind);
					errIm=G_V.rows(even_ind+1)-G_V_out.rows(even_ind+1);
				}
				else
				{
					errRe=G_V-G_V_out;
					//errRe=G_V.rows(0,Nn-1)-G_V_out.rows(0,Nn-1);
					//errIm=G_V.rows(Nn,2*Nn-1)-G_V_out.rows(Nn,2*Nn-1);
				}
			}
			else
			{
				errRe=G_V-G_V_out;
			}
			
			if (cov_diag)
				eigv_ind=wn;
			else
				eigv_ind=linspace<vec>(0,errRe.n_rows-1,errRe.n_rows);
			
			if (NM>0)
			{
				M_out=KM*A;
				M_V_out=KM_V*A;
			}
			
			
			vectors_A.push_back(A);
			vectors_w.push_back(w);
			
			/*
			if (save_spec_func)
			{
				
				
				M_save.zeros(Nw,2);
				M_save.col(0)=w;
				if (!boson || col_Gi>0)
					M_save.submat(1,1,Nw-2,1)=A;
				else
					M_save.submat(0,1,Nw-2,1)=A;
				sprintf(file_name,output_name_format.c_str(),tem,alpha);
				//M_save.save(file_name,raw_ascii);
				
				
				if (!boson || col_Gi>0)
				{
					M_save.zeros(Nn,3);
					uvec even_ind=linspace<uvec>(0,2*Nn-2,Nn);
					M_save.col(0)=wn;
					M_save.col(1)=G_out.rows(even_ind);
					M_save.col(2)=G_out.rows(even_ind+1);
				}
				else
				{
					M_save.zeros(Nn,2);
					M_save.col(0)=wn;
					M_save.col(1)=G_out;
				}
				sprintf(file_name,output_G_format.c_str(),tem,alpha);
				//M_save.save(file_name,raw_ascii);
				
				if (!boson || col_Gi>0)
				{
					if (cov_diag)
					{
						M_save.zeros(Nn,3);
						M_save.col(0)=wn;
						M_save.col(1)=errRe;
						M_save.col(2)=errIm;
					}
					else
					{
						M_save.zeros(2*Nn,2);
						M_save.col(0)=eigv_ind;
						M_save.col(1)=errRe;
					}
				}
				else
				{
					M_save.zeros(Nn,2);
					M_save.col(0)=eigv_ind;
					M_save.col(1)=errRe;
				}
				sprintf(file_name,output_error_format.c_str(),tem,alpha);
				//M_save.save(file_name,raw_ascii);
				
				if (NM>0)
				{
					sprintf(file_name,output_moments_format.c_str(),tem,alpha);
					M_save.zeros(NM,2);
					M_save.col(0)=M;
					M_save.col(1)=M_out;
					//M_save.save(file_name,raw_ascii);
				}
				
			}*/
			
			pow_alpha=pow_alpha-pow_alpha_step;
			alpha_prec=alpha;
			alpha=pow(10,pow_alpha);
		}
		else
		{
			cout<<"minimize(): change in chi^2 is in the wrong direction\n"; //mc: wtf?
			break;
		}
		
		if (ind_alpha_vec>2)
		{
			if (ind_curv0==0)
			{
				lalpha=log10(alpha_vec.rows(0,ind_alpha_vec));
				lchi2=log10(chi2_vec.rows(0,ind_alpha_vec));
				
				j=1;
				la=lalpha(j-1)-lalpha(j);
				lc=lchi2(j-1)-lchi2(j);
				ls=sqrt(pow(fs*la,2)+pow(lc,2));
				while (ls<fs*smooth_length && j<ind_alpha_vec)
				{
					j++;
					la=lalpha(j-1)-lalpha(j);
					lc=lchi2(j-1)-lchi2(j);
					ls=ls+sqrt(pow(fs*la,2)+pow(lc,2));
				}
				if (ls>=fs*smooth_length)	ind_curv0=j;
			}
 
			if (ind_curv0 && ind_curv0<ind_alpha_vec)
			{
				if (ind_curv>=ind_curv0)
				{
					lalpha=log10(alpha_vec.rows(0,ind_alpha_vec));
					lchi2=log10(chi2_vec.rows(0,ind_alpha_vec));
					
					ind_curv++;
					jmin=ind_curv-1;
					la=lalpha(jmin)-lalpha(jmin+1);
					lc=lchi2(jmin)-lchi2(jmin+1);
					ls=sqrt(pow(fs*la,2)+pow(lc,2));
					while (ls<fs*smooth_length && jmin>0)
					{
						jmin--;
						la=lalpha(jmin)-lalpha(jmin+1);
						lc=lchi2(jmin)-lchi2(jmin+1);
						ls=ls+sqrt(pow(fs*la,2)+pow(lc,2));
					}
					jmax=ind_curv+1;
					la=lalpha(jmax-1)-lalpha(jmax);
					lc=lchi2(jmax-1)-lchi2(jmax);
					ls=sqrt(pow(fs*la,2)+pow(lc,2));
					while (ls<fs*smooth_length && jmax<ind_alpha_vec)
					{
						jmax++;
						la=lalpha(jmax-1)-lalpha(jmax);
						lc=lchi2(jmax-1)-lchi2(jmax);
						ls=ls+sqrt(pow(fs*la,2)+pow(lc,2));
					}
					fit_circle_arc(fs*lalpha.rows(jmin,jmax), lchi2.rows(jmin,jmax), arc_params);
					curv_lchi2_lalpha_1(ind_curv)=1.0/arc_params(0);
					
					xc=arc_params(1);
					yc=arc_params(2);
					x1=fs*lalpha(ind_curv);
					y1=lchi2(ind_curv);
					
					dlchi2_lalpha_1(ind_curv)=(xc-x1)/(y1-yc);
				}
				else
				{
					lalpha=log10(alpha_vec.rows(0,ind_alpha_vec));
					lchi2=log10(chi2_vec.rows(0,ind_alpha_vec));
					
					j=ind_alpha_vec-1;
					la=lalpha(j)-lalpha(j+1);
					lc=lchi2(j)-lchi2(j+1);
					ls=sqrt(pow(fs*la,2)+pow(lc,2));
					while (ls<fs*smooth_length && j>ind_curv0)
					{
						j--;
						la=lalpha(j)-lalpha(j+1);
						lc=lchi2(j)-lchi2(j+1);
						ls=ls+sqrt(pow(fs*la,2)+pow(lc,2));
					}
					if (ls>=fs*smooth_length)	ind_curv=j;
				}
			}
			
			if (ind_curv0 && ind_curv>=ind_curv0)
			{
				dlchi2_lalpha_max=dlchi2_lalpha_1.max(ind_max_dchi2_alpha);
				dlchi2_lalpha_min=dlchi2_lalpha_1(ind_curv);
				
				if (alpha<alpha_min)
				{
					if (dlchi2_lalpha_min/dlchi2_lalpha_max>RMAX_dlchi2_lalpha || ind_alpha_vec<Nalpha_min)
					{
						if (!alpha_min_in.size())
						{
							//cout<<"Reducing alpha_min by a factor "<<f_alpha_min<<".\n";
							alpha_min=alpha_min/f_alpha_min;
						}
						else
						{
							cout<<"warning: minimum value of alpha reached, but the calculation does not seem over.\n";
						}
					}
				}
				else if (dlchi2_lalpha_min/dlchi2_lalpha_max<RMAX_dlchi2_lalpha/10 && ind_alpha_vec>Nalpha_min)
				{
					if (!alpha_min_in.size())
					{
						cout<<"dlog(chi2)/dlog(alpha) is below the minimum value. Stopping calculation.\n";
						alpha_min=alpha_prec;
					}
					else if (!alpha_too_small)
					{
						alpha_too_small=true;
						cout<<"warning: minimum value of alpha seems too small.\n";
					}
				}
			}
		}
		
		ind_alpha++;
		ind_alpha_vec++;
	}
	
	if (ind_alpha>Nalpha)
		cout<<"maximum number of values alpha reached\n";
	else if (alpha<alpha_min)
		cout<<"minimum value of alpha reached\n";
}






