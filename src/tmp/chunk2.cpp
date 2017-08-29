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
	
