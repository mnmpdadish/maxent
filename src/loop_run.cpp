#include "maxEnt_data.h"


void OmegaMaxEnt_data::loop_run()
{
	struct stat file_stat;
	char continue_exec;
	string buf;
	save_spec_func=true;
	
	alpha_save_max=DBL_MIN;
	alpha_save_min=DBL_MIN;

	
	read_params=false;
	read_other_params=false;
	cutoff_wn=0;

	if (load_data_file(green_data, fileName)){
		cout<<"data file loaded\n";
	}
	else{
		cout<<"cannot load " << fileName << "\n";
		exit(1);
	}
	
	/// do the preprocessing
	preproc();
	
	
	double pow_alpha0=log10(alpha0), pow_alpha_min=log10(alpha_min);
	
	if (!alpha_opt_max_in.size())
		alpha_opt_max=alpha0;
	
	if (!alpha_opt_min_in.size())
		alpha_opt_min=DBL_MIN;
	
	
	if (w_sample_in.size())
		Nwsamp=w_sample.n_cols;
	else
	{
		if (col_Gi>0)
		{
			double tol_wsamp=1e-2;
			rowvec w_sample_tmp=linspace<rowvec>(wl,wr,Nwsamp);
			int j=0;
			while (j<Nwsamp && w_sample_tmp(j)<0) j++;
			if (j<Nwsamp)
			{
				if (w_sample_tmp(j)/SW<tol_wsamp)
				{
					w_sample_tmp(j)=0;
					w_sample=w_sample_tmp;
				}
				else if (j>0 && fabs(w_sample_tmp(j-1))/SW<tol_wsamp)
				{
					w_sample_tmp(j-1)=0;
					w_sample=w_sample_tmp;
				}
				else
				{
					w_sample.zeros(Nwsamp+1);
					if (j>0)	w_sample.cols(0,j-1)=w_sample_tmp.cols(0,j-1);
					w_sample.cols(j+1,Nwsamp)=w_sample_tmp.cols(j,Nwsamp-1);
					w_sample(j)=0;
				}
			}
			else if (fabs(w_sample_tmp(j-1))/SW<tol_wsamp)
			{
				w_sample_tmp(j-1)=0;
				w_sample=w_sample_tmp;
			}
			else
				w_sample=w_sample_tmp;
			
		}
		else
			w_sample=linspace<rowvec>(0,wr,Nwsamp);
	}
	
	if (col_Gi>0)
		w_sample_ind.ones(Nwsamp);
	else
		w_sample_ind.zeros(Nwsamp);
	
	if (col_Gi>0)
		w_sample_ind=2*w_sample_ind;
	
	int indmax=Nw-3;
	if (col_Gi<=0)	indmax=Nw-2;
	int j;
	for (j=0; j<Nwsamp; j++)
		while (abs(w_sample(j)-w(w_sample_ind(j)+1))<abs(w_sample(j)-w(w_sample_ind(j))) && w_sample_ind(j)<indmax)	w_sample_ind(j)++;
	
	w_sample=trans(w.rows(w_sample_ind));
	
	if (col_Gi>0) w_sample_ind=w_sample_ind-1;
	
	if (initialize_maxent)
	{
		initialize_maxent=false;
		
		if (!alpha_init_in.size()) alpha0=alpha0_default;
		if (!alpha_min_in.size()) alpha_min=alpha_min_default;
		
		
		/*
		if (output_dir_in.size())
		{
			if (stat(output_dir_in.c_str(),&file_stat))
			{
				cout<<"creating output directory: "<<output_dir_in<<endl;
				mkdir(output_dir_in.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
			}
		}
		if (stat(output_dir.c_str(),&file_stat))
		{
			cout<<"creating output directory: "<<output_dir<<endl;
			mkdir(output_dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
		}
		if (stat(output_dir_fin.c_str(),&file_stat))
		{
			cout<<"creating output directory: "<<output_dir_fin<<endl;
			mkdir(output_dir_fin.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
		}
		*/
		
		mat save_mat=zeros<mat>(Nw,2);
		save_mat.col(0)=w;
		if (col_Gi>0)
			save_mat.submat(1,1,Nw-2,1)=default_model;
		else
			save_mat.submat(0,1,Nw-2,1)=default_model;
		string file_name(output_dir);
		file_name+=def_model_output_file_name;
		//save_mat.save(file_name.c_str(),raw_ascii);
		file_name.assign(output_dir_fin);
		file_name+=def_model_output_file_name;
		//save_mat.save(file_name.c_str(),raw_ascii);
		
		A=A0;
		Aprec=A0*ones<rowvec>(NAprec);
		alpha=alpha0;
		pow_alpha_step=pow_alpha_step_init;
		
		alpha_vec.zeros(Nalpha_max);
		chi2_vec.zeros(Nalpha_max);
		S_vec.zeros(Nalpha_max);
		Aw_samp.zeros(Nalpha_max,Nwsamp);
		ind_alpha_vec=0;
		
		lchi2_lalpha_lgth=0;
		ind_curv=0;
		ind_curv0=0;
		curv_lchi2_lalpha_1.zeros(Nalpha_max);
		dlchi2_lalpha_1.zeros(Nalpha_max);
		
		cout<<"\nStarting minimization at initial value of alpha\n";
		
	}
	
	//if (read_params) copy_file(input_params_file_name, "./", output_dir_fin);
	//if (read_other_params) copy_file(other_params_file_name, "./", output_dir_fin);
	
	double pow_alpha=log10(alpha);
	//if (pow_alpha_min<=pow_alpha)
	//	Nalpha_max=(pow_alpha-pow_alpha_min)/pow_alpha_step_min;
	//else
	//	Nalpha_max=500;
	//if (!Nalpha_in.size()) 
	Nalpha=Nalpha_max;
	//Nalpha=500;
	int Nalpha_new=ind_alpha_vec+Nalpha;
	if (alpha_vec.n_rows<Nalpha_new)
	{
		alpha_vec.resize(Nalpha_new);
		chi2_vec.resize(Nalpha_new);
		S_vec.resize(Nalpha_new);
		Aw_samp.resize(Nalpha_new,Nwsamp);
		curv_lchi2_lalpha_1.resize(Nalpha_new);
		dlchi2_lalpha_1.resize(Nalpha_new);
		//MC!
		vectors_A.resize(Nalpha_new);
		vectors_w.resize(Nalpha_new);
	}
	
	handle_gnuplot = gpc_init_image();
	//printf("handle: %d\n", handle_gnuplot);
    
	minimize();
//				pow_alpha_step=pow_alpha_step_min;
//				minimize_increase_alpha();
	
	char alpha_output[100], alpha_output_format[]="alpha: % 1.4e,  Q: % 1.4e,  S: % 1.4e,  chi2: % 1.4e\n";
	//cout << "salut3 " << alpha_vec.n_rows << " " << S_vec.n_rows << " " << chi2_vec.n_rows << " " << ind_alpha_vec << " " <<Nalpha_new << " " <<Nalpha <<" " <<Nalpha_max;
	double Q=chi2_vec(ind_alpha_vec-1)-alpha_vec(ind_alpha_vec-1)*S_vec(ind_alpha_vec-1);
	sprintf(alpha_output,alpha_output_format,alpha_vec(ind_alpha_vec-1),Q,S_vec(ind_alpha_vec-1),chi2_vec(ind_alpha_vec-1));
	cout<<"last output values:\n";
	cout<<alpha_output;
	
	vec lalpha=log10(alpha_vec.rows(0,ind_alpha_vec-1));
	vec lchi2=log10(chi2_vec.rows(0,ind_alpha_vec-1));
	
	char file_name[200];
	sprintf(file_name, chi2_vs_alpha_format.c_str(),tem);
	mat save_mat=join_rows(alpha_vec.rows(0,ind_alpha_vec-1),chi2_vec.rows(0,ind_alpha_vec-1));
	string name(output_dir);
	name.append(file_name);
	//save_mat.save(name.c_str(), raw_ascii);
	name.assign(output_dir_fin);
	name.append(file_name);
	//save_mat.save(name.c_str(), raw_ascii);

	sprintf(file_name, Asamp_vs_alpha_format.c_str(),tem);
	save_mat=join_rows(alpha_vec.rows(0,ind_alpha_vec-1),Aw_samp.rows(0,ind_alpha_vec-1));
	name.assign(output_dir);
	name.append(file_name);
	//save_mat.save(name.c_str(), raw_ascii);
	name.assign(output_dir_fin);
	name.append(file_name);
	//save_mat.save(name.c_str(), raw_ascii);

	vec DM=default_model/exp(1);
	
	if (ind_alpha_vec>3)
	{
		//cout<<"ind_alpha_vec: "<<ind_alpha_vec<<endl;
		
		double fs=f_scale_lalpha_lchi2;
		double smooth_length=chi2_alpha_smooth_range;
		int ind_curv_start, ind_curv_end;
		double la, lc, ls;
		j=1;
		la=lalpha(j-1)-lalpha(j);
		lc=lchi2(j-1)-lchi2(j);
		ls=sqrt(pow(fs*la,2)+pow(lc,2));
		while (ls<fs*smooth_length && j<ind_alpha_vec-1)
		{
			j++;
			la=lalpha(j-1)-lalpha(j);
			lc=lchi2(j-1)-lchi2(j);
			ls=ls+sqrt(pow(fs*la,2)+pow(lc,2));
		}
		ind_curv_start=j;
		//cout<<"ind_curv_start: "<<ind_curv_start<<endl;
		
		j=ind_alpha_vec-2;
		la=lalpha(j)-lalpha(j+1);
		lc=lchi2(j)-lchi2(j+1);
		ls=sqrt(pow(fs*la,2)+pow(lc,2));
		while (ls<fs*smooth_length && j>0)
		{
			j--;
			la=lalpha(j)-lalpha(j+1);
			lc=lchi2(j)-lchi2(j+1);
			ls=ls+sqrt(pow(fs*la,2)+pow(lc,2));
		}
		ind_curv_end=j;
		//cout<<"ind_curv_end: "<<ind_curv_end<<endl;
		
		int Ncurv=ind_curv_end-ind_curv_start+1;
		//cout<<"Ncurv: "<<Ncurv<<endl;
		
		vec dlchi2_lalpha=zeros<vec>(Ncurv);
		vec curv_lchi2_lalpha=zeros<vec>(Ncurv);
		
		vec arc_params;
		int jmin, jmax;
		
		double xc, yc, x1, y1;
		
		for (j=ind_curv_start; j<=ind_curv_end; j++)
		{
			jmin=j-1;
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
			jmax=j+1;
			la=lalpha(jmax-1)-lalpha(jmax);
			lc=lchi2(jmax-1)-lchi2(jmax);
			ls=sqrt(pow(fs*la,2)+pow(lc,2));
			while (ls<fs*smooth_length && jmax<ind_alpha_vec-1)
			{
				jmax++;
				la=lalpha(jmax-1)-lalpha(jmax);
				lc=lchi2(jmax-1)-lchi2(jmax);
				ls=ls+sqrt(pow(fs*la,2)+pow(lc,2));
			}
			fit_circle_arc(fs*lalpha.rows(jmin,jmax), lchi2.rows(jmin,jmax), arc_params);
			curv_lchi2_lalpha(j-ind_curv_start)=1.0/arc_params(0);
			
			xc=arc_params(1);
			yc=arc_params(2);
			x1=fs*lalpha(j);
			y1=lchi2(j);
			
			dlchi2_lalpha(j-ind_curv_start)=(xc-x1)/(y1-yc);
		 
			//x1=fs*lalpha(j-1)-xc;
			//y1=lchi2(j-1)-yc;
			//x2=fs*lalpha(j)-xc;
			//y2=lchi2(j)-yc;
			//v1=sqrt(x1*x1+y1*y1);
			//v2=sqrt(x2*x2+y2*y2);
			//angle(j-1)=acos((x1*x2+y1*y2)/(v1*v2));
			//total_curv_lchi2_lalpha(j-1)=sum(sign(curv_lchi2_lalpha.rows(0,j-1)) % angle.rows(0,j-1));
		}
		
		dlchi2_lalpha_max=dlchi2_lalpha.max();
		dlchi2_lalpha_min=dlchi2_lalpha(Ncurv-1);
		
		uword N_av=5;
		double dlchi2_lalpha_min_av=sum(dlchi2_lalpha.rows(Ncurv-N_av,Ncurv-1))/N_av;
		
		
		int ind_min_curv=0;
		if (alpha_opt_max_in.size())
		{
			while (alpha_vec(ind_min_curv)>alpha_opt_max && ind_min_curv<ind_alpha_vec-1)
				ind_min_curv++;
		}
		
		ind_min_curv=ind_min_curv-ind_curv_start;
		if (ind_min_curv<0) ind_min_curv=0;
		int ind_max_curv=ind_alpha_vec-1;
		while (alpha_vec(ind_max_curv)<alpha_opt_min && ind_max_curv>0)	 ind_max_curv--;
		ind_max_curv=ind_max_curv-ind_curv_start;
		if (ind_max_curv>=Ncurv) ind_max_curv=Ncurv-1;
		if (ind_max_curv<ind_min_curv) ind_max_curv=ind_min_curv;
		
		int ind_alpha_opt, ind_alpha_opt_l, ind_alpha_opt_r;
		
		double max_curv;
		if (ind_min_curv<Ncurv && ind_max_curv<Ncurv) // find ind_alpha_opt
		{
			uword ind_alpha_opt_tmp;
			vec curv=curv_lchi2_lalpha.rows(ind_min_curv,ind_max_curv);
			max_curv=curv.max(ind_alpha_opt_tmp);
			ind_alpha_opt=ind_alpha_opt_tmp+ind_min_curv;
			
			j=ind_alpha_opt;
			while (curv_lchi2_lalpha(j)>max_curv/2 && j>0)  j--;
			ind_alpha_opt_r=j+ind_curv_start;
			j=ind_alpha_opt;
			while (curv_lchi2_lalpha(j)>max_curv/2 && j<Ncurv-1) j++;
			ind_alpha_opt_l=j+ind_curv_start;
			ind_alpha_opt=ind_alpha_opt+ind_curv_start;
		}
		else
		{
			ind_alpha_opt=ind_alpha_vec-1;
		}
		
		vec G_V_out;
		G_V_out=KG_V*A;
		
		vec M_out, M_V_out, errM_out;
		if (NM>0)
		{
			if (!covm_diag) errM=sqrt(COVM.diag());
			M_out=KM*A;
			errM_out=(M-M_out)/errM;
			M_V_out=KM_V*A;
		}
		
		vec errRe, errIm;
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
			errRe=G_V-G_V_out;
		
		vec eigv_ind;
		if (cov_diag)
			eigv_ind=wn;
		else
			eigv_ind=linspace<vec>(0,errRe.n_rows-1,errRe.n_rows);
		
		
		int l;
		vec CRe, CIm;
		if (!boson || col_Gi>0)
		{
			if (cov_diag)
			{
				CRe.zeros(NnC);
				CIm.zeros(NnC);
				for (j=0; j<NnC; j++)
				{
					for (l=0; l<Nn-j; l++)
					{
						CRe(j)=CRe(j)+errRe(l)*errRe(l+j);
						CIm(j)=CIm(j)+errIm(l)*errIm(l+j);
					}
					CRe(j)=CRe(j)/(Nn-j);
					CIm(j)=CIm(j)/(Nn-j);
				}
			}
			else
			{
				CRe.zeros(NnC);
				for (j=0; j<NnC; j++)
				{
					for (l=0; l<2*Nn-j; l++)
					{
						CRe(j)=CRe(j)+errRe(l)*errRe(l+j);
					}
					CRe(j)=CRe(j)/(2*Nn-j);
				}
			}
				
		}
		else
		{
			CRe.zeros(NnC);
			for (j=0; j<NnC; j++)
			{
				for (l=0; l<Nn-j; l++)
					CRe(j)=CRe(j)+errRe(l)*errRe(l+j);
				CRe(j)=CRe(j)/(Nn-j);
			}
		}
		
		
		vec Dn;
		Dn=linspace<vec>(0,NnC-1,NnC);

		//if ( (chi2_vec(ind_alpha_opt)<chi2save || alpha_vec(ind_alpha_opt)<alpha_save) && D2curv_max<0 && ind_alpha_opt_l<Ncurv)
		if (max_curv>0 && (ind_alpha_opt_l-ind_curv_start)<Ncurv-1)
		{
			double f_alpha_save=pow(10,save_alpha_range);
			//cout<<"f_alpha_save: "<<f_alpha_save<<endl;
			if (!alpha_save_max_in.size()) alpha_save_max=f_alpha_save*alpha_vec(ind_alpha_opt);
			if (!alpha_save_min_in.size()) alpha_save_min=alpha_vec(ind_alpha_opt)/f_alpha_save;
			
			cout<<"optimal value of alpha: "<<alpha_vec(ind_alpha_opt)<<endl;
			cout<<"chi2 at optimal alpha: "<<chi2_vec(ind_alpha_opt)<<endl;
			
			//cout << vectors_A.size() << " " << ind_alpha_opt << " " << ind_alpha_opt_r << " " << ind_alpha_opt_l << endl;
			
			vec A_opt, A_opt_l, A_opt_r;
			A_opt.reset();
			A_opt_l.reset();
			A_opt_r.reset();
			
			A_opt   = vectors_A[ind_alpha_opt];
			A_opt_r = vectors_A[ind_alpha_opt_r];
			A_opt_l = vectors_A[ind_alpha_opt_l];
			
			
			//cout << A_opt.size() << " " << w.size() << endl;
			FILE *fileOut = fopen("A.dat","w");
			fprintf(fileOut, "omega          A(omega)        A(omega)         A(omega)\n"); 
			int bosonOffset = 0;
			if (!boson || col_Gi>0) bosonOffset=1;
			for(int ii=bosonOffset; ii<A_opt.size(); ii++) fprintf(fileOut,"% 5.9f % 4.8f % 4.8f % 4.8f\n", w[ii], A_opt[ii-bosonOffset], A_opt_l[ii-bosonOffset], A_opt_r[ii-bosonOffset]); // I still don't know yet why, but the first and last index of all the "w" vector are to be trashed
			fclose(fileOut);


			
		}
		else if (alpha<=alpha_min)
		{
			cout<<"minimum value of alpha reached but optimal spectrum has not been found. Reducing alpha_min by a factor "<<f_alpha_min<<".\n";
			alpha_min=alpha_min/f_alpha_min;
			cout<<"new value of alpha_min: "<<alpha_min<<endl;
		}
			
		if (alpha<=alpha_min && dlchi2_lalpha_min_av/dlchi2_lalpha_max>RMAX_dlchi2_lalpha)
		{
			cout<<"The minimum value of alpha may not be small enough. Reducing alpha_min by a factor "<<f_alpha_min<<".\n";
			alpha_min=alpha_min/f_alpha_min;
			cout<<"new value of alpha_min: "<<alpha_min<<endl;
		}
	}

	return;
}










bool OmegaMaxEnt_data::load_data_file(mat &data_array, string file_name)
{
	string complete_file_name(input_dir);
	complete_file_name+=file_name;
	
	if (data_array.load(complete_file_name.c_str(),auto_detect))
		return true;
	else if (data_array.load(file_name.c_str(),auto_detect))
		return true;
	else
	{
		cout<<"file "<<file_name<<" could not be opened\n";
		return false;
	}
}


