

#include "maxEnt_data.h"

bool polyfit(vec x, vec y, int D, double x0, vec &cfs);
bool polyval(double x0, vec cfs, vec x, vec &y);
void remove_spaces_front(string &str);
void remove_spaces_back(string &str);
void remove_spaces_ends(string &str);
void pascal(int n, imat &P);


bool OmegaMaxEnt_data::set_G_omega_n_fermions()
{
	int j;
	
	signG=1;
	
	wn=green_data.col(0);
	Nn=green_data.n_rows;
	Gr=green_data.col(col_Gr-1);
	Gi=green_data.col(col_Gi-1);
	
	if (static_part_G_in.size()) Gr=Gr-static_part_G;
	
	indG_0=0, indG_f=Nn-1;
	wn_sign_change=false;
	wn_inversed=false;
	if (wn(0)*wn(Nn-1)<0)
	{
		wn_sign_change=true;
		j=1;
		while (wn(j)*wn(Nn-1)<0)
			j=j+1;
		if ((Nn-j)>=j+1)
		{
			indG_0=j;
			indG_f=Nn-1;
		}
		else
		{
			indG_0=0;
			indG_f=j-1;
		}
		wn=wn.rows(indG_0,indG_f);
		Nn=wn.n_rows;
		Gr=Gr.rows(indG_0,indG_f);
		Gi=Gi.rows(indG_0,indG_f);
	}
	if (wn(0)<0)
	{
		wn=abs(wn);
		Gi=-Gi;
	}
	if (wn(0)>wn(Nn-1))
	{
		wn_inversed=true;
		wn=flipud(wn);
		Gr=flipud(Gr);
		Gi=flipud(Gi);
	}
	if (Gi.max()*Gi.min()<0)
	{
		cout<<"error: Im[G] must not change sign\n";
		return false;
	}
	if (Gi.max()>0)
	{
		signG=-1;
	}
	Gr=signG*Gr;
	Gi=signG*Gi;
	
	G.zeros(Nn);
	G.set_real(Gr);
	G.set_imag(Gi);
	
	for (j=0; j<Nn-1; j++)
	{
		if ((wn(j+1)-wn(j))<0)
		{
			cout<<"error: Matsubara frequency is not strictly increasing\n";
			return false;
		}
	}
	
	uvec even_ind=linspace<uvec>(0,2*Nn-2,Nn);
	uvec odd_ind=even_ind+1;
	
	Gchi2.zeros(2*Nn);
	Gchi2.rows(even_ind)=Gr;
	Gchi2.rows(odd_ind)=Gi;
	
	cout<<"Number of Matsubara frequencies in the Green function: "<<Nn<<endl;
	
	double 	wn_min=wn.min();
	double tem_wn=wn_min/PI;
	
	n.zeros(Nn);
	vec ntmp;
	if (tem_in.size()==0)
	{
		if ( (wn(1)-floor(wn(1)))==0 && (wn(2)-floor(wn(2)))==0 )
		{
			cout<<"If the Matsubara frequency are given by index, temperature must be provided.\n";
			return false;
		}
		tem=tem_wn;
		ntmp=(round(wn/(PI*tem))-1)/2;
		for (j=0; j<Nn; j++)	n(j)=ntmp(j);
	}
	else if ( (abs(tem-tem_wn)/tem)>tol_tem )
	{
		if ( wn(1)-floor(wn(1)) )
		{
			ntmp=(round(wn/(PI*tem_wn))-1)/2;
			for (j=0; j<Nn; j++)	n(j)=ntmp(j);
			cout<<"warning: temperature in file "<<data_file_name<<" is different from temperature given in parameter file "<<input_params_file_name<<endl;
			cout<<"provided temperature: "<<tem<<endl;
			cout<<"temperature extracted from first Matsubara frequency: "<<tem_wn<<endl;
		}
		else if ( (wn(2)-floor(wn(2)))==0 )
		{
			ntmp=round(wn);
			for (j=0; j<Nn; j++) n(j)=ntmp(j);
			n=n-n(0);
		}
		else
		{
			cout<<"oops! First Matsubara frequency is an integer but second is not!\n";
			n=linspace<uvec>(0,Nn-1,Nn);
		}
	}
	else
	{
		ntmp=(round(wn/(PI*tem_wn))-1)/2;
		for (j=0; j<Nn; j++) n(j)=ntmp(j);
	}

	
	wn=PI*tem*conv_to<vec>::from(2*n+1);
	
	return true;
}



bool OmegaMaxEnt_data::set_covar_G_omega_n()
{
	int j;
	
	COV.zeros(2*Nn,2*Nn);
	cov_diag=true;
	if (error_file.size())
	{
		if (error_data.n_rows<Nn)
		{
			cout<<"number of lines is too small in file "<<error_file<<endl;
			return false;
		}
		cout<<"error file provided\n";
		errGr=error_data.col(col_errGr-1);
		errGi=error_data.col(col_errGi-1);
		if (wn_sign_change)
		{
			errGr=errGr.rows(indG_0,indG_f);
			errGi=errGi.rows(indG_0,indG_f);
		}
		if (wn_inversed)
		{
			errGr=flipud(errGr);
			errGi=flipud(errGi);
		}
		errG.zeros(2*Nn);
		for (j=0; j<Nn; j++)
		{
			errG(2*j)=errGr(j);
			errG(2*j+1)=errGi(j);
		}
		COV.diag()=square(errG);
	}
	else if ( covar_re_re_file.size() && covar_im_im_file.size() && covar_re_im_file.size() )
	{
		if (CRR.n_rows<Nn || CRR.n_cols<Nn || CII.n_rows<Nn || CII.n_cols<Nn || CRI.n_rows<Nn || CRI.n_cols<Nn)
		{
			cout<<"number of lines and/or columns in covariance file(s) is too small\n";
			return false;
		}
		cout<<"covariance matrix provided\n";
		cov_diag=false;
		uvec even_ind=linspace<uvec>(0,2*Nn-2,Nn);
		uvec odd_ind=even_ind+1;
		COV(even_ind,even_ind)=CRR.submat(0,0,Nn-1,Nn-1);
		COV(odd_ind,odd_ind)=CII.submat(0,0,Nn-1,Nn-1);
		COV(even_ind,odd_ind)=CRI.submat(0,0,Nn-1,Nn-1);
		COV(odd_ind,even_ind)=trans(CRI.submat(0,0,Nn-1,Nn-1));
	}
	else
	{
		cout<<"no errors provided\nusing a constant error\n";
		if (!boson)
		{
			double Gi_max=max(abs(Gi));
			errGi=default_error_G*Gi_max*ones<vec>(Nn);
			errGr=errGi;
			errG=default_error_G*Gi_max*ones<vec>(2*Nn);
		}
		else
		{
			double Gr_max=max(abs(Gr));
			errGr=default_error_G*Gr_max*ones<vec>(Nn);
			errGi=errGr;
			errG=default_error_G*Gr_max*ones<vec>(2*Nn);
		}
		COV.diag()=square(errG);
	}
	
	COV=0.5*(COV+COV.t());
	
	return true;
}







bool OmegaMaxEnt_data::set_moments_fermions()
{
	NM_odd=0;
	
	moments_provided=false;
	covm_diag=true;
	
	M0=1.0;
	if (tau_GF)
		M0=M0t;
	if (M0_in.size())
		M0=stod(M0_in);
	errM0=err_norm*M0;
	
	M.zeros(1);
	M(0)=M0;
	errM.zeros(1);
	errM(0)=errM0;
	NM=1;
	M_ord=linspace<vec>(0,NM-1,NM);
	
	M_even.zeros(1);
	M_even(0)=M0;
	NM_even=1;
	
	std_omega=0;
	
	if (M1_in.size())
	{
		M1=stod(M1_in);
		M1_set=true;
//		NMinput=2;
		moments_provided=true;
		cout<<"moments provided\n";
		M.zeros(2);
		M(0)=M0;
		M(1)=M1;
		NM=2;
		M_ord=linspace<vec>(0,NM-1,NM);
		
		M_odd.zeros(1);
		M_odd(0)=M1;
		NM_odd=1;
		if (!SC_set)
		{
			SC=M1/M0;
			SC_set=true;
		}
		if (M2_in.size())
		{
			M2=stod(M2_in);
			M2_set=true;
			if (M2<=0)
			{
				cout<<"Error: second moment of spectral function must be greater than 0.\n";
				return false;
			}
//			NMinput=3;
			M.zeros(3);
			M(0)=M0;
			M(1)=M1;
			M(2)=M2;
			NM=3;
			M_ord=linspace<vec>(0,NM-1,NM);
			double var_omega=M2/M0-pow(M1/M0,2);
			if (var_omega>0)
			{
				std_omega=sqrt(var_omega);
				if (!SW_set)
				{
					SW=f_SW_std_omega*std_omega;
					SW_set=true;
				}
			}
			else
			{
				cout<<"Error: negative variance. First and/or second moment incorect.\n";
				return false;
			}
			errM2=default_error_M*M2;
			if (errM2_in.size())
				errM2=stod(errM2_in);
		}
		if (M3_in.size())
		{
			M3=stod(M3_in);
			if (M2_in.size())
			{
				M.zeros(4);
				M(0)=M0;
				M(1)=M1;
				M(2)=M2;
				M(3)=M3;
				NM=4;
				M_ord=linspace<vec>(0,NM-1,NM);
			}
			M_odd.zeros(2);
			M_odd(0)=M1;
			M_odd(1)=M3;
			NM_odd=2;
			if (errM3_in.size())
				errM3=stod(errM3_in);
			else if (M2_in.size())
				errM3=default_error_M*pow(std_omega,3);
			else if (SW_in.size())
				errM3=default_error_M*pow(SW,3);
			else
				errM3=default_error_M;
		}
		NM=M.n_rows;
		if (errM1_in.size())
			errM1=stod(errM1_in);
		else if (M2_in.size())
			errM1=default_error_M*std_omega;
		else if (SW_in.size())
			errM1=default_error_M*SW;
		else
			errM1=default_error_M;
		errM.zeros(NM);
		errM(0)=errM0;
		errM(1)=errM1;
		if (M2_in.size())
		{
			errM(2)=errM2;
			if (M3_in.size())
				errM(3)=errM3;
		}
	}
	COVM.zeros(NM,NM);
	COVM.diag()=square(errM);
	
	if ( NM<3 && !SW_in.size() && (!use_grid_params || omega_grid_params.n_cols<3) && !grid_omega_file.size() && !init_spectr_func_file.size() )
	{
		if (!tau_GF)
			cout<<"Not enough moments provided to define the real frequency grid. The program will try to extract moments from the high frequencies of the Green function\n";
		else
			cout<<"Not enough moments provided to define the real frequency grid. The program will try to extract moments from G(tau) around tau=0 and tau=beta.\n";
		eval_moments=true;
	}
	
	if (M2_in.size())
	{
		if (!M1_in.size())
		{
			M2=stod(M2_in);
			M2_set=true;
			if (M2<=0)
			{
				cout<<"Error: second moment of spectral function must be greater than 0.\n";
				return false;
			}
			errM2=default_error_M*M2;
			if (errM2_in.size())
				errM2=stod(errM2_in);
		}
		M_even.zeros(2);
		M_even(0)=M0;
		M_even(1)=M2;
		NM_even=2;
	}
	
	NMinput=NM_even+NM_odd;
	
	return true;
}















bool OmegaMaxEnt_data::set_moments_bosons()
{
	if (!tau_GF) M1n=abs(Gr(0));
	moments_provided=false;
	covm_diag=true;
	std_omega=0;
	NM=0;
	
	if (!tau_GF)
	{
		if (max(abs(Gi))!=0)
		{
			M0=1.0;
			if (M0_in.size())
				M0=stod(M0_in);
			errM0=err_norm*M0;
		}
		else
		{
			M0=0;
			errM0=default_error_M;
			M2=0;
			errM2=default_error_M;
		}
	}
	else
	{
		M0=M0t;
		if (M0_in.size())
			M0=stod(M0_in);
		if (M0)
			errM0=err_norm*M0;
		else
			errM0=default_error_M;
	}
	
	if (!tau_GF)
	{
		if (!SC_set)
		{
			SC=M0/M1n;
			SC_set=true;
		}
	}
	
	if (col_Gi>0)
	{
		M.zeros(1);
		M(0)=M0;
		errM.zeros(1);
		errM(0)=errM0;
		NMinput=1;
		NM=1;
		M_ord=linspace<vec>(0,NM-1,NM);
	}

	if (M1_in.size())
	{
		M1=stod(M1_in);
		M1_set=true;
		if (M1<=0)
		{
			cout<<"error: first moment of spectral function must be greater than 0 for bosons\n";
			return false;
		}
		moments_provided=true;
		cout<<"moments provided\n";
		if (col_Gi>0)
		{
			M.zeros(2);
			M(0)=M0;
			M(1)=M1;
			NM=2;
			M_ord=linspace<vec>(0,NM-1,NM);
		}
		else
		{
			M.zeros(1);
			M(0)=M1;
			NM=1;
			M_ord.zeros(1);
			M_ord(0)=1;
		}
		
		double var_omega=M1/M1n-pow(M0/M1n,2);
		if (var_omega<0)
		{
			cout<<"error: the variance computed from the provided moments is negative\n";
			return false;
		}
		std_omega=sqrt(var_omega);
		if (!SW_set)
		{
			SW=f_SW_std_omega*std_omega;
			SW_set=true;
		}
		if (errM1_in.size())
			errM1=stod(errM1_in);
		else
			errM1=default_error_M*M1;
		if (col_Gi>0)
		{
			if (M2_in.size())
			{
				M2=stod(M2_in);
				M2_set=true;
				M.zeros(3);
				M(0)=M0;
				M(1)=M1;
				M(2)=M2;
				NM=3;
				M_ord=linspace<vec>(0,NM-1,NM);
				if (errM2_in.size())
					errM2=stod(errM2_in);
				else
					errM2=default_error_M*M1n*pow(std_omega,3);
			}
		}
		errM.zeros(NM);
		if (col_Gi>0)
		{
			errM(0)=errM0;
			errM(1)=errM1;
			if (M2_in.size())	errM(2)=errM2;
		}
		else
			errM(0)=errM1;
	}
	COVM=diagmat(square(errM));
	
	if (col_Gi>0)
	{
		if ( NM<2 && !SW_in.size() && (!use_grid_params || omega_grid_params.n_cols<3) && !grid_omega_file.size() && !init_spectr_func_file.size() )
		{
			if (!tau_GF)
				cout<<"note: not enough information provided to define the real frequency grid. The program will try to extract the moments from the high frequencies of the Green function\n";
			else
				cout<<"Not enough moments provided to define the real frequency grid. The program will try to extract moments from G(tau) around tau=0 and tau=beta.\n";
			eval_moments=true;
		}
		
	}
	else
	{
		if ( NM<1 && !SW_in.size() && (!use_grid_params || omega_grid_params.n_cols<3) && !grid_omega_file.size() && !init_spectr_func_file.size() )
		{
			if (!tau_GF)
				cout<<"note: not enough information provided to define the real frequency grid. The program will try to extract the moments from the high frequencies of the Green function\n";
			else
				cout<<"Not enough moments provided to define the real frequency grid. The program will try to extract moments from G(tau) around tau=0 and tau=beta.\n";
			eval_moments=true;
		}
	}
	
	return true;
}



bool OmegaMaxEnt_data::set_covar_Gtau()
{
	cov_diag=false;
	
	if (error_file.size())
	{
		cout<<"error file provided\n";
		if (error_data.n_rows<Ntau)
		{
			cout<<"set_covar_Gtau() error: number of lines in error file is smaller than number of imaginary time steps\n";
			return false;
		}
		errGtau=error_data.col(col_errGtau-1);
		errGtau.resize(Ntau);
		Ctau.zeros(Ntau,Ntau);
		Ctau.diag()=square(errGtau.rows(0,Ntau-1));
		Ctau_all.zeros(Ntau+1,Ntau+1);
		Ctau_all.submat(0,0,Ntau-1,Ntau-1)=Ctau;
		Ctau_all(Ntau,Ntau)=pow(errGtau(0),2);
		
	}
	else if (covar_tau_file.size())
	{
		if (Ctau.n_rows<Ntau)
		{
			cout<<"set_covar_Gtau() error: number of lines in covariance file is smaller than number of imaginary time steps\n";
			return false;
		}
		cout<<"imaginary time covariance matrix provided\n";
		Ctau=Ctau.submat(0,0,Ntau-1,Ntau-1);
		Ctau_all.zeros(Ntau+1,Ntau+1);
		Ctau_all.submat(0,0,Ntau-1,Ntau-1)=Ctau;
		errGtau.reset();
	}
	else
	{
		cout<<"no errors provided\nusing a constant error\n";
		double maxGtau=max(abs(Gtau));
	 	errGtau=default_error_G*maxGtau*ones<vec>(Ntau);
		Ctau.zeros(Ntau,Ntau);
		Ctau.diag()=square(errGtau.rows(0,Ntau-1));
		Ctau_all.zeros(Ntau+1,Ntau+1);
		Ctau_all.submat(0,0,Ntau-1,Ntau-1)=Ctau;
		Ctau_all(Ntau,Ntau)=pow(errGtau(0),2);
	}
	
	COV=Ctau;
	
	Ctau_all.submat(Ntau,0,Ntau,Ntau-1)=-Ctau.row(0);
	Ctau_all.submat(0,Ntau,Ntau-1,Ntau)=-Ctau.col(0);
	Ctau_all(Ntau,Ntau)=Ctau(0,0);
	
	COV=(COV+COV.t())/2.0;
	
	return true;
}






bool OmegaMaxEnt_data::compute_moments_tau_bosons()
{
	cout<<"COMPUTING MOMENTS\n";
	
	int Nv=1;
	int NvN=1;
	int npmin=2;
	int npmax=5;
	int Np=npmax-npmin+1;
	int DNfitmin=0;
	int DNfitmax=Ntau-npmax-1;
	if (16-npmax-1<DNfitmax) DNfitmax=16-npmax-1;
	int NDN=DNfitmax-DNfitmin+1;
	
	mat M0tmp=zeros<mat>(NDN,Np);
	mat M1tmp=zeros<mat>(NDN,Np);
	mat M2tmp=zeros<mat>(NDN,Np);
	
	vec Gtmp, pp;
	int DNfit, np, Nfit;
	for (DNfit=DNfitmin; DNfit<=DNfitmax; DNfit++)
	{
		for (np=npmin; np<=npmax; np++)
		{
			Nfit=np+1+DNfit;
			Gtmp=Gtau.rows(0,Nfit-1)+flipud(Gtau.rows(Ntau-Nfit+1,Ntau));
			polyfit(tau.rows(0,Nfit-1),Gtmp,np,0,pp);
			M1tmp(DNfit-DNfitmin,np-npmin)=pp(np-1);
			
			Gtmp=Gtau.rows(0,Nfit-1)-flipud(Gtau.rows(Ntau-Nfit+1,Ntau));
			polyfit(tau.rows(0,Nfit-1),Gtmp,np,0,pp);
			M0tmp(DNfit-DNfitmin,np-npmin)=-pp(np);
			M2tmp(DNfit-DNfitmin,np-npmin)=-2*pp(np-2);
		}
	}
	
	mat M0m=zeros<mat>(NDN-2*NvN,Np-2*Nv);
	mat M1m=zeros<mat>(NDN-2*NvN,Np-2*Nv);
	mat M2m=zeros<mat>(NDN-2*NvN,Np-2*Nv);
	mat varM0=zeros<mat>(NDN-2*NvN,Np-2*Nv);
	mat varM1=zeros<mat>(NDN-2*NvN,Np-2*Nv);
	mat varM2=zeros<mat>(NDN-2*NvN,Np-2*Nv);
	int j, l;
	for (l=Nv; l<Np-Nv; l++)
	{
		for (j=NvN; j<NDN-NvN; j++)
		{
			M0m(j-NvN,l-Nv)=accu(M0tmp.submat(j-NvN,l-Nv,j+NvN,l+Nv))/((2*Nv+1)*(2*NvN+1));
			M1m(j-NvN,l-Nv)=accu(M1tmp.submat(j-NvN,l-Nv,j+NvN,l+Nv))/((2*Nv+1)*(2*NvN+1));
			M2m(j-NvN,l-Nv)=accu(M2tmp.submat(j-NvN,l-Nv,j+NvN,l+Nv))/((2*Nv+1)*(2*NvN+1));
			varM0(j-NvN,l-Nv)=accu(pow(M0tmp.submat(j-NvN,l-Nv,j+NvN,l+Nv)-M0m(j-NvN,l-Nv),2))/((2*Nv+1)*(2*NvN+1));
			varM1(j-NvN,l-Nv)=accu(pow(M1tmp.submat(j-NvN,l-Nv,j+NvN,l+Nv)-M1m(j-NvN,l-Nv),2))/((2*Nv+1)*(2*NvN+1));
			varM2(j-NvN,l-Nv)=accu(pow(M2tmp.submat(j-NvN,l-Nv,j+NvN,l+Nv)-M2m(j-NvN,l-Nv),2))/((2*Nv+1)*(2*NvN+1));
		}
	}
	
	uword jvmin, lvmin;
	varM0.min(jvmin,lvmin);
	double M0_NP_tmp=M0m(jvmin,lvmin);
	varM1.min(jvmin,lvmin);
	double M1_NP_tmp=M1m(jvmin,lvmin);
	varM2.min(jvmin,lvmin);
	//double M2_NP_tmp=M2m(jvmin,lvmin);
	
	double Wtmp;
	if (M1_NP_tmp/M1n>pow(M0_NP_tmp/M1n,2))
		Wtmp=sqrt(M1_NP_tmp/M1n-pow(M0_NP_tmp/M1n,2));
	else
		Wtmp=abs(M0_NP_tmp/M1n);
	
	int Nfitmax=ceil(FNfitTauW*Ntau*tem/(abs(M0_NP_tmp/M1n)+Wtmp));
	if (Nfitmax>Ntau) Nfitmax=Ntau;
	
	mat X, CG, invCG, AM;
	vec Gchi2tmp, BM, Mtmp;
	npmin=3;
	int Nfitmin=npmin+1;
	int NNfit=Nfitmax-Nfitmin+1;
	mat M0b=zeros<mat>(NNfit,NNfit);
	mat M1b=zeros<mat>(NNfit,NNfit);
	mat M2b=zeros<mat>(NNfit,NNfit);
	mat M3b=zeros<mat>(NNfit,NNfit);
	int p;
	for (Nfit=Nfitmin; Nfit<=Nfitmax; Nfit++)
	{
		for (np=npmin; np<Nfit; np++)
		{
			X=zeros<mat>(Nfit,np+1);
			for (p=0; p<=np; p++)
				X.col(p)=pow(tau.rows(0,Nfit-1),p);
			
			Gchi2tmp=Gtau.rows(0,Nfit-1)+flipud(Gtau.rows(Ntau-Nfit+1,Ntau));
			CG=Ctau_all.submat(0,0,Nfit-1,Nfit-1)+fliplr(Ctau_all.submat(0,Ntau-Nfit+1,Nfit-1,Ntau))+flipud(Ctau_all.submat(Ntau-Nfit+1,0,Ntau,Nfit-1))+flipud(fliplr(Ctau_all.submat(Ntau-Nfit+1,Ntau-Nfit+1,Ntau,Ntau)));
			CG(0,0)=CG(1,1);
			invCG=inv(CG);
//			invCG=inv_sympd(CG);
			AM=(X.t())*invCG*X;
			BM=(X.t())*invCG*Gchi2tmp;
			Mtmp=solve(AM,BM);
			M1b(Nfit-np-1,np-npmin)=Mtmp(1);
			M3b(Nfit-np-1,np-npmin)=6*Mtmp(3);
			
			Gchi2tmp=Gtau.rows(0,Nfit-1)-flipud(Gtau.rows(Ntau-Nfit+1,Ntau));
			CG=Ctau_all.submat(0,0,Nfit-1,Nfit-1)- fliplr(Ctau_all.submat(0,Ntau-Nfit+1,Nfit-1,Ntau))-flipud(Ctau_all.submat(Ntau-Nfit+1,0,Ntau,Nfit-1))+flipud(fliplr(Ctau_all.submat(Ntau-Nfit+1,Ntau-Nfit+1,Ntau,Ntau)));
			invCG=inv(CG);
//			invCG=inv_sympd(CG);
			AM=(X.t())*invCG*X;
			BM=(X.t())*invCG*Gchi2tmp;
			Mtmp=solve(AM,BM);
			M0b(Nfit-np-1,np-npmin)=-Mtmp(0);
			M2b(Nfit-np-1,np-npmin)=-2*Mtmp(2);
			
		}
	}
	
	
	Nv=1;
	NvN=1;
	int jmin=NvN;
	int jmax=NNfit-NvN-2*Nv-1;
	int Nj=jmax-jmin+1;
	int lmin=Nv;
	int lmax=NNfit-2*NvN-Nv-1;
	int Nl=lmax-lmin+1;
	M0m.zeros(Nj,Nl);
	M1m.zeros(Nj,Nl);
	M2m.zeros(Nj,Nl);
	mat M3m=zeros<mat>(Nj,Nl);
	varM0.zeros(Nj,Nl);
	varM1.zeros(Nj,Nl);
	varM2.zeros(Nj,Nl);
	mat varM3=zeros<mat>(Nj,Nl);
	for (j=jmin; j<=jmax; j++)
	{
		for (l=lmin; l<=NNfit-1-j-Nv-NvN; l++)
		{
			M0m(j-jmin,l-lmin)=accu(M0b.submat(j-NvN,l-Nv,j+NvN,l+Nv))/((2*Nv+1)*(2*NvN+1));
			M1m(j-jmin,l-lmin)=accu(M1b.submat(j-NvN,l-Nv,j+NvN,l+Nv))/((2*Nv+1)*(2*NvN+1));
			M2m(j-jmin,l-lmin)=accu(M2b.submat(j-NvN,l-Nv,j+NvN,l+Nv))/((2*Nv+1)*(2*NvN+1));
			M3m(j-jmin,l-lmin)=accu(M3b.submat(j-NvN,l-Nv,j+NvN,l+Nv))/((2*Nv+1)*(2*NvN+1));
			varM0(j-jmin,l-lmin)=accu(pow(M0b.submat(j-NvN,l-Nv,j+NvN,l+Nv)-M0m(j-jmin,l-lmin),2))/((2*Nv+1)*(2*NvN+1));
			varM1(j-jmin,l-lmin)=accu(pow(M1b.submat(j-NvN,l-Nv,j+NvN,l+Nv)-M1m(j-jmin,l-lmin),2))/((2*Nv+1)*(2*NvN+1));
			varM2(j-jmin,l-lmin)=accu(pow(M2b.submat(j-NvN,l-Nv,j+NvN,l+Nv)-M2m(j-jmin,l-lmin),2))/((2*Nv+1)*(2*NvN+1));
			varM3(j-jmin,l-lmin)=accu(pow(M3b.submat(j-NvN,l-Nv,j+NvN,l+Nv)-M3m(j-jmin,l-lmin),2))/((2*Nv+1)*(2*NvN+1));
		}
	}
	
	//	cout<<varM0<<endl;
	
	double varM0max=max(max(varM0));
	double varM1max=max(max(varM1));
	double varM2max=max(max(varM2));
	double varM3max=max(max(varM3));
	
	for (j=1; j<Nj; j++)
	{
		varM0.submat(j,Nl-j,j,Nl-1)=2*varM0max*ones<rowvec>(j);
		varM1.submat(j,Nl-j,j,Nl-1)=2*varM1max*ones<rowvec>(j);
		varM2.submat(j,Nl-j,j,Nl-1)=2*varM2max*ones<rowvec>(j);
		varM3.submat(j,Nl-j,j,Nl-1)=2*varM3max*ones<rowvec>(j);
	}
	
	varM0.min(jvmin,lvmin);
	double M0_N=M0m(jvmin,lvmin);
	varM1.min(jvmin,lvmin);
	double M1_N=M1m(jvmin,lvmin);
	varM2.min(jvmin,lvmin);
	double M2_N=M2m(jvmin,lvmin);
	varM3.min(jvmin,lvmin);
	double M3_N=M3m(jvmin,lvmin);
	
	cout<<"moments determined by polynomial fit to G(tau) at boundaries:\n";
	cout<<"norm: "<<M0_N<<endl;
	cout<<"first moment: "<<M1_N<<endl;
	cout<<"second moment: "<<M2_N<<endl;
	cout<<"third moment: "<<M3_N<<endl;

	varM1.min(jvmin,lvmin);
	np=lvmin+Nv+2;
	Nfit=jvmin+NvN+np+1;
	X=zeros<mat>(Nfit,np+1);
	for (p=0; p<=np; p++)
		X.col(p)=pow(tau.rows(0,Nfit-1),p);
	
	CG=Ctau_all.submat(0,0,Nfit-1,Nfit-1)+fliplr(Ctau_all.submat(0,Ntau-Nfit+1,Nfit-1,Ntau))+flipud(Ctau_all.submat(Ntau-Nfit+1,0,Ntau,Nfit-1))+flipud(fliplr(Ctau_all.submat(Ntau-Nfit+1,Ntau-Nfit+1,Ntau,Ntau)));
	CG(0,0)=CG(1,1);
	invCG=inv(CG);
//	invCG=inv_sympd(CG);
	AM=(X.t())*invCG*X;
	mat invAMp=inv(AM);
	
	CG=Ctau_all.submat(0,0,Nfit-1,Nfit-1)-fliplr(Ctau_all.submat(0,Ntau-Nfit+1,Nfit-1,Ntau))-flipud(Ctau_all.submat(Ntau-Nfit+1,0,Ntau,Nfit-1))+flipud(fliplr(Ctau_all.submat(Ntau-Nfit+1,Ntau-Nfit+1,Ntau,Ntau)));
	invCG=inv(CG);
	AM=(X.t())*invCG*X;
	mat invAMn=inv(AM);
	
	double std_omega_tmp;
	double var_omega=M1_N/M1n-pow(M0_N/M1n,2);
	if (var_omega>0)
		std_omega_tmp=sqrt(var_omega);
	else
	{
		cout<<"Negative variance found during computation of moments.\n";
		return false;
	}
	
	covm_diag=true;
	if (!moments_provided)
	{
		if (col_Gi>0)
		{
			M.zeros(3);
			NM=3;
			M(0)=M0;
			M(1)=M1_N;
			M(2)=M2_N;
			M1=M(1);
			M2=M(2);
			M1_set=true;
			M2_set=true;
			errM.zeros(NM);
			errM(0)=errM0;
			errM(1)=sqrt(invAMp(1,1));
			errM(2)=sqrt(invAMn(2,2));
			errM1=errM(1);
			errM2=errM(2);
			M_ord=linspace<vec>(0,NM-1,NM);
		}
		else
		{
			M.zeros(1);
			NM=1;
			M1=M1_N;
			M(0)=M1;
			errM.zeros(NM);
			errM1=sqrt(invAMp(1,1));
			errM(0)=errM1;
			M_ord.zeros(1);
			M_ord(0)=1;
		}
	}
	else
	{
		if (abs(M1-M1_N)/M1_N>tol_M1)
			cout<<"warning: first moment different from provided one\n";
		
		if (col_Gi>0)
		{
			if (M2_in.size())
			{
				if (abs(M2-M2_N)/pow(std_omega_tmp,2)>tol_M2)
					cout<<"warning: second moment different from provided one\n";
			}
			else
			{
				M2=M2_N;
				M2_set=true;
				errM2=sqrt(invAMn(2,2));
				M.zeros(3);
				NM=3;
				M(0)=M0;
				M(1)=M1;
				M(2)=M2;
				errM.zeros(NM);
				errM(0)=errM0;
				errM(1)=errM1;
				errM(2)=errM2;
				M_ord=linspace<vec>(0,NM-1,NM);
			}
		}
	}
	COVM.zeros(NM,NM);
	COVM.diag()=square(errM);
	
	return true;
}





bool OmegaMaxEnt_data::Fourier_transform_G_tau()
{
	int j;
	dcomplex I(0,1);
	
	cout<<"computing Fourier transform of G(tau)...\n";
	
	Nn=Ntau/2;
	n=linspace<uvec>(0,Nn-1,Nn);
	if (!boson)
		wn=PI*tem*conv_to<vec>::from(2*n+1);
	else
		wn=2*PI*tem*conv_to<vec>::from(n);
	
	vec cfs_LC(4);
	vec LC(2);
	if (!boson)
	{
		LC(0)=M1;
		LC(1)=-M2;
		cfs_LC=ones<vec>(4);
	}
	else
	{
		LC(0)=M1;
		LC(1)=-M2;
		cfs_LC(0)=1;
		cfs_LC(1)=-1;
		cfs_LC(2)=1;
		cfs_LC(3)=-1;
	}
	
	uvec l=linspace<uvec>(0,Ntau-1,Ntau);
	vec p=linspace<vec>(0,Ntau-1,Ntau);
	vec cfs_Gtau;
	spline_coeffs_LC(tau,Gtau,LC,cfs_LC,cfs_Gtau);
	
	G.zeros(Nn);
	Gr.zeros(Nn);
	Gi.zeros(Nn);
	if (boson)
	{
		double a, b, c, d, dtau;
		for (j=0; j<Ntau; j++)
		{
			a=cfs_Gtau(4*j);
			b=cfs_Gtau(4*j+1);
			c=cfs_Gtau(4*j+2);
			d=Gtau(j);
			dtau=tau(j+1)-tau(j);
			Gr(0)=Gr(0)+a*pow(dtau,4)/4+b*pow(dtau,3)/3+c*pow(dtau,2)/2+d*dtau;
		}
	}
	
	cx_vec TF_d3Gtau;
	cx_vec d3Gtau=zeros<cx_vec>(Ntau);
	d3Gtau.set_real(6*cfs_Gtau.rows(4*l));
	if (!boson)
	{
		d3Gtau=exp(I*PI*p/Ntau) % d3Gtau;
		TF_d3Gtau=Ntau*(1.0-exp(I*(2*p+1)*PI/Ntau)) % ifft(d3Gtau);
		G=-I*M0/wn-M1/pow(wn,2)+I*M2/pow(wn,3)+TF_d3Gtau.rows(0,Nn-1)/pow(wn,4);
	}
	else
	{
		TF_d3Gtau=Ntau*(1.0-exp(I*(2*p)*PI/Ntau)) % ifft(d3Gtau);
		G.rows(1,Nn-1)=-I*M0/wn.rows(1,Nn-1)-M1/pow(wn.rows(1,Nn-1),2)+I*M2/pow(wn.rows(1,Nn-1),3)+TF_d3Gtau.rows(1,Nn-1)/pow(wn.rows(1,Nn-1),4);
		G(0)=cx_double(Gr(0),0);
	}
	Gr=real(G);
	Gi=imag(G);
	
	uvec even_ind=linspace<uvec>(0,2*Nn-2,Nn);
	uvec odd_ind=even_ind+1;
	
	if (col_Gi>0)
	{
		Gchi2.zeros(2*Nn);
		Gchi2.rows(even_ind)=Gr;
		Gchi2.rows(odd_ind)=Gi;
	}
	else
	{
		Gchi2=Gr;
	}
	
//	clock_t tc=clock();
	
	cout<<setiosflags(ios::left);
	if (!boson)
	{
		cx_mat Mph=diagmat(exp(I*PI*p/Ntau));
		cx_mat Ctau2=Ctau*Mph;
		cx_mat Ctau_n(Ntau,Ntau);
	
		for (j=0; j<Ntau; j++) Ctau_n.row(j)=ifft(Ctau2.row(j))/tem;
		cx_mat Ctau_n_R=Mph*real(Ctau_n);
		cx_mat Ctau_n_I=Mph*imag(Ctau_n);
		cx_mat Cmn_R(Ntau,Ntau), Cmn_I(Ntau,Ntau);
		for (j=0; j<Nn; j++) Cmn_R.col(j)=ifft(Ctau_n_R.col(j))/tem;
		for (j=0; j<Nn; j++) Cmn_I.col(j)=ifft(Ctau_n_I.col(j))/tem;
		CRR=real(Cmn_R);
		CRI=real(Cmn_I);
		CII=imag(Cmn_I);
	}
	else
	{
		cx_mat Ctau_n(Ntau,Ntau), Ctau_n_R(Ntau,Ntau), Ctau_n_I(Ntau,Ntau), Cmn_R(Ntau,Ntau), Cmn_I(Ntau,Ntau);
	
		for (j=0; j<Ntau; j++) Ctau_n.row(j)=fft(Ctau.row(j))/(Ntau*tem);
		Ctau_n_R.zeros();
		Ctau_n_I.zeros();
		Ctau_n_R.set_real(real(Ctau_n));
		Ctau_n_I.set_real(imag(Ctau_n));
		for (j=0; j<Nn; j++) Cmn_R.col(j)=ifft(Ctau_n_R.col(j))/tem;
		for (j=0; j<Nn; j++) Cmn_I.col(j)=-ifft(Ctau_n_I.col(j))/tem;
		CRR=real(Cmn_R);
		CRI=real(Cmn_I);
		CII=imag(Cmn_I);
		CII(0,0)=CII(1,1);
	}
	
	if (col_Gi>0)
	{
		COV.zeros(2*Nn,2*Nn);
		COV(even_ind,even_ind)=CRR.submat(0,0,Nn-1,Nn-1);
		COV(odd_ind,odd_ind)=CII.submat(0,0,Nn-1,Nn-1);
		COV(even_ind,odd_ind)=CRI.submat(0,0,Nn-1,Nn-1);
		COV(odd_ind,even_ind)=trans(CRI.submat(0,0,Nn-1,Nn-1));
	}
	else
	{
		COV=CRR.submat(0,0,Nn-1,Nn-1);
	}
	
	COV=0.5*(COV+COV.t());
	
	cout<<"Fourier transform computed.\n";
	
	struct stat file_stat;
	
	string output_dir_TF=input_dir;
	output_dir_TF+="Fourier_transformed_data/";

	if (stat(output_dir_TF.c_str(),&file_stat)) mkdir(output_dir_TF.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
	
	mat save_mat_G=zeros<mat>(Nn,3);
	save_mat_G.col(0)=wn;
	save_mat_G.col(1)=Gr;
	save_mat_G.col(2)=Gi;
	string G_file_name("Fourier_transform_G_ascii.dat");
	string complete_file_name_G(output_dir_TF);
	complete_file_name_G+=G_file_name;
	//save_mat_G.save(complete_file_name_G.c_str(),raw_ascii);
	
	G_file_name.assign("Fourier_transform_G.dat");
	complete_file_name_G.assign(output_dir_TF);
	complete_file_name_G+=G_file_name;
	//save_mat_G.save(complete_file_name_G.c_str(),arma_binary);
	
	string cov_file_RR("covar_ReRe.dat");
	string complete_file_name_CRR(output_dir_TF);
	complete_file_name_CRR+=cov_file_RR;
	//CRR.save(complete_file_name_CRR.c_str(),arma_binary);
	
	string cov_file_RI("covar_ReIm.dat");
	string complete_file_name_CRI(output_dir_TF);
	complete_file_name_CRI+=cov_file_RI;
	//CRI.save(complete_file_name_CRI.c_str(),arma_binary);
	
	string cov_file_II("covar_ImIm.dat");
	string complete_file_name_CII(output_dir_TF);
	complete_file_name_CII+=cov_file_II;
	//CII.save(complete_file_name_CII.c_str(),arma_binary);

	jfit=0;
	if (!cutoff_wn_in.size())
	{
		int NC=3;
		
		if (col_Gi>0)
		{
			vec C1b(Nn-2), C2b(Nn-2), C3b(Nn-2), C4b(Nn-2);
			double wn1, wn2, reG1, reG2, imG1, imG2, denom2;
			for (j=1; j<Nn-1; j++)
			{
				wn1=wn(j);
				wn2=wn(j+1);
				reG1=Gr(j);
				reG2=Gr(j+1);
				imG1=Gi(j);
				imG2=Gi(j+1);
				denom2=wn1*wn1 - wn2*wn2;
				C1b(j-1)=-(imG1*pow(wn1,3) - imG2*pow(wn2,3))/denom2;
				C2b(j-1)=-(reG1*pow(wn1,4) - reG2*pow(wn2,4))/denom2;
				C3b(j-1)=((wn1*wn1)*(wn2*wn2)*(-imG1*wn1 + imG2*wn2))/denom2;
				C4b(j-1)=(-reG1*pow(wn1,4)*pow(wn2,2) + reG2*pow(wn1,2)*pow(wn2,4))/denom2;
			}
			
			int p, Nfitmin, jfitmin, jfitmax, NNfit, Nfit;
			jfitmin=2;
			Nfitmin=2*NC+4;
			jfitmax=Nn-Nfitmin;
			NNfit=jfitmax-jfitmin+1;
			
			vec M0v(NNfit), M1v(NNfit), M2v(NNfit), M3v(NNfit);
			mat invCG, A, X, CG;
			vec Mtmp;
			
			for (jfit=jfitmin; jfit<=jfitmax; jfit++)
			{
				Nfit=Nn_fit_max;
				if ((Nn-jfit+1)<Nn_fit_max)
					Nfit=Nn-jfit+1;
				
				X.zeros(2*Nfit,2*NC);
				
				for (j=1; j<=NC; j++)
				{
					for (p=jfit; p<=jfit+Nfit-1; p++)
					{
						X(2*(p-jfit),2*j-1)=pow(-1,j)/pow(wn(p-1),2*j);
						X(2*(p-jfit)+1,2*j-2)=pow(-1,j)/pow(wn(p-1),2*j-1);
					}
				}
				
				CG=COV.submat(2*jfit-2,2*jfit-2,2*(jfit+Nfit-1)-1,2*(jfit+Nfit-1)-1);
				
				invCG=inv(CG);
				//	invCG=inv_sympd(CG);
				A=trans(X)*invCG*X;
				A=0.5*(A+A.t());
				Mtmp=trans(X)*invCG*Gchi2.rows(2*jfit-2,2*(jfit+Nfit-1)-1);
				
				Mtmp=solve(A,Mtmp);
				
				M0v(jfit-jfitmin)=Mtmp(0);
				M1v(jfit-jfitmin)=Mtmp(1);
				M2v(jfit-jfitmin)=Mtmp(2);
				M3v(jfit-jfitmin)=Mtmp(3);
				
			}
			
			int Nv=Nn/16;
			if (Nv<2) Nv=2;
			
			vec varM0(NNfit-2*Nv), varM1(NNfit-2*Nv), varM2(NNfit-2*Nv), varM3(NNfit-2*Nv);
			
			for (j=Nv; j<NNfit-Nv; j++)
			{
				varM0(j-Nv)=var(M0v.rows(j-Nv,j+Nv));
				varM1(j-Nv)=var(M1v.rows(j-Nv,j+Nv));
				varM2(j-Nv)=var(M2v.rows(j-Nv,j+Nv));
				varM3(j-Nv)=var(M3v.rows(j-Nv,j+Nv));
			}
			
			uword j0, j1, j2, j3;
			varM0.min(j0);
			varM1.min(j1);
			varM2.min(j2);
			varM3.min(j3);
			
			j0=j0+Nv;
			j1=j1+Nv;
			j2=j2+Nv;
			j3=j3+Nv;
			
			Mfit.zeros(4);
			Mfit(0)=mean(M0v.rows(j0-Nv,j0+Nv));
			Mfit(1)=mean(M1v.rows(j1-Nv,j1+Nv));
			Mfit(2)=mean(M2v.rows(j2-Nv,j2+Nv));
			Mfit(3)=mean(M3v.rows(j3-Nv,j3+Nv));
			
			j=Nv;
			if (!boson)
			{
				while ((abs(mean(C1b.rows(j-Nv,j+Nv))-Mfit(0))/Mfit(0)>tol_mean_C1 || stddev(C1b.rows(j-Nv,j+Nv))/Mfit(0)>tol_std_C1) && j<Nn-Nv-3)
				{
					j=j+1;
				}
			}
			else
			{
				while ((abs(mean(C2b.rows(j-Nv,j+Nv))-Mfit(1))/Mfit(1)>tol_mean_C1 || stddev(C2b.rows(j-Nv,j+Nv))/Mfit(1)>tol_std_C1) && j<Nn-Nv-3)
				{
					j=j+1;
				}
			}
			jfit=j;
			
			cout<<"frequency range of asymptotic behavior: "<<wn(jfit-1)<<" to "<<wn(Nn-1)<<" (indices "<<jfit-1<<" to "<<Nn-1<<")"<<endl;
		}
		else
		{
			vec C2b(Nn-2), C4b(Nn-2);
			double wn1, wn2, reG1, reG2, denom2;
			for (j=1; j<Nn-1; j++)
			{
				wn1=wn(j);
				wn2=wn(j+1);
				reG1=Gr(j);
				reG2=Gr(j+1);
				denom2=wn1*wn1 - wn2*wn2;
				C2b(j-1)=-(reG1*pow(wn1,4) - reG2*pow(wn2,4))/denom2;
				C4b(j-1)=(-reG1*pow(wn1,4)*pow(wn2,2) + reG2*pow(wn1,2)*pow(wn2,4))/denom2;
			}
			
			int p, Nfitmin, jfitmin, jfitmax, NNfit, Nfit;
			jfitmin=2;
			Nfitmin=2*NC+4;
			jfitmax=Nn-Nfitmin;
			NNfit=jfitmax-jfitmin+1;
			
			vec M1v(NNfit), M3v(NNfit);
			mat invCG, A, X, CG;
			vec Mtmp;
			
			for (jfit=jfitmin; jfit<=jfitmax; jfit++)
			{
				Nfit=Nn_fit_max;
				if ((Nn-jfit+1)<Nn_fit_max)
					Nfit=Nn-jfit+1;
				
				X.zeros(Nfit,NC);
				
				for (j=1; j<=NC; j++)
				{
					for (p=jfit; p<=jfit+Nfit-1; p++)
					{
						X((p-jfit),j-1)=pow(-1,j)/pow(wn(p-1),2*j);
					}
				}
				
				CG=COV.submat(jfit-1,jfit-1,jfit+Nfit-2,jfit+Nfit-2);
				
				invCG=inv(CG);
				//		invCG=inv_sympd(CG);
				A=trans(X)*invCG*X;
				A=0.5*(A+A.t());
				Mtmp=trans(X)*invCG*Gchi2.rows(jfit-1,jfit+Nfit-2);
				
				Mtmp=solve(A,Mtmp);
				
				M1v(jfit-jfitmin)=Mtmp(0);
				M3v(jfit-jfitmin)=Mtmp(1);
				
			}
			
			int Nv=Nn/16;
			if (Nv<2) Nv=2;
			
			vec varM1(NNfit-2*Nv), varM3(NNfit-2*Nv);
			
			for (j=Nv; j<NNfit-Nv; j++)
			{
				varM1(j-Nv)=var(M1v.rows(j-Nv,j+Nv));
				varM3(j-Nv)=var(M3v.rows(j-Nv,j+Nv));
			}
			
			uword j1, j3;
			varM1.min(j1);
			varM3.min(j3);
			
			j1=j1+Nv;
			j3=j3+Nv;
			
			Mfit.zeros(2);
			Mfit(0)=mean(M1v.rows(j1-Nv,j1+Nv));
			Mfit(1)=mean(M3v.rows(j3-Nv,j3+Nv));
			
			j=Nv;
			while ((abs(mean(C2b.rows(j-Nv,j+Nv))-Mfit(0))/Mfit(0)>tol_mean_C1 || stddev(C2b.rows(j-Nv,j+Nv))/Mfit(0)>tol_std_C1) && j<Nn-Nv-3) j++;
			jfit=j;
		
			cout<<"frequency range of asymptotic behavior: "<<wn(jfit-1)<<" to "<<wn(Nn-1)<<" (indices "<<jfit-1<<" to "<<Nn-1<<")"<<endl;
		}
	}
	
	if (Nn-jfit<Nn_as_min)
	{
		jfit=0;
		if (!moments_provided && maxM>0) maxM=0;
	}

	return true;
}






bool OmegaMaxEnt_data::set_wc()
{
	bool use_nu_grid=false;
	
	cout<<setiosflags(ios::left);
	
	if (!step_omega_in.size() && non_uniform_grid && peak_exists)
	{
		double dw_default=0, dw_min;
		
		dw_min=2*dw_peak/R_peak_width_dw;
		
		if (SW_set)
		{
			dw_default=SW/(f_SW_std_omega*Rmin_SW_dw);
		}
		else if (std_omega)
		{
			dw_default=std_omega/Rmin_SW_dw;
		}
		else if (main_spectral_region_set)
		{
			dw_default=(wr-wl)/(f_SW_std_omega*Rmin_SW_dw);
		}
		else if (jfit)
		{
			dw_default=2*wn(jfit)/(R_wncutoff_wr*f_SW_std_omega*Rmin_SW_dw);
		}
		
		if (dw_min<dw_default) use_nu_grid=true;
	}
	
	if (step_omega_in.size() || !non_uniform_grid || !peak_exists || !use_nu_grid)
	{
		double dw;
		
		if (SW_set && SC_set && !main_spectral_region_set)
		{
			wl=SC-SW/2;
			wr=SC+SW/2;
			main_spectral_region_set=true;
		}
		else if (std_omega && SC_set && !main_spectral_region_set)
		{
			SW=f_SW_std_omega*std_omega;
			SW_set=true;
			wl=SC-SW/2;
			wr=SC+SW/2;
			main_spectral_region_set=true;
		}
		else if (main_spectral_region_set)
		{
			SW=wr-wl;
			SW_set=true;
		}
		else
		{
			cout<<"The central part of the grid cannot be defined. Not enough information available.\n";
			return false;
		}
		
//		if (!w_origin_in.size())
		if (!w_origin_set)
		{
			if (SC_set)
				w_origin=SC;
			else
				w_origin=0;
			w_origin_set=true;
		}
		
		if (step_omega_in.size())
		{
			dw=step_omega;
			if (peak_exists)
			{
				if (dw>dw_peak)
				{
					cout<<"warning: step is larger than the estimated width of the peak at low energy\n";
				}
			}
		}
		else
		{
			dw=SW/(f_SW_std_omega*Rmin_SW_dw);
			if (peak_exists)
			{
				if (dw_peak<dw)
					dw=dw_peak;
			}
		}
		
		if (w_origin<wl+dw || w_origin>wr-dw) w_origin=SC;
		
		uint Nwr=round((wr-w_origin)/dw);
		wr=w_origin+Nwr*dw;
		vec wcr=linspace<vec>(w_origin,wr,Nwr+1);
		uint Nwl=round((w_origin-wl)/dw);
		wl=w_origin-Nwl*dw;
		vec wcl=linspace<vec>(w_origin-dw,wl,Nwl);
		wcl=flipud(wcl);
		
		Nwc=Nwl+Nwr+1;
		
		wc.zeros(Nwc);
		wc.rows(0,Nwl-1)=wcl;
		wc.rows(Nwl,Nwc-1)=wcr;
		
		//cout<<"wl: "<<wl<<endl;
		//cout<<"wr: "<<wr<<endl;
		
		//wl=wc(0);
		//wr=wc(Nwc-1);
		dwl=dw;
		dwr=dw;
		wc_exists=true;
		main_spectral_region_set=true;
		
		return true;
	}
	else
	{
		if (!SW_set && !jfit && !std_omega)
		{
			cout<<"The central part of the grid cannot be defined. Not enough information available.\n";
			return false;
		}
		
		double SC_tmp=0;
	//	if (SC_in.size()) SC_tmp=SC;
		
		if (!w_origin_set)
		{
			w_origin=SC_tmp;
			w_origin_set=true;
		}
		
		double wr_tmp, wl_tmp;
		
		if (SW_in.size())
		{
			wr_tmp=SC_tmp+SW/2;
			wl_tmp=SC_tmp-SW/2;
		}
		else if (jfit)
		{
			wr_tmp=wn(jfit)/R_wncutoff_wr;
			wl_tmp=-wr_tmp;
		}
		else if (SW_set)
		{
			wr_tmp=SC_tmp+R_SW_wr*SW/2;
			wl_tmp=SC_tmp-R_SW_wr*SW/2;
		}
		else
		{
			wr_tmp=SC_tmp+R_SW_wr*f_SW_std_omega*std_omega/2;
			wl_tmp=SC_tmp-R_SW_wr*f_SW_std_omega*std_omega/2;
		}
		
		if (SW_set)
		{
			if (wr_tmp<SC_tmp+R_SW_wr*SW/2) wr_tmp=SC_tmp+R_SW_wr*SW/2;
			if (wl_tmp>SC_tmp-R_SW_wr*SW/2) wl_tmp=SC_tmp-R_SW_wr*SW/2;
			double wmax_tmp=SC_tmp+f_w_range*SW/2.0;
			if (wmax_tmp>0)
			{
				if (wr_tmp>wmax_tmp/R_wmax_wr_min) wr_tmp=wmax_tmp/R_wmax_wr_min;
			}
			double wmin_tmp=SC_tmp-f_w_range*SW/2.0;
			if (wmin_tmp<0)
			{
				if (wl_tmp<wmin_tmp/R_wmax_wr_min) wl_tmp=wmin_tmp/R_wmax_wr_min;
			}
		}
		
//		cout<<"wr_tmp: "<<wr_tmp<<endl;
		
		rowvec grid_params_r, grid_params_l;
		
		grid_params_l.zeros(3);
		grid_params_r.zeros(3);
		
//		if (peak_exists)
//		{
			grid_params_l(1)=2*dw_peak/R_peak_width_dw;
			grid_params_r(1)=2*dw_peak/R_peak_width_dw;
			w_origin=0;
//		}

		rowvec params_tmp;
		int j=2;
		grid_params_l(2)=grid_params_l(0)-R_Dw_dw*grid_params_l(1);
		grid_params_r(2)=grid_params_r(0)+R_Dw_dw*grid_params_r(1);
		
		if (SW_set)
		{
			if (grid_params_r(1)>SW/(f_SW_std_omega*Rmin_SW_dw))
			{
				grid_params_l(1)=SW/(f_SW_std_omega*Rmin_SW_dw);
				grid_params_r(1)=SW/(f_SW_std_omega*Rmin_SW_dw);
			}
		}
//		cout<<"wr_tmp: "<<wr_tmp<<endl;
//		cout<<"wl_tmp: "<<wl_tmp<<endl;
		
		while (grid_params_r(j)<wr_tmp)
		{
			params_tmp.zeros(j+3);
			params_tmp.cols(0,j)=grid_params_r;
			grid_params_r=params_tmp;
			j++;
			grid_params_r(j)=2*grid_params_r(j-2);
			//			cout<<grid_params_r(j)<<endl;
			j++;
			grid_params_r(j)=grid_params_r(j-2)+R_Dw_dw*grid_params_r(j-1);
			//			cout<<grid_params_r(j)<<endl;
		}
		j=2;
		while (grid_params_l(j)>wl_tmp)
		{
			params_tmp.zeros(j+3);
			params_tmp.cols(0,j)=grid_params_l;
			grid_params_l=params_tmp;
			j++;
			grid_params_l(j)=2*grid_params_l(j-2);
			//			cout<<grid_params_l(j)<<endl;
			j++;
			grid_params_l(j)=grid_params_l(j-2)-R_Dw_dw*grid_params_l(j-1);
			//			cout<<grid_params_l(j)<<endl;
		}
		
//		cout<<"grid_params_l\n"<<grid_params_l<<endl;
//		cout<<"grid_params_r\n"<<grid_params_r<<endl;
		
		grid_params_l=fliplr(grid_params_l);
		omega_grid_params=join_rows(grid_params_l.cols(0,grid_params_l.n_cols-2),grid_params_r.cols(2,grid_params_r.n_cols-1));
		
//		cout<<"omega_grid_params:\n"<<omega_grid_params<<endl;
		
//		cout<<"SC, SW: "<<setw(20)<<SC<<SW<<endl;
		
		return set_grid_from_params();
	}
}





bool OmegaMaxEnt_data::set_smooth_step_grid()
{
	if (!SC_set || !SW_set || !wc_exists) return false;
	
	int j;
	
	f_width_grid_dens=1;
	
	double Dwh=f_width_grid_dens*SW/2;
	double Dwh2=Dwh*Dwh;
	
	double wmin=SC-f_w_range*SW/2.0;
	int Nwlmax=(wc(0)-wmin)/(wc(1)-wc(0));
	
	vec cfs;
	vec x(3), y(3);
	double x0;
	
	x0=0;
	
	x(0)=wc(0)-Dwh;
	y(0)=2*(wc(1)-wc(0));
	x(1)=wc(1);
	y(1)=wc(1)-wc(0);
	x(2)=wc(2);
	y(2)=wc(2)-wc(1);
	
	polyfit(x, y, 2, x0, cfs);
	
	double w0=-cfs(1)/(2*cfs(0));
	
	if (cfs(0)<0 || (w0>wc(1) && w0<wc(2)))
	{
		x0=wc(0);
		cfs(0)=(wc(1)-wc(0))/Dwh2;
		cfs(1)=0;
		cfs(2)=wc(1)-wc(0);
	}
	
	vec vw(1), vDw(1);
	
	double dwtmp;
	vec wltmp(Nwlmax);
	wltmp(0)=wl;
	j=0;
	while (wltmp(j)>wmin)// && j<Nwlmax-1)
	{
		vw(0)=wltmp(j);
		polyval(x0, cfs, vw, vDw);
		dwtmp=vDw(0);
		wltmp(j+1)=wltmp(j)-dwtmp;
		j++;
	}
	int Nwl=j;
	
	double wmax=SC+f_w_range*SW/2.0;
	int Nwrmax=(wmax-wc(Nwc-1))/(wc(Nwc-1)-wc(Nwc-2));
	
	x0=0;
	
	x(0)=wc(Nwc-3);
	y(0)=wc(Nwc-2)-wc(Nwc-3);
	x(1)=wc(Nwc-2);
	y(1)=wc(Nwc-1)-wc(Nwc-2);
	x(2)=wc(Nwc-1)+Dwh;
	y(2)=2*(wc(Nwc-1)-wc(Nwc-2));
	
	polyfit(x, y, 2, x0, cfs);
	
	w0=-cfs(1)/(2*cfs(0));
	
	if (cfs(0)<0 || (w0>wc(Nwc-3) && w0<wc(Nwc-2)))
	{
		x0=wc(Nwc-1);
		cfs(0)=(wc(Nwc-1)-wc(Nwc-2))/Dwh2;
		cfs(1)=0;
		cfs(2)=wc(Nwc-1)-wc(Nwc-2);
	}
	
	vec wrtmp(Nwrmax);
	wrtmp(0)=wr;
	j=0;
	while (wrtmp(j)<wmax)// && j<Nwrmax-1)
	{
		vw(0)=wrtmp(j);
		polyval(x0, cfs, vw, vDw);
		dwtmp=vDw(0);
		wrtmp(j+1)=wrtmp(j)+dwtmp;
		j++;
	}
	int Nwr=j;
	
	Nw=Nwl+Nwc+Nwr;
	w.zeros(Nw);
	w.rows(0,Nwl-1)=flipud(wltmp.rows(1,Nwl));
	w.rows(Nwl,Nwl+Nwc-1)=wc;
	w.rows(Nwl+Nwc,Nw-1)=wrtmp.rows(1,Nwr);
	w_exists=true;
	Nw_lims.zeros(2);
	Nw_lims(0)=Nwl;
	Nw_lims(1)=Nwl+Nwc-1;
	w0l=SC;
	w0r=SC;
	ws.zeros(2);
	ws(0)=w0l;
	ws(1)=w0r;
	
	return true;
}






bool OmegaMaxEnt_data::set_omega_grid()
{
	if (!SC_set || !SW_set || !wc_exists) return false;
	
	int j;
	
	double wmin=SC-f_w_range*SW/2.0;
	if (wl<0)
	{
		if (wmin>R_wmax_wr_min*wl) wmin=R_wmax_wr_min*wl;
	}
	w0l=wl+sqrt(dwl*(wl-wmin));
	int Nul=ceil((w0l-wl)/dwl);
	w0l=wl+Nul*dwl;
	double dul=dwl/((wl-dwl-w0l)*(wl-w0l));
	wmin=-1.0/dul+w0l;
	ivec ul_int=linspace<ivec>(1,Nul,Nul);
//	vec ul=-dul*ul_int;
	vec ul(Nul);
	for (j=0; j<Nul; j++)
		ul(j)=-dul*ul_int(j);
	
	double wmax=SC+f_w_range*SW/2.0;
	if (wr>0)
	{
		if (wmax<R_wmax_wr_min*wr) wmax=R_wmax_wr_min*wr;
	}
	w0r=wr-sqrt(dwr*(wmax-wr));
	int Nur=ceil((wr-w0r)/dwr);
	w0r=wr-Nur*dwr;
	double dur=dwr/((wr-w0r)*(wr+dwr-w0r));
	wmax=1.0/dur+w0r;
	ivec ur_int=linspace<ivec>(Nur,1,Nur);
//	vec ur=dur*ur_int;
	vec ur(Nur);
	for (j=0; j<Nur; j++)
		ur(j)=dur*ur_int(j);
	
	vec wl_vec=1.0/ul+w0l;
	vec wr_vec=1.0/ur+w0r;
	Nw=Nul+Nur+Nwc;
	w.zeros(Nw);
	w.rows(0,Nul-1)=wl_vec;
	w.rows(Nul,Nul+Nwc-1)=wc;
	w.rows(Nul+Nwc,Nw-1)=wr_vec;
	w_exists=true;
	ws.zeros(2);
	ws(0)=w0l;
	ws(1)=w0r;
	Nw_lims.zeros(2);
	Nw_lims(0)=Nul;
	Nw_lims(1)=Nul+Nwc-1;
	
	Du_constant=true;
	
	return true;
}




bool OmegaMaxEnt_data::truncate_G_omega_n()
{
	double wn_max=wn(ind_cutoff_wn);
//	cout<<"wn_max: "<<wn_max<<endl;
	wn=wn.rows(0,ind_cutoff_wn);
	n=n.rows(0,ind_cutoff_wn);
	Nn=ind_cutoff_wn+1;
	if (Nn>Nn_max) // || Nn>=Nw)
	{
		cout<<"Number of Matsubara frequencies larger than the maximum allowed...\n";
		uvec diff_n=n.rows(1,Nn-1)-n.rows(0,Nn-2);
		if (diff_n.max()==1)
		{
			cout<<"Using a non-uniform Matsubara grid.\n";
			int i;
			int p=ceil(log2(Nn));
			int N0=pow(2,p);
			if (N0+1>Nn_all)
			{
				p=p-1;
				N0=pow(2,p);
			}
			int r=1;
			int N1=N0/pow(2,r);
			int N2=N1/2;
			Nn=N1+r*N2+1;
			n.zeros(Nn);
			n.rows(0,N1-1)=linspace<uvec>(0,N1-1,N1);
			uvec j=linspace<uvec>(0,Nn-N1-1,Nn-N1);
			uvec lj=j/N2;
			for (i=0; i<Nn-N1; i++) n(i+N1)=(j(i) % N2)*pow(2,lj(i)+1) + N1*pow(2,lj(i));
			if (!boson)
				wn=PI*tem*conv_to<vec>::from(2*n+1);
			else
				wn=2*PI*tem*conv_to<vec>::from(n);
			i=Nn-1;
			while (wn(i)>wn_max && i>0) i--;
			Nn=i+1;
			while (Nn>Nn_max) //|| Nn>=Nw)
			{
				r=r+1;
				N1=N0/pow(2,r);
				N2=N1/2;
				Nn=N1+r*N2+1;
				n.zeros(Nn);
				n.rows(0,N1-1)=linspace<uvec>(0,N1-1,N1);
				uvec j=linspace<uvec>(0,Nn-N1-1,Nn-N1);
				uvec lj=j/N2;
				for (i=0; i<Nn-N1; i++) n(i+N1)=(j(i) % N2)*pow(2,lj(i)+1) + N1*pow(2,lj(i));
				if (!boson)
					wn=PI*tem*conv_to<vec>::from(2*n+1);
				else
					wn=2*PI*tem*conv_to<vec>::from(n);
				i=Nn-1;
				while (wn(i)>wn_max && i>0) i--;
				Nn=i+1;
			}
			wn=wn.rows(0,Nn-1);
			n=n.rows(0,Nn-1);
			G=G.rows(n);
			Gr=real(G);
			Gi=imag(G);
			uvec even_ind=linspace<uvec>(0,2*Nn-2,Nn);
			uvec odd_ind=even_ind+1;
			
			Gchi2.zeros(2*Nn);
			Gchi2.rows(even_ind)=Gr;
			Gchi2.rows(odd_ind)=Gi;
			
			uvec ind_R=2*n;
			uvec ind_I=2*n+1;
			
			mat COVtmp(2*Nn,2*Nn);
			COVtmp(even_ind,even_ind)=COV(ind_R,ind_R);
			COVtmp(odd_ind,odd_ind)=COV(ind_I,ind_I);
			COVtmp(even_ind,odd_ind)=COV(ind_R,ind_I);
			COVtmp(odd_ind,even_ind)=COV(ind_I,ind_R);
			COV=COVtmp;
			if (cov_diag)
			{
				errGr=errGr.rows(n);
				errGi=errGi.rows(n);
				errG.zeros(2*Nn);
				errG.rows(even_ind)=errGr;
				errG.rows(odd_ind)=errGi;
			}
		}
		else
		{
			cout<<"Your Matsubara grid is not uniform. Either make it more sparse to reduce further the number of frequencies, or increase \"Nn_max\" in the file \"OmegaMaxEnt_other_params.dat\".\n";
			return false;
		}
	}
	else
	{
		G=G.rows(0,Nn-1);
		Gr=real(G);
		Gi=imag(G);
		Gchi2=Gchi2.rows(0,2*Nn-1);
		COV=COV.submat(0,0,2*Nn-1,2*Nn-1);
		if (cov_diag)
		{
			errG=errG.rows(0,2*Nn-1);
			errGr=errGr.rows(0,Nn-1);
			errGi=errGi.rows(0,Nn-1);
		}
	}
	
	return true;
}






bool OmegaMaxEnt_data::set_default_model()
{
	if (def_model_file.size())
	{
		cout<<"default model provided\n";
		vec w_def=def_data.col(0);
		vec def_m=def_data.col(1);
		if (def_m.min()<0)
		{
			cout<<"set_default_model() error: provided default model has a negative value.\n";
			return false;
		}
		
		uint Nw_def=w_def.n_rows;
		//double wmn=w(0), wmx=w(Nw-1);
		double wmn_def=w_def(0), wmx_def=w_def(Nw_def-1);
		int j1=0, j2=1, j3=2;
		double w0l_def=(2*w_def(j1)*w_def(j3)-w_def(j1)*w_def(j2)-w_def(j2)*w_def(j3))/(w_def(j1)-2*w_def(j2)+w_def(j3));
		j1=Nw_def-3;
		j2=Nw_def-2;
		j3=Nw_def-1;
		double w0r_def=(2*w_def(j1)*w_def(j3)-w_def(j1)*w_def(j2)-w_def(j2)*w_def(j3))/(w_def(j1)-2*w_def(j2)+w_def(j3));
		if (w0l_def<wmn_def || w0r_def>wmx_def)
		{
			//cout<<"user defined default model.\n";
			if (w_def(0)<wl && w_def(Nw_def-1)>wr)
			{
				//int R_wmn_def_SW=4;
				double wmn_def2=w(0), wmx_def2=w(Nw-1);
				//wmn_def2=wl-R_wmn_def_SW*SW;
				//wmx_def2=wr+R_wmn_def_SW*SW;
				int j=0;
				while (w_def(j)<wmn_def2 && j<Nw_def-1)	j++;
				int jl=j-1;
				if (jl<0) jl=0;
				j=Nw_def-1;
				while (w_def(j)>wmx_def2 && j>0) j--;
				int jr=j+1;
				if (jr>Nw_def-1) jr=Nw_def-1;
				
				if (def_m(jl+1)>def_m(jl) && def_m(jr-1)>def_m(jr))
				{
					//w_def=w_def.rows(jl,jr);
					//def_m=def_m.rows(jl,jr);
					//Nw_def=jr-jl+1;
					
					double l21=log(def_m(1)/def_m(0));
					double l31=log(def_m(2)/def_m(0));
					double w1=w_def(0);
					double w2=w_def(1);
					double w3=w_def(2);
					double wcg=(l21*(pow(w1,2)-pow(w3,2))-l31*(pow(w1,2)-pow(w2,2)))/(2*l31*(w2-w1)-2*l21*(w3-w1));
					double C1g=(pow(w1-wcg,2)-pow(w2-wcg,2))/l21;
					double C2g=exp(-pow(w1-wcg,2)/C1g)/def_m(0);
					
					l21=log(def_m(Nw_def-2)/def_m(Nw_def-3));
					l31=log(def_m(Nw_def-1)/def_m(Nw_def-3));
					w1=w_def(Nw_def-3);
					w2=w_def(Nw_def-2);
					w3=w_def(Nw_def-1);
					double wcd=(l21*(pow(w1,2)-pow(w3,2))-l31*(pow(w1,2)-pow(w2,2)))/(2*l31*(w2-w1)-2*l21*(w3-w1));
					double C1d=(pow(w1-wcd,2)-pow(w2-wcd,2))/l21;
					double C2d=exp(-pow(w1-wcd,2)/C1d)/def_m(Nw_def-3);
					
					vec gaussians_params(6);
					gaussians_params(0)=wcg;
					gaussians_params(1)=C1g;
					gaussians_params(2)=C2g;
					gaussians_params(3)=wcd;
					gaussians_params(4)=C1d;
					gaussians_params(5)=C2d;
					
					double dfg=-2*(w_def(0)-wcg)*def_m(0)/C1g;
					double dfd=-2*(w_def(Nw_def-1)-wcd)*def_m(Nw_def-1)/C1d;
					
					uint Nc=4*(Nw_def-1);
					vec coeffs_spline_def(Nc+1);
					coeffs_spline_def(0)=dfg;
					coeffs_spline_def(1)=dfd;
					coeffs_spline_def(Nc)=def_m(Nw_def-1);
					
					if (dfg>=0 && dfd<=0 && C1g>0 && C1d>0 && C2g>0 && C2d>0)
					{
						spline_coeffs_rel(w_def.memptr(),def_m.memptr(),Nw_def,coeffs_spline_def.memptr());
						if (!default_model_val_G(w, w_def, coeffs_spline_def, gaussians_params, default_model))
						{
							return false;
						}
					}
					else
					{
						cout<<"set_default_model(): problem with the user defined default model.\n";
						if (dfg<=0)
							cout<<"The derivative does not have the correct sign at the left boundary\n";
						if (dfd>=0)
							cout<<"The derivative does not have the correct sign at the right boundary\n";
						if ( (C1g<0 || C1d<0 || C2g<0 || C2d<0) && (w_def(1)>w(1) || w_def(Nw_def)<w(Nw)) )
							cout<<"Incorrect parameters found. Unable to extend the model to whole frequency range.\n";
						
						return false;
					}
					
				}
				else if (def_m(jl+1)==def_m(jl) && def_m(jr-1)==def_m(jr))
				{
					default_model.zeros(Nw);
					int l=1;
					for (j=1; j<Nw-1; j++)
					{
						while (l<Nw_def-1 && w(j)>=w_def(l)) l++;
						
						if (w(j)>=w_def(l-1) && w(j)<w_def(l))
							default_model(j)=(w(j)-w_def(l-1))*(def_m(l)-def_m(l-1))/(w_def(l)-w_def(l-1))+def_m(l-1);
						else if (w(j)==w_def(l))
							default_model(j)=def_m(l);
					}
				}
				else
				{
					cout<<"warning: the provided default model is not valid. The function must be decreasing \n";
					return false;
				}
			}
			else
			{
				cout<<"warning: the grid of the user defined default model must extend beyond the frequency range of the spectral function\n";
				return false;
			}
		}
		else
		{
			cout<<"provided default model was generated by this code.\n";
			
			uvec ind_xlims=zeros<uvec>(2);
			vec xs=zeros<vec>(2);
			xs(0)=w0l_def;
			xs(1)=w0r_def;
			int j=0;
			double dw1=w_def(j+1)-w_def(j);
			double dw2=w_def(j+2)-w_def(j+1);
			double rdw=dw2/dw1;
			while (abs(rdw-1.0)>tol_rdw)
			{
				j=j+1;
				dw1=dw2;
				dw2=w_def(j+2)-w_def(j+1);
				rdw=dw2/dw1;
			}
			ind_xlims(0)=j+1;
			j=Nw_def-3;
			dw1=w_def(j+1)-w_def(j);
			dw2=w_def(j+2)-w_def(j+1);
			rdw=dw2/dw1;
			while (abs(rdw-1.0)>tol_rdw)
			{
				j=j-1;
				dw2=dw1;
				dw1=w_def(j+1)-w_def(j);
				rdw=dw2/dw1;
			}
			ind_xlims(1)=j+1;
			
			vec coeffs;
			spline_G_part(w_def, ind_xlims, xs, def_m, coeffs);
			//spline_G_omega_u(w_def, ind_xlims, xs, def_m, coeffs);
			spline_val_G_part(w, w_def, ind_xlims, xs, coeffs, default_model);
			
		}
	}
	else
	{
		if (!default_model_center_in.size())
		{
			if (!boson)
			{
				if (M1_set)
					default_model_center=M1/M0;
				else
					default_model_center=SC;
			}
			else
				default_model_center=M0/M1n;
		}
		if (!default_model_shape_in.size())
			default_model_shape=2.0;
		if (!default_model_width_in.size())
		{
			if (std_omega)
				default_model_width=pow(2.0,1.0/default_model_shape)*std_omega;
			else
				default_model_width=pow(2.0,1.0/default_model_shape)*SW/2;
		}
		general_normal(w, default_model_center, default_model_width, default_model_shape, default_model);
		
	}
	
	default_model=abs(default_model);
	
	return true;
}







bool OmegaMaxEnt_data::set_default_model_chi()
{
	if (def_model_file.size())
	{
		cout<<"default model provided\n";
		vec w_def=def_data.col(0);
		vec def_m=def_data.col(1);
		uint Nw_def=w_def.n_rows;
		if (w_def(0)!=0)
		{
			cout<<"set_default_model_chi(): if \"Im(G) column in data file\" is smaller than 1, the first frequency of the default model must be 0\n";
			return false;
		}
		if (def_m.min()<0)
		{
			cout<<"set_default_model_chi() error: provided default model has a negative value.\n";
			return false;
		}
		//double wmx=w(Nw-1);
		
		double wmx_def=w_def(Nw-1);
		int j1, j2, j3;
		j1=Nw_def-3;
		j2=Nw_def-2;
		j3=Nw_def-1;
		double w0r_def=(2*w_def(j1)*w_def(j3)-w_def(j1)*w_def(j2)-w_def(j2)*w_def(j3))/(w_def(j1)-2*w_def(j2)+w_def(j3));
		if (w0r_def>wmx_def)
		{
			if (w_def(Nw_def-1)>wr)
			{
				int j=Nw_def-1;
				while (w_def(j)>wr && j>0)	 j--;
				int jr=j+1;
				if (jr>Nw_def-1) jr=Nw_def-1;
				
				if (def_m(jr-1)>def_m(jr))
				{
					w_def=w_def.rows(0,jr);
					def_m=def_m.rows(0,jr);
					Nw_def=jr+1;
					
					double l21, l31, w1, w2, w3;
					l21=log(def_m(Nw_def-2)/def_m(Nw_def-3));
					l31=log(def_m(Nw_def-1)/def_m(Nw_def-3));
					w1=w_def(Nw_def-3);
					w2=w_def(Nw_def-2);
					w3=w_def(Nw_def-1);
					double wcd=(l21*(pow(w1,2)-pow(w3,2))-l31*(pow(w1,2)-pow(w2,2)))/(2*l31*(w2-w1)-2*l21*(w3-w1));
					double C1d=(pow(w1-wcd,2)-pow(w2-wcd,2))/l21;
					double C2d=exp(-pow(w1-wcd,2)/C1d)/def_m(Nw_def-3);
					
					vec gaussians_params(3);
					gaussians_params(0)=wcd;
					gaussians_params(1)=C1d;
					gaussians_params(2)=C2d;
					
					double dfd=-2*(w_def(Nw_def-1)-wcd)*def_m(Nw_def-1)/C1d;
					
					uint Nc=4*(Nw_def-1);
					vec coeffs_spline_def(Nc+1);
					coeffs_spline_def(0)=0;
					coeffs_spline_def(1)=dfd;
					coeffs_spline_def(Nc)=def_m(Nw_def-1);
					
					if (dfd<=0 && C1d>0 && C2d>0)
					{
						spline_coeffs_rel(w_def.memptr(),def_m.memptr(),Nw_def,coeffs_spline_def.memptr());
						if (!default_model_val_chi(w, w_def, coeffs_spline_def, gaussians_params, default_model))
						{
							return false;
						}
					}
					else
					{
						cout<<"set_default_model_chi(): problem with the user defined default model.\n";
						if (dfd>=0)
							cout<<"The derivative does not have the correct sign at the right boundary\n";
						if ( (C1d<0 || C2d<0) && w_def(Nw_def)<w(Nw) )
							cout<<"Incorrect parameters found. Unable to extend the model to whole frequency range.\n";
						
						return false;
					}
				}
				else if (def_m(jr-1)==def_m(jr))
				{
					default_model.zeros(Nw);
					default_model(0)=def_m(0);
					int l=1;
					for (j=1; j<Nw-1; j++)
					{
						while (l<Nw_def-1 && w_def(l)<=w(j)) l++;
						
						if (w(j)>=w_def(l-1) && w(j)<w_def(l))
							default_model(j)=(w(j)-w_def(l-1))*(def_m(l)-def_m(l-1))/(w_def(l)-w_def(l-1))+def_m(l-1);
						else if (w(j)==w_def(l))
							default_model(j)=def_m(l);
					}
					//				rowvec dw_tmp(Nw-1);
					//				dw_tmp(0)=w(1)/2;
					//				dw_tmp.cols(1,Nw-2)=trans(w.rows(2,Nw-1)-w.rows(0,Nw-3))/2.0;
					//				mat int_def_m=dw_tmp*(default_model.rows(0,Nw-2))/PI;
					//				default_model=M1n*default_model/int_def_m(0);
				}
				else
				{
					cout<<"warning: the provided default model is not valid\n";
					return false;
				}
			}
			else
			{
				cout<<"warning: the grid of the user defined default model must extend beyond the frequency range of the spectral function\n";
				return false;
			}
		}
		else
		{
			cout<<"provided default model was generated by this code.\n";
			
			uvec ind_xlims(1);
			vec xs(1);
			xs(0)=w0r_def;
			int j;
			double dw1, dw2, rdw;
			j=Nw_def-3;
			dw1=w_def(j+1)-w_def(j);
			dw2=w_def(j+2)-w_def(j+1);
			rdw=dw2/dw1;
			while (abs(rdw-1.0)>tol_rdw)
			{
				j=j-1;
				dw2=dw1;
				dw1=w_def(j+1)-w_def(j);
				rdw=dw2/dw1;
			}
			ind_xlims(0)=j+1;
			
			vec coeffs;
			spline_chi_part(w_def, ind_xlims, xs, def_m, coeffs);
			spline_val_chi_part(w, w_def, ind_xlims, xs, coeffs, default_model);
			
		}
	}
	else
	{
		default_model_center=0;
		if (!default_model_shape_in.size())
			default_model_shape=2.0;
		if (!default_model_width_in.size())
		{
			if (std_omega)
				default_model_width=pow(2.0,1.0/default_model_shape)*std_omega;
			else
				default_model_width=pow(2.0,1.0/default_model_shape)*SW/2;
		}
		general_normal(w, default_model_center, default_model_width, default_model_shape, default_model);
		
	}
	
	default_model=abs(default_model);
	
	return true;
}










bool OmegaMaxEnt_data::Kernel_G_fermions_grid_transf_1()
{
	bool use_HF_exp=true;
	double fg=1.7;
	int pngmax=100;
	double fi=fg;
	double fr=fi;
	int pnmax=50;
	double fd=fg;
	int pndmax=pngmax;
	
	cout<<"defining kernel matrix...\n";
	
	Kcx.zeros(Nn,Nw);
 
	int Nug=Nw_lims(0);
	int Nud=Nw-Nw_lims(1)-1;
	
	vec ug, ud;
	if (Du_constant)
	{
		double dug=dwl/((wl-dwl-w0l)*(wl-w0l));
		vec ug_int=linspace<vec>(1,Nug,Nug);
		ug=-dug*ug_int;
		
		double dud=dwr/((wr-w0r)*(wr+dwr-w0r));
		vec ud_int=linspace<vec>(Nud,1,Nud);
		ud=dud*ud_int;
	}
	else
	{
		ug=1.0/(w.rows(0,Nug-1)-w0l);
		ud=1.0/(w.rows(Nw_lims(1)+1,Nw-1)-w0r);
	}

	mat MM;
	spline_matrix_grid_transf_G_part_1(w, Nw_lims, ws, MM);

	///////////////////////////////////////
	
	int Nint=Nw-1;
	int NCfs=4*Nint;
	
	int Nintg=Nug;
	int NCg=4*Nintg;
	
	mat Pa_g=zeros<mat>(Nintg,NCfs);
	mat Pb_g_r=zeros<mat>(Nintg,NCfs);
	mat Pc_g_r=zeros<mat>(Nintg,NCfs);
	mat Pd_g_r=zeros<mat>(Nintg,NCfs);
	
	int j;
	for (j=0; j<Nintg; j++)
	{
		Pa_g(j,4*j)=1;
		Pb_g_r(j,4*j+1)=1;
		Pc_g_r(j,4*j+2)=1;
		Pd_g_r(j,4*j+3)=1;
	}
	
	vec vDg=((w.rows(1,Nw_lims(0))-w0l)%(w.rows(0,Nw_lims(0)-1)-w0l))/(w.rows(0,Nw_lims(0)-1)-w.rows(1,Nw_lims(0)));
	mat Dg=diagmat(vDg);
	vec vDU=(w.rows(1,Nw_lims(0))-w0l)/(w.rows(0,Nw_lims(0)-1)-w.rows(1,Nw_lims(0)));
	mat DU=diagmat(vDU);
	
	mat Pb_g=Pb_g_r-3*DU*Pa_g;
	mat Pc_g=Pc_g_r+3*pow(DU,2)*Pa_g-2*DU*Pb_g_r;
	mat Pd_g=Pd_g_r-pow(DU,3)*Pa_g+pow(DU,2)*Pb_g_r-DU*Pc_g_r;
	
	Pa_g=pow(Dg,3)*Pa_g;
	Pb_g=pow(Dg,2)*Pb_g;
	Pc_g=Dg*Pc_g;
	
	int	Nintc=Nwc-1;
	
	mat	Pa_c=zeros<mat>(Nintc,NCfs);
	mat	Pb_c=zeros<mat>(Nintc,NCfs);
	mat	Pc_c=zeros<mat>(Nintc,NCfs);
	mat	Pd_c=zeros<mat>(Nintc,NCfs);
	
	for (j=0; j<Nintc; j++)
	{
		Pa_c(j,4*j+NCg)=1;
		Pb_c(j,4*j+1+NCg)=1;
		Pc_c(j,4*j+2+NCg)=1;
		Pd_c(j,4*j+3+NCg)=1;
	}
	
	vec vDc=1.0/(w.rows(Nw_lims(0)+1,Nw_lims(1))-w.rows(Nw_lims(0),Nw_lims(1)-1));
	mat Dc=diagmat(vDc);
	
	Pa_c=pow(Dc,3)*Pa_c;
	Pb_c=pow(Dc,2)*Pb_c;
	Pc_c=Dc*Pc_c;
	
	int NCgc=NCg+4*Nintc;
	
	int Nintd=Nud;
	
	mat Pa_d=zeros<mat>(Nintd,NCfs);
	mat Pb_d_r=zeros<mat>(Nintd,NCfs);
	mat Pc_d_r=zeros<mat>(Nintd,NCfs);
	mat Pd_d_r=zeros<mat>(Nintd,NCfs);
	
	for (j=0; j<Nintd; j++)
	{
		Pa_d(j,4*j+NCgc)=1;
		Pb_d_r(j,4*j+1+NCgc)=1;
		Pc_d_r(j,4*j+2+NCgc)=1;
		Pd_d_r(j,4*j+3+NCgc)=1;
	}
	
	vec vDd=((w.rows(Nw_lims(1)+1,Nw-1)-w0r)%(w.rows(Nw_lims(1),Nw-2)-w0r))/(w.rows(Nw_lims(1)+1,Nw-1)-w.rows(Nw_lims(1),Nw-2));
	mat Dd=diagmat(vDd);
	vDU=(w.rows(Nw_lims(1),Nw-2)-w0r)/(w.rows(Nw_lims(1)+1,Nw-1)-w.rows(Nw_lims(1),Nw-2));
	DU=diagmat(vDU);
	
	mat Pb_d=Pb_d_r-3*DU*Pa_d;
	mat Pc_d=Pc_d_r+3*pow(DU,2)*Pa_d-2*DU*Pb_d_r;
	mat Pd_d=Pd_d_r-pow(DU,3)*Pa_d+pow(DU,2)*Pb_d_r-DU*Pc_d_r;
	
	Pa_d=pow(Dd,3)*Pa_d;
	Pb_d=pow(Dd,2)*Pb_d;
	Pc_d=Dd*Pc_d;
	
	cx_mat Ka_g=zeros<cx_mat>(Nn,Nintg);
	cx_mat Kb_g=zeros<cx_mat>(Nn,Nintg);
	cx_mat Kc_g=zeros<cx_mat>(Nn,Nintg);
	cx_mat Kd_g=zeros<cx_mat>(Nn,Nintg);
	
	rowvec ug2(Nug+1);
	ug2.cols(0,Nug-1)=ug.t();
	ug2(Nug)=1.0/(wl-w0l);
	
	mat Wng=wn*ones<rowvec>(Nintg);
	mat Ug=ones<vec>(Nn)*ug2;
	
	dcomplex i(0,1);
	
	mat atang=atan((Wng % (Ug.cols(1,Nintg)-Ug.cols(0,Nintg-1)))/(1+w0l*(Ug.cols(1,Nintg)+Ug.cols(0,Nintg-1))+(pow(w0l,2)+pow(Wng,2)) % Ug.cols(0,Nintg-1) % Ug.cols(1,Nintg)));
	mat logg=log((1.0+2*w0l*Ug.cols(1,Nintg)+pow(Ug.cols(1,Nintg),2) % (pow(Wng,2)+pow(w0l,2)))/(1.0+2*w0l*Ug.cols(0,Nintg-1)+pow(Ug.cols(0,Nintg-1),2) % (pow(Wng,2)+pow(w0l,2))));
	
	Ka_g=-(Ug.cols(1,Nintg)-Ug.cols(0,Nintg-1))/pow(Wng+i*w0l,2)-i*(pow(Ug.cols(1,Nintg),2)-pow(Ug.cols(0,Nintg-1),2))/(2*(Wng+i*w0l))+atang/pow(Wng+i*w0l,3)+i*logg/(2*pow(Wng+i*w0l,3));
	
	Kb_g=-i*(Ug.cols(1,Nintg)-Ug.cols(0,Nintg-1))/(Wng+i*w0l)+i*atang/pow(Wng+i*w0l,2)-logg/(2*pow(Wng+i*w0l,2));
	
	Kc_g=-atang/(Wng+i*w0l)-i*logg/(2*(Wng+i*w0l));
	
	Kd_g=-i*atang+logg/2-log(Ug.cols(1,Nintg)/Ug.cols(0,Nintg-1));
	
	rowvec wg=trans(w.rows(0,Nug));
	//rowvec wg=1/ug2+w0l;
	mat Wg=ones<vec>(Nn)*wg;
	
	if (use_HF_exp)
	{
		double utmp;
		int jn, p;
		for (j=0; j<Nug; j++)
		{
			utmp=ug2(j);
			jn=0;
			while (abs(wn(jn)*utmp)<fg*abs(1+utmp*w0l) && jn<Nn-1) jn++;
			while (pow(wn(jn),pngmax)==0 && jn<Nn-1) jn++;
			if (jn<Nn-1)
			{
				Ka_g.submat(jn,j,Nn-1,j)=-i*(pow(Ug.submat(jn,j+1,Nn-1,j+1),2)-pow(Ug.submat(jn,j,Nn-1,j),2))/(2*(wn.rows(jn,Nn-1)+i*w0l))-(Ug.submat(jn,j+1,Nn-1,j+1)-Ug.submat(jn,j,Nn-1,j))/pow(wn.rows(jn,Nn-1)+i*w0l,2)+log(Ug.submat(jn,j+1,Nn-1,j+1)/Ug.submat(jn,j,Nn-1,j))/pow(i*wn.rows(jn,Nn-1)-w0l,3);
				Kb_g.submat(jn,j,Nn-1,j)=-i*(Ug.submat(jn,j+1,Nn-1,j+1)-Ug.submat(jn,j,Nn-1,j))/(wn.rows(jn,Nn-1)+i*w0l)+log(Ug.submat(jn,j+1,Nn-1,j+1)/Ug.submat(jn,j,Nn-1,j))/pow(i*wn.rows(jn,Nn-1)-w0l,2);
				Kc_g.submat(jn,j,Nn-1,j)=log(Ug.submat(jn,j+1,Nn-1,j+1)/Ug.submat(jn,j,Nn-1,j))/(i*wn.rows(jn,Nn-1)-w0l);
				Kd_g.submat(jn,j,Nn-1,j)=zeros<cx_vec>(Nn-jn);
				for (p=pngmax; p>=1; p--)
				{
					Ka_g.submat(jn,j,Nn-1,j)=Ka_g.submat(jn,j,Nn-1,j)+(pow(i,p)*pow(-1,p+1)*(pow(Wg.submat(jn,j+1,Nn-1,j+1),p)-pow(Wg.submat(jn,j,Nn-1,j),p))/(p*pow(wn.rows(jn,Nn-1),p)))/pow(i*wn.rows(jn,Nn-1)-w0l,3);
					Kb_g.submat(jn,j,Nn-1,j)=Kb_g.submat(jn,j,Nn-1,j)+(pow(i,p)*pow(-1,p+1)*(pow(Wg.submat(jn,j+1,Nn-1,j+1),p)-pow(Wg.submat(jn,j,Nn-1,j),p))/(p*pow(wn.rows(jn,Nn-1),p)))/pow(i*wn.rows(jn,Nn-1)-w0l,2);
					Kc_g.submat(jn,j,Nn-1,j)=Kc_g.submat(jn,j,Nn-1,j)+(pow(i,p)*pow(-1,p+1)*(pow(Wg.submat(jn,j+1,Nn-1,j+1),p)-pow(Wg.submat(jn,j,Nn-1,j),p))/(p*pow(wn.rows(jn,Nn-1),p)))/(i*wn.rows(jn,Nn-1)-w0l);
					Kd_g.submat(jn,j,Nn-1,j)=Kd_g.submat(jn,j,Nn-1,j)+pow(i,p)*pow(-1,p+1)*(pow(Wg.submat(jn,j+1,Nn-1,j+1),p)-pow(Wg.submat(jn,j,Nn-1,j),p))/(p*pow(wn.rows(jn,Nn-1),p));
				}
			}
		}
	}
	
	Ka_g=-Ka_g/(2*PI);
	Kb_g=-Kb_g/(2*PI);
	Kc_g=-Kc_g/(2*PI);
	Kd_g=-Kd_g/(2*PI);
	
	
	cx_mat Ka_c=zeros<cx_mat>(Nn,Nintc);
	cx_mat Kb_c=zeros<cx_mat>(Nn,Nintc);
	cx_mat Kc_c=zeros<cx_mat>(Nn,Nintc);
	cx_mat Kd_c=zeros<cx_mat>(Nn,Nintc);
	
	mat Wnc=wn*ones<rowvec>(Nintc);
	mat Wc=ones<vec>(Nn)*wc.t();
	
	mat logc=log((pow(Wnc,2)+pow(Wc.cols(1,Nintc),2))/(pow(Wnc,2)+pow(Wc.cols(0,Nintc-1),2)));
	mat atanc2=atan((Wnc % (Wc.cols(1,Nintc)-Wc.cols(0,Nintc-1)))/(Wc.cols(1,Nintc) % Wc.cols(0,Nintc-1)+pow(Wnc,2)));
	cx_mat logc2=logc/2+i*atanc2;
	cx_mat dWn=i*Wnc-Wc.cols(0,Nintc-1);
	mat dWc=Wc.cols(1,Nintc)-Wc.cols(0,Nintc-1);
	
	Ka_c=-pow(dWn,2) % dWc-dWn % pow(dWc,2)/2-pow(dWc,3)/3-pow(dWn,3) % logc2;
	Kb_c=-dWn % dWc-pow(dWc,2)/2-pow(dWn,2) % logc2;
	Kc_c=-dWc-dWn % logc2;
	Kd_c=-logc2;
	
	int Pmax=2*pnmax+4;
	imat MP;
	pascal(Pmax+1,MP);
	
	if (use_HF_exp)
	{
		double wtmp;
		int jni, jnr, p, l;
		for (j=0; j<Nwc-1; j++)
		{
			wtmp=abs(wc(j));
			if (abs(wc(j+1))>wtmp) wtmp=abs(wc(j+1));
			jni=0;
			while (wn(jni)<fi*wtmp && jni<Nn-1) jni++;
			while (pow(wn(jni),2*pnmax+1)==0 && jni<Nn-1) jni++;
			jnr=1;
			while (wn(jnr)<fr*wtmp && jnr<Nn-1) jnr++;
			while (pow(wn(jnr),2*pnmax)==0 && jnr<Nn-1) jnr++;
			if (jni<Nn-1 || jnr<Nn-1)
			{
				double wj=wc(j);
				double dw1=wc(j+1)-wj;
				vec dwp=zeros<vec>(Pmax);
				dwp(0)=dw1;
				dwp(1)=pow(dw1,2)+2*wj*dw1;
				for (p=3; p<=Pmax; p++)
				{
					dwp(p-1)=0;
					for (l=0; l<p; l++)
					{
						dwp(p-1)=dwp(p-1)+MP(l,p-l)*pow(dw1,p-l)*pow(wj,l);
					}
				}
				cx_vec vtmp;
				if (jni<Nn)
				{
					vtmp.zeros(Nn-jni);
					vtmp.set_real(real(Ka_c.submat(jni,j,Nn-1,j)));
					Ka_c.submat(jni,j,Nn-1,j)=vtmp;
					vtmp.set_real(real(Kb_c.submat(jni,j,Nn-1,j)));
					Kb_c.submat(jni,j,Nn-1,j)=vtmp;
					vtmp.set_real(real(Kc_c.submat(jni,j,Nn-1,j)));
					Kc_c.submat(jni,j,Nn-1,j)=vtmp;
					vtmp.set_real(real(Kd_c.submat(jni,j,Nn-1,j)));
					Kd_c.submat(jni,j,Nn-1,j)=vtmp;
					
					for (p=2*pnmax+1; p>=1; p-=2)
					{
						Ka_c.submat(jni,j,Nn-1,j)=Ka_c.submat(jni,j,Nn-1,j) + i*pow(-1,(p-1)/2)*(pow(wj,3)*dwp(p-1)/p - 3.0*pow(wj,2)*dwp(p)/(p+1) + 3.0*wj*dwp(p+1)/(p+2) - dwp(p+2)/(p+3))/pow(wn.rows(jni,Nn-1),p);
						Kb_c.submat(jni,j,Nn-1,j)=Kb_c.submat(jni,j,Nn-1,j) + i*pow(-1,(p-1)/2)*(-pow(wj,2)*dwp(p-1)/p + 2*wj*dwp(p)/(p+1) - dwp(p+1)/(p+2))/pow(wn.rows(jni,Nn-1),p);
						Kc_c.submat(jni,j,Nn-1,j)=Kc_c.submat(jni,j,Nn-1,j) + i*pow(-1,(p-1)/2)*(wj*dwp(p-1)/p - dwp(p)/(p+1))/pow(wn.rows(jni,Nn-1),p);
						Kd_c.submat(jni,j,Nn-1,j)=Kd_c.submat(jni,j,Nn-1,j) - i*pow(-1,(p-1)/2)*dwp(p-1)/(p*pow(wn.rows(jni,Nn-1),p));
					}
				}
				if (jnr<Nn)
				{
					vtmp.zeros(Nn-jnr);
					vtmp.set_imag(imag(Ka_c.submat(jnr,j,Nn-1,j)));
					Ka_c.submat(jnr,j,Nn-1,j)=vtmp;
					vtmp.set_imag(imag(Kb_c.submat(jnr,j,Nn-1,j)));
					Kb_c.submat(jnr,j,Nn-1,j)=vtmp;
					vtmp.set_imag(imag(Kc_c.submat(jnr,j,Nn-1,j)));
					Kc_c.submat(jnr,j,Nn-1,j)=vtmp;
					vtmp.set_imag(imag(Kd_c.submat(jnr,j,Nn-1,j)));
					Kd_c.submat(jnr,j,Nn-1,j)=vtmp;
					for (p=2*pnmax; p>=2; p-=2)
					{
						Ka_c.submat(jnr,j,Nn-1,j)=Ka_c.submat(jnr,j,Nn-1,j) + pow(-1,(p-2)/2)*(pow(wj,3)*dwp(p-1)/p - 3*pow(wj,2)*dwp(p)/(p+1) + 3*wj*dwp(p+1)/(p+2) - dwp(p+2)/(p+3))/pow(wn.rows(jnr,Nn-1),p);
						Kb_c.submat(jnr,j,Nn-1,j)=Kb_c.submat(jnr,j,Nn-1,j) + pow(-1,(p-2)/2)*(-pow(wj,2)*dwp(p-1)/p + 2*wj*dwp(p)/(p+1) - dwp(p+1)/(p+2))/pow(wn.rows(jnr,Nn-1),p);
						Kc_c.submat(jnr,j,Nn-1,j)=Kc_c.submat(jnr,j,Nn-1,j) + pow(-1,(p-2)/2)*(wj*dwp(p-1)/p - dwp(p)/(p+1))/pow(wn.rows(jnr,Nn-1),p);
						Kd_c.submat(jnr,j,Nn-1,j)=Kd_c.submat(jnr,j,Nn-1,j) - pow(-1,(p-2)/2)*dwp(p-1)/(p*pow(wn.rows(jnr,Nn-1),p));
					}
				}
			}
		}
	}
	
	Ka_c=Ka_c/(2*PI);
	Kb_c=Kb_c/(2*PI);
	Kc_c=Kc_c/(2*PI);
	Kd_c=Kd_c/(2*PI);
	
	cx_mat Ka_d=zeros<cx_mat>(Nn,Nintd);
	cx_mat Kb_d=zeros<cx_mat>(Nn,Nintd);
	cx_mat Kc_d=zeros<cx_mat>(Nn,Nintd);
	cx_mat Kd_d=zeros<cx_mat>(Nn,Nintd);
	
	rowvec ud2(Nud+1);
	ud2.cols(1,Nud)=ud.t();
	ud2(0)=1.0/(wr-w0r);
	
	mat Wnd=wn*ones<rowvec>(Nintd);
	mat Ud=ones<vec>(Nn)*ud2;
	
	mat atand=atan((Wnd % (Ud.cols(1,Nintd)-Ud.cols(0,Nintd-1)))/(1.0+w0r*(Ud.cols(1,Nintd)+Ud.cols(0,Nintd-1))+(pow(w0r,2)+pow(Wnd,2)) % Ud.cols(0,Nintd-1) % Ud.cols(1,Nintd)));
	mat logd=log((1+2*w0r*Ud.cols(1,Nintd)+pow(Ud.cols(1,Nintd),2) % (pow(Wnd,2)+pow(w0r,2)))/(1+2*w0r*Ud.cols(0,Nintd-1)+pow(Ud.cols(0,Nintd-1),2) % (pow(Wnd,2)+pow(w0r,2))));
	
	Ka_d=-(Ud.cols(1,Nintd)-Ud.cols(0,Nintd-1))/pow(Wnd+i*w0r,2) - i*(pow(Ud.cols(1,Nintd),2)-pow(Ud.cols(0,Nintd-1),2))/(2*(Wnd+i*w0r))+atand/pow(Wnd+i*w0r,3) +i*logd/(2*pow(Wnd+i*w0r,3));
	
	Kb_d=-i*(Ud.cols(1,Nintd)-Ud.cols(0,Nintd-1))/(Wnd+i*w0r)+i*atand/pow(Wnd+i*w0r,2)-logd/(2*pow(Wnd+i*w0r,2));
	
	Kc_d=-atand/(Wnd+i*w0r)-i*logd/(2*(Wnd+i*w0r));
	
	Kd_d=-i*atand-log(Ud.cols(1,Nintd)/Ud.cols(0,Nintd-1))+logd/2;
	
	//rowvec wd=1/ud2+w0r;
	rowvec wd=trans(w.rows(Nw_lims(1),Nw-1));
	mat Wd=ones<vec>(Nn)*wd;
	
	if (use_HF_exp)
	{
		double utmp;
		int jn, p;
		for (j=0; j<Nud; j++)
		{
			utmp=ud2(j);
			jn=0;
			while (abs(wn(jn)*utmp)<fd*abs(1+utmp*w0r) && jn<Nn-1) jn++;
			while (pow(wn(jn),pndmax) && jn<Nn-1) jn++;
			if (jn<Nn-1)
			{
				Ka_d.submat(jn,j,Nn-1,j)=-(Ud.submat(jn,j+1,Nn-1,j+1)-Ud.submat(jn,j,Nn-1,j))/pow(wn.rows(jn,Nn-1)+i*w0r,2) - i*(pow(Ud.submat(jn,j+1,Nn-1,j+1),2)-pow(Ud.submat(jn,j,Nn-1,j),2))/(2*(wn.rows(jn,Nn-1)+i*w0r))+log(Ud.submat(jn,j+1,Nn-1,j+1)/Ud.submat(jn,j,Nn-1,j))/pow(i*wn.rows(jn,Nn-1)-w0r,3);
				Kb_d.submat(jn,j,Nn-1,j)=-i*(Ud.submat(jn,j+1,Nn-1,j+1)-Ud.submat(jn,j,Nn-1,j))/(wn.rows(jn,Nn-1)+i*w0r)+log(Ud.submat(jn,j+1,Nn-1,j+1)/Ud.submat(jn,j,Nn-1,j))/pow(i*wn.rows(jn,Nn-1)-w0r,2);
				Kc_d.submat(jn,j,Nn-1,j)=log(Ud.submat(jn,j+1,Nn-1,j+1)/Ud.submat(jn,j,Nn-1,j))/(i*wn.rows(jn,Nn-1)-w0r);
				Kd_d.submat(jn,j,Nn-1,j)=zeros<cx_vec>(Nn-jn);
				for (p=pndmax; p>=1; p--)
				{
					Ka_d.submat(jn,j,Nn-1,j)=Ka_d.submat(jn,j,Nn-1,j)+(pow(i,p)*pow(-1,p+1)*(pow(Wd.submat(jn,j+1,Nn-1,j+1),p)-pow(Wd.submat(jn,j,Nn-1,j),p))/(p*pow(wn.rows(jn,Nn-1),p)))/pow(i*wn.rows(jn,Nn-1)-w0r,3);
					Kb_d.submat(jn,j,Nn-1,j)=Kb_d.submat(jn,j,Nn-1,j)+(pow(i,p)*pow(-1,p+1)*(pow(Wd.submat(jn,j+1,Nn-1,j+1),p)-pow(Wd.submat(jn,j,Nn-1,j),p))/(p*pow(wn.rows(jn,Nn-1),p)))/pow(i*wn.rows(jn,Nn-1)-w0r,2);
					Kc_d.submat(jn,j,Nn-1,j)=Kc_d.submat(jn,j,Nn-1,j)+(pow(i,p)*pow(-1,p+1)*(pow(Wd.submat(jn,j+1,Nn-1,j+1),p)-pow(Wd.submat(jn,j,Nn-1,j),p))/(p*pow(wn.rows(jn,Nn-1),p)))/(i*wn.rows(jn,Nn-1)-w0r);
					Kd_d.submat(jn,j,Nn-1,j)=Kd_d.submat(jn,j,Nn-1,j)+pow(i,p)*pow(-1,p+1)*(pow(Wd.submat(jn,j+1,Nn-1,j+1),p)-pow(Wd.submat(jn,j,Nn-1,j),p))/(p*pow(wn.rows(jn,Nn-1),p));
				}
			}
		}
	}
	
	Ka_d=-Ka_d/(2*PI);
	Kb_d=-Kb_d/(2*PI);
	Kc_d=-Kc_d/(2*PI);
	Kd_d=-Kd_d/(2*PI);
	
	cx_mat KG=(Ka_g*Pa_g+Kb_g*Pb_g+Kc_g*Pc_g+Kd_g*Pd_g)*MM;
	cx_mat KC=(Ka_c*Pa_c+Kb_c*Pb_c+Kc_c*Pc_c+Kd_c*Pd_c)*MM;
	cx_mat KD=(Ka_d*Pa_d+Kb_d*Pb_d+Kc_d*Pc_d+Kd_d*Pd_d)*MM;
	
	Kcx=KG+KC+KD;
	
	

	rowvec Knorm_a_g=zeros<rowvec>(Nintg);
	rowvec Knorm_b_g=zeros<rowvec>(Nintg);
	rowvec Knorm_c_g=zeros<rowvec>(Nintg);
	rowvec Knorm_d_g=zeros<rowvec>(Nintg);
	
	Knorm_a_g=-(pow(ug2.cols(1,Nug),2)-pow(ug2.cols(0,Nug-1),2))/2;
	Knorm_b_g=-(ug2.cols(1,Nug)-ug2.cols(0,Nug-1));
	Knorm_c_g=-log(ug2.cols(1,Nug)/ug2.cols(0,Nug-1));
	Knorm_d_g=1.0/ug2.cols(1,Nug)-1.0/ug2.cols(0,Nug-1);
	
	rowvec KM0g=(Knorm_a_g*Pa_g+Knorm_b_g*Pb_g+Knorm_c_g*Pc_g+Knorm_d_g*Pd_g)*MM/(2*PI);
	
	rowvec Knorm_a_c=zeros<rowvec>(Nintc);
	rowvec Knorm_b_c=zeros<rowvec>(Nintc);
	rowvec Knorm_c_c=zeros<rowvec>(Nintc);
	rowvec Knorm_d_c=zeros<rowvec>(Nintc);
	
	Knorm_a_c=trans(pow(wc.rows(1,Nwc-1)-wc.rows(0,Nwc-2),4)/4);
	Knorm_b_c=trans(pow(wc.rows(1,Nwc-1)-wc.rows(0,Nwc-2),3)/3);
	Knorm_c_c=trans(pow(wc.rows(1,Nwc-1)-wc.rows(0,Nwc-2),2)/2);
	Knorm_d_c=trans(wc.rows(1,Nwc-1)-wc.rows(0,Nwc-2));
	
	rowvec KM0c=(Knorm_a_c*Pa_c+Knorm_b_c*Pb_c+Knorm_c_c*Pc_c+Knorm_d_c*Pd_c)*MM/(2*PI);
	
	rowvec Knorm_a_d=zeros<rowvec>(Nintd);
	rowvec Knorm_b_d=zeros<rowvec>(Nintd);
	rowvec Knorm_c_d=zeros<rowvec>(Nintd);
	rowvec Knorm_d_d=zeros<rowvec>(Nintd);
	
	Knorm_a_d=-(pow(ud2.cols(1,Nud),2)-pow(ud2.cols(0,Nud-1),2))/2;
	Knorm_b_d=-(ud2.cols(1,Nud)-ud2.cols(0,Nud-1));
	Knorm_c_d=-log(ud2.cols(1,Nud)/ud2.cols(0,Nud-1));
	Knorm_d_d=1.0/ud2.cols(1,Nud)-1.0/ud2.cols(0,Nud-1);
	
	rowvec KM0d=(Knorm_a_d*Pa_d+Knorm_b_d*Pb_d+Knorm_c_d*Pc_d+Knorm_d_d*Pd_d)*MM/(2*PI);
	
	rowvec KM0=KM0g+KM0c+KM0d;
	
	
	rowvec KM1_a_g=zeros<rowvec>(Nintg);
	rowvec KM1_b_g=zeros<rowvec>(Nintg);
	rowvec KM1_c_g=zeros<rowvec>(Nintg);
	rowvec KM1_d_g=zeros<rowvec>(Nintg);
	
	KM1_a_g=Knorm_b_g;
	KM1_b_g=Knorm_c_g;
	KM1_c_g=Knorm_d_g;
	KM1_d_g=(1.0/pow(ug2.cols(1,Nug),2)-1.0/pow(ug2.cols(0,Nug-1),2))/2;
	
	rowvec KM1g_tmp=(KM1_a_g*Pa_g+KM1_b_g*Pb_g+KM1_c_g*Pc_g+KM1_d_g*Pd_g)*MM/(2*PI);
	
	rowvec KM1g=w0l*KM0g+KM1g_tmp;
	
	mat Wjc=diagmat(wc.rows(0,Nintc-1));
	
	rowvec KM1_a_c_tmp=trans(pow(wc.rows(1,Nwc-1)-wc.rows(0,Nwc-2),5)/5.0);
	
	rowvec KM1_a_c=KM1_a_c_tmp+Knorm_a_c*Wjc;
	rowvec KM1_b_c=Knorm_a_c+Knorm_b_c*Wjc;
	rowvec KM1_c_c=Knorm_b_c+Knorm_c_c*Wjc;
	rowvec KM1_d_c=Knorm_c_c+Knorm_d_c*Wjc;
	
	rowvec KM1c=(KM1_a_c*Pa_c+KM1_b_c*Pb_c+KM1_c_c*Pc_c+KM1_d_c*Pd_c)*MM/(2*PI);
	
	rowvec KM1_a_d=zeros<rowvec>(Nintd);
	rowvec KM1_b_d=zeros<rowvec>(Nintd);
	rowvec KM1_c_d=zeros<rowvec>(Nintd);
	rowvec KM1_d_d=zeros<rowvec>(Nintd);
	
	KM1_a_d=Knorm_b_d;
	KM1_b_d=Knorm_c_d;
	KM1_c_d=Knorm_d_d;
	KM1_d_d=(1.0/pow(ud2.cols(1,Nud),2)-1.0/pow(ud2.cols(0,Nud-1),2))/2;
	
	rowvec KM1d_tmp=(KM1_a_d*Pa_d+KM1_b_d*Pb_d+KM1_c_d*Pc_d+KM1_d_d*Pd_d)*MM/(2*PI);
	
	rowvec KM1d=w0r*KM0d+KM1d_tmp;
	
	rowvec KM1=KM1g+KM1c+KM1d;
	
	rowvec KM2_a_g=zeros<rowvec>(Nintg);
	rowvec KM2_b_g=zeros<rowvec>(Nintg);
	rowvec KM2_c_g=zeros<rowvec>(Nintg);
	rowvec KM2_d_g=zeros<rowvec>(Nintg);
	
	KM2_a_g=Knorm_c_g;
	KM2_b_g=Knorm_d_g;
	KM2_c_g=KM1_d_g;
	KM2_d_g=(1.0/pow(ug2.cols(1,Nug),3)-1.0/pow(ug2.cols(0,Nug-1),3))/3;
	
	rowvec KM2g_tmp=(KM2_a_g*Pa_g+KM2_b_g*Pb_g+KM2_c_g*Pc_g+KM2_d_g*Pd_g)*MM/(2*PI);
	
	rowvec KM2g=pow(w0l,2)*KM0g+2*w0l*KM1g_tmp+KM2g_tmp;
	
	rowvec KM2_a_c_tmp=trans(pow(wc.rows(1,Nwc-1)-wc.rows(0,Nwc-2),6)/6);
	
	rowvec KM2_a_c=KM2_a_c_tmp+2*KM1_a_c_tmp*Wjc+Knorm_a_c*pow(Wjc,2);
	rowvec KM2_b_c=KM1_a_c_tmp+2*Knorm_a_c*Wjc+Knorm_b_c*pow(Wjc,2);
	rowvec KM2_c_c=Knorm_a_c+2*Knorm_b_c*Wjc+Knorm_c_c*pow(Wjc,2);
	rowvec KM2_d_c=Knorm_b_c+2*Knorm_c_c*Wjc+Knorm_d_c*pow(Wjc,2);
	
	rowvec KM2c=(KM2_a_c*Pa_c+KM2_b_c*Pb_c+KM2_c_c*Pc_c+KM2_d_c*Pd_c)*MM/(2*PI);
	
	rowvec KM2_a_d=zeros<rowvec>(Nintd);
	rowvec KM2_b_d=zeros<rowvec>(Nintd);
	rowvec KM2_c_d=zeros<rowvec>(Nintd);
	rowvec KM2_d_d=zeros<rowvec>(Nintd);
	
	KM2_a_d=Knorm_c_d;
	KM2_b_d=Knorm_d_d;
	KM2_c_d=KM1_d_d;
	KM2_d_d=(1.0/pow(ud2.cols(1,Nud),3)-1.0/pow(ud2.cols(0,Nud-1),3))/3;
	
	rowvec KM2d_tmp=(KM2_a_d*Pa_d+KM2_b_d*Pb_d+KM2_c_d*Pc_d+KM2_d_d*Pd_d)*MM/(2*PI);
	
	rowvec KM2d=pow(w0r,2)*KM0d+2*w0r*KM1d_tmp+KM2d_tmp;
	
	rowvec KM2=KM2g+KM2c+KM2d;
	
	rowvec KM3_a_g=zeros<rowvec>(Nintg);
	rowvec KM3_b_g=zeros<rowvec>(Nintg);
	rowvec KM3_c_g=zeros<rowvec>(Nintg);
	rowvec KM3_d_g=zeros<rowvec>(Nintg);
	
	KM3_a_g=Knorm_d_g;
	KM3_b_g=KM1_d_g;
	KM3_c_g=KM2_d_g;
	KM3_d_g=(1.0/pow(ug2.cols(1,Nug),4)-1.0/pow(ug2.cols(0,Nug-1),4))/4;
	
	rowvec KM3g_tmp=(KM3_a_g*Pa_g+KM3_b_g*Pb_g+KM3_c_g*Pc_g+KM3_d_g*Pd_g)*MM/(2*PI);
	
	rowvec KM3g=pow(w0l,3)*KM0g+3*pow(w0l,2)*KM1g_tmp+3*w0l*KM2g_tmp+KM3g_tmp;
	
	rowvec KM3_a_c_tmp=trans(pow(wc.rows(1,Nwc-1)-wc.rows(0,Nwc-2),7)/7);
	
	rowvec KM3_a_c=KM3_a_c_tmp+3*KM2_a_c_tmp*Wjc+3*KM1_a_c_tmp*pow(Wjc,2)+Knorm_a_c*pow(Wjc,3);
	rowvec KM3_b_c=KM2_a_c_tmp+3*KM1_a_c_tmp*Wjc+3*Knorm_a_c*pow(Wjc,2)+Knorm_b_c*pow(Wjc,3);
	rowvec KM3_c_c=KM1_a_c_tmp+3*Knorm_a_c*Wjc+3*Knorm_b_c*pow(Wjc,2)+Knorm_c_c*pow(Wjc,3);
	rowvec KM3_d_c=Knorm_a_c+3*Knorm_b_c*Wjc+3*Knorm_c_c*pow(Wjc,2)+Knorm_d_c*pow(Wjc,3);
	
	rowvec KM3c=(KM3_a_c*Pa_c+KM3_b_c*Pb_c+KM3_c_c*Pc_c+KM3_d_c*Pd_c)*MM/(2*PI);
	
	rowvec KM3_a_d=zeros<rowvec>(Nintd);
	rowvec KM3_b_d=zeros<rowvec>(Nintd);
	rowvec KM3_c_d=zeros<rowvec>(Nintd);
	rowvec KM3_d_d=zeros<rowvec>(Nintd);
	
	KM3_a_d=Knorm_d_d;
	KM3_b_d=KM1_d_d;
	KM3_c_d=KM2_d_d;
	KM3_d_d=(1.0/pow(ud2.cols(1,Nud),4)-1.0/pow(ud2.cols(0,Nud-1),4))/4;
	
	rowvec KM3d_tmp=(KM3_a_d*Pa_d+KM3_b_d*Pb_d+KM3_c_d*Pc_d+KM3_d_d*Pd_d)*MM/(2*PI);
	
	rowvec KM3d=pow(w0r,3)*KM0d+3*pow(w0r,2)*KM1d_tmp+3*w0r*KM2d_tmp+KM3d_tmp;
	
	rowvec KM3=KM3g+KM3c+KM3d;
	
	KM.zeros(4,Nw);
	 
	KM.row(0)=KM0;
	KM.row(1)=KM1;
	KM.row(2)=KM2;
	KM.row(3)=KM3;
	
	K.zeros(2*Nn,Nw);
	uvec even_ind=linspace<uvec>(0,2*Nn-2,Nn);
	K.rows(even_ind)=real(Kcx);
	K.rows(even_ind+1)=imag(Kcx);
	
	cout<<"kernel matrix defined.\n";
	
	/*
	 //check moments
	 
	 rowvec x0={-2, 0, 3};
	 rowvec s0={0.8, 0.3, 1};
	 rowvec wgt={1,1,1};
	 
	 vec test_A;
	 sum_gaussians(w, x0, s0, wgt, test_A);
	
 	cout<<KM*test_A<<endl;
	 */
	
	return true;
}






bool OmegaMaxEnt_data::diagonalize_covariance()
{
	mat VM, WM, VG, WG;
	
	if (NM>0)
	{
		if (covm_diag)
		{
			VM.eye(NM,NM);
			WM=diagmat(1.0/errM);
		}
		else
		{
			int ind0_M=1;
			if (boson)	ind0_M=0;
			vec COVM_eig;
			mat VMtmp;
			if (!eig_sym(COVM_eig,VMtmp,COVM.submat(ind0_M,ind0_M,NM-1,NM-1)))
			{
				cout<<"diagonalize_covariance(): diagonalization of moments covariance matrix failed\n";
				return false;
			}
			if (COVM_eig.min()<=0)
			{
				cout<<"error: the moments covariance matrix has non-positive eigenvalues\n";
				return false;
			}
			VM.zeros(NM,NM);
			if (!boson) VM(0,0)=1;
			VM.submat(ind0_M,ind0_M,NM-1,NM-1)=VMtmp;
			vec errM_eig(NM);
			if (!boson) errM_eig(0)=sqrt(COVM(0,0));
			errM_eig.rows(ind0_M,NM-1)=sqrt(COVM_eig);
			WM=diagmat(1.0/errM_eig);
			
//			cout<<"COVM:\n"<<COVM<<endl;
//			cout<<"VM.t()*COVM*VM:\n"<<VM.t()*COVM*VM<<endl;
//			cout<<"COVM_eig:\n"<<COVM_eig<<endl;
		}
		
		if (!boson && covm_diag && M1_in.size() && M3_in.size() && !M2_in.size() && !eval_moments && maxM>2)
		{
			NM=3;
			M_ord.zeros(NM);
			M_ord(0)=0;
			M_ord(1)=1;
			M_ord(2)=3;
			M.zeros(NM);
			M(0)=M0;
			M(1)=M1;
			M(2)=M3;
			errM.zeros(NM);
			errM(0)=errM0;
			errM(1)=errM1;
			errM(2)=errM3;
			WM=diagmat(1.0/errM);
			M_V=WM*M;
			uvec indM={0,1,3};
			KM_V=WM*KM.rows(indM);
		}
		else
		{
			KM=KM.rows(0,NM-1);
			M_V=WM*VM.t()*M;
			KM_V=WM*VM.t()*KM;
		}
		
		KM_V=KM_V.cols(ind0,Nw-2);
		KM=KM.cols(ind0,Nw-2);
	}

	if (cov_diag)
	{
		VG.eye(2*Nn,2*Nn);
		WG=diagmat(1.0/errG);
	}
	else
	{
		vec COVG_eig;
		mat VGtmp;
		if (!eig_sym(COVG_eig,VGtmp,COV))
		{
			cout<<"diagonalize_covariance(): diagonalization of covariance matrix failed\n";
			return false;
		}
		mat PV=eye<mat>(2*Nn,2*Nn);
		PV=flipud(PV);
		COVG_eig=flipud(COVG_eig);
		VG=VGtmp*PV;
		vec GVtmp=VG.t()*Gchi2;
		mat SGV=diagmat(sign(GVtmp));
		SGV=SGV-abs(SGV)+eye<mat>(2*Nn,2*Nn);
		VG=VG*SGV;
		WG=diagmat(1.0/sqrt(COVG_eig));
		if (COVG_eig.min()<=0)
		{
			cout<<"error: the covariance matrix has non-positive eigenvalues\n";
			return false;
		}
	}
	
	G_V=WG*VG.t()*Gchi2;
	KG_V=WG*VG.t()*K.cols(ind0,Nw-2);
	K=K.cols(ind0,Nw-2);
	
	if (NM>0)
	{
		GM=join_vert(M_V, G_V);
		KGM=join_vert(KM_V, KG_V);
	}
	else
	{
		GM=G_V;
		KGM=KG_V;
	}
	NGM=GM.n_rows;
	
	cout<<"number of terms in chi2: "<<NGM<<endl;
	
	return true;
}









bool OmegaMaxEnt_data::set_G_omega_n_bosons()
{
	int j;
	
	signG=1;
	
	wn=green_data.col(0);
	Nn=green_data.n_rows;
	Gr=green_data.col(col_Gr-1);
	if (col_Gi>0)
		Gi=green_data.col(col_Gi-1);
	else
		Gi.zeros(Nn);
	
	if (static_part_G_in.size()) Gr=Gr-static_part_G;
	
	indG_0=0, indG_f=Nn-1;
	wn_sign_change=false;
	wn_inversed=false;
	if (wn(0)*wn(Nn-1)<0)
	{
		wn_sign_change=true;
		j=1;
		while (wn(j)*wn(Nn-1)<0) j++;
		if ((Nn-j)>=j+1)
		{
			indG_0=j;
			indG_f=Nn-1;
		}
		else
		{
			indG_0=0;
			indG_f=j-1;
		}
		wn=wn.rows(indG_0,indG_f);
		Nn=wn.n_rows;
		Gr=Gr.rows(indG_0,indG_f);
		Gi=Gi.rows(indG_0,indG_f);
	}
	if (wn(1)<0 )
	{
		wn=abs(wn);
		Gi=-Gi;
	}
	if (wn(0)>wn(Nn-1))
	{
		wn_inversed=true;
		wn=flipud(wn);
		Gr=flipud(Gr);
		Gi=flipud(Gi);
	}
	if (Gr.max()*Gr.min()<0)
	{
		cout<<"error: Re[G] must not change sign.\n";
		return false;
	}
	if (Gr.max()>0) signG=-1;
	Gr=signG*Gr;
	Gi=signG*Gi;
	
	G.zeros(Nn);
	G.set_real(Gr);
	G.set_imag(Gi);
	
	for (j=0; j<Nn-1; j++)
	{
		if ((wn(j+1)-wn(j))<0)
		{
			cout<<"error: Matsubara frequency is not strictly increasing\n";
			return false;
		}
	}
	
	if (col_Gi>0)
	{
		uvec even_ind=linspace<uvec>(0,2*Nn-2,Nn);
		uvec odd_ind=even_ind+1;
		
		Gchi2.zeros(2*Nn);
		Gchi2.rows(even_ind)=Gr;
		Gchi2.rows(odd_ind)=Gi;
	}
	else
	{
		Gchi2=Gr;
	}
	
	cout<<"Number of Matsubara frequencies in the Green function: "<<Nn<<endl;
	
	double tem_wn=wn(1)/(2*PI);
	
	n.zeros(Nn);
	vec ntmp;
	if (tem_in.size()==0)
	{
		if (  (wn(1)-floor(wn(1)))==0  && (wn(2)-floor(wn(2)))==0 )
		{
			cout<<"If the Matsubara frequency are given by index, temperature must be provided.\n";
			return false;
		}
		tem=tem_wn;
		ntmp=round(wn/(2*PI*tem));
		for (j=0; j<Nn; j++)	n(j)=ntmp(j);
	}
	else if ( (abs(tem-tem_wn)/tem)>tol_tem )
	{
		if ( wn(1)-floor(wn(1)) )
		{
			ntmp=round(wn/(2*PI*tem_wn));
			for (j=0; j<Nn; j++)	n(j)=ntmp(j);
			cout<<"warning: temperature in file "<<data_file_name<<" is different from temperature given in parameter file "<<input_params_file_name<<endl;
			cout<<"provided temperature: "<<tem<<endl;
			cout<<"temperature extracted from first finite Matsubara frequency: "<<tem_wn<<endl;
		}
		else if ( (wn(2)-floor(wn(2)))==0 )
		{
			ntmp=round(wn);
			for (j=0; j<Nn; j++) n(j)=ntmp(j);
		}
		else
		{
			cout<<"oops! First finite Matsubara frequency is an integer but second is not!\n";
			n=linspace<uvec>(0,Nn-1,Nn);
		}
	}
	else
	{
		ntmp=round(wn/(2*PI*tem_wn));
		for (j=0; j<Nn; j++) n(j)=ntmp(j);
	}
	
	wn=2*PI*tem*conv_to<vec>::from(n);
	
	return true;
}










bool OmegaMaxEnt_data::compute_moments_omega_n()
{
	int j, NC=3;
	//bool sol_found;
	
	cout<<"COMPUTING MOMENTS\n";
	
	vec C1b(Nn-2), C2b(Nn-2), C3b(Nn-2), C4b(Nn-2);
	double wn1, wn2, reG1, reG2, imG1, imG2, denom2;
	for (j=1; j<Nn-1; j++)
	{
		wn1=wn(j);
		wn2=wn(j+1);
		reG1=Gr(j);
		reG2=Gr(j+1);
		imG1=Gi(j);
		imG2=Gi(j+1);
		denom2=wn1*wn1 - wn2*wn2;
		C1b(j-1)=-(imG1*pow(wn1,3) - imG2*pow(wn2,3))/denom2;
		C2b(j-1)=-(reG1*pow(wn1,4) - reG2*pow(wn2,4))/denom2;
		C3b(j-1)=((wn1*wn1)*(wn2*wn2)*(-imG1*wn1 + imG2*wn2))/denom2;
		C4b(j-1)=(-reG1*pow(wn1,4)*pow(wn2,2) + reG2*pow(wn1,2)*pow(wn2,4))/denom2;
	}
	
	int p, Nfitmin, jfitmin, jfitmax, NNfit, Nfit;
	jfitmin=2;
	Nfitmin=2*NC+4;
	jfitmax=Nn-Nfitmin;
	NNfit=jfitmax-jfitmin+1;
	
	//int jfitmax_fin=jfitmax;
	
	vec M0v(NNfit), M1v(NNfit), M2v(NNfit), M3v(NNfit);
	mat invCG, A, X, CG;
	vec Mtmp;
	
//	char UPLO='U';
//	int NA=2*NC, NRHS=1, INFO;
	
//	mat LC, invLC;
	for (jfit=jfitmin; jfit<=jfitmax; jfit++)
	{
		Nfit=Nn_fit_max;
		if ((Nn-jfit+1)<Nn_fit_max)
			Nfit=Nn-jfit+1;
		
		X.zeros(2*Nfit,2*NC);
		
		for (j=1; j<=NC; j++)
		{
			for (p=jfit; p<=jfit+Nfit-1; p++)
			{
				X(2*(p-jfit),2*j-1)=pow(-1,j)/pow(wn(p-1),2*j);
				X(2*(p-jfit)+1,2*j-2)=pow(-1,j)/pow(wn(p-1),2*j-1);
			}
		}
		
		CG=COV.submat(2*jfit-2,2*jfit-2,2*(jfit+Nfit-1)-1,2*(jfit+Nfit-1)-1);

//		if (jfit==jfitmin) cout<<CG.submat(0,0,10,10);
	
//		LC=chol(CG);
//		invLC=inv(LC);
//		invCG=invLC*invLC.t();
		
		invCG=inv(CG);
		//	invCG=inv_sympd(CG);
		A=trans(X)*invCG*X;
		A=0.5*(A+A.t());
		Mtmp=trans(X)*invCG*Gchi2.rows(2*jfit-2,2*(jfit+Nfit-1)-1);
		
		//		dposv_(&UPLO, &NA, &NRHS, A.memptr(), &NA, Mtmp.memptr(), &NA, &INFO);
		Mtmp=solve(A,Mtmp);
//		sol_found=solve(Mtmp,A,Mtmp);
		
//		if (jfitmax_fin==jfitmax && (!sol_found || rcond(A)<EPSILON))
//		{
//			jfitmax_fin=jfit;
//			//break;
//		}
		
		M0v(jfit-jfitmin)=Mtmp(0);
		M1v(jfit-jfitmin)=Mtmp(1);
		M2v(jfit-jfitmin)=Mtmp(2);
		M3v(jfit-jfitmin)=Mtmp(3);
	}
	
	int Nv=Nn/16;
	if (Nv<2) Nv=2;
	
//	NNfit=jfitmax_fin-jfitmin;
	
	vec varM0(NNfit-2*Nv), varM1(NNfit-2*Nv), varM2(NNfit-2*Nv), varM3(NNfit-2*Nv);
	
	for (j=Nv; j<NNfit-Nv; j++)
	{
		varM0(j-Nv)=var(M0v.rows(j-Nv,j+Nv));
		varM1(j-Nv)=var(M1v.rows(j-Nv,j+Nv));
		varM2(j-Nv)=var(M2v.rows(j-Nv,j+Nv));
		varM3(j-Nv)=var(M3v.rows(j-Nv,j+Nv));
	}
	
	uword j0, j1, j2, j3;
	varM0.min(j0);
	varM1.min(j1);
	varM2.min(j2);
	varM3.min(j3);
	
	j0=j0+Nv;
	j1=j1+Nv;
	j2=j2+Nv;
	j3=j3+Nv;
	
	Mfit.zeros(4);
	Mfit(0)=mean(M0v.rows(j0-Nv,j0+Nv));
	Mfit(1)=mean(M1v.rows(j1-Nv,j1+Nv));
	Mfit(2)=mean(M2v.rows(j2-Nv,j2+Nv));
	Mfit(3)=mean(M3v.rows(j3-Nv,j3+Nv));
	
	
	int jfit0;
	j=Nv;
	if (!boson)
	{
		while ((abs(mean(C1b.rows(j-Nv,j+Nv))-Mfit(0))/Mfit(0)>tol_mean_C1 || stddev(C1b.rows(j-Nv,j+Nv))/Mfit(0)>tol_std_C1) && j<Nn-Nv-3)
		{
			j=j+1;
		}
	}
	else
	{
		while ((abs(mean(C2b.rows(j-Nv,j+Nv))-Mfit(1))/Mfit(1)>tol_mean_C1 || stddev(C2b.rows(j-Nv,j+Nv))/Mfit(1)>tol_std_C1) && j<Nn-Nv-3)
		{
			j=j+1;
		}
	}
	jfit0=j;
	
	if (!boson)
		jfit=j0+jfitmin-1;
	else
		jfit=j1+jfitmin-1;
	Nfit=Nn_fit_fin;
	if (Nfit>Nn-jfit+1)
		Nfit=Nn-jfit+1;
	
	X.zeros(2*Nfit,2*NC);
	
	for (j=1; j<=NC; j++)
	{
		for (p=jfit; p<=jfit+Nfit-1; p++)
		{
			X(2*(p-jfit),2*j-1)=pow(-1,j)/pow(wn(p-1),2*j);
			X(2*(p-jfit)+1,2*j-2)=pow(-1,j)/pow(wn(p-1),2*j-1);
		}
	}
	
	CG=COV.submat(2*jfit-2,2*jfit-2,2*(jfit+Nfit-1)-1,2*(jfit+Nfit-1)-1);

//	LC=chol(CG);
//	invLC=inv(LC);
//	invCG=invLC*invLC.t();
	
	invCG=inv(CG);
	//	invCG=inv_sympd(CG);
	A=trans(X)*invCG*X;
	A=0.5*(A+A.t());
	mat COVMtmp=inv(A);
	COVMfit=COVMtmp.submat(0,0,3,3);
	
	
	jfit=jfit0;
	
	cout<<"frequency range of asymptotic behavior: "<<wn(jfit-1)<<" to "<<wn(Nn-1)<<" (indices "<<jfit-1<<" to "<<Nn-1<<")"<<endl;
	
	cout<<"norm extracted from high frequencies of G: "<<Mfit(0)<<endl;
	cout<<"1st moment extracted from high frequencies: "<<Mfit(1)<<endl;
	cout<<"2nd moment extracted from high frequencies: "<<Mfit(2)<<endl;
	cout<<"3rd moment extracted from high frequencies: "<<Mfit(3)<<endl;
	
	double std_omega_tmp;
	double var_omega;
	if (!boson)
	{
		if (M0)
		{
			if (abs(Mfit(0)-M0)/Mfit(0)>tol_norm)
			{
				if (M0_in.size())
					cout<<"warning: norm of spectral function is different from provided one.\n";
				else
					cout<<"warning: spectral function is not normalized.\n";
			}
		}
		
		var_omega=Mfit(2)/Mfit(0)-pow(Mfit(1)/Mfit(0),2);
		
		if (var_omega>0)
			std_omega_tmp=sqrt(var_omega);
		else
		{
			cout<<"Negative variance found during computation of moments.\n";
			return false;
		}
		
		if (!moments_provided)
		{
			M=Mfit;
			NM=M.n_rows;
			COVM.zeros(NM,NM);
			COVM(0,0)=pow(errM0,2);
			COVM.submat(1,1,NM-1,NM-1)=COVMfit.submat(1,1,NM-1,NM-1);
			//		COVM=COVMfit;
			M(0)=M0;
			M1=M(1);
			M2=M(2);
			M3=M(3);
			covm_diag=false;
			M1_set=true;
			M2_set=true;
		}
		else
		{
			if (abs(M1-Mfit(1))/std_omega_tmp>tol_M1)
				cout<<"warning: first moment different from provided one\n";
			if (M2_in.size())
			{
				if (abs(M2-Mfit(2))/Mfit(2)>tol_M2)
					cout<<"warning: second moment different from provided one\n";
				
				if (M3_in.size())
				{
					if (abs(M3-Mfit(3))/pow(std_omega_tmp,3)>tol_M3)
					{
						cout<<"warning: third moment different from provided one\n";
					}
				}
				else
				{
					M=Mfit;
					M(0)=M0;
					M(1)=M1;
					M(2)=M2;
					M3=M(3);
					COVMtmp=COVM.submat(0,0,2,2);
					COVM.zeros(4,4);
					COVM.submat(0,0,2,2)=COVMtmp;
					COVM(3,3)=COVMfit(3,3);
				}
			}
			else if (M3_in.size())
			{
				M=Mfit;
				M(0)=M0;
				M(1)=M1;
				M2=M(2);
				M2_set=true;
				M(3)=M3;
				COVMtmp=COVM.submat(0,0,1,1);
				COVM.zeros(4,4);
				COVM.submat(0,0,1,1)=COVMtmp;
				if (errM3_in.size())
				{
					COVM(2,2)=COVMfit(2,2);
					COVM(3,3)=pow(errM3,2);
				}
				else
				{
					COVM.submat(2,2,3,3)=COVMfit.submat(2,2,3,3);
					covm_diag=false;
				}
			}
			else
			{
				M=Mfit;
				M(0)=M0;
				M(1)=M1;
				M2=M(2);
				M2_set=true;
				M3=M(3);
				COVMtmp=COVM.submat(0,0,1,1);
				COVM.zeros(4,4);
				COVM.submat(0,0,1,1)=COVMtmp;
				COVM.submat(2,2,3,3)=COVMfit.submat(2,2,3,3);
				covm_diag=false;
			}
		}
		
		NM=M.n_rows;
		M_ord=linspace<vec>(0,NM-1,NM);
		if (errM.n_rows<NM)
		{
			vec errM_tmp=errM;
			errM.zeros(NM);
			errM.rows(0,errM_tmp.n_rows-1)=errM_tmp;
			COVMtmp=COVM.submat(errM_tmp.n_rows,errM_tmp.n_rows,NM-1,NM-1);
			errM.rows(errM_tmp.n_rows,NM-1)=sqrt(COVMtmp.diag());
		}
		
		if (!std_omega)
		{
			var_omega=M2/M0-pow(M1/M0,2);
			std_omega=sqrt(var_omega);
		}
		if (!SC_set)
		{
			SC=M1/M0;
			SC_set=true;
		}
		if (!SW_set)
		{
			SW=f_SW_std_omega*std_omega;
			SW_set=true;
		}
		if (M(2)<0)
		{
			M=M.rows(0,1);
			COVM=COVM.submat(0,0,1,1);
			NM=2;
		}
	}
	else
	{
		if (M0_in.size())
		{
			if (M0)
			{
				if (abs((Mfit(0)-M0)/Mfit(0))>tol_norm)
					cout<<"warning: norm of spectral function is different from provided one.\n";
			}
		}
		
		var_omega=Mfit(1)/M1n-pow(Mfit(0)/M1n,2);
		if (var_omega>0)
			std_omega_tmp=sqrt(var_omega);
		else
		{
			cout<<"Negative variance found during computation of moments.\n";
			return false;
		}
		
		if (!moments_provided)
		{
			M=Mfit.rows(0,2);
			NM=3;
			COVM=COVMfit.submat(0,0,NM-1,NM-1);
			if (M0_in.size())
			{
				M(0)=M0;
				COVM.row(0)=zeros<rowvec>(NM);
				COVM.col(0)=zeros<vec>(NM);
				COVM(0,0)=pow(errM0,2);
			}
			else if (M0)
			{
				if (abs((Mfit(0)-M0)/Mfit(0))<tol_norm)
				{
					M(0)=M0;
					COVM.row(0)=zeros<rowvec>(NM);
					COVM.col(0)=zeros<vec>(NM);
					COVM(0,0)=pow(errM0,2);
				}
				else
					M0=M(0);
			}
			else
			{
				M(0)=M0;
				COVM.row(0)=zeros<rowvec>(NM);
				COVM.col(0)=zeros<vec>(NM);
				COVM(0,0)=pow(errM0,2);
			}
			M1=M(1);
			M2=M(2);
			covm_diag=false;
			M1_set=true;
			M2_set=true;
		}
		else
		{
			if (abs((M1-Mfit(1))/Mfit(1))>tol_M1)
				cout<<"warning: first moment different from provided one\n";
			if (M2_in.size())
			{
				if (M2)
				{
					if (abs((M2-Mfit(2))/Mfit(2))>tol_M2)
						cout<<"warning: second moment different from provided one\n";
				}
			}
			else
			{
				M=Mfit.rows(0,2);
				NM=3;
				M(0)=M0;
				M(1)=M1;
				M2=M(2);
				M2_set=true;
				vec errMtmp=errM;
				errM.zeros(NM);
				errM.rows(0,1)=errMtmp;
				errM(2)=sqrt(COVMfit(2,2));
				COVM=diagmat(square(errM));
			}
		}
		M_ord=linspace<vec>(0,NM-1,NM);
		
		if (!std_omega)
		{
			var_omega=M1/M1n-pow(M0/M1n,2);
			if (var_omega>0)
				std_omega=sqrt(var_omega);
			else
			{
				cout<<"Negative variance found during computation of moments.\n";
				return false;
			}
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
	}
	
	if (Nn-jfit<Nn_as_min)
	{
		jfit=0;
		if (!moments_provided && maxM>0) maxM=0;
	}
	
	return true;
}
















bool OmegaMaxEnt_data::test_low_energy_peak_bosons()
{
	
	peak_exists=false;

		cout<<"Looking for a peak in the spectral function at low energy...\n";
		
		int nmax=n(Nn-1);
		
		int NCnmin=1;
		int NCnmax=3;
		int NNCn=NCnmax-NCnmin+1;
		ivec NCn=linspace<ivec>(NCnmin,NCnmax,NNCn);
		
		int NCpmin=1;
		int NCpmax=15;
		if (NCpmax>(nmax-NCn.max()-2))
		{
			NCpmax=nmax-NCn.max()-2;
		}
		int NNCp=NCpmax-NCpmin+1;
		ivec NCp=linspace<ivec>(NCpmin,NCpmax,NNCp);
		
		int DNwn=2;
		
		int p0_min=2;
		int p0_max=2;
		int Np=p0_max-p0_min+1;
		ivec p0=linspace<ivec>(p0_min,p0_max,Np);
		
		imat NCmin(Np,2);
		
		vec M0_inc_opt(Np,fill::zeros);
		vec M1_pk_opt(Np,fill::zeros);
		vec M2_pk_opt(Np,fill::zeros);
		vec varM2_pk_opt(Np,fill::zeros);
		
		mat X, CG, P, invCG, AM, BM, AMP, BMP, MP, Mtmp, chi2tmp;
		mat M0_inc(NNCp,NNCn), M1_pk(NNCp,NNCn), M2_pk(NNCp,NNCn), chi2_pk(NNCp,NNCn);
		int j, l, m, p, Nfit;
		uword q;
		vec Gtmp, Glftmp, diffG;
		rowvec maxX;
		//double vartmp;
		
		for (q=0; q<Np; q++)
		{
			M0_inc.zeros();
			M1_pk.zeros();
			M2_pk.zeros();
			chi2_pk.zeros();
			
			for (l=0; l<NNCp; l++)
			{
				for (m=0; m<NNCn; m++)
				{
					Nfit=NCp(l)+NCn(m)+1+DNwn;
					
					X.zeros(2*Nfit,2*(NCp(l)+NCn(m)+1));
					
					for (j=NCp(l); j>=-NCn(m); j--)
					{
						for (p=p0(q); p<=Nfit+p0(q)-1; p++)
						{
							X(2*(p-p0(q)),2*NCp(l)-2*j+1)=pow(wn(p-1),2*j);
							X(2*(p-p0(q))+1,2*NCp(l)-2*j)=pow(wn(p-1),2*j+1);
						}
					}
					
					Gtmp=Gchi2.rows(2*p0(q)-2,2*(Nfit+p0(q))-3);
					CG=COV.submat(2*p0(q)-2,2*p0(q)-2,2*(Nfit+p0(q))-3,2*(Nfit+p0(q))-3);
					
					invCG=inv(CG);
					//	invCG=inv_sympd(CG);
					AM=(X.t())*invCG*X;
					BM=(X.t())*invCG*Gtmp;
					//					BM2=BM;
					//					AM2=AM;
					
					Mtmp=solve(AM,BM);
					
					
					Glftmp=X*Mtmp;
					
					diffG=Gtmp-Glftmp;
					
					chi2tmp=((diffG.t())*invCG*diffG)/(2.0*Nfit);

					chi2_pk(l,m)=chi2tmp(0,0);
					M0_inc(l,m)=-Mtmp(2*NCp(l)+1);
					M1_pk(l,m)=-Mtmp(2*NCp(l)+2);
					M2_pk(l,m)=-Mtmp(2*NCp(l)+3);
				}
			}
			
			uint NvM=2;
			uint NvarM=NNCp-2*NvM;
			mat varM2_pk(NvarM,NNCn,fill::zeros);
			for (j=NvM; j<NNCp-NvM; j++)
			{
				varM2_pk.row(j-NvM)=var(M2_pk.rows(j-NvM,j+NvM));
			}
			
			uword indpM0, indnM0;
			
			double varM2_pk_min=varM2_pk.min(indpM0,indnM0);
			
			l=indpM0+NvM;
			m=indnM0;
			
			NCmin(q,0)=l;
			NCmin(q,1)=m;
			
			M0_inc_opt(q)=M0_inc(l,m);
			M1_pk_opt(q)=M1_pk(l,m);
			M2_pk_opt(q)=M2_pk(l,m);
			varM2_pk_opt(q)=varM2_pk_min;
			
		}
		
		double varM2_min=varM2_pk_opt.min(q);
		
		double peak_weight=-Gr(0)-M0_inc_opt(q);
		double peak_center=M1_pk_opt(q)/peak_weight;
		double var_peak=M2_pk_opt(q)/peak_weight-pow(peak_center,2);
		double peak_width=0;
		if (var_peak>0)		peak_width=sqrt(var_peak);
		//double m2_lf_min=M2_pk_opt(q);
	
		if (peak_width>100*EPSILON && M0_inc_opt(q)>100*EPSILON && varM2_min/M2_pk_opt(q)<varM2_peak_max && peak_weight>peak_weight_min*M1n)
		{
			peak_exists=true;
			
			l=NCmin(q,0);
			m=NCmin(q,1);
			
			Nfit=NCp(l)+NCn(m)+1+DNwn;
			
			X.zeros(2*Nfit,2*(NCp(l)+NCn(m)+1));
			
			for (j=NCp(l); j>=-NCn(m); j--)
			{
				for (p=p0(q); p<=Nfit+p0(q)-1; p++)
				{
					X(2*(p-p0(q)),2*NCp(l)-2*j+1)=pow(wn(p-1),2*j);
					X(2*(p-p0(q))+1,2*NCp(l)-2*j)=pow(wn(p-1),2*j+1);
				}
			}
			
			Gtmp=Gchi2.rows(2*p0(q)-2,2*(Nfit+p0(q))-3);
			CG=COV.submat(2*p0(q)-2,2*p0(q)-2,2*(Nfit+p0(q))-3,2*(Nfit+p0(q))-3);
			
			invCG=inv(CG);
			//	invCG=inv_sympd(CG);
			AM=(X.t())*invCG*X;
			AM=0.5*(AM.t()+AM);
			BM=(X.t())*invCG*Gtmp;
			
			Mtmp=solve(AM,BM);
						
			Glftmp=X*Mtmp;
			
			mat COVpeak=inv(AM);
			
			//double err_norm_peak=sqrt(COVpeak(2*NCp(l)+1,2*NCp(l)+1));
			//double err_M1_peak=sqrt(COVpeak(2*NCp(l)+2,2*NCp(l)+2));
			//double err_M2_peak=sqrt(COVpeak(2*NCp(l)+3,2*NCp(l)+3));
			//double err_peak_position=err_M1_peak/peak_weight-err_norm_peak*peak_center/peak_weight;
			//double err_std_peak=(-err_norm_peak*m2_lf_min/peak_weight+2*err_norm_peak*peak_center*peak_center/peak_weight+err_M2_peak/peak_weight-2*err_M1_peak*peak_center/peak_weight)/(2*peak_width);
			
			cout<<"Peak detected\n";
			cout<<"peak width: "<<peak_width<<endl;
			cout<<"peak weight: "<<peak_weight<<endl;
			
			dw_peak=peak_width/2.0;
		}
		else
		{
			cout<<"no peak found\n";
		}
//	}
	
	return peak_exists;
}





bool OmegaMaxEnt_data::set_initial_spectrum()
{
	Nw_lims.zeros(2);
	ws.zeros(2);
	
	bool A0_loaded=false;

	cout<<"initial spectral function provided\n";
	
	w=Aw_data.col(0);
	Nw=w.n_rows;
	A0=Aw_data.col(1);
	A0=A0.rows(1,Nw-2);
	double wmin=w(0), wmax=w(Nw-1);
	int j1=0, j2=1, j3=2;
	w0l=(2*w(j1)*w(j3)-w(j1)*w(j2)-w(j2)*w(j3))/(w(j1)-2*w(j2)+w(j3));
	j1=Nw-3;
	j2=Nw-2;
	j3=Nw-1;
	w0r=(2*w(j1)*w(j3)-w(j1)*w(j2)-w(j2)*w(j3))/(w(j1)-2*w(j2)+w(j3));
	if (w0l<wmin || w0r>wmax)
	{
		cout<<"the real frequency grid in file "<<init_spectr_func_file<<" was not generated by this code\n";
		cout<<"the spectrum in that file cannot be used as the initial one in the calculation\n";
	}
	else
	{
		int j=0;
		double dw1=w(j+1)-w(j);
		double dw2=w(j+2)-w(j+1);
		double rdw=dw2/dw1;
		while (abs(rdw-1.0)>tol_rdw)
		{
			j=j+1;
			dw1=dw2;
			dw2=w(j+2)-w(j+1);
			rdw=dw2/dw1;
		}
		Nw_lims(0)=j+1;
//		int Nu_l=j+1;
		dwl=dw2;
		wl=w(j+1);
		j=Nw-3;
		dw1=w(j+1)-w(j);
		dw2=w(j+2)-w(j+1);
		rdw=dw2/dw1;
		while (abs(rdw-1.0)>tol_rdw)
		{
			j=j-1;
			dw2=dw1;
			dw1=w(j+1)-w(j);
			rdw=dw2/dw1;
		}
		Nw_lims(1)=j+1;
//		int Nu_r=Nw-j-2;
		dwr=dw1;
		wr=w(j+1);
		wc=w.rows(Nw_lims(0),Nw_lims(1));
		A0_loaded=true;
		Nwc=wc.n_rows;
//		if (!SW_set)
//		{
			SW=wr-wl;
			SW_set=true;
//		}
		if (!SC_set)
		{
			SC=(wl+wr)/2.0;
			SC_set=true;
		}
		ws(0)=w0l;
		ws(1)=w0r;
		
		w_exists=true;
		wc_exists=true;
		main_spectral_region_set=true;
		
//		if (!w_origin_in.size()) w_origin=SC;
		if (!w_origin_set)
		{
			w_origin=SC;
			w_origin_set=true;
		}
		j=0;
		while (abs(w_origin-wc(j+1))<abs(w_origin-wc(j)) && j<Nwc-2) j=j+1;
		w_origin=wc(j);
		
		Du_constant=true;
	}

	return A0_loaded;
}









bool OmegaMaxEnt_data::set_initial_spectrum_chi()
{
	Nw_lims.zeros(1);
	ws.zeros(1);
	
	bool A0_loaded=false;
	
	cout<<"initial spectral function provided\n";
	
	w=Aw_data.col(0);
	Nw=w.n_rows;
	A0=Aw_data.col(1);
	A0=A0.rows(0,Nw-2);
	double wmax=w(Nw-1);
	if (w(0)!=0)
	{
		cout<<"set_initial_spectrum_chi(): if \"Im(G) column in data file\" is smaller than 1, the first frequency in initial spectrum file must be equal to 0\n";
		return false;
	}
	
	int j1, j2, j3;
	j1=Nw-3;
	j2=Nw-2;
	j3=Nw-1;
	w0r=(2*w(j1)*w(j3)-w(j1)*w(j2)-w(j2)*w(j3))/(w(j1)-2*w(j2)+w(j3));
	if (w0r>wmax)
	{
		cout<<"the real frequency grid in file "<<init_spectr_func_file<<" was not generated by this code\n";
		cout<<"the spectrum in that file cannot be used as the initial one in the calculation\n";
	}
	else
	{
		int j;
		double dw1, dw2, rdw;
		j=Nw-3;
		dw1=w(j+1)-w(j);
		dw2=w(j+2)-w(j+1);
		rdw=dw2/dw1;
		while (abs(rdw-1.0)>tol_rdw)
		{
			j=j-1;
			dw2=dw1;
			dw1=w(j+1)-w(j);
			rdw=dw2/dw1;
		}
		Nw_lims(0)=j+1;
		dwr=dw1;
		wr=w(j+1);
		wc=w.rows(0,Nw_lims(0));
		A0_loaded=true;
		Nwc=wc.n_rows;
		if (!SW_set)
		{
			SW=2*wr;
			SW_set=true;
		}
		if (!SC_set)
		{
			SC=0;
			SC_set=true;
		}
		ws(0)=w0r;
		
		w_exists=true;
		wc_exists=true;
		main_spectral_region_set=true;
		Du_constant=true;
	}
	
	return A0_loaded;
}






bool OmegaMaxEnt_data::set_grid_omega_from_file()
{
	bool grid_set=false;
	
	int j, jl, jr;
	cout<<"real frequency grid file provided\n";
	
	Nw=grid_w_data.n_rows;
	vec grid_w=grid_w_data.col(0);
	if (Nw<Nw_min)
	{
		cout<<"warning: size of privided real frequency grid must be larger than "<<Nw_min<<endl;
	}
	double wmin=grid_w(0), wmax=grid_w(Nw-1);
	int j1=0, j2=1, j3=2;
	w0l=(2*grid_w(j1)*grid_w(j3)-grid_w(j1)*grid_w(j2)-grid_w(j2)*grid_w(j3))/(grid_w(j1)-2*grid_w(j2)+grid_w(j3));
	j1=Nw-3;
	j2=Nw-2;
	j3=Nw-1;
	w0r=(2*grid_w(j1)*grid_w(j3)-grid_w(j1)*grid_w(j2)-grid_w(j2)*grid_w(j3))/(grid_w(j1)-2*grid_w(j2)+grid_w(j3));
	if (w0l<wmin || w0r>wmax)
	{
		cout<<"The real frequency grid in file "<<grid_omega_file<<" was not generated by this code,\n";
		if (!main_spectral_region_set)
		{
			wl=grid_w(0)+EPSILON;
			wr=grid_w(Nw-1)-EPSILON;
		}
		
		cout<<"only the part in the main spectral region will be used.\n";
		
		if (grid_w(0)<wl && grid_w(Nw-1)>wr)
		{
			j=0;
			while (grid_w(j)<wl && j<Nw-1)
			{
				j=j+1;
			}
			jl=j;
			if ((grid_w(j)-wl)>(wl-grid_w(j-1)))
			{
				jl=j-1;
			}
			wl=grid_w(jl);
			dwl=grid_w(jl+1)-grid_w(jl);
			while (grid_w(j)<wr && j<Nw-1)
			{
				j=j+1;
			}
			jr=j;
			if ((grid_w(j)-wr)>(wr-grid_w(j-1)))
			{
				jr=j-1;
			}
			wr=grid_w(jr);
			main_spectral_region_set=true;
			dwr=grid_w(jr)-grid_w(jr-1);
			wc=grid_w.rows(jl,jr);
			wc_exists=true;
			Nwc=wc.n_rows;
			
//			if (!SW_set)
//			{
				SW=wr-wl;
				SW_set=true;
//			}
			if (!SC_set)
			{
				SC=(wl+wr)/2.0;
				SC_set=true;
			}
			grid_set=true;
//			if (!w_origin_in.size())
			if (!w_origin_set)
			{
				w_origin=SC;
				w_origin_set=true;
			}
			j=0;
			while (abs(w_origin-wc(j+1))<abs(w_origin-wc(j)) && j<Nwc-2) j=j+1;
			w_origin=wc(j);
		}
		else
		{
			cout<<"The provided real frequency grid must extend beyond the spectral function frequency range.\n";
			cout<<"The grid cannot be used.\n";
		}
	}
	else
	{
		Nw_lims.zeros(2);
		ws.zeros(2);
		ws(0)=w0l;
		ws(1)=w0r;
		w=grid_w;
		int j=0;
		double dw1=w(j+1)-w(j);
		double dw2=w(j+2)-w(j+1);
		double rdw=dw2/dw1;
		while (abs(rdw-1.0)>tol_rdw)
		{
			j=j+1;
			dw1=dw2;
			dw2=w(j+2)-w(j+1);
			rdw=dw2/dw1;
		}
		Nw_lims(0)=j+1;
//		int Nu_l=j+1;
		dwl=dw2;
		wl=w(j+1);
		j=Nw-3;
		dw1=w(j+1)-w(j);
		dw2=w(j+2)-w(j+1);
		rdw=dw2/dw1;
		while (abs(rdw-1.0)>tol_rdw)
		{
			j=j-1;
			dw2=dw1;
			dw1=w(j+1)-w(j);
			rdw=dw2/dw1;
		}
		Nw_lims(1)=j+1;
//		int Nu_r=Nw-j-2;
		dwr=dw1;
		wr=w(j+1);
		wc=w.rows(Nw_lims(0),Nw_lims(1));
		Nwc=wc.n_rows;
		main_spectral_region_set=true;
//		if (!SW_set)
//		{
			SW=wr-wl;
			SW_set=true;
//		}
		if (!SC_set)
		{
			SC=(wl+wr)/2.0;
			SC_set=true;
		}
		w_exists=true;
		wc_exists=true;
		grid_set=true;
//		if (!w_origin_in.size()) w_origin=SC;
		if (!w_origin_set)
		{
			w_origin=SC;
			w_origin_set=true;
		}
		j=0;
		while (abs(w_origin-wc(j+1))<abs(w_origin-wc(j)) && j<Nwc-2) j=j+1;
		w_origin=wc(j);
		
		Du_constant=true;
	}
	
	return grid_set;
}










bool OmegaMaxEnt_data::set_grid_omega_from_file_chi()
{
	bool grid_set=false;
	
	int j, jr;
	cout<<"real frequency grid file provided\n";
	
	Nw=grid_w_data.n_rows;
	vec grid_w=grid_w_data.col(0);
	if (Nw<Nw_min)
	{
		cout<<"warning: size of privided real frequency grid must be larger than "<<Nw_min<<endl;
	}
	double wmax=grid_w(Nw-1);
	int j1, j2, j3;
	j1=Nw-3;
	j2=Nw-2;
	j3=Nw-1;
	w0r=(2*grid_w(j1)*grid_w(j3)-grid_w(j1)*grid_w(j2)-grid_w(j2)*grid_w(j3))/(grid_w(j1)-2*grid_w(j2)+grid_w(j3));
	if (w0r>wmax)
	{
		cout<<"The real frequency grid in file "<<grid_omega_file<<" was not generated by this code,\n";
		if (!main_spectral_region_set)
			wr=grid_w(Nw-1)-EPSILON;
		
		cout<<"only the part in the main spectral region will be used.\n";
		
		if (grid_w(0)<=0 && grid_w(Nw-1)>wr)
		{
			j=0;
			while (grid_w(j)<wr && j<Nw-1)  j++;
			jr=j;
			if ((grid_w(j)-wr)>(wr-grid_w(j-1)))  jr=j-1;
			wr=grid_w(jr);
			main_spectral_region_set=true;
			dwr=grid_w(jr)-grid_w(jr-1);
			wc=grid_w.rows(0,jr);
			wc_exists=true;
			Nwc=wc.n_rows;
			
			if (!SW_set)
			{
				SW=2*wr;
				SW_set=true;
			}
			if (!SC_set)
			{
				SC=0;
				SC_set=true;
			}
			grid_set=true;
		}
		else
		{
			cout<<"The provided real frequency grid must extend beyond the spectral function frequency range.\n";
			cout<<"The grid cannot be used.\n";
		}
	}
	else
	{
		Nw_lims.zeros(1);
		ws.zeros(1);
		ws(0)=w0r;
		w=grid_w;
		int j;
		double dw1, dw2, rdw;
		j=Nw-3;
		dw1=w(j+1)-w(j);
		dw2=w(j+2)-w(j+1);
		rdw=dw2/dw1;
		while (abs(rdw-1.0)>tol_rdw)
		{
			j=j-1;
			dw2=dw1;
			dw1=w(j+1)-w(j);
			rdw=dw2/dw1;
		}
		Nw_lims(0)=j+1;
		dwr=dw1;
		wr=w(j+1);
		wc=w.rows(0,Nw_lims(0));
		Nwc=wc.n_rows;
		main_spectral_region_set=true;
		if (!SW_set)
		{
			SW=2*wr;
			SW_set=true;
		}
		if (!SC_set)
		{
			SC=0;
			SC_set=true;
		}
		w_exists=true;
		wc_exists=true;
		grid_set=true;
		Du_constant=true;
	}
	
	return grid_set;
}














bool OmegaMaxEnt_data::set_grid_from_params_chi()
{
	
	int j;
	bool grid_set=false;
	cout<<"grid parameters provided\n";
	uint Npar_grid=omega_grid_params.n_cols;
	
	if (omega_grid_params(0)!=0)
	{
		cout<<"set_grid_from_params_chi(): if \"Im(G) column in data file\" is smaller than 1, the first parameter in grid parameters must be 0\n";
		return false;
	}
	
	if ( Npar_grid % 2 )
	{
		uint Nint=(Npar_grid-1)/2;
		uint Nlims=Nint+1;
		uvec odd_ind=linspace<uvec>(1,Npar_grid-2,Nint);
		rowvec dw_par=omega_grid_params.cols(odd_ind);
		bool Rdw_too_large=false, Rdw_too_small=false;
		for (j=0; j<Nint-1; j++)
		{
			if (dw_par(j+1)/dw_par(j)>Rdw_max) Rdw_too_large=true;
			if (dw_par(j+1)/dw_par(j)<(1.0/Rdw_max)) Rdw_too_small=true;
		}
		if (Rdw_too_large || Rdw_too_small)
		{
			cout<<"The ratio of frequency steps in adjacent intervals of \"grid parameters\" is too large. Spurious oscillations may appear in the spectral function. The maximum ratio recommended is "<<Rdw_max<<endl;
		}
		uvec even_ind=linspace<uvec>(0,Npar_grid-1,Nlims);
		rowvec wlims_par=omega_grid_params.cols(even_ind);
		
		//		cout<<"Dw: "<<dw_par<<endl;
		//		cout<<"wlims: "<<wlims_par<<endl;
		bool ordered=true;
		for (j=1; j<Nlims; j++)
		{
			if (wlims_par(j)<wlims_par(j-1)) ordered=false;
		}
		if (ordered)
		{
			bool grid_par_ok=true;
			bool grid_step_pos=true;
			for (j=1; j<Nlims; j++)
			{
				if (dw_par(j-1)>=(wlims_par(j)-wlims_par(j-1))/Rmin_Dw_dw) grid_par_ok=false;
				if (dw_par(j-1)<0)
				{
					grid_par_ok=false;
					grid_step_pos=false;
				}
			}
			if (grid_par_ok)
			{
				double w0tmp=(wlims_par(0)+wlims_par(Nlims-1))/2.0;
				j=0;
				while (w0tmp>=wlims_par(j+1))
				{
					j=j+1;
				}
				double w0_par=0;
				vec R(2);
				R(0)=RW_grid;
				R(1)=RWD_grid;
				if (non_uniform_frequency_grid(dw_par, wlims_par, w0_par, R, wc))
				{
					wc_exists=true;
					Nwc=wc.n_rows;
					wr=wc(Nwc-1);
					main_spectral_region_set=true;
					dwr=wc(Nwc-1)-wc(Nwc-2);
					if (!SW_set)
					{
						SW=2*wr;
						SW_set=true;
					}
					if (!SC_set)
					{
						SC=0;
						SC_set=true;
					}
					grid_set=true;
				}
				else
				{
					cout<<"parameterized grid definition failed\n";
				}
			}
			else
			{
				if (grid_step_pos)
					cout<<"Error: the frequency step in \"grid parameters\" is too large in at least one interval. The step must be at least "<<Rmin_Dw_dw<<" times smaller than the interval size.\n";
				else
					cout<<"Error: negative step value found in \"grid parameters\"\n";
			}
		}
		else
		{
			cout<<"Error: interval boundaries in \"grid parameters\" are not strictly increasing values\n";
		}
	}
	else
	{
		cout<<"Error: \"grid parameters\" must have an odd number of elements\n";
	}
	
//	cout<<"grid_set: "<<grid_set<<endl;
	
	return grid_set;
}











bool OmegaMaxEnt_data::set_grid_from_params()
{
	int j;
	bool grid_set=false;
	cout<<"grid parameters provided\n";
	uint Npar_grid=omega_grid_params.n_cols;
	if ( Npar_grid % 2 )
	{
		uint Nint=(Npar_grid-1)/2;
		uint Nlims=Nint+1;
		uvec odd_ind=linspace<uvec>(1,Npar_grid-2,Nint);
		rowvec dw_par=omega_grid_params.cols(odd_ind);
		bool Rdw_too_large=false, Rdw_too_small=false;
		for (j=0; j<Nint-1; j++)
		{
			if (dw_par(j+1)/dw_par(j)>Rdw_max) Rdw_too_large=true;
			if (dw_par(j+1)/dw_par(j)<(1.0/Rdw_max)) Rdw_too_small=true;
		}
		if (Rdw_too_large || Rdw_too_small)
		{
			cout<<"The ratio of frequency steps in adjacent intervals of \"grid parameters\" is too large. Spurious oscillations may appear in the spectral function. The maximum ratio recommended is "<<Rdw_max<<endl;
		}
		uvec even_ind=linspace<uvec>(0,Npar_grid-1,Nlims);
		rowvec wlims_par=omega_grid_params.cols(even_ind);
		
//		cout<<"Dw: "<<dw_par<<endl;
//		cout<<"wlims: "<<wlims_par<<endl;
		bool ordered=true;
		for (j=1; j<Nlims; j++)
		{
			if (wlims_par(j)<wlims_par(j-1)) ordered=false;
		}
		if (ordered)
		{
			bool grid_par_ok=true;
			bool grid_step_pos=true;
			for (j=1; j<Nlims; j++)
			{
				if (dw_par(j-1)>=(wlims_par(j)-wlims_par(j-1))/Rmin_Dw_dw) grid_par_ok=false;
				if (dw_par(j-1)<0)
				{
					grid_par_ok=false;
					grid_step_pos=false;
				}
			}
			if (grid_par_ok)
			{
				double w0tmp=(wlims_par(0)+wlims_par(Nlims-1))/2.0;
				j=0;
				while (w0tmp>=wlims_par(j+1))
				{
					j=j+1;
				}
				double w0_par=(wlims_par(j)+wlims_par(j+1))/2.0;
//				if (w_origin_in.size())
				if (w_origin_set)
				{
					w0_par=w_origin;
				}
				else
				{
					w_origin=w0_par;
					w_origin_set=true;
				}
				vec R(2);
				R(0)=RW_grid;
				R(1)=RWD_grid;
				if (non_uniform_frequency_grid(dw_par, wlims_par, w0_par, R, wc))
				{
					wc_exists=true;
					Nwc=wc.n_rows;
					wl=wc(0);
					wr=wc(Nwc-1);
					main_spectral_region_set=true;
					dwl=wc(1)-wc(0);
					dwr=wc(Nwc-1)-wc(Nwc-2);
					if (!SW_set)
					{
						SW=wr-wl;
						SW_set=true;
					}
					if (!SC_set)
					{
						SC=(wl+wr)/2;
						SC_set=true;
					}
					grid_set=true;
				}
				else
				{
					cout<<"parameterized grid definition failed\n";
				}
			}
			else
			{
				if (grid_step_pos)
					cout<<"Error: the frequency step in \"grid parameters\" is too large in at least one interval. The step must be at least "<<Rmin_Dw_dw<<" times smaller than the interval size.\n";
				else
					cout<<"Error: negative step value found in \"grid parameters\"\n";
			}
		}
		else
		{
			cout<<"Error: interval boundaries in \"grid parameters\" are not strictly increasing values\n";
		}
	}
	else
	{
		cout<<"Error: \"grid parameters\" must have an odd number of elements\n";
	}
	
	return grid_set;
}















bool OmegaMaxEnt_data::set_wc_chi()
{
	bool use_nu_grid=false;
	
	if (!step_omega_in.size() && non_uniform_grid && peak_exists)
	{
		double dw_default=0, dw_min;
		
		dw_min=2*dw_peak/R_peak_width_dw;
		
		if (SW_set)
		{
			dw_default=SW/(f_SW_std_omega*Rmin_SW_dw);
		}
		else if (std_omega)
		{
			dw_default=std_omega/Rmin_SW_dw;
		}
		else if (main_spectral_region_set)
		{
			dw_default=2*wr/(f_SW_std_omega*Rmin_SW_dw);
		}
		else if (jfit)
		{
			dw_default=2*wn(jfit)/(R_wncutoff_wr*f_SW_std_omega*Rmin_SW_dw);
		}
		
		if (dw_min<dw_default) use_nu_grid=true;
	}
	
	
	if (step_omega_in.size() || !non_uniform_grid || !peak_exists || !use_nu_grid)
	{
		//cout<<"uniform grid\n";
		
		double dw;
		
		if (SW_set && !main_spectral_region_set)
		{
			wr=SW/2;
			main_spectral_region_set=true;
		}
		else if (std_omega && !main_spectral_region_set)
		{
			SW=f_SW_std_omega*std_omega;
			SW_set=true;
			wr=SW/2;
			main_spectral_region_set=true;
		}
		else if (main_spectral_region_set)
		{
			SW=2*wr;
			SW_set=true;
		}
		else
		{
			cout<<"The central part of the grid cannot be defined. Not enough information available.\n";
			return false;
		}
		
		if (step_omega_in.size())	dw=step_omega;
		else dw=SW/(f_SW_std_omega*Rmin_SW_dw);
		
		if (peak_exists)
		{
			if (dw>dw_peak)
			{
				cout<<"warning: step is larger than the estimated width of the peak at low energy\n";
			}
		}
		
		wl=0;
		Nwc=round(wr/dw);
		wr=Nwc*dw;
		Nwc=Nwc+1;
		wc=linspace<vec>(0,wr,Nwc);
		
		dwr=dw;
		wc_exists=true;
		main_spectral_region_set=true;
		
		return true;
	}
	else
	{
		if (!SW_set && !jfit && !std_omega)
		{
			cout<<"The central part of the grid cannot be defined. Not enough information available.\n";
			return false;
		}
		
		//double R_peak_width_dw=10, R_wncutoff_wr=10, R_Dw_dw=40, R_SW_wr=1, R_wmax_wr_min=3;
		double wr_tmp;
		
		omega_grid_params.zeros(3);
		
		if (SW_in.size())
		{
			wr_tmp=SW/2;
		}
		else if (jfit)
		{
			wr_tmp=wn(jfit)/R_wncutoff_wr;
		}
		else if (SW_set)
		{
			wr_tmp=R_SW_wr*SW/2;
		}
		else
		{
			wr_tmp=R_SW_wr*f_SW_std_omega*std_omega/2;
		}
		
		if (SW_set)
		{
			double wmax_tmp=f_w_range*SW/2.0;
			if (wr_tmp>wmax_tmp/R_wmax_wr_min) wr_tmp=wmax_tmp/R_wmax_wr_min;
		}
		
//		if (peak_exists)
//		{
			omega_grid_params(1)=2*dw_peak/R_peak_width_dw;
//		}
//		else
//		{
//			omega_grid_params(1)=wr_tmp/(f_SW_std_omega*Rmin_SW_dw);
//		}
		
		rowvec params_tmp;
		int j=2;
		omega_grid_params(2)=R_Dw_dw*omega_grid_params(1);
		
		if (SW_set)
		{
			if (omega_grid_params(1)>SW/(f_SW_std_omega*Rmin_SW_dw))	omega_grid_params(1)=SW/(f_SW_std_omega*Rmin_SW_dw);
		}
//		cout<<"wr_tmp: "<<wr_tmp<<endl;
		
//		cout<<omega_grid_params(0)<<endl;
//		cout<<omega_grid_params(1)<<endl;
//		cout<<omega_grid_params(2)<<endl;
		while (omega_grid_params(j)<wr_tmp)
		{
			params_tmp.zeros(j+3);
			params_tmp.cols(0,j)=omega_grid_params;
			omega_grid_params=params_tmp;
			j++;
			omega_grid_params(j)=2*omega_grid_params(j-2);
//			cout<<omega_grid_params(j)<<endl;
			j++;
			omega_grid_params(j)=omega_grid_params(j-2)+R_Dw_dw*omega_grid_params(j-1);
//			cout<<omega_grid_params(j)<<endl;
		}
		
		return set_grid_from_params_chi();
	}
}























bool OmegaMaxEnt_data::Kernel_G_bosons()
{
	//bool use_HF_exp=true;
	double fg=1.7;
	double fi=fg;
	int pnmax=100;
	
	cout<<"defining kernel matrix...\n";
	
	Kcx.zeros(Nn,Nw);
	KM.zeros(3,Nw);
 
	int Nug=Nw_lims(0);
	int Nud=Nw-Nw_lims(1)-1;
 
	vec ud, ug;
	if (Du_constant)
	{
		double dug=dwl/((wl-dwl-w0l)*(wl-w0l));
		vec ug_int=linspace<vec>(1,Nug,Nug);
		ug=-dug*ug_int;
		
		double dud=dwr/((wr-w0r)*(wr+dwr-w0r));
		vec ud_int=linspace<vec>(Nud,1,Nud);
		ud=dud*ud_int;
	}
	else
	{
		ug=1.0/(w.rows(0,Nug-1)-w0l);
		ud=1.0/(w.rows(Nw_lims(1)+1,Nw-1)-w0r);
	}
 
	mat MM;
	spline_matrix_G_part(w, Nw_lims, ws, MM);
	
	int Ncfs0=MM.n_rows;
	
	mat MC(Ncfs0+Nw,Nw);
	
	MC.submat(0,0,Ncfs0-1,Nw-1)=MM;
	MC.submat(Ncfs0,0,Ncfs0+Nw-1,Nw-1)=eye<mat>(Nw,Nw);
	
	int Nint=Nw+1;
	int Nintg=Nug+1;
	
	mat Pa_g=zeros<mat>(Nintg,3*Nint-2+Nw);
	mat Pb_g_r=zeros<mat>(Nintg,3*Nint-2+Nw);
	mat Pc_g_r=zeros<mat>(Nintg,3*Nint-2+Nw);
	mat Pd_g_r=zeros<mat>(Nintg,3*Nint-2+Nw);
	
	Pa_g(0,0)=1;
	Pb_g_r(0,1)=1;
	
	int j;
	for (j=2; j<=Nintg; j++)
	{
		Pa_g(j-1,3*j-4)=1;
		Pb_g_r(j-1,3*j-3)=1;
		Pc_g_r(j-1,3*j-2)=1;
		Pd_g_r(j-1,j+3*Nint-4)=1;
	}
	
	mat U=zeros<mat>(Nug+1,Nug+1);
	vec ug1=zeros<vec>(Nug+1);
	ug1.rows(1,Nug)=ug;
	
	U.diag()=ug1;
	
	mat Pb_g=Pb_g_r-3*U*Pa_g;
	mat Pc_g=Pc_g_r+3*pow(U,2)*Pa_g-2*U*Pb_g_r;
	mat Pd_g=Pd_g_r-pow(U,3)*Pa_g+pow(U,2)*Pb_g_r-U*Pc_g_r;
	
	int	Nintc=Nwc-1;
	
	mat	Pa_c=zeros<mat>(Nintc,3*Nint-2+Nw);
	mat	Pb_c_r=zeros<mat>(Nintc,3*Nint-2+Nw);
	mat	Pc_c_r=zeros<mat>(Nintc,3*Nint-2+Nw);
	mat	Pd_c_r=zeros<mat>(Nintc,3*Nint-2+Nw);
	
	for (j=1; j<=Nintc; j++)
	{
		Pa_c(j-1,3*j-4+3*Nintg)=1;
		Pb_c_r(j-1,3*j-3+3*Nintg)=1;
		Pc_c_r(j-1,3*j-2+3*Nintg)=1;
		Pd_c_r(j-1,j+3*Nint-3+Nug)=1;
	}
	
	mat W=diagmat(wc.rows(0,Nwc-2));
	
	mat Pb_c=Pb_c_r-3*W*Pa_c;
	mat Pc_c=Pc_c_r+3*pow(W,2)*Pa_c-2*W*Pb_c_r;
	mat Pd_c=Pd_c_r-pow(W,3)*Pa_c+pow(W,2)*Pb_c_r-W*Pc_c_r;
	
	int Nintd=Nud+1;
	
	mat Pa_d=zeros<mat>(Nintd,3*Nint-2+Nw);
	mat Pb_d_r=zeros<mat>(Nintd,3*Nint-2+Nw);
	mat Pc_d_r=zeros<mat>(Nintd,3*Nint-2+Nw);
	mat Pd_d_r=zeros<mat>(Nintd,3*Nint-2+Nw);
	
	for (j=1; j<Nintd; j++)
	{
		Pa_d(j-1,3*j-4+3*Nintg+3*Nintc)=1;
		Pb_d_r(j-1,3*j-3+3*Nintg+3*Nintc)=1;
		Pc_d_r(j-1,3*j-2+3*Nintg+3*Nintc)=1;
		Pd_d_r(j-1,j+3*Nint-3+Nug+Nwc)=1;
	}
	
	j=Nintd;
	Pa_d(j-1,3*j-4+3*Nintg+3*Nintc)=1;
	Pb_d_r(j-1,3*j-3+3*Nintg+3*Nintc)=1;
	
	U.zeros(Nud+1,Nud+1);
	vec ud1=zeros<vec>(Nud+1);
	ud1.rows(0,Nud-1)=ud;
	U.diag()=ud1;
	
	mat Pb_d=Pb_d_r-3*U*Pa_d;
	mat Pc_d=Pc_d_r+3*pow(U,2)*Pa_d-2*U*Pb_d_r;
	mat Pd_d=Pd_d_r-pow(U,3)*Pa_d+pow(U,2)*Pb_d_r-U*Pc_d_r;

	cx_mat Ka_g=zeros<cx_mat>(Nn,Nintg);
 	cx_mat Kb_g=zeros<cx_mat>(Nn,Nintg);
 	cx_mat Kc_g=zeros<cx_mat>(Nn,Nintg);
	cx_mat Kd_g=zeros<cx_mat>(Nn,Nintg);
 
 	rowvec ug2(Nug+1);
	ug2.cols(0,Nug-1)=ug.t();
	ug2(Nug)=1.0/(wl-w0l);
	
	mat Wng=wn.rows(1,Nn-1)*ones<rowvec>(Nintg-1);
	mat Ug=ones<vec>(Nn-1)*ug2;
	
	dcomplex i(0,1);
 
	rowvec wg=1/ug2+w0l;
	mat Wg=ones<vec>(Nn-1)*wg;
	
	cx_rowvec cx_vec_tmp(Nug);
	cx_vec_tmp.zeros();

	rowvec vec_tmp=-(pow(ug2.cols(1,Nintg-1),2)-pow(ug2.cols(0,Nintg-2),2))/(4*PI);
	cx_vec_tmp.set_real(vec_tmp);
	Ka_g.submat(0,1,0,Nintg-1)=cx_vec_tmp;
	vec_tmp=-(ug2.cols(1,Nintg-1)-ug2.cols(0,Nintg-2))/(2*PI);
	cx_vec_tmp.set_real(vec_tmp);
	Kb_g.submat(0,1,0,Nintg-1)=cx_vec_tmp;
	vec_tmp=-log(ug2.cols(1,Nintg-1)/ug2.cols(0,Nintg-2))/(2*PI);
	cx_vec_tmp.set_real(vec_tmp);
	Kc_g.submat(0,1,0,Nintg-1)=cx_vec_tmp;
	vec_tmp=trans(w.rows(1,Nintg-1)-w.rows(0,Nintg-2))/(2*PI);
	cx_vec_tmp.set_real(vec_tmp);
	Kd_g.submat(0,1,0,Nintg-1)=cx_vec_tmp;
	
// 	Ka_g.submat(0,1,0,Nintg-1)=-(pow(ug2.cols(1,Nintg-1),2)-pow(ug2.cols(0,Nintg-2),2))/(4*PI);
// 	Kb_g.submat(0,1,0,Nintg-1)=-(ug2.cols(1,Nintg-1)-ug2.cols(0,Nintg-2))/(2*PI);
// 	Kc_g.submat(0,1,0,Nintg-1)=-log(ug2.cols(1,Nintg-1)/ug2.cols(0,Nintg-2))/(2*PI);
// 	Kd_g.submat(0,1,0,Nintg-1)=(w.rows(1,Nintg-1)-w.rows(0,Nintg-2))/(2*PI);


 	mat atang=atan((Wng % (Ug.cols(1,Nintg-1)-Ug.cols(0,Nintg-2)))/(1+w0l*(Ug.cols(1,Nintg-1)+Ug.cols(0,Nintg-2))+(pow(w0l,2)+pow(Wng,2)) % Ug.cols(0,Nintg-2) % Ug.cols(1,Nintg-1)));
	mat logg=log((1+2*w0l*Ug.cols(1,Nintg-1)+pow(Ug.cols(1,Nintg-1),2) % (pow(Wng,2)+pow(w0l,2)))/(1+2*w0l*Ug.cols(0,Nintg-2)+pow(Ug.cols(0,Nintg-2),2) % (pow(Wng,2)+pow(w0l,2))));
	
 	Ka_g.submat(1,1,Nn-1,Nintg-1)=(-i*(2*Wng % (Ug.cols(1,Nintg-1)-Ug.cols(0,Nintg-2)))/pow(Wng + i*w0l,2) - i*w0l*(pow(Ug.cols(1,Nintg-1),2)-pow(Ug.cols(0,Nintg-2),2))/(Wng + i*w0l) +  i*(2*Wng % atang)/pow(Wng + i*w0l,3) - (Wng % logg)/pow(Wng + i*w0l,3) )/(4*PI);
 	Kb_g.submat(1,1,Nn-1,Nintg-1)=(-i*2.0*w0l*(Ug.cols(1,Nintg-1)-Ug.cols(0,Nintg-2))/(Wng + i*w0l) -(2*Wng % atang)/pow(Wng + i*w0l,2) -i*(Wng % logg)/pow(Wng + i*w0l,2) )/(4*PI);
 	Kc_g.submat(1,1,Nn-1,Nintg-1)=( -2*log(Ug.cols(1,Nintg-1)/Ug.cols(0,Nintg-2)) - i*(2*Wng % atang)/(Wng + i*w0l) + (Wng % logg)/(Wng + i*w0l) )/(4*PI);
 	Kd_g.submat(1,1,Nn-1,Nintg-1)=( 2*(Wg.cols(1,Nintg-1)-Wg.cols(0,Nintg-2)) -i*2.0*Wng % log(Ug.cols(1,Nintg-1)/Ug.cols(0,Nintg-2)) + 2*Wng % atang + i*Wng % logg )/(4*PI);
	

	cx_mat Ka_c=zeros<cx_mat>(Nn,Nintc);
 	cx_mat Kb_c=zeros<cx_mat>(Nn,Nintc);
 	cx_mat Kc_c=zeros<cx_mat>(Nn,Nintc);
 	cx_mat Kd_c=zeros<cx_mat>(Nn,Nintc);
 
	cx_vec_tmp.zeros(Nintc);
	
	vec_tmp=-trans(pow(wc.rows(1,Nwc-1),4)-pow(wc.rows(0,Nwc-2),4))/(8*PI);
	cx_vec_tmp.set_real(vec_tmp);
	Ka_c.row(0)=cx_vec_tmp;
	vec_tmp=-trans(pow(wc.rows(1,Nwc-1),3)-pow(wc.rows(0,Nwc-2),3))/(6*PI);
	cx_vec_tmp.set_real(vec_tmp);
	Kb_c.row(0)=cx_vec_tmp;
	vec_tmp=-trans(pow(wc.rows(1,Nwc-1),2)-pow(wc.rows(0,Nwc-2),2))/(4*PI);
	cx_vec_tmp.set_real(vec_tmp);
	Kc_c.row(0)=cx_vec_tmp;
	vec_tmp=-trans(wc.rows(1,Nwc-1)-wc.rows(0,Nwc-2))/(2*PI);
	cx_vec_tmp.set_real(vec_tmp);
	Kd_c.row(0)=cx_vec_tmp;
	
//	Ka_c.row(0)=-(pow(wc.rows(1,Nwc-1),4)-pow(wc.rows(0,Nwc-2),4))/(8*PI);
//	Kb_c.row(0)=-(pow(wc.rows(1,Nwc-1),3)-pow(wc.rows(0,Nwc-2),3))/(6*PI);
//	Kc_c.row(0)=-(pow(wc.rows(1,Nwc-1),2)-pow(wc.rows(0,Nwc-2),2))/(4*PI);
//	Kd_c.row(0)=-(wc.rows(1,Nwc-1)-wc.rows(0,Nwc-2))/(2*PI);
 
	mat Wnc=wn.rows(1,Nn-1)*ones<rowvec>(Nintc);
	mat Wc=ones<vec>(Nn-1)*wc.t();
 
 	mat logc=log((pow(Wnc,2)+pow(Wc.cols(1,Nintc),2))/(pow(Wnc,2)+pow(Wc.cols(0,Nintc-1),2)));
 	mat atanc=atan((Wnc % (Wc.cols(0,Nintc-1)-Wc.cols(1,Nintc)))/(Wc.cols(1,Nintc) % Wc.cols(0,Nintc-1)+pow(Wnc,2)));
	
	mat dWc=Wc.cols(1,Nintc)-Wc.cols(0,Nintc-1);
	mat dWc2=pow(Wc.cols(1,Nintc),2)-pow(Wc.cols(0,Nintc-1),2);
	mat dWc3=pow(Wc.cols(1,Nintc),3)-pow(Wc.cols(0,Nintc-1),3);
	mat Wnc2=pow(Wnc,2);
	mat Wnc3=pow(Wnc,3);
	mat Wnc4=pow(Wnc,4);
 
 	Ka_c.rows(1,Nn-1)=-i*( -Wnc3 % dWc + i*(Wnc2 % dWc2)/2 + (Wnc % dWc3)/3 -i*(pow(Wc.cols(1,Nintc),4)-pow(Wc.cols(0,Nintc-1),4))/4 - Wnc4 % atanc - i*Wnc4 % logc/2 )/(2*PI);
	Kb_c.rows(1,Nn-1)=-i*( i*Wnc2 % dWc + (Wnc % dWc2)/2 - i*dWc3/3 + i*Wnc3 % atanc - Wnc3 % logc/2 )/(2*PI);
	Kc_c.rows(1,Nn-1)=-i*( Wnc % dWc - i*dWc2/2 + Wnc2 % atanc + i*Wnc2 % logc/2 )/(2*PI);
 	Kd_c.rows(1,Nn-1)=-i*( -i*dWc -i*Wnc % atanc + Wnc % logc/2 )/(2*PI);

 
 	int Pmax=pnmax+4;
	imat MP;
	pascal(Pmax+1,MP);
 
	double wtmp, wj, dw1;
	vec dwp;
	int p,l, jni;
	for (j=0; j<Nwc-1; j++)
	{
		wtmp=abs(wc(j+1));
		if (abs(wc(j))>wtmp) wtmp=abs(wc(j));
 		jni=1;
		while (wn(jni)<fi*wtmp && jni<Nn)	jni++;
		while (pow(wn(jni),pnmax)==0 && jni<Nn)		jni++;
		
		if (jni<Nn-1)
		{
			wj=wc(j);
			dw1=wc(j+1)-wj;
			dwp.zeros(Pmax);
			dwp(0)=dw1;
 			dwp(1)=pow(dw1,2)+2*wj*dw1;
			for (p=3; p<=Pmax; p++)
			{
				dwp(p-1)=0;
				for (l=0; l<p; l++)
				{
					dwp(p-1)=dwp(p-1)+MP(l,p-l)*pow(dw1,p-l)*pow(wj,l);
				}
			}
 			if (jni<Nn-1)
 			{
	 			Ka_c.submat(jni,j,Nn-1,j)=zeros<cx_vec>(Nn-jni);
 				Kb_c.submat(jni,j,Nn-1,j)=zeros<cx_vec>(Nn-jni);
 				Kc_c.submat(jni,j,Nn-1,j)=zeros<cx_vec>(Nn-jni);
 				Kd_c.submat(jni,j,Nn-1,j)=zeros<cx_vec>(Nn-jni);
				for (p=pnmax; p>=1; p--)
 				{
	 				Ka_c.submat(jni,j,Nn-1,j)=Ka_c.submat(jni,j,Nn-1,j) + pow(i,-p)*dwp(p+3)/(2*PI*(p+4)*pow(wn.rows(jni,Nn-1),p));
 					Kb_c.submat(jni,j,Nn-1,j)=Kb_c.submat(jni,j,Nn-1,j) + pow(i,-p)*dwp(p+2)/(2*PI*(p+3)*pow(wn.rows(jni,Nn-1),p));
 					Kc_c.submat(jni,j,Nn-1,j)=Kc_c.submat(jni,j,Nn-1,j) + pow(i,-p)*dwp(p+1)/(2*PI*(p+2)*pow(wn.rows(jni,Nn-1),p));
 					Kd_c.submat(jni,j,Nn-1,j)=Kd_c.submat(jni,j,Nn-1,j) + pow(i,-p)*dwp(p)/(2*PI*(p+1)*pow(wn.rows(jni,Nn-1),p));
 				}
 			}
		}
	}


	cx_mat Ka_d=zeros<cx_mat>(Nn,Nintd);
 	cx_mat Kb_d=zeros<cx_mat>(Nn,Nintd);
 	cx_mat Kc_d=zeros<cx_mat>(Nn,Nintd);
 	cx_mat Kd_d=zeros<cx_mat>(Nn,Nintd);
	
	rowvec ud2(Nud+1);
	ud2.cols(1,Nud)=ud.t();
	ud2(0)=1.0/(wr-w0r);
	
	mat Wnd=wn.rows(1,Nn-1)*ones<rowvec>(Nintd-1);
	mat Ud=ones<vec>(Nn-1)*ud2;
	
	rowvec wd=1/ud2+w0r;
	mat Wd=ones<vec>(Nn-1)*wd;
	
	cx_vec_tmp.zeros(Nud);


	vec_tmp=-(pow(ud2.cols(1,Nintd-1),2)-pow(ud2.cols(0,Nintd-2),2))/(4*PI);
	cx_vec_tmp.set_real(vec_tmp);
	Ka_d.submat(0,0,0,Nintd-2)=cx_vec_tmp;
	vec_tmp=-(ud2.cols(1,Nintd-1)-ud2.cols(0,Nintd-2))/(2*PI);
	cx_vec_tmp.set_real(vec_tmp);
	Kb_d.submat(0,0,0,Nintd-2)=cx_vec_tmp;
	vec_tmp=-log(ud2.cols(1,Nintd-1)/ud2.cols(0,Nintd-2))/(2*PI);
	cx_vec_tmp.set_real(vec_tmp);
	Kc_d.submat(0,0,0,Nintd-2)=cx_vec_tmp;
	vec_tmp=(wd.cols(1,Nintd-1)-wd.cols(0,Nintd-2))/(2*PI);
	cx_vec_tmp.set_real(vec_tmp);
	Kd_d.submat(0,0,0,Nintd-2)=cx_vec_tmp;
	
//	Ka_d.submat(0,0,0,Nintd-2)=-(pow(ud2.cols(1,Nintd-1),2)-pow(ud2.cols(0,Nintd-2),2))/(4*PI);
// 	Kb_d.submat(0,0,0,Nintd-2)=-(ud2.cols(1,Nintd-1)-ud2.cols(0,Nintd-2))/(2*PI);
// 	Kc_d.submat(0,0,0,Nintd-2)=-log(ud2.cols(1,Nintd-1)/ud2.cols(0,Nintd-2))/(2*PI);
//	Kd_d.submat(0,0,0,Nintd-2)=(wd.rows(1,Nintd-1)-wd.rows(0,Nintd-2))/(2*PI);


	mat atand=atan((Wnd % (Ud.cols(1,Nintd-1)-Ud.cols(0,Nintd-2)))/(1+w0r*(Ud.cols(1,Nintd-1)+Ud.cols(0,Nintd-2))+(pow(w0r,2)+pow(Wnd,2)) % Ud.cols(0,Nintd-2) % Ud.cols(1,Nintd-1)));
 	mat logd=log((1+2*w0r*Ud.cols(1,Nintd-1)+pow(Ud.cols(1,Nintd-1),2) % (pow(Wnd,2)+pow(w0r,2)))/(1+2*w0r*Ud.cols(0,Nintd-2)+pow(Ud.cols(0,Nintd-2),2) % (pow(Wnd,2)+pow(w0r,2))));
	
	Ka_d.submat(1,0,Nn-1,Nintd-2)=(-i*(2*Wnd % (Ud.cols(1,Nintd-1)-Ud.cols(0,Nintd-2)))/pow(Wnd + i*w0r,2) - i*(w0r*(pow(Ud.cols(1,Nintd-1),2)-pow(Ud.cols(0,Nintd-2),2)))/(Wnd + i*w0r) +  i*(2*Wnd % atand)/pow(Wnd + i*w0r,3)  - (Wnd % logd)/pow(Wnd + i*w0r,3) )/(4*PI);
 	Kb_d.submat(1,0,Nn-1,Nintd-2)=(-i*(2*w0r*(Ud.cols(1,Nintd-1)-Ud.cols(0,Nintd-2)))/(Wnd + i*w0r) -(2*Wnd % atand)/pow(Wnd + i*w0r,2) -i*(Wnd % logd)/pow(Wnd + i*w0r,2) )/(4*PI);
 	Kc_d.submat(1,0,Nn-1,Nintd-2)=( -2*log(Ud.cols(1,Nintd-1)/Ud.cols(0,Nintd-2)) -i*(2*Wnd % atand)/(Wnd + i*w0r) + (Wnd % logd)/(Wnd + i*w0r) )/(4*PI);
 	Kd_d.submat(1,0,Nn-1,Nintd-2)=( 2*(Wd.cols(1,Nintd-1)-Wd.cols(0,Nintd-2)) -i*2.0*Wnd % log(Ud.cols(1,Nintd-1)/Ud.cols(0,Nintd-2)) + 2*Wnd % atand +i*Wnd % logd )/(4*PI);
	
	
	cx_mat KG=-(Ka_g*Pa_g+Kb_g*Pb_g+Kc_g*Pc_g+Kd_g*Pd_g)*MC;
 	cx_mat KC=(Ka_c*Pa_c+Kb_c*Pb_c+Kc_c*Pc_c+Kd_c*Pd_c)*MC;
 	cx_mat KD=-(Ka_d*Pa_d+Kb_d*Pb_d+Kc_d*Pc_d+Kd_d*Pd_d)*MC;
 
 	Kcx=KG+KC+KD;
	

	rowvec Knorm_a_g=zeros<rowvec>(Nintg);
	rowvec Knorm_b_g=zeros<rowvec>(Nintg);
	rowvec Knorm_c_g=zeros<rowvec>(Nintg);
	rowvec Knorm_d_g=zeros<rowvec>(Nintg);
	
	Knorm_a_g(0)=-pow(ug(0),2)/2;
 	Knorm_b_g(0)=-ug(0);
	
	Knorm_a_g.cols(1,Nintg-1)=-(pow(ug2.cols(1,Nug),2)-pow(ug2.cols(0,Nug-1),2))/2;
 	Knorm_b_g.cols(1,Nintg-1)=-(ug2.cols(1,Nug)-ug2.cols(0,Nug-1));
 	Knorm_c_g.cols(1,Nintg-1)=-log(ug2.cols(1,Nug)/ug2.cols(0,Nug-1));
 	Knorm_d_g.cols(1,Nintg-1)=1.0/ug2.cols(1,Nug)-1.0/ug2.cols(0,Nug-1);
	
	rowvec KM0g=(Knorm_a_g*Pa_g+Knorm_b_g*Pb_g+Knorm_c_g*Pc_g+Knorm_d_g*Pd_g)*MC/(2*PI);
	
	rowvec Knorm_a_c_r=zeros<rowvec>(Nintc);
	rowvec Knorm_b_c_r=zeros<rowvec>(Nintc);
	rowvec Knorm_c_c_r=zeros<rowvec>(Nintc);
	rowvec Knorm_d_c_r=zeros<rowvec>(Nintc);
	
	Knorm_a_c_r=trans(pow(wc.rows(1,Nwc-1)-wc.rows(0,Nwc-2),4)/4);
	Knorm_b_c_r=trans(pow(wc.rows(1,Nwc-1)-wc.rows(0,Nwc-2),3)/3);
	Knorm_c_c_r=trans(pow(wc.rows(1,Nwc-1)-wc.rows(0,Nwc-2),2)/2);
	Knorm_d_c_r=trans(wc.rows(1,Nwc-1)-wc.rows(0,Nwc-2));
	
	rowvec KM0c=(Knorm_a_c_r*Pa_c+Knorm_b_c_r*Pb_c_r+Knorm_c_c_r*Pc_c_r+Knorm_d_c_r*Pd_c_r)*MC/(2*PI);
	
	rowvec Knorm_a_d=zeros<rowvec>(Nintd);
	rowvec Knorm_b_d=zeros<rowvec>(Nintd);
	rowvec Knorm_c_d=zeros<rowvec>(Nintd);
	rowvec Knorm_d_d=zeros<rowvec>(Nintd);
	
	Knorm_a_d.cols(0,Nintd-2)=-(pow(ud2.cols(1,Nud),2)-pow(ud2.cols(0,Nud-1),2))/2;
	Knorm_b_d.cols(0,Nintd-2)=-(ud2.cols(1,Nud)-ud2.cols(0,Nud-1));
	Knorm_c_d.cols(0,Nintd-2)=-log(ud2.cols(1,Nud)/ud2.cols(0,Nud-1));
	Knorm_d_d.cols(0,Nintd-2)=1.0/ud2.cols(1,Nud)-1.0/ud2.cols(0,Nud-1);
	
	Knorm_a_d(Nintd-1)=pow(ud2(Nud-1),2)/2;
	Knorm_b_d(Nintd-1)=ud2(Nud-1);
	
	rowvec KM0d=(Knorm_a_d*Pa_d+Knorm_b_d*Pb_d+Knorm_c_d*Pc_d+Knorm_d_d*Pd_d)*MC/(2*PI);
	
	rowvec KM0=KM0g+KM0c+KM0d;
	

	rowvec KM1_a_g=zeros<rowvec>(Nintg);
	rowvec KM1_b_g=zeros<rowvec>(Nintg);
	rowvec KM1_c_g=zeros<rowvec>(Nintg);
	rowvec KM1_d_g=zeros<rowvec>(Nintg);
	
	KM1_a_g.cols(1,Nintg-1)=Knorm_b_g.cols(1,Nintg-1);
	KM1_b_g.cols(1,Nintg-1)=Knorm_c_g.cols(1,Nintg-1);
	KM1_c_g.cols(1,Nintg-1)=Knorm_d_g.cols(1,Nintg-1);
	KM1_d_g.cols(1,Nintg-1)=(1.0/pow(ug2.cols(1,Nug),2)-1.0/pow(ug2.cols(0,Nug-1),2))/2;
	
	rowvec KM1g_tmp=(KM1_a_g*Pa_g+KM1_b_g*Pb_g+KM1_c_g*Pc_g+KM1_d_g*Pd_g)*MC/(2*PI);
	
	rowvec KM1g=w0l*KM0g+KM1g_tmp;
	
	mat Wjc=diagmat(wc.rows(0,Nintc-1));
	
	rowvec KM1_a_c_tmp=trans(pow(wc.rows(1,Nwc-1)-wc.rows(0,Nwc-2),5)/5.0);
	
	rowvec KM1_a_c=KM1_a_c_tmp+Knorm_a_c_r*Wjc;
	rowvec KM1_b_c=Knorm_a_c_r+Knorm_b_c_r*Wjc;
	rowvec KM1_c_c=Knorm_b_c_r+Knorm_c_c_r*Wjc;
	rowvec KM1_d_c=Knorm_c_c_r+Knorm_d_c_r*Wjc;
	
	rowvec KM1c=(KM1_a_c*Pa_c+KM1_b_c*Pb_c_r+KM1_c_c*Pc_c_r+KM1_d_c*Pd_c_r)*MC/(2*PI);
	
	rowvec KM1_a_d=zeros<rowvec>(Nintd);
	rowvec KM1_b_d=zeros<rowvec>(Nintd);
	rowvec KM1_c_d=zeros<rowvec>(Nintd);
	rowvec KM1_d_d=zeros<rowvec>(Nintd);
	
	KM1_a_d.cols(0,Nintd-2)=Knorm_b_d.cols(0,Nintd-2);
	KM1_b_d.cols(0,Nintd-2)=Knorm_c_d.cols(0,Nintd-2);
	KM1_c_d.cols(0,Nintd-2)=Knorm_d_d.cols(0,Nintd-2);
	KM1_d_d.cols(0,Nintd-2)=(1.0/pow(ud2.cols(1,Nud),2)-1.0/pow(ud2.cols(0,Nud-1),2))/2;
	
	rowvec KM1d_tmp=(KM1_a_d*Pa_d+KM1_b_d*Pb_d+KM1_c_d*Pc_d+KM1_d_d*Pd_d)*MC/(2*PI);
	
	rowvec KM1d=w0r*KM0d+KM1d_tmp;
	
	rowvec KM1=KM1g+KM1c+KM1d;
	
	
	rowvec KM2_a_g=zeros<rowvec>(Nintg);
	rowvec KM2_b_g=zeros<rowvec>(Nintg);
	rowvec KM2_c_g=zeros<rowvec>(Nintg);
	rowvec KM2_d_g=zeros<rowvec>(Nintg);
	
	KM2_a_g.cols(1,Nintg-1)=Knorm_c_g.cols(1,Nintg-1);
	KM2_b_g.cols(1,Nintg-1)=Knorm_d_g.cols(1,Nintg-1);
	KM2_c_g.cols(1,Nintg-1)=KM1_d_g.cols(1,Nintg-1);
	KM2_d_g.cols(1,Nintg-1)=(1.0/pow(ug2.cols(1,Nug),3)-1.0/pow(ug2.cols(0,Nug-1),3))/3;
	
	rowvec KM2g_tmp=(KM2_a_g*Pa_g+KM2_b_g*Pb_g+KM2_c_g*Pc_g+KM2_d_g*Pd_g)*MC/(2*PI);
	
	rowvec KM2g=pow(w0l,2)*KM0g+2*w0l*KM1g_tmp+KM2g_tmp;
	
	rowvec KM2_a_c_tmp=trans(pow(wc.rows(1,Nwc-1)-wc.rows(0,Nwc-2),6)/6);
	
	rowvec KM2_a_c=KM2_a_c_tmp+2*KM1_a_c_tmp*Wjc+Knorm_a_c_r*pow(Wjc,2);
	rowvec KM2_b_c=KM1_a_c_tmp+2*Knorm_a_c_r*Wjc+Knorm_b_c_r*pow(Wjc,2);
	rowvec KM2_c_c=Knorm_a_c_r+2*Knorm_b_c_r*Wjc+Knorm_c_c_r*pow(Wjc,2);
	rowvec KM2_d_c=Knorm_b_c_r+2*Knorm_c_c_r*Wjc+Knorm_d_c_r*pow(Wjc,2);
	
	rowvec KM2c=(KM2_a_c*Pa_c+KM2_b_c*Pb_c_r+KM2_c_c*Pc_c_r+KM2_d_c*Pd_c_r)*MC/(2*PI);
	
	rowvec KM2_a_d=zeros<rowvec>(Nintd);
	rowvec KM2_b_d=zeros<rowvec>(Nintd);
	rowvec KM2_c_d=zeros<rowvec>(Nintd);
	rowvec KM2_d_d=zeros<rowvec>(Nintd);
	
	KM2_a_d.cols(0,Nintd-2)=Knorm_c_d.cols(0,Nintd-2);
	KM2_b_d.cols(0,Nintd-2)=Knorm_d_d.cols(0,Nintd-2);
	KM2_c_d.cols(0,Nintd-2)=KM1_d_d.cols(0,Nintd-2);
	KM2_d_d.cols(0,Nintd-2)=(1.0/pow(ud2.cols(1,Nud),3)-1.0/pow(ud2.cols(0,Nud-1),3))/3;
	
	rowvec KM2d_tmp=(KM2_a_d*Pa_d+KM2_b_d*Pb_d+KM2_c_d*Pc_d+KM2_d_d*Pd_d)*MC/(2*PI);
	
	rowvec KM2d=pow(w0r,2)*KM0d+2*w0r*KM1d_tmp+KM2d_tmp;
	
	rowvec KM2=KM2g+KM2c+KM2d;
	
	rowvec KM3_a_g=zeros<rowvec>(Nintg);
	rowvec KM3_b_g=zeros<rowvec>(Nintg);
	rowvec KM3_c_g=zeros<rowvec>(Nintg);
	rowvec KM3_d_g=zeros<rowvec>(Nintg);
	
	KM3_a_g.cols(1,Nintg-1)=Knorm_d_g.cols(1,Nintg-1);
	KM3_b_g.cols(1,Nintg-1)=KM1_d_g.cols(1,Nintg-1);
	KM3_c_g.cols(1,Nintg-1)=KM2_d_g.cols(1,Nintg-1);
	KM3_d_g.cols(1,Nintg-1)=(1.0/pow(ug2.cols(1,Nug),4)-1.0/pow(ug2.cols(0,Nug-1),4))/4;
	
	rowvec KM3g_tmp=(KM3_a_g*Pa_g+KM3_b_g*Pb_g+KM3_c_g*Pc_g+KM3_d_g*Pd_g)*MC/(2*PI);
	
	rowvec KM3g=pow(w0l,3)*KM0g+3*pow(w0l,2)*KM1g_tmp+3*w0l*KM2g_tmp+KM3g_tmp;
	
	rowvec KM3_a_c_tmp=trans(pow(wc.rows(1,Nwc-1)-wc.rows(0,Nwc-2),7)/7);
	
	rowvec KM3_a_c=KM3_a_c_tmp+3*KM2_a_c_tmp*Wjc+3*KM1_a_c_tmp*pow(Wjc,2)+Knorm_a_c_r*pow(Wjc,3);
	rowvec KM3_b_c=KM2_a_c_tmp+3*KM1_a_c_tmp*Wjc+3*Knorm_a_c_r*pow(Wjc,2)+Knorm_b_c_r*pow(Wjc,3);
	rowvec KM3_c_c=KM1_a_c_tmp+3*Knorm_a_c_r*Wjc+3*Knorm_b_c_r*pow(Wjc,2)+Knorm_c_c_r*pow(Wjc,3);
	rowvec KM3_d_c=Knorm_a_c_r+3*Knorm_b_c_r*Wjc+3*Knorm_c_c_r*pow(Wjc,2)+Knorm_d_c_r*pow(Wjc,3);
	
	rowvec KM3c=(KM3_a_c*Pa_c+KM3_b_c*Pb_c_r+KM3_c_c*Pc_c_r+KM3_d_c*Pd_c_r)*MC/(2*PI);
	
	rowvec KM3_a_d=zeros<rowvec>(Nintd);
	rowvec KM3_b_d=zeros<rowvec>(Nintd);
	rowvec KM3_c_d=zeros<rowvec>(Nintd);
	rowvec KM3_d_d=zeros<rowvec>(Nintd);
	
	KM3_a_d.cols(0,Nintd-2)=Knorm_d_d.cols(0,Nintd-2);
	KM3_b_d.cols(0,Nintd-2)=KM1_d_d.cols(0,Nintd-2);
	KM3_c_d.cols(0,Nintd-2)=KM2_d_d.cols(0,Nintd-2);
	KM3_d_d.cols(0,Nintd-2)=(1.0/pow(ud2.cols(1,Nud),4)-1.0/pow(ud2.cols(0,Nud-1),4))/4;
	
	rowvec KM3d_tmp=(KM3_a_d*Pa_d+KM3_b_d*Pb_d+KM3_c_d*Pc_d+KM3_d_d*Pd_d)*MC/(2*PI);
	
	rowvec KM3d=pow(w0r,3)*KM0d+3*pow(w0r,2)*KM1d_tmp+3*w0r*KM2d_tmp+KM3d_tmp;
	
	rowvec KM3=KM3g+KM3c+KM3d;
	
	KM.row(0)=KM1;
	KM.row(1)=KM2;
	KM.row(2)=KM3;
	
	K.zeros(2*Nn,Nw);
	uvec even_ind=linspace<uvec>(0,2*Nn-2,Nn);
	K.rows(even_ind)=real(Kcx);
	K.rows(even_ind+1)=imag(Kcx);

   ///////////////////////////////
	
	return true;
}










bool OmegaMaxEnt_data::set_covar_chi_omega_n()
{
	//int j;
	
	COV.zeros(Nn,Nn);
	cov_diag=true;
	if (error_file.size())
	{
		if (error_data.n_rows<Nn)
		{
			cout<<"number of lines is too small in file "<<error_file<<endl;
			return false;
		}
		cout<<"error file provided\n";
		errGr=error_data.col(col_errGr-1);
		if (wn_sign_change)
		{
			errGr=errGr.rows(indG_0,indG_f);
		}
		if (wn_inversed)
		{
			errGr=flipud(errGr);
		}
		errG=errGr;
		COV.diag()=square(errG);
	}
	else if ( covar_re_re_file.size() )
	{
		if (CRR.n_rows<Nn || CRR.n_cols<Nn)
		{
			cout<<"number of lines and/or columns in covariance file(s) is too small\n";
			return false;
		}
		cout<<"covariance matrix provided\n";
		cov_diag=false;
		COV=CRR.submat(0,0,Nn-1,Nn-1);
	}
	else
	{
		cout<<"no errors provided\nusing a constant error\n";
		
		double Gr_max=max(abs(Gr));
		errGr=default_error_G*Gr_max*ones<vec>(Nn);
		errG=errGr;
		COV.diag()=square(errG);
	}
	
	COV=0.5*(COV+COV.t());
	
	return true;
}












bool OmegaMaxEnt_data::compute_moments_chi_omega_n()
{
	int j, NC=3;
	
	printf("COMPUTING MOMENTS\n");
	
	vec C2b(Nn-2), C4b(Nn-2);
	double wn1, wn2, reG1, reG2, denom2;
	for (j=1; j<Nn-1; j++)
	{
		wn1=wn(j);
		wn2=wn(j+1);
		reG1=Gr(j);
		reG2=Gr(j+1);
		denom2=wn1*wn1 - wn2*wn2;
		C2b(j-1)=-(reG1*pow(wn1,4) - reG2*pow(wn2,4))/denom2;
		C4b(j-1)=(-reG1*pow(wn1,4)*pow(wn2,2) + reG2*pow(wn1,2)*pow(wn2,4))/denom2;
	}
	
	int p, Nfitmin, jfitmin, jfitmax, NNfit, Nfit;
	jfitmin=2;
	Nfitmin=2*NC+4;
	jfitmax=Nn-Nfitmin;
	NNfit=jfitmax-jfitmin+1;
	
	vec M1v(NNfit), M3v(NNfit);
	mat invCG, A, X, CG;
	vec Mtmp;
	
	//	char UPLO='U';
	//	int NA=NC, NRHS=1, INFO;
	
//	mat LC, invLC;
	for (jfit=jfitmin; jfit<=jfitmax; jfit++)
	{
		Nfit=Nn_fit_max;
		if ((Nn-jfit+1)<Nn_fit_max)
			Nfit=Nn-jfit+1;
		
		X.zeros(Nfit,NC);
		
		for (j=1; j<=NC; j++)
		{
			for (p=jfit; p<=jfit+Nfit-1; p++)
			{
				X((p-jfit),j-1)=pow(-1,j)/pow(wn(p-1),2*j);
			}
		}
		
		CG=COV.submat(jfit-1,jfit-1,jfit+Nfit-2,jfit+Nfit-2);
		
		//		if (jfit==jfitmin) cout<<CG.submat(0,0,10,10);
		
		//		LC=chol(CG);
		//		invLC=inv(LC);
		//		invCG=invLC*invLC.t();
		
		invCG=inv(CG);
//		invCG=inv_sympd(CG);
		A=trans(X)*invCG*X;
		A=0.5*(A+A.t());
		Mtmp=trans(X)*invCG*Gchi2.rows(jfit-1,jfit+Nfit-2);
		
		//		dposv_(&UPLO, &NA, &NRHS, A.memptr(), &NA, Mtmp.memptr(), &NA, &INFO);
		Mtmp=solve(A,Mtmp);
		
		M1v(jfit-jfitmin)=Mtmp(0);
		M3v(jfit-jfitmin)=Mtmp(1);
		
	}
	
	int Nv=Nn/16;
	if (Nv<2) Nv=2;
	
	vec varM1(NNfit-2*Nv), varM3(NNfit-2*Nv);
	
	for (j=Nv; j<NNfit-Nv; j++)
	{
		varM1(j-Nv)=var(M1v.rows(j-Nv,j+Nv));
		varM3(j-Nv)=var(M3v.rows(j-Nv,j+Nv));
	}
	
	uword j1, j3;
	varM1.min(j1);
	varM3.min(j3);
	
	j1=j1+Nv;
	j3=j3+Nv;
	
	Mfit.zeros(2);
	Mfit(0)=mean(M1v.rows(j1-Nv,j1+Nv));
	Mfit(1)=mean(M3v.rows(j3-Nv,j3+Nv));
	
	
	int jfit0;
	j=Nv;
	while ((abs(mean(C2b.rows(j-Nv,j+Nv))-Mfit(0))/Mfit(0)>tol_mean_C1 || stddev(C2b.rows(j-Nv,j+Nv))/Mfit(0)>tol_std_C1) && j<Nn-Nv-3) j++;
	jfit0=j;
	
	jfit=j1+jfitmin-1;
	Nfit=Nn_fit_fin;
	if (Nfit>Nn-jfit+1)
		Nfit=Nn-jfit+1;
	
	X.zeros(Nfit,NC);
	
	for (j=1; j<=NC; j++)
	{
		for (p=jfit; p<=jfit+Nfit-1; p++)
		{
			X((p-jfit),j-1)=pow(-1,j)/pow(wn(p-1),2*j);
		}
	}
	
	CG=COV.submat(jfit-1,jfit-1,jfit+Nfit-2,jfit+Nfit-2);
	
	//	LC=chol(CG);
	//	invLC=inv(LC);
	//	invCG=invLC*invLC.t();
	
	invCG=inv(CG);
//		invCG=inv_sympd(CG);
	A=trans(X)*invCG*X;
	A=0.5*(A+A.t());
	mat COVMtmp=inv(A);
	COVMfit=COVMtmp.submat(0,0,0,0);
	
	jfit=jfit0;
	
	cout<<"frequency range of asymptotic behavior: "<<wn(jfit-1)<<" to "<<wn(Nn-1)<<" (indices "<<jfit-1<<" to "<<Nn-1<<")"<<endl;
	
	cout<<"1st moment extracted from high frequencies: "<<Mfit(0)<<endl;
	cout<<"3rd moment extracted from high frequencies: "<<Mfit(1)<<endl;
	
	double std_omega_tmp;
	double var_omega;
		
	var_omega=Mfit(0)/M1n;
	if (var_omega>0)
		std_omega_tmp=sqrt(var_omega);
	else
	{
		cout<<"Negative variance found during computation of moments.\n";
		return false;
	}
		
	if (!moments_provided)
	{
		M=Mfit.rows(0,0);
		NM=1;
		COVM=COVMfit.submat(0,0,0,0);
		M1=M(0);
		covm_diag=true;
		M1_set=true;
		M_ord.zeros(1);
		M_ord(0)=1;
		errM.zeros(NM);
		errM(0)=sqrt(COVM(0,0));
	}
	else if (abs((M1-Mfit(0))/Mfit(0))>tol_M1)
		cout<<"warning: first moment different from provided one\n";
	
	if (!std_omega)
	{
		var_omega=M1/M1n;
		if (var_omega>0)
			std_omega=sqrt(var_omega);
		else
		{
			cout<<"Negative variance found during computation of moments.\n";
			return false;
		}
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
	
	if (Nn-jfit<Nn_as_min)
	{
		jfit=0;
		if (!moments_provided) maxM=0;
	}
	
	return true;
}










bool OmegaMaxEnt_data::test_low_energy_peak_chi()
{
//	if (Nwn_test_metal>Nn-2)	Nwn_test_metal=Nn-2;
	
	peak_exists=false;
//	bool test_peak=false;
	
//	vec D2G=(Gr.rows(0,Nn-3)-2*Gr.rows(1,Nn-2)+Gr.rows(2,Nn-1));
//	if (abs(D2G(0)/Gr(1))>R_d2G_chi_peak && Nn>20) test_peak=true;
	
//	if (test_peak)
//	{
//		cout<<"Re(G) has a jump at low frequency. Looking for a peak in the spectral function at low energy...\n";
		cout<<"Looking for a peak in the spectral function at low energy...\n";
		
		int nmax=n(Nn-1);
		
		int NCnmin=1;
		int NCnmax=3;
		int NNCn=NCnmax-NCnmin+1;
		ivec NCn=linspace<ivec>(NCnmin,NCnmax,NNCn);
		
		int NCpmin=1;
		int NCpmax=15;
		if (NCpmax>(nmax-NCn.max()-2))
		{
			NCpmax=nmax-NCn.max()-2;
		}
		int NNCp=NCpmax-NCpmin+1;
		ivec NCp=linspace<ivec>(NCpmin,NCpmax,NNCp);
		
		int DNwn=2;
		
		int p0_min=2;
		int p0_max=2;
		int Np=p0_max-p0_min+1;
		ivec p0=linspace<ivec>(p0_min,p0_max,Np);
		
		imat NCmin(Np,2);
		
		vec M0_inc_opt(Np,fill::zeros);
		vec M2_pk_opt(Np,fill::zeros);
		vec varM2_pk_opt(Np,fill::zeros);
		
		mat X, CG, P, invCG, AM, BM, AMP, BMP, MP, Mtmp, chi2tmp;
		mat M0_inc(NNCp,NNCn), M2_pk(NNCp,NNCn), chi2_pk(NNCp,NNCn);
		int j, l, m, p, Nfit;
		uword q;
		vec Gtmp, Glftmp, diffG;
		rowvec maxX;
		//double vartmp;
		
		for (q=0; q<Np; q++)
		{
			M0_inc.zeros();
			M2_pk.zeros();
			chi2_pk.zeros();
			
			for (l=0; l<NNCp; l++)
			{
				for (m=0; m<NNCn; m++)
				{
					Nfit=NCp(l)+NCn(m)+1+DNwn;
					
					X.zeros(Nfit,NCp(l)+NCn(m)+1);
					
					for (j=NCp(l); j>=-NCn(m); j--)
					{
						for (p=p0(q); p<=Nfit+p0(q)-1; p++)
						{
							X(p-p0(q),NCp(l)-j)=pow(wn(p-1),2*j);
						}
					}
					
					Gtmp=Gchi2.rows(p0(q)-1,Nfit+p0(q)-2);
					CG=COV.submat(p0(q)-1,p0(q)-1,Nfit+p0(q)-2,Nfit+p0(q)-2);
					
					invCG=inv(CG);
					//	invCG=inv_sympd(CG);
					AM=(X.t())*invCG*X;
					BM=(X.t())*invCG*Gtmp;
					
					Mtmp=solve(AM,BM);
					
					Glftmp=X*Mtmp;
					
					diffG=Gtmp-Glftmp;
					
					chi2tmp=((diffG.t())*invCG*diffG)/Nfit;
					
					chi2_pk(l,m)=chi2tmp(0,0);
					M0_inc(l,m)=-Mtmp(NCp(l));
					M2_pk(l,m)=-Mtmp(NCp(l)+1);
				}
			}
			
			uint NvM=2;
			uint NvarM=NNCp-2*NvM;
			mat varM2_pk(NvarM,NNCn,fill::zeros);
			for (j=NvM; j<NNCp-NvM; j++)
				varM2_pk.row(j-NvM)=var(M2_pk.rows(j-NvM,j+NvM));
			
			uword indpM0, indnM0;
			
			double varM2_pk_min=varM2_pk.min(indpM0,indnM0);
			
			l=indpM0+NvM;
			m=indnM0;
			
			NCmin(q,0)=l;
			NCmin(q,1)=m;
			
			M0_inc_opt(q)=M0_inc(l,m);
			M2_pk_opt(q)=M2_pk(l,m);
			varM2_pk_opt(q)=varM2_pk_min;
		}
		
		double varM2_min=varM2_pk_opt.min(q);
		
		double peak_weight=-Gr(0)-M0_inc_opt(q);
		double var_peak=M2_pk_opt(q)/peak_weight;
		double peak_width=0;
		if (var_peak>0)	peak_width=sqrt(var_peak);
		//double m2_lf_min=M2_pk_opt(q);
		
		if (peak_width>100*EPSILON && varM2_min/M2_pk_opt(q)<varM2_peak_max && peak_weight>peak_weight_min*M1n && var_peak>0)
		{
			peak_exists=true;
			
			l=NCmin(q,0);
			m=NCmin(q,1);
			
			Nfit=NCp(l)+NCn(m)+1+DNwn;
			
			X.zeros(Nfit,NCp(l)+NCn(m)+1);
			
			for (j=NCp(l); j>=-NCn(m); j--)
			{
				for (p=p0(q); p<=Nfit+p0(q)-1; p++)
				{
					X(p-p0(q),NCp(l)-j)=pow(wn(p-1),2*j);
				}
			}
			
			Gtmp=Gchi2.rows(p0(q)-1,Nfit+p0(q)-2);
			CG=COV.submat(p0(q)-1,p0(q)-1,Nfit+p0(q)-2,Nfit+p0(q)-2);
			
			invCG=inv(CG);
			//	invCG=inv_sympd(CG);
			AM=(X.t())*invCG*X;
			AM=0.5*(AM.t()+AM);
			BM=(X.t())*invCG*Gtmp;
			
			Mtmp=solve(AM,BM);
			
			Glftmp=X*Mtmp;
			
			mat COVpeak=inv(AM);
			
			//double err_norm_peak=sqrt(COVpeak(NCp(l),NCp(l)));
			//double err_M2_peak=sqrt(COVpeak(NCp(l)+1,NCp(l)+1));
			//double err_std_peak=(-err_norm_peak*m2_lf_min/peak_weight+err_M2_peak/peak_weight)/(2*peak_width);
			
			cout<<"Peak detected\n";
			cout<<"peak width: "<<peak_width<<endl;
			//cout<<"error on width: "<<err_std_peak<<endl;
			cout<<"peak weight: "<<peak_weight<<endl;
			//cout<<"error on weight: "<<err_norm_peak<<endl;
			
			dw_peak=peak_width/2.0;
		}
		else
		{
			cout<<"no peak found\n";
		}
//	}
	
	return peak_exists;
}








bool OmegaMaxEnt_data::set_omega_grid_chi()
{
	if (!SW_set || !wc_exists) return false;
	
	int j;
	
	double wmax=f_w_range*SW/2.0;
	if (wr>0)
	{
		if (wmax<R_wmax_wr_min*wr) wmax=R_wmax_wr_min*wr;
	}
	w0r=wr-sqrt(dwr*(wmax-wr));
	int Nur=ceil((wr-w0r)/dwr);
	w0r=wr-Nur*dwr;
	double dur=dwr/((wr-w0r)*(wr+dwr-w0r));
	wmax=1.0/dur+w0r;
	ivec ur_int=linspace<ivec>(Nur,1,Nur);
	vec ur(Nur);
	for (j=0; j<Nur; j++)
		ur(j)=dur*ur_int(j);
	
	vec wr_vec=1.0/ur+w0r;
	Nw=Nur+Nwc;
	w.zeros(Nw);
	w.rows(0,Nwc-1)=wc;
	w.rows(Nwc,Nw-1)=wr_vec;
	w_exists=true;
	ws.zeros(1);
	ws(0)=w0r;
	Nw_lims.zeros(1);
	Nw_lims(0)=Nwc-1;
	
	Du_constant=true;
	
	return true;
}













bool OmegaMaxEnt_data::truncate_chi_omega_n()
{
	double wn_max=wn(ind_cutoff_wn);
	wn=wn.rows(0,ind_cutoff_wn);
	n=n.rows(0,ind_cutoff_wn);
	Nn=ind_cutoff_wn+1;
	if (Nn>Nn_max) // || Nn>=Nw)
	{
		cout<<"Number of Matsubara frequencies larger than the maximum allowed...\n";
		uvec diff_n=n.rows(1,Nn-1)-n.rows(0,Nn-2);
		if (diff_n.max()==1)
		{
			cout<<"Using a non-uniform Matsubara grid.\n";
			int i;
			int p=ceil(log2(Nn));
			int N0=pow(2,p);
			if (N0+1>Nn_all)
			{
				p=p-1;
				N0=pow(2,p);
			}
			int r=1;
			int N1=N0/pow(2,r);
			int N2=N1/2;
			Nn=N1+r*N2+1;
			n.zeros(Nn);
			n.rows(0,N1-1)=linspace<uvec>(0,N1-1,N1);
			uvec j=linspace<uvec>(0,Nn-N1-1,Nn-N1);
			uvec lj=j/N2;
			for (i=0; i<Nn-N1; i++) n(i+N1)=(j(i) % N2)*pow(2,lj(i)+1) + N1*pow(2,lj(i));
			wn=2*PI*tem*conv_to<vec>::from(n);
			i=Nn-1;
			while (wn(i)>wn_max && i>0) i--;
			Nn=i+1;
			while (Nn>Nn_max) //|| Nn>=Nw)
			{
				r=r+1;
				N1=N0/pow(2,r);
				N2=N1/2;
				Nn=N1+r*N2+1;
				n.zeros(Nn);
				n.rows(0,N1-1)=linspace<uvec>(0,N1-1,N1);
				uvec j=linspace<uvec>(0,Nn-N1-1,Nn-N1);
				uvec lj=j/N2;
				for (i=0; i<Nn-N1; i++) n(i+N1)=(j(i) % N2)*pow(2,lj(i)+1) + N1*pow(2,lj(i));
				wn=2*PI*tem*conv_to<vec>::from(n);
				i=Nn-1;
				while (wn(i)>wn_max && i>0) i--;
				Nn=i+1;
			}
			wn=wn.rows(0,Nn-1);
			n=n.rows(0,Nn-1);
			G=G.rows(n);
			Gr=real(G);
			Gi=imag(G);
			
			Gchi2=Gr;
			COV=COV(n,n);
			
			if (cov_diag)
			{
				errGr=errGr.rows(n);
				errG=errGr;
			}
			
		}
		else
		{
			cout<<"Your Matsubara grid is not uniform. Either make it more sparse to reduce further the number of frequencies, or increase \"Nn_max\" in the file \"OmegaMaxEnt_other_params.dat\".\n";
			return false;
		}
	}
	else
	{
		G=G.rows(0,Nn-1);
		Gr=real(G);
		Gi=imag(G);
		Gchi2=Gr;
		COV=COV.submat(0,0,Nn-1,Nn-1);
		if (cov_diag)
		{
			errGr=errG.rows(0,Nn-1);
			errG=errGr;
		}
	}
	
	return true;
}














bool OmegaMaxEnt_data::Kernel_chi()
{
	dcomplex i(0,1);
	
	double fg=1.7;
	double fi=fg;
	int pnmax=100;
	
	cout<<"defining kernel matrix...\n";
	
	Kcx.zeros(Nn,Nw);
	KM.zeros(1,Nw);
 
	int Nud=Nw-Nw_lims(0)-1;
 
	vec ud;
	if (Du_constant)
	{
		double dud=dwr/((wr-w0r)*(wr+dwr-w0r));
		vec ud_int=linspace<vec>(Nud,1,Nud);
		ud=dud*ud_int;
	}
	else
	{
		ud=1.0/(w.rows(Nw_lims(0)+1,Nw-1)-w0r);
	}

	mat MM;
	spline_matrix_chi(w, Nw_lims, ws, MM);

	int Nint=Nw;
	int	Nintc=Nwc-1;

	int Ncfs0=3*Nint-1;

	mat MC(Ncfs0+Nw,Nw);
	
	MC.submat(0,0,Ncfs0-1,Nw-1)=MM;
	MC.submat(Ncfs0,0,Ncfs0+Nw-1,Nw-1)=eye<mat>(Nw,Nw);

	mat	Pa_c=zeros<mat>(Nintc,Ncfs0+Nw);
	mat	Pb_c_r=zeros<mat>(Nintc,Ncfs0+Nw);
	mat	Pc_c_r=zeros<mat>(Nintc,Ncfs0+Nw);
	mat	Pd_c_r=zeros<mat>(Nintc,Ncfs0+Nw);
	
	int j;
	for (j=0; j<Nintc; j++)
	{
		Pa_c(j,3*j)=1;
		Pb_c_r(j,3*j+1)=1;
		Pc_c_r(j,3*j+2)=1;
		Pd_c_r(j,j+Ncfs0)=1;
	}

	mat W=diagmat(wc.rows(0,Nwc-2));
	
	mat Pb_c=Pb_c_r-3*W*Pa_c;
	mat Pc_c=Pc_c_r+3*pow(W,2)*Pa_c-2*W*Pb_c_r;
	mat Pd_c=Pd_c_r-pow(W,3)*Pa_c+pow(W,2)*Pb_c_r-W*Pc_c_r;
	
	int Nintd=Nud+1;
	
	mat Pa_d=zeros<mat>(Nintd,Ncfs0+Nw);
	mat Pb_d_r=zeros<mat>(Nintd,Ncfs0+Nw);
	mat Pc_d_r=zeros<mat>(Nintd,Ncfs0+Nw);
	mat Pd_d_r=zeros<mat>(Nintd,Ncfs0+Nw);
	
	for (j=0; j<Nud; j++)
	{
		Pa_d(j,3*j+3*Nintc)=1;
		Pb_d_r(j,3*j+1+3*Nintc)=1;
		Pc_d_r(j,3*j+2+3*Nintc)=1;
		Pd_d_r(j,j+Ncfs0+Nwc)=1;
	}
	
	j=Nud;
	Pa_d(j,3*j+3*Nintc)=1;
	Pb_d_r(j,3*j+1+3*Nintc)=1;
	
	mat U=zeros<mat>(Nud+1,Nud+1);
	vec ud1=zeros<vec>(Nud+1);
	ud1.rows(0,Nud-1)=ud;
	U.diag()=ud1;
	
	mat Pb_d=Pb_d_r-3*U*Pa_d;
	mat Pc_d=Pc_d_r+3*pow(U,2)*Pa_d-2*U*Pb_d_r;
	mat Pd_d=Pd_d_r-pow(U,3)*Pa_d+pow(U,2)*Pb_d_r-U*Pc_d_r;

	mat Ka_c=zeros<mat>(Nn,Nintc);
	mat Kb_c=zeros<mat>(Nn,Nintc);
	mat Kc_c=zeros<mat>(Nn,Nintc);
	mat Kd_c=zeros<mat>(Nn,Nintc);
 
	Ka_c.row(0)=-trans(pow(wc.rows(1,Nwc-1),4)-pow(wc.rows(0,Nwc-2),4))/(8*PI);
	Kb_c.row(0)=-trans(pow(wc.rows(1,Nwc-1),3)-pow(wc.rows(0,Nwc-2),3))/(6*PI);
	Kc_c.row(0)=-trans(pow(wc.rows(1,Nwc-1),2)-pow(wc.rows(0,Nwc-2),2))/(4*PI);
	Kd_c.row(0)=-trans(wc.rows(1,Nwc-1)-wc.rows(0,Nwc-2))/(2*PI);
	
//	cx_vec_tmp.zeros(Nintc);
//	vec_tmp=-trans(pow(wc.rows(1,Nwc-1),4)-pow(wc.rows(0,Nwc-2),4))/(8*PI);
//	cx_vec_tmp.set_real(vec_tmp);
//	Ka_c.row(0)=cx_vec_tmp;
//	vec_tmp=-trans(pow(wc.rows(1,Nwc-1),3)-pow(wc.rows(0,Nwc-2),3))/(6*PI);
//	cx_vec_tmp.set_real(vec_tmp);
//	Kb_c.row(0)=cx_vec_tmp;
//	vec_tmp=-trans(pow(wc.rows(1,Nwc-1),2)-pow(wc.rows(0,Nwc-2),2))/(4*PI);
//	cx_vec_tmp.set_real(vec_tmp);
//	Kc_c.row(0)=cx_vec_tmp;
//	vec_tmp=-trans(wc.rows(1,Nwc-1)-wc.rows(0,Nwc-2))/(2*PI);
//	cx_vec_tmp.set_real(vec_tmp);
//	Kd_c.row(0)=cx_vec_tmp;
 
	mat Wnc=wn.rows(1,Nn-1)*ones<rowvec>(Nintc);
	mat Wc=ones<vec>(Nn-1)*wc.t();
 
	mat logc=log((pow(Wnc,2)+pow(Wc.cols(1,Nintc),2))/(pow(Wnc,2)+pow(Wc.cols(0,Nintc-1),2)));
	mat atanc=atan((Wnc % (Wc.cols(0,Nintc-1)-Wc.cols(1,Nintc)))/(Wc.cols(1,Nintc) % Wc.cols(0,Nintc-1)+pow(Wnc,2)));
	
	mat dWc=Wc.cols(1,Nintc)-Wc.cols(0,Nintc-1);
	mat dWc2=pow(Wc.cols(1,Nintc),2)-pow(Wc.cols(0,Nintc-1),2);
	mat dWc3=pow(Wc.cols(1,Nintc),3)-pow(Wc.cols(0,Nintc-1),3);
	mat Wnc2=pow(Wnc,2);
	mat Wnc3=pow(Wnc,3);
	mat Wnc4=pow(Wnc,4);
 
	Ka_c.rows(1,Nn-1)=real(-i*( -Wnc3 % dWc + i*(Wnc2 % dWc2)/2 + (Wnc % dWc3)/3 -i*(pow(Wc.cols(1,Nintc),4)-pow(Wc.cols(0,Nintc-1),4))/4 - Wnc4 % atanc - i*Wnc4 % logc/2 ))/(2*PI);
	Kb_c.rows(1,Nn-1)=real(-i*( i*Wnc2 % dWc + (Wnc % dWc2)/2 - i*dWc3/3 + i*Wnc3 % atanc - Wnc3 % logc/2 ))/(2*PI);
	Kc_c.rows(1,Nn-1)=real(-i*( Wnc % dWc - i*dWc2/2 + Wnc2 % atanc + i*Wnc2 % logc/2 ))/(2*PI);
	Kd_c.rows(1,Nn-1)=real(-i*( -i*dWc -i*Wnc % atanc + Wnc % logc/2 ))/(2*PI);
	
 
	int Pmax=pnmax+4;
	imat MP;
	pascal(Pmax+1,MP);
	int pnmax2=pnmax/2;
 
	double wtmp, wj, dw1;
	vec dwp;
	int p,l, jni;
	for (j=0; j<Nwc-1; j++)
	{
		wtmp=abs(wc(j+1));
		if (abs(wc(j))>wtmp) wtmp=abs(wc(j));
		jni=1;
		while (wn(jni)<fi*wtmp && jni<Nn)	jni++;
		while (pow(wn(jni),pnmax)==0 && jni<Nn)		jni++;
		
		if (jni<Nn-1)
		{
			wj=wc(j);
			dw1=wc(j+1)-wj;
			dwp.zeros(Pmax);
			dwp(0)=dw1;
			dwp(1)=pow(dw1,2)+2*wj*dw1;
			for (p=3; p<=Pmax; p++)
			{
				dwp(p-1)=0;
				for (l=0; l<p; l++)
				{
					dwp(p-1)=dwp(p-1)+MP(l,p-l)*pow(dw1,p-l)*pow(wj,l);
				}
			}
			if (jni<Nn-1)
			{
				Ka_c.submat(jni,j,Nn-1,j)=zeros<vec>(Nn-jni);
				Kb_c.submat(jni,j,Nn-1,j)=zeros<vec>(Nn-jni);
				Kc_c.submat(jni,j,Nn-1,j)=zeros<vec>(Nn-jni);
				Kd_c.submat(jni,j,Nn-1,j)=zeros<vec>(Nn-jni);
//				for (p=pnmax; p>=1; p--)
				for (l=pnmax2; l>=1; l--)
				{
					p=2*l;
					Ka_c.submat(jni,j,Nn-1,j)=Ka_c.submat(jni,j,Nn-1,j) + pow(-1,l)*dwp(p+3)/(2*PI*(p+4)*pow(wn.rows(jni,Nn-1),p));
					Kb_c.submat(jni,j,Nn-1,j)=Kb_c.submat(jni,j,Nn-1,j) + pow(-1,l)*dwp(p+2)/(2*PI*(p+3)*pow(wn.rows(jni,Nn-1),p));
					Kc_c.submat(jni,j,Nn-1,j)=Kc_c.submat(jni,j,Nn-1,j) + pow(-1,l)*dwp(p+1)/(2*PI*(p+2)*pow(wn.rows(jni,Nn-1),p));
					Kd_c.submat(jni,j,Nn-1,j)=Kd_c.submat(jni,j,Nn-1,j) + pow(-1,l)*dwp(p)/(2*PI*(p+1)*pow(wn.rows(jni,Nn-1),p));
				}
			}
		}
	}
	
	mat Ka_d=zeros<mat>(Nn,Nintd);
	mat Kb_d=zeros<mat>(Nn,Nintd);
	mat Kc_d=zeros<mat>(Nn,Nintd);
	mat Kd_d=zeros<mat>(Nn,Nintd);
	
	rowvec ud2(Nud+1);
	ud2.cols(1,Nud)=ud.t();
	ud2(0)=1.0/(wr-w0r);
	
	mat Wnd=wn.rows(1,Nn-1)*ones<rowvec>(Nintd-1);
	mat Ud=ones<vec>(Nn-1)*ud2;
	
	rowvec wd=1/ud2+w0r;
	mat Wd=ones<vec>(Nn-1)*wd;

	Ka_d.submat(0,0,0,Nintd-2)=-(pow(ud2.cols(1,Nintd-1),2)-pow(ud2.cols(0,Nintd-2),2))/(4*PI);
	Kb_d.submat(0,0,0,Nintd-2)=-(ud2.cols(1,Nintd-1)-ud2.cols(0,Nintd-2))/(2*PI);
	Kc_d.submat(0,0,0,Nintd-2)=-log(ud2.cols(1,Nintd-1)/ud2.cols(0,Nintd-2))/(2*PI);
	Kd_d.submat(0,0,0,Nintd-2)=(wd.cols(1,Nintd-1)-wd.cols(0,Nintd-2))/(2*PI);
	
//	cx_vec_tmp.zeros(Nud);
//	vec_tmp=-(pow(ud2.cols(1,Nintd-1),2)-pow(ud2.cols(0,Nintd-2),2))/(4*PI);
//	cx_vec_tmp.set_real(vec_tmp);
//	Ka_d.submat(0,0,0,Nintd-2)=cx_vec_tmp;
//	vec_tmp=-(ud2.cols(1,Nintd-1)-ud2.cols(0,Nintd-2))/(2*PI);
//	cx_vec_tmp.set_real(vec_tmp);
//	Kb_d.submat(0,0,0,Nintd-2)=cx_vec_tmp;
//	vec_tmp=-log(ud2.cols(1,Nintd-1)/ud2.cols(0,Nintd-2))/(2*PI);
//	cx_vec_tmp.set_real(vec_tmp);
//	Kc_d.submat(0,0,0,Nintd-2)=cx_vec_tmp;
//	vec_tmp=(wd.cols(1,Nintd-1)-wd.cols(0,Nintd-2))/(2*PI);
//	cx_vec_tmp.set_real(vec_tmp);
//	Kd_d.submat(0,0,0,Nintd-2)=cx_vec_tmp;
	
	
	mat atand=atan((Wnd % (Ud.cols(1,Nintd-1)-Ud.cols(0,Nintd-2)))/(1+w0r*(Ud.cols(1,Nintd-1)+Ud.cols(0,Nintd-2))+(pow(w0r,2)+pow(Wnd,2)) % Ud.cols(0,Nintd-2) % Ud.cols(1,Nintd-1)));
	mat logd=log((1+2*w0r*Ud.cols(1,Nintd-1)+pow(Ud.cols(1,Nintd-1),2) % (pow(Wnd,2)+pow(w0r,2)))/(1+2*w0r*Ud.cols(0,Nintd-2)+pow(Ud.cols(0,Nintd-2),2) % (pow(Wnd,2)+pow(w0r,2))));
	
	Ka_d.submat(1,0,Nn-1,Nintd-2)=real(-i*(2*Wnd % (Ud.cols(1,Nintd-1)-Ud.cols(0,Nintd-2)))/pow(Wnd + i*w0r,2) - i*(w0r*(pow(Ud.cols(1,Nintd-1),2)-pow(Ud.cols(0,Nintd-2),2)))/(Wnd + i*w0r) +  i*(2*Wnd % atand)/pow(Wnd + i*w0r,3)  - (Wnd % logd)/pow(Wnd + i*w0r,3) )/(4*PI);
	Kb_d.submat(1,0,Nn-1,Nintd-2)=real(-i*(2*w0r*(Ud.cols(1,Nintd-1)-Ud.cols(0,Nintd-2)))/(Wnd + i*w0r) -(2*Wnd % atand)/pow(Wnd + i*w0r,2) -i*(Wnd % logd)/pow(Wnd + i*w0r,2) )/(4*PI);
	Kc_d.submat(1,0,Nn-1,Nintd-2)=real( -2*log(Ud.cols(1,Nintd-1)/Ud.cols(0,Nintd-2)) -i*(2*Wnd % atand)/(Wnd + i*w0r) + (Wnd % logd)/(Wnd + i*w0r) )/(4*PI);
	Kd_d.submat(1,0,Nn-1,Nintd-2)=real( 2*(Wd.cols(1,Nintd-1)-Wd.cols(0,Nintd-2)) -i*2.0*Wnd % log(Ud.cols(1,Nintd-1)/Ud.cols(0,Nintd-2)) + 2*Wnd % atand +i*Wnd % logd )/(4*PI);
	
	
	mat KC=(Ka_c*Pa_c+Kb_c*Pb_c+Kc_c*Pc_c+Kd_c*Pd_c)*MC;
	mat KD=-(Ka_d*Pa_d+Kb_d*Pb_d+Kc_d*Pc_d+Kd_d*Pd_d)*MC;
 
	Kcx.set_real(KC+KD);
	Kcx=2*Kcx;
	
	rowvec Knorm_a_c_r=zeros<rowvec>(Nintc);
	rowvec Knorm_b_c_r=zeros<rowvec>(Nintc);
	rowvec Knorm_c_c_r=zeros<rowvec>(Nintc);
	rowvec Knorm_d_c_r=zeros<rowvec>(Nintc);
	
	Knorm_a_c_r=trans(pow(wc.rows(1,Nwc-1)-wc.rows(0,Nwc-2),4)/4);
	Knorm_b_c_r=trans(pow(wc.rows(1,Nwc-1)-wc.rows(0,Nwc-2),3)/3);
	Knorm_c_c_r=trans(pow(wc.rows(1,Nwc-1)-wc.rows(0,Nwc-2),2)/2);
	Knorm_d_c_r=trans(wc.rows(1,Nwc-1)-wc.rows(0,Nwc-2));
	
	rowvec KM0c=(Knorm_a_c_r*Pa_c+Knorm_b_c_r*Pb_c_r+Knorm_c_c_r*Pc_c_r+Knorm_d_c_r*Pd_c_r)*MC/(2*PI);
	
	rowvec Knorm_a_d=zeros<rowvec>(Nintd);
	rowvec Knorm_b_d=zeros<rowvec>(Nintd);
	rowvec Knorm_c_d=zeros<rowvec>(Nintd);
	rowvec Knorm_d_d=zeros<rowvec>(Nintd);
	
	Knorm_a_d.cols(0,Nintd-2)=-(pow(ud2.cols(1,Nud),2)-pow(ud2.cols(0,Nud-1),2))/2;
	Knorm_b_d.cols(0,Nintd-2)=-(ud2.cols(1,Nud)-ud2.cols(0,Nud-1));
	Knorm_c_d.cols(0,Nintd-2)=-log(ud2.cols(1,Nud)/ud2.cols(0,Nud-1));
	Knorm_d_d.cols(0,Nintd-2)=1.0/ud2.cols(1,Nud)-1.0/ud2.cols(0,Nud-1);
	
	Knorm_a_d(Nintd-1)=pow(ud2(Nud-1),2)/2;
	Knorm_b_d(Nintd-1)=ud2(Nud-1);
	
	rowvec KM0d=(Knorm_a_d*Pa_d+Knorm_b_d*Pb_d+Knorm_c_d*Pc_d+Knorm_d_d*Pd_d)*MC/(2*PI);
	
	rowvec KM0=KM0c+KM0d;
	
	mat Wjc=diagmat(wc.rows(0,Nintc-1));
	
	rowvec KM1_a_c_tmp=trans(pow(wc.rows(1,Nwc-1)-wc.rows(0,Nwc-2),5)/5.0);
	
	rowvec KM1_a_c=KM1_a_c_tmp+Knorm_a_c_r*Wjc;
	rowvec KM1_b_c=Knorm_a_c_r+Knorm_b_c_r*Wjc;
	rowvec KM1_c_c=Knorm_b_c_r+Knorm_c_c_r*Wjc;
	rowvec KM1_d_c=Knorm_c_c_r+Knorm_d_c_r*Wjc;
	
	rowvec KM1c=(KM1_a_c*Pa_c+KM1_b_c*Pb_c_r+KM1_c_c*Pc_c_r+KM1_d_c*Pd_c_r)*MC/(2*PI);
	
	rowvec KM1_a_d=zeros<rowvec>(Nintd);
	rowvec KM1_b_d=zeros<rowvec>(Nintd);
	rowvec KM1_c_d=zeros<rowvec>(Nintd);
	rowvec KM1_d_d=zeros<rowvec>(Nintd);
	
	KM1_a_d.cols(0,Nintd-2)=Knorm_b_d.cols(0,Nintd-2);
	KM1_b_d.cols(0,Nintd-2)=Knorm_c_d.cols(0,Nintd-2);
	KM1_c_d.cols(0,Nintd-2)=Knorm_d_d.cols(0,Nintd-2);
	KM1_d_d.cols(0,Nintd-2)=(1.0/pow(ud2.cols(1,Nud),2)-1.0/pow(ud2.cols(0,Nud-1),2))/2;
	
	rowvec KM1d_tmp=(KM1_a_d*Pa_d+KM1_b_d*Pb_d+KM1_c_d*Pc_d+KM1_d_d*Pd_d)*MC/(2*PI);
	
	rowvec KM1d=w0r*KM0d+KM1d_tmp;
	
	rowvec KM1=KM1c+KM1d;
	
	rowvec KM2_a_c_tmp=trans(pow(wc.rows(1,Nwc-1)-wc.rows(0,Nwc-2),6)/6);
	
	rowvec KM2_a_c=KM2_a_c_tmp+2*KM1_a_c_tmp*Wjc+Knorm_a_c_r*pow(Wjc,2);
	rowvec KM2_b_c=KM1_a_c_tmp+2*Knorm_a_c_r*Wjc+Knorm_b_c_r*pow(Wjc,2);
	rowvec KM2_c_c=Knorm_a_c_r+2*Knorm_b_c_r*Wjc+Knorm_c_c_r*pow(Wjc,2);
	rowvec KM2_d_c=Knorm_b_c_r+2*Knorm_c_c_r*Wjc+Knorm_d_c_r*pow(Wjc,2);
	
	rowvec KM2c=(KM2_a_c*Pa_c+KM2_b_c*Pb_c_r+KM2_c_c*Pc_c_r+KM2_d_c*Pd_c_r)*MC/(2*PI);
	
	rowvec KM2_a_d=zeros<rowvec>(Nintd);
	rowvec KM2_b_d=zeros<rowvec>(Nintd);
	rowvec KM2_c_d=zeros<rowvec>(Nintd);
	rowvec KM2_d_d=zeros<rowvec>(Nintd);
	
	KM2_a_d.cols(0,Nintd-2)=Knorm_c_d.cols(0,Nintd-2);
	KM2_b_d.cols(0,Nintd-2)=Knorm_d_d.cols(0,Nintd-2);
	KM2_c_d.cols(0,Nintd-2)=KM1_d_d.cols(0,Nintd-2);
	KM2_d_d.cols(0,Nintd-2)=(1.0/pow(ud2.cols(1,Nud),3)-1.0/pow(ud2.cols(0,Nud-1),3))/3;
	
	rowvec KM2d_tmp=(KM2_a_d*Pa_d+KM2_b_d*Pb_d+KM2_c_d*Pc_d+KM2_d_d*Pd_d)*MC/(2*PI);
	
	rowvec KM2d=pow(w0r,2)*KM0d+2*w0r*KM1d_tmp+KM2d_tmp;
	
	rowvec KM2=KM2c+KM2d;
	
	KM.row(0)=2*KM2;
	K=real(Kcx);
	
	cout<<"kernel matrix defined.\n";


	
	return true;
}









bool OmegaMaxEnt_data::diagonalize_covariance_chi()
{
	mat VM, WM, VG, WG;
	
	if (NM>0)
	{
		if (covm_diag)
		{
			VM.eye(NM,NM);
			WM=diagmat(1.0/errM);
		}
		else
		{
			vec COVM_eig;
			mat VMtmp;
			if (!eig_sym(COVM_eig,VMtmp,COVM.submat(0,0,NM-1,NM-1)))
			{
				cout<<"diagonalize_covariance(): diagonalization of moments covariance matrix failed\n";
				return false;
			}
			if (COVM_eig.min()<=0)
			{
				cout<<"error: the moments covariance matrix has non-positive eigenvalues\n";
				return false;
			}
			VM.zeros(NM,NM);
			VM.submat(0,0,NM-1,NM-1)=VMtmp;
			vec errM_eig(NM);
			errM_eig.rows(0,NM-1)=sqrt(COVM_eig);
			WM=diagmat(1.0/errM_eig);
		}
		
		
		KM=KM.rows(0,NM-1);
		M_V=WM*VM.t()*M;
		KM_V=WM*VM.t()*KM;
		
		KM_V=KM_V.cols(ind0,Nw-2);
		KM=KM.cols(ind0,Nw-2);
	}
	
	if (cov_diag)
	{
		VG.eye(Nn,Nn);
		WG=diagmat(1.0/errGr);
	}
	else
	{
		vec COVG_eig;
		mat VGtmp;
		if (!eig_sym(COVG_eig,VGtmp,COV))
		{
			cout<<"diagonalize_covariance(): diagonalization of covariance matrix failed\n";
			return false;
		}
		mat PV=eye<mat>(Nn,Nn);
		PV=flipud(PV);
		COVG_eig=flipud(COVG_eig);
		VG=VGtmp*PV;
		vec GVtmp=VG.t()*Gchi2;
		mat SGV=diagmat(sign(GVtmp));
		SGV=SGV-abs(SGV)+eye<mat>(Nn,Nn);
		VG=VG*SGV;
		WG=diagmat(1.0/sqrt(COVG_eig));
		if (COVG_eig.min()<=0)
		{
			cout<<"error: the covariance matrix has non-positive eigenvalues\n";
			return false;
		}
	}
	
	G_V=WG*VG.t()*Gchi2;
	KG_V=WG*VG.t()*K.cols(ind0,Nw-2);
	K=K.cols(ind0,Nw-2);
	
	if (NM>0)
	{
		GM=join_vert(M_V, G_V);
		KGM=join_vert(KM_V, KG_V);
	}
	else
	{
		GM=G_V;
		KGM=KG_V;
	}
	NGM=GM.n_rows;
	
	cout<<"number of terms in chi2: "<<NGM<<endl;
	
	return true;
}





bool OmegaMaxEnt_data::spline_val_chi_part(vec x, vec x0, uvec ind_xlims, vec xs, vec coeffs, vec &sv)
{
	int Nx=x.n_rows;
	int Nx0=x0.n_rows;
	
	int Ncfs=4*Nx0;
	if (coeffs.n_rows<Ncfs)
	{
		cout<<"spline_val_chi_part(): number of elements in vectors \"coeffs\" and \"x0\" are not consistant\n";
		return false;
	}
	
	//	int Nc=ind_xlims(0);
	//	int Nd=Nx0-ind_xlims(1);
	
	double xr=x0(ind_xlims(0));
	
	double x0r=xs(0);
	
	double a,b,c,d, Dx, Du;
	
	sv.zeros(Nx);
	
	int j, l;
    for (j=0; j<Nx; j++)
	{
		if (x(j)>=0 && x(j)<xr)
		{
			l=0;
			while (x(j)>=x0(l+1) && l<ind_xlims(0)-1) l++;
			a=coeffs(4*l);
			b=coeffs(4*l+1);
			c=coeffs(4*l+2);
			d=coeffs(4*l+3);
			Dx=x(j)-x0(l);
			sv(j)=a*pow(Dx,3)+b*pow(Dx,2)+c*Dx+d;
		}
		else if (x(j)>xr)
		{
			l=ind_xlims(0);
			while (x(j)>x0(l+1) && l<Nx0-1) l++;
			a=coeffs(4*l);
			b=coeffs(4*l+1);
			c=coeffs(4*l+2);
			d=coeffs(4*l+3);
			if (l<Nx0-1)
				Du=(x0(l+1)-x(j))/((x(j)-x0r)*(x0(l+1)-x0r));
			else
				Du=1.0/(x(j)-x0r);
			sv(j)=a*pow(Du,3)+b*pow(Du,2)+c*Du+d;
		}
		else
		{
			l=ind_xlims(0)-1;
			a=coeffs(4*l);
			b=coeffs(4*l+1);
			c=coeffs(4*l+2);
			d=coeffs(4*l+3);
			Dx=x(j)-x0(l);
			sv(j)=a*pow(Dx,3)+b*pow(Dx,2)+c*Dx+d;
		}
	}
	
	return true;
}











bool OmegaMaxEnt_data::sum_gaussians_chi(vec x, rowvec x0, rowvec s0, rowvec wgt, vec &F)
{
	int Npks=x0.n_cols;
	
	F.zeros(x.n_rows);

	double W=sum(wgt.cols(0,Npks-1)), s2;
	
	int j;
	for (j=0; j<Npks; j++)
	{
		s2=pow(s0(j),2);
		F=F+wgt(j)*(exp(-pow(x-x0(j),2)/(2*s2)) + exp(-pow(x+x0(j),2)/(2*s2)))/(2*s0(j));
	}
	
	F=sqrt(2*PI)*F/W;
	
	return true;
}

















bool OmegaMaxEnt_data::spline_chi_part(vec x, uvec ind_xlims, vec xs, vec F, vec &coeffs)
{
	int j;
	
	int Nx=x.n_rows;
	
	int Nc=ind_xlims(0);
	int Nd=Nx-ind_xlims(0);
	
	double x0d=xs(0);
	
	//solve the spline in the interval [0,wr]
	int NSc=3*Nc-1;
	
	mat Aspl=zeros<mat>(NSc,NSc);
	vec Dc=zeros<vec>(NSc);
	
	double Dx=x(1)-x(0);
	
	Aspl(0,0)=pow(Dx,3);
	Aspl(0,1)=pow(Dx,2);
	Dc(0)=F(1)-F(0);
	
	Aspl(1,0)=3*pow(Dx,2);
	Aspl(1,1)=2*Dx;
	Aspl(1,4)=-1;
	
	Aspl(2,0)=6*Dx;
	Aspl(2,1)=2;
	Aspl(2,3)=-2;
	
	for (j=2; j<Nc; j++)
	{
		Dx=x(j)-x(j-1);
	 
		Dc(3*j-3)=F(j)-F(j-1);
		
		Aspl(3*j-3,3*j-4)=pow(Dx,3);
		Aspl(3*j-3,3*j-3)=pow(Dx,2);
		Aspl(3*j-3,3*j-2)=Dx;
	 
		Aspl(3*j-2,3*j-4)=3*pow(Dx,2);
		Aspl(3*j-2,3*j-3)=2*Dx;
		Aspl(3*j-2,3*j-2)=1;
		Aspl(3*j-2,3*j+1)=-1;
	 
		Aspl(3*j-1,3*j-4)=6*Dx;
		Aspl(3*j-1,3*j-3)=2;
		Aspl(3*j-1,3*j)=-2;
	}
	
	j=Nc;
	Dx=x(j)-x(j-1);
	Dc(3*j-3)=F(j)-F(j-1);
	Aspl(3*j-3,3*j-4)=pow(Dx,3);
	Aspl(3*j-3,3*j-3)=pow(Dx,2);
	Aspl(3*j-3,3*j-2)=Dx;
	
	Dc(3*j-2)=(F(j+1)-F(j-1))/(x(ind_xlims(0)+1)-x(ind_xlims(0)-1));
	Aspl(3*j-2,3*j-4)=3*pow(Dx,2);
	Aspl(3*j-2,3*j-3)=2*Dx;
	Aspl(3*j-2,3*j-2)=1;
	
	vec Cc=solve(Aspl,Dc);
/*
	mat T=zeros<mat>(NSc,Nx);
	for (j=0; j<Nc; j++)
	{
		T(3*j,j)=-1;
		T(3*j,j+1)=1;
	}
	T(3*Nc-2,ind_xlims(0)-1)=-1/(x(ind_xlims(0)+1)-x(ind_xlims(0)-1));
	T(3*Nc-2,ind_xlims(0)+1)=1/(x(ind_xlims(0)+1)-x(ind_xlims(0)-1));
	
	Cc=Cc*T;
	
	mat Cc2=zeros<mat>(NSc+1,Nx);
	
	Cc2.rows(0,1)=Cc.rows(0,1);
	Cc2.rows(3,NSc)=Cc.rows(2,NSc-1);
*/
	
	//solve the spline in the interval [wr,inf
	int NSd=3*Nd-1;
	Aspl.zeros(NSd,NSd);
	vec Dd=zeros<vec>(NSd);
	
	double u=1/(x(ind_xlims(0))-x0d)-1/(x(ind_xlims(0)+1)-x0d);
	double ud=1/(x(ind_xlims(0))-x0d);
	
	Aspl(0,0)=-3*pow(ud,2)*pow(u,2);
	Aspl(0,1)=-2*pow(ud,2)*u;
	Aspl(0,2)=-pow(ud,2);
	
	j=ind_xlims(0);
	Dd(0)=(F(j+1)-F(j-1))/(x(j+1)-x(j-1));
	
	Aspl(1,0)=pow(u,3);
	Aspl(1,1)=pow(u,2);
	Aspl(1,2)=u;
	Dd(1)=F(j)-F(j+1);
	
	for (j=1; j<Nd-1; j++)
	{
		u=1/(x(ind_xlims(0)+j)-x0d)-1/(x(ind_xlims(0)+j+1)-x0d);
		
		Aspl(3*j-1,3*j-1)=-1;
		Aspl(3*j-1,3*j)=3*pow(u,2);
		Aspl(3*j-1,3*j+1)=2*u;
		Aspl(3*j-1,3*j+2)=1;
		
		Aspl(3*j,3*j-2)=-2;
		Aspl(3*j,3*j)=6*u;
		Aspl(3*j,3*j+1)=2;
	 
		Dd(3*j+1)=F(ind_xlims(0)+j)-F(ind_xlims(0)+j+1);
		Aspl(3*j+1,3*j)=pow(u,3);
		Aspl(3*j+1,3*j+1)=pow(u,2);
		Aspl(3*j+1,3*j+2)=u;
	}
	
	j=Nd-1;
	u=1/(x(ind_xlims(0)+j)-x0d);
	
	Aspl(3*j-1,3*j-1)=-1;
	Aspl(3*j-1,3*j)=3*pow(u,2);
	Aspl(3*j-1,3*j+1)=2*u;
	
	Aspl(3*j,3*j-2)=-2;
	Aspl(3*j,3*j)=6*u;
	Aspl(3*j,3*j+1)=2;
	
	Dd(3*j+1)=F(ind_xlims(0)+j);
	Aspl(3*j+1,3*j)=pow(u,3);
	Aspl(3*j+1,3*j+1)=pow(u,2);
	
	vec Cd=solve(Aspl,Dd);
	
	coeffs.zeros(4*Nx);
	coeffs(0)=Cc(0);
	coeffs(1)=Cc(1);
	coeffs(3)=F(0);
	uvec ind_tmp=linspace<uvec>(1,Nc-1,Nc-1);
	coeffs.rows(4*ind_tmp)=Cc(3*ind_tmp-1);
	coeffs.rows(4*ind_tmp+1)=Cc(3*ind_tmp);
	coeffs.rows(4*ind_tmp+2)=Cc(3*ind_tmp+1);
	coeffs.rows(4*ind_tmp+3)=F.rows(ind_tmp);
	
	ind_tmp=linspace<uvec>(0,Nd-2,Nd-1);
	coeffs.rows(4*ind_tmp+4*Nc)=Cd(3*ind_tmp);
	coeffs.rows(4*ind_tmp+1+4*Nc)=Cd(3*ind_tmp+1);
	coeffs.rows(4*ind_tmp+2+4*Nc)=Cd(3*ind_tmp+2);
	coeffs.rows(4*ind_tmp+3+4*Nc)=F.rows(ind_xlims(0)+ind_tmp+1);
	coeffs(4*Nx-4)=Cd(3*Nd-3);
	coeffs(4*Nx-3)=Cd(3*Nd-2);

/*
	T.zeros(NSd,Nx);
	T(0,ind_xlims(0)-1)=-1/(x(ind_xlims(0)+1)-x(ind_xlims(0)-1));
	T(0,ind_xlims(0)+1)=1/(x(ind_xlims(0)+1)-x(ind_xlims(0)-1));
	for (j=1; j<Nd; j++)
	{
		T(3*j-2,ind_xlims(0)+j-1)=1;
		T(3*j-2,ind_xlims(0)+j)=-1;
	}
	j=Nd;
	T(3*j-2,ind_xlims(0)+j-1)=1;
	
	Cd=Cd*T;
	
	int NS_tot=NSc+NSd+1;
	
	Mspl.zeros(NS_tot,Nx);
	Mspl.rows(0,NSc)=Cc2;
	Mspl.rows(NSc+1,NS_tot-1)=Cd;
*/
	return true;
}






bool OmegaMaxEnt_data::spline_matrix_chi(vec x, uvec ind_xlims, vec xs, mat &Mspl)
{
	int j;
	
	int Nx=x.n_rows;
	
	int Nc=ind_xlims(0);
	int Nd=Nx-ind_xlims(0);
	
	double x0d=xs(0);
	
	//solve the spline in the interval [0,wr]
	int NSc=3*Nc-1;
	
	mat Aspl=zeros<mat>(NSc,NSc);
	
	double Dx=x(1)-x(0);
	
	Aspl(0,0)=pow(Dx,3);
	Aspl(0,1)=pow(Dx,2);
	
	Aspl(1,0)=3*pow(Dx,2);
	Aspl(1,1)=2*Dx;
	Aspl(1,4)=-1;
	
	Aspl(2,0)=6*Dx;
	Aspl(2,1)=2;
	Aspl(2,3)=-2;
	
	for (j=2; j<Nc; j++)
	{
		Dx=x(j)-x(j-1);
	 
		Aspl(3*j-3,3*j-4)=pow(Dx,3);
		Aspl(3*j-3,3*j-3)=pow(Dx,2);
		Aspl(3*j-3,3*j-2)=Dx;
	 
		Aspl(3*j-2,3*j-4)=3*pow(Dx,2);
		Aspl(3*j-2,3*j-3)=2*Dx;
		Aspl(3*j-2,3*j-2)=1;
		Aspl(3*j-2,3*j+1)=-1;
	 
		Aspl(3*j-1,3*j-4)=6*Dx;
		Aspl(3*j-1,3*j-3)=2;
		Aspl(3*j-1,3*j)=-2;
	}
	
	j=Nc;
	Dx=x(j)-x(j-1);
	Aspl(3*j-3,3*j-4)=pow(Dx,3);
	Aspl(3*j-3,3*j-3)=pow(Dx,2);
	Aspl(3*j-3,3*j-2)=Dx;
	
	Aspl(3*j-2,3*j-4)=3*pow(Dx,2);
	Aspl(3*j-2,3*j-3)=2*Dx;
	Aspl(3*j-2,3*j-2)=1;
	
	mat D(NSc,NSc);
	D.eye();
	mat Cc=solve(Aspl,D);
	
	mat T=zeros<mat>(NSc,Nx);
	for (j=0; j<Nc; j++)
	{
		T(3*j,j)=-1;
		T(3*j,j+1)=1;
	}
	T(3*Nc-2,ind_xlims(0)-1)=-1/(x(ind_xlims(0)+1)-x(ind_xlims(0)-1));
	T(3*Nc-2,ind_xlims(0)+1)=1/(x(ind_xlims(0)+1)-x(ind_xlims(0)-1));
	
	Cc=Cc*T;
	
	mat Cc2=zeros<mat>(NSc+1,Nx);
	
	Cc2.rows(0,1)=Cc.rows(0,1);
	Cc2.rows(3,NSc)=Cc.rows(2,NSc-1);
	
	//solve the spline in the interval [wr,inf
	int NSd=3*Nd-1;
	Aspl.zeros(NSd,NSd);
	
	double u=1/(x(ind_xlims(0))-x0d)-1/(x(ind_xlims(0)+1)-x0d);
	double ud=1/(x(ind_xlims(0))-x0d);
	
	Aspl(0,0)=-3*pow(ud,2)*pow(u,2);
	Aspl(0,1)=-2*pow(ud,2)*u;
	Aspl(0,2)=-pow(ud,2);
	
	Aspl(1,0)=pow(u,3);
	Aspl(1,1)=pow(u,2);
	Aspl(1,2)=u;
	
	for (j=1; j<Nd-1; j++)
	{
		u=1/(x(ind_xlims(0)+j)-x0d)-1/(x(ind_xlims(0)+j+1)-x0d);
		
		Aspl(3*j-1,3*j-1)=-1;
		Aspl(3*j-1,3*j)=3*pow(u,2);
		Aspl(3*j-1,3*j+1)=2*u;
		Aspl(3*j-1,3*j+2)=1;
		
		Aspl(3*j,3*j-2)=-2;
		Aspl(3*j,3*j)=6*u;
		Aspl(3*j,3*j+1)=2;
	 
		Aspl(3*j+1,3*j)=pow(u,3);
		Aspl(3*j+1,3*j+1)=pow(u,2);
		Aspl(3*j+1,3*j+2)=u;
	}
	
	j=Nd-1;
	u=1/(x(ind_xlims(0)+j)-x0d);
	
	Aspl(3*j-1,3*j-1)=-1;
	Aspl(3*j-1,3*j)=3*pow(u,2);
	Aspl(3*j-1,3*j+1)=2*u;
	
	Aspl(3*j,3*j-2)=-2;
	Aspl(3*j,3*j)=6*u;
	Aspl(3*j,3*j+1)=2;
	
	Aspl(3*j+1,3*j)=pow(u,3);
	Aspl(3*j+1,3*j+1)=pow(u,2);
	
	D.eye(NSd,NSd);
	mat Cd=solve(Aspl,D);
	
	T.zeros(NSd,Nx);
	T(0,ind_xlims(0)-1)=-1/(x(ind_xlims(0)+1)-x(ind_xlims(0)-1));
	T(0,ind_xlims(0)+1)=1/(x(ind_xlims(0)+1)-x(ind_xlims(0)-1));
	for (j=1; j<Nd; j++)
	{
		T(3*j-2,ind_xlims(0)+j-1)=1;
		T(3*j-2,ind_xlims(0)+j)=-1;
	}
	j=Nd;
	T(3*j-2,ind_xlims(0)+j-1)=1;
	
	Cd=Cd*T;
	
	int NS_tot=NSc+NSd+1;
	
	Mspl.zeros(NS_tot,Nx);
	Mspl.rows(0,NSc)=Cc2;
	Mspl.rows(NSc+1,NS_tot-1)=Cd;
	
	return true;
}







bool OmegaMaxEnt_data::default_model_val_chi(vec x, vec x0, vec coeffs, vec gaussians_params, vec &dm)
{
	int Nx=x.n_rows;
	int Nx0=x0.n_rows;
	
	double wcd=gaussians_params(0);
	double C1d=gaussians_params(1);
	double C2d=gaussians_params(2);
	
	if (x0(0)!=0)
	{
		cout<<"default_model_val_chi(): the first value of x0 must be 0\n";
		return false;
	}
	
	vec diff_x=x.rows(1,Nx-1)-x.rows(0,Nx-2);
	if (diff_x.min()<=0)
	{
		cout<<"default_model_val_chi(): values in position vector are not strictly increasing\n";
		return false;
	}
	
	dm.zeros(Nx);
	
	int j=0;
	while (x(j)<x0(Nx0-1) && j<Nx-1) j++;
	int jr=j-1;
	int Nm=jr+1;
	int Nr=Nx-Nm;
	
	if (Nm)
	{
		vec vm;
		spline_val(x.rows(0,jr), x0, coeffs, vm);
		dm.rows(0,jr)=vm;
	}
	if (Nr)
	{
		vec vr=exp(-pow(x.rows(jr+1,Nx-1)-wcd,2)/C1d)/C2d;
		dm.rows(jr+1,Nx-1)=vr;
	}
	
	return true;
}
























bool OmegaMaxEnt_data::compute_moments_tau_fermions()
{
	//cout<<"COMPUTING MOMENTS with compute_moments_tau_fermions()\n";
	cout<<"COMPUTING MOMENTS\n";
	
	int Nv=1;
	int NvN=1;
	int npmin=2;
	int npmax=5;
	int Np=npmax-npmin+1;
	int DNfitmin=0;
	int DNfitmax=Ntau-npmax-1;
	if (16-npmax-1<DNfitmax) DNfitmax=16-npmax-1;
	int NDN=DNfitmax-DNfitmin+1;
	
	mat M0tmp=zeros<mat>(NDN,Np);
	mat M1tmp=zeros<mat>(NDN,Np);
	mat M2tmp=zeros<mat>(NDN,Np);
	
	vec Gtmp, pp;
	int DNfit, np, Nfit;
	for (DNfit=DNfitmin; DNfit<=DNfitmax; DNfit++)
	{
		for (np=npmin; np<=npmax; np++)
		{
			Nfit=np+1+DNfit;
			Gtmp=Gtau.rows(0,Nfit-1)+flipud(Gtau.rows(Ntau-Nfit+1,Ntau));
			polyfit(tau.rows(0,Nfit-1),Gtmp,np,0,pp);
			M0tmp(DNfit-DNfitmin,np-npmin)=-pp(np);
			M2tmp(DNfit-DNfitmin,np-npmin)=-2*pp(np-2);
			
			Gtmp=Gtau.rows(0,Nfit-1)-flipud(Gtau.rows(Ntau-Nfit+1,Ntau));
			polyfit(tau.rows(0,Nfit-1),Gtmp,np,0,pp);
			M1tmp(DNfit-DNfitmin,np-npmin)=pp(np-1);
		}
	}
	
	//cout<<"1\n";
	
	mat M0m=zeros<mat>(NDN-2*NvN,Np-2*Nv);
	mat M1m=zeros<mat>(NDN-2*NvN,Np-2*Nv);
	mat M2m=zeros<mat>(NDN-2*NvN,Np-2*Nv);
	mat varM0=zeros<mat>(NDN-2*NvN,Np-2*Nv);
	mat varM1=zeros<mat>(NDN-2*NvN,Np-2*Nv);
	mat varM2=zeros<mat>(NDN-2*NvN,Np-2*Nv);
	int j, l;
	for (l=Nv; l<Np-Nv; l++)
	{
		for (j=NvN; j<NDN-NvN; j++)
		{
			M0m(j-NvN,l-Nv)=accu(M0tmp.submat(j-NvN,l-Nv,j+NvN,l+Nv))/((2*Nv+1)*(2*NvN+1));
			M1m(j-NvN,l-Nv)=accu(M1tmp.submat(j-NvN,l-Nv,j+NvN,l+Nv))/((2*Nv+1)*(2*NvN+1));
			M2m(j-NvN,l-Nv)=accu(M2tmp.submat(j-NvN,l-Nv,j+NvN,l+Nv))/((2*Nv+1)*(2*NvN+1));
			varM0(j-NvN,l-Nv)=accu(pow(M0tmp.submat(j-NvN,l-Nv,j+NvN,l+Nv)-M0m(j-NvN,l-Nv),2))/((2*Nv+1)*(2*NvN+1));
			varM1(j-NvN,l-Nv)=accu(pow(M1tmp.submat(j-NvN,l-Nv,j+NvN,l+Nv)-M1m(j-NvN,l-Nv),2))/((2*Nv+1)*(2*NvN+1));
			varM2(j-NvN,l-Nv)=accu(pow(M2tmp.submat(j-NvN,l-Nv,j+NvN,l+Nv)-M2m(j-NvN,l-Nv),2))/((2*Nv+1)*(2*NvN+1));
		}
	}
	
	uword jvmin, lvmin;
	varM0.min(jvmin,lvmin);
	double M0_NP_tmp=M0m(jvmin,lvmin);
	varM1.min(jvmin,lvmin);
	double M1_NP_tmp=M1m(jvmin,lvmin);
	varM2.min(jvmin,lvmin);
	double M2_NP_tmp=M2m(jvmin,lvmin);
	
	double Wtmp;
	if (M2_NP_tmp/M0_NP_tmp > pow(M1_NP_tmp/M0_NP_tmp,2))
		Wtmp=sqrt(M2_NP_tmp/M0_NP_tmp-pow(M1_NP_tmp/M0_NP_tmp,2));
	else
		Wtmp=abs(M1_NP_tmp/M0_NP_tmp);
	

	int Nfitmax=ceil(FNfitTauW*Ntau*tem/(abs(M1_NP_tmp/M0_NP_tmp)+Wtmp));
	if (Nfitmax>Ntau) Nfitmax=Ntau;
	
	//cout<<"2\n";

	mat X, CG, invCG, AM;
	vec Gchi2tmp, BM, Mtmp;
	npmin=3;
	int Nfitmin=npmin+1;
	int NNfit=Nfitmax-Nfitmin+1;
	mat M0b=zeros<mat>(NNfit,NNfit);
	mat M1b=zeros<mat>(NNfit,NNfit);
	mat M2b=zeros<mat>(NNfit,NNfit);
	mat M3b=zeros<mat>(NNfit,NNfit);
	int p;
	for (Nfit=Nfitmin; Nfit<=Nfitmax; Nfit++)
	{
		for (np=npmin; np<Nfit; np++)
		{
			X=zeros<mat>(Nfit,np+1);
			for (p=0; p<=np; p++)
				X.col(p)=pow(tau.rows(0,Nfit-1),p);
			
			Gchi2tmp=Gtau.rows(0,Nfit-1)+flipud(Gtau.rows(Ntau-Nfit+1,Ntau));
			CG=Ctau_all.submat(0,0,Nfit-1,Nfit-1)+fliplr(Ctau_all.submat(0,Ntau-Nfit+1,Nfit-1,Ntau))+flipud(Ctau_all.submat(Ntau-Nfit+1,0,Ntau,Nfit-1))+flipud(fliplr(Ctau_all.submat(Ntau-Nfit+1,Ntau-Nfit+1,Ntau,Ntau)));
			CG(0,0)=CG(1,1);
			invCG=inv(CG);
			//	invCG=inv_sympd(CG);
			AM=(X.t())*invCG*X;
			BM=(X.t())*invCG*Gchi2tmp;
			Mtmp=solve(AM,BM);
			M0b(Nfit-np-1,np-npmin)=-Mtmp(0);
			M2b(Nfit-np-1,np-npmin)=-2*Mtmp(2);
		 
		 	Gchi2tmp=Gtau.rows(0,Nfit-1)-flipud(Gtau.rows(Ntau-Nfit+1,Ntau));
		 	CG=Ctau_all.submat(0,0,Nfit-1,Nfit-1)- fliplr(Ctau_all.submat(0,Ntau-Nfit+1,Nfit-1,Ntau))-flipud(Ctau_all.submat(Ntau-Nfit+1,0,Ntau,Nfit-1))+flipud(fliplr(Ctau_all.submat(Ntau-Nfit+1,Ntau-Nfit+1,Ntau,Ntau)));
		 	invCG=inv(CG);
			//	invCG=inv_sympd(CG);
		 	AM=(X.t())*invCG*X;
			BM=(X.t())*invCG*Gchi2tmp;
			Mtmp=solve(AM,BM);
			M1b(Nfit-np-1,np-npmin)=Mtmp(1);
			M3b(Nfit-np-1,np-npmin)=6*Mtmp(3);
 
		}
	}

	//cout<<"3\n";

	Nv=1;
	NvN=1;
	int jmin=NvN;
	int jmax=NNfit-NvN-2*Nv-1;
	int Nj=jmax-jmin+1;
	int lmin=Nv;
	int lmax=NNfit-2*NvN-Nv-1;
	int Nl=lmax-lmin+1;
	M0m.zeros(Nj,Nl);
	M1m.zeros(Nj,Nl);
	M2m.zeros(Nj,Nl);
	mat M3m=zeros<mat>(Nj,Nl);
	varM0.zeros(Nj,Nl);
	varM1.zeros(Nj,Nl);
	varM2.zeros(Nj,Nl);
	mat varM3=zeros<mat>(Nj,Nl);
	for (j=jmin; j<=jmax; j++)
	{
		for (l=lmin; l<=NNfit-1-j-Nv-NvN; l++)
		{
			M0m(j-jmin,l-lmin)=accu(M0b.submat(j-NvN,l-Nv,j+NvN,l+Nv))/((2*Nv+1)*(2*NvN+1));
			M1m(j-jmin,l-lmin)=accu(M1b.submat(j-NvN,l-Nv,j+NvN,l+Nv))/((2*Nv+1)*(2*NvN+1));
			M2m(j-jmin,l-lmin)=accu(M2b.submat(j-NvN,l-Nv,j+NvN,l+Nv))/((2*Nv+1)*(2*NvN+1));
			M3m(j-jmin,l-lmin)=accu(M3b.submat(j-NvN,l-Nv,j+NvN,l+Nv))/((2*Nv+1)*(2*NvN+1));
			varM0(j-jmin,l-lmin)=accu(pow(M0b.submat(j-NvN,l-Nv,j+NvN,l+Nv)-M0m(j-jmin,l-lmin),2))/((2*Nv+1)*(2*NvN+1));
			varM1(j-jmin,l-lmin)=accu(pow(M1b.submat(j-NvN,l-Nv,j+NvN,l+Nv)-M1m(j-jmin,l-lmin),2))/((2*Nv+1)*(2*NvN+1));
			varM2(j-jmin,l-lmin)=accu(pow(M2b.submat(j-NvN,l-Nv,j+NvN,l+Nv)-M2m(j-jmin,l-lmin),2))/((2*Nv+1)*(2*NvN+1));
			varM3(j-jmin,l-lmin)=accu(pow(M3b.submat(j-NvN,l-Nv,j+NvN,l+Nv)-M3m(j-jmin,l-lmin),2))/((2*Nv+1)*(2*NvN+1));
		}
	}
	
//	cout<<varM0<<endl;

	double varM0max=max(max(varM0));
	double varM1max=max(max(varM1));
	double varM2max=max(max(varM2));
	double varM3max=max(max(varM3));

	for (j=1; j<Nj; j++)
	{
		varM0.submat(j,Nl-j,j,Nl-1)=2*varM0max*ones<rowvec>(j);
		varM1.submat(j,Nl-j,j,Nl-1)=2*varM1max*ones<rowvec>(j);
		varM2.submat(j,Nl-j,j,Nl-1)=2*varM2max*ones<rowvec>(j);
		varM3.submat(j,Nl-j,j,Nl-1)=2*varM3max*ones<rowvec>(j);
	}
	
//	cout<<varM0<<endl;
	
	//cout<<"4\n";

	varM0.min(jvmin,lvmin);
	double M0_N=M0m(jvmin,lvmin);
	varM1.min(jvmin,lvmin);
	double M1_N=M1m(jvmin,lvmin);
	varM2.min(jvmin,lvmin);
	double M2_N=M2m(jvmin,lvmin);
	varM3.min(jvmin,lvmin);
	double M3_N=M3m(jvmin,lvmin);
	
	cout<<"moments determined by polynomial fit to G(tau) at boundaries:\n";
	cout<<"norm: "<<M0_N<<endl;
	cout<<"first moment: "<<M1_N<<endl;
	cout<<"second moment: "<<M2_N<<endl;
	cout<<"third moment: "<<M3_N<<endl;
	
	varM1.min(jvmin,lvmin);
	np=lvmin+Nv+2;
	Nfit=jvmin+NvN+np+1;
	X=zeros<mat>(Nfit,np+1);
	for (p=0; p<=np; p++)
		X.col(p)=pow(tau.rows(0,Nfit-1),p);
	
	CG=Ctau_all.submat(0,0,Nfit-1,Nfit-1)+fliplr(Ctau_all.submat(0,Ntau-Nfit+1,Nfit-1,Ntau))+flipud(Ctau_all.submat(Ntau-Nfit+1,0,Ntau,Nfit-1))+flipud(fliplr(Ctau_all.submat(Ntau-Nfit+1,Ntau-Nfit+1,Ntau,Ntau)));
	CG(0,0)=CG(1,1);
	invCG=inv(CG);
	//	invCG=inv_sympd(CG);
	AM=(X.t())*invCG*X;
	mat invAMp=inv(AM);
	
	CG=Ctau_all.submat(0,0,Nfit-1,Nfit-1)-fliplr(Ctau_all.submat(0,Ntau-Nfit+1,Nfit-1,Ntau))-flipud(Ctau_all.submat(Ntau-Nfit+1,0,Ntau,Nfit-1))+flipud(fliplr(Ctau_all.submat(Ntau-Nfit+1,Ntau-Nfit+1,Ntau,Ntau)));
	invCG=inv(CG);
	//	invCG=inv_sympd(CG);
	AM=(X.t())*invCG*X;
	mat invAMn=inv(AM);
	
	double std_omega_tmp;
	double var_omega=M2_N/M0_N-pow(M1_N/M0_N,2);
	if (var_omega>0)
		std_omega_tmp=sqrt(var_omega);
	else
	{
		cout<<"Negative variance found during computation of moments.\n";
		return false;
	}
	
	covm_diag=true;
	if (!moments_provided)
	{
		M.zeros(4);
		M(0)=M0;
		M(1)=M1_N;
		M(2)=M2_N;
		M(3)=M3_N;
		M1=M(1);
		M2=M(2);
		M3=M(3);
		NM=4;
		M1_set=true;
		M2_set=true;
		errM.zeros(NM);
		errM(0)=errM0;
		errM(1)=sqrt(invAMn(1,1));
		errM(2)=sqrt(invAMp(2,2));
		errM(3)=sqrt(invAMn(3,3));
		errM1=errM(1);
		errM2=errM(2);
		errM3=errM(3);
		M_ord=linspace<vec>(0,NM-1,NM);
	}
	else
	{
		if (M0)
		{
			if (abs(M0_N-M0)/M0_N>tol_norm)
			{
				if (M0_in.size())
					cout<<"warning: norm of spectral function is different from provided one.\n";
				else
				{
					cout<<"warning: spectral function is not normalized.\n";
					cout<<"Use parameter \"norm of spectral function:\" in subsection DATA PARAMETERS to provide a norm different from 1.\n";
				}
			}
		}
		
		if (abs(M1-M1_N)/std_omega_tmp>tol_M1)
			cout<<"warning: first moment different from provided one\n";
		
		if (M2_in.size())
		{
			if (abs(M2-M2_N)/M2_N>tol_M2)
				cout<<"warning: second moment different from provided one\n";
			
			if (M3_in.size())
			{
				if (abs(M3-M3_N)/pow(std_omega_tmp,3)>tol_M3)
				{
					cout<<"warning: third moment different from provided one\n";
				}
			}
			else
			{
				M3=M3_N;
				M.zeros(4);
				M(0)=M0;
				M(1)=M1;
				M(2)=M2;
				M(3)=M3;
				errM3=sqrt(invAMn(3,3));
				NM=4;
				errM.zeros(NM);
				errM(0)=errM0;
				errM(1)=errM1;
				errM(2)=errM2;
				errM(3)=errM3;
				M_ord=linspace<vec>(0,NM-1,NM);
			}
		}
		else if (M3_in.size())
		{
			M2=M2_N;
			M2_set=true;
			M.zeros(4);
			M(0)=M0;
			M(1)=M1;
			M(2)=M2;
			M(3)=M3;
			errM2=sqrt(invAMp(2,2));
			NM=4;
			errM.zeros(NM);
			errM(0)=errM0;
			errM(1)=errM1;
			errM(2)=errM2;
			errM(3)=errM3;
			M_ord=linspace<vec>(0,NM-1,NM);
		}
		else
		{
			M2=M2_N;
			M2_set=true;
			M3=M3_N;
			errM2=sqrt(invAMp(2,2));
			errM3=sqrt(invAMn(3,3));
			M.zeros(4);
			M(0)=M0;
			M(1)=M1;
			M(2)=M2;
			M(3)=M3;
			NM=4;
			errM.zeros(NM);
			errM(0)=errM0;
			errM(1)=errM1;
			errM(2)=errM2;
			errM(3)=errM3;
			M_ord=linspace<vec>(0,NM-1,NM);
		}
	}
	COVM.zeros(NM,NM);
	COVM.diag()=square(errM);
	
	if (!std_omega)
	{
		var_omega=M2/M0-pow(M1/M0,2);
		std_omega=sqrt(var_omega);
	}
	
	if (!SC_set)
	{
		SC=M1/M0;
		SC_set=true;
	}
	if (!SW_set)
	{
		SW=f_SW_std_omega*std_omega;
		SW_set=true;
	}
	if (M(2)<0)
	{
		M=M.rows(0,1);
		COVM=COVM.submat(0,0,1,1);
		NM=2;
	}
	
	return true;
}







bool OmegaMaxEnt_data::fit_circle_arc(vec x, vec y, vec &arc_params)
{
	
	arc_params.set_size(3);
	arc_params(0)=INF;
	arc_params(1)=NAN;
	arc_params(2)=NAN;
	
	double tol_r=1e-10;
	int Niter_max=20;
	
	int N=x.n_rows;
	
	if (N<3)
	{
		cout<<"fit_circle_arc(): number of points must be larger of equal to 3\n";
		return false;
	}
	
	mat X(N,3);
	X.zeros();
	
	X.col(0)=pow(x-x(0),2);
	X.col(1)=x-x(0);
	X.col(2)=ones<vec>(N);
	
	mat A=X.t()*X;
	vec B=X.t()*y;
	
	vec cfs=solve(A,B);
	
	int jmid=floor(N/2);
	
	double Dx;
	double x0=x(jmid-1);
	Dx=x0-x(0);
	double y0=cfs(0)*pow(Dx,2)+cfs(1)*Dx+cfs(2);
	double x1=x(jmid);
	Dx=x1-x(0);
	double y1=cfs(0)*pow(Dx,2)+cfs(1)*Dx+cfs(2);
	double x2=x(jmid+1);
	Dx=x2-x(0);
	double y2=cfs(0)*pow(Dx,2)+cfs(1)*Dx+cfs(2);
	
	double denom=2*(x2*(y1-y0) + x1*(y0-y2) + x0*(y2-y1));
	
	if (denom==0)
	{
		arc_params(0)=INF;
		arc_params(1)=NAN;
		arc_params(2)=NAN;
		return true;
	}
	
	double xc_init=(pow(x2,2)*(y1-y0) + pow(x1,2)*(y0-y2) + pow(x0,2)*(y2-y1) + (y1-y0)*(y0-y2)*(y1-y2))/denom;
	double yc_init=-(pow(y2,2)*(x1-x0) + pow(y1,2)*(x0-x2) + pow(y0,2)*(x2-x1) + (x1-x0)*(x2-x0)*(x2-x1))/denom;
	
//	double r2_init=pow(x1-xc_init,2)+pow(y1-yc_init,2);
//	double r_init=sqrt(r2_init);
	
	double xc=xc_init;
	double yc=yc_init;
	
	vec rj=sqrt(pow(x-xc,2)+pow(y-yc,2));
	double r=sum(rj)/N;
	vec drj=r-rj;
	
	vec dr2=pow(drj,2);
	double chi2=sum(dr2);
	
	double xc_m=2*x1-xc;
	double yc_m=2*y1-yc;
	
	vec rjm=sqrt(pow(x-xc_m,2)+pow(y-yc_m,2));
	double rm=sum(rjm)/N;
	vec drjm=rm-rjm;
	dr2=pow(drjm,2);
	double chi2m=sum(dr2);
	
	if (chi2m<chi2)
	{
		xc=xc_m;
		yc=yc_m;
		rj=rjm;
		r=rm;
		drj=drjm;
		chi2=chi2m;
	}
	
	vec Drjxc=(xc-x)/rj;
	vec Drjyc=(yc-y)/rj;
	double Drxc=sum(Drjxc)/N;
	double Dryc=sum(Drjyc)/N;
	
	vec D2rjxc=1.0/rj-pow(x-xc,2)/pow(rj,3);
	vec D2rjyc=1.0/rj-pow(y-yc,2)/pow(rj,3);
	vec D2rjxcyc=-((x-xc) % (y-yc))/pow(rj,3);
	double D2rxc=sum(D2rjxc)/N;
	double D2ryc=sum(D2rjyc)/N;
	double D2rxcyc=sum(D2rjxcyc)/N;
	
	double Dchi2xc=sum(drj % (Drxc-Drjxc));
	double Dchi2yc=sum(drj % (Dryc-Drjyc));
	
	vec gradchi2(2);
	gradchi2(0)=-Dchi2xc;
	gradchi2(1)=-Dchi2yc;
	
	double D2chi2xc=sum( pow(Drxc-Drjxc,2) + drj % (D2rxc-D2rjxc) );
	double D2chi2xcyc=sum( (Drxc-Drjxc) % (Dryc-Drjyc) + drj % (D2rxcyc-D2rjxcyc) );
	double D2chi2yc=sum( pow(Dryc-Drjyc,2) + drj % (D2ryc-D2rjyc) );
	
	vec chi2_vec_arc(Niter_max), r_vec(Niter_max), xc_vec(Niter_max), yc_vec(Niter_max);
	chi2_vec_arc(0)=chi2;
	r_vec(0)=r;
	xc_vec(0)=xc;
	yc_vec(0)=yc;
	
	mat H(2,2);
	H(0,0)=D2chi2xc;
	H(0,1)=D2chi2xcyc;
	H(1,0)=D2chi2xcyc;
	H(1,1)=D2chi2yc;
	
	vec p=solve(H,gradchi2);
	
	
	double xc_prec=xc;
	double yc_prec=yc;
	double r_prec=r;
	
	xc=xc+p(0);
	yc=yc+p(1);
	
	rj=sqrt(pow(x-xc,2)+pow(y-yc,2));
	r=sum(rj)/N;
	drj=r-rj;
	
	double chi2_prec=chi2;
	dr2=pow(drj,2);
	chi2=sum(dr2);
	
	chi2_vec_arc(1)=chi2;
	r_vec(1)=r;
	xc_vec(1)=xc;
	yc_vec(1)=yc;
	
	int j=1;
	while (abs(r-r_prec)/r>tol_r && j<Niter_max-1 && chi2<chi2_prec)
	{
		Drjxc=(xc-x)/rj;
		Drjyc=(yc-y)/rj;
		Drxc=sum(Drjxc)/N;
		Dryc=sum(Drjyc)/N;
		
		D2rjxc=1.0/rj-pow(x-xc,2)/pow(rj,3);
		D2rjyc=1.0/rj-pow(y-yc,2)/pow(rj,3);
		D2rjxcyc=-((x-xc) % (y-yc))/pow(rj,3);
		D2rxc=sum(D2rjxc)/N;
		D2ryc=sum(D2rjyc)/N;
		D2rxcyc=sum(D2rjxcyc)/N;
		
		Dchi2xc=sum(drj % (Drxc-Drjxc));
		Dchi2yc=sum(drj % (Dryc-Drjyc));
		
		gradchi2(0)=-Dchi2xc;
		gradchi2(1)=-Dchi2yc;
		
		D2chi2xc=sum( pow(Drxc-Drjxc,2) + drj % (D2rxc-D2rjxc) );
		D2chi2xcyc=sum( (Drxc-Drjxc) % (Dryc-Drjyc) + drj % (D2rxcyc-D2rjxcyc) );
		D2chi2yc=sum( pow(Dryc-Drjyc,2) + drj % (D2ryc-D2rjyc) );
		
		H(0,0)=D2chi2xc;
		H(0,1)=D2chi2xcyc;
		H(1,0)=D2chi2xcyc;
		H(1,1)=D2chi2yc;
		
		p=solve(H,gradchi2);
		
		xc_prec=xc;
		yc_prec=yc;
	 	r_prec=r;
	 
	 	xc=xc+p(0);
	 	yc=yc+p(1);
	 
		rj=sqrt(pow(x-xc,2)+pow(y-yc,2));
		r=sum(rj)/N;
		drj=r-rj;
		
		chi2_prec=chi2;
		dr2=pow(drj,2);
		chi2=sum(dr2);
	 
	 	xc_m=2*x1-xc;
	 	yc_m=2*y1-yc;
	 
	 	rjm=sqrt(pow(x-xc_m,2)+pow(y-yc_m,2));
	 	rm=sum(rjm)/N;
	 	drjm=rm-rjm;
	 
		dr2=pow(drjm,2);
	 	chi2m=sum(dr2);
	 
		if (chi2m<chi2)
		{
			xc=xc_m;
			yc=yc_m;
			rj=rjm;
			r=rm;
			drj=drjm;
			chi2=chi2m;
		}
	 
		j++;
		
		chi2_vec_arc(j)=chi2;
		r_vec(j)=r;
		xc_vec(j)=xc;
		yc_vec(j)=yc;
	}
	//int Nvec=j+1;
	
	if (chi2>chi2_prec)
	{
		xc=xc_prec;
		yc=yc_prec;
		r=r_prec;
	}
	
	double r2=pow(r,2);
	vec ypl=zeros<vec>(N);
	vec dr2xc2=r2-pow(x-xc,2);
	vec xpl=x;
	double sgn=0;
	for (j=0; j<N; j++)
	{
		if (dr2xc2(j)<0)
		{
		 	if (x(j)>xc)
				xpl(j)=xc+r;
		 	else
				xpl(j)=xc-r;
		}
		if (y(j)<yc)
		{
			ypl(j)=yc-sqrt(r2-pow(xpl(j)-xc,2));
			sgn=sgn+1;
		}
		else
		{
			ypl(j)=yc+sqrt(r2-pow(xpl(j)-xc,2));
			if (y(j)>yc)
				sgn=sgn-1;
		}
	}
	
	if (sgn) r=sgn*r/abs(sgn);
	
	arc_params(0)=r;
	arc_params(1)=xc;
	arc_params(2)=yc;
	
	return true;
}






bool OmegaMaxEnt_data::Kernel_G_fermions_grid_transf_omega()
{
	bool use_HF_exp=true;
	double fg=1.7;
	//int pngmax=100;
	double fi=fg;
	double fr=fi;
	int pnmax=50;
	//double fd=fg;
	//int pndmax=pngmax;
	
	cout<<"defining kernel matrix...\n";
	
	Kcx.zeros(Nn,Nw);
	KM.zeros(4,Nw);
	
	mat MM;
	spline_matrix_grid_transf(w, MM);
	
	
	dcomplex i(0,1);
	int Nint=Nw-1;
	
	cx_mat Ka=zeros<cx_mat>(Nn,Nint);
	cx_mat Kb=zeros<cx_mat>(Nn,Nint);
	cx_mat Kc=zeros<cx_mat>(Nn,Nint);
	cx_mat Kd=zeros<cx_mat>(Nn,Nint);
	
	mat Wn=wn*ones<rowvec>(Nint);
	mat W=ones<vec>(Nn)*w.t();
	vec vdw=w.rows(1,Nw-1)-w.rows(0,Nw-2);
	
	mat logc=log((pow(Wn,2)+pow(W.cols(1,Nint),2))/(pow(Wn,2)+pow(W.cols(0,Nint-1),2)));
	mat atanc=atan((Wn % (W.cols(1,Nint)-W.cols(0,Nint-1)))/(W.cols(1,Nint) % W.cols(0,Nint-1)+pow(Wn,2)));
	cx_mat logc2=logc/2+i*atanc;
	cx_mat dWn=i*Wn-W.cols(0,Nint-1);
	mat dW=W.cols(1,Nint)-W.cols(0,Nint-1);
	cx_mat Rdwn_dw=dWn/dW;
	
	Ka=-pow(Rdwn_dw,2)-Rdwn_dw/2-1.0/3.0-pow(Rdwn_dw,3) % logc2;
	Kb=-Rdwn_dw-0.5-pow(Rdwn_dw,2) % logc2;
	Kc=-1.0 - Rdwn_dw % logc2;
	Kd=-logc2;
	
	int Pmax=2*pnmax+4;
	imat MP;
	pascal(Pmax+1,MP);
	
	int j;
	if (use_HF_exp)
	{
		double wtmp;
		int jni, jnr, p, l;
		for (j=0; j<Nw-1; j++)
		{
			wtmp=abs(w(j));
			if (abs(w(j+1))>wtmp) wtmp=abs(w(j+1));
			jni=0;
			while (wn(jni)<fi*wtmp && jni<Nn-1) jni++;
			while (pow(wn(jni),2*pnmax+1)==0 && jni<Nn-1) jni++;
			jnr=1;
			while (wn(jnr)<fr*wtmp && jnr<Nn-1) jnr++;
			while (pow(wn(jnr),2*pnmax)==0 && jnr<Nn-1) jnr++;
			if (jni<Nn-1 || jnr<Nn-1)
			{
				double wj=w(j);
				double dw1=w(j+1)-wj;
				vec dwp=zeros<vec>(Pmax);
				dwp(0)=dw1;
				dwp(1)=pow(dw1,2)+2*wj*dw1;
				for (p=3; p<=Pmax; p++)
				{
					dwp(p-1)=0;
					for (l=0; l<p; l++)
					{
						dwp(p-1)=dwp(p-1)+MP(l,p-l)*pow(dw1,p-l)*pow(wj,l);
					}
				}
				cx_vec vtmp;
				if (jni<Nn)
				{
					vtmp.zeros(Nn-jni);
					vtmp.set_real(real(Ka.submat(jni,j,Nn-1,j)));
					Ka.submat(jni,j,Nn-1,j)=vtmp;
					vtmp.set_real(real(Kb.submat(jni,j,Nn-1,j)));
					Kb.submat(jni,j,Nn-1,j)=vtmp;
					vtmp.set_real(real(Kc.submat(jni,j,Nn-1,j)));
					Kc.submat(jni,j,Nn-1,j)=vtmp;
					vtmp.set_real(real(Kd.submat(jni,j,Nn-1,j)));
					Kd.submat(jni,j,Nn-1,j)=vtmp;
					
					for (p=2*pnmax+1; p>=1; p-=2)
					{
						Ka.submat(jni,j,Nn-1,j)=Ka.submat(jni,j,Nn-1,j) + (i*pow(-1,(p-1)/2)*(pow(wj,3)*dwp(p-1)/p - 3.0*pow(wj,2)*dwp(p)/(p+1) + 3.0*wj*dwp(p+1)/(p+2) - dwp(p+2)/(p+3))/pow(wn.rows(jni,Nn-1),p))/pow(dw1,3);
						Kb.submat(jni,j,Nn-1,j)=Kb.submat(jni,j,Nn-1,j) + (i*pow(-1,(p-1)/2)*(-pow(wj,2)*dwp(p-1)/p + 2*wj*dwp(p)/(p+1) - dwp(p+1)/(p+2))/pow(wn.rows(jni,Nn-1),p))/pow(dw1,2);
						Kc.submat(jni,j,Nn-1,j)=Kc.submat(jni,j,Nn-1,j) + (i*pow(-1,(p-1)/2)*(wj*dwp(p-1)/p - dwp(p)/(p+1))/pow(wn.rows(jni,Nn-1),p))/dw1;
						Kd.submat(jni,j,Nn-1,j)=Kd.submat(jni,j,Nn-1,j) - i*pow(-1,(p-1)/2)*dwp(p-1)/(p*pow(wn.rows(jni,Nn-1),p));
					}
				}
				if (jnr<Nn)
				{
					vtmp.zeros(Nn-jnr);
					vtmp.set_imag(imag(Ka.submat(jnr,j,Nn-1,j)));
					Ka.submat(jnr,j,Nn-1,j)=vtmp;
					vtmp.set_imag(imag(Kb.submat(jnr,j,Nn-1,j)));
					Kb.submat(jnr,j,Nn-1,j)=vtmp;
					vtmp.set_imag(imag(Kc.submat(jnr,j,Nn-1,j)));
					Kc.submat(jnr,j,Nn-1,j)=vtmp;
					vtmp.set_imag(imag(Kd.submat(jnr,j,Nn-1,j)));
					Kd.submat(jnr,j,Nn-1,j)=vtmp;
					for (p=2*pnmax; p>=2; p-=2)
					{
						Ka.submat(jnr,j,Nn-1,j)=Ka.submat(jnr,j,Nn-1,j) + (pow(-1,(p-2)/2)*(pow(wj,3)*dwp(p-1)/p - 3*pow(wj,2)*dwp(p)/(p+1) + 3*wj*dwp(p+1)/(p+2) - dwp(p+2)/(p+3))/pow(wn.rows(jnr,Nn-1),p))/pow(dw1,3);
						Kb.submat(jnr,j,Nn-1,j)=Kb.submat(jnr,j,Nn-1,j) + (pow(-1,(p-2)/2)*(-pow(wj,2)*dwp(p-1)/p + 2*wj*dwp(p)/(p+1) - dwp(p+1)/(p+2))/pow(wn.rows(jnr,Nn-1),p))/pow(dw1,2);
						Kc.submat(jnr,j,Nn-1,j)=Kc.submat(jnr,j,Nn-1,j) + (pow(-1,(p-2)/2)*(wj*dwp(p-1)/p - dwp(p)/(p+1))/pow(wn.rows(jnr,Nn-1),p))/dw1;
						Kd.submat(jnr,j,Nn-1,j)=Kd.submat(jnr,j,Nn-1,j) - pow(-1,(p-2)/2)*dwp(p-1)/(p*pow(wn.rows(jnr,Nn-1),p));
					}
				}
			}
		}
	}
	
	Ka=Ka/(2*PI);
	Kb=Kb/(2*PI);
	Kc=Kc/(2*PI);
	Kd=Kd/(2*PI);
	
	mat	Pa=zeros<mat>(Nint,4*Nint);
	mat	Pb=zeros<mat>(Nint,4*Nint);
	mat	Pc=zeros<mat>(Nint,4*Nint);
	mat	Pd=zeros<mat>(Nint,4*Nint);
	
	for (j=0; j<Nint; j++)
	{
		Pa(j,4*j)=1;
		Pb(j,4*j+1)=1;
		Pc(j,4*j+2)=1;
		Pd(j,4*j+3)=1;
	}
	
	Kcx=(Ka*Pa+Kb*Pb+Kc*Pc+Kd*Pd)*MM;
	
	uvec even_ind=linspace<uvec>(0,2*Nn-2,Nn);
	K.zeros(2*Nn,Nw);
	K.rows(even_ind)=real(Kcx);
	K.rows(even_ind+1)=imag(Kcx);
	
	//////////////////////
	
	rowvec Knorm_a=zeros<rowvec>(Nint);
	rowvec Knorm_b=zeros<rowvec>(Nint);
	rowvec Knorm_c=zeros<rowvec>(Nint);
	rowvec Knorm_d=zeros<rowvec>(Nint);
	
	Knorm_a=trans(pow(w.rows(1,Nw-1)-w.rows(0,Nw-2),4)/4);
	Knorm_b=trans(pow(w.rows(1,Nw-1)-w.rows(0,Nw-2),3)/3);
	Knorm_c=trans(pow(w.rows(1,Nw-1)-w.rows(0,Nw-2),2)/2);
	Knorm_d=trans(w.rows(1,Nw-1)-w.rows(0,Nw-2));
	
	rowvec KM0=(Knorm_a*Pa+Knorm_b*Pb+Knorm_c*Pc+Knorm_d*Pd)*MM/(2*PI);
	
	mat Wjc=diagmat(w.rows(0,Nint-1));
	
	rowvec KM1_a_c_tmp=trans(pow(w.rows(1,Nw-1)-w.rows(0,Nw-2),5)/5.0);
	
	rowvec KM1_a_c=KM1_a_c_tmp+Knorm_a*Wjc;
	rowvec KM1_b_c=Knorm_a+Knorm_b*Wjc;
	rowvec KM1_c_c=Knorm_b+Knorm_c*Wjc;
	rowvec KM1_d_c=Knorm_c+Knorm_d*Wjc;
	
	rowvec KM1=(KM1_a_c*Pa+KM1_b_c*Pb+KM1_c_c*Pc+KM1_d_c*Pd)*MM/(2*PI);
	
	rowvec KM2_a_c_tmp=trans(pow(w.rows(1,Nw-1)-w.rows(0,Nw-2),6)/6);
	
	rowvec KM2_a_c=KM2_a_c_tmp+2*KM1_a_c_tmp*Wjc+Knorm_a*pow(Wjc,2);
	rowvec KM2_b_c=KM1_a_c_tmp+2*Knorm_a*Wjc+Knorm_b*pow(Wjc,2);
	rowvec KM2_c_c=Knorm_a+2*Knorm_b*Wjc+Knorm_c*pow(Wjc,2);
	rowvec KM2_d_c=Knorm_b+2*Knorm_c*Wjc+Knorm_d*pow(Wjc,2);
	
	rowvec KM2=(KM2_a_c*Pa+KM2_b_c*Pb+KM2_c_c*Pc+KM2_d_c*Pd)*MM/(2*PI);
	
	rowvec KM3_a_c_tmp=trans(pow(w.rows(1,Nw-1)-w.rows(0,Nw-2),7)/7);
	
	rowvec KM3_a_c=KM3_a_c_tmp+3*KM2_a_c_tmp*Wjc+3*KM1_a_c_tmp*pow(Wjc,2)+Knorm_a*pow(Wjc,3);
	rowvec KM3_b_c=KM2_a_c_tmp+3*KM1_a_c_tmp*Wjc+3*Knorm_a*pow(Wjc,2)+Knorm_b*pow(Wjc,3);
	rowvec KM3_c_c=KM1_a_c_tmp+3*Knorm_a*Wjc+3*Knorm_b*pow(Wjc,2)+Knorm_c*pow(Wjc,3);
	rowvec KM3_d_c=Knorm_a+3*Knorm_b*Wjc+3*Knorm_c*pow(Wjc,2)+Knorm_d*pow(Wjc,3);
	
	rowvec KM3=(KM3_a_c*Pa+KM3_b_c*Pb+KM3_c_c*Pc+KM3_d_c*Pd)*MM/(2*PI);
	
	KM.row(0)=KM0;
	KM.row(1)=KM1;
	KM.row(2)=KM2;
	KM.row(3)=KM3;
	
	/*
	//test moments
	rowvec x0={-2, 0, 3};
	rowvec s0={0.8, 0.3, 1};
	rowvec wgt={1,1,1};
	
	vec test_A;
	sum_gaussians(w, x0, s0, wgt, test_A);
	
	cout<<KM*test_A<<endl;
	*/
	
	return true;
}




bool OmegaMaxEnt_data::Kernel_G_fermions_grid_transf()
{
	bool use_HF_exp=true;
	double fg=1.7;
	int pngmax=100;
	double fi=fg;
	double fr=fi;
	int pnmax=50;
	double fd=fg;
	int pndmax=pngmax;
	
	cout<<"defining kernel matrix...\n";
	
	Kcx.zeros(Nn,Nw);
 
	int Nug=Nw_lims(0);
	int Nud=Nw-Nw_lims(1)-1;
	
	vec ug, ud;
	if (Du_constant)
	{
		double dug=dwl/((wl-dwl-w0l)*(wl-w0l));
		vec ug_int=linspace<vec>(1,Nug,Nug);
		ug=-dug*ug_int;
		
		double dud=dwr/((wr-w0r)*(wr+dwr-w0r));
		vec ud_int=linspace<vec>(Nud,1,Nud);
		ud=dud*ud_int;
	}
	else
	{
		ug=1.0/(w.rows(0,Nug-1)-w0l);
		ud=1.0/(w.rows(Nw_lims(1)+1,Nw-1)-w0r);
	}
	
	mat MM;
	spline_matrix_grid_transf_G_part(w, Nw_lims, ws, MM);
	
   //////////////////////////////////
	
	int Nint=Nw-1;
	int NCfs=4*Nint;
	
	int Nintg=Nug;
	int NCg=4*Nintg;
	
	mat Pa_g=zeros<mat>(Nintg,NCfs);
	mat Pb_g_r=zeros<mat>(Nintg,NCfs);
	mat Pc_g_r=zeros<mat>(Nintg,NCfs);
	mat Pd_g_r=zeros<mat>(Nintg,NCfs);
	
	int j;
	for (j=0; j<Nintg; j++)
	{
		Pa_g(j,4*j)=1;
		Pb_g_r(j,4*j+1)=1;
		Pc_g_r(j,4*j+2)=1;
		Pd_g_r(j,4*j+3)=1;
	}
	
	vec vDg=((w.rows(1,Nw_lims(0))-w0l)%(w.rows(0,Nw_lims(0)-1)-w0l))/(w.rows(0,Nw_lims(0)-1)-w.rows(1,Nw_lims(0)));
	mat Dg=diagmat(vDg);
	vec vDU=(w.rows(1,Nw_lims(0))-w0l)/(w.rows(0,Nw_lims(0)-1)-w.rows(1,Nw_lims(0)));
	mat DU=diagmat(vDU);
	
	mat Pb_g=Pb_g_r-3*DU*Pa_g;
	mat Pc_g=Pc_g_r+3*pow(DU,2)*Pa_g-2*DU*Pb_g_r;
	mat Pd_g=Pd_g_r-pow(DU,3)*Pa_g+pow(DU,2)*Pb_g_r-DU*Pc_g_r;
	
	Pa_g=pow(Dg,3)*Pa_g;
	Pb_g=pow(Dg,2)*Pb_g;
	Pc_g=Dg*Pc_g;
	
	int	Nintc=Nwc-1;
	
	mat	Pa_c=zeros<mat>(Nintc,NCfs);
	mat	Pb_c=zeros<mat>(Nintc,NCfs);
	mat	Pc_c=zeros<mat>(Nintc,NCfs);
	mat	Pd_c=zeros<mat>(Nintc,NCfs);
	
	for (j=0; j<Nintc; j++)
	{
		Pa_c(j,4*j+NCg)=1;
		Pb_c(j,4*j+1+NCg)=1;
		Pc_c(j,4*j+2+NCg)=1;
		Pd_c(j,4*j+3+NCg)=1;
	}
	
	vec vDc=1.0/(w.rows(Nw_lims(0)+1,Nw_lims(1))-w.rows(Nw_lims(0),Nw_lims(1)-1));
	mat Dc=diagmat(vDc);
	
	Pa_c=pow(Dc,3)*Pa_c;
	Pb_c=pow(Dc,2)*Pb_c;
	Pc_c=Dc*Pc_c;
	
	int NCgc=NCg+4*Nintc;
	
	int Nintd=Nud;
	
	mat Pa_d=zeros<mat>(Nintd,NCfs);
	mat Pb_d_r=zeros<mat>(Nintd,NCfs);
	mat Pc_d_r=zeros<mat>(Nintd,NCfs);
	mat Pd_d_r=zeros<mat>(Nintd,NCfs);
	
	for (j=0; j<Nintd; j++)
	{
		Pa_d(j,4*j+NCgc)=1;
		Pb_d_r(j,4*j+1+NCgc)=1;
		Pc_d_r(j,4*j+2+NCgc)=1;
		Pd_d_r(j,4*j+3+NCgc)=1;
	}
	
	vec vDd=((w.rows(Nw_lims(1)+1,Nw-1)-w0r)%(w.rows(Nw_lims(1),Nw-2)-w0r))/(w.rows(Nw_lims(1),Nw-2)-w.rows(Nw_lims(1)+1,Nw-1));
	mat Dd=diagmat(vDd);
	vDU=(w.rows(Nw_lims(1)+1,Nw-1)-w0r)/(w.rows(Nw_lims(1),Nw-2)-w.rows(Nw_lims(1)+1,Nw-1));
	DU=diagmat(vDU);
	
	mat Pb_d=Pb_d_r-3*DU*Pa_d;
	mat Pc_d=Pc_d_r+3*pow(DU,2)*Pa_d-2*DU*Pb_d_r;
	mat Pd_d=Pd_d_r-pow(DU,3)*Pa_d+pow(DU,2)*Pb_d_r-DU*Pc_d_r;
	
	Pa_d=pow(Dd,3)*Pa_d;
	Pb_d=pow(Dd,2)*Pb_d;
	Pc_d=Dd*Pc_d;
	
	cx_mat Ka_g=zeros<cx_mat>(Nn,Nintg);
	cx_mat Kb_g=zeros<cx_mat>(Nn,Nintg);
	cx_mat Kc_g=zeros<cx_mat>(Nn,Nintg);
	cx_mat Kd_g=zeros<cx_mat>(Nn,Nintg);
	
	rowvec ug2(Nug+1);
	ug2.cols(0,Nug-1)=ug.t();
	ug2(Nug)=1.0/(wl-w0l);
	
	mat Wng=wn*ones<rowvec>(Nintg);
	mat Ug=ones<vec>(Nn)*ug2;
	
	dcomplex i(0,1);
	
	mat atang=atan((Wng % (Ug.cols(1,Nintg)-Ug.cols(0,Nintg-1)))/(1+w0l*(Ug.cols(1,Nintg)+Ug.cols(0,Nintg-1))+(pow(w0l,2)+pow(Wng,2)) % Ug.cols(0,Nintg-1) % Ug.cols(1,Nintg)));
	mat logg=log((1.0+2*w0l*Ug.cols(1,Nintg)+pow(Ug.cols(1,Nintg),2) % (pow(Wng,2)+pow(w0l,2)))/(1.0+2*w0l*Ug.cols(0,Nintg-1)+pow(Ug.cols(0,Nintg-1),2) % (pow(Wng,2)+pow(w0l,2))));
	
	Ka_g=-(Ug.cols(1,Nintg)-Ug.cols(0,Nintg-1))/pow(Wng+i*w0l,2)-i*(pow(Ug.cols(1,Nintg),2)-pow(Ug.cols(0,Nintg-1),2))/(2*(Wng+i*w0l))+atang/pow(Wng+i*w0l,3)+i*logg/(2*pow(Wng+i*w0l,3));
	
	Kb_g=-i*(Ug.cols(1,Nintg)-Ug.cols(0,Nintg-1))/(Wng+i*w0l)+i*atang/pow(Wng+i*w0l,2)-logg/(2*pow(Wng+i*w0l,2));
	
	Kc_g=-atang/(Wng+i*w0l)-i*logg/(2*(Wng+i*w0l));
	
	Kd_g=-i*atang+logg/2-log(Ug.cols(1,Nintg)/Ug.cols(0,Nintg-1));
	
	rowvec wg=trans(w.rows(0,Nug));
	//rowvec wg=1/ug2+w0l;
	mat Wg=ones<vec>(Nn)*wg;
	
	if (use_HF_exp)
	{
		double utmp;
		int jn, p;
		for (j=0; j<Nug; j++)
		{
			utmp=ug2(j);
			jn=0;
			while (abs(wn(jn)*utmp)<fg*abs(1+utmp*w0l) && jn<Nn-1) jn++;
			while (pow(wn(jn),pngmax)==0 && jn<Nn-1) jn++;
			if (jn<Nn-1)
			{
				Ka_g.submat(jn,j,Nn-1,j)=-i*(pow(Ug.submat(jn,j+1,Nn-1,j+1),2)-pow(Ug.submat(jn,j,Nn-1,j),2))/(2*(wn.rows(jn,Nn-1)+i*w0l))-(Ug.submat(jn,j+1,Nn-1,j+1)-Ug.submat(jn,j,Nn-1,j))/pow(wn.rows(jn,Nn-1)+i*w0l,2)+log(Ug.submat(jn,j+1,Nn-1,j+1)/Ug.submat(jn,j,Nn-1,j))/pow(i*wn.rows(jn,Nn-1)-w0l,3);
				Kb_g.submat(jn,j,Nn-1,j)=-i*(Ug.submat(jn,j+1,Nn-1,j+1)-Ug.submat(jn,j,Nn-1,j))/(wn.rows(jn,Nn-1)+i*w0l)+log(Ug.submat(jn,j+1,Nn-1,j+1)/Ug.submat(jn,j,Nn-1,j))/pow(i*wn.rows(jn,Nn-1)-w0l,2);
				Kc_g.submat(jn,j,Nn-1,j)=log(Ug.submat(jn,j+1,Nn-1,j+1)/Ug.submat(jn,j,Nn-1,j))/(i*wn.rows(jn,Nn-1)-w0l);
				Kd_g.submat(jn,j,Nn-1,j)=zeros<cx_vec>(Nn-jn);
				for (p=pngmax; p>=1; p--)
				{
					Ka_g.submat(jn,j,Nn-1,j)=Ka_g.submat(jn,j,Nn-1,j)+(pow(i,p)*pow(-1,p+1)*(pow(Wg.submat(jn,j+1,Nn-1,j+1),p)-pow(Wg.submat(jn,j,Nn-1,j),p))/(p*pow(wn.rows(jn,Nn-1),p)))/pow(i*wn.rows(jn,Nn-1)-w0l,3);
					Kb_g.submat(jn,j,Nn-1,j)=Kb_g.submat(jn,j,Nn-1,j)+(pow(i,p)*pow(-1,p+1)*(pow(Wg.submat(jn,j+1,Nn-1,j+1),p)-pow(Wg.submat(jn,j,Nn-1,j),p))/(p*pow(wn.rows(jn,Nn-1),p)))/pow(i*wn.rows(jn,Nn-1)-w0l,2);
					Kc_g.submat(jn,j,Nn-1,j)=Kc_g.submat(jn,j,Nn-1,j)+(pow(i,p)*pow(-1,p+1)*(pow(Wg.submat(jn,j+1,Nn-1,j+1),p)-pow(Wg.submat(jn,j,Nn-1,j),p))/(p*pow(wn.rows(jn,Nn-1),p)))/(i*wn.rows(jn,Nn-1)-w0l);
					Kd_g.submat(jn,j,Nn-1,j)=Kd_g.submat(jn,j,Nn-1,j)+pow(i,p)*pow(-1,p+1)*(pow(Wg.submat(jn,j+1,Nn-1,j+1),p)-pow(Wg.submat(jn,j,Nn-1,j),p))/(p*pow(wn.rows(jn,Nn-1),p));
				}
			}
		}
	}
	
	Ka_g=-Ka_g/(2*PI);
	Kb_g=-Kb_g/(2*PI);
	Kc_g=-Kc_g/(2*PI);
	Kd_g=-Kd_g/(2*PI);
	
	
	cx_mat Ka_c=zeros<cx_mat>(Nn,Nintc);
	cx_mat Kb_c=zeros<cx_mat>(Nn,Nintc);
	cx_mat Kc_c=zeros<cx_mat>(Nn,Nintc);
	cx_mat Kd_c=zeros<cx_mat>(Nn,Nintc);
	
	mat Wnc=wn*ones<rowvec>(Nintc);
	mat Wc=ones<vec>(Nn)*wc.t();
	
	mat logc=log((pow(Wnc,2)+pow(Wc.cols(1,Nintc),2))/(pow(Wnc,2)+pow(Wc.cols(0,Nintc-1),2)));
	mat atanc2=atan((Wnc % (Wc.cols(1,Nintc)-Wc.cols(0,Nintc-1)))/(Wc.cols(1,Nintc) % Wc.cols(0,Nintc-1)+pow(Wnc,2)));
	cx_mat logc2=logc/2+i*atanc2;
	cx_mat dWn=i*Wnc-Wc.cols(0,Nintc-1);
	mat dWc=Wc.cols(1,Nintc)-Wc.cols(0,Nintc-1);
	
	Ka_c=-pow(dWn,2) % dWc-dWn % pow(dWc,2)/2-pow(dWc,3)/3-pow(dWn,3) % logc2;
	Kb_c=-dWn % dWc-pow(dWc,2)/2-pow(dWn,2) % logc2;
	Kc_c=-dWc-dWn % logc2;
	Kd_c=-logc2;
	
	int Pmax=2*pnmax+4;
	imat MP;
	pascal(Pmax+1,MP);
	
	if (use_HF_exp)
	{
		double wtmp;
		int jni, jnr, p, l;
		for (j=0; j<Nwc-1; j++)
		{
			wtmp=abs(wc(j));
			if (abs(wc(j+1))>wtmp) wtmp=abs(wc(j+1));
			jni=0;
			while (wn(jni)<fi*wtmp && jni<Nn-1) jni++;
			while (pow(wn(jni),2*pnmax+1)==0 && jni<Nn-1) jni++;
			jnr=1;
			while (wn(jnr)<fr*wtmp && jnr<Nn-1) jnr++;
			while (pow(wn(jnr),2*pnmax)==0 && jnr<Nn-1) jnr++;
			if (jni<Nn-1 || jnr<Nn-1)
			{
				double wj=wc(j);
				double dw1=wc(j+1)-wj;
				vec dwp=zeros<vec>(Pmax);
				dwp(0)=dw1;
				dwp(1)=pow(dw1,2)+2*wj*dw1;
				for (p=3; p<=Pmax; p++)
				{
					dwp(p-1)=0;
					for (l=0; l<p; l++)
					{
						dwp(p-1)=dwp(p-1)+MP(l,p-l)*pow(dw1,p-l)*pow(wj,l);
					}
				}
				cx_vec vtmp;
				if (jni<Nn)
				{
					vtmp.zeros(Nn-jni);
					vtmp.set_real(real(Ka_c.submat(jni,j,Nn-1,j)));
					Ka_c.submat(jni,j,Nn-1,j)=vtmp;
					vtmp.set_real(real(Kb_c.submat(jni,j,Nn-1,j)));
					Kb_c.submat(jni,j,Nn-1,j)=vtmp;
					vtmp.set_real(real(Kc_c.submat(jni,j,Nn-1,j)));
					Kc_c.submat(jni,j,Nn-1,j)=vtmp;
					vtmp.set_real(real(Kd_c.submat(jni,j,Nn-1,j)));
					Kd_c.submat(jni,j,Nn-1,j)=vtmp;
					
					for (p=2*pnmax+1; p>=1; p-=2)
					{
						Ka_c.submat(jni,j,Nn-1,j)=Ka_c.submat(jni,j,Nn-1,j) + i*pow(-1,(p-1)/2)*(pow(wj,3)*dwp(p-1)/p - 3.0*pow(wj,2)*dwp(p)/(p+1) + 3.0*wj*dwp(p+1)/(p+2) - dwp(p+2)/(p+3))/pow(wn.rows(jni,Nn-1),p);
						Kb_c.submat(jni,j,Nn-1,j)=Kb_c.submat(jni,j,Nn-1,j) + i*pow(-1,(p-1)/2)*(-pow(wj,2)*dwp(p-1)/p + 2*wj*dwp(p)/(p+1) - dwp(p+1)/(p+2))/pow(wn.rows(jni,Nn-1),p);
						Kc_c.submat(jni,j,Nn-1,j)=Kc_c.submat(jni,j,Nn-1,j) + i*pow(-1,(p-1)/2)*(wj*dwp(p-1)/p - dwp(p)/(p+1))/pow(wn.rows(jni,Nn-1),p);
						Kd_c.submat(jni,j,Nn-1,j)=Kd_c.submat(jni,j,Nn-1,j) - i*pow(-1,(p-1)/2)*dwp(p-1)/(p*pow(wn.rows(jni,Nn-1),p));
					}
				}
				if (jnr<Nn)
				{
					vtmp.zeros(Nn-jnr);
					vtmp.set_imag(imag(Ka_c.submat(jnr,j,Nn-1,j)));
					Ka_c.submat(jnr,j,Nn-1,j)=vtmp;
					vtmp.set_imag(imag(Kb_c.submat(jnr,j,Nn-1,j)));
					Kb_c.submat(jnr,j,Nn-1,j)=vtmp;
					vtmp.set_imag(imag(Kc_c.submat(jnr,j,Nn-1,j)));
					Kc_c.submat(jnr,j,Nn-1,j)=vtmp;
					vtmp.set_imag(imag(Kd_c.submat(jnr,j,Nn-1,j)));
					Kd_c.submat(jnr,j,Nn-1,j)=vtmp;
					for (p=2*pnmax; p>=2; p-=2)
					{
						Ka_c.submat(jnr,j,Nn-1,j)=Ka_c.submat(jnr,j,Nn-1,j) + pow(-1,(p-2)/2)*(pow(wj,3)*dwp(p-1)/p - 3*pow(wj,2)*dwp(p)/(p+1) + 3*wj*dwp(p+1)/(p+2) - dwp(p+2)/(p+3))/pow(wn.rows(jnr,Nn-1),p);
						Kb_c.submat(jnr,j,Nn-1,j)=Kb_c.submat(jnr,j,Nn-1,j) + pow(-1,(p-2)/2)*(-pow(wj,2)*dwp(p-1)/p + 2*wj*dwp(p)/(p+1) - dwp(p+1)/(p+2))/pow(wn.rows(jnr,Nn-1),p);
						Kc_c.submat(jnr,j,Nn-1,j)=Kc_c.submat(jnr,j,Nn-1,j) + pow(-1,(p-2)/2)*(wj*dwp(p-1)/p - dwp(p)/(p+1))/pow(wn.rows(jnr,Nn-1),p);
						Kd_c.submat(jnr,j,Nn-1,j)=Kd_c.submat(jnr,j,Nn-1,j) - pow(-1,(p-2)/2)*dwp(p-1)/(p*pow(wn.rows(jnr,Nn-1),p));
					}
				}
			}
		}
	}
	
	Ka_c=Ka_c/(2*PI);
	Kb_c=Kb_c/(2*PI);
	Kc_c=Kc_c/(2*PI);
	Kd_c=Kd_c/(2*PI);
	
	cx_mat Ka_d=zeros<cx_mat>(Nn,Nintd);
	cx_mat Kb_d=zeros<cx_mat>(Nn,Nintd);
	cx_mat Kc_d=zeros<cx_mat>(Nn,Nintd);
	cx_mat Kd_d=zeros<cx_mat>(Nn,Nintd);
	
	rowvec ud2(Nud+1);
	ud2.cols(1,Nud)=ud.t();
	ud2(0)=1.0/(wr-w0r);
	
	mat Wnd=wn*ones<rowvec>(Nintd);
	mat Ud=ones<vec>(Nn)*ud2;
	
	mat atand=atan((Wnd % (Ud.cols(1,Nintd)-Ud.cols(0,Nintd-1)))/(1.0+w0r*(Ud.cols(1,Nintd)+Ud.cols(0,Nintd-1))+(pow(w0r,2)+pow(Wnd,2)) % Ud.cols(0,Nintd-1) % Ud.cols(1,Nintd)));
	mat logd=log((1+2*w0r*Ud.cols(1,Nintd)+pow(Ud.cols(1,Nintd),2) % (pow(Wnd,2)+pow(w0r,2)))/(1+2*w0r*Ud.cols(0,Nintd-1)+pow(Ud.cols(0,Nintd-1),2) % (pow(Wnd,2)+pow(w0r,2))));
	
	Ka_d=-(Ud.cols(1,Nintd)-Ud.cols(0,Nintd-1))/pow(Wnd+i*w0r,2) - i*(pow(Ud.cols(1,Nintd),2)-pow(Ud.cols(0,Nintd-1),2))/(2*(Wnd+i*w0r))+atand/pow(Wnd+i*w0r,3) +i*logd/(2*pow(Wnd+i*w0r,3));
	
	Kb_d=-i*(Ud.cols(1,Nintd)-Ud.cols(0,Nintd-1))/(Wnd+i*w0r)+i*atand/pow(Wnd+i*w0r,2)-logd/(2*pow(Wnd+i*w0r,2));
	
	Kc_d=-atand/(Wnd+i*w0r)-i*logd/(2*(Wnd+i*w0r));
	
	Kd_d=-i*atand-log(Ud.cols(1,Nintd)/Ud.cols(0,Nintd-1))+logd/2;
	
	//rowvec wd=1/ud2+w0r;
	rowvec wd=trans(w.rows(Nw_lims(1),Nw-1));
	mat Wd=ones<vec>(Nn)*wd;
	
	if (use_HF_exp)
	{
		double utmp;
		int jn, p;
		for (j=0; j<Nud; j++)
		{
			utmp=ud2(j);
			jn=0;
			while (abs(wn(jn)*utmp)<fd*abs(1+utmp*w0r) && jn<Nn-1) jn++;
			while (pow(wn(jn),pndmax) && jn<Nn-1) jn++;
			if (jn<Nn-1)
			{
				Ka_d.submat(jn,j,Nn-1,j)=-(Ud.submat(jn,j+1,Nn-1,j+1)-Ud.submat(jn,j,Nn-1,j))/pow(wn.rows(jn,Nn-1)+i*w0r,2) - i*(pow(Ud.submat(jn,j+1,Nn-1,j+1),2)-pow(Ud.submat(jn,j,Nn-1,j),2))/(2*(wn.rows(jn,Nn-1)+i*w0r))+log(Ud.submat(jn,j+1,Nn-1,j+1)/Ud.submat(jn,j,Nn-1,j))/pow(i*wn.rows(jn,Nn-1)-w0r,3);
				Kb_d.submat(jn,j,Nn-1,j)=-i*(Ud.submat(jn,j+1,Nn-1,j+1)-Ud.submat(jn,j,Nn-1,j))/(wn.rows(jn,Nn-1)+i*w0r)+log(Ud.submat(jn,j+1,Nn-1,j+1)/Ud.submat(jn,j,Nn-1,j))/pow(i*wn.rows(jn,Nn-1)-w0r,2);
				Kc_d.submat(jn,j,Nn-1,j)=log(Ud.submat(jn,j+1,Nn-1,j+1)/Ud.submat(jn,j,Nn-1,j))/(i*wn.rows(jn,Nn-1)-w0r);
				Kd_d.submat(jn,j,Nn-1,j)=zeros<cx_vec>(Nn-jn);
				for (p=pndmax; p>=1; p--)
				{
					Ka_d.submat(jn,j,Nn-1,j)=Ka_d.submat(jn,j,Nn-1,j)+(pow(i,p)*pow(-1,p+1)*(pow(Wd.submat(jn,j+1,Nn-1,j+1),p)-pow(Wd.submat(jn,j,Nn-1,j),p))/(p*pow(wn.rows(jn,Nn-1),p)))/pow(i*wn.rows(jn,Nn-1)-w0r,3);
					Kb_d.submat(jn,j,Nn-1,j)=Kb_d.submat(jn,j,Nn-1,j)+(pow(i,p)*pow(-1,p+1)*(pow(Wd.submat(jn,j+1,Nn-1,j+1),p)-pow(Wd.submat(jn,j,Nn-1,j),p))/(p*pow(wn.rows(jn,Nn-1),p)))/pow(i*wn.rows(jn,Nn-1)-w0r,2);
					Kc_d.submat(jn,j,Nn-1,j)=Kc_d.submat(jn,j,Nn-1,j)+(pow(i,p)*pow(-1,p+1)*(pow(Wd.submat(jn,j+1,Nn-1,j+1),p)-pow(Wd.submat(jn,j,Nn-1,j),p))/(p*pow(wn.rows(jn,Nn-1),p)))/(i*wn.rows(jn,Nn-1)-w0r);
					Kd_d.submat(jn,j,Nn-1,j)=Kd_d.submat(jn,j,Nn-1,j)+pow(i,p)*pow(-1,p+1)*(pow(Wd.submat(jn,j+1,Nn-1,j+1),p)-pow(Wd.submat(jn,j,Nn-1,j),p))/(p*pow(wn.rows(jn,Nn-1),p));
				}
			}
		}
	}
	
	Ka_d=-Ka_d/(2*PI);
	Kb_d=-Kb_d/(2*PI);
	Kc_d=-Kc_d/(2*PI);
	Kd_d=-Kd_d/(2*PI);
	
	cx_mat KG=(Ka_g*Pa_g+Kb_g*Pb_g+Kc_g*Pc_g+Kd_g*Pd_g)*MM;
	cx_mat KC=(Ka_c*Pa_c+Kb_c*Pb_c+Kc_c*Pc_c+Kd_c*Pd_c)*MM;
	cx_mat KD=(Ka_d*Pa_d+Kb_d*Pb_d+Kc_d*Pc_d+Kd_d*Pd_d)*MM;
	
	Kcx=KG+KC+KD;
	
	//////////////////////////////////////
	
	
	rowvec Knorm_a_g=zeros<rowvec>(Nintg);
	rowvec Knorm_b_g=zeros<rowvec>(Nintg);
	rowvec Knorm_c_g=zeros<rowvec>(Nintg);
	rowvec Knorm_d_g=zeros<rowvec>(Nintg);
	
	Knorm_a_g=-(pow(ug2.cols(1,Nug),2)-pow(ug2.cols(0,Nug-1),2))/2;
	Knorm_b_g=-(ug2.cols(1,Nug)-ug2.cols(0,Nug-1));
	Knorm_c_g=-log(ug2.cols(1,Nug)/ug2.cols(0,Nug-1));
	Knorm_d_g=1.0/ug2.cols(1,Nug)-1.0/ug2.cols(0,Nug-1);
	
	rowvec KM0g=(Knorm_a_g*Pa_g+Knorm_b_g*Pb_g+Knorm_c_g*Pc_g+Knorm_d_g*Pd_g)*MM/(2*PI);
	
	rowvec Knorm_a_c=zeros<rowvec>(Nintc);
	rowvec Knorm_b_c=zeros<rowvec>(Nintc);
	rowvec Knorm_c_c=zeros<rowvec>(Nintc);
	rowvec Knorm_d_c=zeros<rowvec>(Nintc);
	
	Knorm_a_c=trans(pow(wc.rows(1,Nwc-1)-wc.rows(0,Nwc-2),4)/4);
	Knorm_b_c=trans(pow(wc.rows(1,Nwc-1)-wc.rows(0,Nwc-2),3)/3);
	Knorm_c_c=trans(pow(wc.rows(1,Nwc-1)-wc.rows(0,Nwc-2),2)/2);
	Knorm_d_c=trans(wc.rows(1,Nwc-1)-wc.rows(0,Nwc-2));
	
	rowvec KM0c=(Knorm_a_c*Pa_c+Knorm_b_c*Pb_c+Knorm_c_c*Pc_c+Knorm_d_c*Pd_c)*MM/(2*PI);
	
	rowvec Knorm_a_d=zeros<rowvec>(Nintd);
	rowvec Knorm_b_d=zeros<rowvec>(Nintd);
	rowvec Knorm_c_d=zeros<rowvec>(Nintd);
	rowvec Knorm_d_d=zeros<rowvec>(Nintd);
	
	Knorm_a_d=-(pow(ud2.cols(1,Nud),2)-pow(ud2.cols(0,Nud-1),2))/2;
	Knorm_b_d=-(ud2.cols(1,Nud)-ud2.cols(0,Nud-1));
	Knorm_c_d=-log(ud2.cols(1,Nud)/ud2.cols(0,Nud-1));
	Knorm_d_d=1.0/ud2.cols(1,Nud)-1.0/ud2.cols(0,Nud-1);
	
	rowvec KM0d=(Knorm_a_d*Pa_d+Knorm_b_d*Pb_d+Knorm_c_d*Pc_d+Knorm_d_d*Pd_d)*MM/(2*PI);
	
	rowvec KM0=KM0g+KM0c+KM0d;
	
	
	rowvec KM1_a_g=zeros<rowvec>(Nintg);
	rowvec KM1_b_g=zeros<rowvec>(Nintg);
	rowvec KM1_c_g=zeros<rowvec>(Nintg);
	rowvec KM1_d_g=zeros<rowvec>(Nintg);
	
	KM1_a_g=Knorm_b_g;
	KM1_b_g=Knorm_c_g;
	KM1_c_g=Knorm_d_g;
	KM1_d_g=(1.0/pow(ug2.cols(1,Nug),2)-1.0/pow(ug2.cols(0,Nug-1),2))/2;
	
	rowvec KM1g_tmp=(KM1_a_g*Pa_g+KM1_b_g*Pb_g+KM1_c_g*Pc_g+KM1_d_g*Pd_g)*MM/(2*PI);
	
	rowvec KM1g=w0l*KM0g+KM1g_tmp;
	
	mat Wjc=diagmat(wc.rows(0,Nintc-1));
	
	rowvec KM1_a_c_tmp=trans(pow(wc.rows(1,Nwc-1)-wc.rows(0,Nwc-2),5)/5.0);
	
	rowvec KM1_a_c=KM1_a_c_tmp+Knorm_a_c*Wjc;
	rowvec KM1_b_c=Knorm_a_c+Knorm_b_c*Wjc;
	rowvec KM1_c_c=Knorm_b_c+Knorm_c_c*Wjc;
	rowvec KM1_d_c=Knorm_c_c+Knorm_d_c*Wjc;
	
	rowvec KM1c=(KM1_a_c*Pa_c+KM1_b_c*Pb_c+KM1_c_c*Pc_c+KM1_d_c*Pd_c)*MM/(2*PI);
	
	rowvec KM1_a_d=zeros<rowvec>(Nintd);
	rowvec KM1_b_d=zeros<rowvec>(Nintd);
	rowvec KM1_c_d=zeros<rowvec>(Nintd);
	rowvec KM1_d_d=zeros<rowvec>(Nintd);
	
	KM1_a_d=Knorm_b_d;
	KM1_b_d=Knorm_c_d;
	KM1_c_d=Knorm_d_d;
	KM1_d_d=(1.0/pow(ud2.cols(1,Nud),2)-1.0/pow(ud2.cols(0,Nud-1),2))/2;
	
	rowvec KM1d_tmp=(KM1_a_d*Pa_d+KM1_b_d*Pb_d+KM1_c_d*Pc_d+KM1_d_d*Pd_d)*MM/(2*PI);
	
	rowvec KM1d=w0r*KM0d+KM1d_tmp;
	
	rowvec KM1=KM1g+KM1c+KM1d;
	
	rowvec KM2_a_g=zeros<rowvec>(Nintg);
	rowvec KM2_b_g=zeros<rowvec>(Nintg);
	rowvec KM2_c_g=zeros<rowvec>(Nintg);
	rowvec KM2_d_g=zeros<rowvec>(Nintg);
	
	KM2_a_g=Knorm_c_g;
	KM2_b_g=Knorm_d_g;
	KM2_c_g=KM1_d_g;
	KM2_d_g=(1.0/pow(ug2.cols(1,Nug),3)-1.0/pow(ug2.cols(0,Nug-1),3))/3;
	
	rowvec KM2g_tmp=(KM2_a_g*Pa_g+KM2_b_g*Pb_g+KM2_c_g*Pc_g+KM2_d_g*Pd_g)*MM/(2*PI);
	
	rowvec KM2g=pow(w0l,2)*KM0g+2*w0l*KM1g_tmp+KM2g_tmp;
	
	rowvec KM2_a_c_tmp=trans(pow(wc.rows(1,Nwc-1)-wc.rows(0,Nwc-2),6)/6);
	
	rowvec KM2_a_c=KM2_a_c_tmp+2*KM1_a_c_tmp*Wjc+Knorm_a_c*pow(Wjc,2);
	rowvec KM2_b_c=KM1_a_c_tmp+2*Knorm_a_c*Wjc+Knorm_b_c*pow(Wjc,2);
	rowvec KM2_c_c=Knorm_a_c+2*Knorm_b_c*Wjc+Knorm_c_c*pow(Wjc,2);
	rowvec KM2_d_c=Knorm_b_c+2*Knorm_c_c*Wjc+Knorm_d_c*pow(Wjc,2);
	
	rowvec KM2c=(KM2_a_c*Pa_c+KM2_b_c*Pb_c+KM2_c_c*Pc_c+KM2_d_c*Pd_c)*MM/(2*PI);
	
	rowvec KM2_a_d=zeros<rowvec>(Nintd);
	rowvec KM2_b_d=zeros<rowvec>(Nintd);
	rowvec KM2_c_d=zeros<rowvec>(Nintd);
	rowvec KM2_d_d=zeros<rowvec>(Nintd);
	
	KM2_a_d=Knorm_c_d;
	KM2_b_d=Knorm_d_d;
	KM2_c_d=KM1_d_d;
	KM2_d_d=(1.0/pow(ud2.cols(1,Nud),3)-1.0/pow(ud2.cols(0,Nud-1),3))/3;
	
	rowvec KM2d_tmp=(KM2_a_d*Pa_d+KM2_b_d*Pb_d+KM2_c_d*Pc_d+KM2_d_d*Pd_d)*MM/(2*PI);
	
	rowvec KM2d=pow(w0r,2)*KM0d+2*w0r*KM1d_tmp+KM2d_tmp;
	
	rowvec KM2=KM2g+KM2c+KM2d;
	
	rowvec KM3_a_g=zeros<rowvec>(Nintg);
	rowvec KM3_b_g=zeros<rowvec>(Nintg);
	rowvec KM3_c_g=zeros<rowvec>(Nintg);
	rowvec KM3_d_g=zeros<rowvec>(Nintg);
	
	KM3_a_g=Knorm_d_g;
	KM3_b_g=KM1_d_g;
	KM3_c_g=KM2_d_g;
	KM3_d_g=(1.0/pow(ug2.cols(1,Nug),4)-1.0/pow(ug2.cols(0,Nug-1),4))/4;
	
	rowvec KM3g_tmp=(KM3_a_g*Pa_g+KM3_b_g*Pb_g+KM3_c_g*Pc_g+KM3_d_g*Pd_g)*MM/(2*PI);
	
	rowvec KM3g=pow(w0l,3)*KM0g+3*pow(w0l,2)*KM1g_tmp+3*w0l*KM2g_tmp+KM3g_tmp;
	
	rowvec KM3_a_c_tmp=trans(pow(wc.rows(1,Nwc-1)-wc.rows(0,Nwc-2),7)/7);
	
	rowvec KM3_a_c=KM3_a_c_tmp+3*KM2_a_c_tmp*Wjc+3*KM1_a_c_tmp*pow(Wjc,2)+Knorm_a_c*pow(Wjc,3);
	rowvec KM3_b_c=KM2_a_c_tmp+3*KM1_a_c_tmp*Wjc+3*Knorm_a_c*pow(Wjc,2)+Knorm_b_c*pow(Wjc,3);
	rowvec KM3_c_c=KM1_a_c_tmp+3*Knorm_a_c*Wjc+3*Knorm_b_c*pow(Wjc,2)+Knorm_c_c*pow(Wjc,3);
	rowvec KM3_d_c=Knorm_a_c+3*Knorm_b_c*Wjc+3*Knorm_c_c*pow(Wjc,2)+Knorm_d_c*pow(Wjc,3);
	
	rowvec KM3c=(KM3_a_c*Pa_c+KM3_b_c*Pb_c+KM3_c_c*Pc_c+KM3_d_c*Pd_c)*MM/(2*PI);
	
	rowvec KM3_a_d=zeros<rowvec>(Nintd);
	rowvec KM3_b_d=zeros<rowvec>(Nintd);
	rowvec KM3_c_d=zeros<rowvec>(Nintd);
	rowvec KM3_d_d=zeros<rowvec>(Nintd);
	
	KM3_a_d=Knorm_d_d;
	KM3_b_d=KM1_d_d;
	KM3_c_d=KM2_d_d;
	KM3_d_d=(1.0/pow(ud2.cols(1,Nud),4)-1.0/pow(ud2.cols(0,Nud-1),4))/4;
	
	rowvec KM3d_tmp=(KM3_a_d*Pa_d+KM3_b_d*Pb_d+KM3_c_d*Pc_d+KM3_d_d*Pd_d)*MM/(2*PI);
	
	rowvec KM3d=pow(w0r,3)*KM0d+3*pow(w0r,2)*KM1d_tmp+3*w0r*KM2d_tmp+KM3d_tmp;
	
	rowvec KM3=KM3g+KM3c+KM3d;
	
	KM.zeros(4,Nw);
	
	KM.row(0)=KM0;
	KM.row(1)=KM1;
	KM.row(2)=KM2;
	KM.row(3)=KM3;
	
	K.zeros(2*Nn,Nw);
	uvec even_ind=linspace<uvec>(0,2*Nn-2,Nn);
	K.rows(even_ind)=real(Kcx);
	K.rows(even_ind+1)=imag(Kcx);
	
	/*
	 //check moments
	 rowvec x0={-2, 0, 3};
	 rowvec s0={0.8, 0.3, 1};
	 rowvec wgt={1,1,1};
	 
	 vec test_A;
	 sum_gaussians(w, x0, s0, wgt, test_A);
	 
 	cout<<KM*test_A<<endl;
	*/
	
	return true;
}







bool OmegaMaxEnt_data::Kernel_G_fermions()
{
	bool use_HF_exp=true;
	double fg=1.7;
	int pngmax=100;
	double fi=fg;
	double fr=fi;
	int pnmax=50;
	double fd=fg;
	int pndmax=pngmax;
	
	cout<<"defining kernel matrix...\n";
	
	Kcx.zeros(Nn,Nw);
 	KM.zeros(4,Nw);
 
	int Nug=Nw_lims(0);
 	int Nud=Nw-Nw_lims(1)-1;
 
	vec ug, ud;
	if (Du_constant)
	{
 		double dug=dwl/((wl-dwl-w0l)*(wl-w0l));
		vec ug_int=linspace<vec>(1,Nug,Nug);
		ug=-dug*ug_int;
		
 		double dud=dwr/((wr-w0r)*(wr+dwr-w0r));
		vec ud_int=linspace<vec>(Nud,1,Nud);
		ud=dud*ud_int;
	}
	else
	{
		ug=1.0/(w.rows(0,Nug-1)-w0l);
		ud=1.0/(w.rows(Nw_lims(1)+1,Nw-1)-w0r);
	}
 
	mat MM;
	spline_matrix_G_part(w, Nw_lims, ws, MM);
	
	int Ncfs0=MM.n_rows;


	mat MC(Ncfs0+Nw,Nw);
	
	MC.submat(0,0,Ncfs0-1,Nw-1)=MM;
	MC.submat(Ncfs0,0,Ncfs0+Nw-1,Nw-1)=eye<mat>(Nw,Nw);
	 
	int Nint=Nw+1;
	int Nintg=Nug+1;

	mat Pa_g=zeros<mat>(Nintg,3*Nint-2+Nw);
	mat Pb_g_r=zeros<mat>(Nintg,3*Nint-2+Nw);
	mat Pc_g_r=zeros<mat>(Nintg,3*Nint-2+Nw);
	mat Pd_g_r=zeros<mat>(Nintg,3*Nint-2+Nw);

	Pa_g(0,0)=1;
	Pb_g_r(0,1)=1;
	
	int j;
	for (j=2; j<=Nintg; j++)
	{
		Pa_g(j-1,3*j-4)=1;
		Pb_g_r(j-1,3*j-3)=1;
		Pc_g_r(j-1,3*j-2)=1;
		Pd_g_r(j-1,j+3*Nint-4)=1;
	}
	
	mat U=zeros<mat>(Nug+1,Nug+1);
	vec ug1=zeros<vec>(Nug+1);
	ug1.rows(1,Nug)=ug;
	
	U.diag()=ug1;
		
	mat Pb_g=Pb_g_r-3*U*Pa_g;
	mat Pc_g=Pc_g_r+3*pow(U,2)*Pa_g-2*U*Pb_g_r;
	mat Pd_g=Pd_g_r-pow(U,3)*Pa_g+pow(U,2)*Pb_g_r-U*Pc_g_r;
	
	int	Nintc=Nwc-1;
		
	mat	Pa_c=zeros<mat>(Nintc,3*Nint-2+Nw);
	mat	Pb_c_r=zeros<mat>(Nintc,3*Nint-2+Nw);
	mat	Pc_c_r=zeros<mat>(Nintc,3*Nint-2+Nw);
	mat	Pd_c_r=zeros<mat>(Nintc,3*Nint-2+Nw);
		
	for (j=1; j<=Nintc; j++)
	{
	 	Pa_c(j-1,3*j-4+3*Nintg)=1;
		Pb_c_r(j-1,3*j-3+3*Nintg)=1;
		Pc_c_r(j-1,3*j-2+3*Nintg)=1;
		Pd_c_r(j-1,j+3*Nint-3+Nug)=1;
	}
	
	int Nintd=Nud+1;
	
	mat Pa_d=zeros<mat>(Nintd,3*Nint-2+Nw);
	mat Pb_d_r=zeros<mat>(Nintd,3*Nint-2+Nw);
	mat Pc_d_r=zeros<mat>(Nintd,3*Nint-2+Nw);
	mat Pd_d_r=zeros<mat>(Nintd,3*Nint-2+Nw);
	
	for (j=1; j<Nintd; j++)
	{
		Pa_d(j-1,3*j-4+3*Nintg+3*Nintc)=1;
		Pb_d_r(j-1,3*j-3+3*Nintg+3*Nintc)=1;
		Pc_d_r(j-1,3*j-2+3*Nintg+3*Nintc)=1;
		Pd_d_r(j-1,j+3*Nint-3+Nug+Nwc)=1;
	}
	
	j=Nintd;
	Pa_d(j-1,3*j-4+3*Nintg+3*Nintc)=1;
	Pb_d_r(j-1,3*j-3+3*Nintg+3*Nintc)=1;
	
	U.zeros(Nud+1,Nud+1);
	vec ud1=zeros<vec>(Nud+1);
	ud1.rows(0,Nud-1)=ud;
	U.diag()=ud1;
	
	mat Pb_d=Pb_d_r-3*U*Pa_d;
	mat Pc_d=Pc_d_r+3*pow(U,2)*Pa_d-2*U*Pb_d_r;
	mat Pd_d=Pd_d_r-pow(U,3)*Pa_d+pow(U,2)*Pb_d_r-U*Pc_d_r;

// check projectors
//	vec a=Pa_g*MC*test_A;
//	vec b=Pb_g_r*MC*test_A;
//	vec c=Pc_g_r*MC*test_A;
//	vec d=Pd_g_r*MC*test_A;
//	ind_tmp=linspace<uvec>(0,Nug,Nug+1);
//	vec a=Pa_c*MC*test_A;
//	vec b=Pb_c_r*MC*test_A;
//	vec c=Pc_c_r*MC*test_A;
//	vec d=Pd_c_r*MC*test_A;
//	ind_tmp=linspace<uvec>(Nug+1,Nug+Nintc,Nintc);
//	vec a=Pa_d*MC*test_A;
//	vec b=Pb_d_r*MC*test_A;
//	vec c=Pc_d_r*MC*test_A;
//	vec d=Pd_d_r*MC*test_A;
//	ind_tmp=linspace<uvec>(Nug+Nintc+1,Nug+Nintc+Nud+1,Nud+1);
	
//	cout<<"a:\n"<<a<<endl;
//	cout<<"diff(a):\n"<<a-coeffs1.rows(4*ind_tmp)<<endl;
//	cout<<"b:\n"<<b<<endl;
//	cout<<"diff(b):\n"<<b-coeffs1.rows(4*ind_tmp+1)<<endl;
//	cout<<"c:\n"<<c<<endl;
//	cout<<"diff(c):\n"<<c-coeffs1.rows(4*ind_tmp+2)<<endl;
//	cout<<"d:\n"<<d<<endl;
//	cout<<"diff(d):\n"<<d-coeffs1.rows(4*ind_tmp+3)<<endl;
	

	cx_mat Ka_g=zeros<cx_mat>(Nn,Nintg);
	cx_mat Kb_g=zeros<cx_mat>(Nn,Nintg);
	cx_mat Kc_g=zeros<cx_mat>(Nn,Nintg);
	cx_mat Kd_g=zeros<cx_mat>(Nn,Nintg);
	
	rowvec ug2(Nug+1);
	ug2.cols(0,Nug-1)=ug.t();
	ug2(Nug)=1.0/(wl-w0l);
	
	mat Wng=wn*ones<rowvec>(Nintg-1);
	mat Ug=ones<vec>(Nn)*ug2;
	
	dcomplex i(0,1);
	
	Ka_g.col(0)=-ug(0)/pow(wn+i*w0l,2)-i*pow(ug(0),2)/(2*(wn+i*w0l))+atan((wn*ug(0))/(1.0+ug(0)*w0l))/pow(wn+i*w0l,3)+i*log(1.0+2*w0l*ug(0)+pow(ug(0),2)*(pow(wn,2)+pow(w0l,2)))/(2*pow(wn+i*w0l,3));
	Kb_g.col(0)=-i*ug(0)/(wn+i*w0l)+i*atan((wn*ug(0))/(1+ug(0)*w0l))/pow(wn+i*w0l,2)-log(1.0+2*w0l*ug(0)+pow(ug(0),2)*(pow(wn,2)+pow(w0l,2)))/(2*pow(wn+i*w0l,2));

	mat atang=atan((Wng % (Ug.cols(1,Nintg-1)-Ug.cols(0,Nintg-2)))/(1+w0l*(Ug.cols(1,Nintg-1)+Ug.cols(0,Nintg-2))+(pow(w0l,2)+pow(Wng,2)) % Ug.cols(0,Nintg-2) % Ug.cols(1,Nintg-1)));
	mat logg=log((1.0+2*w0l*Ug.cols(1,Nintg-1)+pow(Ug.cols(1,Nintg-1),2) % (pow(Wng,2)+pow(w0l,2)))/(1.0+2*w0l*Ug.cols(0,Nintg-2)+pow(Ug.cols(0,Nintg-2),2) % (pow(Wng,2)+pow(w0l,2))));
	 
	Ka_g.cols(1,Nintg-1)=-(Ug.cols(1,Nintg-1)-Ug.cols(0,Nintg-2))/pow(Wng+i*w0l,2)-i*(pow(Ug.cols(1,Nintg-1),2)-pow(Ug.cols(0,Nintg-2),2))/(2*(Wng+i*w0l))+atang/pow(Wng+i*w0l,3)+i*logg/(2*pow(Wng+i*w0l,3));
	 
	Kb_g.cols(1,Nintg-1)=-i*(Ug.cols(1,Nintg-1)-Ug.cols(0,Nintg-2))/(Wng+i*w0l)+i*atang/pow(Wng+i*w0l,2)-logg/(2*pow(Wng+i*w0l,2));
	 
	Kc_g.cols(1,Nintg-1)=-atang/(Wng+i*w0l)-i*logg/(2*(Wng+i*w0l));
	 
	Kd_g.cols(1,Nintg-1)=-i*atang+logg/2-log(Ug.cols(1,Nintg-1)/Ug.cols(0,Nintg-2));

	
	rowvec wg=1/ug2+w0l;
	mat Wg=ones<vec>(Nn)*wg;
	 
	if (use_HF_exp)
	{
		double utmp;
		int jn, p;
		for (j=1; j<=Nug; j++)
		{
			utmp=ug2(j-1);
			jn=0;
			while (abs(wn(jn)*utmp)<fg*abs(1+utmp*w0l) && jn<Nn-1) jn++;
			while (pow(wn(jn),pngmax)==0 && jn<Nn-1) jn++;
			if (jn<Nn-1)
			{
				Ka_g.submat(jn,j,Nn-1,j)=-i*(pow(Ug.submat(jn,j,Nn-1,j),2)-pow(Ug.submat(jn,j-1,Nn-1,j-1),2))/(2*(wn.rows(jn,Nn-1)+i*w0l))-(Ug.submat(jn,j,Nn-1,j)-Ug.submat(jn,j-1,Nn-1,j-1))/pow(wn.rows(jn,Nn-1)+i*w0l,2)+log(Ug.submat(jn,j,Nn-1,j)/Ug.submat(jn,j-1,Nn-1,j-1))/pow(i*wn.rows(jn,Nn-1)-w0l,3);
				Kb_g.submat(jn,j,Nn-1,j)=-i*(Ug.submat(jn,j,Nn-1,j)-Ug.submat(jn,j-1,Nn-1,j-1))/(wn.rows(jn,Nn-1)+i*w0l)+log(Ug.submat(jn,j,Nn-1,j)/Ug.submat(jn,j-1,Nn-1,j-1))/pow(i*wn.rows(jn,Nn-1)-w0l,2);
				Kc_g.submat(jn,j,Nn-1,j)=log(Ug.submat(jn,j,Nn-1,j)/Ug.submat(jn,j-1,Nn-1,j-1))/(i*wn.rows(jn,Nn-1)-w0l);
				Kd_g.submat(jn,j,Nn-1,j)=zeros<cx_vec>(Nn-jn);
				for (p=pngmax; p>=1; p--)
				{
					Ka_g.submat(jn,j,Nn-1,j)=Ka_g.submat(jn,j,Nn-1,j)+(pow(i,p)*pow(-1,p+1)*(pow(Wg.submat(jn,j,Nn-1,j),p)-pow(Wg.submat(jn,j-1,Nn-1,j-1),p))/(p*pow(wn.rows(jn,Nn-1),p)))/pow(i*wn.rows(jn,Nn-1)-w0l,3);
					Kb_g.submat(jn,j,Nn-1,j)=Kb_g.submat(jn,j,Nn-1,j)+(pow(i,p)*pow(-1,p+1)*(pow(Wg.submat(jn,j,Nn-1,j),p)-pow(Wg.submat(jn,j-1,Nn-1,j-1),p))/(p*pow(wn.rows(jn,Nn-1),p)))/pow(i*wn.rows(jn,Nn-1)-w0l,2);
					Kc_g.submat(jn,j,Nn-1,j)=Kc_g.submat(jn,j,Nn-1,j)+(pow(i,p)*pow(-1,p+1)*(pow(Wg.submat(jn,j,Nn-1,j),p)-pow(Wg.submat(jn,j-1,Nn-1,j-1),p))/(p*pow(wn.rows(jn,Nn-1),p)))/(i*wn.rows(jn,Nn-1)-w0l);
					Kd_g.submat(jn,j,Nn-1,j)=Kd_g.submat(jn,j,Nn-1,j)+pow(i,p)*pow(-1,p+1)*(pow(Wg.submat(jn,j,Nn-1,j),p)-pow(Wg.submat(jn,j-1,Nn-1,j-1),p))/(p*pow(wn.rows(jn,Nn-1),p));
				}
			}
		}
	}
	 
	Ka_g=-Ka_g/(2*PI);
	Kb_g=-Kb_g/(2*PI);
	Kc_g=-Kc_g/(2*PI);
	Kd_g=-Kd_g/(2*PI);
	 

	cx_mat Ka_c_r=zeros<cx_mat>(Nn,Nintc);
	cx_mat Kb_c_r=zeros<cx_mat>(Nn,Nintc);
	cx_mat Kc_c_r=zeros<cx_mat>(Nn,Nintc);
	cx_mat Kd_c_r=zeros<cx_mat>(Nn,Nintc);

	mat Wnc=wn*ones<rowvec>(Nintc);
	mat Wc=ones<vec>(Nn)*wc.t();
	 
	mat logc=log((pow(Wnc,2)+pow(Wc.cols(1,Nintc),2))/(pow(Wnc,2)+pow(Wc.cols(0,Nintc-1),2)));
	mat atanc2=atan((Wnc % (Wc.cols(1,Nintc)-Wc.cols(0,Nintc-1)))/(Wc.cols(1,Nintc) % Wc.cols(0,Nintc-1)+pow(Wnc,2)));
	cx_mat logc2=logc/2+i*atanc2;
	cx_mat dWn=i*Wnc-Wc.cols(0,Nintc-1);
	mat dWc=Wc.cols(1,Nintc)-Wc.cols(0,Nintc-1);
	 
	Ka_c_r=-pow(dWn,2) % dWc-dWn % pow(dWc,2)/2-pow(dWc,3)/3-pow(dWn,3) % logc2;
	Kb_c_r=-dWn % dWc-pow(dWc,2)/2-pow(dWn,2) % logc2;
	Kc_c_r=-dWc-dWn % logc2;
	Kd_c_r=-logc2;
	
	int Pmax=2*pnmax+4;
	imat MP;
	pascal(Pmax+1,MP);

	if (use_HF_exp)
	{
		double wtmp;
		int jni, jnr, p, l;
		for (j=0; j<Nwc-1; j++)
		{
			wtmp=abs(wc(j));
			if (abs(wc(j+1))>wtmp) wtmp=abs(wc(j+1));
			jni=0;
			while (wn(jni)<fi*wtmp && jni<Nn-1) jni++;
			while (pow(wn(jni),2*pnmax+1)==0 && jni<Nn-1) jni++;
			jnr=1;
			while (wn(jnr)<fr*wtmp && jnr<Nn-1) jnr++;
			while (pow(wn(jnr),2*pnmax)==0 && jnr<Nn-1) jnr++;
			if (jni<Nn-1 || jnr<Nn-1)
			{
				double wj=wc(j);
				double dw1=wc(j+1)-wj;
				vec dwp=zeros<vec>(Pmax);
				dwp(0)=dw1;
				dwp(1)=pow(dw1,2)+2*wj*dw1;
				for (p=3; p<=Pmax; p++)
				{
					dwp(p-1)=0;
					for (l=0; l<p; l++)
					{
						dwp(p-1)=dwp(p-1)+MP(l,p-l)*pow(dw1,p-l)*pow(wj,l);
					}
				}
				cx_vec vtmp;
				if (jni<Nn)
				{
					vtmp.zeros(Nn-jni);
					vtmp.set_real(real(Ka_c_r.submat(jni,j,Nn-1,j)));
					Ka_c_r.submat(jni,j,Nn-1,j)=vtmp;
					vtmp.set_real(real(Kb_c_r.submat(jni,j,Nn-1,j)));
					Kb_c_r.submat(jni,j,Nn-1,j)=vtmp;
					vtmp.set_real(real(Kc_c_r.submat(jni,j,Nn-1,j)));
					Kc_c_r.submat(jni,j,Nn-1,j)=vtmp;
					vtmp.set_real(real(Kd_c_r.submat(jni,j,Nn-1,j)));
					Kd_c_r.submat(jni,j,Nn-1,j)=vtmp;

					for (p=2*pnmax+1; p>=1; p-=2)
					{
						Ka_c_r.submat(jni,j,Nn-1,j)=Ka_c_r.submat(jni,j,Nn-1,j) + i*pow(-1,(p-1)/2)*(pow(wj,3)*dwp(p-1)/p - 3.0*pow(wj,2)*dwp(p)/(p+1) + 3.0*wj*dwp(p+1)/(p+2) - dwp(p+2)/(p+3))/pow(wn.rows(jni,Nn-1),p);
						Kb_c_r.submat(jni,j,Nn-1,j)=Kb_c_r.submat(jni,j,Nn-1,j) + i*pow(-1,(p-1)/2)*(-pow(wj,2)*dwp(p-1)/p + 2*wj*dwp(p)/(p+1) - dwp(p+1)/(p+2))/pow(wn.rows(jni,Nn-1),p);
						Kc_c_r.submat(jni,j,Nn-1,j)=Kc_c_r.submat(jni,j,Nn-1,j) + i*pow(-1,(p-1)/2)*(wj*dwp(p-1)/p - dwp(p)/(p+1))/pow(wn.rows(jni,Nn-1),p);
						Kd_c_r.submat(jni,j,Nn-1,j)=Kd_c_r.submat(jni,j,Nn-1,j) - i*pow(-1,(p-1)/2)*dwp(p-1)/(p*pow(wn.rows(jni,Nn-1),p));
					}
				}
				if (jnr<Nn)
				{
					vtmp.zeros(Nn-jnr);
					vtmp.set_imag(imag(Ka_c_r.submat(jnr,j,Nn-1,j)));
					Ka_c_r.submat(jnr,j,Nn-1,j)=vtmp;
					vtmp.set_imag(imag(Kb_c_r.submat(jnr,j,Nn-1,j)));
					Kb_c_r.submat(jnr,j,Nn-1,j)=vtmp;
					vtmp.set_imag(imag(Kc_c_r.submat(jnr,j,Nn-1,j)));
					Kc_c_r.submat(jnr,j,Nn-1,j)=vtmp;
					vtmp.set_imag(imag(Kd_c_r.submat(jnr,j,Nn-1,j)));
					Kd_c_r.submat(jnr,j,Nn-1,j)=vtmp;
					for (p=2*pnmax; p>=2; p-=2)
					{
						Ka_c_r.submat(jnr,j,Nn-1,j)=Ka_c_r.submat(jnr,j,Nn-1,j) + pow(-1,(p-2)/2)*(pow(wj,3)*dwp(p-1)/p - 3*pow(wj,2)*dwp(p)/(p+1) + 3*wj*dwp(p+1)/(p+2) - dwp(p+2)/(p+3))/pow(wn.rows(jnr,Nn-1),p);
						Kb_c_r.submat(jnr,j,Nn-1,j)=Kb_c_r.submat(jnr,j,Nn-1,j) + pow(-1,(p-2)/2)*(-pow(wj,2)*dwp(p-1)/p + 2*wj*dwp(p)/(p+1) - dwp(p+1)/(p+2))/pow(wn.rows(jnr,Nn-1),p);
						Kc_c_r.submat(jnr,j,Nn-1,j)=Kc_c_r.submat(jnr,j,Nn-1,j) + pow(-1,(p-2)/2)*(wj*dwp(p-1)/p - dwp(p)/(p+1))/pow(wn.rows(jnr,Nn-1),p);
						Kd_c_r.submat(jnr,j,Nn-1,j)=Kd_c_r.submat(jnr,j,Nn-1,j) - pow(-1,(p-2)/2)*dwp(p-1)/(p*pow(wn.rows(jnr,Nn-1),p));
					}
				}
			}
	 	}
	}
	
	Ka_c_r=Ka_c_r/(2*PI);
	Kb_c_r=Kb_c_r/(2*PI);
	Kc_c_r=Kc_c_r/(2*PI);
	Kd_c_r=Kd_c_r/(2*PI);

	cx_mat Ka_d=zeros<cx_mat>(Nn,Nintd);
	cx_mat Kb_d=zeros<cx_mat>(Nn,Nintd);
	cx_mat Kc_d=zeros<cx_mat>(Nn,Nintd);
	cx_mat Kd_d=zeros<cx_mat>(Nn,Nintd);
	
	rowvec ud2(Nud+1);
	ud2.cols(1,Nud)=ud.t();
	ud2(0)=1.0/(wr-w0r);
	 
	mat Wnd=wn*ones<rowvec>(Nintd-1);
	mat Ud=ones<vec>(Nn)*ud2;

	mat atand=atan((Wnd % (Ud.cols(1,Nintd-1)-Ud.cols(0,Nintd-2)))/(1.0+w0r*(Ud.cols(1,Nintd-1)+Ud.cols(0,Nintd-2))+(pow(w0r,2)+pow(Wnd,2)) % Ud.cols(0,Nintd-2) % Ud.cols(1,Nintd-1)));
	mat logd=log((1+2*w0r*Ud.cols(1,Nintd-1)+pow(Ud.cols(1,Nintd-1),2) % (pow(Wnd,2)+pow(w0r,2)))/(1+2*w0r*Ud.cols(0,Nintd-2)+pow(Ud.cols(0,Nintd-2),2) % (pow(Wnd,2)+pow(w0r,2))));
	 
	Ka_d.cols(0,Nintd-2)=-(Ud.cols(1,Nintd-1)-Ud.cols(0,Nintd-2))/pow(Wnd+i*w0r,2) - i*(pow(Ud.cols(1,Nintd-1),2)-pow(Ud.cols(0,Nintd-2),2))/(2*(Wnd+i*w0r))+atand/pow(Wnd+i*w0r,3) +i*logd/(2*pow(Wnd+i*w0r,3));
	 
	Kb_d.cols(0,Nintd-2)=-i*(Ud.cols(1,Nintd-1)-Ud.cols(0,Nintd-2))/(Wnd+i*w0r)+i*atand/pow(Wnd+i*w0r,2)-logd/(2*pow(Wnd+i*w0r,2));
	 
	Kc_d.cols(0,Nintd-2)=-atand/(Wnd+i*w0r)-i*logd/(2*(Wnd+i*w0r));
	 
	Kd_d.cols(0,Nintd-2)=-i*atand-log(Ud.cols(1,Nintd-1)/Ud.cols(0,Nintd-2))+logd/2;
	 
	Ka_d.col(Nintd-1)=ud(Nud-1)/pow(wn+i*w0r,2)+i*pow(ud(Nud-1),2)/(2*(wn+i*w0r))-atan((ud(Nud-1)*wn)/(1+ud(Nud-1)*w0r))/pow(wn+i*w0r,3)-i*log(1+2*ud(Nud-1)*w0r+pow(ud(Nud-1),2)*(pow(wn,2)+pow(w0r,2)))/(2*pow(wn+i*w0r,3));
	Kb_d.col(Nintd-1)=i*ud(Nud-1)/(wn+i*w0r)-i*atan((ud(Nud-1)*wn)/(1+ud(Nud-1)*w0r))/pow(wn+i*w0r,2)+log(1+2*ud(Nud-1)*w0r+pow(ud(Nud-1),2)*(pow(wn,2)+pow(w0r,2)))/(2*pow(wn+i*w0r,2));
	 
	rowvec wd=1/ud2+w0r;
	mat Wd=ones<vec>(Nn)*wd;
	 
	if (use_HF_exp)
	{
		double utmp;
		int jn, p;
		for (j=0; j<Nud; j++)
		{
		 	utmp=ud2(j);
			jn=0;
			while (abs(wn(jn)*utmp)<fd*abs(1+utmp*w0r) && jn<Nn-1) jn++;
			while (pow(wn(jn),pndmax) && jn<Nn-1) jn++;
			if (jn<Nn-1)
		 	{
				Ka_d.submat(jn,j,Nn-1,j)=-(Ud.submat(jn,j+1,Nn-1,j+1)-Ud.submat(jn,j,Nn-1,j))/pow(wn.rows(jn,Nn-1)+i*w0r,2) - i*(pow(Ud.submat(jn,j+1,Nn-1,j+1),2)-pow(Ud.submat(jn,j,Nn-1,j),2))/(2*(wn.rows(jn,Nn-1)+i*w0r))+log(Ud.submat(jn,j+1,Nn-1,j+1)/Ud.submat(jn,j,Nn-1,j))/pow(i*wn.rows(jn,Nn-1)-w0r,3);
				Kb_d.submat(jn,j,Nn-1,j)=-i*(Ud.submat(jn,j+1,Nn-1,j+1)-Ud.submat(jn,j,Nn-1,j))/(wn.rows(jn,Nn-1)+i*w0r)+log(Ud.submat(jn,j+1,Nn-1,j+1)/Ud.submat(jn,j,Nn-1,j))/pow(i*wn.rows(jn,Nn-1)-w0r,2);
				Kc_d.submat(jn,j,Nn-1,j)=log(Ud.submat(jn,j+1,Nn-1,j+1)/Ud.submat(jn,j,Nn-1,j))/(i*wn.rows(jn,Nn-1)-w0r);
				Kd_d.submat(jn,j,Nn-1,j)=zeros<cx_vec>(Nn-jn);
				for (p=pndmax; p>=1; p--)
			 	{
					Ka_d.submat(jn,j,Nn-1,j)=Ka_d.submat(jn,j,Nn-1,j)+(pow(i,p)*pow(-1,p+1)*(pow(Wd.submat(jn,j+1,Nn-1,j+1),p)-pow(Wd.submat(jn,j,Nn-1,j),p))/(p*pow(wn.rows(jn,Nn-1),p)))/pow(i*wn.rows(jn,Nn-1)-w0r,3);
					Kb_d.submat(jn,j,Nn-1,j)=Kb_d.submat(jn,j,Nn-1,j)+(pow(i,p)*pow(-1,p+1)*(pow(Wd.submat(jn,j+1,Nn-1,j+1),p)-pow(Wd.submat(jn,j,Nn-1,j),p))/(p*pow(wn.rows(jn,Nn-1),p)))/pow(i*wn.rows(jn,Nn-1)-w0r,2);
					Kc_d.submat(jn,j,Nn-1,j)=Kc_d.submat(jn,j,Nn-1,j)+(pow(i,p)*pow(-1,p+1)*(pow(Wd.submat(jn,j+1,Nn-1,j+1),p)-pow(Wd.submat(jn,j,Nn-1,j),p))/(p*pow(wn.rows(jn,Nn-1),p)))/(i*wn.rows(jn,Nn-1)-w0r);
					Kd_d.submat(jn,j,Nn-1,j)=Kd_d.submat(jn,j,Nn-1,j)+pow(i,p)*pow(-1,p+1)*(pow(Wd.submat(jn,j+1,Nn-1,j+1),p)-pow(Wd.submat(jn,j,Nn-1,j),p))/(p*pow(wn.rows(jn,Nn-1),p));
			 	}
		 	}
		}
	}
	
	Ka_d=-Ka_d/(2*PI);
	Kb_d=-Kb_d/(2*PI);
	Kc_d=-Kc_d/(2*PI);
	Kd_d=-Kd_d/(2*PI);
	

	cx_mat KG=(Ka_g*Pa_g+Kb_g*Pb_g+Kc_g*Pc_g+Kd_g*Pd_g)*MC;
	cx_mat KC=(Ka_c_r*Pa_c+Kb_c_r*Pb_c_r+Kc_c_r*Pc_c_r+Kd_c_r*Pd_c_r)*MC;
	cx_mat KD=(Ka_d*Pa_d+Kb_d*Pb_d+Kc_d*Pc_d+Kd_d*Pd_d)*MC;
	 
	Kcx=KG+KC+KD;


	rowvec Knorm_a_g=zeros<rowvec>(Nintg);
	rowvec Knorm_b_g=zeros<rowvec>(Nintg);
	rowvec Knorm_c_g=zeros<rowvec>(Nintg);
	rowvec Knorm_d_g=zeros<rowvec>(Nintg);
	 
	Knorm_a_g(0)=-pow(ug(0),2)/2;
	Knorm_b_g(0)=-ug(0);
	 
	Knorm_a_g.cols(1,Nintg-1)=-(pow(ug2.cols(1,Nug),2)-pow(ug2.cols(0,Nug-1),2))/2;
	Knorm_b_g.cols(1,Nintg-1)=-(ug2.cols(1,Nug)-ug2.cols(0,Nug-1));
	Knorm_c_g.cols(1,Nintg-1)=-log(ug2.cols(1,Nug)/ug2.cols(0,Nug-1));
	Knorm_d_g.cols(1,Nintg-1)=1.0/ug2.cols(1,Nug)-1.0/ug2.cols(0,Nug-1);
	 
	rowvec KM0g=(Knorm_a_g*Pa_g+Knorm_b_g*Pb_g+Knorm_c_g*Pc_g+Knorm_d_g*Pd_g)*MC/(2*PI);

	rowvec Knorm_a_c_r=zeros<rowvec>(Nintc);
	rowvec Knorm_b_c_r=zeros<rowvec>(Nintc);
	rowvec Knorm_c_c_r=zeros<rowvec>(Nintc);
	rowvec Knorm_d_c_r=zeros<rowvec>(Nintc);
	 
	Knorm_a_c_r=trans(pow(wc.rows(1,Nwc-1)-wc.rows(0,Nwc-2),4)/4);
	Knorm_b_c_r=trans(pow(wc.rows(1,Nwc-1)-wc.rows(0,Nwc-2),3)/3);
	Knorm_c_c_r=trans(pow(wc.rows(1,Nwc-1)-wc.rows(0,Nwc-2),2)/2);
	Knorm_d_c_r=trans(wc.rows(1,Nwc-1)-wc.rows(0,Nwc-2));
	 
	rowvec KM0c=(Knorm_a_c_r*Pa_c+Knorm_b_c_r*Pb_c_r+Knorm_c_c_r*Pc_c_r+Knorm_d_c_r*Pd_c_r)*MC/(2*PI);

	rowvec Knorm_a_d=zeros<rowvec>(Nintd);
	rowvec Knorm_b_d=zeros<rowvec>(Nintd);
	rowvec Knorm_c_d=zeros<rowvec>(Nintd);
	rowvec Knorm_d_d=zeros<rowvec>(Nintd);
	
	Knorm_a_d.cols(0,Nintd-2)=-(pow(ud2.cols(1,Nud),2)-pow(ud2.cols(0,Nud-1),2))/2;
	Knorm_b_d.cols(0,Nintd-2)=-(ud2.cols(1,Nud)-ud2.cols(0,Nud-1));
	Knorm_c_d.cols(0,Nintd-2)=-log(ud2.cols(1,Nud)/ud2.cols(0,Nud-1));
	Knorm_d_d.cols(0,Nintd-2)=1.0/ud2.cols(1,Nud)-1.0/ud2.cols(0,Nud-1);
	
	Knorm_a_d(Nintd-1)=pow(ud2(Nud-1),2)/2;
	Knorm_b_d(Nintd-1)=ud2(Nud-1);
	
	rowvec KM0d=(Knorm_a_d*Pa_d+Knorm_b_d*Pb_d+Knorm_c_d*Pc_d+Knorm_d_d*Pd_d)*MC/(2*PI);
	
	rowvec KM0=KM0g+KM0c+KM0d;
	

	rowvec KM1_a_g=zeros<rowvec>(Nintg);
	rowvec KM1_b_g=zeros<rowvec>(Nintg);
	rowvec KM1_c_g=zeros<rowvec>(Nintg);
	rowvec KM1_d_g=zeros<rowvec>(Nintg);
	
	KM1_a_g.cols(1,Nintg-1)=Knorm_b_g.cols(1,Nintg-1);
	KM1_b_g.cols(1,Nintg-1)=Knorm_c_g.cols(1,Nintg-1);
	KM1_c_g.cols(1,Nintg-1)=Knorm_d_g.cols(1,Nintg-1);
	KM1_d_g.cols(1,Nintg-1)=(1.0/pow(ug2.cols(1,Nug),2)-1.0/pow(ug2.cols(0,Nug-1),2))/2;
	
	rowvec KM1g_tmp=(KM1_a_g*Pa_g+KM1_b_g*Pb_g+KM1_c_g*Pc_g+KM1_d_g*Pd_g)*MC/(2*PI);
	
	rowvec KM1g=w0l*KM0g+KM1g_tmp;

	mat Wjc=diagmat(wc.rows(0,Nintc-1));
	
	rowvec KM1_a_c_tmp=trans(pow(wc.rows(1,Nwc-1)-wc.rows(0,Nwc-2),5)/5.0);
	
	rowvec KM1_a_c=KM1_a_c_tmp+Knorm_a_c_r*Wjc;
	rowvec KM1_b_c=Knorm_a_c_r+Knorm_b_c_r*Wjc;
	rowvec KM1_c_c=Knorm_b_c_r+Knorm_c_c_r*Wjc;
	rowvec KM1_d_c=Knorm_c_c_r+Knorm_d_c_r*Wjc;
	
	rowvec KM1c=(KM1_a_c*Pa_c+KM1_b_c*Pb_c_r+KM1_c_c*Pc_c_r+KM1_d_c*Pd_c_r)*MC/(2*PI);
	 
	rowvec KM1_a_d=zeros<rowvec>(Nintd);
	rowvec KM1_b_d=zeros<rowvec>(Nintd);
	rowvec KM1_c_d=zeros<rowvec>(Nintd);
	rowvec KM1_d_d=zeros<rowvec>(Nintd);
	 
	KM1_a_d.cols(0,Nintd-2)=Knorm_b_d.cols(0,Nintd-2);
	KM1_b_d.cols(0,Nintd-2)=Knorm_c_d.cols(0,Nintd-2);
	KM1_c_d.cols(0,Nintd-2)=Knorm_d_d.cols(0,Nintd-2);
	KM1_d_d.cols(0,Nintd-2)=(1.0/pow(ud2.cols(1,Nud),2)-1.0/pow(ud2.cols(0,Nud-1),2))/2;
	
	rowvec KM1d_tmp=(KM1_a_d*Pa_d+KM1_b_d*Pb_d+KM1_c_d*Pc_d+KM1_d_d*Pd_d)*MC/(2*PI);
	 
	rowvec KM1d=w0r*KM0d+KM1d_tmp;
	 
	rowvec KM1=KM1g+KM1c+KM1d;


	rowvec KM2_a_g=zeros<rowvec>(Nintg);
	rowvec KM2_b_g=zeros<rowvec>(Nintg);
	rowvec KM2_c_g=zeros<rowvec>(Nintg);
	rowvec KM2_d_g=zeros<rowvec>(Nintg);
	 
	KM2_a_g.cols(1,Nintg-1)=Knorm_c_g.cols(1,Nintg-1);
	KM2_b_g.cols(1,Nintg-1)=Knorm_d_g.cols(1,Nintg-1);
	KM2_c_g.cols(1,Nintg-1)=KM1_d_g.cols(1,Nintg-1);
	KM2_d_g.cols(1,Nintg-1)=(1.0/pow(ug2.cols(1,Nug),3)-1.0/pow(ug2.cols(0,Nug-1),3))/3;
	 
	rowvec KM2g_tmp=(KM2_a_g*Pa_g+KM2_b_g*Pb_g+KM2_c_g*Pc_g+KM2_d_g*Pd_g)*MC/(2*PI);
	 
	rowvec KM2g=pow(w0l,2)*KM0g+2*w0l*KM1g_tmp+KM2g_tmp;
	 
	rowvec KM2_a_c_tmp=trans(pow(wc.rows(1,Nwc-1)-wc.rows(0,Nwc-2),6)/6);
	 
	rowvec KM2_a_c=KM2_a_c_tmp+2*KM1_a_c_tmp*Wjc+Knorm_a_c_r*pow(Wjc,2);
	rowvec KM2_b_c=KM1_a_c_tmp+2*Knorm_a_c_r*Wjc+Knorm_b_c_r*pow(Wjc,2);
	rowvec KM2_c_c=Knorm_a_c_r+2*Knorm_b_c_r*Wjc+Knorm_c_c_r*pow(Wjc,2);
	rowvec KM2_d_c=Knorm_b_c_r+2*Knorm_c_c_r*Wjc+Knorm_d_c_r*pow(Wjc,2);
	 
	rowvec KM2c=(KM2_a_c*Pa_c+KM2_b_c*Pb_c_r+KM2_c_c*Pc_c_r+KM2_d_c*Pd_c_r)*MC/(2*PI);
	 
	rowvec KM2_a_d=zeros<rowvec>(Nintd);
	rowvec KM2_b_d=zeros<rowvec>(Nintd);
	rowvec KM2_c_d=zeros<rowvec>(Nintd);
	rowvec KM2_d_d=zeros<rowvec>(Nintd);
	 
	KM2_a_d.cols(0,Nintd-2)=Knorm_c_d.cols(0,Nintd-2);
	KM2_b_d.cols(0,Nintd-2)=Knorm_d_d.cols(0,Nintd-2);
	KM2_c_d.cols(0,Nintd-2)=KM1_d_d.cols(0,Nintd-2);
	KM2_d_d.cols(0,Nintd-2)=(1.0/pow(ud2.cols(1,Nud),3)-1.0/pow(ud2.cols(0,Nud-1),3))/3;
	 
	rowvec KM2d_tmp=(KM2_a_d*Pa_d+KM2_b_d*Pb_d+KM2_c_d*Pc_d+KM2_d_d*Pd_d)*MC/(2*PI);
	 
	rowvec KM2d=pow(w0r,2)*KM0d+2*w0r*KM1d_tmp+KM2d_tmp;
	 
	rowvec KM2=KM2g+KM2c+KM2d;

	rowvec KM3_a_g=zeros<rowvec>(Nintg);
	rowvec KM3_b_g=zeros<rowvec>(Nintg);
	rowvec KM3_c_g=zeros<rowvec>(Nintg);
	rowvec KM3_d_g=zeros<rowvec>(Nintg);
	 
	KM3_a_g.cols(1,Nintg-1)=Knorm_d_g.cols(1,Nintg-1);
	KM3_b_g.cols(1,Nintg-1)=KM1_d_g.cols(1,Nintg-1);
	KM3_c_g.cols(1,Nintg-1)=KM2_d_g.cols(1,Nintg-1);
	KM3_d_g.cols(1,Nintg-1)=(1.0/pow(ug2.cols(1,Nug),4)-1.0/pow(ug2.cols(0,Nug-1),4))/4;
	 
	rowvec KM3g_tmp=(KM3_a_g*Pa_g+KM3_b_g*Pb_g+KM3_c_g*Pc_g+KM3_d_g*Pd_g)*MC/(2*PI);
	 
	rowvec KM3g=pow(w0l,3)*KM0g+3*pow(w0l,2)*KM1g_tmp+3*w0l*KM2g_tmp+KM3g_tmp;
	 
	rowvec KM3_a_c_tmp=trans(pow(wc.rows(1,Nwc-1)-wc.rows(0,Nwc-2),7)/7);
	 
	rowvec KM3_a_c=KM3_a_c_tmp+3*KM2_a_c_tmp*Wjc+3*KM1_a_c_tmp*pow(Wjc,2)+Knorm_a_c_r*pow(Wjc,3);
	rowvec KM3_b_c=KM2_a_c_tmp+3*KM1_a_c_tmp*Wjc+3*Knorm_a_c_r*pow(Wjc,2)+Knorm_b_c_r*pow(Wjc,3);
	rowvec KM3_c_c=KM1_a_c_tmp+3*Knorm_a_c_r*Wjc+3*Knorm_b_c_r*pow(Wjc,2)+Knorm_c_c_r*pow(Wjc,3);
	rowvec KM3_d_c=Knorm_a_c_r+3*Knorm_b_c_r*Wjc+3*Knorm_c_c_r*pow(Wjc,2)+Knorm_d_c_r*pow(Wjc,3);
	 
	rowvec KM3c=(KM3_a_c*Pa_c+KM3_b_c*Pb_c_r+KM3_c_c*Pc_c_r+KM3_d_c*Pd_c_r)*MC/(2*PI);
	 
	rowvec KM3_a_d=zeros<rowvec>(Nintd);
	rowvec KM3_b_d=zeros<rowvec>(Nintd);
	rowvec KM3_c_d=zeros<rowvec>(Nintd);
	rowvec KM3_d_d=zeros<rowvec>(Nintd);
	 
	KM3_a_d.cols(0,Nintd-2)=Knorm_d_d.cols(0,Nintd-2);
	KM3_b_d.cols(0,Nintd-2)=KM1_d_d.cols(0,Nintd-2);
	KM3_c_d.cols(0,Nintd-2)=KM2_d_d.cols(0,Nintd-2);
	KM3_d_d.cols(0,Nintd-2)=(1.0/pow(ud2.cols(1,Nud),4)-1.0/pow(ud2.cols(0,Nud-1),4))/4;
	 
	rowvec KM3d_tmp=(KM3_a_d*Pa_d+KM3_b_d*Pb_d+KM3_c_d*Pc_d+KM3_d_d*Pd_d)*MC/(2*PI);
	 
	rowvec KM3d=pow(w0r,3)*KM0d+3*pow(w0r,2)*KM1d_tmp+3*w0r*KM2d_tmp+KM3d_tmp;
	 
	rowvec KM3=KM3g+KM3c+KM3d;

	KM.row(0)=KM0;
	KM.row(1)=KM1;
	KM.row(2)=KM2;
	KM.row(3)=KM3;
	
	K.zeros(2*Nn,Nw);
	uvec even_ind=linspace<uvec>(0,2*Nn-2,Nn);
	K.rows(even_ind)=real(Kcx);
	K.rows(even_ind+1)=imag(Kcx);

	
//	cout<<"kernel matrix defined.\n";
	
	return true;
}



bool OmegaMaxEnt_data::set_A_ref()
{
	if (Aref_data.n_cols<2)
	{
		cout<<"number of columns in file "<<A_ref_file<<" is smaller than 2\n";
		return false;
	}
	
	w_ref=Aref_data.col(0);
	A_ref=Aref_data.col(1);
	
	return true;
}




double OmegaMaxEnt_data::spline_val_G_part_int(double x, void *par[])
{
	vec *x0=reinterpret_cast<vec*>(par[0]);
	uvec *ind_xlims=reinterpret_cast<uvec*>(par[1]);
	vec *xs=reinterpret_cast<vec*>(par[2]);
	vec *coeffs=reinterpret_cast<vec*>(par[3]);
	int *p=reinterpret_cast<int*>(par[4]);
	
	return pow(x,p[0])*spline_val_G_part(x, *x0, *ind_xlims, *xs, *coeffs);
}

double OmegaMaxEnt_data::spline_val_G_part(double x, vec &x0, uvec &ind_xlims, vec &xs, vec &coeffs)
{
	int Nx0=x0.n_rows;
	
	double xl=x0(ind_xlims(0));
	double xr=x0(ind_xlims(1));
	
	double x0l=xs(0);
	double x0r=xs(1);
	
	double a,b,c,d, Dx, Du, sv;
	
	int l;
	sv=0;
	
	//if (x>=x0(0) && x<=x0(Nx0-1))
	//{
		if (x<xl)
		{
			l=0;
			while (x>=x0(l) && l<ind_xlims(0)) l++;
			a=coeffs(4*l);
			b=coeffs(4*l+1);
			c=coeffs(4*l+2);
			d=coeffs(4*l+3);
			if (l>0)
				Du=(x0(l-1)-x)/((x-x0l)*(x0(l-1)-x0l));
			else
				Du=1.0/(x-x0l);
			sv=a*pow(Du,3)+b*pow(Du,2)+c*Du+d;
		}
		else if (x>=xl && x<xr)
		{
			l=ind_xlims(0)+1;
			while (x>=x0(l) && l<ind_xlims(1)) l++;
			a=coeffs(4*l);
			b=coeffs(4*l+1);
			c=coeffs(4*l+2);
			d=coeffs(4*l+3);
			Dx=x-x0(l-1);
			sv=a*pow(Dx,3)+b*pow(Dx,2)+c*Dx+d;
		}
		else
		{
			l=ind_xlims(1)+1;
			while (x>x0(l) && l<Nx0-1) l++;
			if (x>x0(Nx0-1)) l=Nx0;
			a=coeffs(4*l);
			b=coeffs(4*l+1);
			c=coeffs(4*l+2);
			d=coeffs(4*l+3);
			if (l<Nx0)
				Du=(x0(l)-x)/((x-x0r)*(x0(l)-x0r));
			else
				Du=1.0/(x-x0r);
			sv=a*pow(Du,3)+b*pow(Du,2)+c*Du+d;
		}
	//}
	
	return sv;
}

bool OmegaMaxEnt_data::spline_val_G_part_grid_transf(vec x, vec x0, uvec ind_xlims, vec xs, vec coeffs, vec &sv)
{
	int Nx=x.n_rows;
	int Nx0=x0.n_rows;
	
	int Ncfs=4*(Nx0-1);
	if (coeffs.n_rows<Ncfs)
	{
		cout<<"spline_val_G_part_grid_transf(): number of elements in vectors \"coeffs\" and \"x0\" are not consistant\n";
		return false;
	}
	
	double xl=x0(ind_xlims(0));
	double xr=x0(ind_xlims(1));
	
	double x0l=xs(0);
	double x0r=xs(1);
	
	double a,b,c,d, Dx, Du;
	
	sv.zeros(Nx);
	
	int j, l;
	for (j=0; j<Nx; j++)
	{
		if (x(j)>=x0(0) && x(j)<=x0(Nx0-1))
		{
			if (x(j)<xl)
			{
				l=0;
				while (x(j)>=x0(l+1) && l<ind_xlims(0)-1) l++;
				a=coeffs(4*l);
				b=coeffs(4*l+1);
				c=coeffs(4*l+2);
				d=coeffs(4*l+3);
				
				Du=((x0(l+1)-x0l)/(x(j)-x0l))*((x0(l)-x(j))/(x0(l)-x0(l+1)));
				
				sv(j)=a*pow(Du,3)+b*pow(Du,2)+c*Du+d;
			}
			else if (x(j)>=xl && x(j)<=xr)
			{
				if (x(j)<xr)
				{
					l=ind_xlims(0);
					while (x(j)>=x0(l+1) && l<ind_xlims(1)-1) l++;
					a=coeffs(4*l);
					b=coeffs(4*l+1);
					c=coeffs(4*l+2);
					d=coeffs(4*l+3);
					Dx=(x(j)-x0(l))/(x0(l+1)-x0(l));
					sv(j)=a*pow(Dx,3)+b*pow(Dx,2)+c*Dx+d;
				}
				else
				{
					l=ind_xlims(1)-1;
					a=coeffs(4*l);
					b=coeffs(4*l+1);
					c=coeffs(4*l+2);
					d=coeffs(4*l+3);
					sv(j)=a+b+c+d;
				}
			}
			else
			{
				l=ind_xlims(1);
				while (x(j)>x0(l+1) && l<Nx0-2) l++;
				a=coeffs(4*l);
				b=coeffs(4*l+1);
				c=coeffs(4*l+2);
				d=coeffs(4*l+3);
				Du=((x0(l+1)-x0r)/(x(j)-x0r))*((x0(l)-x(j))/(x0(l)-x0(l+1)));
				sv(j)=a*pow(Du,3)+b*pow(Du,2)+c*Du+d;
			}
		}
	}
	
	return true;
}

bool OmegaMaxEnt_data::spline_val_G_part_grid_transf_1(vec x, vec x0, uvec ind_xlims, vec xs, vec coeffs, vec &sv)
{
	int Nx=x.n_rows;
	int Nx0=x0.n_rows;
	
	int Ncfs=4*(Nx0-1);
	if (coeffs.n_rows<Ncfs)
	{
		cout<<"spline_val_G_part_grid_transf(): number of elements in vectors \"coeffs\" and \"x0\" are not consistant\n";
		return false;
	}
	
	double xl=x0(ind_xlims(0));
	double xr=x0(ind_xlims(1));
	
	double x0l=xs(0);
	double x0r=xs(1);
	
	double a,b,c,d, Dx, Du;
	
	sv.zeros(Nx);
	
	int j, l;
	for (j=0; j<Nx; j++)
	{
		if (x(j)>=x0(0) && x(j)<=x0(Nx0-1))
		{
			if (x(j)<xl)
			{
				l=0;
				while (x(j)>=x0(l+1) && l<ind_xlims(0)-1) l++;
				a=coeffs(4*l);
				b=coeffs(4*l+1);
				c=coeffs(4*l+2);
				d=coeffs(4*l+3);
				
				Du=((x0(l+1)-x0l)/(x(j)-x0l))*((x0(l)-x(j))/(x0(l)-x0(l+1)));
				
				sv(j)=a*pow(Du,3)+b*pow(Du,2)+c*Du+d;
			}
			else if (x(j)>=xl && x(j)<=xr)
			{
				if (x(j)<xr)
				{
					l=ind_xlims(0);
					while (x(j)>=x0(l+1) && l<ind_xlims(1)-1) l++;
					a=coeffs(4*l);
					b=coeffs(4*l+1);
					c=coeffs(4*l+2);
					d=coeffs(4*l+3);
					Dx=(x(j)-x0(l))/(x0(l+1)-x0(l));
					sv(j)=a*pow(Dx,3)+b*pow(Dx,2)+c*Dx+d;
				}
				else
				{
					l=ind_xlims(1)-1;
					a=coeffs(4*l);
					b=coeffs(4*l+1);
					c=coeffs(4*l+2);
					d=coeffs(4*l+3);
					sv(j)=a+b+c+d;
				}
			}
			else
			{
				l=ind_xlims(1);
				while (x(j)>x0(l+1) && l<Nx0-2) l++;
				a=coeffs(4*l);
				b=coeffs(4*l+1);
				c=coeffs(4*l+2);
				d=coeffs(4*l+3);
				Du=((x0(l)-x0r)/(x(j)-x0r))*((x0(l+1)-x(j))/(x0(l+1)-x0(l)));
				sv(j)=a*pow(Du,3)+b*pow(Du,2)+c*Du+d;
			}
		}
	}
	
	return true;
}

bool OmegaMaxEnt_data::spline_val_G_part(vec x, vec x0, uvec ind_xlims, vec xs, vec coeffs, vec &sv)
{
	int Nx=x.n_rows;
	int Nx0=x0.n_rows;
	
	int Ncfs=4*(Nx0+1);
	if (coeffs.n_rows<Ncfs)
	{
		cout<<"spline_val_G_part(): number of elements in vectors \"coeffs\" and \"x0\" are not consistant\n";
		return false;
	}
	
//	int Ng=ind_xlims(0)+1;
//	int Nc=ind_xlims(1)-ind_xlims(0);
//	int Nd=Nx0-ind_xlims(1);
	
	double xl=x0(ind_xlims(0));
	double xr=x0(ind_xlims(1));
	
	double x0l=xs(0);
	double x0r=xs(1);
	
	double a,b,c,d, Dx, Du;
	
	sv.zeros(Nx);
	
	int j, l;
	for (j=0; j<Nx; j++)
	{
		//if (x(j)>=x0(0) && x(j)<=x0(Nx0-1))
		//{
			if (x(j)<xl)
			{
				l=0;
				while (x(j)>=x0(l) && l<ind_xlims(0)) l++;
				a=coeffs(4*l);
				b=coeffs(4*l+1);
				c=coeffs(4*l+2);
				d=coeffs(4*l+3);
				if (l>0)
					Du=(x0(l-1)-x(j))/((x(j)-x0l)*(x0(l-1)-x0l));
				else
					Du=1.0/(x(j)-x0l);
				sv(j)=a*pow(Du,3)+b*pow(Du,2)+c*Du+d;
			}
			else if (x(j)>=xl && x(j)<xr)
			{
				l=ind_xlims(0)+1;
				while (x(j)>=x0(l) && l<ind_xlims(1)) l++;
				a=coeffs(4*l);
				b=coeffs(4*l+1);
				c=coeffs(4*l+2);
				d=coeffs(4*l+3);
				Dx=x(j)-x0(l-1);
				sv(j)=a*pow(Dx,3)+b*pow(Dx,2)+c*Dx+d;
			}
			else
			{
				l=ind_xlims(1)+1;
				while (x(j)>x0(l) && l<Nx0-1) l++;
				if (x(j)>x0(Nx0-1)) l=Nx0;
				a=coeffs(4*l);
				b=coeffs(4*l+1);
				c=coeffs(4*l+2);
				d=coeffs(4*l+3);
				if (l<Nx0)
					Du=(x0(l)-x(j))/((x(j)-x0r)*(x0(l)-x0r));
				else
					Du=1.0/(x(j)-x0r);
				sv(j)=a*pow(Du,3)+b*pow(Du,2)+c*Du+d;
			}
		//}
	}
	
	return true;
}

bool OmegaMaxEnt_data::spline_G_omega_u(vec x, uvec ind_xlims, vec xs, vec F, vec &coeffs)
{
	int Nx=x.n_rows;
	
	if (F.n_rows!=Nx)
	{
		cout<<"spline_G() error: function vector and position vector do not have the same size.\n";
		return false;
	}
	
	//int Ng=ind_xlims(0)+1;
	//int Nc=ind_xlims(1)-ind_xlims(0);
	//int Nd=Nx-ind_xlims(1);
	int Nint=Nx+1;
	
	double x0g=xs(0);
	double x0d=xs(1);
	
	//	 solve the spline in the interval -inf,wl]
	int NS=3*Nint-2;
	mat Aspl=zeros<mat>(NS,NS);
	
	vec D=zeros<vec>(NS);
	D(0)=F(0);
	
	double u=1.0/(x(0)-x0g);
	
	Aspl(0,0)=pow(u,3);
	Aspl(0,1)=pow(u,2);
	Aspl(1,0)=3*pow(u,2);
	Aspl(1,1)=2*u;
	Aspl(1,4)=-1;
	Aspl(2,0)=6*u;
	Aspl(2,1)=2;
	Aspl(2,3)=-2;
	
	int j;
	for (j=1; j<ind_xlims(0); j++)
	{
		u=1/(x(j)-x0g)-1/(x(j-1)-x0g);
	 
		D(3*j)=F(j)-F(j-1);
	 
		Aspl(3*j,3*j-1)=pow(u,3);
		Aspl(3*j,3*j)=pow(u,2);
		Aspl(3*j,3*j+1)=u;
	 
		Aspl(3*j+1,3*j-1)=3*pow(u,2);
		Aspl(3*j+1,3*j)=2*u;
		Aspl(3*j+1,3*j+1)=1;
		Aspl(3*j+1,3*j+4)=-1;
	 
		Aspl(3*j+2,3*j-1)=6*u;
		Aspl(3*j+2,3*j)=2;
		Aspl(3*j+2,3*j+3)=-2;
	}
	
	j=ind_xlims(0);
	
	u=1/(x(j)-x0g)-1/(x(j-1)-x0g);
	double ug=1/(x(j)-x0g);
	
	D(3*j)=F(j)-F(j-1);
	Aspl(3*j,3*j-1)=pow(u,3);
	Aspl(3*j,3*j)=pow(u,2);
	Aspl(3*j,3*j+1)=u;
	
	Aspl(3*j+1,3*j-1)=-3*pow(ug,2)*pow(u,2);
	Aspl(3*j+1,3*j)=-2*pow(ug,2)*u;
	Aspl(3*j+1,3*j+1)=-pow(ug,2);
	Aspl(3*j+1,3*j+4)=-1;
	
	Aspl(3*j+2,3*j-1)=6*pow(ug,3)*pow(u,2)+6*pow(ug,4)*u;
	Aspl(3*j+2,3*j)=4*pow(ug,3)*u+2*pow(ug,4);
	Aspl(3*j+2,3*j+1)=2*pow(ug,3);
	Aspl(3*j+2,3*j+3)=-2;
	
	double Dx;
	for (j=ind_xlims(0)+1; j<ind_xlims(1); j++)
	{
		Dx=x(j)-x(j-1);
	 
		D(3*j)=F(j)-F(j-1);
	 
		Aspl(3*j,3*j-1)=pow(Dx,3);
		Aspl(3*j,3*j)=pow(Dx,2);
		Aspl(3*j,3*j+1)=Dx;
	 
		Aspl(3*j+1,3*j-1)=3*pow(Dx,2);
		Aspl(3*j+1,3*j)=2*Dx;
		Aspl(3*j+1,3*j+1)=1;
		Aspl(3*j+1,3*j+4)=-1;
	 
		Aspl(3*j+2,3*j-1)=6*Dx;
		Aspl(3*j+2,3*j)=2;
		Aspl(3*j+2,3*j+3)=-2;
	}
	
	j=ind_xlims(1);
	
	Dx=x(j)-x(j-1);
	
	D(3*j)=F(j)-F(j-1);
	Aspl(3*j,3*j-1)=pow(Dx,3);
	Aspl(3*j,3*j)=pow(Dx,2);
	Aspl(3*j,3*j+1)=Dx;
	
	u=1/(x(j)-x0d)-1/(x(j+1)-x0d);
	double ud=1/(x(j)-x0d);
	
	Aspl(3*j+1,3*j-1)=3*pow(Dx,2);
	Aspl(3*j+1,3*j)=2*Dx;
	Aspl(3*j+1,3*j+1)=1;
	Aspl(3*j+1,3*j+2)=3*pow(ud,2)*pow(u,2);
	Aspl(3*j+1,3*j+3)=2*pow(ud,2)*u;
	Aspl(3*j+1,3*j+4)=pow(ud,2);
	
	Aspl(3*j+2,3*j-1)=6*Dx;
	Aspl(3*j+2,3*j)=2;
	Aspl(3*j+2,3*j+2)=-6*pow(ud,3)*pow(u,2)-6*pow(ud,4)*u;
	Aspl(3*j+2,3*j+3)=-4*pow(ud,3)*u-2*pow(ud,4);
	Aspl(3*j+2,3*j+4)=-2*pow(ud,3);
	
	D(3*j+3)=F(j)-F(j+1);
	Aspl(3*j+3,3*j+2)=pow(u,3);
	Aspl(3*j+3,3*j+3)=pow(u,2);
	Aspl(3*j+3,3*j+4)=u;
	
	for (j=ind_xlims(1)+1; j<Nx-1; j++)
	{
		u=1/(x(j)-x0d)-1/(x(j+1)-x0d);
		
		Aspl(3*j+1,3*j+1)=-1;
		Aspl(3*j+1,3*j+2)=3*pow(u,2);
		Aspl(3*j+1,3*j+3)=2*u;
		Aspl(3*j+1,3*j+4)=1;
		
		Aspl(3*j+2,3*j)=-2;
		Aspl(3*j+2,3*j+2)=6*u;
		Aspl(3*j+2,3*j+3)=2;
		
		D(3*j+3)=F(j)-F(j+1);
		Aspl(3*j+3,3*j+2)=pow(u,3);
		Aspl(3*j+3,3*j+3)=pow(u,2);
		Aspl(3*j+3,3*j+4)=u;
	}
	
	j=Nx-1;
	u=1/(x(j)-x0d);
	
	Aspl(3*j+1,3*j+1)=-1;
	Aspl(3*j+1,3*j+2)=3*pow(u,2);
	Aspl(3*j+1,3*j+3)=2*u;
	
	Aspl(3*j+2,3*j)=-2;
	Aspl(3*j+2,3*j+2)=6*u;
	Aspl(3*j+2,3*j+3)=2;
	
	D(3*j+3)=F(j);
	Aspl(3*j+3,3*j+2)=pow(u,3);
	Aspl(3*j+3,3*j+3)=pow(u,2);
	
	vec C=solve(Aspl,D);
	
	coeffs.zeros(4*(Nx+1));
	
	coeffs(0)=C(0);
	coeffs(1)=C(1);
	
	uvec ind_tmp=linspace<uvec>(1,Nx-1,Nx-1);
	coeffs.rows(4*ind_tmp)=C(3*ind_tmp-1);
	coeffs.rows(4*ind_tmp+1)=C(3*ind_tmp);
	coeffs.rows(4*ind_tmp+2)=C(3*ind_tmp+1);
	
	ind_tmp=linspace<uvec>(1,ind_xlims(1),ind_xlims(1));
	
	coeffs.rows(4*ind_tmp+3)=F.rows(ind_tmp-1);
	
	ind_tmp=linspace<uvec>(ind_xlims(1)+1,Nx-1,Nx-ind_xlims(1)-1);
	coeffs.rows(4*ind_tmp+3)=F.rows(ind_tmp);
	
	coeffs(4*Nx)=C(3*Nx-1);
	coeffs(4*Nx+1)=C(3*Nx);
	
	return true;
}

bool OmegaMaxEnt_data::spline_G_part(vec x, uvec ind_xlims, vec xs, vec F, vec &coeffs)
{
	int Nx=x.n_rows;
	
	if (F.n_rows!=Nx)
	{
		cout<<"spline_G_part() error: function vector and position vector do not have the same size.\n";
		return false;
	}
	
	int Ng=ind_xlims(0)+1;
	int Nc=ind_xlims(1)-ind_xlims(0);
	int Nd=Nx-ind_xlims(1);
	
	double x0g=xs(0);
	double x0d=xs(1);
	
	//	 solve the spline in the interval -inf,wl]
	int NSg=3*Ng-1;
	mat Aspl=zeros<mat>(NSg,NSg);
	
	vec Dg=zeros<vec>(NSg);
	Dg(0)=F(0);
	
	double u=1.0/(x(0)-x0g);
	
	Aspl(0,0)=pow(u,3);
	Aspl(0,1)=pow(u,2);
	Aspl(1,0)=3*pow(u,2);
	Aspl(1,1)=2*u;
	Aspl(1,4)=-1;
	Aspl(2,0)=6*u;
	Aspl(2,1)=2;
	Aspl(2,3)=-2;
	
	int j;
	for (j=2; j<Ng; j++)
	{
		u=1/(x(j-1)-x0g)-1/(x(j-2)-x0g);
	 
		Dg(3*j-3)=F(j-1)-F(j-2);
	 
		Aspl(3*j-3,3*j-4)=pow(u,3);
		Aspl(3*j-3,3*j-3)=pow(u,2);
		Aspl(3*j-3,3*j-2)=u;
	 
		Aspl(3*j-2,3*j-4)=3*pow(u,2);
		Aspl(3*j-2,3*j-3)=2*u;
		Aspl(3*j-2,3*j-2)=1;
		Aspl(3*j-2,3*j+1)=-1;
	 
		Aspl(3*j-1,3*j-4)=6*u;
		Aspl(3*j-1,3*j-3)=2;
		Aspl(3*j-1,3*j)=-2;
	}
	
	j=ind_xlims(0)+1;
	Dg(3*j-3)=F(j-1)-F(j-2);
	u=1/(x(j-1)-x0g)-1/(x(j-2)-x0g);
	double ug=1/(x(j-1)-x0g);
	
	Aspl(3*j-3,3*j-4)=pow(u,3);
	Aspl(3*j-3,3*j-3)=pow(u,2);
	Aspl(3*j-3,3*j-2)=u;
	
	Aspl(3*j-2,3*j-4)=-3*pow(ug,2)*pow(u,2);
	Aspl(3*j-2,3*j-3)=-2*pow(ug,2)*u;
	Aspl(3*j-2,3*j-2)=-pow(ug,2);
	
	Dg(3*j-2)=(F(j)-F(j-2))/(x(j)-x(j-2));
	
	vec Cg=solve(Aspl,Dg);
/*
	mat T(NSg,Nx);
	T.zeros();
	T(0,0)=1;
	for (j=2; j<=Ng; j++)
	{
	 T(3*j-3,j-2)=-1;
	 T(3*j-3,j-1)=1;
	}
	T(3*Ng-2,Ng-2)=-1/(x(Ng)-x(Ng-2));
	T(3*Ng-2,Ng)=1/(x(Ng)-x(Ng-2));
	
	Cg=Cg*T;
*/
	
	//solve the spline in the interval [wl,wr]
	int NSc=3*Nc-1;
	Aspl.zeros(NSc,NSc);
	
	double Dx=x(Ng)-x(Ng-1);
	
	vec Dc=zeros<vec>(NSc);
	Dc(0)=F(Ng)-Dx*(F(Ng)-F(Ng-2))/(x(Ng)-x(Ng-2))-F(Ng-1);
	
	Aspl(0,0)=pow(Dx,3);
	Aspl(0,1)=pow(Dx,2);
	
	Dc(1)=-(F(Ng)-F(Ng-2))/(x(Ng)-x(Ng-2));
	
	Aspl(1,0)=3*pow(Dx,2);
	Aspl(1,1)=2*Dx;
	Aspl(1,4)=-1;
	
	Aspl(2,0)=6*Dx;
	Aspl(2,1)=2;
	Aspl(2,3)=-2;
	
	for (j=2; j<Nc; j++)
	{
		Dx=x(Ng+j-1)-x(Ng+j-2);
	 
		Dc(3*j-3)=F(Ng+j-1)-F(Ng+j-2);
	 
		Aspl(3*j-3,3*j-4)=pow(Dx,3);
		Aspl(3*j-3,3*j-3)=pow(Dx,2);
		Aspl(3*j-3,3*j-2)=Dx;
	 
		Aspl(3*j-2,3*j-4)=3*pow(Dx,2);
		Aspl(3*j-2,3*j-3)=2*Dx;
		Aspl(3*j-2,3*j-2)=1;
		Aspl(3*j-2,3*j+1)=-1;
	 
		Aspl(3*j-1,3*j-4)=6*Dx;
		Aspl(3*j-1,3*j-3)=2;
		Aspl(3*j-1,3*j)=-2;
	}
	
	j=Nc;
	Dx=x(Ng+j-1)-x(Ng+j-2);
	Dc(3*j-3)=F(Ng+j-1)-F(Ng+j-2);
	Aspl(3*j-3,3*j-4)=pow(Dx,3);
	Aspl(3*j-3,3*j-3)=pow(Dx,2);
	Aspl(3*j-3,3*j-2)=Dx;
	
	Dc(3*j-2)=(F(Ng+j)-F(Ng+j-2))/(x(Ng+j)-x(Ng+j-2));
	Aspl(3*j-2,3*j-4)=3*pow(Dx,2);
	Aspl(3*j-2,3*j-3)=2*Dx;
	Aspl(3*j-2,3*j-2)=1;
	
	vec Cc=solve(Aspl,Dc);
	
/*
	T.zeros(NSc,Nx);
	T(0,Ng-2)=(x(Ng)-x(Ng-1))/(x(Ng)-x(Ng-2));
	T(0,Ng-1)=-1;
	T(0,Ng)=1-(x(Ng)-x(Ng-1))/(x(Ng)-x(Ng-2));
	T(1,Ng-2)=1/(x(Ng)-x(Ng-2));
	T(1,Ng)=-1/(x(Ng)-x(Ng-2));
	for (j=2; j<=Nc; j++)
	{
		T(3*j-3,Ng+j-2)=-1;
		T(3*j-3,Ng+j-1)=1;
	}
	T(3*Nc-2,ind_xlims(1)-1)=-1/(x(ind_xlims(1)+1)-x(ind_xlims(1)-1));
	T(3*Nc-2,ind_xlims(1)+1)=1/(x(ind_xlims(1)+1)-x(ind_xlims(1)-1));
	
	Cc=Cc*T;
	
	mat Cc2(NSc+1,Nx);
	Cc2.zeros();
	
	Cc2.rows(0,1)=Cc.rows(0,1);
	Cc2.rows(3,NSc)=Cc.rows(2,NSc-1);
	Cc2(2,Ng-2)=-1/(x(Ng)-x(Ng-2));
	Cc2(2,Ng)=1/(x(Ng)-x(Ng-2));
*/
	
	//solve the spline in the interval [wr,inf
	int NSd=3*Nd-1;
	Aspl.zeros(NSd,NSd);
	vec Dd=zeros<vec>(NSd);
	
	u=1/(x(ind_xlims(1))-x0d)-1/(x(ind_xlims(1)+1)-x0d);
	double ud=1/(x(ind_xlims(1))-x0d);
	
	Dd(0)=(F(ind_xlims(1)+1)-F(ind_xlims(1)-1))/(x(ind_xlims(1)+1)-x(ind_xlims(1)-1));
	Aspl(0,0)=-3*pow(ud,2)*pow(u,2);
	Aspl(0,1)=-2*pow(ud,2)*u;
	Aspl(0,2)=-pow(ud,2);
	
	Dd(1)=F(ind_xlims(1))-F(ind_xlims(1)+1);
	Aspl(1,0)=pow(u,3);
	Aspl(1,1)=pow(u,2);
	Aspl(1,2)=u;
	
	for (j=1; j<Nd-1; j++)
	{
		u=1/(x(ind_xlims(1)+j)-x0d)-1/(x(ind_xlims(1)+j+1)-x0d);
		
		Aspl(3*j-1,3*j-1)=-1;
		Aspl(3*j-1,3*j)=3*pow(u,2);
		Aspl(3*j-1,3*j+1)=2*u;
		Aspl(3*j-1,3*j+2)=1;
		
		Aspl(3*j,3*j-2)=-2;
		Aspl(3*j,3*j)=6*u;
		Aspl(3*j,3*j+1)=2;
	 
		Dd(3*j+1)=F(ind_xlims(1)+j)-F(ind_xlims(1)+j+1);
		Aspl(3*j+1,3*j)=pow(u,3);
		Aspl(3*j+1,3*j+1)=pow(u,2);
		Aspl(3*j+1,3*j+2)=u;
	}
	
	j=Nd-1;
	u=1/(x(ind_xlims(1)+j)-x0d);
	
	Aspl(3*j-1,3*j-1)=-1;
	Aspl(3*j-1,3*j)=3*pow(u,2);
	Aspl(3*j-1,3*j+1)=2*u;
	
	Aspl(3*j,3*j-2)=-2;
	Aspl(3*j,3*j)=6*u;
	Aspl(3*j,3*j+1)=2;
	
	Dd(3*j+1)=F(ind_xlims(1)+j);
	Aspl(3*j+1,3*j)=pow(u,3);
	Aspl(3*j+1,3*j+1)=pow(u,2);
	
	vec Cd=solve(Aspl,Dd);
	
	int Nint_tot=Nx+1;
	int NS_tot=4*Nint_tot;
	coeffs.zeros(NS_tot);
	
	coeffs(0)=Cg(0);
	coeffs(1)=Cg(1);
	
	uvec ind_tmp=linspace<uvec>(1,Ng-1,Ng-1);
	coeffs.rows(4*ind_tmp)=Cg(3*ind_tmp-1);
	coeffs.rows(4*ind_tmp+1)=Cg(3*ind_tmp);
	coeffs.rows(4*ind_tmp+2)=Cg(3*ind_tmp+1);
	coeffs.rows(4*ind_tmp+3)=F.rows(0,ind_xlims(0)-1);
	
	coeffs(4*Ng)=Cc(0);
	coeffs(4*Ng+1)=Cc(1);
	coeffs(4*Ng+2)=(F(Ng)-F(Ng-2))/(x(Ng)-x(Ng-2));
	coeffs(4*Ng+3)=F(ind_xlims(0));
	
	ind_tmp=linspace<uvec>(1,Nc-1,Nc-1);
	coeffs.rows(4*ind_tmp+4*Ng)=Cc(3*ind_tmp-1);
	coeffs.rows(4*ind_tmp+4*Ng+1)=Cc(3*ind_tmp);
	coeffs.rows(4*ind_tmp+4*Ng+2)=Cc(3*ind_tmp+1);
	coeffs.rows(4*ind_tmp+4*Ng+3)=F.rows(ind_xlims(0)+1,ind_xlims(1)-1);
	
	ind_tmp=linspace<uvec>(0,Nd-2,Nd-1);
	coeffs.rows(4*ind_tmp+4*(Ng+Nc))=Cd(3*ind_tmp);
	coeffs.rows(4*ind_tmp+4*(Ng+Nc)+1)=Cd(3*ind_tmp+1);
	coeffs.rows(4*ind_tmp+4*(Ng+Nc)+2)=Cd(3*ind_tmp+2);
	coeffs.rows(4*ind_tmp+4*(Ng+Nc)+3)=F.rows(ind_xlims(1)+1,Nx-1);
	
	coeffs(NS_tot-4)=Cd(3*Nd-3);
	coeffs(NS_tot-3)=Cd(3*Nd-2);
	
/*
	T.zeros(NSd,Nx);
	T(0,ind_xlims(1)-1)=-1/(x(ind_xlims(1)+1)-x(ind_xlims(1)-1));
	T(0,ind_xlims(1)+1)=1/(x(ind_xlims(1)+1)-x(ind_xlims(1)-1));
	for (j=1; j<Nd; j++)
	{
		T(3*j-2,ind_xlims(1)+j-1)=1;
		T(3*j-2,ind_xlims(1)+j)=-1;
	}
	j=Nd;
	T(3*j-2,ind_xlims(1)+j-1)=1;
	
	Cd=Cd*T;
	
	int NS_tot=NSg+NSc+NSd+1;
	
	Mspl.zeros(NS_tot,Nx);
	Mspl.rows(0,NSg-1)=Cg;
	Mspl.rows(NSg,NSg+NSc)=Cc2;
	Mspl.rows(NSg+NSc+1,NS_tot-1)=Cd;
*/
	return true;
}

bool OmegaMaxEnt_data::spline_val_grid_transf(vec x, vec x0, vec coeffs, vec &s)
{
	int Nx=x.n_rows;
	int Nx0=x0.n_rows;
	int Ns=Nx0-1;
	int Nc=4*Ns;
	
	vec D=1.0/(x0.rows(1,Nx0-1)-x0.rows(0,Nx0-2));
	
	if (coeffs.n_rows<Nc)
	{
		cout<<"spline_val_grid_transf(): either the coefficients vector is incomplete or the position vector is too long\n";
		return false;
	}
	
	s.zeros(Nx);
	
	double a, b, c, d, Dx;
	
	int j, l;
	for (j=0; j<Nx; j++)
	{
		if (x(j)>=x0(0) && x(j)<x0(Nx0-1))
		{
			l=0;
			while ( x(j)>=x0(l+1) && l<Ns-1) l++;
			a=coeffs(4*l);
			b=coeffs(4*l+1);
			c=coeffs(4*l+2);
			d=coeffs(4*l+3);
			
			Dx=D(l)*(x(j)-x0(l));
			s(j)=a*pow(Dx,3)+b*pow(Dx,2)+c*Dx+d;
		}
		else if (x(j)==x0(Nx0-1))
		{
			if (coeffs.n_rows>Nc)
			{
				s(j)=coeffs[Nc];
			}
			else
			{
				l=Ns-1;
				a=coeffs(4*l);
				b=coeffs(4*l+1);
				c=coeffs(4*l+2);
				d=coeffs(4*l+3);
				
				s(j)=a+b+c+d;
			}
		}
		else
		{
			cout<<"spline_val_grid_transf(): the provided position vector has a value outside the boundaries of the spline.\n";
			return false;
		}
	}
	
	return true;
}

bool OmegaMaxEnt_data::spline_matrix_grid_transf(vec w0, mat &M)
{
	int N=w0.n_rows;
	vec D=1.0/(w0.rows(1,N-1)-w0.rows(0,N-2));
	
	int Nint=N-1;
	
	mat B=zeros<mat>(3*Nint-1,3*Nint-1);
	//vec vTl(3*Nint-1), vTr(3*Nint-1);
	mat Ps=zeros<mat>(3*Nint-1, N);
	mat Pg=zeros<mat>(4*Nint,4*Nint-1);
	
	//double N1p5=pow(1.0*Nint,1.5), N0p5=pow(1.0*Nint,0.5), Nm0p5=pow(1.0*Nint,-0.5);
	
	B(0,0)=1;
	B(0,1)=1;
	B(1,0)=3;
	B(1,1)=2;
	B(1,4)=-D(1)/D(0);
	B(2,0)=6;
	B(2,1)=2;
	B(2,3)=-2*pow(D(1)/D(0),2);
	
//	vTl(0)=N1p5;
//	vTl(1)=N0p5;
//	vTl(2)=Nm0p5;
	
//	vTr(0)=N1p5;
//	vTr(1)=N0p5;
	
	Ps(0,0)=-1;
	Ps(0,1)=1;
	
	Pg(0,0)=1;
	Pg(1,1)=1;
	Pg(3,3*Nint-1)=1;
	
	int j;
	for (j=1; j<Nint-1; j++)
	{
		B(3*j,3*j-1)=1;
		B(3*j,3*j)=1;
		B(3*j,3*j+1)=1;
		B(3*j+1,3*j-1)=3;
		B(3*j+1,3*j)=2;
		B(3*j+1,3*j+1)=1;
		B(3*j+1,3*j+4)=-D(j+1)/D(j);
		B(3*j+2,3*j-1)=6;
		B(3*j+2,3*j)=2;
		B(3*j+2,3*j+3)=-2*pow(D(j+1)/D(j),2);
		
//		vTl(3*j)=N1p5;
//		vTl(3*j+1)=N0p5;
//		vTl(3*j+2)=Nm0p5;
		
//		vTr(3*j-1)=N1p5;
//		vTr(3*j)=N0p5;
//		vTr(3*j+1)=Nm0p5;
		
		Ps(3*j,j)=-1;
		Ps(3*j,j+1)=1;
		
		Pg(4*j,3*j-1)=1;
		Pg(4*j+1,3*j)=1;
		Pg(4*j+2,3*j+1)=1;
		Pg(4*j+3,3*Nint-1+j)=1;
	}
	
	j=Nint-1;
	B(3*j,3*j-1)=1;
	B(3*j,3*j)=1;
	B(3*j,3*j+1)=1;
	B(3*j+1,3*j-1)=3;
	B(3*j+1,3*j)=2;
	B(3*j+1,3*j+1)=1;
	
//	vTl(3*j)=N1p5;
//	vTl(3*j+1)=N0p5;
	
//	vTr(3*j-1)=N1p5;
//	vTr(3*j)=N0p5;
//	vTr(3*j+1)=Nm0p5;
	
	Ps(3*j,j)=-1;
	Ps(3*j,j+1)=1;
	
	Pg(4*j,3*j-1)=1;
	Pg(4*j+1,3*j)=1;
	Pg(4*j+2,3*j+1)=1;
	Pg(4*j+3,3*Nint-1+j)=1;
	
//	mat Tl=diagmat(vTl);
//	mat Tr=diagmat(vTr);
	
	mat IB=eye(3*Nint-1,3*Nint-1);
	mat invB=solve(B,IB);
	
	mat IA=eye(N,N);
	mat PA=IA.submat(0,0,N-2,N-1);
	
	mat Lg=join_vert(invB*Ps,PA);
	
//	mat Lg=join_vert(Tr*invB*Tl*Ps,PA);
	
	M=Pg*Lg;
	
	return true;
}

bool OmegaMaxEnt_data::spline_matrix_grid_transf_G_part(vec x, uvec ind_xlims, vec xs, mat &M)
{
	int Nx=x.n_rows;
	
	int Ng=ind_xlims(0);
	int Nc=ind_xlims(1)-ind_xlims(0);
	int Nd=Nx-ind_xlims(1)-1;
	
	double x0g=xs(0);
	double x0d=xs(1);
	
	//	 solve the spline in the interval -inf,wl]
	//vec Dg=-((x.rows(1,ind_xlims(0))-x0g)%(x.rows(0,ind_xlims(0)-1)-x0g))/(x.rows(1,ind_xlims(0))-x.rows(0,ind_xlims(0)-1));
	
	vec RDg=((x.rows(2,ind_xlims(0))-x0g)/(x.rows(0,ind_xlims(0)-2)-x0g))%((x.rows(1,ind_xlims(0)-1)-x.rows(0,ind_xlims(0)-2))/(x.rows(2,ind_xlims(0))-x.rows(1,ind_xlims(0)-1)));
	
	int NCg=3*Ng-1;
	mat B=zeros<mat>(NCg,NCg);
	mat Ps=zeros<mat>(NCg, Nx);
	mat Pg=zeros<mat>(4*Ng,4*Ng-1);
	
	B(0,0)=1;
	B(0,1)=1;
	B(1,0)=3;
	B(1,1)=2;
	B(1,4)=-RDg(0);
	B(2,0)=6;
	B(2,1)=2;
	B(2,3)=-2*pow(RDg(0),2);
	
	Ps(0,0)=-1;
	Ps(0,1)=1;
	
	Pg(0,0)=1;
	Pg(1,1)=1;
	Pg(3,NCg)=1;
	
	int j;
	for (j=1; j<Ng-1; j++)
	{
		B(3*j,3*j-1)=1;
		B(3*j,3*j)=1;
		B(3*j,3*j+1)=1;
		B(3*j+1,3*j-1)=3;
		B(3*j+1,3*j)=2;
		B(3*j+1,3*j+1)=1;
		B(3*j+1,3*j+4)=-RDg(j);
		B(3*j+2,3*j-1)=6;
		B(3*j+2,3*j)=2;
		B(3*j+2,3*j+3)=-2*pow(RDg(j),2);
		
		Ps(3*j,j)=-1;
		Ps(3*j,j+1)=1;
		
		Pg(4*j,3*j-1)=1;
		Pg(4*j+1,3*j)=1;
		Pg(4*j+2,3*j+1)=1;
		Pg(4*j+3,NCg+j)=1;
	}
	
	j=Ng-1;
	B(3*j,3*j-1)=1;
	B(3*j,3*j)=1;
	B(3*j,3*j+1)=1;
	B(3*j+1,3*j-1)=3;
	B(3*j+1,3*j)=2;
	B(3*j+1,3*j+1)=1;
	
	double fdAg=((x(j+1)-x0g)/(x(j)-x0g))*((x(j+1)-x(j))/(x(j+2)-x(j)));
	
	Ps(3*j,j)=-1;
	Ps(3*j,j+1)=1;
	Ps(3*j+1,j)=-fdAg;
	Ps(3*j+1,j+2)=fdAg;
	
	Pg(4*j,3*j-1)=1;
	Pg(4*j+1,3*j)=1;
	Pg(4*j+2,3*j+1)=1;
	Pg(4*j+3,NCg+j)=1;
	
	mat IB=eye(NCg,NCg);
	mat invB=solve(B,IB);
	
	mat IA=eye(Nx,Nx);
	mat PA=IA.submat(0,0,Ng-1,Nx-1);
	
	mat Lg=join_vert(invB*Ps,PA);
	
	mat Mg=Pg*Lg;
	
	// solve in the interval [wl,wr]
	//vec Dc=1.0/(x.rows(ind_xlims(0)+1,ind_xlims(1))-x.rows(ind_xlims(0),ind_xlims(1)-1));
	
	vec RDc=(x.rows(ind_xlims(0)+1,ind_xlims(1)-1)-x.rows(ind_xlims(0),ind_xlims(1)-2))/(x.rows(ind_xlims(0)+2,ind_xlims(1))-x.rows(ind_xlims(0)+1,ind_xlims(1)-1));
	
	int NCc=3*Nc-1;
	
	B.zeros(NCc,NCc);
	Ps.zeros(NCc, Nx);
	Pg.zeros(4*Nc,4*Nc);
	
	B(0,0)=1;
	B(0,1)=1;
	B(1,0)=3;
	B(1,1)=2;
	B(1,4)=-RDc(0);
	B(2,0)=6;
	B(2,1)=2;
	B(2,3)=-2*pow(RDc(0),2);
	
	Ps(0,Ng-1)=(x(Ng+1)-x(Ng))/(x(Ng+1)-x(Ng-1));
	Ps(0,Ng)=-1;
	Ps(0,Ng+1)=1-(x(Ng+1)-x(Ng))/(x(Ng+1)-x(Ng-1));
	Ps(1,Ng-1)=(x(Ng+1)-x(Ng))/(x(Ng+1)-x(Ng-1));
	Ps(1,Ng+1)=-(x(Ng+1)-x(Ng))/(x(Ng+1)-x(Ng-1));
	
	Pg(0,0)=1;
	Pg(1,1)=1;
	Pg(2,NCc)=-(x(Ng+1)-x(Ng))/(x(Ng+1)-x(Ng-1));
	Pg(2,NCc+2)=(x(Ng+1)-x(Ng))/(x(Ng+1)-x(Ng-1));
	Pg(3,NCc+1)=1;
	
	for (j=1; j<Nc-1; j++)
	{
		B(3*j,3*j-1)=1;
		B(3*j,3*j)=1;
		B(3*j,3*j+1)=1;
		B(3*j+1,3*j-1)=3;
		B(3*j+1,3*j)=2;
		B(3*j+1,3*j+1)=1;
		B(3*j+1,3*j+4)=-RDc(j);
		B(3*j+2,3*j-1)=6;
		B(3*j+2,3*j)=2;
		B(3*j+2,3*j+3)=-2*pow(RDc(j),2);
		
		Ps(3*j,Ng+j)=-1;
		Ps(3*j,Ng+j+1)=1;
		
		Pg(4*j,3*j-1)=1;
		Pg(4*j+1,3*j)=1;
		Pg(4*j+2,3*j+1)=1;
		Pg(4*j+3,NCc+j+1)=1;
	}
	
	j=Nc-1;
	B(3*j,3*j-1)=1;
	B(3*j,3*j)=1;
	B(3*j,3*j+1)=1;
	B(3*j+1,3*j-1)=3;
	B(3*j+1,3*j)=2;
	B(3*j+1,3*j+1)=1;
	
	Ps(3*j,Ng+j)=-1;
	Ps(3*j,Ng+j+1)=1;
	Ps(3*j+1,Ng+j)=-(x(ind_xlims(1))-x(ind_xlims(1)-1))/(x(ind_xlims(1)+1)-x(ind_xlims(1)-1));
	Ps(3*j+1,Ng+j+2)=(x(ind_xlims(1))-x(ind_xlims(1)-1))/(x(ind_xlims(1)+1)-x(ind_xlims(1)-1));
	
	Pg(4*j,3*j-1)=1;
	Pg(4*j+1,3*j)=1;
	Pg(4*j+2,3*j+1)=1;
	Pg(4*j+3,NCc+j+1)=1;
	
	IB.eye(NCc,NCc);
	invB=solve(B,IB);
	
	PA=IA.submat(Ng-1,0,ind_xlims(1)-1,Nx-1);
	
	Lg=join_vert(invB*Ps,PA);
	
	mat Mc=Pg*Lg;
	
	int NCd=3*Nd-1;
	
	//vec RDd=((x.rows(ind_xlims(1)+2,Nx-1)-x.rows(ind_xlims(1)+1,Nx-2))/(x.rows(ind_xlims(1)+1,Nx-2)-x.rows(ind_xlims(1),Nx-3)))%((x.rows(ind_xlims(1),Nx-3)-x0d)/(x.rows(ind_xlims(1)+2,Nx-1)-x0d));
	
	double fdAd=((x(ind_xlims(1)+1)-x(ind_xlims(1)))/(x(ind_xlims(1)+1)-x(ind_xlims(1)-1)))*((x(ind_xlims(1))-x0d)/(x(ind_xlims(1)+1)-x0d));
	
	vec RDd=((x.rows(ind_xlims(1)+2,Nx-1)-x0d)/(x.rows(ind_xlims(1),Nx-3)-x0d))%((x.rows(ind_xlims(1),Nx-3)-x.rows(ind_xlims(1)+1,Nx-2))/(x.rows(ind_xlims(1)+1,Nx-2)-x.rows(ind_xlims(1)+2,Nx-1)));
	
	B.zeros(NCd,NCd);
	Ps.zeros(NCd, Nx);
	Pg.zeros(4*Nd,4*Nd);
	
	B(0,0)=1;
	B(0,1)=1;
	B(1,0)=3;
	B(1,1)=2;
	B(1,4)=-RDd(0);
	B(2,0)=6;
	B(2,1)=2;
	B(2,3)=-2*pow(RDd(0),2);
	
	Ps(0,ind_xlims(1)-1)=fdAd;
	Ps(0,ind_xlims(1))=-1;
	Ps(0,ind_xlims(1)+1)=1-fdAd;
	Ps(1,ind_xlims(1)-1)=fdAd;
	Ps(1,ind_xlims(1)+1)=-fdAd;
	
	Pg(0,0)=1;
	Pg(1,1)=1;
	Pg(2,NCd)=-fdAd;
	Pg(2,NCd+2)=fdAd;
	Pg(3,NCd+1)=1;
	
	for (j=1; j<Nd-1; j++)
	{
		B(3*j,3*j-1)=1;
		B(3*j,3*j)=1;
		B(3*j,3*j+1)=1;
		B(3*j+1,3*j-1)=3;
		B(3*j+1,3*j)=2;
		B(3*j+1,3*j+1)=1;
		B(3*j+1,3*j+4)=-RDd(j);
		B(3*j+2,3*j-1)=6;
		B(3*j+2,3*j)=2;
		B(3*j+2,3*j+3)=-2*pow(RDd(j),2);
		
		Ps(3*j,ind_xlims(1)+j)=-1;
		Ps(3*j,ind_xlims(1)+j+1)=1;
		
		Pg(4*j,3*j-1)=1;
		Pg(4*j+1,3*j)=1;
		Pg(4*j+2,3*j+1)=1;
		Pg(4*j+3,NCd+j+1)=1;
	}
	
	j=Nd-1;
	B(3*j,3*j-1)=1;
	B(3*j,3*j)=1;
	B(3*j,3*j+1)=1;
	B(3*j+1,3*j-1)=3;
	B(3*j+1,3*j)=2;
	B(3*j+1,3*j+1)=1;
	
	Ps(3*j,ind_xlims(1)+j)=-1;
	Ps(3*j,ind_xlims(1)+j+1)=1;
	
	Pg(4*j,3*j-1)=1;
	Pg(4*j+1,3*j)=1;
	Pg(4*j+2,3*j+1)=1;
	Pg(4*j+3,NCd+j+1)=1;
	
	IB.eye(NCd,NCd);
	invB=solve(B,IB);
	
	PA=IA.submat(ind_xlims(1)-1,0,Nx-2,Nx-1);
	
	Lg=join_vert(invB*Ps,PA);
	
	mat Md=Pg*Lg;
	
	/*
	Pg.zeros(4*Nd,4*Nd-1);
	 
	B(0,0)=3;
	B(0,1)=2;
	B(0,2)=1;
	B(1,0)=1;
	B(1,1)=1;
	B(1,2)=1;
	
	Ps(0,ind_xlims(1)-1)=-fdAd;
	Ps(0,ind_xlims(1)+1)=fdAd;
	Ps(1,ind_xlims(1))=1;
	Ps(1,ind_xlims(1)+1)=-1;
	
	Pg(0,0)=1;
	Pg(1,1)=1;
	Pg(2,2)=1;
	Pg(3,NCd)=1;
	
	for (j=1; j<Nd-1; j++)
	{
		B(3*j-1,3*j-1)=-RDd(j-1);
		B(3*j-1,3*j)=3;
		B(3*j-1,3*j+1)=2;
		B(3*j-1,3*j+2)=1;
		
		B(3*j,3*j-2)=-2*pow(RDd(j-1),2);
		B(3*j,3*j)=6;
		B(3*j,3*j+1)=2;
		
		B(3*j+1,3*j)=1;
		B(3*j+1,3*j+1)=1;
		B(3*j+1,3*j+2)=1;
		
		Ps(3*j+1,ind_xlims(1)+j)=1;
		Ps(3*j+1,ind_xlims(1)+j+1)=-1;
		
		Pg(4*j,3*j)=1;
		Pg(4*j+1,3*j+1)=1;
		Pg(4*j+2,3*j+2)=1;
		Pg(4*j+3,NCd+j)=1;
	}
	
	j=Nd-1;
	B(3*j-1,3*j-1)=-RDd(j-1);
	B(3*j-1,3*j)=3;
	B(3*j-1,3*j+1)=2;
	
	B(3*j,3*j-2)=-2*pow(RDd(j-1),2);
	B(3*j,3*j)=6;
	B(3*j,3*j+1)=2;
	
	B(3*j+1,3*j)=1;
	B(3*j+1,3*j+1)=1;
	
	Ps(3*j+1,ind_xlims(1)+j)=1;
	Ps(3*j+1,ind_xlims(1)+j+1)=-1;
	
	Pg(4*j,3*j)=1;
	Pg(4*j+1,3*j+1)=1;
	Pg(4*j+3,NCd+j)=1;
	 
	IB.eye(NCd,NCd);
	invB=solve(B,IB);
	
	PA=IA.submat(ind_xlims(1)+1,0,Nx-1,Nx-1);
	
	Lg=join_vert(invB*Ps,PA);
	
	mat Md=Pg*Lg;
	*/
	 
	M=join_vert(Mg,Mc);
	M=join_vert(M,Md);
	
	return true;
}

bool OmegaMaxEnt_data::spline_matrix_grid_transf_G_part_1(vec x, uvec ind_xlims, vec xs, mat &M)
{
	int Nx=x.n_rows;
	
	int Ng=ind_xlims(0);
	int Nc=ind_xlims(1)-ind_xlims(0);
	int Nd=Nx-ind_xlims(1)-1;
	
	double x0g=xs(0);
	double x0d=xs(1);
	
	//	 solve the spline in the interval -inf,wl]
	//vec Dg=-((x.rows(1,ind_xlims(0))-x0g)%(x.rows(0,ind_xlims(0)-1)-x0g))/(x.rows(1,ind_xlims(0))-x.rows(0,ind_xlims(0)-1));
	
	vec RDg=((x.rows(2,ind_xlims(0))-x0g)/(x.rows(0,ind_xlims(0)-2)-x0g))%((x.rows(1,ind_xlims(0)-1)-x.rows(0,ind_xlims(0)-2))/(x.rows(2,ind_xlims(0))-x.rows(1,ind_xlims(0)-1)));
	
	int NCg=3*Ng-1;
	mat B=zeros<mat>(NCg,NCg);
	mat Ps=zeros<mat>(NCg, Nx);
	mat Pg=zeros<mat>(4*Ng,4*Ng-1);
	
	B(0,0)=1;
	B(0,1)=1;
	B(1,0)=3;
	B(1,1)=2;
	B(1,4)=-RDg(0);
	B(2,0)=6;
	B(2,1)=2;
	B(2,3)=-2*pow(RDg(0),2);
	
	Ps(0,0)=-1;
	Ps(0,1)=1;
	
	Pg(0,0)=1;
	Pg(1,1)=1;
	Pg(3,NCg)=1;
	
	int j;
	for (j=1; j<Ng-1; j++)
	{
		B(3*j,3*j-1)=1;
		B(3*j,3*j)=1;
		B(3*j,3*j+1)=1;
		B(3*j+1,3*j-1)=3;
		B(3*j+1,3*j)=2;
		B(3*j+1,3*j+1)=1;
		B(3*j+1,3*j+4)=-RDg(j);
		B(3*j+2,3*j-1)=6;
		B(3*j+2,3*j)=2;
		B(3*j+2,3*j+3)=-2*pow(RDg(j),2);
		
		Ps(3*j,j)=-1;
		Ps(3*j,j+1)=1;
		
		Pg(4*j,3*j-1)=1;
		Pg(4*j+1,3*j)=1;
		Pg(4*j+2,3*j+1)=1;
		Pg(4*j+3,NCg+j)=1;
	}
	
	j=Ng-1;
	B(3*j,3*j-1)=1;
	B(3*j,3*j)=1;
	B(3*j,3*j+1)=1;
	B(3*j+1,3*j-1)=3;
	B(3*j+1,3*j)=2;
	B(3*j+1,3*j+1)=1;
	
	double fdAg=((x(j+1)-x0g)/(x(j)-x0g))*((x(j+1)-x(j))/(x(j+2)-x(j)));
	
	Ps(3*j,j)=-1;
	Ps(3*j,j+1)=1;
	Ps(3*j+1,j)=-fdAg;
	Ps(3*j+1,j+2)=fdAg;
	
	Pg(4*j,3*j-1)=1;
	Pg(4*j+1,3*j)=1;
	Pg(4*j+2,3*j+1)=1;
	Pg(4*j+3,NCg+j)=1;
	
	mat IB=eye(NCg,NCg);
	mat invB=solve(B,IB);
	
	mat IA=eye(Nx,Nx);
	mat PA=IA.submat(0,0,Ng-1,Nx-1);
	
	mat Lg=join_vert(invB*Ps,PA);
	
	mat Mg=Pg*Lg;
	
	// solve in the interval [wl,wr]
	//vec Dc=1.0/(x.rows(ind_xlims(0)+1,ind_xlims(1))-x.rows(ind_xlims(0),ind_xlims(1)-1));
	
	vec RDc=(x.rows(ind_xlims(0)+1,ind_xlims(1)-1)-x.rows(ind_xlims(0),ind_xlims(1)-2))/(x.rows(ind_xlims(0)+2,ind_xlims(1))-x.rows(ind_xlims(0)+1,ind_xlims(1)-1));
	
	int NCc=3*Nc-1;
	
	B.zeros(NCc,NCc);
	Ps.zeros(NCc, Nx);
	Pg.zeros(4*Nc,4*Nc);
	
	B(0,0)=1;
	B(0,1)=1;
	B(1,0)=3;
	B(1,1)=2;
	B(1,4)=-RDc(0);
	B(2,0)=6;
	B(2,1)=2;
	B(2,3)=-2*pow(RDc(0),2);
	
	Ps(0,Ng-1)=(x(Ng+1)-x(Ng))/(x(Ng+1)-x(Ng-1));
	Ps(0,Ng)=-1;
	Ps(0,Ng+1)=1-(x(Ng+1)-x(Ng))/(x(Ng+1)-x(Ng-1));
	Ps(1,Ng-1)=(x(Ng+1)-x(Ng))/(x(Ng+1)-x(Ng-1));
	Ps(1,Ng+1)=-(x(Ng+1)-x(Ng))/(x(Ng+1)-x(Ng-1));
	
	Pg(0,0)=1;
	Pg(1,1)=1;
	Pg(2,NCc)=-(x(Ng+1)-x(Ng))/(x(Ng+1)-x(Ng-1));
	Pg(2,NCc+2)=(x(Ng+1)-x(Ng))/(x(Ng+1)-x(Ng-1));
	Pg(3,NCc+1)=1;
	
	for (j=1; j<Nc-1; j++)
	{
		B(3*j,3*j-1)=1;
		B(3*j,3*j)=1;
		B(3*j,3*j+1)=1;
		B(3*j+1,3*j-1)=3;
		B(3*j+1,3*j)=2;
		B(3*j+1,3*j+1)=1;
		B(3*j+1,3*j+4)=-RDc(j);
		B(3*j+2,3*j-1)=6;
		B(3*j+2,3*j)=2;
		B(3*j+2,3*j+3)=-2*pow(RDc(j),2);
		
		Ps(3*j,Ng+j)=-1;
		Ps(3*j,Ng+j+1)=1;
		
		Pg(4*j,3*j-1)=1;
		Pg(4*j+1,3*j)=1;
		Pg(4*j+2,3*j+1)=1;
		Pg(4*j+3,NCc+j+1)=1;
	}
	
	j=Nc-1;
	B(3*j,3*j-1)=1;
	B(3*j,3*j)=1;
	B(3*j,3*j+1)=1;
	B(3*j+1,3*j-1)=3;
	B(3*j+1,3*j)=2;
	B(3*j+1,3*j+1)=1;
	
	Ps(3*j,Ng+j)=-1;
	Ps(3*j,Ng+j+1)=1;
	Ps(3*j+1,Ng+j)=-(x(ind_xlims(1))-x(ind_xlims(1)-1))/(x(ind_xlims(1)+1)-x(ind_xlims(1)-1));
	Ps(3*j+1,Ng+j+2)=(x(ind_xlims(1))-x(ind_xlims(1)-1))/(x(ind_xlims(1)+1)-x(ind_xlims(1)-1));
	
	Pg(4*j,3*j-1)=1;
	Pg(4*j+1,3*j)=1;
	Pg(4*j+2,3*j+1)=1;
	Pg(4*j+3,NCc+j+1)=1;
	
	IB.eye(NCc,NCc);
	invB=solve(B,IB);
	
	PA=IA.submat(Ng-1,0,ind_xlims(1)-1,Nx-1);
	
	Lg=join_vert(invB*Ps,PA);
	
	mat Mc=Pg*Lg;
	
	int NCd=3*Nd-1;
	
	//vec Dd=((x.rows(ind_xlims(1),Nx-2)-x0d)%(x.rows(ind_xlims(1)+1,Nx-1)-x0d))/(x.rows(ind_xlims(1)+1,Nx-1)-x.rows(ind_xlims(1),Nx-2));
	
	vec RDd=((x.rows(ind_xlims(1)+2,Nx-1)-x.rows(ind_xlims(1)+1,Nx-2))/(x.rows(ind_xlims(1)+1,Nx-2)-x.rows(ind_xlims(1),Nx-3)))%((x.rows(ind_xlims(1),Nx-3)-x0d)/(x.rows(ind_xlims(1)+2,Nx-1)-x0d));
	
	double fdAd=((x(ind_xlims(1))-x0d)/(x(ind_xlims(1)+1)-x0d))*((x(ind_xlims(1))-x(ind_xlims(1)+1))/(x(ind_xlims(1)+1)-x(ind_xlims(1)-1)));
	
	B.zeros(NCd,NCd);
	Ps.zeros(NCd, Nx);
	Pg.zeros(4*Nd,4*Nd-1);
	
	B(0,0)=3;
	B(0,1)=2;
	B(0,2)=1;
	B(1,0)=1;
	B(1,1)=1;
	B(1,2)=1;
	
	Ps(0,ind_xlims(1)-1)=-fdAd;
	Ps(0,ind_xlims(1)+1)=fdAd;
	Ps(1,ind_xlims(1))=1;
	Ps(1,ind_xlims(1)+1)=-1;
	
	Pg(0,0)=1;
	Pg(1,1)=1;
	Pg(2,2)=1;
	Pg(3,NCd)=1;
	
	for (j=1; j<Nd-1; j++)
	{
		B(3*j-1,3*j-1)=-RDd(j-1);
		B(3*j-1,3*j)=3;
		B(3*j-1,3*j+1)=2;
		B(3*j-1,3*j+2)=1;
		
		B(3*j,3*j-2)=-2*pow(RDd(j-1),2);
		B(3*j,3*j)=6;
		B(3*j,3*j+1)=2;
		
		B(3*j+1,3*j)=1;
		B(3*j+1,3*j+1)=1;
		B(3*j+1,3*j+2)=1;
		
		Ps(3*j+1,ind_xlims(1)+j)=1;
		Ps(3*j+1,ind_xlims(1)+j+1)=-1;
		
		Pg(4*j,3*j)=1;
		Pg(4*j+1,3*j+1)=1;
		Pg(4*j+2,3*j+2)=1;
		Pg(4*j+3,NCd+j)=1;
	}
	
	j=Nd-1;
	B(3*j-1,3*j-1)=-RDd(j-1);
	B(3*j-1,3*j)=3;
	B(3*j-1,3*j+1)=2;
	
	B(3*j,3*j-2)=-2*pow(RDd(j-1),2);
	B(3*j,3*j)=6;
	B(3*j,3*j+1)=2;
	
	B(3*j+1,3*j)=1;
	B(3*j+1,3*j+1)=1;
	
	Ps(3*j+1,ind_xlims(1)+j)=1;
	Ps(3*j+1,ind_xlims(1)+j+1)=-1;
	
	Pg(4*j,3*j)=1;
	Pg(4*j+1,3*j+1)=1;
	Pg(4*j+3,NCd+j)=1;
	
	IB.eye(NCd,NCd);
	invB=solve(B,IB);
	
	PA=IA.submat(ind_xlims(1)+1,0,Nx-1,Nx-1);
	
	Lg=join_vert(invB*Ps,PA);
	
	mat Md=Pg*Lg;
	
	M=join_vert(Mg,Mc);
	M=join_vert(M,Md);
	
	return true;
}

bool OmegaMaxEnt_data::spline_matrix_G_part(vec x, uvec ind_xlims, vec xs, mat &Mspl)
{
	int Nx=x.n_rows;
	
	int Ng=ind_xlims(0)+1;
	int Nc=ind_xlims(1)-ind_xlims(0);
	int Nd=Nx-ind_xlims(1);
	
	double x0g=xs(0);
	double x0d=xs(1);
	
	
	//	 solve the spline in the interval -inf,wl]
	int NSg=3*Ng-1;
	mat Aspl=zeros<mat>(NSg,NSg);
	
	// Dg=zeros(NSg,1);
	// Dg(1)=S(1);
	
	double u=1.0/(x(0)-x0g);
	
	Aspl(0,0)=pow(u,3);
	Aspl(0,1)=pow(u,2);
	Aspl(1,0)=3*pow(u,2);
	Aspl(1,1)=2*u;
	Aspl(1,4)=-1;
	Aspl(2,0)=6*u;
	Aspl(2,1)=2;
	Aspl(2,3)=-2;
	
	int j;
	for (j=2; j<Ng; j++)
	{
		u=1/(x(j-1)-x0g)-1/(x(j-2)-x0g);
	 
	 //    Dg(3*j-2)=S(j)-S(j-1);
	 
		Aspl(3*j-3,3*j-4)=pow(u,3);
		Aspl(3*j-3,3*j-3)=pow(u,2);
		Aspl(3*j-3,3*j-2)=u;
	 
		Aspl(3*j-2,3*j-4)=3*pow(u,2);
		Aspl(3*j-2,3*j-3)=2*u;
		Aspl(3*j-2,3*j-2)=1;
		Aspl(3*j-2,3*j+1)=-1;
	 
		Aspl(3*j-1,3*j-4)=6*u;
		Aspl(3*j-1,3*j-3)=2;
		Aspl(3*j-1,3*j)=-2;
	}
	
	j=ind_xlims(0)+1;
	//Dg(3*j-2)=S(j)-S(j-1);
	u=1/(x(j-1)-x0g)-1/(x(j-2)-x0g);
	double ug=1/(x(j-1)-x0g);
	
	Aspl(3*j-3,3*j-4)=pow(u,3);
	Aspl(3*j-3,3*j-3)=pow(u,2);
	Aspl(3*j-3,3*j-2)=u;
	
	Aspl(3*j-2,3*j-4)=-3*pow(ug,2)*pow(u,2);
	Aspl(3*j-2,3*j-3)=-2*pow(ug,2)*u;
	Aspl(3*j-2,3*j-2)=-pow(ug,2);
	
	//Dg(3*j-1)=(S(j+1)-S(j-1))/(x(j+1)-x(j-1));
	//Cg=Ag\Dg;
	
	mat D(NSg,NSg);
	D.eye();
	
	mat Cg=solve(Aspl,D);
	
	mat T(NSg,Nx);
	T.zeros();
	T(0,0)=1;
	for (j=2; j<=Ng; j++)
	{
		T(3*j-3,j-2)=-1;
		T(3*j-3,j-1)=1;
	}
	T(3*Ng-2,Ng-2)=-1/(x(Ng)-x(Ng-2));
	T(3*Ng-2,Ng)=1/(x(Ng)-x(Ng-2));
	
	Cg=Cg*T;
	
	//solve the spline in the interval [wl,wr]
	int NSc=3*Nc-1;
	Aspl.zeros(NSc,NSc);
	
	//Dc=zeros(NSc,1);
	//Dc(1)=S(Ng+1)-(x(Ng+1)-x(Ng))*(S(Ng+1)-S(Ng-1))/(x(Ng+1)-x(Ng-1))-S(Ng);
	
	double Dx=x(Ng)-x(Ng-1);
	
	Aspl(0,0)=pow(Dx,3);
	Aspl(0,1)=pow(Dx,2);
	
	//Dc(2)=-(S(Ng+1)-S(Ng-1))/(x(Ng+1)-x(Ng-1));
	
	Aspl(1,0)=3*pow(Dx,2);
	Aspl(1,1)=2*Dx;
	Aspl(1,4)=-1;
	
	Aspl(2,0)=6*Dx;
	Aspl(2,1)=2;
	Aspl(2,3)=-2;
	
	for (j=2; j<Nc; j++)
	{
		Dx=x(Ng+j-1)-x(Ng+j-2);
	 
	 //    Dc(3*j-2)=S(Ng+j)-S(Ng+j-1);
	 
		Aspl(3*j-3,3*j-4)=pow(Dx,3);
		Aspl(3*j-3,3*j-3)=pow(Dx,2);
		Aspl(3*j-3,3*j-2)=Dx;
	 
		Aspl(3*j-2,3*j-4)=3*pow(Dx,2);
		Aspl(3*j-2,3*j-3)=2*Dx;
		Aspl(3*j-2,3*j-2)=1;
		Aspl(3*j-2,3*j+1)=-1;
	 
		Aspl(3*j-1,3*j-4)=6*Dx;
		Aspl(3*j-1,3*j-3)=2;
		Aspl(3*j-1,3*j)=-2;
	}
	
	j=Nc;
	Dx=x(Ng+j-1)-x(Ng+j-2);
	//Dc(3*j-2)=S(Ng+j)-S(Ng+j-1);
	Aspl(3*j-3,3*j-4)=pow(Dx,3);
	Aspl(3*j-3,3*j-3)=pow(Dx,2);
	Aspl(3*j-3,3*j-2)=Dx;
	
	//Dc(3*j-1)=(S(Ng+j+1)-S(Ng+j-1))/(x(Ng+j+1)-x(Ng+j-1));
	Aspl(3*j-2,3*j-4)=3*pow(Dx,2);
	Aspl(3*j-2,3*j-3)=2*Dx;
	Aspl(3*j-2,3*j-2)=1;
	
	//Cc=A\Dc;
	
	D.eye(NSc,NSc);
	mat Cc=solve(Aspl,D);
	
	T.zeros(NSc,Nx);
	T(0,Ng-2)=(x(Ng)-x(Ng-1))/(x(Ng)-x(Ng-2));
	T(0,Ng-1)=-1;
	T(0,Ng)=1-(x(Ng)-x(Ng-1))/(x(Ng)-x(Ng-2));
	T(1,Ng-2)=1/(x(Ng)-x(Ng-2));
	T(1,Ng)=-1/(x(Ng)-x(Ng-2));
	for (j=2; j<=Nc; j++)
	{
		T(3*j-3,Ng+j-2)=-1;
		T(3*j-3,Ng+j-1)=1;
	}
	T(3*Nc-2,ind_xlims(1)-1)=-1/(x(ind_xlims(1)+1)-x(ind_xlims(1)-1));
	T(3*Nc-2,ind_xlims(1)+1)=1/(x(ind_xlims(1)+1)-x(ind_xlims(1)-1));
	
	Cc=Cc*T;
	
	mat Cc2(NSc+1,Nx);
	Cc2.zeros();
	
	Cc2.rows(0,1)=Cc.rows(0,1);
	Cc2.rows(3,NSc)=Cc.rows(2,NSc-1);
	Cc2(2,Ng-2)=-1/(x(Ng)-x(Ng-2));
	Cc2(2,Ng)=1/(x(Ng)-x(Ng-2));
	
	//solve the spline in the interval [wr,inf
	int NSd=3*Nd-1;
	Aspl.zeros(NSd,NSd);
	//Dd=zeros(NSd,1);
	
	u=1/(x(ind_xlims(1))-x0d)-1/(x(ind_xlims(1)+1)-x0d);
	double ud=1/(x(ind_xlims(1))-x0d);
	
	//Dd(1)=(S(ind_xlims(2)+1)-S(ind_xlims(2)-1))/(x(ind_xlims(2)+1)-x(ind_xlims(2)-1));
	Aspl(0,0)=-3*pow(ud,2)*pow(u,2);
	Aspl(0,1)=-2*pow(ud,2)*u;
	Aspl(0,2)=-pow(ud,2);
	
	//Dd(2)=S(ind_xlims(2))-S(ind_xlims(2)+1);
	Aspl(1,0)=pow(u,3);
	Aspl(1,1)=pow(u,2);
	Aspl(1,2)=u;
	
	for (j=1; j<Nd-1; j++)
	{
		u=1/(x(ind_xlims(1)+j)-x0d)-1/(x(ind_xlims(1)+j+1)-x0d);
		
		Aspl(3*j-1,3*j-1)=-1;
		Aspl(3*j-1,3*j)=3*pow(u,2);
		Aspl(3*j-1,3*j+1)=2*u;
		Aspl(3*j-1,3*j+2)=1;
		
		Aspl(3*j,3*j-2)=-2;
		Aspl(3*j,3*j)=6*u;
		Aspl(3*j,3*j+1)=2;
	 
	 //    Dd(3*j+2)=S(ind_xlims(2)+j)-S(ind_xlims(2)+j+1);
		Aspl(3*j+1,3*j)=pow(u,3);
		Aspl(3*j+1,3*j+1)=pow(u,2);
		Aspl(3*j+1,3*j+2)=u;
	}
	
	j=Nd-1;
	u=1/(x(ind_xlims(1)+j)-x0d);
	
	Aspl(3*j-1,3*j-1)=-1;
	Aspl(3*j-1,3*j)=3*pow(u,2);
	Aspl(3*j-1,3*j+1)=2*u;
	
	Aspl(3*j,3*j-2)=-2;
	Aspl(3*j,3*j)=6*u;
	Aspl(3*j,3*j+1)=2;
	
	//Dd(3*j+2)=S(ind_xlims(2)+j);
	Aspl(3*j+1,3*j)=pow(u,3);
	Aspl(3*j+1,3*j+1)=pow(u,2);
	
	//Cd=A\Dd;
	
	D.eye(NSd,NSd);
	mat Cd=solve(Aspl,D);
	
	T.zeros(NSd,Nx);
	T(0,ind_xlims(1)-1)=-1/(x(ind_xlims(1)+1)-x(ind_xlims(1)-1));
	T(0,ind_xlims(1)+1)=1/(x(ind_xlims(1)+1)-x(ind_xlims(1)-1));
	for (j=1; j<Nd; j++)
	{
		T(3*j-2,ind_xlims(1)+j-1)=1;
		T(3*j-2,ind_xlims(1)+j)=-1;
	}
	j=Nd;
	T(3*j-2,ind_xlims(1)+j-1)=1;
	
	Cd=Cd*T;
	
	int NS_tot=NSg+NSc+NSd+1;
	
	Mspl.zeros(NS_tot,Nx);
	Mspl.rows(0,NSg-1)=Cg;
	Mspl.rows(NSg,NSg+NSc)=Cc2;
	Mspl.rows(NSg+NSc+1,NS_tot-1)=Cd;
	
	return true;
}




bool OmegaMaxEnt_data::general_normal(vec x, double x0, double s0, double p, vec &F)
{
	if (p<=0 || s0<=0)
	{
		if (p<0)
			cout<<"general_normal(): power must be greater than 0\n";
		else
			cout<<"general_normal(): width must be greater than 0\n";
		return false;
	}
	
	F=2*PI*p*exp(-pow(abs((x-x0)/s0),p))/(2*s0*tgamma(1.0/p));
	
	return true;
}



double OmegaMaxEnt_data::general_normal_val(double x, void *par[])
{
	double *x0=reinterpret_cast<double*>(par[0]);
	double *s0=reinterpret_cast<double*>(par[1]);
	double *p=reinterpret_cast<double*>(par[2]);
	
	if (p[0]<=0 || s0[0]<=0)
	{
		if (p[0]<0)
			cout<<"general_normal(): power must be greater than 0\n";
		else
			cout<<"general_normal(): width must be greater than 0\n";
		return false;
	}
	
	return 2*PI*p[0]*exp(-pow(abs((x-x0[0])/s0[0]),p[0]))/(2*s0[0]*tgamma(1.0/p[0]));
}

bool OmegaMaxEnt_data::sum_gaussians(vec x, rowvec x0, rowvec s0, rowvec wgt, vec &F)
{
	int Nx=x.n_rows;
	int Np=x0.n_cols;
	
	if (s0.n_cols<Np || wgt.n_cols<Np)
	{
		cout<<"sum_gaussians(): third or fourth parameter has too few elements\n";
		return false;
	}
	
	F.zeros(Nx);

	double wgtTot=0;

	for (int j=0; j<Np; j++)
	{
		F=F+wgt(j)*exp(-pow(x-x0(j),2)/(2*pow(s0(j),2)))/s0(j);
		wgtTot=wgtTot+wgt(j);
	}

	F=sqrt(2.0*PI)*F/wgtTot;
	
	return true;
}

bool OmegaMaxEnt_data::default_model_val_G(vec x, vec x0, vec coeffs, vec gaussians_params, vec &dm)
{
	int Nx=x.n_rows;
	int Nx0=x0.n_rows;

	double wcg=gaussians_params(0);
	double C1g=gaussians_params(1);
	double C2g=gaussians_params(2);
	double wcd=gaussians_params(3);
	double C1d=gaussians_params(4);
	double C2d=gaussians_params(5);
	
	vec diff_x=x.rows(1,Nx-1)-x.rows(0,Nx-2);
	if (diff_x.min()<=0)
	{
		cout<<"default_model_val_G(): values in position vector are not strictly increasing\n";
		return false;
	}
	
	dm.zeros(Nx);
	
	int j=0;
	while (x(j)<x0(0) && j<Nx-1) j++;
	int jl=j;
	while (x(j)<x0(Nx0-1) && j<Nx-1) j++;
	int jr=j-1;
	int Nm=jr-jl+1, Nl=jl;
	int Nr=Nx-Nm-Nl;
	
	if (Nl)
	{
		vec vl=exp(-pow(x.rows(0,jl-1)-wcg,2)/C1g)/C2g;
		dm.rows(0,jl-1)=vl;
	}
	if (Nm)
	{
		vec vm;
		spline_val(x.rows(jl,jr), x0, coeffs, vm);
		dm.rows(jl,jr)=vm;
	}
	if (Nr)
	{
		vec vr=exp(-pow(x.rows(jr+1,Nx-1)-wcd,2)/C1d)/C2d;
		dm.rows(jr+1,Nx-1)=vr;
	}
	
	return true;
}

bool OmegaMaxEnt_data::spline_val(vec x, vec x0, vec coeffs, vec &s)
{
	int Nx=x.n_rows;
	int Nx0=x0.n_rows;
	int Ns=Nx0-1;
	int Nc=4*Ns;
	
	if (coeffs.n_rows<Nc)
	{
		cout<<"spline_val(): either the coefficients vector is incomplete or the position vector is too long\n";
		return false;
	}

	s.zeros(Nx);
	
	double a, b, c, d, Dx;

	int j, l;
	for (j=0; j<Nx; j++)
	{
		if (x(j)>=x0(0) && x(j)<x0(Nx0-1))
		{
			l=0;
			while ( x(j)>=x0(l+1) && l<Ns-1) l++;
			a=coeffs(4*l);
			b=coeffs(4*l+1);
			c=coeffs(4*l+2);
			d=coeffs(4*l+3);
			
			Dx=x(j)-x0(l);
			s(j)=a*pow(Dx,3)+b*pow(Dx,2)+c*Dx+d;
		}
		else if (x(j)==x0(Nx0-1))
		{
			if (coeffs.n_rows>Nc)
			{
				s(j)=coeffs[Nc];
			}
			else
			{
				l=Ns-1;
				a=coeffs(4*l);
				b=coeffs(4*l+1);
				c=coeffs(4*l+2);
				d=coeffs(4*l+3);
				
				Dx=x(j)-x0(l);
				s(j)=a*pow(Dx,3)+b*pow(Dx,2)+c*Dx+d;
			}
		}
		else
		{
			cout<<"spline_val(): the provided position vector has a value outside the boundaries of the spline.\n";
			return false;
		}
	}
	
	return true;
}

//!Compute the coefficients of the cubic spline for V(x) known at N0 positions x0, size of coeffs must be 4*(N0-1).
//!upon entry, coeffs[0] and coeffs[1] must contain the derivatives of V(x) at x0[0] and x[N0-1],
//!upon exit, the form of coeffs is {a_0, b_0, c_0, d_0, ... a_(N0-2), b_(N0-2), c_(N0-2), d_(N0-2)}.
//!the spline values are given by S_i(x)=a_i(x-x0[i])^3+b_i(x-x0[i])^2+c_i(x-x0[i])+d_i
void OmegaMaxEnt_data::spline_coeffs_rel(double *x0, double *V, int N0, double *coeffs)
{
	int j;
	
	int NS=N0-1;
	int N=3*NS-1;
	int KL=3;
	int KU=2;
	int NA=2*KL+KU+1;
	int SA=NA*N;
	
	double *A=new double[SA];
	for (j=0; j<SA; j++) A[j]=0;
	int *P=new int[N];
	double *coeffs_tmp=new double[N];
	for (j=0; j<N; j++) coeffs_tmp[j]=0;
	
	double dV1=coeffs[0], dVN0=coeffs[1];
	double x;
	
	x=x0[1]-x0[0];
	
	A[KL+KU]=(double) (x*x*x);
	A[KL+KU+1]=(double) (-6.0*x);
	A[KL+KU+2]=(double) (-3.0*x*x);
	
	A[NA+KL+KU-1]=(double) (x*x);
	A[NA+KL+KU]=(double) (-2.0);
	A[NA+KL+KU+1]=(double) (-2.0*x);
	
	coeffs_tmp[0]=(double) (V[1]-V[0]-dV1*x);
	coeffs_tmp[2]=(double)dV1;
	
	for (j=1; j<NS-1; j++)
	{
		x=x0[j+1]-x0[j];
		
		A[(3*j-1)*NA+KL+KU+1]=(double) (x*x*x);
		A[(3*j-1)*NA+KL+KU+2]=(double) (-6.0*x);
		A[(3*j-1)*NA+KL+KU+3]=(double) (-3.0*x*x);
		
		A[3*j*NA+KL]=2.0;
		A[3*j*NA+KL+KU]=(double) (x*x);
		A[3*j*NA+KL+KU+1]=(double) (-2.0);
		A[3*j*NA+KL+KU+2]=(double) (-2.0*x);
		
		A[(3*j+1)*NA+KL]=1.0;
		A[(3*j+1)*NA+KL+1]=(double) x;
		A[(3*j+1)*NA+KL+KU+1]=-1.0;
		
		coeffs_tmp[3*j]=(double) (V[j+1]-V[j]);
	}
	
	j=NS-1;
	x=x0[j+1]-x0[j];
	
	A[(3*j-1)*NA+KL+KU+1]=(double) (x*x*x);
	A[(3*j-1)*NA+KL+KU+2]=(double) (3.0*x*x);
	
	A[3*j*NA+KL]=2.0;
	A[3*j*NA+KL+KU]=(double) (x*x);
	A[3*j*NA+KL+KU+1]=(double) (2.0*x);
	
	A[(3*j+1)*NA+KL]=1.0;
	A[(3*j+1)*NA+KL+1]=(double) x;
	A[(3*j+1)*NA+KL+KU]=1.0;
	
	coeffs_tmp[3*j]=(double)(V[j+1]-V[j]);
	coeffs_tmp[3*j+1]=(double)dVN0;
	
	int NRHS=1;
	int INFO=0;
	dgbsv_(&N, &KL, &KU, &NRHS, A, &NA, P, coeffs_tmp, &N, &INFO );
	
	coeffs[0]=coeffs_tmp[0];
	coeffs[1]=coeffs_tmp[1];
	coeffs[2]=dV1;
	coeffs[3]=V[0];
	for (j=1; j<NS; j++)
	{
		coeffs[4*j]=coeffs_tmp[3*j-1];
		coeffs[4*j+1]=coeffs_tmp[3*j];
		coeffs[4*j+2]=coeffs_tmp[3*j+1];
		coeffs[4*j+3]=V[j];
	}
	
	if (INFO)	cout<<"spline_coeffs(): INFO:  "<<INFO<<'\n';
	
	delete [] A;
	delete [] P;
	delete [] coeffs_tmp;
}






bool OmegaMaxEnt_data::set_exp_step_grid()
{
	double dw;

	double f_width_exp_grid=4;
	
	if (SW_set && SC_set && !main_spectral_region_set)
	{
		wl=SC-SW/2;
		wr=SC+SW/2;
		main_spectral_region_set=true;
	}
	else if (std_omega && SC_set && !main_spectral_region_set)
	{
		SW=f_SW_std_omega*std_omega;
		SW_set=true;
		wl=SC-SW/2;
		wr=SC+SW/2;
		main_spectral_region_set=true;
	}
	else if (!main_spectral_region_set)
	{
		cout<<"The central part of the grid cannot be defined. Not enough information available.\n";
		return false;
	}
//	if (!w_origin_in.size())
	if (!w_origin_set)
	{
		w_origin=SC;
		w_origin_set=true;
	}
	
	double Dwh=f_width_exp_grid*SW/2;
	
	if (step_omega_in.size())
	{
		dw=step_omega;
		if (peak_exists)
		{
			if (dw>dw_peak)
			{
				cout<<"warning: step is larger than the estimated width of the peak at low energy\n";
			}
		}
	}
	else
	{
		dw=SW/(f_SW_std_omega*Rmin_SW_dw);
	//	dw=Dwh/(2*Nw_min);
		if (peak_exists)
		{
			if (dw_peak<dw)
				dw=dw_peak;
		}
	}
	
	if (w_origin<wl+dw) w_origin=wl+dw;
	else if (w_origin>wr-dw) w_origin=wr-dw;
	
	int j;
	
	int Nwmax=SW/dw;
	
	double w0=SC;
	
	double dwtmp;
	vec wltmp(Nwmax);
	wltmp(0)=w_origin;
	j=0;
	while (wltmp(j)>wl)// && j<Nwmax-1)
	{
		dwtmp=dw*exp(pow((w0-wltmp(j))/Dwh,2)/2);
		//if (abs(wltmp(j)-dwtmp-wl)<abs(wltmp(j)-wl))
		//{
			wltmp(j+1)=wltmp(j)-dwtmp;
			j++;
		//}
		//else
		//	wl=wltmp(j);
	}
	int Nwcl=j+1;
	
	vec wrtmp(Nwmax);
	wrtmp(0)=w_origin;
	j=0;
	while (wrtmp(j)<wr)// && j<Nwmax-1)
	{
		dwtmp=dw*exp(pow((w0-wrtmp(j))/Dwh,2)/2);
		//if (abs(wrtmp(j)+dwtmp-wr)<abs(wrtmp(j)-wr))
		//{
			wrtmp(j+1)=wrtmp(j)+dwtmp;
			j++;
		//}
		//else
		//	wr=wrtmp(j);
	}
	int Nwcr=j+1;
	
	wc=join_vert(flipud(wltmp.rows(1,Nwcl-1)),wrtmp.rows(0,Nwcr-1));
	Nwc=wc.n_rows;
	wl=wc(0);
	wr=wc(Nwc-1);
	wc_exists=true;
	
	
	return true;
}

bool OmegaMaxEnt_data::set_quad_step_grid()
{
	double dw;
	
	f_width_grid_dens=1;
	
	if (SW_set && SC_set && !main_spectral_region_set)
	{
		wl=SC-SW/2;
		wr=SC+SW/2;
		main_spectral_region_set=true;
	}
	else if (std_omega && SC_set && !main_spectral_region_set)
	{
		SW=f_SW_std_omega*std_omega;
		SW_set=true;
		wl=SC-SW/2;
		wr=SC+SW/2;
		main_spectral_region_set=true;
	}
	else if (!main_spectral_region_set)
	{
		cout<<"The central part of the grid cannot be defined. Not enough information available.\n";
		return false;
	}
//	if (!w_origin_in.size())
	if (!w_origin_set)
	{
		w_origin=SC;
		w_origin_set=true;
	}
	
	double Dwh=f_width_grid_dens*SW/2;
	
	if (step_omega_in.size())
	{
		dw=step_omega;
		if (peak_exists)
		{
			if (dw>dw_peak)
			{
				cout<<"warning: step is larger than the estimated width of the peak at low energy\n";
			}
		}
	}
	else
	{
		dw=SW/(f_SW_std_omega*Rmin_SW_dw);
		//dw=Dwh/(2*Nw_min);
		if (peak_exists)
		{
			if (dw_peak<dw)
				dw=dw_peak;
		}
	}
	
	if (w_origin<wl+dw) w_origin=wl+dw;
	else if (w_origin>wr-dw) w_origin=wr-dw;
	
	int j;
	
	int Nwmax=6*Dwh/dw;
	
	double Dwh2=Dwh*Dwh;
	double w0=w_origin;
	
	double dwtmp;
	vec wltmp(Nwmax);
	wltmp(0)=w_origin;
	j=0;
	while (wltmp(j)>wl && j<Nwmax-1)
	{
		dwtmp=dw*(1+pow(w0-wltmp(j),2)/Dwh2);
		wltmp(j+1)=wltmp(j)-dwtmp;
		j++;
	}
	int Nwcl=j+1;
	
	vec wrtmp(Nwmax);
	wrtmp(0)=w_origin;
	j=0;
	while (wrtmp(j)<wr && j<Nwmax-1)
	{
		dwtmp=dw*(1+pow(wrtmp(j)-w0,2)/Dwh2);
		wrtmp(j+1)=wrtmp(j)+dwtmp;
		j++;
	}
	int Nwcr=j+1;
	
	wc=join_vert(flipud(wltmp.rows(1,Nwcl-1)),wrtmp.rows(0,Nwcr-1));
	Nwc=wc.n_rows;
	wl=wc(0);
	wr=wc(Nwc-1);
	
	double wmin=SC-f_w_range*SW/2.0;
	
	j=Nwcl-1;
	while (wltmp(j)>wmin && j+1-Nwcl<Nwc)
	{
		dwtmp=dw*(1+pow(w0-wltmp(j),2)/Dwh2);
		wltmp(j+1)=wltmp(j)-dwtmp;
		j++;
	}
	int Nwl=j+1;
	
	double wmax=SC+f_w_range*SW/2.0;
	
	j=Nwcr-1;
	while (wrtmp(j)<wmax && j+1-Nwcr<Nwc)
	{
		dwtmp=dw*(1+pow(wrtmp(j)-w0,2)/Dwh2);
		wrtmp(j+1)=wrtmp(j)+dwtmp;
		j++;
	}
	int Nwr=j+1;
	
	w=join_vert(flipud(wltmp.rows(1,Nwl-1)),wrtmp.rows(0,Nwr-1));
	Nw=w.n_rows;
	wc_exists=true;
	w_exists=true;
	Nw_lims.zeros(2);
	Nw_lims(0)=Nwl-Nwcl;
	Nw_lims(1)=Nw_lims(0)+Nwc-1;
	
	w0l=SC;
	w0r=SC;
	ws.zeros(2);
	ws(0)=w0l;
	ws(1)=w0r;
	
	return true;
}







bool OmegaMaxEnt_data::non_uniform_frequency_grid(rowvec w_steps, rowvec wlims, double w0, vec R, vec &grid)
{
	bool ordered=true, grid_set=true;
	
	uint Ndw=w_steps.n_cols;
	uint Nlims=wlims.n_cols;
	
	if (Nlims<=Ndw)
		Ndw=Nlims-1;
	else if (Nlims>Ndw+1)
		Nlims=Ndw+1;
	
	double RW=R(0);
	double RWD=R(1);
	
	int j;
	for (j=1; j<Nlims; j++)
	{
		if (wlims(j)<wlims(j-1)) ordered=false;
	}
	
	if (ordered )
	{
		double D, wc_j, dwtmp;
		int max_grid_size=ceil((wlims(Nlims-1)-wlims(0))/w_steps.min());
		int j0, l, Nm, Np;
		if (w0>wlims(0) && w0<wlims(Nlims-1))
		{
			j=0;
			while (w0>=wlims(j+1))	j=j+1;
			j0=j;
			double wstart=w0;
//			double wstart=(wlims(j)+wlims(j+1))/2.0;
			vec wp(max_grid_size);
			wp(0)=wstart;
			l=0;
			for (j=j0; j<Ndw-1; j++)
			{
				D=(wlims(j+2)-wlims(j))/(RW*RWD);
				wc_j=wlims(j+1)+(wlims(j+2)-2*wlims(j+1)+wlims(j))/(2*RW);
				while (wp(l)<(wlims(j+2)+wlims(j+1))/2.0)
				{
					dwtmp=w_steps(j+1)+(w_steps(j)-w_steps(j+1))/(exp((wp(l)-wc_j)/D)+1.0);
					wp(l+1)=wp(l)+dwtmp;
					l=l+1;
				}
			}
			while (wp(l)<wlims(Nlims-1))
			{
				wp(l+1)=wp(l)+w_steps(Ndw-1);
				l=l+1;
			}
			Np=l+1;
			vec wm(max_grid_size);
			wm(0)=wstart;
			l=0;
			for (j=j0-1; j>=0; j--)
			{
				D=(wlims(j+2)-wlims(j))/(RW*RWD);
				wc_j=wlims(j+1)+(wlims(j+2)-2*wlims(j+1)+wlims(j))/(2*RW);
				while (wm(l)>(wlims(j)+wlims(j+1))/2)
				{
					dwtmp=w_steps(j+1)+(w_steps(j)-w_steps(j+1))/(exp((wm(l)-wc_j)/D)+1.0);
					wm(l+1)=wm(l)-dwtmp;
					l=l+1;
				}
			}
			while (wm(l)>wlims(0))
			{
				wm(l+1)=wm(l)-w_steps(0);
				l=l+1;
			}
			Nm=l+1;
			grid.zeros(Nm+Np-1);
			grid.rows(0,Nm-1)=flipud(wm.rows(0,Nm-1));
			grid.rows(Nm,Nm+Np-2)=wp.rows(1,Np-1);
		}
		else if (w0==wlims(0))
		{
			vec wp(max_grid_size);
			wp(0)=w0;
			j0=0;
			l=0;
			for (j=j0; j<Ndw-1; j++)
			{
				D=(wlims(j+2)-wlims(j))/(RW*RWD);
				wc_j=wlims(j+1)+(wlims(j+2)-2*wlims(j+1)+wlims(j))/(2*RW);
				while (wp(l)<(wlims(j+2)+wlims(j+1))/2.0)
				{
					dwtmp=w_steps(j+1)+(w_steps(j)-w_steps(j+1))/(exp((wp(l)-wc_j)/D)+1.0);
					wp(l+1)=wp(l)+dwtmp;
					l=l+1;
				}
			}
			while (wp(l)<wlims(Nlims-1))
			{
				wp(l+1)=wp(l)+w_steps(Ndw-1);
				l=l+1;
			}
			Np=l+1;
			grid=wp.rows(0,Np-1);
		}
		else if (w0==wlims(Nlims-1))
		{
			j0=Nlims-2;
			vec wm(max_grid_size);
			wm(0)=w0;
			l=0;
			for (j=j0-1; j<=0; j--)
			{
				D=(wlims(j+2)-wlims(j))/(RW*RWD);
				wc_j=wlims(j+1)+(wlims(j+2)-2*wlims(j+1)+wlims(j))/(2*RW);
				while (wm(l)>(wlims(j)+wlims(j+1))/2.0)
				{
					dwtmp=w_steps(j+1)+(w_steps(j)-w_steps(j+1))/(exp((wm(l)-wc_j)/D)+1.0);
					wm(l+1)=wm(l)-dwtmp;
					l=l+1;
				}
			}
			while (wm(l)>wlims(0))
			{
				wm(l+1)=wm(l)-w_steps(0);
				l=l+1;
			}
			Nm=l+1;
			grid=flipud(wm.rows(0,Nm-1));
		}
		else
		{
			cout<<"grid origin is outside grid boundaries\n";
			grid_set=false;
		}
	}
	else
	{
		cout<<"non_uniform_frequency_grid(): interval boundaries are not strictly increasing values\n";
		grid_set=false;
	}
	
	return grid_set;
}






bool OmegaMaxEnt_data::test_low_energy_peak_fermions()
{
//	if (Nwn_test_metal>Nn-2)	Nwn_test_metal=Nn-2;
	
	peak_exists=false;
//	bool test_peak=false;
	
//	vec D2Gi=Gi.rows(0,Nn-3)-2.0*Gi.rows(1,Nn-2)+Gi.rows(2,Nn-1);
//	if (max(D2Gi.rows(0,Nwn_test_metal-1))<0 && Nn>20) test_peak=true;

//	if (test_peak)
//	{
//		cout<<"Im(G) seems to diverge at low Matsubara frequency. Looking for a sharp peak in the spectral function at low energy...\n";
		cout<<"Looking for a peak in the spectral function at low energy...\n";
		
		int nmax=n(Nn-1);
		
		int NCnmin=2;
		int NCnmax=3;
		int NNCn=NCnmax-NCnmin+1;
		ivec NCn=linspace<ivec>(NCnmin,NCnmax,NNCn);
		
		int NCpmin=1;
		int NCpmax=15;
		if (NCpmax>(nmax-NCn.max()-2))
		{
			NCpmax=nmax-NCn.max()-2;
		}
		int NNCp=NCpmax-NCpmin+1;
		ivec NCp=linspace<ivec>(NCpmin,NCpmax,NNCp);
		
		int DNwn=2;
		
		int p0_min=1;
		int p0_max=1;
		int Np=p0_max-p0_min+1;
		ivec p0=linspace<ivec>(p0_min,p0_max,Np);
		
		imat NCmin(Np,2);
		
		vec norm_lf_min(Np,fill::zeros);
		vec std_peak_min(Np,fill::zeros);
		vec mean_omega_lf_min(Np,fill::zeros);
		vec mean_omega2_lf_min(Np,fill::zeros);
		vec chi2_lf_min(Np,fill::zeros);
		vec varM2_q(Np,fill::zeros);
		vec varM0_q(Np,fill::zeros);
		
		mat X, CG, P, invCG, AM, BM, AMP, BMP, MP, Mtmp, chi2tmp;
		mat norm_lf(NNCp,NNCn), std_peak(NNCp,NNCn), mean_omega_lf(NNCp,NNCn), mean_omega2_lf(NNCp,NNCn), chi2_lf(NNCp,NNCn), M2_lf(NNCp,NNCn);
		int j, l, m, p, Nfit;
		uword q;
		vec Gtmp, Glftmp, diffG;
		rowvec maxX;
		double vartmp;
		
//		char UPLO='U';
//		int NA, NRHS=1, INFO;
		
//		mat test_M, AM2;
//		vec BM2;


		for (q=0; q<Np; q++)
		{
			norm_lf.zeros();
			std_peak.zeros();
			mean_omega_lf.zeros();
			mean_omega2_lf.zeros();
			chi2_lf.zeros();
			M2_lf.zeros();
			

			for (l=0; l<NNCp; l++)
			{
				for (m=0; m<NNCn; m++)
				{
					Nfit=NCp(l)+NCn(m)+1+DNwn;
					
					X.zeros(2*Nfit,2*(NCp(l)+NCn(m)+1));
					
					for (j=NCp(l); j>=-NCn(m); j--)
					{
						for (p=p0(q); p<=Nfit+p0(q)-1; p++)
						{
							X(2*(p-p0(q)),2*NCp(l)-2*j+1)=pow(-1,j)*pow(wn(p-1),2*j);
							X(2*(p-p0(q))+1,2*NCp(l)-2*j)=pow(-1,j)*pow(wn(p-1),2*j+1);
						}
					}

					Gtmp=Gchi2.rows(2*p0(q)-2,2*(Nfit+p0(q))-3);
					CG=COV.submat(2*p0(q)-2,2*p0(q)-2,2*(Nfit+p0(q))-3,2*(Nfit+p0(q))-3);

					invCG=inv(CG);
					//	invCG=inv_sympd(CG);
					AM=(X.t())*invCG*X;
					BM=(X.t())*invCG*Gtmp;
//					BM2=BM;
//					AM2=AM;
					
					Mtmp=solve(AM,BM);
					
//					cout<<l<<" "<<m<<endl;
					
//					NA=AM.n_rows;
//					dposv_(&UPLO, &NA, &NRHS, AM.memptr(), &NA, BM.memptr(), &NA, &INFO);
//					Mtmp=BM;
					
//					dposv_(&UPLO, &NA, &NRHS, AM2.memptr(), &NA, BM2.memptr(), &NA, &INFO);
//					test_M.zeros(NA,2);
//					test_M.col(0)=Mtmp;
//					test_M.col(1)=BM2;
//					cout<<test_M<<endl;

//					maxX=max(abs(X),0);
//					P=diagmat(1.0/maxX);
//					AMP=(P*(AM*P)+(P*AM)*P)/2;
//					BMP=P*BM;
//					dposv_(&UPLO, &NA, &NRHS, AMP.memptr(), &NA, BMP.memptr(), &NA, &INFO);
//					MP=BMP;
//					MP=solve(AMP,BMP);
//					Mtmp=P*MP;
					
					Glftmp=X*Mtmp;
					
					diffG=Gtmp-Glftmp;
					
					chi2tmp=((diffG.t())*invCG*diffG)/(2.0*Nfit);
					
					chi2_lf(l,m)=chi2tmp(0,0);

					norm_lf(l,m)=Mtmp(2*NCp(l)+2);
					mean_omega_lf(l,m)=Mtmp(2*NCp(l)+3)/norm_lf(l,m);
					mean_omega2_lf(l,m)=Mtmp(2*NCp(l)+4)/norm_lf(l,m);
					vartmp=mean_omega2_lf(l,m)-pow(mean_omega_lf(l,m),2);
					if (vartmp>0)
						std_peak(l,m)=sqrt(vartmp);
					M2_lf(l,m)=Mtmp(2*NCp(l)+4);

//					cout<<setw(30)<<"norm_lf: "<<setw(20)<<norm_lf(l,m)<<"	"<<BM2(2*NCp(l)+2)<<endl;
//					cout<<setw(30)<<"mean_omega_lf: "<<setw(20)<<mean_omega_lf(l,m)<<"	"<<BM2(2*NCp(l)+3)/BM2(2*NCp(l)+2)<<endl;
//					cout<<setw(30)<<"mean_omega2_lf: "<<setw(20)<<mean_omega2_lf(l,m)<<"	"<<BM2(2*NCp(l)+4)/BM2(2*NCp(l)+2)<<endl;
	
				}
			}

			for (l=0; l<NNCp; l++)
			{
				for (m=0; m<NNCn; m++)
				{
					if (std_peak(l,m)==0)
					{
						chi2_lf(l,m)=max(max(chi2_lf));
					}
				}
			}
			
//			cout<<"std_peak:\n"<<std_peak<<endl;
//			cout<<"norm_lf:\n"<<norm_lf<<endl;
//			cout<<"mean_omega_lf:\n"<<mean_omega_lf<<endl;
			
			uint NvM=2;
			uint NvarM=NNCp-2*NvM;
			mat varM0(NvarM,NNCn,fill::zeros);
			mat varM2(NvarM,NNCn,fill::zeros);
			for (j=NvM; j<NNCp-NvM; j++)
			{
				varM0.row(j-NvM)=var(norm_lf.rows(j-NvM,j+NvM));
				varM2.row(j-NvM)=var(M2_lf.rows(j-NvM,j+NvM));
			}
			
			uword indpM0, indnM0;
			
			double varM0_min=varM0.min(indpM0,indnM0);
			double varM2_min=varM2.min();
			
			l=indpM0+NvM;
			m=indnM0;
			
			NCmin(q,0)=l;
			NCmin(q,1)=m;
			
			norm_lf_min(q)=norm_lf(l,m);
			mean_omega_lf_min(q)=mean_omega_lf(l,m);
			mean_omega2_lf_min(q)=mean_omega2_lf(l,m);
			std_peak_min(q)=std_peak(l,m);
			chi2_lf_min(q)=chi2_lf(l,m);
			varM0_q(q)=varM0_min;
			varM2_q(q)=varM2_min;
			
		}

		double varM0_min=varM0_q.min(q);
		
		double peak_weight=norm_lf_min(q);
		//double peak_center=mean_omega_lf_min(q);
		double peak_width=std_peak_min(q);
//		double chi2_min=chi2_lf_min(q);
		//double m2_lf_min=mean_omega2_lf_min(q);
		
		if (peak_width>100*EPSILON && sqrt(varM0_min)/peak_weight<std_norm_peak_max && peak_weight>peak_weight_min*M0)
		{
			peak_exists=true;
			
			l=NCmin(q,0);
			m=NCmin(q,1);
			
			Nfit=NCp(l)+NCn(m)+1+DNwn;
			
			X.zeros(2*Nfit,2*(NCp(l)+NCn(m)+1));
			
			for (j=NCp(l); j>=-NCn(m); j--)
			{
				for (p=p0(q); p<=Nfit+p0(q)-1; p++)
				{
					X(2*(p-p0(q)),2*NCp(l)-2*j+1)=pow(-1,j)*pow(wn(p-1),2*j);
					X(2*(p-p0(q))+1,2*NCp(l)-2*j)=pow(-1,j)*pow(wn(p-1),2*j+1);
				}
			}
			
			Gtmp=Gchi2.rows(2*p0(q)-2,2*(Nfit+p0(q))-3);
			CG=COV.submat(2*p0(q)-2,2*p0(q)-2,2*(Nfit+p0(q))-3,2*(Nfit+p0(q))-3);
			
			invCG=inv(CG);
			//	invCG=inv_sympd(CG);
			AM=(X.t())*invCG*X;
			AM=0.5*(AM.t()+AM);
			BM=(X.t())*invCG*Gtmp;
			
			Mtmp=solve(AM,BM);
			
//			NA=AM.n_rows;
//			dposv_(&UPLO, &NA, &NRHS, AM.memptr(), &NA, BM.memptr(), &NA, &INFO);
//			Mtmp=BM;
			
//			rowvec maxX=max(abs(X),0);
//			mat P=diagmat(1.0/maxX);
//			AMP=(P*(AM*P)+(P*AM)*P)/2.0;
//			BMP=P*BM;
//			dposv_(&UPLO, &NA, &NRHS, AMP.memptr(), &NA, BMP.memptr(), &NA, &INFO);
//			MP=BMP;
//			MP=solve(AMP,BMP);
//			Mtmp=P*MP;
			
			Glftmp=X*Mtmp;
			
			mat COVpeak=inv(AM);
			
//			invAMP=inv(AMP);
//			COVpeak=((P*invAMP)*P+P*(invAMP*P))/2;
			
			/*
			double err_norm_peak=sqrt(COVpeak(2*NCp(l)+2,2*NCp(l)+2));
			double err_M1_peak=sqrt(COVpeak(2*NCp(l)+3,2*NCp(l)+3));
			double err_M2_peak=sqrt(COVpeak(2*NCp(l)+4,2*NCp(l)+4));
			double err_peak_position=err_M1_peak/peak_weight-err_norm_peak*peak_center/peak_weight;
			double err_std_peak=(-err_norm_peak*m2_lf_min/peak_weight+2*err_norm_peak*peak_center*peak_center/peak_weight+err_M2_peak/peak_weight-2*err_M1_peak*peak_center/peak_weight)/(2*peak_width);
			*/
			cout<<"Peak detected\n";
			cout<<"peak width: "<<peak_width<<endl;
			//cout<<"error on width: "<<err_std_peak<<endl;
			cout<<"peak weight: "<<peak_weight<<endl;
			//cout<<"error on weight: "<<err_norm_peak<<endl;
			//cout<<"peak position: "<<peak_center<<endl;
			//cout<<"error on position: "<<err_peak_position<<endl;
			
			dw_peak=peak_width/2.0;
		}
		else
		{
			cout<<"no peak found\n";
		}
//	}

	return peak_exists;
}




void OmegaMaxEnt_data::init_variables()
{
	NAprec=5;
	
    input_dir_prec=NULL;
    data_file_name_prec=NULL;
    boson_prec=NULL;
    tau_GF_prec=NULL;
    tem_prec=NULL;
//    signG_prec=NULL;
    M0_prec=NULL;
    M1_prec=NULL;
    errM1_prec=NULL;
    M2_prec=NULL;
    errM2_prec=NULL;
    M3_prec=NULL;
    errM3_prec=NULL;
	static_part_G_prec=NULL;
    col_Gr_prec=NULL;
    col_Gi_prec=NULL;
    error_file_prec=NULL;
    col_errGr_prec=NULL;
    col_errGi_prec=NULL;
    covar_re_re_file_prec=NULL;
    covar_im_im_file_prec=NULL;
    covar_re_im_file_prec=NULL;
    col_Gtau_prec=NULL;
    col_errGtau_prec=NULL;
    covar_tau_file_prec=NULL;
    cutoff_wn_prec=NULL;
    SW_prec=NULL;
    SC_prec=NULL;
    w_origin_prec=NULL;
    step_omega_prec=NULL;
    grid_omega_file_prec=NULL;
    use_grid_params_prec=NULL;
    omega_grid_params_prec=NULL;
    eval_moments_prec=NULL;
    maxM_prec=NULL;
//    peaked_default_model_prec=NULL;
    def_model_file_prec=NULL;
    init_spectr_func_file_prec=NULL;
//    interp_type_prec=NULL;
	output_dir_prec=NULL;
    w_sample_prec=NULL;
    alpha_init_prec=NULL;
	default_model_center_prec=NULL;
	default_model_width_prec=NULL;
	default_model_shape_prec=NULL;
	A_ref_file_prec=NULL;
	non_uniform_grid_prec=NULL;
    time_params_file=NULL;
    time_other_params_file=NULL;
}

void OmegaMaxEnt_data::free_variables()
{
	
    if (input_dir_prec) delete input_dir_prec;
    if (data_file_name_prec) delete data_file_name_prec;
    if (boson_prec) delete boson_prec;
    if (tau_GF_prec) delete tau_GF_prec;
    if (tem_prec) delete tem_prec;
//    if (signG_prec) delete signG_prec;
    if (M0_prec) delete M0_prec;
    if (M1_prec) delete M1_prec;
    if (errM1_prec) delete errM1_prec;
    if (M2_prec) delete M2_prec;
    if (errM2_prec) delete errM2_prec;
    if (M3_prec) delete M3_prec;
    if (errM3_prec) delete errM3_prec;
	if (static_part_G_prec) delete static_part_G_prec;
    if (col_Gr_prec) delete col_Gr_prec;
    if (col_Gi_prec) delete col_Gi_prec;
    if (error_file_prec) delete error_file_prec;
    if (col_errGr_prec) delete col_errGr_prec;
    if (col_errGi_prec) delete col_errGi_prec;
    if (covar_re_re_file_prec) delete covar_re_re_file_prec;
    if (covar_im_im_file_prec) delete covar_im_im_file_prec;
    if (covar_re_im_file_prec) delete covar_re_im_file_prec;
    if (col_Gtau_prec) delete col_Gtau_prec;
    if (col_errGtau_prec) delete col_errGtau_prec;
    if (covar_tau_file_prec) delete covar_tau_file_prec;
    if (cutoff_wn_prec) delete cutoff_wn_prec;
    if (SW_prec) delete SW_prec;
    if (SC_prec) delete SC_prec;
    if (w_origin_prec) delete w_origin_prec;
    if (step_omega_prec) delete step_omega_prec;
    if (grid_omega_file_prec) delete grid_omega_file_prec;
    if (use_grid_params_prec) delete use_grid_params_prec;
    if (omega_grid_params_prec) delete omega_grid_params_prec;
    if (eval_moments_prec) delete eval_moments_prec;
    if (maxM_prec) delete maxM_prec;
//    if (peaked_default_model_prec) delete peaked_default_model_prec;
    if (def_model_file_prec) delete def_model_file_prec;
    if (init_spectr_func_file_prec) delete init_spectr_func_file_prec;
//    if (interp_type_prec) delete interp_type_prec;
	if (output_dir_prec) delete output_dir_prec;
    if (w_sample_prec) delete w_sample_prec;
    if (alpha_init_prec) delete alpha_init_prec;
	if (default_model_center_prec) delete default_model_center_prec;
	if (default_model_width_prec) delete default_model_width_prec;
	if (default_model_shape_prec) delete default_model_shape_prec;
	if (A_ref_file_prec) delete A_ref_file_prec;
	 
    if (time_params_file) delete time_params_file;
    if (time_other_params_file) delete time_other_params_file;
}


bool OmegaMaxEnt_data::create_default_input_params_file()
{
	initialize=true;
	
    ofstream file(input_params_file_name);
    if (file)
    {
        file<<"data file:\n\nOPTIONAL PREPROCESSING TIME PARAMETERS\n\nDATA PARAMETERS\n";
        for (auto item: Data_params)
        {
            file<<item.second<<'\n';
        }
        file<<"\nINPUT FILES PARAMETERS\n";
        for (auto item: Input_files_params)
        {
            file<<item.second<<'\n';
        }
        file<<"\nFREQUENCY GRID PARAMETERS\n";
        for (auto item: Grid_params)
        {
            file<<item.second<<'\n';
        }
        file<<"\nCOMPUTATION OPTIONS\n";
        for (auto item: Preproc_comp_params)
        {
            file<<item.second<<'\n';
        }
        file<<"\nPREPROCESSING EXECUTION OPTIONS\n";
        for (auto item: Preproc_exec_params)
        {
            file<<item.second<<'\n';
        }
        file<<"\n\nOPTIONAL MINIMIZATION TIME PARAMETERS\n\nOUTPUT FILES PARAMETERS\n";
        for (auto item: Output_files_params)
        {
            file<<item.second<<'\n';
        }
        file<<"\nCOMPUTATION PARAMETERS\n";
        for (auto item: Optim_comp_params)
        {
            file<<item.second<<'\n';
        }
        file<<"\nMINIMIZATION EXECUTION OPTIONS\n";
        for (auto item: Optim_exec_params)
        {
            file<<item.second<<'\n';
        }
        file<<"\nDISPLAY OPTIONS\n";
        for (auto item: Optim_displ_params)
        {
            file<<item.second<<'\n';
        }
        file.close();
		
		cout<<"The file "<<input_params_file_name<<" was created.\n";
		cout<<"Fill it according to the instructions given in the User Guide before restarting the program.\n";
    }
    else
    {
        cout<<"cannot create file "<<input_params_file_name<<'\n';
        return false;
    }
    
    return true;
}



bool polyval(double x0, vec cfs, vec x, vec &y)
{
	int D=cfs.n_rows-1;
	vec Dx=x-x0;
	y=cfs(0)*pow(Dx,D);
	for (int j=1; j<D-1; j++)
		y=y+cfs(j)*pow(Dx,D-j);
	
	y=y+cfs(D-1)*Dx+cfs(D);
	
	return y.is_finite();
}

bool polyfit(vec x, vec y, int D, double x0, vec &cfs)
{
	int N=x.n_rows;
	
	mat X(N,D+1);
	X.zeros();
	
	int j;
	
	vec Dx=x-x0;
	for (j=0; j<D-1; j++)
	{
		X.col(j)=pow(Dx,D-j);
	}
	X.col(D-1)=Dx;
	X.col(D)=ones<vec>(N);
	
	mat A=X.t()*X;
	vec B=X.t()*y;
	
	cfs=solve(A,B);
	
	return cfs.is_finite();
}

void pascal(int n, imat &P)
{
	P.zeros(n,n);
	P.col(0)=ones<ivec>(n);
	P.row(0)=ones<irowvec>(n);
	int j,l;
	for (j=1; j<n; j++)
		for (l=1; l<n; l++)
			P(j,l)=P(j,l-1)+P(j-1,l);
}

void remove_spaces_front(string &str)
{
    int j=0;
    while (str[j]==' ' || str[j]=='\t') j++;
    str=str.substr(j);
}

void remove_spaces_back(string &str)
{
    int j=str.size()-1;
    while (str[j]==' ' || str[j]=='\t') j--;
    str.resize(j+1);
}

void remove_spaces_ends(string &str)
{
	remove_spaces_front(str);
	remove_spaces_back(str);
}



