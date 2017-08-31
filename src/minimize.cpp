
#include "maxEnt_data.h"

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
		
		
		/*
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
		*/
		
		double tmp_mean_int_dA=1e9; //arbitrairy large number (to ensure one while iteration at least)
		double tmp_mean_int_dA_prec=1e10; //arbitrairy large number (to ensure one while iteration at least)
		
		iter_dA=1;
		while ( (tmp_mean_int_dA>tol_int_dA || A1min<0) && iter_dA<Niter_dA_max && (tmp_mean_int_dA<tmp_mean_int_dA_prec))
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
			
			mean_int_dA=abs(dA1.t())*dwS;  // MC: this is the integral of the difference dA1
			tmp_mean_int_dA_prec=tmp_mean_int_dA;
			tmp_mean_int_dA=mean_int_dA(0);
			
			printf("% 4.8f  % 4.8f \n", tmp_mean_int_dA, mean_int_dA_prec);
			
			
			//if (mean_int_dA(0)<mean_int_dA_prec) //no need for this test, do the correction either ways
			//{
		   A1=A1+dA1;
			A1min=min(A1-Amin);
			//dA1.t().print();
			//}
			
			ind_Anul=find(A1==0);
			if (ind_Anul.n_rows)
			{
				A1.rows(ind_Anul)=Amin.rows(ind_Anul);
			}
			
			iter_dA++;
		}
		

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
			vec logA=log(A/default_model);
			ind_An=find(A/default_model<rADchange);
			if (ind_An.n_rows)
			{
				logA.rows(ind_An)=log(rADmin)*ones<vec>(ind_An.n_rows);
			}
			//Aw_samp.row(ind_alpha_vec)=trans(A(w_sample_ind));
			
			G_out=K*A;
			G_V_out=KG_V*A;

			int bosonOffset = 0;
			if (!boson || col_Gi>0)
			{
         	bosonOffset=1;			
				if (cov_diag)
				{
					uvec even_ind=linspace<uvec>(0,2*Nn-2,Nn);
					errRe=G_V.rows(even_ind)-G_V_out.rows(even_ind);
					errIm=G_V.rows(even_ind+1)-G_V_out.rows(even_ind+1);
				}
				else
				{
					errRe=G_V-G_V_out;
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

			gpc_plot_image(handle_gnuplot, vectors_A, w, ind_alpha-1, bosonOffset);
		   
			pow_alpha=pow_alpha-pow_alpha_step;
			alpha_prec=alpha;
			alpha=pow(10,pow_alpha);
		}
		else
		{
			cout<<"minimize(): change in chi^2 is in the wrong direction\n"; //mc: wtf?
			break;
		}
		
		/*
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
		*/
		
		ind_alpha++;
		ind_alpha_vec++;
	}
	
	if (ind_alpha>Nalpha)
		cout<<"maximum number of values alpha reached\n";
	else if (alpha<alpha_min)
		cout<<"minimum value of alpha reached\n";
}

