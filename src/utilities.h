/*
 file utilities.h
 deinition of class "utilities" containing multipurpose numerical routines. The best way to use it is to derive a class from it.
 */
#pragma once


//start of includeDef.h
#ifdef USE_MPI
#include "fftw3-mpi.h"
#include "mpi.h"
#endif

//#include "fftw3.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdio>
#include <complex>
#include <string>
#include <cstring>
#include <ctime>
#include <limits>
#include <stdio.h>
#include <sys/stat.h>

#ifndef PI
#define PI acos((double)-1.0)
#endif

#define EPSILON numeric_limits<double>::epsilon()
#define DBL_MIN numeric_limits<double>::min()
#define DBL_MAX numeric_limits<double>::max()
#define INF numeric_limits<double>::infinity()

using namespace std;

typedef complex<double> dcomplex;
typedef unsigned int  uint;
//end of includeDef.h



#include "armadillo"

using namespace arma;


extern "C"
{
	// routines LAPACK  (descriptions sur http://www.netlib.org/lapack )
	// resout le systeme d'equations AX=B ou A est tridiagonal, reels double precision 	
	void dgtsv_(int *N, int *NRHS, double *DL, double *D, double *DU, double *B, int *LDB, int *INFO );
	// resout le systeme d'equations AX=B ou A est tridiagonal, complex double precision 	
	void zgtsv_(int *N, int *NRHS, dcomplex *DL, dcomplex *D, dcomplex *DU, dcomplex *B, int *LDB, int *INFO );
	// resout le systeme d'equations AX=B ou A a plusieurs diagonales, reals double precision
	void dgbsv_(int *N, int *KL, int *KU, int *NRHS, double *AB, int *LDAB, int *IPIV, double *B, int *LDB, int *INFO );
	// routine avancee pour resoudre le systeme d'equations AX=B ou A a plusieurs diagonales, reals double precision
	void dgbsvx_( char *FACT, char *TRANS, int *N, int *KL, int *KU, int *NRHS, double *AB, int *LDAB, double *AFB,
				 int *LDAFB, int *IPIV, char *EQUED, double *R, double *C, double *B, int *LDB, double *X, int *LDX,
				 double *RCOND, double *FERR, double *BERR, double *WORK, int *IWORK, int *INFO );
	// resout le systeme d'equations AX=B ou A a plusieurs diagonales, complex double precision
	void zgbsv_(int *N, int *KL, int *KU, int *NRHS, dcomplex *AB, int *LDAB, int *IPIV, dcomplex *B, int *LDB, int *INFO );

	//computes the solution to system of linear equations A * X = B ( http://www.netlib.org/lapack/explore-html/d8/d72/dgesv_8f_source.html )
	void dgesv_(int *N, int *NRHS, double *A, int *LDA, int *IPIV, double *B, int *LDB, int *INFO );
	// computes the solution to a real system of linear equations A * X = B, where A is an N-by-N symmetric positive definite matrix and X and B are N-by-NRHS matrices. ( http://www.netlib.org/lapack/explore-html/d9/d6f/dposv_8f_source.html )
	void dposv_(char *UPLO, int *N, int *NRHS, double *A, int *LDA, double *B, int *LDB, int *INFO );
	//computes the singular value decomposition ( http://www.netlib.org/lapack/explore-html/d8/d2d/dgesvd_8f_source.html )
	void dgesvd_(char *JOBU, char *JOBVT, int *M, int *N, double *A, int *LDA, double *S, double *U, int *LDU, double *VT, int *LDVT, double *WORK, int *LWORK, int *INFO );
	//computes the singular value decomposition. If singular vectors are desired, uses a divide and conquer algorithm. ( http://www.netlib.org/lapack/explore-html/db/db4/dgesdd_8f_source.html )
	void dgesdd_( char *JOBZ, int *M, int *N, double *A, int *LDA, double *S, double *U, int *LDU, double *VT, int *LDVT, double *WORK, int *LWORK, int *IWORK, int *INFO );
	//computes the eigenvalues and, optionally, the left and/or right eigenvectors for SY matrices ( http://www.netlib.org/lapack/explore-html/dd/d4c/dsyev_8f_source.html )
	void dsyev_( char *JOBZ, char *UPLO, int *N, double *A, int *LDA, double *W, double *WORK, int *LWORK, int *INFO );
	//solves a general Gauss-Markov linear model (GLM) problem. ( http://www.netlib.org/lapack/explore-html/d3/df4/dggglm_8f_source.html )
	void dggglm_(int *N, int *M, int *P, double *A, int *LDA, double *B, int *LDB, double *D, double *X, double *Y, double *WORK, int *LWORK, int *INFO );
	//computes the solution to system of linear equations A * X = B, where A is a band matrix ( http://www.netlib.org/lapack/explore-html/dd/dc2/dgbsv_8f_source.html )
	void dgbsv_(int *N, int *KL, int *KU, int *NRHS, double *AB, int *LDAB, int *IPIV, double *B, int *LDB, int *INFO );
}



extern "C++" {
	
//! Note: to use the routines that take a pointer to a function as a parameter, use the static_cast() function to cast the pointer to a function of your class as a function of this class. Ex: fctPtr1 Ptr=static_cast<fctPtr1> (&myclass::myfunc);
class utilities
{
 public:
	typedef double (utilities::*fctPtr1) (double, void*[]);
	typedef dcomplex (utilities::*cx_fctPtr1) (double, void*[]);
	
	typedef double (utilities::*fctPtr) (double, double[]);
	typedef double (utilities::*IntPtr)(double, double, int[], double[], fctPtr, double[2]);
	
	typedef dcomplex (utilities::*complexFctPtr) (double, double[]);
	typedef dcomplex (utilities::*complexIntPtr)(double, double, int[], double[], complexFctPtr, double[2]);
	
						
//@{
//! Constructor
	utilities(){};
//! Destructor
	~utilities(){};
	
//! root finding routine	
	bool find_zero(fctPtr func, double init[], double params[], double root[]);
	
//! root finding routine	
	bool find_zero(fctPtr func, double init[], double params[], double root[], double lims[]);	
		
//@{ Adaptive quadratic 1D and 2D integration routines for real functions
//! 1D integation routine (also used by 2D routines). In direct use for 1D functions, last two arguments should be NULL
	double quadInteg(fctPtr, double lims[2], double tol, int nbEval[], double params[], IntPtr, double limx[2]);
//! 1D subinterval integration step
	double quadStep(fctPtr, double lims[2], double vals[3], double tol, int nbEval[], double hmin, double params[], 
					IntPtr, double limx[2]);
//! 2D integration
	double quadInteg2D(fctPtr, double limx[2], double limy[2], double tol, int nbEval[], double params[]);
//! Inner (x directions) integral in 2D routine
	double innerInteg(double y, double tol, int nbEval[], double params[], fctPtr, double limx[2]);
	
	//! complex 1D integration routine whith a pointer array as parameter
	dcomplex cx_quadInteg1D(cx_fctPtr1, double lims[2], double tol, int nbEval[], void *params[]);
	//! 1D subinterval complex integration step
	dcomplex cx_quadStep1D(cx_fctPtr1, double lims[2], dcomplex vals[3], double tol, int nbEval[], double hmin, void *params[]);
	
//! another 1D integration routine whith a pointer array as parameter
	double quadInteg1D(fctPtr1, double lims[2], double tol, int nbEval[], void *params[]);
//! 1D subinterval integration step
	double quadStep1D(fctPtr1, double lims[2], double vals[3], double tol, int nbEval[], double hmin, void *params[]);
//@}
	
	//! compute the coefficients of a cubic spline using linear combinations of the first and second derivatives as the two mandatory additionnal constraints to the spline. x0 is the position vector, F the function vector, LC contains the two right-hand side values of the constraints, coeffs_LC(0) and coeffs_LC(1) are the coefficients of the first derivatives at the left and right boundary and coeffs_LC(2) and coeffs_LC(3) are the coefficients of the second derivatives at the same boundaries.
	void spline_coeffs_LC(vec x0, vec F, vec LC, vec coeffs_LC, vec &coeffs);
	
	//! compute the coefficients of a clamped spline.
	void spline_coefficients(double *vec_coeff,const double *xx,const double *yy,const double *fp,const int size);
	
	//! evaluate the spline at x	
	double Interp_1d_spline(const double x,const double *coeff,const double *xx,const int size);
	
	//! locate the index of the interval where is located
	int locate(const double x,const double *xx,const int size);
	
	//! Integrate the spline
	double oned_spline_integral(const double *coeff,const double *xx,const double *yy,const int size);

	//! compute the coefficients of a cubic spline for V(x) of the form S_j(x)=a_j*(x-x_j)^3+b_j*(x-x_j)^2+c_j*(x-x_j)+d_j, where j is the interval's index
	//coeffs[0] and coeffs[1] contains the derivatives of V(x) at x_1 et x_N0, the size of x0 and V is N0 and the size of coeffs is 4*(N0-1)
	void spline_coeffs_rel(double *x0, double *V, int N0, double *coeffs);
	
	//! compute the coefficients of a cubic spline for V(x) of the form S_j(x)=a_j*x^3+b_j*x^2+c_j*x+d_j, where j is the interval's index
	//coeffs[0] and coeffs[1] contains the derivatives of V(x) at x_1 et x_N0, the size of x0 and V is N0 and the size of coeffs is 4*(N0-1)
	void spline_coeffs(double *x0, double *V, int N0, double *coeffs);
	
	//! return the value at x of the cubic spline computed with spline_coeffs()
	double spline(double x, double *x0, double *V, int N0, double *coeffs);

	//! pade approximant for an input value x after coef is computed with pade_cont_frac_coef (or pade_cont_frac_coef_rec),
	//! n: number of input points in the Pade approximant, x0: vector of the input points
	dcomplex pade(dcomplex x, int n, dcomplex *x0, dcomplex *coef);
	
	//! recursive formula in the Pad? coefficients calculation
	dcomplex pade_recursion(int indFunc, int indx, dcomplex *func, dcomplex *x, dcomplex *coef);
	
	//! compute recursively the coefficients in the continued fraction representation of the Pade aproximant,
	//! func: vector of values of the function, x: vector of the input points, N: size of 'func' and 'x',
	//! coefficients are returned in 'coef'
	void pade_cont_frac_coef_rec(dcomplex *func, dcomplex *x, int N, dcomplex *coef) { for (int j=0; j<N; j++) coef[j]=pade_recursion(j, j, func, x, coef); }
	
	//! compute the coefficients in the continued fraction representation of the Pade aproximant,
	// func: vector of values of the function, x: vector of the input points, N: size of 'func' and 'x',
	// coefficients are returned in 'coef'
	inline void pade_cont_frac_coef(dcomplex *func, dcomplex *x, int N, dcomplex *coef);
	
	//! find the minimum of a function
	bool find_min(fctPtr func, double init[], double params[], double min[]);	
	
	//! absolute value for complex number (modulus)
	double fabs(dcomplex a){ return sqrt( real(a)*real(a) + imag(a)*imag(a) ); }
	
	//! find the roots of the polynomial Ax^2+Bx+C with high precision
	void find_roots_2pol(dcomplex coef[], dcomplex roots[]);
	
	//! find the roots of the polynomial Ax^3+Bx^2+Cx+D
	void find_roots_3pol(dcomplex coef[], dcomplex roots[]);
	
	//! find the roots of the polynomial Ax^4+Bx^3+Cx^2+Dx+F
	void find_roots_4pol(dcomplex coef[], dcomplex roots[]);

	//! test if a vector contains NaN
	bool contains_NaN(double *v, long int size)
	{
		bool NaN_present=false;
		for (long int j=0; j<size; j++)
			if (v[j]!=v[j]) 
			{
				NaN_present=true;
				break;
			}
		return NaN_present;
	}
	
	bool copy_file(const string fname, const string in_dir, const string out_dir);
	
		
};
	
} /* extern "C++" */


			
