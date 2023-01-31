/*********************************************************************
We assume that the eigenstrain is dilatational and the elastic constants 
are cubic, and hence shortened the lengthy summations, including only 
the non-zero terms 
**********************************************************************/

#include<stdio.h>
#include<stdlib.h>

#include<complex.h>
#include<fftw3.h>

#include "../headers/functions.h"
#include "../headers/nrutil.h"

void calculate_hom_strain(int n_x, int n_y, int n_z, 
double delta_x, double delta_y, double delta_z, double Ceff[4][4][4][4], 
double DeltaC[4][4][4][4], double S[4][4][4][4], 
double **sig_zero_exp, fftw_complex *comp, double ave_comp,
fftw_complex *eps_star11, fftw_complex *eps_star12,
fftw_complex *eps_star22, fftw_complex *eps_star13, fftw_complex *eps_star33,
fftw_complex *eps_star23, double **sig_app, double **E, 
int n_alpha, double *alpha){

int i1,i2,i3;
int J;

double A11=0.0;
double A22=0.0;
double A33=0.0;
double A12=0.0;
double A21=0.0;
double A13=0.0;
double A31=0.0;
double A23=0.0;
double A32=0.0;
double epsstr11,epsstr22,epsstr33,epsstr12,epsstr13,epsstr23;

double a;
double inv_nxnynz;

inv_nxnynz = 1.0/(n_x*n_y*n_z);

/* Calculate the expectation value for the periodic strain */

for(i1=0; i1<n_x; ++i1){
for(i2=0; i2<n_y; ++i2){
for(i3=0; i3<n_z; ++i3){
J = i3 + n_z*(i2+n_y*i1);
epsstr11 = creal(eps_star11[J]);
epsstr22 = creal(eps_star22[J]);
epsstr33 = creal(eps_star33[J]);
epsstr12 = creal(eps_star12[J]);
epsstr13 = creal(eps_star13[J]);
epsstr23 = creal(eps_star23[J]);
a = alpha[(int) (n_alpha*creal(comp[J]))];

A11 = A11 
	+ (Ceff[1][1][1][1]+DeltaC[1][1][1][1]*a)*epsstr11
	+ (Ceff[1][1][2][2]+DeltaC[1][1][2][2]*a)*epsstr22
	+ (Ceff[1][1][3][3]+DeltaC[1][1][3][3]*a)*epsstr33
	+ (Ceff[1][1][1][2]+DeltaC[1][1][1][2]*a)*epsstr12
	+ (Ceff[1][1][2][1]+DeltaC[1][1][2][1]*a)*epsstr12
	+ (Ceff[1][1][1][3]+DeltaC[1][1][1][3]*a)*epsstr13
	+ (Ceff[1][1][3][1]+DeltaC[1][1][3][1]*a)*epsstr13
	+ (Ceff[1][1][2][3]+DeltaC[1][1][2][3]*a)*epsstr23
	+ (Ceff[1][1][3][2]+DeltaC[1][1][3][2]*a)*epsstr23;

A22 = A22 
	+ (Ceff[2][2][1][1]+DeltaC[2][2][1][1]*a)*epsstr11
	+ (Ceff[2][2][2][2]+DeltaC[2][2][2][2]*a)*epsstr22 
	+ (Ceff[2][2][3][3]+DeltaC[2][2][3][3]*a)*epsstr33 
	+ (Ceff[2][2][1][2]+DeltaC[2][2][1][2]*a)*epsstr12
	+ (Ceff[2][2][2][1]+DeltaC[2][2][2][1]*a)*epsstr12
	+ (Ceff[2][2][1][3]+DeltaC[2][2][1][3]*a)*epsstr13
	+ (Ceff[2][2][3][1]+DeltaC[2][2][3][1]*a)*epsstr13
	+ (Ceff[2][2][2][3]+DeltaC[2][2][2][3]*a)*epsstr23
	+ (Ceff[2][2][3][2]+DeltaC[2][2][3][2]*a)*epsstr23;

A33 = A33 
	+ (Ceff[3][3][1][1]+DeltaC[3][3][1][1]*a)*epsstr11
	+ (Ceff[3][3][2][2]+DeltaC[3][3][2][2]*a)*epsstr22 
	+ (Ceff[3][3][3][3]+DeltaC[3][3][3][3]*a)*epsstr33 
	+ (Ceff[3][3][1][2]+DeltaC[3][3][1][2]*a)*epsstr12
	+ (Ceff[3][3][2][1]+DeltaC[3][3][2][1]*a)*epsstr12
	+ (Ceff[3][3][1][3]+DeltaC[3][3][1][3]*a)*epsstr13
	+ (Ceff[3][3][3][1]+DeltaC[3][3][3][1]*a)*epsstr13
	+ (Ceff[3][3][2][3]+DeltaC[3][3][2][3]*a)*epsstr23
	+ (Ceff[3][3][3][2]+DeltaC[3][3][3][2]*a)*epsstr23;

A12 = A12 
	+ (Ceff[1][2][1][2]+DeltaC[1][2][1][2]*a
	+ Ceff[1][2][2][1]+DeltaC[1][2][2][1]*a)*epsstr12;

A13 = A13 
	+ (Ceff[1][3][1][3]+DeltaC[1][3][1][3]*a
	+ Ceff[1][3][3][1]+DeltaC[1][3][3][1]*a)*epsstr13;

A23 = A23 
	+ (Ceff[2][3][2][3]+DeltaC[2][3][2][3]*a
	+ Ceff[2][3][3][2]+DeltaC[2][3][3][2]*a)*epsstr23;

}}}

/* Calculate the 'effective' stress */

A11 = sig_zero_exp[1][1] + sig_app[1][1] - A11*inv_nxnynz;
A22 = sig_zero_exp[2][2] + sig_app[2][2] - A22*inv_nxnynz;
A33 = sig_zero_exp[3][3] + sig_app[3][3] - A33*inv_nxnynz;
A12 = sig_app[1][2] - A12*inv_nxnynz;
A13 = sig_app[1][3] - A13*inv_nxnynz;
A23 = sig_app[2][3] - A23*inv_nxnynz;
A21 = A12;
A31 = A13;
A32 = A23;

/* Multiply the 'effective' stress by the 'effective' compliance to
 * obtain the homogeneous strain */

E[1][1] = S[1][1][1][1]*A11 + S[1][1][2][2]*A22 + S[1][1][3][3]*A33;
E[2][2] = S[2][2][1][1]*A11 + S[2][2][2][2]*A22 + S[2][2][3][3]*A33;
E[3][3] = S[3][3][1][1]*A11 + S[3][3][2][2]*A22 + S[3][3][3][3]*A33;
E[1][2] = (S[1][2][1][2] + S[1][2][2][1])*A12;
E[2][1] = E[1][2];
E[1][3] = (S[1][3][1][3] + S[1][3][3][1])*A13;
E[3][1] = E[1][3];
E[2][3] = (S[2][3][2][3] + S[2][3][3][2])*A23;
E[3][2] = E[2][3];

}
