/*********************************************************************
In an attempt to avoid unnecessary algebra, we have simplified the
summations to include only the non-zero terms while assuming that the 
elastic constants are cubic and the eigenstrain is dilatational
**********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include <fftw3.h>

#include "../headers/functions.h"
#include "../headers/nrutil.h"

void refine_u(int n_x, int n_y, int n_z, int half_nx, int half_ny, int half_nz,
double delta_x, double delta_y, double delta_z, double **epsilon_T,
double delta_kx, double delta_ky, double delta_kz, fftw_complex *comp,
double ave_comp, double Ceff[4][4][4][4], double DeltaC[4][4][4][4], 
double **sigma_zero, fftw_complex *u1_old, fftw_complex *u2_old, fftw_complex *u3_old, 
fftw_complex *u1_new, fftw_complex *u2_new, fftw_complex *u3_new, double **Del_sigma_T, 
double **sig_app, double **E, int MAXITS, double MAXERR,
int n_alpha, int n_beta, double *alpha, double *beta,
double *Omega11, double *Omega12, double *Omega21,
 double *Omega22, double *Omega13, double *Omega31, double *Omega33, double *Omega23, double *Omega32,
fftw_plan planF, fftw_plan planB){

FILE *fpw;
FILE *fp1;

#ifndef DEBUG
FILE *fp2, *fp3;
#endif

#ifdef DETAILED_DEBUG
int i5,i6;
double IDENT[4][4][4][4];
double C[4][4][4][4];
double sum1;
#endif

int i1,i2,i3,i4,i5;

int ITS=0;
double ERR=1.0;
double tmp;

double *kx, *ky, *kz;

double S[4][4][4][4];
double **sig_zero_exp;
double zeta11, zeta12, zeta21, zeta22, zeta13, zeta31, zeta33, zeta23, zeta32;

int nxnynz, J;
double inv_nxnynz;

sig_zero_exp = dmatrix(1,3,1,3);

fftw_complex gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7, gamma8;

fftw_complex *u1_tmp, *u2_tmp, *u3_tmp;
fftw_complex *chi11, *chi12, *chi21, *chi22, *chi13, *chi31, *chi33, *chi23, *chi32;
fftw_complex *Alpha, *Beta, *AlpBeta;
fftw_complex *eps_star11, *eps_star22, *eps_star33, *eps_star12, *eps_star13, *eps_star23;

double a;

/* Allocate memory */

nxnynz = n_x*n_y*n_z;
inv_nxnynz = 1.0/nxnynz;
u1_tmp = fftw_malloc(nxnynz* sizeof(fftw_complex));
u2_tmp = fftw_malloc(nxnynz* sizeof(fftw_complex));
u3_tmp = fftw_malloc(nxnynz* sizeof(fftw_complex));
chi11 = fftw_malloc(nxnynz* sizeof(fftw_complex));
chi12 = fftw_malloc(nxnynz* sizeof(fftw_complex));
chi21 = fftw_malloc(nxnynz* sizeof(fftw_complex));
chi22 = fftw_malloc(nxnynz* sizeof(fftw_complex));
chi13 = fftw_malloc(nxnynz* sizeof(fftw_complex));
chi31 = fftw_malloc(nxnynz* sizeof(fftw_complex));
chi33 = fftw_malloc(nxnynz* sizeof(fftw_complex));
chi23 = fftw_malloc(nxnynz* sizeof(fftw_complex));
chi32 = fftw_malloc(nxnynz* sizeof(fftw_complex));
Alpha = fftw_malloc(nxnynz* sizeof(fftw_complex));
Beta = fftw_malloc(nxnynz* sizeof(fftw_complex));
AlpBeta = fftw_malloc(nxnynz* sizeof(fftw_complex));
eps_star11 = fftw_malloc(nxnynz* sizeof(fftw_complex));
eps_star22 = fftw_malloc(nxnynz* sizeof(fftw_complex));
eps_star33 = fftw_malloc(nxnynz* sizeof(fftw_complex));
eps_star12 = fftw_malloc(nxnynz* sizeof(fftw_complex));
eps_star13 = fftw_malloc(nxnynz* sizeof(fftw_complex));
eps_star23 = fftw_malloc(nxnynz* sizeof(fftw_complex));

/* Generate the array of reciprocal vector components */

kx = (double *) malloc((size_t) n_x* sizeof(double));
ky = (double *) malloc((size_t) n_y* sizeof(double));
kz = (double *) malloc((size_t) n_z* sizeof(double));

for(i1=0; i1 < n_x; ++i1){
if(i1 < half_nx) kx[i1] = i1*delta_kx;
else kx[i1] = (i1-n_x)*delta_kx;
}
for(i2=0; i2 < n_y; ++i2){
if(i2 < half_ny) ky[i2] = i2*delta_ky;
else ky[i2] = (i2-n_y)*delta_ky;
}
for(i3=0; i3 < n_z; ++i3){
if(i3 < half_nz) kz[i3] = i3*delta_kz;
else kz[i3] = (i3-n_z)*delta_kz;
}

/* Calculate the derivative of displacement */

for(i1=0; i1 < n_x; ++i1){
for(i2=0; i2 < n_y; ++i2){
for(i3=0; i3 < n_z; ++i3){
	J = i3 + n_z*(i2+n_y*i1);
  	chi11[J] = _Complex_I*kx[i1]*u1_old[J];
	chi22[J] = _Complex_I*ky[i2]*u2_old[J];
	chi33[J] = _Complex_I*kz[i3]*u3_old[J];

	chi12[J] = _Complex_I*ky[i2]*u1_old[J];
	chi21[J] = _Complex_I*kx[i1]*u2_old[J];

	chi13[J] = _Complex_I*kz[i3]*u1_old[J];
	chi31[J] = _Complex_I*kx[i1]*u3_old[J];

	chi23[J] = _Complex_I*kz[i3]*u2_old[J];
	chi32[J] = _Complex_I*ky[i2]*u3_old[J];

	eps_star12[J] = 0.5*(chi12[J]+chi21[J]);
	eps_star13[J] = 0.5*(chi13[J]+chi31[J]);
	eps_star23[J] = 0.5*(chi23[J]+chi32[J]);
}}}

fftw_execute_dft(planB,chi11,chi11);
fftw_execute_dft(planB,chi22,chi22);
fftw_execute_dft(planB,chi33,chi33);
fftw_execute_dft(planB,chi12,chi12);
fftw_execute_dft(planB,chi21,chi21);
fftw_execute_dft(planB,chi13,chi13);
fftw_execute_dft(planB,chi31,chi31);
fftw_execute_dft(planB,chi23,chi23);
fftw_execute_dft(planB,chi32,chi32);

fftw_execute_dft(planB,eps_star12,eps_star12);
fftw_execute_dft(planB,eps_star13,eps_star13);
fftw_execute_dft(planB,eps_star23,eps_star23);

/* Multiply the derivative of displacement by alpha */

for(i1=0; i1<n_x; ++i1){
for(i2=0; i2<n_y; ++i2){
for(i3=0; i3<n_z; ++i3){

	J = i3 + n_z*(i2+n_y*i1);
	chi11[J] = chi11[J]*inv_nxnynz;
	chi12[J] = chi12[J]*inv_nxnynz;
	chi21[J] = chi21[J]*inv_nxnynz;
	chi22[J] = chi22[J]*inv_nxnynz;
	chi13[J] = chi13[J]*inv_nxnynz;
	chi31[J] = chi31[J]*inv_nxnynz;
	chi33[J] = chi33[J]*inv_nxnynz;
	chi23[J] = chi23[J]*inv_nxnynz;
	chi32[J] = chi32[J]*inv_nxnynz;

	eps_star11[J] = chi11[J];
	eps_star22[J] = chi22[J];
	eps_star33[J] = chi33[J];
	eps_star12[J] = eps_star12[J]*inv_nxnynz;
	eps_star13[J] = eps_star13[J]*inv_nxnynz;
	eps_star23[J] = eps_star23[J]*inv_nxnynz;


	a = alpha[(int) (n_alpha*creal(comp[J]))];


	chi11[J] = a*chi11[J];
	chi12[J] = a*chi12[J];
	chi21[J] = a*chi21[J];
	chi22[J] = a*chi22[J];
	chi13[J] = a*chi13[J];
	chi31[J] = a*chi31[J];
	chi33[J] = a*chi33[J];
	chi23[J] = a*chi23[J];
	chi32[J] = a*chi32[J];
}}}

fftw_execute_dft(planF,chi11,chi11);
fftw_execute_dft(planF,chi12,chi12);
fftw_execute_dft(planF,chi21,chi21);
fftw_execute_dft(planF,chi22,chi22);
fftw_execute_dft(planF,chi13,chi13);
fftw_execute_dft(planF,chi31,chi31);
fftw_execute_dft(planF,chi33,chi33);
fftw_execute_dft(planF,chi23,chi23);
fftw_execute_dft(planF,chi32,chi32);

/* Calculate the displacements, alpha, beta and alpha*beta */

calculate_S_exp(n_x,n_y,n_z,delta_x,delta_y,delta_z,comp,ave_comp,Ceff,DeltaC,S,
n_alpha,alpha);

#ifdef DETAILED_DEBUG

for(i1=1; i1<4; ++i1){
for(i2=1; i2<4; ++i2){
for(i3=1; i3<4; ++i3){
for(i4=1; i4<4; ++i4){
IDENT[i1][i2][i3][i4]=0.0;
}}}}

sum1 = 0.0;
for(i1=0; i1<n_x; ++i1){
for(i2=0; i2<n_y; ++i2){
for(i3=0; i3<n_z; ++i3){
J = i3 + n_z*(i2+n_y*i1);
sum1 = sum1 + alpha[(int) (n_alpha*creal(comp[J]))];
}}}
sum1 = sum1/(n_x*n_y*n_z);

C[1][1][1][1] = Ceff[1][1][1][1] + DeltaC[1][1][1][1]*sum1;
C[1][2][1][2] = Ceff[1][2][1][2] + DeltaC[1][2][1][2]*sum1;
C[1][1][2][2] = Ceff[1][1][2][2] + DeltaC[1][1][2][2]*sum1;
/* c11 */
C[2][2][2][2] = C[1][1][1][1];
C[3][3][3][3] = C[1][1][1][1];
/* c12 */
C[2][2][1][1] = C[1][1][2][2];
C[1][1][3][3] = C[3][3][1][1] = C[1][1][2][2];
C[2][2][3][3] = C[3][3][2][2] = C[1][1][2][2];
/* c44 */
C[1][2][2][1] = C[2][1][1][2] = C[2][1][2][1]  =
C[1][2][1][2];
C[1][3][1][3] = C[1][3][3][1] = C[3][1][1][3] = C[3][1][3][1]  =
C[1][2][1][2];
C[2][3][2][3] = C[2][3][3][2] = C[3][2][2][3] = C[3][2][3][2]  =
C[1][2][1][2];

for(i1=1; i1<4; ++i1){
for(i2=1; i2<4; ++i2){
for(i3=1; i3<4; ++i3){
for(i4=1; i4<4; ++i4){
for(i5=1; i5<4; ++i5){
for(i6=1; i6<4; ++i6){
IDENT[i1][i2][i3][i4] = IDENT[i1][i2][i3][i4] 
+ C[i1][i2][i5][i6]*S[i5][i6][i3][i4];
}}}}}}

for(i1=1; i1<4; ++i1){
for(i2=1; i2<4; ++i2){
for(i3=1; i3<4; ++i3){
for(i4=1; i4<4; ++i4){
printf("%d %d %d %d %le %le %le\n",i1,i2,i3,i4,
C[i1][i2][i3][i4],S[i1][i2][i3][i4],IDENT[i1][i2][i3][i4]);
}}}}

#endif

calculate_sig_zero_exp(n_x,n_y,n_z,delta_x,delta_y,delta_z,
comp,ave_comp,Ceff,DeltaC,epsilon_T,sig_zero_exp,n_alpha,n_beta,
alpha,beta);

calculate_hom_strain(n_x,n_y,n_z,delta_x,delta_y,delta_z,
Ceff,DeltaC,S,sig_zero_exp,comp,ave_comp,eps_star11,eps_star12,
eps_star22,eps_star13,eps_star33,eps_star23,sig_app,E,n_alpha,alpha);

fftw_execute_dft(planB,u1_old,u1_old);
fftw_execute_dft(planB,u2_old,u2_old);
fftw_execute_dft(planB,u3_old,u3_old);

for(i1=0; i1<n_x; ++i1){
for(i2=0; i2<n_y; ++i2){
for(i3=0; i3<n_z; ++i3){
	J = i3 + n_z*(i2+n_y*i1);
	u1_old[J] = u1_old[J]*inv_nxnynz;
	u2_old[J] = u2_old[J]*inv_nxnynz;
	u3_old[J] = u3_old[J]*inv_nxnynz;
	
	u1_tmp[J] = u1_old[J];
	u2_tmp[J] = u2_old[J];
	u3_tmp[J] = u3_old[J];
	
	Alpha[J] = alpha[(int) (n_alpha*creal(comp[J]))];
	Beta[J] = beta[(int) (n_beta*creal(comp[J]))];
	AlpBeta[J] = Alpha[J]*Beta[J];
}}}
fftw_execute_dft(planF,Alpha,Alpha);
fftw_execute_dft(planF,Beta,Beta);
fftw_execute_dft(planF,AlpBeta,AlpBeta);

/* Calculte zeta = DeltaC*Hom_strain */

zeta11 =  DeltaC[1][1][1][1]*E[1][1]
	+ DeltaC[1][1][2][2]*E[2][2]
	+ DeltaC[1][1][3][3]*E[3][3];

zeta22 =  DeltaC[2][2][1][1]*E[1][1]
	+ DeltaC[2][2][2][2]*E[2][2]
	+ DeltaC[2][2][3][3]*E[3][3];

zeta33 =  DeltaC[3][3][1][1]*E[1][1]
	+ DeltaC[3][3][2][2]*E[2][2]
	+ DeltaC[3][3][3][3]*E[3][3];

zeta12 =  DeltaC[1][2][1][2]*E[1][2]
	+ DeltaC[1][2][2][1]*E[2][1];

zeta21 =  DeltaC[2][1][1][2]*E[1][2]
	+ DeltaC[2][1][2][1]*E[2][1];

zeta13 =  DeltaC[1][3][1][3]*E[1][3]
	+ DeltaC[1][3][3][1]*E[3][1];

zeta31 =  DeltaC[3][1][1][3]*E[1][3]
	+ DeltaC[3][1][3][1]*E[3][1];

zeta23 =  DeltaC[2][3][2][3]*E[2][3]
	+ DeltaC[2][3][3][2]*E[3][2];

zeta32 =  DeltaC[3][2][2][3]*E[2][3]
	+ DeltaC[3][2][3][2]*E[3][2];

/* Refine the solution */

while(ITS != MAXITS && ERR > MAXERR ){

	for(i1=0; i1 < n_x; ++i1){
	for(i2=0; i2 < n_y; ++i2){
	for(i3=0; i3 < n_z; ++i3){
	
		
		J = i3 + n_z*(i2+n_y*i1);

		gamma0 = DeltaC[1][1][1][1]*chi11[J]
							+ DeltaC[1][1][2][2]*chi22[J]
							+ DeltaC[1][1][3][3]*chi33[J];
		gamma1 = DeltaC[1][2][1][2]*chi12[J]
							+ DeltaC[1][2][2][1]*chi21[J];
		gamma2 = DeltaC[2][1][1][2]*chi12[J]
							+ DeltaC[2][1][2][1]*chi21[J];
		gamma3 = DeltaC[2][2][1][1]*chi11[J]
							+ DeltaC[2][2][2][2]*chi22[J]
							+ DeltaC[2][2][3][3]*chi33[J];

		gamma4 = DeltaC[1][3][1][3]*chi13[J]
							+ DeltaC[1][3][3][1]*chi31[J];
		gamma5 = DeltaC[3][1][1][3]*chi13[J]
							+ DeltaC[3][1][3][1]*chi31[J];
		gamma6 = DeltaC[3][3][1][1]*chi11[J]
							+ DeltaC[3][3][2][2]*chi22[J]
							+ DeltaC[3][3][3][3]*chi33[J];
		gamma7 = DeltaC[2][3][2][3]*chi23[J]
							+ DeltaC[2][3][3][2]*chi32[J];
		gamma8 = DeltaC[3][2][2][3]*chi23[J]
							+ DeltaC[3][2][3][2]*chi32[J];

		u1_new[J] =	- _Complex_I*Omega11[J]*kx[i1]
		*(sigma_zero[1][1]*Beta[J] - zeta11*Alpha[J] + Del_sigma_T[1][1]*AlpBeta[J] - gamma0)

				- _Complex_I*Omega21[J]*kx[i1]
		*(sigma_zero[2][1]*Beta[J] - zeta21*Alpha[J] + Del_sigma_T[2][1]*AlpBeta[J] - gamma2)

				- _Complex_I*Omega11[J]*ky[i2]
		*(sigma_zero[1][2]*Beta[J] - zeta12*Alpha[J] + Del_sigma_T[1][2]*AlpBeta[J] - gamma1)

				- _Complex_I*Omega21[J]*ky[i2]
		*(sigma_zero[2][2]*Beta[J] - zeta22*Alpha[J] + Del_sigma_T[2][2]*AlpBeta[J] - gamma3)

				- _Complex_I*Omega31[J]*kx[i1]
		*(sigma_zero[3][1]*Beta[J] - zeta31*Alpha[J] + Del_sigma_T[3][1]*AlpBeta[J] - gamma5)

				- _Complex_I*Omega11[J]*kz[i3]
		*(sigma_zero[1][3]*Beta[J] - zeta13*Alpha[J] + Del_sigma_T[1][3]*AlpBeta[J] - gamma4)

				- _Complex_I*Omega31[J]*kz[i3]
		*(sigma_zero[3][3]*Beta[J] - zeta33*Alpha[J] + Del_sigma_T[3][3]*AlpBeta[J] - gamma6);


		u2_new[J] =  	- _Complex_I*Omega12[J]*kx[i1]
		*(sigma_zero[1][1]*Beta[J] - zeta11*Alpha[J] + Del_sigma_T[1][1]*AlpBeta[J] - gamma0)
		
				- _Complex_I*Omega22[J]*kx[i1]
		*(sigma_zero[2][1]*Beta[J] - zeta21*Alpha[J] + Del_sigma_T[2][1]*AlpBeta[J] - gamma2)
		
				- _Complex_I*Omega12[J]*ky[i2]
		*(sigma_zero[1][2]*Beta[J] - zeta12*Alpha[J] + Del_sigma_T[1][2]*AlpBeta[J] - gamma1)
		
				- _Complex_I*Omega22[J]*ky[i2]
		*(sigma_zero[2][2]*Beta[J] - zeta22*Alpha[J] + Del_sigma_T[2][2]*AlpBeta[J] - gamma3)

				- _Complex_I*Omega32[J]*ky[i2]
		*(sigma_zero[3][2]*Beta[J] - zeta32*Alpha[J] + Del_sigma_T[3][2]*AlpBeta[J] - gamma8)
		
				- _Complex_I*Omega22[J]*kz[i3]
		*(sigma_zero[2][3]*Beta[J] - zeta23*Alpha[J] + Del_sigma_T[2][3]*AlpBeta[J] - gamma7)
		
				- _Complex_I*Omega32[J]*kz[i3]
		*(sigma_zero[3][3]*Beta[J] - zeta33*Alpha[J] + Del_sigma_T[3][3]*AlpBeta[J] - gamma6);


		u3_new[J] =  	- _Complex_I*Omega13[J]*kx[i1]
		*(sigma_zero[1][1]*Beta[J] - zeta11*Alpha[J] + Del_sigma_T[1][1]*AlpBeta[J] - gamma0)
		
				- _Complex_I*Omega33[J]*kx[i1]
		*(sigma_zero[3][1]*Beta[J] - zeta31*Alpha[J] + Del_sigma_T[3][1]*AlpBeta[J] - gamma5)
		
				- _Complex_I*Omega13[J]*kz[i3]
		*(sigma_zero[1][3]*Beta[J] - zeta13*Alpha[J] + Del_sigma_T[1][3]*AlpBeta[J] - gamma4)
		
				- _Complex_I*Omega33[J]*kz[i3]
		*(sigma_zero[3][3]*Beta[J] - zeta33*Alpha[J] + Del_sigma_T[3][3]*AlpBeta[J] - gamma6)

				- _Complex_I*Omega23[J]*kz[i3]
		*(sigma_zero[2][3]*Beta[J] - zeta23*Alpha[J] + Del_sigma_T[2][3]*AlpBeta[J] - gamma7)
		
				- _Complex_I*Omega33[J]*ky[i2]
		*(sigma_zero[3][2]*Beta[J] - zeta32*Alpha[J] + Del_sigma_T[3][2]*AlpBeta[J] - gamma8)
		
				- _Complex_I*Omega23[J]*ky[i2]
		*(sigma_zero[2][2]*Beta[J] - zeta22*Alpha[J] + Del_sigma_T[2][2]*AlpBeta[J] - gamma3);
		
	chi11[J] = _Complex_I*kx[i1]*u1_new[J];
	chi12[J] = _Complex_I*ky[i2]*u1_new[J];
	chi21[J] = _Complex_I*kx[i1]*u2_new[J];
	chi22[J] = _Complex_I*ky[i2]*u2_new[J];
	chi13[J] = _Complex_I*kz[i3]*u1_new[J];
	chi31[J] = _Complex_I*kx[i1]*u3_new[J];
	chi33[J] = _Complex_I*kz[i3]*u3_new[J];
	chi23[J] = _Complex_I*kz[i3]*u2_new[J];
	chi32[J] = _Complex_I*ky[i2]*u3_new[J];
	eps_star12[J] = 0.5*(chi12[J]+chi21[J]);
	eps_star13[J] = 0.5*(chi13[J]+chi31[J]);
	eps_star23[J] = 0.5*(chi23[J]+chi32[J]);
	
	}}}

	fftw_execute_dft(planB,chi11,chi11);
	fftw_execute_dft(planB,chi12,chi12);
	fftw_execute_dft(planB,chi21,chi21);
	fftw_execute_dft(planB,chi22,chi22);
	fftw_execute_dft(planB,chi13,chi13);
	fftw_execute_dft(planB,chi31,chi31);
	fftw_execute_dft(planB,chi33,chi33);
	fftw_execute_dft(planB,chi23,chi23);
	fftw_execute_dft(planB,chi32,chi32);
	fftw_execute_dft(planB,eps_star12,eps_star12);
	fftw_execute_dft(planB,eps_star13,eps_star13);
	fftw_execute_dft(planB,eps_star23,eps_star23);
	for(i1=0; i1<n_x; ++i1){
	for(i2=0; i2<n_y; ++i2){
	for(i3=0; i3<n_z; ++i3){
	J = i3 + n_z*(i2+n_y*i1);
	chi11[J] = chi11[J]*inv_nxnynz;  
	chi12[J] = chi12[J]*inv_nxnynz;
	chi21[J] = chi21[J]*inv_nxnynz;
	chi22[J] = chi22[J]*inv_nxnynz;
	chi13[J] = chi13[J]*inv_nxnynz;
	chi31[J] = chi31[J]*inv_nxnynz;
	chi33[J] = chi33[J]*inv_nxnynz;
	chi23[J] = chi23[J]*inv_nxnynz;
	chi32[J] = chi32[J]*inv_nxnynz;

	eps_star11[J] = chi11[J];
	eps_star22[J] = chi22[J];
	eps_star33[J] = chi33[J];
	eps_star12[J] = eps_star12[J]*inv_nxnynz;
	eps_star13[J] = eps_star13[J]*inv_nxnynz;
	eps_star23[J] = eps_star23[J]*inv_nxnynz;
	a = alpha[(int) (n_alpha*creal(comp[J]))];
	chi11[J] = a*chi11[J];  
	chi12[J] = a*chi12[J];
	chi21[J] = a*chi21[J];
	chi22[J] = a*chi22[J];
	chi13[J] = a*chi13[J];
	chi31[J] = a*chi31[J];
	chi33[J] = a*chi33[J];
	chi23[J] = a*chi23[J];
	chi32[J] = a*chi32[J];
	}}}
	fftw_execute_dft(planF,chi11,chi11);
	fftw_execute_dft(planF,chi12,chi12);
	fftw_execute_dft(planF,chi21,chi21);
	fftw_execute_dft(planF,chi22,chi22);
	fftw_execute_dft(planF,chi13,chi13);
	fftw_execute_dft(planF,chi31,chi31);
	fftw_execute_dft(planF,chi33,chi33);
	fftw_execute_dft(planF,chi23,chi23);
	fftw_execute_dft(planF,chi32,chi32);

	for(i3=0; i3<n_x; ++i3){
	for(i4=0; i4<n_y; ++i4){
	for(i5=0; i5<n_z; ++i5){
		J = i5 + n_z*(i4+n_y*i3);
		u1_tmp[J] = u1_old[J];
		u2_tmp[J] = u2_old[J];
		u3_tmp[J] = u3_old[J];
		u1_old[J] = u1_new[J];
		u2_old[J] = u2_new[J];
		u3_old[J] = u3_new[J];
	}}}

	calculate_hom_strain(n_x,n_y,n_z,delta_x,delta_y,delta_z,Ceff,
	DeltaC,S,sig_zero_exp,comp,ave_comp,eps_star11,eps_star12,
	eps_star22,eps_star13,eps_star33,eps_star23,
	sig_app,E,n_alpha,alpha);

/* Calculte zeta = DeltaC*Hom_strain with the new homogeneous strain
 * components */

zeta11 =  DeltaC[1][1][1][1]*E[1][1]
	+ DeltaC[1][1][2][2]*E[2][2]
	+ DeltaC[1][1][3][3]*E[3][3];

zeta22 =  DeltaC[2][2][1][1]*E[1][1]
	+ DeltaC[2][2][2][2]*E[2][2]
	+ DeltaC[2][2][3][3]*E[3][3];

zeta33 =  DeltaC[3][3][1][1]*E[1][1]
	+ DeltaC[3][3][2][2]*E[2][2]
	+ DeltaC[3][3][3][3]*E[3][3];

zeta12 =  DeltaC[1][2][1][2]*E[1][2]
	+ DeltaC[1][2][2][1]*E[2][1];

zeta21 =  DeltaC[2][1][1][2]*E[1][2]
	+ DeltaC[2][1][2][1]*E[2][1];

zeta13 =  DeltaC[1][3][1][3]*E[1][3]
	+ DeltaC[1][3][3][1]*E[3][1];

zeta31 =  DeltaC[3][1][1][3]*E[1][3]
	+ DeltaC[3][1][3][1]*E[3][1];

zeta23 =  DeltaC[2][3][2][3]*E[2][3]
	+ DeltaC[2][3][3][2]*E[3][2];

zeta32 =  DeltaC[3][2][2][3]*E[2][3]
	+ DeltaC[3][2][3][2]*E[3][2];

	fftw_execute_dft(planB,u1_old,u1_old);
	fftw_execute_dft(planB,u2_old,u2_old);
	fftw_execute_dft(planB,u3_old,u3_old);

/* Calculate the error */

	ERR = 0.0;
	for(i1=0; i1<n_x; ++i1){
	for(i2=0; i2<n_y; ++i2){
	for(i3=0; i3<n_z; ++i3){
		J = i3 + n_z*(i2+n_y*i1);	
		u1_old[J] = u1_old[J]*inv_nxnynz;  
		u2_old[J] = u2_old[J]*inv_nxnynz; 
		u3_old[J] = u3_old[J]*inv_nxnynz; 
		tmp  = creal((u1_old[J]-u1_tmp[J])
								*(u1_old[J]-u1_tmp[J])
								+(u2_old[J]-u2_tmp[J])
								*(u2_old[J]-u2_tmp[J])
								+(u3_old[J]-u3_tmp[J])
								*(u3_old[J]-u3_tmp[J]));
		ERR = ERR + tmp;
	}}}

#ifndef DEBUG
	ERR = sqrt(ERR);
	ERR = ERR*inv_nxnynz;
	if(ERR > 1.0){
	printf("The ERR is %le after %d iterations\n",ERR,ITS);
	printf("It cannot be more than unity. Hence exiting\n");
	exit(0);
	}
#endif

#ifdef DEBUG
	ERR = sqrt(ERR);
	ERR = ERR*inv_nxnynz;
	printf("%d %le\n",ITS,ERR);
	if(ERR > 1.0){
	printf("The ERR is %le after %d iterations\n",ERR,ITS);
	printf("It cannot be more than unity. Hence discontinuing.\n");
	printf("I am not exiting because I am in the DEBUG mode\n");
	printf("I will go through till the end of the program\n");
	printf("The results, however, might not make any sense\n");
	ERR = 1.0e-16;
	}
#endif

#ifdef APPROX_CALC
	ERR = sqrt(ERR);
	ERR = ERR*inv_nxnynz;
	printf("%d %le\n",ITS,ERR);
	if(ERR > 1.0){
	printf("The ERR is %le after %d iterations\n",ERR,ITS);
	printf("It cannot be more than unity. Hence discontinuing.\n");
	printf("I am not exiting because I am in the APPROX_CALC mode\n");
	printf("I will go through till the end of the program\n");
	printf("The results, however, might not make any sense\n");
	ERR = 1.0e-16;
	}
#endif
	
	ITS = ITS + 1;
}

for(i1=0; i1<n_x; ++i1){
for(i2=0; i2<n_y; ++i2){
for(i3=0; i3<n_z; ++i3){
	J = i3 + n_z*(i2+n_y*i1);
	u1_old[J] = u1_new[J];
	u2_old[J] = u2_new[J];
	u3_old[J] = u3_new[J];
}}}

fpw = fopen("../output/itnos_err","a");
fprintf(fpw,"%d %le\n",ITS-1,ERR);
fclose(fpw);

#ifndef APPROX_CALC
#ifndef DEBUG 
if(ITS == MAXITS){
	printf("No convergence after %d iterations\n",MAXITS);
	printf("Exiting.\n"); 
	fp1 = fopen("../output/unconverged/README","w");
	fprintf(fp1,"MAXITS:\n%d\n",MAXITS);
	fprintf(fp1,"n_x:\n%d\nn_y:\n%d\n",n_x,n_y);
	fprintf(fp1,"E[1][1]:\n%le\n",E[1][1]);
	fprintf(fp1,"E[1][2]:\n%le\n",E[1][2]);
	fprintf(fp1,"E[2][1]:\n%le\n",E[2][1]);
	fprintf(fp1,"E[2][2]:\n%le\n",E[2][2]);
	fclose(fp1);
	fp1 = fopen("../output/unconverged/u1","w");
	fp2 = fopen("../output/unconverged/u2","w");
	fp3 = fopen("../output/unconverged/u3","w");
	for(i1=0; i1<n_x; ++i1){
	for(i2=0; i2<n_y; ++i2){
	for(i3=0; i3<n_z; ++i3){
	J = i3 + n_z*(i2+n_y*i1);
		fprintf(fp1,"%le %le\n",creal(u1_new[J]),cimag(u1_new[J]));
		fprintf(fp2,"%le %le\n",creal(u2_new[J]),cimag(u2_new[J]));
		fprintf(fp3,"%le %le\n",creal(u3_new[J]),cimag(u3_new[J]));
	}}}
	fclose(fp1);
	fclose(fp2);
	fclose(fp3);
	/* exit(0); */
}
#endif
#endif

#ifdef DEBUG
fp1 = fopen("output/itnos_err.dat","w");
fprintf(fp1,"%d %le\n",ITS-1,ERR);
fclose(fp1);
#endif

#ifdef APPROX_CALC
fp1 = fopen("../output/itnos_err","a");
fprintf(fp1,"%d %le\n",ITS-1,ERR);
fclose(fp1);
#endif

/* Free all the allocated memory */

free(kx);
free(ky);
free(kz);
fftw_free(u1_tmp);
fftw_free(u2_tmp);
fftw_free(u3_tmp);
fftw_free(chi11);
fftw_free(chi12);
fftw_free(chi21);
fftw_free(chi22);
fftw_free(chi13);
fftw_free(chi31);
fftw_free(chi33);
fftw_free(chi23);
fftw_free(chi32);
fftw_free(Alpha);
fftw_free(Beta);
fftw_free(AlpBeta);
fftw_free(eps_star11);
fftw_free(eps_star22);
fftw_free(eps_star33);
fftw_free(eps_star12);
fftw_free(eps_star13);
fftw_free(eps_star23);

free_dmatrix(sig_zero_exp,1,3,1,3);
}
