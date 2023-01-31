/*************************************************************
We previously assumed that the eigenstrain is dilatational and the elastic
constants are cubic, and hence shortened the lengthy summations, 
including only the non-zero terms. However, now the sig_zero
expectation values are calculated without making any such assumption
and hence skipping terms.
M P Gururajan and N Davessar, 1 April, 20
***************************************************************/

#include<complex.h>
#include<fftw3.h>

#include "../headers/functions.h"

void calculate_sig_zero_exp(int n_x, int n_y, int n_z, double delta_x, 
double delta_y, double delta_z, fftw_complex *comp, double ave_comp, 
double Ceff[4][4][4][4], double DeltaC[4][4][4][4], double **epsilon_T,
double **sig_zero_exp, int n_alpha, int n_beta, double *alpha, 
double *beta){

int i5,i6,i7,J;
double c;
double a, b, d, e, f, g, h;

/* Calculate the 11, 22 and 12 component of 'expectation' value of eigenstrain */

sig_zero_exp[1][1] = 0.0;
sig_zero_exp[2][2] = 0.0;
sig_zero_exp[3][3] = 0.0;
sig_zero_exp[1][2] = 0.0;
sig_zero_exp[1][3] = 0.0;
sig_zero_exp[2][3] = 0.0;
for(i5=0; i5 < n_x; ++i5){
for(i6=0; i6 < n_y; ++i6){
for(i7=0; i7 < n_z; ++i7){
J = i7 + n_z*(i6+n_y*i5);
c = creal(comp[J]);
a = alpha[(int) (n_alpha*c)];
b = beta[(int) (n_beta*c)] * epsilon_T[1][1];
d = beta[(int) (n_beta*c)] * epsilon_T[2][2];
e = beta[(int) (n_beta*c)] * epsilon_T[3][3];
f = beta[(int) (n_beta*c)] * epsilon_T[1][2];
g = beta[(int) (n_beta*c)] * epsilon_T[1][3];
h = beta[(int) (n_beta*c)] * epsilon_T[2][3];
sig_zero_exp[1][1] = sig_zero_exp[1][1] 
+ (Ceff[1][1][1][1]+DeltaC[1][1][1][1]*a)*b 
+ (Ceff[1][1][2][2]+DeltaC[1][1][2][2]*a)*d 
+ (Ceff[1][1][3][3]+DeltaC[1][1][3][3]*a)*e 
+ (Ceff[1][1][1][2]+DeltaC[1][1][1][2]*a)*f
+ (Ceff[1][1][2][1]+DeltaC[1][1][2][1]*a)*f
+ (Ceff[1][1][1][3]+DeltaC[1][1][1][3]*a)*g
+ (Ceff[1][1][3][1]+DeltaC[1][1][3][1]*a)*g
+ (Ceff[1][1][2][3]+DeltaC[1][1][2][3]*a)*h
+ (Ceff[1][1][3][2]+DeltaC[1][1][3][2]*a)*h;

sig_zero_exp[2][2] = sig_zero_exp[2][2] 
+ (Ceff[2][2][1][1]+DeltaC[2][2][1][1]*a)*b 
+ (Ceff[2][2][2][2]+DeltaC[2][2][2][2]*a)*d 
+ (Ceff[2][2][3][3]+DeltaC[2][2][3][3]*a)*e 
+ (Ceff[2][2][1][2]+DeltaC[2][2][1][2]*a)*f
+ (Ceff[2][2][2][1]+DeltaC[2][2][2][1]*a)*f
+ (Ceff[2][2][1][3]+DeltaC[2][2][1][3]*a)*g
+ (Ceff[2][2][3][1]+DeltaC[2][2][3][1]*a)*g
+ (Ceff[2][2][2][3]+DeltaC[2][2][2][3]*a)*h
+ (Ceff[2][2][3][2]+DeltaC[2][2][3][2]*a)*h;

sig_zero_exp[3][3] = sig_zero_exp[3][3] 
+ (Ceff[3][3][1][1]+DeltaC[3][3][1][1]*a)*b 
+ (Ceff[3][3][2][2]+DeltaC[3][3][2][2]*a)*d 
+ (Ceff[3][3][3][3]+DeltaC[3][3][3][3]*a)*e 
+ (Ceff[3][3][1][2]+DeltaC[3][3][1][2]*a)*f
+ (Ceff[3][3][2][1]+DeltaC[3][3][2][1]*a)*f
+ (Ceff[3][3][1][3]+DeltaC[3][3][1][3]*a)*g
+ (Ceff[3][3][3][1]+DeltaC[3][3][3][1]*a)*g
+ (Ceff[3][3][2][3]+DeltaC[3][3][2][3]*a)*h
+ (Ceff[3][3][3][2]+DeltaC[3][3][3][2]*a)*h;

sig_zero_exp[1][2] = sig_zero_exp[1][2] 
+ (Ceff[1][2][1][1]+DeltaC[1][2][1][1]*a)*b 
+ (Ceff[1][2][2][2]+DeltaC[1][2][2][2]*a)*d 
+ (Ceff[1][2][3][3]+DeltaC[1][2][3][3]*a)*e 
+ (Ceff[1][2][1][2]+DeltaC[1][2][1][2]*a)*f
+ (Ceff[1][2][2][1]+DeltaC[1][2][2][1]*a)*f
+ (Ceff[1][2][1][3]+DeltaC[1][2][1][3]*a)*g
+ (Ceff[1][2][3][1]+DeltaC[1][2][3][1]*a)*g
+ (Ceff[1][2][2][3]+DeltaC[1][2][2][3]*a)*h
+ (Ceff[1][2][3][2]+DeltaC[1][2][3][2]*a)*h;

sig_zero_exp[1][3] = sig_zero_exp[1][3] 
+ (Ceff[1][3][1][1]+DeltaC[1][3][1][1]*a)*b 
+ (Ceff[1][3][2][2]+DeltaC[1][3][2][2]*a)*d 
+ (Ceff[1][3][3][3]+DeltaC[1][3][3][3]*a)*e 
+ (Ceff[1][3][1][2]+DeltaC[1][3][1][2]*a)*f
+ (Ceff[1][3][2][1]+DeltaC[1][3][2][1]*a)*f
+ (Ceff[1][3][1][3]+DeltaC[1][3][1][3]*a)*g
+ (Ceff[1][3][3][1]+DeltaC[1][3][3][1]*a)*g
+ (Ceff[1][3][2][3]+DeltaC[1][3][2][3]*a)*h
+ (Ceff[1][3][3][2]+DeltaC[1][3][3][2]*a)*h;

sig_zero_exp[2][3] = sig_zero_exp[2][3] 
+ (Ceff[2][3][1][1]+DeltaC[2][3][1][1]*a)*b 
+ (Ceff[2][3][2][2]+DeltaC[2][3][2][2]*a)*d 
+ (Ceff[2][3][3][3]+DeltaC[2][3][3][3]*a)*e 
+ (Ceff[2][3][1][2]+DeltaC[2][3][1][2]*a)*f
+ (Ceff[2][3][2][1]+DeltaC[2][3][2][1]*a)*f
+ (Ceff[2][3][1][3]+DeltaC[2][3][1][3]*a)*g
+ (Ceff[2][3][3][1]+DeltaC[2][3][3][1]*a)*g
+ (Ceff[2][3][2][3]+DeltaC[2][3][2][3]*a)*h
+ (Ceff[2][3][3][2]+DeltaC[2][3][3][2]*a)*h;
}}}

sig_zero_exp[1][1] = sig_zero_exp[1][1]/(n_x*n_y*n_z);
sig_zero_exp[2][2] = sig_zero_exp[2][2]/(n_x*n_y*n_z);
sig_zero_exp[3][3] = sig_zero_exp[3][3]/(n_x*n_y*n_z);
sig_zero_exp[1][2] = sig_zero_exp[1][2]/(n_x*n_y*n_z);
sig_zero_exp[1][3] = sig_zero_exp[1][3]/(n_x*n_y*n_z);
sig_zero_exp[2][3] = sig_zero_exp[2][3]/(n_x*n_y*n_z);
sig_zero_exp[2][1] = sig_zero_exp[1][2];
sig_zero_exp[3][1] = sig_zero_exp[1][3];
sig_zero_exp[3][2] = sig_zero_exp[2][3];

}
