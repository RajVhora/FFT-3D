/****************************************************************
The compliance tensor for the plane strain problem is given in the
following paper (Equation 5.11 - p. 3093):

I. Schmidt and D. Gross

Directional coarsening in Ni-base superalloys: nalytical results for
an elasticity-based model

Proc. R. Soc. Lond. A (1999), 455, p. 3085-3106.

However, note that the factors multiplying s12 and s44 are not two and
four (as in 3D compliance tensor) for the S_{1122} and S_{1212},
but unity and two, respectively.
******************************************************************/

#include<complex.h>
#include<fftw3.h>

#include "../headers/functions.h"

void calculate_S_exp(int n_x, int n_y, int n_z, double delta_x, double delta_y, double delta_z,
fftw_complex *comp, double ave_comp, double Ceff[4][4][4][4], 
double DeltaC[4][4][4][4], double S[4][4][4][4], int n_alpha, 
double *alpha){

int i1,i2,i3,i4;
double C[3][3][3][3];
double A, nu;
double sum1;

/* Calculate 'effective' modulus */

sum1 = 0.0;
for(i1=0; i1<n_x; ++i1){
for(i2=0; i2<n_y; ++i2){
for(i3=0; i3<n_z; ++i3){
sum1 = sum1 + alpha[(int) (n_alpha*creal(comp[i3 + n_z*(i2+n_y*i1)]))];
}}}
sum1 = sum1/(n_x*n_y);

C[1][1][1][1] = Ceff[1][1][1][1] + DeltaC[1][1][1][1]*sum1;
C[1][2][1][2] = Ceff[1][2][1][2] + DeltaC[1][2][1][2]*sum1;
C[1][1][2][2] = Ceff[1][1][2][2] + DeltaC[1][1][2][2]*sum1;

/* Get 'effective' Poisson's ratio and Zener anisotropy parameter */

nu = 0.5*C[1][1][2][2]/(C[1][1][2][2]+C[1][2][1][2]);
A = 2.0*C[1][2][1][2]/(C[1][1][1][1]-C[1][1][2][2]);

/* Initialise the compliance tensor */

for(i1=1;i1<4;++i1){
for(i2=1;i2<4;++i2){
for(i3=1;i3<4;++i3){
for(i4=1;i4<4;++i4){
S[i1][i2][i3][i4] = 0.0;
}}}}

/* Get the non-zero compliance components */

S[1][1][1][1] = (C[1][1][1][1] + C[1][1][2][2])
		/((C[1][1][1][1] - C[1][1][2][2])*(C[1][1][1][1] + 2.*C[1][1][2][2]));  
S[1][1][2][2] = -(C[1][1][2][2])
		/((C[1][1][1][1] - C[1][1][2][2])*(C[1][1][1][1] + 2.*C[1][1][2][2]));
S[1][2][1][2] = (1.0)/(C[1][2][1][2]);

/* c11 */
S[2][2][2][2] = S[1][1][1][1];
S[3][3][3][3] = S[1][1][1][1];
/* c12 */
S[2][2][1][1] = S[1][1][2][2];
S[1][1][3][3] = S[3][3][1][1] = S[1][1][2][2];
S[2][2][3][3] = S[3][3][2][2] = S[1][1][2][2];
/* c44 */
S[1][2][2][1] = S[2][1][1][2] = S[2][1][2][1]  =
S[1][2][1][2];
S[1][3][1][3] = S[1][3][3][1] = S[3][1][1][3] = S[3][1][3][1]  =
S[1][2][1][2];
S[2][3][2][3] = S[2][3][3][2] = S[3][2][2][3] = S[3][2][3][2]  =
S[1][2][1][2];

}
