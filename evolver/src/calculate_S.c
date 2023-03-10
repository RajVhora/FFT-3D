/****************************************************************
The compliance tensor for the plane strain problem is given in the
following paper (Equation 5.11 - p. 3093):

I. Schmidt and D. Gross

Directional coarsening in Ni-base superalloys: analytical results for
an elasticity-based model 

Proc. R. Soc. Lond. A (1999), 455, p. 3085-3106.

However, note that the factors multiplying s12 and s44 are not two and
four (as in 3D compliance tensor) for the S_{1122} and S_{1212},
but unity and two, respectively.
******************************************************************/

#include<complex.h>
#include<fftw3.h>

#include "../headers/functions.h"

void calculate_S(double Ceff[4][4][4][4], double S[4][4][4][4]){

int i1,i2,i3,i4;
double A, nu;

/* nu and A are the Poisson's ratio and anisotropy parameter
 * repsectively */

nu = 0.5*Ceff[1][1][2][2]/(Ceff[1][1][2][2]+Ceff[1][2][1][2]);
A = 2.0*Ceff[1][2][1][2]/(Ceff[1][1][1][1]-Ceff[1][1][2][2]);

/* Initialise the compliance tensor */

for(i1=1;i1<4;++i1){
for(i2=1;i2<4;++i2){
for(i3=1;i3<4;++i3){
for(i4=1;i4<4;++i4){
S[i1][i2][i3][i4] = 0.0;
}}}}

/* Calculate the non-zero compliance tensor components */

S[1][1][1][1] = (Ceff[1][1][1][1] + Ceff[1][1][2][2])
		/((Ceff[1][1][1][1] - Ceff[1][1][2][2])*(Ceff[1][1][1][1] + 2.*Ceff[1][1][2][2])); 
S[1][1][2][2] = -(Ceff[1][1][2][2])
		/((Ceff[1][1][1][1] - Ceff[1][1][2][2])*(Ceff[1][1][1][1] + 2.*Ceff[1][1][2][2]));
S[1][2][1][2] = (1.0)/(Ceff[1][2][1][2]);
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
