/********************************************************************
This routine returns the zeroth-order displacement vector
corresponding to the zeroth order solution - The displacement vector
is returned in the Fourier space, and not in real space
*********************************************************************/

#include <complex.h>
#include <fftw3.h>

#include "../headers/functions.h"
#include "../headers/nrutil.h"

void calculate_uzero(int n_x, int n_y, int n_z, int half_nx, int half_ny, int half_nz,
double delta_kx, double delta_ky, double delta_kz, double ave_comp, fftw_complex *comp, 
double Ceff[4][4][4][4], double **sigma_zero, fftw_complex *u1_old, 
fftw_complex *u2_old, fftw_complex *u3_old, int n_beta, double *beta, double *Omega11,
double *Omega12, double *Omega21, double *Omega22, double *Omega13, double *Omega31, double *Omega33,
 double *Omega23, double *Omega32, fftw_plan planF){

int i1,i2,i3,i4,i5;
double n[4];
int J;
fftw_complex *Beta;

/* Calculate the Fourier transform of the eigenstrain interpolation
 * function */

Beta = fftw_malloc(n_x*n_y*n_z* sizeof(fftw_complex));
for(i3=0; i3<n_x; ++i3){
for(i4=0; i4<n_y; ++i4){
for(i5=0; i5<n_z; ++i5){
	J = i5 + n_z*(i4+n_y*i3);
	u1_old[J] = 0.0;
	u2_old[J] = 0.0;
	u3_old[J] = 0.0;
	Beta[J] = beta[(int) (n_beta*creal(comp[J]))] + _Complex_I*0.0;
}}}
fftw_execute_dft(planF,Beta,Beta);

/* Calculate the displacements */

for(i1=0; i1 < n_x; ++i1){
if(i1 < half_nx) n[1] =  i1*delta_kx;
else n[1] = (i1-n_x)*delta_kx;
for(i2=0; i2 < n_y; ++i2){
if(i2 < half_ny) n[2] =  i2*delta_ky;
else n[2] = (i2-n_y)*delta_ky;
for(i3=0; i3 < n_z; ++i3){
if(i3 < half_nz) n[3] =  i3*delta_kz;
else n[3] = (i3-n_z)*delta_kz;
	J = i3 + n_z*(i2+n_y*i1);
	for(i4=1; i4<4; ++i4){
		u1_old[J] = u1_old[J] 
		- _Complex_I*Omega11[J]*n[i4]*sigma_zero[1][i4]*Beta[J]
		- _Complex_I*Omega21[J]*n[i4]*sigma_zero[2][i4]*Beta[J]
		- _Complex_I*Omega31[J]*n[i4]*sigma_zero[3][i4]*Beta[J];
		u2_old[J] = u2_old[J] 
		- _Complex_I*Omega12[J]*n[i4]*sigma_zero[1][i4]*Beta[J]
		- _Complex_I*Omega22[J]*n[i4]*sigma_zero[2][i4]*Beta[J]
		- _Complex_I*Omega32[J]*n[i4]*sigma_zero[3][i4]*Beta[J];
		u3_old[J] = u3_old[J] 
		- _Complex_I*Omega13[J]*n[i4]*sigma_zero[1][i4]*Beta[J]
		- _Complex_I*Omega23[J]*n[i4]*sigma_zero[2][i4]*Beta[J]
		- _Complex_I*Omega33[J]*n[i4]*sigma_zero[3][i4]*Beta[J];
	}
}}}

fftw_free(Beta);
}
