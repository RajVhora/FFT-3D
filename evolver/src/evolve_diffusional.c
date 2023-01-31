#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <complex.h>

#include <fftw3.h>

#include "../headers/functions.h"
#include "../headers/nrutil.h"

void evolve_diffusional(int n_x, int n_y, int n_z, double delta_x, double delta_y, double delta_z, 
double kappa, double A, double delta_t, int time_steps, double t0,
fftw_complex *comp, double ave_comp, int STEPS, int plan_indicator, 
int SNA, int NOISEADDER, double noise_str){

FILE *fpd;
FILE *fptmp;

int INDEX=0;
int i1,i2,i3;
int J;

int half_nx, half_ny, half_nz;
double kx, ky, kz, delta_kx, delta_ky, delta_kz;
double k2, k4;
double inv_denom;
int nxnynz;
double inv_nxnynz;
char NAME[50];
double *c;
double realc;

size_t tmp;
unsigned FLAG;

fftw_complex *g;
fftw_plan planF, planB;

nxnynz = n_x*n_y*n_z;
inv_nxnynz = 1.0/nxnynz;

/* Allocate memory to g */

g = fftw_malloc(nxnynz* sizeof(fftw_complex));

/* Decide on the FFTW plan determination */

if(plan_indicator == 0){
FLAG = FFTW_ESTIMATE;
}
else if(plan_indicator == 1){
FLAG = FFTW_MEASURE;
}
else if(plan_indicator == 2){
FLAG = FFTW_PATIENT;
}
else if(plan_indicator == 3){
FLAG = FFTW_EXHAUSTIVE;
}
else{
printf("Unable to determine the plan type\n");
printf("Exiting from evolve_diffusional.c\n");
exit(0);
}

/* Create two plans - for forward and backward transforms*/

planF =
fftw_plan_dft_3d(n_x,n_y,n_z,g,g,FFTW_FORWARD,FLAG);
planB =
fftw_plan_dft_3d(n_x,n_y,n_z,g,g,FFTW_BACKWARD,FLAG);

/* Write the initial configuration */

c = (double *) malloc((size_t) nxnynz* sizeof(double));
for(i1=0; i1 < n_x; ++i1){
for(i2=0; i2 < n_y; ++i2){
for(i3=0; i3 < n_z; ++i3){
J = i3 + n_z*(i2+n_y*i1);
c[J] = creal(comp[J]);
}}}

sprintf(NAME,"../output/data/time%d.dat",
	(int) (INDEX + t0));
fpd = fopen(NAME,"w");
tmp = fwrite(&c[0],sizeof(double),(size_t) nxnynz,fpd);
fclose(fpd);

/* Evolve */

half_nx = (int) n_x/2;
half_ny = (int) n_y/2;
half_nz = (int) n_z/2;

delta_kx = (2.0*M_PI)/(n_x*delta_x);
delta_ky = (2.0*M_PI)/(n_y*delta_y);
delta_kz = (2.0*M_PI)/(n_z*delta_z);

/* Calculate g and its Fourier transform */

fptmp = fopen("../output/noise_info","w");

for(INDEX=0; INDEX < time_steps+1; ++INDEX){

for(i1=0; i1<n_x; ++i1){
for(i2=0; i2<n_y; ++i2){
for(i3=0; i3 < n_z; ++i3){
J = i3 + n_z*(i2+n_y*i1);
realc = creal(comp[J]);
g[J] = 2.0*A*realc*(1.0-realc)*(1.0-2.0*realc) + 0.0*_Complex_I;
}}}
fftw_execute(planF);
fftw_execute_dft(planF,comp,comp);

/* Evolve the composition profile */

for(i1=0; i1 < n_x; ++i1){
	if(i1 < half_nx) kx = i1*delta_kx;
	else kx = (i1-n_x)*delta_kx;
for(i2=0; i2 < n_y; ++i2){
	if(i2 < half_ny) ky = i2*delta_ky;
	else ky = (i2-n_y)*delta_ky;
for(i3=0; i3 < n_z; ++i3){
	if(i3 < half_nz) kz = i3*delta_kz;
	else kz = (i3-n_z)*delta_kz;

J = i3 + n_z*(i2+n_y*i1);

k2 = kx*kx + ky*ky + kz*kz;
k4 = k2*k2;

inv_denom = 1.0 + 2.0*kappa*k4*delta_t;
inv_denom = 1.0/inv_denom;

comp[J] = inv_denom*( comp[J] - k2*delta_t*g[J] );

}}}

/* Get the composition back to real space */

fftw_execute_dft(planB,comp,comp);

for(i1=0; i1<n_x; ++i1){
for(i2=0; i2<n_y; ++i2){
for(i3=0; i3<n_z; ++i3){
J = i3 + n_z*(i2+n_y*i1);
comp[J] = comp[J]*inv_nxnynz;
}}}

/* Once in a while check that the composition is within bounds,
 * and write the output **/

if(INDEX != 0 && INDEX%STEPS ==0){
	for(i1=0; i1 < n_x; ++i1){
	for(i2=0; i2 < n_y; ++i2){
	for(i3=0; i3 < n_z; ++i3){
	J = i3 + n_z*(i2+n_y*i1);
		c[J] = creal(comp[J]);
		if(c[J] <= -0.1 || c[J] >= 1.1){
		printf("The composition goes out of bounds\n");
		printf("Exiting\n");
		exit(0);
		}
	}}}
	
	sprintf(NAME,"../output/data/time%d.dat",
		(int) (INDEX + t0) );
	fpd = fopen(NAME,"w");
	tmp = fwrite(&c[0],sizeof(double),(size_t) nxnynz,fpd);
	fclose(fpd);

}


if(INDEX != 0 && INDEX%NOISEADDER == 0 && SNA == 1){
	fprintf(fptmp,"Step %d: Noise added\n",INDEX);
	for(i1=0; i1 < n_x; ++i1){
	for(i2=0; i2 < n_y; ++i2){
	for(i3=0; i3 < n_z; ++i3){
	J = i3 + n_z*(i2+n_y*i1);
		c[J] = creal(comp[J]);
	}}}
	add_noise(n_x,n_y,n_z,c,noise_str);
}

}
fclose(fptmp);

/* Free the variables and destroy the plans */

fftw_free(g);
fftw_destroy_plan(planF);
fftw_destroy_plan(planB);
free(c);
}
