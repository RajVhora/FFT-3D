#include <math.h>
#include <float.h>

#include "../headers/nrutil.h"

void calculate_Omega(int n_x, int n_y, int n_z, int half_nx, int half_ny, int half_nz,
double *kx, double *ky, double *kz, double Ceff[4][4][4][4], double *Omega11, 
double *Omega12, double*Omega21, double *Omega22,
double *Omega13, double*Omega31, double *Omega33,
double *Omega23, double*Omega32){

int i5,i6,i7;

double Omega_inv11, Omega_inv12, Omega_inv21, Omega_inv22, Omega_inv13, Omega_inv31, Omega_inv33, Omega_inv23, Omega_inv32;
double det_Omega_inv;
double n[4];
int J;

/* First step is to calculate the inverse of the acoustic tensor. Then
 * we invert it - We assume plane elasticity,  and hence matrix inversion
 * is straightforward */

for(i5=0; i5<n_x; ++i5){
	 n[1] = kx[i5];
for(i6=0; i6<n_y; ++i6){
	 n[2] = ky[i6];
for(i7=0; i7<n_z; ++i7){
	 n[3] = kz[i7];
		J = i7 + n_z*(i6+n_y*i5);
		Omega_inv11 = Ceff[1][1][1][1]*n[1]*n[1]+
			Ceff[1][1][2][1]*n[1]*n[2]+
			Ceff[1][2][1][1]*n[2]*n[1]+
			Ceff[1][2][2][1]*n[2]*n[2]+

			Ceff[1][1][3][1]*n[1]*n[3]+
			Ceff[1][3][1][1]*n[3]*n[1]+
			Ceff[1][3][3][1]*n[3]*n[3]+

			Ceff[1][2][3][1]*n[2]*n[3]+
			Ceff[1][3][2][1]*n[3]*n[2];

		Omega_inv12 = Ceff[1][1][1][2]*n[1]*n[1]+
			Ceff[1][1][2][2]*n[1]*n[2]+
			Ceff[1][2][1][2]*n[2]*n[1]+
			Ceff[1][2][2][2]*n[2]*n[2]+

			Ceff[1][1][3][2]*n[1]*n[3]+
			Ceff[1][3][1][2]*n[3]*n[1]+
			Ceff[1][3][3][2]*n[3]*n[3]+

			Ceff[1][2][3][2]*n[2]*n[3]+
			Ceff[1][3][2][2]*n[3]*n[2];

		Omega_inv21 = Ceff[2][1][1][1]*n[1]*n[1]+
			Ceff[2][1][2][1]*n[1]*n[2]+
			Ceff[2][2][1][1]*n[2]*n[1]+
			Ceff[2][2][2][1]*n[2]*n[2]+

			Ceff[2][1][3][1]*n[1]*n[3]+
			Ceff[2][3][1][1]*n[3]*n[1]+
			Ceff[2][3][3][1]*n[3]*n[3]+

			Ceff[2][2][3][1]*n[2]*n[3]+
			Ceff[2][3][2][1]*n[3]*n[2];

		Omega_inv22 = Ceff[2][1][1][2]*n[1]*n[1]+
			Ceff[2][1][2][2]*n[1]*n[2]+
			Ceff[2][2][1][2]*n[2]*n[1]+
			Ceff[2][2][2][2]*n[2]*n[2]+

			Ceff[2][1][3][2]*n[1]*n[3]+
			Ceff[2][3][1][2]*n[3]*n[1]+
			Ceff[2][3][3][2]*n[3]*n[3]+

			Ceff[2][2][3][2]*n[2]*n[3]+
			Ceff[2][3][2][2]*n[3]*n[2];

		Omega_inv13 = Ceff[1][1][1][3]*n[1]*n[1]+
			Ceff[1][1][2][3]*n[1]*n[2]+
			Ceff[1][2][1][3]*n[2]*n[1]+
			Ceff[1][2][2][3]*n[2]*n[2]+

			Ceff[1][1][3][3]*n[1]*n[3]+
			Ceff[1][3][1][3]*n[3]*n[1]+
			Ceff[1][3][3][3]*n[3]*n[3]+

			Ceff[1][2][3][3]*n[2]*n[3]+
			Ceff[1][3][2][3]*n[3]*n[2];

		Omega_inv31 = Ceff[3][1][1][1]*n[1]*n[1]+
			Ceff[3][1][2][1]*n[1]*n[2]+
			Ceff[3][2][1][1]*n[2]*n[1]+
			Ceff[3][2][2][1]*n[2]*n[2]+

			Ceff[3][1][3][1]*n[1]*n[3]+
			Ceff[3][3][1][1]*n[3]*n[1]+
			Ceff[3][3][3][1]*n[3]*n[3]+

			Ceff[3][2][3][1]*n[2]*n[3]+
			Ceff[3][3][2][1]*n[3]*n[2];

		Omega_inv33 = Ceff[3][1][1][3]*n[1]*n[1]+
			Ceff[3][1][2][3]*n[1]*n[2]+
			Ceff[3][2][1][3]*n[2]*n[1]+
			Ceff[3][2][2][3]*n[2]*n[2]+

			Ceff[3][1][3][3]*n[1]*n[3]+
			Ceff[3][3][1][3]*n[3]*n[1]+
			Ceff[3][3][3][3]*n[3]*n[3]+

			Ceff[3][2][3][3]*n[2]*n[3]+
			Ceff[3][3][2][3]*n[3]*n[2];

		Omega_inv23 = Ceff[2][1][1][3]*n[1]*n[1]+
			Ceff[2][1][2][3]*n[1]*n[2]+
			Ceff[2][2][1][3]*n[2]*n[1]+
			Ceff[2][2][2][3]*n[2]*n[2]+

			Ceff[2][1][3][3]*n[1]*n[3]+
			Ceff[2][3][1][3]*n[3]*n[1]+
			Ceff[2][3][3][3]*n[3]*n[3]+

			Ceff[2][2][3][3]*n[2]*n[3]+
			Ceff[2][3][2][3]*n[3]*n[2];

		Omega_inv32 = Ceff[3][1][1][2]*n[1]*n[1]+
			Ceff[3][1][2][2]*n[1]*n[2]+
			Ceff[3][2][1][2]*n[2]*n[1]+
			Ceff[3][2][2][2]*n[2]*n[2]+

			Ceff[3][1][3][2]*n[1]*n[3]+
			Ceff[3][3][1][2]*n[3]*n[1]+
			Ceff[3][3][3][2]*n[3]*n[3]+

			Ceff[3][2][3][2]*n[2]*n[3]+
			Ceff[3][3][2][2]*n[3]*n[2];
			
	det_Omega_inv = Omega_inv11*(Omega_inv22*Omega_inv33 - Omega_inv23*Omega_inv32)-
			Omega_inv12*(Omega_inv21*Omega_inv33 - Omega_inv23*Omega_inv31)+
			Omega_inv13*(Omega_inv21*Omega_inv32 - Omega_inv22*Omega_inv31);
	if(det_Omega_inv != 0.0){
		Omega11[J] = (Omega_inv22*Omega_inv33 - Omega_inv23*Omega_inv32)/det_Omega_inv;  
		Omega22[J] = (Omega_inv11*Omega_inv33 - Omega_inv13*Omega_inv31)/det_Omega_inv;  
		Omega33[J] = (Omega_inv11*Omega_inv22 - Omega_inv12*Omega_inv21)/det_Omega_inv;  

		Omega12[J] = (Omega_inv13*Omega_inv32 - Omega_inv12*Omega_inv33)/det_Omega_inv;  
		Omega21[J] = (Omega_inv23*Omega_inv31 - Omega_inv21*Omega_inv33)/det_Omega_inv;  

		Omega13[J] = (Omega_inv12*Omega_inv23 - Omega_inv13*Omega_inv22)/det_Omega_inv;  
		Omega31[J] = (Omega_inv21*Omega_inv32 - Omega_inv22*Omega_inv31)/det_Omega_inv;

		Omega23[J] = (Omega_inv13*Omega_inv21 - Omega_inv11*Omega_inv23)/det_Omega_inv;  
		Omega32[J] = (Omega_inv12*Omega_inv31 - Omega_inv11*Omega_inv32)/det_Omega_inv;  
	}
	else{
		Omega11[J] = (Omega_inv22*Omega_inv33 - Omega_inv23*Omega_inv32);  
		Omega22[J] = (Omega_inv11*Omega_inv33 - Omega_inv13*Omega_inv31);  
		Omega33[J] = (Omega_inv11*Omega_inv22 - Omega_inv12*Omega_inv21);  

		Omega12[J] = (Omega_inv13*Omega_inv32 - Omega_inv12*Omega_inv33);  
		Omega21[J] = (Omega_inv23*Omega_inv31 - Omega_inv21*Omega_inv33);  

		Omega13[J] = (Omega_inv12*Omega_inv23 - Omega_inv13*Omega_inv22);  
		Omega31[J] = (Omega_inv21*Omega_inv32 - Omega_inv22*Omega_inv31);

		Omega23[J] = (Omega_inv13*Omega_inv21 - Omega_inv11*Omega_inv23);  
		Omega32[J] = (Omega_inv12*Omega_inv31 - Omega_inv11*Omega_inv32); 
	}

}}}

}
