#include<math.h>

void thinfilm(int n_x, int n_y, int n_z, int h, double comp_precipitate,
double comp_matrix, double *comp){

int i1,i2,i3;
int halfh;
double n_x_by2,n_y_by2,n_z_by2;
double tmp1,tmp2,tmp3,vol,d;


halfh = (int) (h/2);
vol = n_x*((n_y + halfh) - (n_y - halfh))*n_z;
for(i1=0; i1< n_x; ++i1){
for(i2=0; i2< n_y; ++i2){
for(i3=0; i3< n_z; ++i3){


	if( i1 >=0 && i1 <= n_x && i2 >= (n_y/2 -halfh) && i2 <= (n_y/2 +halfh) && i3 >=0  && i3 <=n_z )
		comp[i3 + n_z*(i2+n_y*i1)] = comp_precipitate;
	
	else
		comp[i3 + n_z*(i2+n_y*i1)] = comp_matrix;		
}}}




}
