#include<math.h>

void circular(int n_x, int n_y, int n_z, int R, double comp_precipitate,
double comp_matrix, double *comp){

int i1,i2,i3;
double a,b,c,d;
double tmp1,tmp2,tmp3;
double n_x_by2,n_y_by2,n_z_by2;

R = R*R;
a = (double) R;
b = (double) R;
c = (double) R;
n_x_by2 = (int) (n_x/2) - 1;
n_y_by2 = (int) (n_y/2) - 1;
n_z_by2 = (int) (n_z/2) - 1;
for(i1=0; i1< n_x; ++i1){
for(i2=0; i2< n_y; ++i2){
for(i3=0; i3< n_z; ++i3){
	tmp1 = 1.0*i1 - n_x/2.0;
	tmp1 = tmp1*tmp1/a;
	tmp2 = 1.0*i2 - n_y/2.0;
	tmp2 = tmp2*tmp2/b;
	tmp3 = 1.0*i3 - n_z/2.0;
	tmp3 = tmp3*tmp3/c;
	d = tmp1 + tmp2 + tmp3;

	if( d <= 1.0 )
		comp[i3 + n_z*(i2+n_y*i1)] = comp_precipitate;
	else
		comp[i3 + n_z*(i2+n_y*i1)] = comp_matrix;
}}}

}
