
void spinodal(int n_x, int n_y, int n_z, double ave_comp, double *comp){

int i1, i2, i3;

for(i1=0; i1< n_x; ++i1){
for(i2=0; i2< n_y; ++i2){
for(i3=0; i3< n_z; ++i3){
comp[i3 + n_z*(i2+n_y*i1)] = ave_comp;
}}}
	
}
