/** This is a machine generated C-file. **/
#include<stdio.h>
#include<stdlib.h>

extern void generate_psfile(char *NAME,char *name,
int n_x,int n_y,int N);

int main(void){
FILE *fp;
if( (fp=fopen("../../../postprocess/ps/time0_eta1.ps","r")) == NULL){
printf("Generating the psfile ../../../postprocess/ps/time0_eta1.ps\n");
generate_psfile("../../data/time0_eta1.dat","../../../postprocess/ps/time0_eta1.ps",512,512,512);
}
if( (fp=fopen("../../../postprocess/ps/time0_eta2.ps","r")) == NULL){
printf("Generating the psfile ../../../postprocess/ps/time0_eta2.ps\n");
generate_psfile("../../data/time0_eta2.dat","../../../postprocess/ps/time0_eta2.ps",512,512,512);
}
if( (fp=fopen("../../../postprocess/ps/time100_eta1.ps","r")) == NULL){
printf("Generating the psfile ../../../postprocess/ps/time100_eta1.ps\n");
generate_psfile("../../data/time100_eta1.dat","../../../postprocess/ps/time100_eta1.ps",512,512,512);
}
if( (fp=fopen("../../../postprocess/ps/time100_eta2.ps","r")) == NULL){
printf("Generating the psfile ../../../postprocess/ps/time100_eta2.ps\n");
generate_psfile("../../data/time100_eta2.dat","../../../postprocess/ps/time100_eta2.ps",512,512,512);
}
if( (fp=fopen("../../../postprocess/ps/time100_eta1.ps","r")) == NULL){
printf("Generating the psfile ../../../postprocess/ps/time100_eta1.ps\n");
generate_psfile("../../data/time100_eta1.dat","../../../postprocess/ps/time100_eta1.ps",512,512,512);
}
if( (fp=fopen("../../../postprocess/ps/time100_eta2.ps","r")) == NULL){
printf("Generating the psfile ../../../postprocess/ps/time100_eta2.ps\n");
generate_psfile("../../data/time100_eta2.dat","../../../postprocess/ps/time100_eta2.ps",512,512,512);
}
if( (fp=fopen("../../../postprocess/ps/time200_eta1.ps","r")) == NULL){
printf("Generating the psfile ../../../postprocess/ps/time200_eta1.ps\n");
generate_psfile("../../data/time200_eta1.dat","../../../postprocess/ps/time200_eta1.ps",512,512,512);
}
if( (fp=fopen("../../../postprocess/ps/time200_eta2.ps","r")) == NULL){
printf("Generating the psfile ../../../postprocess/ps/time200_eta2.ps\n");
generate_psfile("../../data/time200_eta2.dat","../../../postprocess/ps/time200_eta2.ps",512,512,512);
}
if( (fp=fopen("../../../postprocess/ps/time200_eta1.ps","r")) == NULL){
printf("Generating the psfile ../../../postprocess/ps/time200_eta1.ps\n");
generate_psfile("../../data/time200_eta1.dat","../../../postprocess/ps/time200_eta1.ps",512,512,512);
}
if( (fp=fopen("../../../postprocess/ps/time200_eta2.ps","r")) == NULL){
printf("Generating the psfile ../../../postprocess/ps/time200_eta2.ps\n");
generate_psfile("../../data/time200_eta2.dat","../../../postprocess/ps/time200_eta2.ps",512,512,512);
}
if( (fp=fopen("../../../postprocess/ps/time300_eta1.ps","r")) == NULL){
printf("Generating the psfile ../../../postprocess/ps/time300_eta1.ps\n");
generate_psfile("../../data/time300_eta1.dat","../../../postprocess/ps/time300_eta1.ps",512,512,512);
}
if( (fp=fopen("../../../postprocess/ps/time300_eta2.ps","r")) == NULL){
printf("Generating the psfile ../../../postprocess/ps/time300_eta2.ps\n");
generate_psfile("../../data/time300_eta2.dat","../../../postprocess/ps/time300_eta2.ps",512,512,512);
}
return(0);
}
