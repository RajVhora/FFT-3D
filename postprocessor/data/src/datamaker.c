/** This is a machine generated C-file. **/
#include<stdio.h>
#include<stdlib.h>

extern void generate_datafile(char *NAME,char *name,
int n_x,int n_y,int N);

int main(void){
FILE *fp;
if( (fp=fopen("../../../postprocess/data/time0_eta1.data","r")) == NULL){
printf("Generating the datafile ../../../postprocess/data/time0_eta1.data\n");
generate_datafile("../../data/time0_eta1.dat","../../../postprocess/data/time0_eta1.data",512,512,512);
}
if( (fp=fopen("../../../postprocess/data/time0_eta2.data","r")) == NULL){
printf("Generating the datafile ../../../postprocess/data/time0_eta2.data\n");
generate_datafile("../../data/time0_eta2.dat","../../../postprocess/data/time0_eta2.data",512,512,512);
}
if( (fp=fopen("../../../postprocess/data/time100_eta1.data","r")) == NULL){
printf("Generating the datafile ../../../postprocess/data/time100_eta1.data\n");
generate_datafile("../../data/time100_eta1.dat","../../../postprocess/data/time100_eta1.data",512,512,512);
}
if( (fp=fopen("../../../postprocess/data/time100_eta2.data","r")) == NULL){
printf("Generating the datafile ../../../postprocess/data/time100_eta2.data\n");
generate_datafile("../../data/time100_eta2.dat","../../../postprocess/data/time100_eta2.data",512,512,512);
}
if( (fp=fopen("../../../postprocess/data/time100_eta1.data","r")) == NULL){
printf("Generating the datafile ../../../postprocess/data/time100_eta1.data\n");
generate_datafile("../../data/time100_eta1.dat","../../../postprocess/data/time100_eta1.data",512,512,512);
}
if( (fp=fopen("../../../postprocess/data/time100_eta2.data","r")) == NULL){
printf("Generating the datafile ../../../postprocess/data/time100_eta2.data\n");
generate_datafile("../../data/time100_eta2.dat","../../../postprocess/data/time100_eta2.data",512,512,512);
}
if( (fp=fopen("../../../postprocess/data/time200_eta1.data","r")) == NULL){
printf("Generating the datafile ../../../postprocess/data/time200_eta1.data\n");
generate_datafile("../../data/time200_eta1.dat","../../../postprocess/data/time200_eta1.data",512,512,512);
}
if( (fp=fopen("../../../postprocess/data/time200_eta2.data","r")) == NULL){
printf("Generating the datafile ../../../postprocess/data/time200_eta2.data\n");
generate_datafile("../../data/time200_eta2.dat","../../../postprocess/data/time200_eta2.data",512,512,512);
}
if( (fp=fopen("../../../postprocess/data/time200_eta1.data","r")) == NULL){
printf("Generating the datafile ../../../postprocess/data/time200_eta1.data\n");
generate_datafile("../../data/time200_eta1.dat","../../../postprocess/data/time200_eta1.data",512,512,512);
}
if( (fp=fopen("../../../postprocess/data/time200_eta2.data","r")) == NULL){
printf("Generating the datafile ../../../postprocess/data/time200_eta2.data\n");
generate_datafile("../../data/time200_eta2.dat","../../../postprocess/data/time200_eta2.data",512,512,512);
}
if( (fp=fopen("../../../postprocess/data/time300_eta1.data","r")) == NULL){
printf("Generating the datafile ../../../postprocess/data/time300_eta1.data\n");
generate_datafile("../../data/time300_eta1.dat","../../../postprocess/data/time300_eta1.data",512,512,512);
}
if( (fp=fopen("../../../postprocess/data/time300_eta2.data","r")) == NULL){
printf("Generating the datafile ../../../postprocess/data/time300_eta2.data\n");
generate_datafile("../../data/time300_eta2.dat","../../../postprocess/data/time300_eta2.data",512,512,512);
}
return(0);
}
