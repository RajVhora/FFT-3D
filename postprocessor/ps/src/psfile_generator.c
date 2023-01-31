#include<stdio.h>
#include<stdlib.h>

void generate_psfile(char *NAME, char *name, int n_x, int n_y, int n_z, int N){

FILE *fp;
int i,j,k;
int temp;
size_t tmp;
double *c;

c=(double *) malloc((size_t)n_x*n_y*n_z*sizeof(double));

/** Let us open the data file and read the binary data **/

if( (fp = fopen(NAME,"r")) == NULL){
printf("Unable to open the data file to read.\n");
printf("Exiting from the functions generate_psfile%s\n",name);
exit(0);
}
else{
fp = fopen(NAME,"r");
}
tmp = fread(&c[0],sizeof(double),(size_t) n_x*n_y*n_z,fp);
fclose(fp);

/* Let us open the ps file to be written (if it isn't there already) */

if( (fp = fopen(name,"w")) == NULL){
printf("Unable to open the ps file to write.Exiting\n");
exit(0);
}
else{
fp = fopen(name,"w");
}

/** The preamble to the ps file **/
/**
fprintf(fp,"%%!PS-Adobe-3.0 EPSF-2.0\n");
**/
fprintf(fp,"%%!PS-Adobe-2.0 EPSF-2.0\n");
fprintf(fp,"%%%%Generated by C-program output\n");
fprintf(fp,"%%%%Creator: Guru (inherited from Saswata)\n");
fprintf(fp,"%%%%Title: Microstructure %s\n",NAME);
fprintf(fp,"%%%%BoundingBox: 0 0 %d %d\n",n_x,n_y);
if(n_x <= 512) fprintf(fp,"%d 0 translate\n",n_x);
else fprintf(fp,"512 0 translate\n");
fprintf(fp,"90 rotate\n");
if(n_x<=512) fprintf(fp,"%d %d scale\n",n_x,n_y);
else fprintf(fp,"512 512 scale\n");
fprintf(fp,"/picstr %d string def\n",n_x);
fprintf(fp,"%d %d 8 [%d 0 0 -%d 0 %d]\n",n_y,n_x,n_x,n_y,n_x);
fprintf(fp,"{currentfile picstr readhexstring pop} image\n");

/** Let us write the ps file proper **/

for(i=0; i < n_x; ++i){
	for(j=0; j < n_y; ++j){
	for(k=0; k < n_z; ++k){
		temp = (int) (255.0*c[k + n_z*(j+n_y*i)]); //>>>> check
		if(temp < 0) temp = 0;
		else if(temp > 255) temp = 255;
		fprintf(fp,"%02x",temp);
	}
	fprintf(fp,"\n");
	}
	fprintf(fp,"\n");
}

fprintf(fp,"grestore\n");
fprintf(fp,"showpage\n");
fprintf(fp,"%%%%EOF\n");
fclose(fp);

free(c);
}

