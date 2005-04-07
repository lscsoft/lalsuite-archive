#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "rastermagic.h"
#include "grid.h"

float resolution;
SKY_GRID *grid=NULL;
RGBPic *pic=NULL;
float *data=NULL;

/* bunch of variables that we don't need, but rastermagic and grid want */
FILE *LOG=NULL;
double band_axis[3]={1.0, 0.0, 0.0};
char *output_dir="";
char args_info[1024*65]; /* dummy */

void *do_alloc(long a, long b)
{
void *r;
int i=0;
r=calloc(a,b);
while(r==NULL){
        fprintf(stderr,"Could not allocate %ld chunks of %ld bytes each (%ld bytes total)\n",a,b,a*b);
        if(i>10)exit(-1);
        sleep(10);
        r=calloc(a,b);
        i++;
        }
return r;
}

int clear_name_png(char *name)
{
return 1;
}

int main(int argc, char *argv[])
{
FILE *fin;

LOG=stderr;

if(argc<4){
	fprintf(stderr, "Usage: %s resolution file.dat skymap.png\n", argv[0]);
	return -1;
	}
resolution=atof(argv[1]);

fin=fopen(argv[2], "r");
if(fin==NULL){
	perror("Cannot read file:");
	return -1;
	}

grid=make_sin_theta_grid(resolution);

data=do_alloc(grid->npoints, sizeof(*data));
fread(data, sizeof(*data), grid->npoints, fin);

pic=make_RGBPic(grid->max_n_ra+140, grid->max_n_dec);
plot_grid_f(pic, grid, data, 1);

RGBPic_dump_png(argv[3], pic);
return 0;
}
