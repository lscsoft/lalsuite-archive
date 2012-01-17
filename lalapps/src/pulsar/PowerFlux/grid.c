/*
*  Copyright (C) 2007 Vladimir Dergachev
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "global.h"
#include "grid.h"
#include "util.h"
#include "dataset.h"
#include "jobs.h"

extern FILE *LOG;

extern double average_det_velocity[3];
extern double band_axis_norm;
extern double band_axis[3];
extern double spindown;

SKY_GRID_TYPE spherical_distance(SKY_GRID_TYPE ra0, SKY_GRID_TYPE dec0,
			  SKY_GRID_TYPE ra1, SKY_GRID_TYPE dec1)
{
SKY_GRID_TYPE ds;
ds=acos(sin(dec0)*sin(dec1)+cos(dec0)*cos(dec1)*cos(ra0-ra1));
return ds;
}

SKY_GRID_TYPE fast_spherical_distance(SKY_GRID_TYPE ra0, SKY_GRID_TYPE dec0,
			  SKY_GRID_TYPE ra1, SKY_GRID_TYPE dec1)
{
SKY_GRID_TYPE ds;
ds=1.0-(sin(dec0)*sin(dec1)+cos(dec0)*cos(dec1)*cos(ra0-ra1));
return ds;
}

/* Precompute values that are used later */
void precompute_values(SKY_GRID *grid)
{
long i, k;
SKY_GRID_TYPE e2,e3,e4,e5;

for(i=0;i<GRID_E_COUNT;i++) {
	if(grid->e[i]==NULL)grid->e[i]=do_alloc(grid->npoints, sizeof(SKY_GRID_TYPE));
	}

for(k=0;k<grid->npoints;k++){
		e2=cos(M_PI_2-grid->latitude[k]);
		e3=sin(M_PI_2-grid->latitude[k]);
		e4=cos(grid->longitude[k]);
		e5=sin(grid->longitude[k]);
		/* unit vector */
		grid->e[0][k]=e3*e4;
		grid->e[1][k]=e3*e5;
		grid->e[2][k]=e2;
		/* other useful values */
		grid->e[3][k]=e3;
		grid->e[4][k]=e4;
		grid->e[5][k]=e5;
		/* these values are needed for regression of plus and cross */
		grid->e[6][k]=e4*e5;
		grid->e[7][k]=e3*e4;
		grid->e[8][k]=e3*e5;
		
		grid->e[9][k]=e2*e2*e4;
		grid->e[10][k]=e2*e2*e5;
		grid->e[11][k]=e3*e3*e4;
		grid->e[12][k]=e3*e3*e5;
		grid->e[13][k]=e2*e3*e4;
		grid->e[14][k]=e2*e3*e5;

		grid->e[15][k]=e2*e4*e4;
		grid->e[16][k]=e2*e5*e5;
		grid->e[17][k]=e3*e4*e4;
		grid->e[18][k]=e3*e5*e5;
		grid->e[19][k]=e2*e4*e5;
		
		grid->e[20][k]=e2*e2*e4*e4;
		grid->e[21][k]=e2*e2*e5*e5;
		grid->e[22][k]=e3*e3*e4*e4;
		grid->e[23][k]=e3*e3*e5*e5;

		grid->e[24][k]=e2*e2*e4*e5;
		grid->e[25][k]=e3*e3*e4*e5;
	}

}

void free_values(SKY_GRID *grid)
{
int i;
for(i=0;i<GRID_E_COUNT;i++) {
	if(grid->e[i]!=NULL) {
		free(grid->e[i]);
		grid->e[i]=NULL;
		}
	}
}

SKY_GRID *make_arcsin_grid(long num_ra, long num_dec)
{
SKY_GRID *grid;
RECT_SKY_GRID_PRIV *priv;
long i,j,k;
SKY_GRID_TYPE a,b;

/* set all up */
grid=do_alloc(1,sizeof(*grid));
grid->npoints=num_ra*num_dec;
grid->max_n_dec=num_dec;
grid->max_n_ra=num_ra;
grid->name="arcsin rectangular";
grid->latitude=do_alloc(grid->npoints, sizeof(SKY_GRID_TYPE));
grid->longitude=do_alloc(grid->npoints, sizeof(SKY_GRID_TYPE));

grid->max_x=num_dec;
grid->max_y=num_ra;
grid->x=do_alloc(grid->npoints, sizeof(*(grid->x)));
grid->y=do_alloc(grid->npoints, sizeof(*(grid->y)));

grid->nbands_size=2000;
grid->nbands=0;
grid->band_name=do_alloc(grid->nbands_size, sizeof(*grid->band_name));

grid->band=do_alloc(grid->npoints, sizeof(*grid->band));
grid->band_f=do_alloc(grid->npoints, sizeof(*grid->band_f));

for(i=0;i<GRID_E_COUNT;i++)
	grid->e[i]=do_alloc(grid->npoints, sizeof(SKY_GRID_TYPE));
priv=do_alloc(1, sizeof(*priv));
priv->num_ra=num_ra;
priv->num_dec=num_dec;
grid->grid_priv=priv;

/* fill in the coordinates */
for(i=0;i<num_ra;i++){
	a=M_PI*(1.0+2.0*i)/num_ra;
	for(j=0;j<num_dec;j++){
		k=i*num_dec+j;
		b=asin(-1.0+(1.0+2.0*j)/num_dec);
		grid->latitude[k]=b;
		grid->longitude[k]=a;
		grid->band[k]=-1;

		grid->x[k]=i;
		grid->y[k]=j;
		}
	}
return grid;
}

SKY_GRID *make_rect_grid(long num_ra, long num_dec)
{
SKY_GRID *grid;
RECT_SKY_GRID_PRIV *priv;
long i,j,k;
SKY_GRID_TYPE a,b;

/* set all up */
grid=do_alloc(1,sizeof(*grid));
grid->npoints=num_ra*num_dec;
grid->max_n_dec=num_dec;
grid->max_n_ra=num_ra;
grid->name="plain rectangular";
grid->latitude=do_alloc(grid->npoints, sizeof(SKY_GRID_TYPE));
grid->longitude=do_alloc(grid->npoints, sizeof(SKY_GRID_TYPE));

grid->max_x=num_dec;
grid->max_y=num_ra;
grid->x=do_alloc(grid->npoints, sizeof(*(grid->x)));
grid->y=do_alloc(grid->npoints, sizeof(*(grid->y)));

grid->nbands_size=2000;
grid->nbands=0;
grid->band_name=do_alloc(grid->nbands_size, sizeof(*grid->band_name));
grid->band=do_alloc(grid->npoints, sizeof(*grid->band));
grid->band_f=do_alloc(grid->npoints, sizeof(*grid->band_f));
for(i=0;i<GRID_E_COUNT;i++)
	grid->e[i]=do_alloc(grid->npoints, sizeof(SKY_GRID_TYPE));
priv=do_alloc(1, sizeof(*priv));
priv->num_ra=num_ra;
priv->num_dec=num_dec;
grid->grid_priv=priv;

/* fill in the coordinates */
for(i=0;i<num_ra;i++){
	a=M_PI*(1.0+2.0*i)/num_ra;
	for(j=0;j<num_dec;j++){
		k=i*num_dec+j;
		b=M_PI_2*(-1.0+(1.0+2.0*j)/num_dec);
		grid->latitude[k]=b;
		grid->longitude[k]=a;
		grid->band[k]=-1;

		grid->x[k]=i;
		grid->y[k]=j;
		}
	}
return grid;
}

/* This function is meant for targeted search in a small disc centered on (RA, DEC) */
/* focus-* options are required ! */
SKY_GRID *make_targeted_rect_grid(SKY_GRID_TYPE ra, SKY_GRID_TYPE dec, SKY_GRID_TYPE radius, long num_dec)
{
SKY_GRID *grid;
TARGETED_RECT_SKY_GRID_PRIV *priv;
long i,j,k;
double a,b, e1[3], e2[3];

num_dec|=1; /* Make sure num_dec is odd */


grid=do_alloc(1,sizeof(*grid));
grid->npoints=num_dec*num_dec;
grid->max_n_dec=num_dec;
grid->max_n_ra=num_dec;
grid->name="targeted rectangular";
grid->latitude=do_alloc(grid->npoints, sizeof(SKY_GRID_TYPE));
grid->longitude=do_alloc(grid->npoints, sizeof(SKY_GRID_TYPE));

grid->max_x=num_dec;
grid->max_y=num_dec;
grid->x=do_alloc(grid->npoints, sizeof(*(grid->x)));
grid->y=do_alloc(grid->npoints, sizeof(*(grid->y)));

grid->nbands_size=2000;
grid->nbands=0;
grid->band_name=do_alloc(grid->nbands_size, sizeof(*grid->band_name));
grid->band=do_alloc(grid->npoints, sizeof(*grid->band));
grid->band_f=do_alloc(grid->npoints, sizeof(*grid->band_f));
for(i=0;i<GRID_E_COUNT;i++)
	grid->e[i]=do_alloc(grid->npoints, sizeof(SKY_GRID_TYPE));
priv=do_alloc(1, sizeof(*priv));
priv->ra=ra;
priv->dec=dec;
priv->radius=radius;
priv->num_dec=num_dec;
grid->grid_priv=priv;

/* fill in the coordinates */
for(i=0;i<num_dec;i++){
	a=radius*(i*2.0-num_dec+1)/(num_dec-1);
	for(j=0;j<num_dec;j++){
		k=i*num_dec+j;
		b=radius*(j*2.0-num_dec+1)/(num_dec-1);

		/* transition to standard coordinates around point (RA, DEC) */

		/* (0, 0) -> (1, 0, 0) */
		e1[0]=cos(b)*cos(a);
		e1[1]=cos(b)*sin(a);
		e1[2]=sin(b);

		/* rotate by DEC around Oy */

		e2[0]=e1[0]*cos(dec)-e1[2]*sin(dec);
		e2[1]=e1[1];
		e2[2]=e1[0]*sin(dec)+e1[2]*cos(dec);

		/* rotate by RA around 0z */

		e1[0]=e2[0]*cos(ra)-e2[1]*sin(ra);
		e1[1]=e2[0]*sin(ra)+e2[1]*cos(ra);
		e1[2]=e2[2];


		grid->latitude[k]=asin(e1[2]);
		grid->longitude[k]=atan2(e1[1], e1[0]);
		/* Fixup (0, 0, 1) vector which would produce NaNs */
		if(e1[0]*e1[0]+e1[1]*e1[1]<=0)grid->longitude[k]=0.0;

		/* make sure right ascension is positive as in other grids */
		if(grid->longitude[k]<0.0)grid->longitude[k]+=2*M_PI;

		//fprintf(stderr, "a=%.5f b=%.5f e=(%.5f, %.5f, %.5f) (%.5f, %.5f)\n", a, b, e1[0], e1[1], e1[2], grid->longitude[k], grid->latitude[k]);

		grid->band[k]=-1;

		grid->x[k]=i;
		grid->y[k]=j;
		}
	}
return grid;
}


SKY_GRID *make_sin_theta_grid(SKY_GRID_TYPE resolution)
{
SKY_GRID *grid;
SIN_THETA_SKY_GRID_PRIV *priv;
long i,j,k;
SKY_GRID_TYPE a,b;

/* set all up */
grid=do_alloc(1,sizeof(*grid));
grid->name="sin theta";
priv=do_alloc(1, sizeof(*priv));
priv->num_dec=ceil(M_PI/resolution);
priv->resolution=resolution;
priv->num_ra=do_alloc(priv->num_dec,sizeof(*priv->num_ra));
/* the total number of points is variable... */
grid->max_n_dec=priv->num_dec;
grid->max_n_ra=0;
grid->npoints=0;
for(i=0;i<priv->num_dec;i++){
	a=M_PI_2*(-1.0+(1.0+2.0*i)/priv->num_dec);
	priv->num_ra[i]=ceil(2.0*M_PI*cos(a)/resolution);
	/* make it always contain odd number of points.. makes plotting easier */
	priv->num_ra[i]|=1;
	if(priv->num_ra[i]>grid->max_n_ra)grid->max_n_ra=priv->num_ra[i];
	grid->npoints+=priv->num_ra[i];
	}

grid->latitude=do_alloc(grid->npoints, sizeof(SKY_GRID_TYPE));
grid->longitude=do_alloc(grid->npoints, sizeof(SKY_GRID_TYPE));

grid->max_x=2*(grid->max_n_ra>>1)+1;
grid->max_y=grid->max_n_dec;
grid->x=do_alloc(grid->npoints, sizeof(*(grid->x)));
grid->y=do_alloc(grid->npoints, sizeof(*(grid->y)));

grid->nbands_size=2000;
grid->nbands=0;
grid->band_name=do_alloc(grid->nbands_size, sizeof(*grid->band_name));

grid->band=do_alloc(grid->npoints, sizeof(*grid->band));
grid->band_f=do_alloc(grid->npoints, sizeof(*grid->band_f));

for(i=0;i<GRID_E_COUNT;i++)
	grid->e[i]=do_alloc(grid->npoints, sizeof(SKY_GRID_TYPE));
grid->grid_priv=priv;

/* fill in the coordinates */
k=0;
for(i=0;i<priv->num_dec;i++) {
	a=M_PI_2*(-1.0+(1.0+2.0*i)/priv->num_dec);
	for(j=0;j<priv->num_ra[i];j++){
		b=M_PI*(1.0+2.0*j)/priv->num_ra[i];
		grid->latitude[k]=a;
		grid->longitude[k]=b;
		grid->band[k]=-1;

		grid->x[k]=(grid->max_x>>1)-(priv->num_ra[i]>>1)+j;
		grid->y[k]=i;

		k++;
		}
	}
return grid;
}

long find_sin_theta_closest(SKY_GRID *grid, float RA, float DEC)
{
int i,j_start,j_stop,k, k_start, k_stop, best_i;
SIN_THETA_SKY_GRID_PRIV *priv=grid->grid_priv;
SKY_GRID_TYPE ds, best_ds=10.0;

k=floor((priv->num_dec*(DEC+M_PI_2))/M_PI+0.5);
k_start=k-1;
k_stop=k+1;
if(k_start<0)k_start=0;
if(k_stop>=priv->num_dec)k_stop=priv->num_dec-1;
if(k_stop<k_start)k_start=k_stop;

/* Find 3 lines worth of points to search .. Brute force, but reliable */
for(i=0,j_start=0;i<k_start;i++)j_start+=priv->num_ra[i];
for(j_stop=j_start;i<=k_stop;i++)j_stop+=priv->num_ra[i];

/* Search them for closest point */
best_i=-1;
for(i=j_start;i<j_stop;i++){
	ds=spherical_distance(RA, DEC, grid->longitude[i], grid->latitude[i]);
	if((best_i<0) || (ds<best_ds)){
		best_i=i;
		best_ds=ds;
		}
	}
return best_i;
}

/* This reduces the grid by eliminating band=-1 points */
SKY_GRID * reduced_grid(SKY_GRID *g)
{
SKY_GRID *grid;
REDUCED_SKY_GRID_PRIV *priv;
int i,k;

/* set all up */
grid=do_alloc(1,sizeof(*grid));
grid->npoints=0;
for(i=0;i<g->npoints;i++)if(g->band[i]>=0)grid->npoints++;
grid->max_n_dec=g->max_n_dec;
grid->max_n_ra=g->max_n_ra;
grid->name="reduced";
grid->latitude=do_alloc(grid->npoints, sizeof(SKY_GRID_TYPE));
grid->longitude=do_alloc(grid->npoints, sizeof(SKY_GRID_TYPE));

grid->max_x=g->max_x;
grid->max_y=g->max_y;
grid->x=do_alloc(grid->npoints, sizeof(*(grid->x)));
grid->y=do_alloc(grid->npoints, sizeof(*(grid->y)));

grid->nbands_size=2000;
grid->nbands=g->nbands;
grid->band_name=do_alloc(grid->nbands_size, sizeof(*grid->band_name));
for(i=0;i<grid->nbands;i++)grid->band_name[i]=strdup(g->band_name[i]);

grid->band=do_alloc(grid->npoints, sizeof(*grid->band));
grid->band_f=do_alloc(grid->npoints, sizeof(*grid->band_f));
for(i=0;i<GRID_E_COUNT;i++)
	grid->e[i]=do_alloc(grid->npoints, sizeof(SKY_GRID_TYPE));
priv=do_alloc(1, sizeof(*priv));
priv->original_index=do_alloc(grid->npoints, sizeof(*(priv->original_index)));
grid->grid_priv=priv;

/* fill in the points */
k=0;
for(i=0;i<g->npoints;i++) {
	if(g->band[i]<0)continue;
	priv->original_index[k]=i;
	grid->latitude[k]=g->latitude[i];
	grid->longitude[k]=g->longitude[i];
	grid->band[k]=g->band[i];
	grid->band_f[k]=g->band_f[i];
	grid->x[k]=g->x[i];
	grid->y[k]=g->y[i];
	k++;
	}

return grid;
}

void free_grid(SKY_GRID *grid)
{
long i;
/* free private fields */
if(!strcmp(grid->name,"sin theta")) {
	SIN_THETA_SKY_GRID_PRIV *priv;
	priv=grid->grid_priv;
	free(priv->num_ra);
	free(priv);	
	} else
if(!strcmp(grid->name,"reduced")) {
	REDUCED_SKY_GRID_PRIV *priv;
	priv=grid->grid_priv;
	free(priv->original_index);
	free(priv);	
	} else
if(!strcmp(grid->name,"plain rectangular")) {
	free(grid->grid_priv);
	} else
if(!strcmp(grid->name,"arcsin")) {
	free(grid->grid_priv);
	} else {
	if(grid->grid_priv!=NULL) {
		fprintf(stderr,"** Unknown grid type \"%s\", possible memory leak when freeing private structure\n", grid->name);
		fprintf(LOG,"** Unknown grid type \"%s\", possible memory leak when freeing private structure\n", grid->name);
		free(grid->grid_priv);
		}
	}
free(grid->latitude);
free(grid->longitude);
free(grid->x);
free(grid->y);
free(grid->band);
free(grid->band_f);
for(i=0;i<grid->nbands;i++)free(grid->band_name[i]);
free(grid->band_name);
for(i=0;i<GRID_E_COUNT;i++) {
	if(grid->e[i]!=NULL)free(grid->e[i]);
	}
free(grid);
}

void compute_list_map(SKY_SUPERGRID *sg)
{
long i,j,k,n;
sg->max_npatch=1;
for(i=0;i<sg->super_grid->npoints;i++){
	k=sg->reverse_map[i];
	j=sg->first_map[k];
	if(j==i)continue; /* do not mark very first item - we already know what it is */
	for(n=2;sg->list_map[j]>=0;j=sg->list_map[j],n++);
	sg->list_map[j]=i;
	if(n>sg->max_npatch)sg->max_npatch=n;
	}
}

void print_grid_statistics(FILE *file, char *prefix, SKY_GRID *grid)
{
long *count;
long masked_count;
int i;

count=aligned_alloca(grid->nbands*sizeof(*count));

for(i=0;i<grid->nbands;i++)count[i]=0;
masked_count=0;

for(i=0;i<grid->npoints;i++){
	if(grid->band[i]<0) {
		masked_count++;
		continue;
		}
	count[grid->band[i]]++;	
	}
for(i=0;i<grid->nbands;i++){
	fprintf(file, "%sgrid_points: %d \"%s\" %ld\n", prefix, i, grid->band_name[i], count[i]);
	}
fprintf(file, "%smasked_points: %ld\n", prefix, masked_count);
}

int add_band(SKY_GRID *sky_grid, char *name, int length)
{
int i;
if(length<0)length=strlen(name);
/* empty name matches nothing */
if(length==0)return (-1);

for(i=0;i<sky_grid->nbands;i++) {
	if(!strncmp(sky_grid->band_name[i], name, length) && strlen(sky_grid->band_name[i])==length)return i;
	}
/* band does not exist, add it */
if(sky_grid->nbands>=sky_grid->nbands_size) {
	/* one can put expanding code here, but I do not anticipate large number of bands, so make it static */
	fprintf(stderr, "*** ERROR: run out of band name array space\n");
	exit(-1);
	}

sky_grid->band_name[sky_grid->nbands]=do_alloc(length+1, 1);
memcpy(sky_grid->band_name[sky_grid->nbands], name, length);
sky_grid->band_name[sky_grid->nbands][length]=0;
sky_grid->nbands++;
return(sky_grid->nbands-1);
}


void angle_assign_bands(SKY_GRID *grid, int n_bands)
{
int i,k;
int *band_id;
char s[30];
SKY_GRID_TYPE angle, proj, x,y,z;

band_id=aligned_alloca(n_bands*sizeof(int));
for(i=0;i<n_bands;i++) {
	sprintf(s, "Angle_%d", i);
	band_id[i]=add_band(grid, s, -1);
	}

for(i=0;i<grid->npoints;i++) {
	/* convert into 3d */
	x=cos(grid->longitude[i])*cos(grid->latitude[i]);
	y=sin(grid->longitude[i])*cos(grid->latitude[i]);
	z=sin(grid->latitude[i]);

	proj=x*band_axis[0]+y*band_axis[1]+z*band_axis[2];
	if(proj<-1.0)proj=-1.0;
	if(proj>1.0)proj=1.0;
	
	angle=acosf(proj);
	
	k=floor((angle/M_PI)*n_bands);
	if(k<0)k=0;
	if(k>=n_bands)k=n_bands-1;
	grid->band[i]=band_id[k];
	grid->band_f[i]=band_id[k];
	}
}

void S_assign_bands(SKY_GRID *grid, int n_bands, double large_S, double spindown, double frequency)
{
int i,k;
double S;
int *band_id;
char s[30];
SKY_GRID_TYPE x,y,z;

band_id=aligned_alloca(n_bands*sizeof(int));
for(i=0;i<n_bands;i++) {
	sprintf(s, "S_%d", i);
	band_id[i]=add_band(grid, s, -1);
	}

for(i=0;i<grid->npoints;i++) {
	/* convert into 3d */
	x=cos(grid->longitude[i])*cos(grid->latitude[i]);
	y=sin(grid->longitude[i])*cos(grid->latitude[i]);
	z=sin(grid->latitude[i]);

 	S=spindown+
		band_axis_norm*frequency*(x*band_axis[0]+
		y*band_axis[1]+
		z*band_axis[2]);

	S=fabs(S);
	
	if(S>=large_S){
		grid->band[i]=0;
		grid->band_f[i]=0;
		continue;
		}
	k=n_bands-floor(S*(n_bands-1)/large_S)-1;
	if(k>=n_bands)k=n_bands-1;
	if(k<1)k=1;

	grid->band[i]=band_id[k];
	grid->band_f[i]=band_id[k];
	}
}

void mask_far_points(SKY_GRID *grid, SKY_GRID_TYPE ra, SKY_GRID_TYPE dec, SKY_GRID_TYPE radius)
{
int i;
SKY_GRID_TYPE ds;
for(i=0;i<grid->npoints;i++){
	ds=acos(sin(grid->latitude[i])*sin(dec)+
		cos(grid->latitude[i])*cos(dec)*
		cos(grid->longitude[i]-ra));
	if(ds>radius){
		grid->band[i]=-1;
		grid->band_f[i]=-1;
		}
	}
}

void mask_small_cos(SKY_GRID *grid, SKY_GRID_TYPE x, SKY_GRID_TYPE y, SKY_GRID_TYPE z, SKY_GRID_TYPE cos_level)
{
int i;
SKY_GRID_TYPE ds,x0,y0,z0;
ds=sqrt(x*x+y*y+z*z);
x0=x/ds;
y0=y/ds;
z0=z/ds;

for(i=0;i<grid->npoints;i++){
	ds=fabs(grid->e[0][i]*x0+
	   grid->e[1][i]*y0+
	   grid->e[2][i]*z0);
	if(ds<cos_level){
		grid->band[i]=-1;
		grid->band_f[i]=-1;
		}
	}
}

void mark_closest(SKY_GRID *sky_grid, int band_to, int band_from, SKY_GRID_TYPE ra, SKY_GRID_TYPE dec)
{
SKY_GRID_TYPE ds, ds_min;
int i,k;

k=-1;
ds_min=100;

for(i=0;i<sky_grid->npoints;i++) {
	if((band_from>=0) && sky_grid->band[i]!=band_from)continue;

	ds=acos(sin(sky_grid->latitude[i])*sin(dec)+
		cos(sky_grid->latitude[i])*cos(dec)*
		cos(sky_grid->longitude[i]-ra));

	if((ds<ds_min)) {
		ds_min=ds;
		k=i;
		}
	}

if( k>=0 ) {
	sky_grid->band[k]=band_to;
	sky_grid->band_f[k]=band_to;
	}
}

typedef struct {
	SKY_GRID *sky_grid;
	int band_to;
	int band_from;
	SKY_GRID_TYPE ra;
	SKY_GRID_TYPE dec;
	float ref_spindown;
	float weight_ratio_level;
	float bin_tolerance;
	float spindown_tolerance;
	int remainder;
	int divisor;
	} SIGNAL_SWEEP_DATA;

void signal_sweep_cruncher(int thread_id, SIGNAL_SWEEP_DATA *ssd)
{
int i;

for(i=0;i< ssd->sky_grid->npoints;i++) {
	if((i % ssd->divisor)!=ssd->remainder)continue;

	if((ssd->band_from>=0) && ssd->sky_grid->band[i]!=ssd->band_from)continue;


	if(effective_weight_ratio(ssd->sky_grid->longitude[i], ssd->sky_grid->latitude[i], ssd->ra, ssd->dec, ssd->ref_spindown, ssd->bin_tolerance, ssd->spindown_tolerance)> ssd->weight_ratio_level) {
		ssd->sky_grid->band[i]=ssd->band_to;
		ssd->sky_grid->band_f[i]=ssd->band_to;
		}
	}
}

void signal_sweep(SKY_GRID *sky_grid, int band_to, int band_from, SKY_GRID_TYPE ra, SKY_GRID_TYPE dec, float ref_spindown, float weight_ratio_level, float bin_tolerance, float spindown_tolerance)
{
int i;
int n_units;
SIGNAL_SWEEP_DATA *units;

n_units=get_max_threads()*2;
units=do_alloc(n_units, sizeof(*units));

for(i=0;i<n_units;i++) {
	units[i].remainder=i;
	units[i].divisor=n_units;
	units[i].sky_grid=sky_grid;
	units[i].band_to=band_to;
	units[i].band_from=band_from;
	units[i].ra=ra;
	units[i].dec=dec;
	units[i].ref_spindown=ref_spindown;
	units[i].weight_ratio_level=weight_ratio_level;
	units[i].bin_tolerance=bin_tolerance;
	units[i].spindown_tolerance=spindown_tolerance;
	
	submit_job(signal_sweep_cruncher, &(units[i]));
	}
while(do_single_job(-1));

wait_for_all_done();

free(units);
units=NULL;
}

void stationary_sweep(SKY_GRID *sky_grid, int band_to, int band_from, float weight_ratio_level, float tolerance)
{
int i;

for(i=0;i<sky_grid->npoints;i++) {
	if((band_from>=0) && sky_grid->band[i]!=band_from)continue;

	if(stationary_effective_weight_ratio(sky_grid->longitude[i], sky_grid->latitude[i], tolerance)>weight_ratio_level) {
		sky_grid->band[i]=band_to;
		sky_grid->band_f[i]=band_to;
		}
	}
}

void process_band_definition_line(SKY_GRID *sky_grid, char *line, int length)
{
int ai,aj, i;
SKY_GRID_TYPE ra, dec, radius, ds;
int band_to, band_from;

/* skip whitespace in the beginning */
while(((*line)==' ') || ((*line)=='\t'))line++;
/* skip comments */
if((*line)=='#')return;
/* skip empty lines */
if((*line)=='\n')return;
if((*line)=='\r')return;
if((*line)==0)return;

/* General format of the command:
	command band_to band_from [other_args]
*/

locate_arg(line, length, 1, &ai, &aj);
band_to=add_band(sky_grid, &(line[ai]), aj-ai);

locate_arg(line, length, 2, &ai, &aj);
band_from=add_band(sky_grid, &(line[ai]), aj-ai);

if(!strncasecmp(line, "disk", 4)) {
	int count;
	
	locate_arg(line, length, 3, &ai, &aj);
	sscanf(&(line[ai]), "%g", &ra);

	locate_arg(line, length, 4, &ai, &aj);
	sscanf(&(line[ai]), "%g", &dec);

	locate_arg(line, length, 5, &ai, &aj);
	sscanf(&(line[ai]), "%g", &radius);

	fprintf(stderr, "Marking disk (%d <- %d) around (%g, %g) with radius %g\n", band_to, band_from, ra, dec, radius);
	fprintf(LOG, "Marking disk (%d <- %d) around (%g, %g) with radius %g\n", band_to, band_from, ra, dec, radius);

	count=0;
	/* mark disk */
	for(i=0;i<sky_grid->npoints;i++){
		if((band_from>=0) && sky_grid->band[i]!=band_from)continue;

		ds=acos(sin(sky_grid->latitude[i])*sin(dec)+
			cos(sky_grid->latitude[i])*cos(dec)*
			cos(sky_grid->longitude[i]-ra));
		if(ds<radius) {
			count++;
			sky_grid->band[i]=band_to;
			sky_grid->band_f[i]=band_to;
			}
		}

	/* if the grid was too coarse and the radius too small just mark the closest point */
	if(count<1)mark_closest(sky_grid, band_to, band_from, ra, dec);
	} else
if(!strncasecmp(line, "band", 4)) {
	SKY_GRID_TYPE x0, y0, z0, level1, level2;
	locate_arg(line, length, 3, &ai, &aj);
	sscanf(&(line[ai]), "%g", &ra);

	locate_arg(line, length, 4, &ai, &aj);
	sscanf(&(line[ai]), "%g", &dec);

	locate_arg(line, length, 5, &ai, &aj);
	sscanf(&(line[ai]), "%g", &level1);

	locate_arg(line, length, 6, &ai, &aj);
	sscanf(&(line[ai]), "%g", &level2);

	fprintf(stderr, "Marking band (%d <- %d) around (%g, %g) with cos in [%g, %g]\n", band_to, band_from, ra, dec, level1, level2);
	fprintf(LOG, "Marking band (%d <- %d) around (%g, %g) with cos in [%g, %g]\n", band_to, band_from, ra, dec, level1, level2);

	x0=cos(ra)*sin(M_PI_2-dec);
	y0=sin(ra)*sin(M_PI_2-dec);
	z0=cos(M_PI_2-dec);

	for(i=0;i<sky_grid->npoints;i++) {
		if((band_from>=0) && sky_grid->band[i]!=band_from)continue;

		ds=cos(sky_grid->latitude[i])*cos(sky_grid->longitude[i])*x0+
			cos(sky_grid->latitude[i])*sin(sky_grid->longitude[i])*y0+
			sin(sky_grid->latitude[i])*z0;

		if((ds>level1) && (ds<=level2)) {
			sky_grid->band[i]=band_to;
			sky_grid->band_f[i]=band_to;
			}
		}
	} else 
if(!strncasecmp(line, "closest", 7)) {

	locate_arg(line, length, 3, &ai, &aj);
	sscanf(&(line[ai]), "%g", &ra);

	locate_arg(line, length, 4, &ai, &aj);
	sscanf(&(line[ai]), "%g", &dec);

	fprintf(stderr, "Marking point (%d <- %d) closest to (%g, %g)\n", band_to, band_from, ra, dec);
	fprintf(LOG, "Marking point (%d <- %d) closest to (%g, %g)\n", band_to, band_from, ra, dec);

	mark_closest(sky_grid, band_to, band_from, ra, dec);
	} else
if(!strncasecmp(line, "response", 8)) {
	float weight_ratio_level, bin_tolerance, spindown_tolerance;
	locate_arg(line, length, 3, &ai, &aj);
	sscanf(&(line[ai]), "%g", &ra);

	locate_arg(line, length, 4, &ai, &aj);
	sscanf(&(line[ai]), "%g", &dec);

	locate_arg(line, length, 5, &ai, &aj);
	sscanf(&(line[ai]), "%g", &weight_ratio_level);

	locate_arg(line, length, 6, &ai, &aj);
	sscanf(&(line[ai]), "%g", &bin_tolerance);

	locate_arg(line, length, 7, &ai, &aj);
	sscanf(&(line[ai]), "%g", &spindown_tolerance);

	fprintf(stderr, "Marking points (%d <- %d) swept by (%g, %g) weight_ratio=%g bin_width=%g spindown_width=%g\n", band_to, band_from, ra, dec, weight_ratio_level, bin_tolerance, spindown_tolerance);
	fprintf(LOG, "Marking points (%d <- %d) swept by (%g, %g) weight_ratio=%g bin_width=%g spindown_width=%g\n", band_to, band_from, ra, dec, weight_ratio_level, bin_tolerance, spindown_tolerance);

	signal_sweep(sky_grid, band_to, band_from, ra, dec, spindown, weight_ratio_level, bin_tolerance, spindown_tolerance);
	} else
if(!strncasecmp(line, "echo_response", 13)) {
	float weight_ratio_level, bin_tolerance, ref_spindown, spindown_tolerance;
	locate_arg(line, length, 3, &ai, &aj);
	sscanf(&(line[ai]), "%g", &ra);

	locate_arg(line, length, 4, &ai, &aj);
	sscanf(&(line[ai]), "%g", &dec);

	locate_arg(line, length, 5, &ai, &aj);
	sscanf(&(line[ai]), "%g", &ref_spindown);

	locate_arg(line, length, 6, &ai, &aj);
	sscanf(&(line[ai]), "%g", &weight_ratio_level);

	locate_arg(line, length, 7, &ai, &aj);
	sscanf(&(line[ai]), "%g", &bin_tolerance);

	locate_arg(line, length, 8, &ai, &aj);
	sscanf(&(line[ai]), "%g", &spindown_tolerance);

	fprintf(stderr, "Marking points (%d <- %d) swept by (%g, %g) spindown=%g weight_ratio=%g bin_width=%g spindown_width=%g\n", band_to, band_from, ra, dec, ref_spindown, weight_ratio_level, bin_tolerance, spindown_tolerance);
	fprintf(LOG, "Marking points (%d <- %d) swept by (%g, %g) spindown=%g weight_ratio=%g bin_width=%g spindown_width=%g\n", band_to, band_from, ra, dec, ref_spindown, weight_ratio_level, bin_tolerance, spindown_tolerance);

	signal_sweep(sky_grid, band_to, band_from, ra, dec, ref_spindown, weight_ratio_level, bin_tolerance, spindown_tolerance);
	} else
if(!strncasecmp(line, "line_response", 13)) {
	float weight_ratio_level, bin_tolerance;

	locate_arg(line, length, 3, &ai, &aj);
	sscanf(&(line[ai]), "%g", &weight_ratio_level);

	locate_arg(line, length, 4, &ai, &aj);
	sscanf(&(line[ai]), "%g", &bin_tolerance);

	fprintf(stderr, "Marking points (%d <- %d) swept by lines weight_ratio=%g bin_width=%g\n", band_to, band_from, weight_ratio_level, bin_tolerance);
	fprintf(LOG, "Marking points (%d <- %d) swept by lines weight_ratio=%g bin_width=%g\n", band_to, band_from, weight_ratio_level, bin_tolerance);

	stationary_sweep(sky_grid, band_to, band_from, weight_ratio_level, bin_tolerance);
	} else
	{
	fprintf(stderr, "*** UNKNOWN masking command \"%s\"\n", line);
	}
}

void process_marks(SKY_GRID *sky_grid, char *s, int length)
{
int ai, aj;
ai=0;
aj=0;
while(aj<length) {
	ai=aj;
	while(s[aj] && s[aj]!='\n' && (aj<length))aj++;
	process_band_definition_line(sky_grid, &(s[ai]), aj-ai);
	aj++;
	}
}

void propagate_far_points_to_super_grid(SKY_GRID *grid, SKY_SUPERGRID *super_grid)
{
long k, offset, pi;
for(pi=0;pi<grid->npoints;pi++){
	if(grid->band[pi]>=0)continue;
	for(k=0,offset=super_grid->first_map[pi];offset>=0;offset=super_grid->list_map[offset],k++){
		super_grid->super_grid->band[offset]=-1;
		super_grid->super_grid->band_f[offset]=-1;
		}
	}
}

void propagate_far_points_from_super_grid(SKY_GRID *grid, SKY_SUPERGRID *super_grid)
{
long k, offset, pi;
int nonzero;
for(pi=0;pi<grid->npoints;pi++){
	nonzero=0;
	for(k=0,offset=super_grid->first_map[pi];offset>=0;offset=super_grid->list_map[offset],k++){
		if(super_grid->super_grid->band[offset]>=0){
			nonzero=1;
			break;
			}
		}
	if(!nonzero){
		grid->band[pi]=-1;
		grid->band_f[pi]=-1;
		}
	}
}

SKY_SUPERGRID *make_rect_supergrid(SKY_GRID *grid, int ra_factor, int dec_factor)
{
SKY_SUPERGRID *sg;
RECT_SKY_GRID_PRIV *priv;
long i,j,k, ra_start, dec_start, di,dj;
if(strcmp(grid->name,"plain rectangular") && 
   strcmp(grid->name, "arcsin")){
   	fprintf(stderr,"** Internal error: cannot make rectangular supergrid from %s\n", grid->name);
	exit(-1);
   	}
priv=grid->grid_priv;

sg=do_alloc(1, sizeof(*sg));
if(!strcmp(grid->name,"arcsin")){
	sg->super_grid=make_arcsin_grid(priv->num_ra*ra_factor, priv->num_dec*dec_factor);
	} else {
	sg->super_grid=make_rect_grid(priv->num_ra*ra_factor, priv->num_dec*dec_factor);
	}

sg->subgrid_npoints=grid->npoints;
sg->first_map=do_alloc(grid->npoints, sizeof(*sg->first_map));
sg->reverse_map=do_alloc(sg->super_grid->npoints, sizeof(*sg->reverse_map));
sg->list_map=do_alloc(sg->super_grid->npoints, sizeof(*sg->list_map));
//sg->max_npatch=ra_factor*dec_factor;

/* clear the arrays */
for(i=0;i<grid->npoints;i++)sg->first_map[i]=-1;
for(i=0;i<sg->super_grid->npoints;i++){
	sg->reverse_map[i]=-1;
	sg->list_map[i]=-1;
	}

ra_start=(ra_factor+1)/2;
dec_start=(dec_factor+1)/2;

k=0;
for(i=ra_start;i<priv->num_ra*ra_factor;i+=ra_factor)
	for(j=dec_start;j<priv->num_dec*dec_factor;j+=dec_factor){
		sg->first_map[k]=i*priv->num_dec*dec_factor+j;
		for(di=0;di<ra_factor;di++)
			for(dj=0;dj<dec_factor;dj++){
				sg->reverse_map[(i+di-ra_start)*priv->num_dec*dec_factor+j+dj-dec_start]=k;
				}
		k++;
		}

compute_list_map(sg);
return sg;
}

SKY_SUPERGRID *make_targeted_rect_supergrid(SKY_GRID *grid, int factor)
{
SKY_SUPERGRID *sg;
TARGETED_RECT_SKY_GRID_PRIV *priv;
TARGETED_RECT_SKY_GRID_PRIV *sg_priv;
long i,j,k, ra_start, dec_start, pi, pj, pk, di,dj, best_di, best_dj, best_pk;
double ds, best_ds;
if(strcmp(grid->name,"targeted rectangular")){
   	fprintf(stderr,"** Internal error: cannot make rectangular supergrid from %s\n", grid->name);
	exit(-1);
   	}
priv=grid->grid_priv;

sg=do_alloc(1, sizeof(*sg));
sg->super_grid=make_targeted_rect_grid(priv->ra, priv->dec, priv->radius, (priv->num_dec-1)*factor+1);

sg->subgrid_npoints=grid->npoints;
sg->first_map=do_alloc(grid->npoints, sizeof(*sg->first_map));
sg->reverse_map=do_alloc(sg->super_grid->npoints, sizeof(*sg->reverse_map));
sg->list_map=do_alloc(sg->super_grid->npoints, sizeof(*sg->list_map));
//sg->max_npatch=ra_factor*dec_factor;
sg_priv=sg->super_grid->grid_priv;

/* clear the arrays */
for(i=0;i<grid->npoints;i++)sg->first_map[i]=-1;
for(i=0;i<sg->super_grid->npoints;i++){
	sg->reverse_map[i]=-1;
	sg->list_map[i]=-1;
	}

ra_start=(factor)/2;
dec_start=(factor)/2;

for(i=0;i<sg_priv->num_dec;i++)
	for(j=0;j<sg_priv->num_dec;j++) {
		pi=i/factor;
		pj=j/factor;

		k=i*sg_priv->num_dec+j;
		
		best_ds=100.0;
		best_di=-2;
		best_dj=-2;
		best_pk=-1;

		for(di=-2;di<=2;di++)
			for(dj=-2;dj<=2;dj++) {
				if(pi+di<0)continue;
				if(pi+di>=priv->num_dec)continue;
				if(pj+dj<0)continue;
				if(pj+dj>=priv->num_dec)continue;
				
				pk=(pi+di)*priv->num_dec+pj+dj;
				
				ds=spherical_distance(grid->longitude[pk], grid->latitude[pk],
					sg->super_grid->longitude[k], sg->super_grid->latitude[k]);
				if(ds<best_ds) {
					best_ds=ds;
					best_pk=pk;
					best_di=di;
					best_dj=dj;
					}
				}
		if(best_pk<0) {
			fprintf(stderr, "***** INTERNAL ERROR: targeted grid did not find closest patch\n");
			exit(-1);
			}
		sg->first_map[best_pk]=k;
		sg->reverse_map[k]=best_pk;
		}

compute_list_map(sg);
/*for(i=0;i<sg->super_grid->npoints;i++){
	 if(i % (sg_priv->num_dec)==0)fprintf(stderr, "\n");
	 fprintf(stderr, "%d ", sg->list_map[i]);
	}
for(i=0;i<grid->npoints;i++){
	 if(i % (priv->num_dec)==0)fprintf(stderr, "\n");
	 fprintf(stderr, "%d ", sg->first_map[i]);
	}
for(i=0;i<sg->super_grid->npoints;i++){
	 if(i % (sg_priv->num_dec)==0)fprintf(stderr, "\n");
	 fprintf(stderr, "%d ", sg->reverse_map[i]);
	}
for(i=0;i<sg->super_grid->npoints;i++){
	 if(i % (priv->num_dec*factor)==0)fprintf(stderr, "\n");
	 fprintf(stderr, "(%.2f, %.2f) ", sg->super_grid->longitude[i], sg->super_grid->latitude[i]);
	}
fprintf(stderr, "\n");*/
return sg;
}

SKY_SUPERGRID *make_sin_theta_supergrid(SKY_GRID *grid, int factor)
{
SKY_SUPERGRID *sg;
SIN_THETA_SKY_GRID_PRIV *priv;
int i,j,k, pi, pk, ra_pk;
SKY_GRID_TYPE ds, best_ds;
if(strcmp(grid->name,"sin theta")){
   	fprintf(stderr,"** Internal error: cannot make sin theta supergrid from %s\n", grid->name);
	exit(-1);
   	}
priv=grid->grid_priv;

sg=do_alloc(1, sizeof(*sg));
sg->super_grid=make_sin_theta_grid(priv->resolution/factor);

sg->subgrid_npoints=grid->npoints;
sg->first_map=do_alloc(grid->npoints, sizeof(*sg->first_map));
sg->reverse_map=do_alloc(sg->super_grid->npoints, sizeof(*sg->reverse_map));
sg->list_map=do_alloc(sg->super_grid->npoints, sizeof(*sg->list_map));

fprintf(stderr,"npoints=%d super grid npoints=%d\n", grid->npoints, sg->super_grid->npoints);
/* clear the arrays */
for(i=0;i<grid->npoints;i++)sg->first_map[i]=-1;
for(i=0;i<sg->super_grid->npoints;i++){
	sg->reverse_map[i]=-1;
	sg->list_map[i]=-1;
	}

k=0;
pk=0;
ra_pk=0;
pi=0;
sg->first_map[0]=0;
sg->first_map[grid->npoints-1]=sg->super_grid->npoints-1;
sg->reverse_map[sg->super_grid->npoints-1]=grid->npoints-1;
sg->reverse_map[0]=0;
for(k=1;k<sg->super_grid->npoints-1;k++){
	#if 0
	fprintf(stderr, "patch=(%.2f,%.2f) super=(%.2f,%.2f)\n", 
		grid->longitude[pk],grid->latitude[pk],
		sg->super_grid->longitude[k], sg->super_grid->latitude[k]
		);
	#endif
	
	#if 0  /* older and really fast way.. is it right ?? */
	if((sg->super_grid->longitude[k]<sg->super_grid->longitude[k-1])){
			/* crossing RA=0 boundary */
	//		fprintf(stderr, "@");
			//pk=pk-priv->num_ra[pi]+1;
			pk=ra_pk;
			} else
	if((pk+1<grid->npoints)&& (grid->longitude[pk+1]>=grid->longitude[pk]) && 
		(sg->super_grid->longitude[k]-grid->longitude[pk]>=
			grid->longitude[pk+1]-sg->super_grid->longitude[k])){
	//		fprintf(stderr, "+");
			pk++;
			}
	if((pi+1<priv->num_dec) && (sg->super_grid->latitude[k]-grid->latitude[pk]>=
		grid->latitude[ra_pk+priv->num_ra[pi]]-sg->super_grid->latitude[k])){
	//		fprintf(stderr, "#");
			ra_pk+=priv->num_ra[pi];
			pk=ra_pk;
			pi++;
			}
	//fprintf(stderr, "k=%d pk=%d pi=%d\n", k, pk, pi);
	sg->reverse_map[k]=pk;	
	sg->first_map[pk]=k; /* not the most efficient, but it works */

        #else /* a good deal slower, but much surer */
	
	if((sg->super_grid->longitude[k]<sg->super_grid->longitude[k-1])){
			/* crossing RA=0 boundary */
	//		fprintf(stderr, "@");
			//pk=pk-priv->num_ra[pi]+1;
			pk=ra_pk;
			}
			else
	if((pk+1<grid->npoints)&& (grid->longitude[pk+1]>=grid->longitude[pk]) && 
		(sg->super_grid->longitude[k]-grid->longitude[pk]>=
			grid->longitude[pk+1]-sg->super_grid->longitude[k])){
	//		fprintf(stderr, "+");
			pk++;
			}
	if((pi+1<priv->num_dec) && (sg->super_grid->latitude[k]-grid->latitude[pk]>=
		grid->latitude[ra_pk+priv->num_ra[pi]]-sg->super_grid->latitude[k])){
	//		fprintf(stderr, "#");
			ra_pk+=priv->num_ra[pi];
			pk=ra_pk;
			pi++;
			}
	

	best_ds=10.0;
	j=pk;
	
	i=pk-2*grid->max_n_ra-1;
	if(i<0)i=0;
	for(;(i<(pk+2*grid->max_n_ra+1)) && (i<grid->npoints);i++){
		/* Try approximate comparison first */
		ds=fabs(grid->longitude[i]-sg->super_grid->longitude[k]);
		/* check that we are not far enough in RA */
		/* The (ds<1.0) is to check that we are not jumping 2*PI */
		if((cos(grid->latitude[i])*ds>best_ds)&&(ds<6.0))continue;
		ds=spherical_distance(grid->longitude[i], grid->latitude[i],
				sg->super_grid->longitude[k], sg->super_grid->latitude[k]);
		if(ds<best_ds){
			j=i;
			best_ds=ds;
			}
		}

	//if(pk!=j)fprintf(stderr, "k=%d pk=%d pi=%d j=%d  %d\n", k, pk, pi, j, pk-j);
	sg->reverse_map[k]=j;	
	sg->first_map[j]=k; /* not the most efficient, but it works */

	#endif

	}
compute_list_map(sg);
return sg;
}

SKY_SUPERGRID * reduced_supergrid(SKY_SUPERGRID *sg0)
{
SKY_SUPERGRID *sg;
REDUCED_SKY_GRID_PRIV *priv;
int i,k;
int *new_index;

sg=do_alloc(1, sizeof(*sg));
sg->super_grid=reduced_grid(sg0->super_grid);
priv=sg->super_grid->grid_priv;

sg->subgrid_npoints=sg0->subgrid_npoints;
sg->first_map=do_alloc(sg->subgrid_npoints, sizeof(*sg->first_map));
sg->reverse_map=do_alloc(sg->super_grid->npoints, sizeof(*sg->reverse_map));
sg->list_map=do_alloc(sg->super_grid->npoints, sizeof(*sg->list_map));
sg->max_npatch=sg0->max_npatch;

fprintf(stderr,"reduced npoints=%d to npoints=%d\n", sg0->super_grid->npoints, sg->super_grid->npoints);

new_index=do_alloc(sg0->super_grid->npoints, sizeof(*new_index));
for(i=0;i<sg0->super_grid->npoints;i++)new_index[i]=-1;
for(i=0;i<sg->super_grid->npoints;i++)new_index[priv->original_index[i]]=i;

for(i=0;i<sg->subgrid_npoints;i++) {
	k=sg0->first_map[i];
	if(k<0) {
		sg->first_map[i]=k;
		} else {
		while( (sg0->super_grid->band[k]<0) && (sg0->list_map[k]>=0))k=sg0->list_map[k];
		sg->first_map[i]=new_index[k]; /* will be -1 if band[k]<0 */
		}

	}

for(i=0;i<sg->super_grid->npoints;i++) {
	sg->reverse_map[i]=sg0->reverse_map[priv->original_index[i]];

	k=sg0->list_map[priv->original_index[i]];
	while( (k>=0) && (sg0->super_grid->band[k]<0))k=sg0->list_map[k];
	if(k<0)sg->list_map[i]=-1;
		else sg->list_map[i]=new_index[k];
	}
free(new_index);
/*for(i=0;i<sg->super_grid->npoints;i++){
	 //if(i % (priv->num_dec*factor)==0)fprintf(stderr, "\n");
	 fprintf(stderr, "%d ", sg->list_map[i]);
	}
fprintf(stderr, "\n");
fprintf(stderr, "\n");
for(i=0;i<sg->subgrid_npoints;i++){
	 //if(i % (priv->num_dec)==0)fprintf(stderr, "\n");
	 fprintf(stderr, "%d ", sg->first_map[i]);
	}
fprintf(stderr, "\n");
*/
return(sg);
}

void free_supergrid(SKY_SUPERGRID *sg)
{
free(sg->first_map);
free(sg->list_map);
free(sg->reverse_map);
free_grid(sg->super_grid);
sg->first_map=NULL;
sg->list_map=NULL;
sg->reverse_map=NULL;
free(sg);
}


void rotate_xz(SKY_GRID_TYPE RA_in, SKY_GRID_TYPE DEC_in, 
			SKY_GRID_TYPE * RA_out, SKY_GRID_TYPE * DEC_out, 
			SKY_GRID_TYPE angle)
{
SKY_GRID_TYPE x,y,z,x2,y2,z2;

/* convert into 3d */
x=cos(RA_in)*cos(DEC_in);
y=sin(RA_in)*cos(DEC_in);
z=sin(DEC_in);

x2=cos(angle)*x+sin(angle)*z;
y2=y;
z2=-sin(angle)*x+cos(angle)*z;

*DEC_out=atan2f(z2, sqrt(x2*x2+y2*y2));
*RA_out=atan2f(y2, x2);
if(*RA_out <0) *RA_out+=2*M_PI;
//fprintf(stderr,"%f %f --> %f %f\n", RA_in, DEC_in, *RA_out, *DEC_out);
}

void rotate_xy(SKY_GRID_TYPE RA_in, SKY_GRID_TYPE DEC_in, 
			SKY_GRID_TYPE * RA_out, SKY_GRID_TYPE * DEC_out, 
			SKY_GRID_TYPE angle)
{
SKY_GRID_TYPE x,y,z,x2,y2,z2;

/* convert into 3d */
x=cos(RA_in)*cos(DEC_in);
y=sin(RA_in)*cos(DEC_in);
z=sin(DEC_in);

x2=cos(angle)*x-sin(angle)*y;
y2=sin(angle)*x+cos(angle)*y;
z2=z;

*DEC_out=atan2f(z2, sqrt(x2*x2+y2*y2));
*RA_out=atan2f(y2, x2);
if(*RA_out <0) *RA_out+=2*M_PI;
//fprintf(stderr,"%f %f --> %f %f\n", RA_in, DEC_in, *RA_out, *DEC_out);
}

void rotate_grid_xz(SKY_GRID *grid, SKY_GRID_TYPE angle)
{
long i;
for(i=0;i<grid->npoints;i++){
	rotate_xz((grid->longitude[i]), (grid->latitude[i]),
		&(grid->longitude[i]), &(grid->latitude[i]),
		angle);
	}
}

void rotate_grid_xy(SKY_GRID *grid, SKY_GRID_TYPE angle)
{
long i;
for(i=0;i<grid->npoints;i++){
	rotate_xy((grid->longitude[i]), (grid->latitude[i]),
		&(grid->longitude[i]), &(grid->latitude[i]),
		angle);
	}
}
