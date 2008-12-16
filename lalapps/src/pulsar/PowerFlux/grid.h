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

#ifndef __GRID_H__
#define __GRID_H__

typedef float SKY_GRID_TYPE;

SKY_GRID_TYPE spherical_distance(SKY_GRID_TYPE ra0, SKY_GRID_TYPE dec0,
			  SKY_GRID_TYPE ra1, SKY_GRID_TYPE dec1);

/* this saves a call to acos() and is different
   from above by a monotonic function
   It can be used for distance comparisons */
SKY_GRID_TYPE fast_spherical_distance(SKY_GRID_TYPE ra0, SKY_GRID_TYPE dec0,
			  SKY_GRID_TYPE ra1, SKY_GRID_TYPE dec1);


#define GRID_E_COUNT	26
#define GRID_FIT_START	2
#define GRID_FIT_COUNT  24

typedef struct {
	int npoints; /* total number of points in this grid */
	int max_n_ra;
	int max_n_dec;
	char *name; 
	SKY_GRID_TYPE *latitude;
	SKY_GRID_TYPE *longitude;

	char **band_name;
	int nbands;
	int nbands_size;

	/* for plotting - coordinates of points to mark */
	int max_x;
	int max_y;
	int *x;
	int *y;

	int *band;
	float *band_f;  /* just a convenience data.. */
	SKY_GRID_TYPE *e[GRID_E_COUNT];  /* 3d coordinates: 
				0,1,2 - coordinates of the unit vector
				2  - cos(M_PI_2-latitude)
				
				others - precomputed values:
				3  - sin(M_PI_2-latitude)
				4  - cos(longitude)
				5  - sin(longitude)
				
				2 through 22 are used in regression of
				plus and cross response values
				*/
	void *grid_priv; /* private */
	} SKY_GRID;

typedef struct {
	int num_ra;
	int num_dec;
	} RECT_SKY_GRID_PRIV;	

typedef struct {
	SKY_GRID_TYPE resolution;
	int num_dec;
	int *num_ra;
	} SIN_THETA_SKY_GRID_PRIV;

typedef struct {
	int *original_index;
	} REDUCED_SKY_GRID_PRIV;	

typedef struct {
	SKY_GRID *super_grid;		/* larger grid */
	int subgrid_npoints;
	int *first_map;      /* these are indices of subgrid points */
	int *reverse_map;    /* reverse map from grid to nearest subgrid point */
	int *list_map;    /* these indices form lists of nearest neighbours */
	int max_npatch;    /* maximum number of fine points in a patch */
	} SKY_SUPERGRID;

SKY_GRID_TYPE spherical_distance(SKY_GRID_TYPE ra0, SKY_GRID_TYPE dec0,
			  SKY_GRID_TYPE ra1, SKY_GRID_TYPE dec1);
			  
SKY_GRID_TYPE fast_spherical_distance(SKY_GRID_TYPE ra0, SKY_GRID_TYPE dec0,
			  SKY_GRID_TYPE ra1, SKY_GRID_TYPE dec1);

SKY_GRID *make_arcsin_grid(long num_ra, long num_dec);
SKY_GRID *make_rect_grid(long num_ra, long num_dec);
SKY_GRID *make_sin_theta_grid(SKY_GRID_TYPE resolution);
/* This reduces the grid by eliminating band=-1 points */
SKY_GRID * reduced_grid(SKY_GRID *g);

void precompute_values(SKY_GRID *grid);
void free_values(SKY_GRID *grid);

void free_grid(SKY_GRID *grid);

long find_sin_theta_closest(SKY_GRID *grid, float RA, float DEC);

SKY_SUPERGRID *make_rect_supergrid(SKY_GRID *grid, int ra_factor, int dec_factor);
SKY_SUPERGRID *make_sin_theta_supergrid(SKY_GRID *grid, int factor);
SKY_SUPERGRID * reduced_supergrid(SKY_SUPERGRID *sg0);
void free_supergrid(SKY_SUPERGRID *sg);

int add_band(SKY_GRID *sky_grid, char *name, int length);
void print_grid_statistics(FILE *file, char *prefix, SKY_GRID *grid);
void angle_assign_bands(SKY_GRID *grid, int n_bands);
void S_assign_bands(SKY_GRID *grid, int n_bands, double large_S, double s, double f);
void mask_far_points(SKY_GRID *grid, SKY_GRID_TYPE ra, SKY_GRID_TYPE dec, SKY_GRID_TYPE radius);
void mask_small_cos(SKY_GRID *grid, SKY_GRID_TYPE x, SKY_GRID_TYPE y, SKY_GRID_TYPE z, SKY_GRID_TYPE cos_level);
void process_marks(SKY_GRID *sky_grid, char *s, int length);

void propagate_far_points_to_super_grid(SKY_GRID *grid, SKY_SUPERGRID *super_grid);
void propagate_far_points_from_super_grid(SKY_GRID *grid, SKY_SUPERGRID *super_grid);

void rotate_grid_xz(SKY_GRID *grid, SKY_GRID_TYPE angle);
void rotate_grid_xy(SKY_GRID *grid, SKY_GRID_TYPE angle);

#endif
