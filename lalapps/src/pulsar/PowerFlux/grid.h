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
	long npoints; /* total number of points in this grid */
	long max_n_ra;
	long max_n_dec;
	char *name; 
	SKY_GRID_TYPE *latitude;
	SKY_GRID_TYPE *longitude;
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
	long num_ra;
	long num_dec;
	} RECT_SKY_GRID_PRIV;	

typedef struct {
	SKY_GRID_TYPE resolution;
	long num_dec;
	long *num_ra;
	} SIN_THETA_SKY_GRID_PRIV;

typedef struct {
	SKY_GRID *super_grid;		/* larger grid */
	long *first_map;      /* these are indices of subgrid points */
	long *reverse_map;    /* reverse map from grid to nearest subgrid point */
	long *list_map;    /* these indices form lists of nearest neighbours */
	long max_npatch;    /* maximum number of fine points in a patch */
	} SKY_SUPERGRID;

SKY_GRID *make_arcsin_grid(long num_ra, long num_dec);
SKY_GRID *make_rect_grid(long num_ra, long num_dec);
SKY_GRID *make_sin_theta_grid(SKY_GRID_TYPE resolution);
void free_grid(SKY_GRID *grid);

SKY_SUPERGRID *make_rect_supergrid(SKY_GRID *grid, int ra_factor, int dec_factor);
SKY_SUPERGRID *make_sin_theta_supergrid(SKY_GRID *grid, int factor);

void assign_bands(SKY_GRID *grid, int n_bands);
void mask_far_points(SKY_GRID *grid, SKY_GRID_TYPE ra, SKY_GRID_TYPE dec, SKY_GRID_TYPE radius);
void mask_small_cos(SKY_GRID *grid, SKY_GRID_TYPE x, SKY_GRID_TYPE y, SKY_GRID_TYPE z, SKY_GRID_TYPE cos_level);

void propagate_far_points_to_super_grid(SKY_GRID *grid, SKY_SUPERGRID *super_grid);
void propagate_far_points_from_super_grid(SKY_GRID *grid, SKY_SUPERGRID *super_grid);

void rotate_grid_xz(SKY_GRID *grid, SKY_GRID_TYPE angle);
void rotate_grid_xy(SKY_GRID *grid, SKY_GRID_TYPE angle);

#endif
