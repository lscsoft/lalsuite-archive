#ifndef __POLARIZATION_H__
#define __POLARIZATION_H__

#include "grid.h"
#include "global.h"

typedef struct {
	float  plus_proj;
	float  cross_proj;
	
	SKY_GRID_TYPE *AM_coeffs;
	float *patch_CutOff;
	
	char *name;
	/* these arrays have stored_fine_bins*useful_bins entries */
	SUM_TYPE *fine_grid_sum;
	SUM_TYPE *fine_grid_sq_sum;
	SUM_TYPE *fine_grid_weight;
	COUNT_TYPE *fine_grid_count;
	
	/* results */
	/* these arrays have find_grid->npoints entries */
	SUM_TYPE *total_weight;
	SUM_TYPE *max_sub_weight;
	COUNT_TYPE *total_count;
	SUM_TYPE *max_dx;
	SUM_TYPE *M_map;
	SUM_TYPE *S_map;
	SUM_TYPE *max_upper_limit;
	SUM_TYPE *max_lower_limit;
	SUM_TYPE *freq_map;
	SUM_TYPE *cor1;
	SUM_TYPE *cor2;
	} POLARIZATION;

static float inline AM_response(int segment, SKY_GRID *grid, int point, float *coeff)
{
float a;
int ii;
a=0.0;
for(ii=0;ii<GRID_FIT_COUNT;ii++)
	a+=coeff[segment*GRID_FIT_COUNT+ii]*grid->e[ii+GRID_FIT_START][point];
#if 0 /* just for testing */
return a;
#else
return (a*a);
#endif
}


void init_polarizations(void);
void allocate_polarization_arrays(void);
void clear_accumulation_arrays(void);



#endif
