#ifndef __POLARIZATION_H__
#define __POLARIZATION_H__

#include "grid.h"
#include "global.h"

typedef struct S_POLARIZATION {
	float orientation;
	float  plus_proj;
	float  cross_proj;
	struct S_POLARIZATION *conjugate; /* plus for cross and cross for plus */
	/* For linear polarizations plus_factor is always 1.0 and
	   cross factor is always 0.0 */
	float plus_factor;
	float cross_factor;
	
	SKY_GRID_TYPE *AM_coeffs;
	float *patch_CutOff;
	
	char *name;
	/* these arrays have stored_fine_bins*useful_bins entries */
	SUM_TYPE *fine_grid_sum;
	SUM_TYPE *fine_grid_sq_sum;
	SUM_TYPE *fine_grid_weight;
	COUNT_TYPE *fine_grid_count;
	
	/* results - sky maps, projection along frequency */
	/* these arrays have fine_grid->npoints entries */
	struct {
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
		SUM_TYPE *ks_test;
		SUM_TYPE *ks_count;
		
		/* universality coefficients - these are needed to 
		   translate from several sampled linear polarizations 
		   to limits on circular and arbitrary linear polarizations.
		   
		   They can also be used to set limits on arbitrary pulsar
		   signal 
		   
		   WARNING: they are *not* applicable for signals with many
		   vetoed bins
		   
		   */
		   
		SUM_TYPE *beta1;
		SUM_TYPE *beta2;
		} skymap;
	
	/* results - spectral plots, projection along sky bands */
	/* these arrays have useful_bins*args_info.dec_bands_arg entries */
	struct {
		SUM_TYPE *max_upper_limit;
		SUM_TYPE *ul_dec;
		SUM_TYPE *ul_ra;
		SUM_TYPE *max_dx;
		SUM_TYPE *dx_dec;
		SUM_TYPE *dx_ra;
		SUM_TYPE *max_mask_ratio;
		} spectral_plot;
	} POLARIZATION;

extern int no_am_response;

static float inline F_plus(int segment, SKY_GRID *grid, int point, float *coeff)
{
float a;
int ii;
if(no_am_response)return 1.0;
a=0.0;
for(ii=0;ii<GRID_FIT_COUNT;ii++)
	a+=coeff[segment*GRID_FIT_COUNT+ii]*grid->e[ii+GRID_FIT_START][point];
return a;
}

static float inline AM_response(int segment, SKY_GRID *grid, int point, float *coeff)
{
float a;
a=F_plus(segment, grid, point, coeff);
return (a*a);
}


void init_polarizations(void);
void allocate_polarization_arrays(void);
void clear_accumulation_arrays(void);

#endif
