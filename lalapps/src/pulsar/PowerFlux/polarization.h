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
	
	SKY_GRID_TYPE * __restrict__ AM_coeffs;
	float *patch_CutOff;
	
	char *name;

	} POLARIZATION;

typedef struct {
	/* intermediate results: this array has stored_fine_bins entries */
	SUM_TYPE *total_weight;
	/* intermediate results: these arrays have stored_fine_bins*useful_bins entries */
	SUM_TYPE *fine_grid_sum;
	SUM_TYPE *fine_grid_sq_sum;
	SUM_TYPE *fine_grid_weight;
	COUNT_TYPE *fine_grid_count;

		/* double variants of the same structure - for intermediate accumulation */
	double *total_weight_d;
	double *fine_grid_sum_d;

	} ACCUMULATION_ARRAYS;


typedef struct {
	char *name;

	float orientation;
	float plus_proj;
	float cross_proj;

	float plus_factor;
	float cross_factor;

	int conjugate;

	/* results - sky maps, projection along frequency */
	/* these arrays have fine_grid->npoints entries */
	struct {
		SUM_TYPE *total_weight;
		SUM_TYPE *max_sub_weight;
		SUM_TYPE *total_count;
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

	} POLARIZATION_RESULTS;

extern int no_am_response;

inline static void precompute_am_constants(float *e, float ra, float dec)
{
float e2, e3, e4, e5;

e2=cos(M_PI_2-dec);
e3=sin(M_PI_2-dec);
e4=cos(ra);
e5=sin(ra);
/* unit vector */
e[0]=e3*e4;
e[1]=e3*e5;
e[2]=e2;
/* other useful values */
e[3]=e3;
e[4]=e4;
e[5]=e5;
/* these values are needed for regression of plus and cross */
e[6]=e4*e5;
e[7]=e3*e4;
e[8]=e3*e5;

e[9]=e2*e2*e4;
e[10]=e2*e2*e5;
e[11]=e3*e3*e4;
e[12]=e3*e3*e5;
e[13]=e2*e3*e4;
e[14]=e2*e3*e5;

e[15]=e2*e4*e4;
e[16]=e2*e5*e5;
e[17]=e3*e4*e4;
e[18]=e3*e5*e5;
e[19]=e2*e4*e5;

e[20]=e2*e2*e4*e4;
e[21]=e2*e2*e5*e5;
e[22]=e3*e3*e4*e4;
e[23]=e3*e3*e5*e5;

e[24]=e2*e2*e4*e5;
e[25]=e3*e3*e4*e5;
}

static float inline F_plus_coeff(int segment, float *e, float *coeff)
{
float a;
int ii;
if(no_am_response)return 1.0;
a=0.0;
for(ii=0;ii<GRID_FIT_COUNT;ii++)
	a+=coeff[segment*GRID_FIT_COUNT+ii]*e[ii+GRID_FIT_START];
return a;
}

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


void init_polarizations0(void);
void init_polarizations1(POLARIZATION *polarizations, SKY_GRID_TYPE *AM_coeffs_plus, SKY_GRID_TYPE *AM_coeffs_cross, long AM_coeffs_size);
void allocate_polarization_arrays(void);
void free_polarization_arrays(void);
void clear_polarization_arrays(void);

ACCUMULATION_ARRAYS *new_accumulation_arrays(void);
void free_accumulation_arrays(ACCUMULATION_ARRAYS *);
ACCUMULATION_ARRAYS *get_thread_accumulation_arrays(int thread_id);
void clear_accumulation_arrays(ACCUMULATION_ARRAYS *);
void update_d_accumulation_arrays(ACCUMULATION_ARRAYS *r);
void finalize_accumulation_arrays(ACCUMULATION_ARRAYS *r);

/* How to use:

      ar=new_accumulation_arrays()

      loop {
          clear_accumulation_arrays()
               loop {
                   compute and add 
                   periodically update_d_accumulation_arrays()
                    }
          finalize_accumulation_arrays()
          Use results
          }
     free_accumulation_arrays() 
*/


#endif
