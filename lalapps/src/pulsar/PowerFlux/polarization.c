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
#include <math.h>
#include <string.h>

#include "global.h"
#include "polarization.h"
#include "cmdline.h"
#include "jobs.h"


extern struct gengetopt_args_info args_info;
extern FILE *LOG;

extern SKY_GRID *fine_grid, *patch_grid;

extern int stored_fine_bins;
extern int useful_bins;

int ntotal_polarizations=-1, nlinear_polarizations=-1;

POLARIZATION_RESULTS *polarization_results=NULL;

int max_cruncher_threads=-2;
ACCUMULATION_ARRAYS **aar=NULL;


void init_polarizations0()
{
int i, nderived=1;
float a;


nlinear_polarizations=args_info.nlinear_polarizations_arg;;
ntotal_polarizations=nlinear_polarizations+nderived;

fprintf(LOG,"nlinear_polarizations: %d\n", nlinear_polarizations);
fprintf(LOG,"ntotal_polarizations : %d\n", ntotal_polarizations);
if(nlinear_polarizations<2){
	fprintf(stderr,"*** ERROR: number of polarizations is less than 2, aborting\n");
	exit(-1);
	}
if(nlinear_polarizations & 1){
	fprintf(stderr,"*** ERROR: number of polarizations is not even, aborting\n");
	exit(-1);
	}

polarization_results=do_alloc(ntotal_polarizations, sizeof(*polarization_results));

polarization_results[0].orientation=0;
polarization_results[0].name="plus";
polarization_results[0].plus_proj=1.0;
polarization_results[0].cross_proj=0.0;
polarization_results[0].plus_factor=1.0;
polarization_results[0].cross_factor=0.0;
polarization_results[0].conjugate=1;
fprintf(stderr,"\t%s 0.0\n",polarization_results[0].name);

polarization_results[1].orientation=M_PI/4.0;
polarization_results[1].name="cross";
polarization_results[1].plus_proj=0.0;
polarization_results[1].cross_proj=1.0;
polarization_results[1].plus_factor=1.0;
polarization_results[1].cross_factor=0.0;
polarization_results[1].conjugate=0;
fprintf(stderr,"\t%s %g\n",polarization_results[1].name, M_PI/4.0);

for(i=2;i<nlinear_polarizations;i++){
	polarization_results[i].name=do_alloc(16,sizeof(char));
	a=(i>>1)*M_PI/(2.0*nlinear_polarizations);
	if(i & 1){
		a+=M_PI/4.0;
		snprintf(polarization_results[i].name,16,"pi_%d_%d",((i+nlinear_polarizations)>>1), 2*nlinear_polarizations);
		} else {
		snprintf(polarization_results[i].name,16,"pi_%d_%d", (i>>1), 2*nlinear_polarizations);
		}
	polarization_results[i].orientation=a;
	polarization_results[i].plus_proj=cos(2*a);
	polarization_results[i].cross_proj=sin(2*a);
	polarization_results[i].plus_factor=1.0;
	polarization_results[i].cross_factor=0.0;
	polarization_results[i].conjugate=i^1;
	fprintf(stderr,"\t%s a=%g\n",polarization_results[i].name, a);
	}

polarization_results[nlinear_polarizations+0].orientation=0;
polarization_results[nlinear_polarizations+0].name="circular";
polarization_results[nlinear_polarizations+0].plus_proj=1.0;
polarization_results[nlinear_polarizations+0].cross_proj=0.0;
polarization_results[nlinear_polarizations+0].plus_factor=1.0;
polarization_results[nlinear_polarizations+0].cross_factor=1.0;
polarization_results[nlinear_polarizations+0].conjugate=1;
fprintf(stderr,"\t%s %f %f\n",polarization_results[nlinear_polarizations+0].name,
	polarization_results[nlinear_polarizations+0].plus_factor, 
	polarization_results[nlinear_polarizations+0].cross_factor);

}

void init_polarizations1(POLARIZATION *polarizations, SKY_GRID_TYPE *AM_coeffs_plus, SKY_GRID_TYPE *AM_coeffs_cross, long AM_coeffs_size)
{
int i, k;

for(i=0;i<ntotal_polarizations;i++) {
	fprintf(stderr,"\t%s\n",polarization_results[i].name);
	polarizations[i].orientation=polarization_results[i].orientation;
	polarizations[i].name=polarization_results[i].name;
	polarizations[i].plus_proj=polarization_results[i].plus_proj;
	polarizations[i].cross_proj=polarization_results[i].cross_proj;
	polarizations[i].plus_factor=polarization_results[i].plus_factor;
	polarizations[i].cross_factor=polarization_results[i].cross_factor;
	polarizations[i].conjugate=&(polarizations[polarization_results[i].conjugate]);
	
	polarizations[i].AM_coeffs=do_alloc(AM_coeffs_size, sizeof(*(polarizations[i].AM_coeffs)));
	for(k=0;k<AM_coeffs_size;k++) {
		polarizations[i].AM_coeffs[k]=AM_coeffs_plus[k]*polarizations[i].plus_proj+
			AM_coeffs_cross[k]*polarizations[i].cross_proj;
		}
	polarizations[i].patch_CutOff=do_alloc(patch_grid->npoints,sizeof(*polarizations[i].patch_CutOff));
	}
}

void allocate_polarization_arrays(void)
{
int i;
	
fprintf(stderr, "Skymap arrays size: %.1f MB\n", ntotal_polarizations*(11.0+2.0*args_info.compute_betas_arg)*fine_grid->npoints*sizeof(SUM_TYPE)/(1024.0*1024.0));
fprintf(LOG, "Skymap arrays size: %f MB\n", ntotal_polarizations*(11.0+2.0*args_info.compute_betas_arg)*fine_grid->npoints*sizeof(SUM_TYPE)/(1024.0*1024.0));

fprintf(stderr, "Spectral plot arrays size: %.1f KB\n", ntotal_polarizations*7.0*useful_bins*fine_grid->nbands*sizeof(SUM_TYPE)/1024.0);
fprintf(LOG, "Spectral plot arrays size: %f KB\n", ntotal_polarizations*7.0*useful_bins*fine_grid->nbands*sizeof(SUM_TYPE)/1024.0);
		
for(i=0;i<ntotal_polarizations;i++) {

	#define SECONDARY  args_info.no_secondary_skymaps_arg ? NULL :

	#ifdef WEIGHTED_SUM
	polarization_results[i].skymap.total_weight=do_alloc(fine_grid->npoints,sizeof(*polarization_results[i].skymap.total_weight));
	polarization_results[i].skymap.total_count= SECONDARY do_alloc(fine_grid->npoints,sizeof(*polarization_results[i].skymap.total_count));
	#else
	polarization_results[i].skymap.total_count=do_alloc(fine_grid->npoints,sizeof(*polarization_results[i].skymap.total_count));
	#endif

	/* Output arrrays - skymaps*/
	polarization_results[i].skymap.max_sub_weight=do_alloc(fine_grid->npoints, sizeof(SUM_TYPE));
	polarization_results[i].skymap.max_dx=do_alloc(fine_grid->npoints, sizeof(SUM_TYPE));
	polarization_results[i].skymap.M_map= SECONDARY do_alloc(fine_grid->npoints, sizeof(SUM_TYPE));
	polarization_results[i].skymap.S_map= SECONDARY do_alloc(fine_grid->npoints, sizeof(SUM_TYPE));
	polarization_results[i].skymap.max_upper_limit=do_alloc(fine_grid->npoints, sizeof(SUM_TYPE));
	polarization_results[i].skymap.max_lower_limit=do_alloc(fine_grid->npoints, sizeof(SUM_TYPE));
	polarization_results[i].skymap.freq_map=do_alloc(fine_grid->npoints, sizeof(SUM_TYPE));
	polarization_results[i].skymap.cor1= SECONDARY do_alloc(fine_grid->npoints, sizeof(SUM_TYPE));
	polarization_results[i].skymap.cor2= SECONDARY do_alloc(fine_grid->npoints, sizeof(SUM_TYPE));
	polarization_results[i].skymap.ks_test=do_alloc(fine_grid->npoints, sizeof(SUM_TYPE));
	polarization_results[i].skymap.ks_count= SECONDARY do_alloc(fine_grid->npoints, sizeof(SUM_TYPE));

	if(args_info.compute_betas_arg){
		polarization_results[i].skymap.beta1=do_alloc(fine_grid->npoints, sizeof(SUM_TYPE));
		polarization_results[i].skymap.beta2=do_alloc(fine_grid->npoints, sizeof(SUM_TYPE));
		} else {
		polarization_results[i].skymap.beta1=NULL;
		polarization_results[i].skymap.beta2=NULL;
		}

	/* Output arrays - spectral plots */
	polarization_results[i].spectral_plot.max_upper_limit=do_alloc(useful_bins*fine_grid->nbands, sizeof(SUM_TYPE));
	polarization_results[i].spectral_plot.ul_dec=do_alloc(useful_bins*fine_grid->nbands, sizeof(SUM_TYPE));
	polarization_results[i].spectral_plot.ul_ra=do_alloc(useful_bins*fine_grid->nbands, sizeof(SUM_TYPE));
	polarization_results[i].spectral_plot.max_dx=do_alloc(useful_bins*fine_grid->nbands, sizeof(SUM_TYPE));
	polarization_results[i].spectral_plot.dx_dec=do_alloc(useful_bins*fine_grid->nbands, sizeof(SUM_TYPE));
	polarization_results[i].spectral_plot.dx_ra=do_alloc(useful_bins*fine_grid->nbands, sizeof(SUM_TYPE));
	polarization_results[i].spectral_plot.max_mask_ratio=do_alloc(useful_bins*fine_grid->nbands, sizeof(SUM_TYPE));

	}

max_cruncher_threads=get_max_threads();

aar=do_alloc(max_cruncher_threads, sizeof(*aar));
for(i=0;i<max_cruncher_threads;i++)aar[i]=new_accumulation_arrays();
}

void free_polarization_arrays(void)
{
int i;
#define FREE(var) { if(var!=NULL)free(var); var=NULL; }
for(i=0;i<ntotal_polarizations;i++){

	#ifdef WEIGHTED_SUM
	FREE(polarization_results[i].skymap.total_weight);
	FREE(polarization_results[i].skymap.total_count);
	#else
	FREE(polarization_results[i].skymap.total_count);
	#endif

	/* Output arrrays - skymaps*/
	FREE(polarization_results[i].skymap.max_sub_weight);
	FREE(polarization_results[i].skymap.max_dx);
	FREE(polarization_results[i].skymap.M_map);
	FREE(polarization_results[i].skymap.S_map);
	FREE(polarization_results[i].skymap.max_upper_limit);
	FREE(polarization_results[i].skymap.max_lower_limit);
	FREE(polarization_results[i].skymap.freq_map);
	FREE(polarization_results[i].skymap.cor1);
	FREE(polarization_results[i].skymap.cor2);
	FREE(polarization_results[i].skymap.ks_test);
	FREE(polarization_results[i].skymap.ks_count);

	if(args_info.compute_betas_arg){
		FREE(polarization_results[i].skymap.beta1);
		FREE(polarization_results[i].skymap.beta2);
		} else {
		polarization_results[i].skymap.beta1=NULL;
		polarization_results[i].skymap.beta2=NULL;
		}

	/* Output arrays - spectral plots */
	FREE(polarization_results[i].spectral_plot.max_upper_limit);
	FREE(polarization_results[i].spectral_plot.ul_dec);
	FREE(polarization_results[i].spectral_plot.ul_ra);
	FREE(polarization_results[i].spectral_plot.max_dx);
	FREE(polarization_results[i].spectral_plot.dx_dec);
	FREE(polarization_results[i].spectral_plot.dx_ra);
	FREE(polarization_results[i].spectral_plot.max_mask_ratio);
	}

for(i=0;i<max_cruncher_threads;i++) {
	free_accumulation_arrays(aar[i]);
	aar[i]=NULL;
	}
free(aar);
aar=NULL;

max_cruncher_threads=-2;
}

void clear_polarization_arrays(void)
{
int i,k;
for(i=0;i<ntotal_polarizations;i++){
	for(k=0;k<fine_grid->npoints;k++){
		polarization_results[i].skymap.max_dx[k]=-1.0;
		if(!args_info.no_secondary_skymaps_arg)polarization_results[i].skymap.M_map[k]=-1.0;
		if(!args_info.no_secondary_skymaps_arg)polarization_results[i].skymap.S_map[k]=-1.0;
		polarization_results[i].skymap.max_upper_limit[k]=-1.0;
		polarization_results[i].skymap.max_lower_limit[k]=-1.0;
		polarization_results[i].skymap.freq_map[k]=-1.0;
		if(!args_info.no_secondary_skymaps_arg)polarization_results[i].skymap.cor1[k]=-1.0;
		if(!args_info.no_secondary_skymaps_arg)polarization_results[i].skymap.cor2[k]=-1.0;
		polarization_results[i].skymap.ks_test[k]=-1.0;
		if(!args_info.no_secondary_skymaps_arg)polarization_results[i].skymap.ks_count[k]=-1.0;
		#ifdef WEIGHTED_SUM
		polarization_results[i].skymap.total_weight[k]=0.0;
		if(!args_info.no_secondary_skymaps_arg)polarization_results[i].skymap.total_count[k]=0;
		#else
		polarization_results[i].skymap.total_count[k]=0;
		#endif

		if(args_info.compute_betas_arg){
			polarization_results[i].skymap.beta1[k]=0.0;
			polarization_results[i].skymap.beta2[k]=0.0;
			}
		}

	for(k=0;k<useful_bins*fine_grid->nbands;k++){
		polarization_results[i].spectral_plot.max_upper_limit[k]=-1.0;
		polarization_results[i].spectral_plot.ul_dec[k]=-10.0;
		polarization_results[i].spectral_plot.ul_ra[k]=-10.0;
		polarization_results[i].spectral_plot.max_dx[k]=-1.0;
		polarization_results[i].spectral_plot.dx_dec[k]=-10.0;
		polarization_results[i].spectral_plot.dx_ra[k]=-10.0;
		polarization_results[i].spectral_plot.max_mask_ratio[k]=-1.0;
		}
	}
}

ACCUMULATION_ARRAYS * new_accumulation_arrays(void)
{
ACCUMULATION_ARRAYS *r;
long total,i;

r=do_alloc(ntotal_polarizations, sizeof(*r));


total=ntotal_polarizations*sizeof(*r[0].fine_grid_sum);

#ifdef COMPUTE_SIGMA
total+=ntotal_polarizations*sizeof(*r[0].fine_grid_sq_sum);
#endif

#ifdef WEIGHTED_SUM
total+=ntotal_polarizations*sizeof(*r[0].fine_grid_sum);
total+=ntotal_polarizations*sizeof(*r[0].fine_grid_count);
#else
total+=ntotal_polarizations*sizeof(*r[0].fine_grid_count);
#endif

if(args_info.compute_betas_arg){
	fprintf(LOG,"Compute betas: true\n");
	} else {
	fprintf(LOG,"Compute betas: false\n");
	}


fprintf(stderr, "Allocating accumulation arrays: %.1f KB\n", 
	(stored_fine_bins*useful_bins*total)/(1024.0));
fprintf(LOG, "Accumulation set size: %f KB\n", 
	(stored_fine_bins*useful_bins*total)/(1024.0));


for(i=0;i<ntotal_polarizations;i++){
	/* Accumulation arrays */
	r[i].total_weight=do_alloc(stored_fine_bins,sizeof(*r[i].total_weight));
	r[i].total_weight_d=do_alloc(stored_fine_bins,sizeof(*r[i].total_weight_d));

	r[i].fine_grid_sum=do_alloc(stored_fine_bins*useful_bins,sizeof(*r[i].fine_grid_sum));
	r[i].fine_grid_sum_d=do_alloc(stored_fine_bins*useful_bins,sizeof(*r[i].fine_grid_sum_d));
	#ifdef COMPUTE_SIGMA
	r[i].fine_grid_sq_sum=do_alloc(stored_fine_bins*useful_bins,sizeof(*r[i].fine_grid_sq_sum));
	#endif


	#ifdef WEIGHTED_SUM
	r[i].fine_grid_weight=do_alloc(stored_fine_bins*useful_bins,sizeof(*r[i].fine_grid_weight));
	r[i].fine_grid_count=do_alloc(stored_fine_bins*useful_bins,sizeof(*r[i].fine_grid_count));
	#else
	r[i].fine_grid_count=do_alloc(stored_fine_bins*useful_bins,sizeof(*r[i].fine_grid_count));
	#endif
	}
return(r);
}

void free_accumulation_arrays(ACCUMULATION_ARRAYS *r)
{
int i;
#define FREE(var) { if(var!=NULL)free(var); var=NULL; }
for(i=0;i<ntotal_polarizations;i++){
	/* Accumulation arrays */
	FREE(r[i].total_weight);
	FREE(r[i].total_weight_d);
	FREE(r[i].fine_grid_sum);
	FREE(r[i].fine_grid_sum_d);
	#ifdef COMPUTE_SIGMA
	FREE(r[i].fine_grid_sq_sum);
	#endif


	#ifdef WEIGHTED_SUM
	FREE(r[i].fine_grid_weight);
	FREE(r[i].fine_grid_count);
	#else
	FREE(r[i].fine_grid_count);
	#endif
	}
free(r);
}

ACCUMULATION_ARRAYS *get_thread_accumulation_arrays(int thread_id)
{
/* Reserve -1 thread_id for main program, get_max_threads always returns one more than have been spawned */
thread_id++;
if(thread_id<0) {
	fprintf(stderr, "Aieee ! - thread id is smaller than -1\n");
	exit(-1);
 	}
if(thread_id>=max_cruncher_threads) {
	fprintf(stderr, "Aieee ! - thread id is too large (%d), maximum is %d\n", thread_id, max_cruncher_threads);
	exit(-1);
 	}
return(aar[thread_id]);
}

void clear_accumulation_arrays(ACCUMULATION_ARRAYS *r)
{
long i,k;

for(k=0;k<ntotal_polarizations;k++){
	for(i=0;i<stored_fine_bins;i++){
		r[k].total_weight[i]=0.0;
		r[k].total_weight_d[i]=0.0;
		}

	for(i=0;i<stored_fine_bins*useful_bins;i++){
		r[k].fine_grid_sum[i]=0.0;
		r[k].fine_grid_sum_d[i]=0.0;
		
		#ifdef COMPUTE_SIGMA
		r[k].fine_grid_sq_sum[i]=0.0;
		#endif

		#ifdef WEIGHTED_SUM	
		r[k].fine_grid_weight[i]=0.0;
		r[k].fine_grid_count[i]=0;
		#else	
		r[k].fine_grid_count[i]=0;
		#endif
		}
	}
}

void update_d_accumulation_arrays(ACCUMULATION_ARRAYS *r)
{
long i,k;

for(k=0;k<ntotal_polarizations;k++){
	for(i=0;i<stored_fine_bins;i++){
		r[k].total_weight_d[i]+=r[k].total_weight[i];
		r[k].total_weight[i]=0.0;
		}

	for(i=0;i<stored_fine_bins*useful_bins;i++){
		r[k].fine_grid_sum_d[i]+=r[k].fine_grid_sum[i];
		r[k].fine_grid_sum[i]=0.0;
		}
	}
}

void finalize_accumulation_arrays(ACCUMULATION_ARRAYS *r)
{
long i,k;

for(k=0;k<ntotal_polarizations;k++){
	for(i=0;i<stored_fine_bins;i++){
		r[k].total_weight_d[i]+=r[k].total_weight[i];
		r[k].total_weight[i]=r[k].total_weight_d[i];
		}

	for(i=0;i<stored_fine_bins*useful_bins;i++){
		r[k].fine_grid_sum_d[i]+=r[k].fine_grid_sum[i];
		r[k].fine_grid_sum[i]=r[k].fine_grid_sum_d[i];
		}
	}
}
