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
/* We need this define to get NAN values */
#define __USE_ISOC99
#include <math.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/time.h>
#include <time.h>
#include <string.h>

#include "global.h"
#include "rastermagic.h"
#include "cmdline.h"
#include "hookup.h"
#include "grid.h"
#include "polarization.h"
#include "statistics.h"
#include "dataset.h"
#include "candidates.h"
#include "util.h"
#include "jobs.h"

extern POLARIZATION_RESULTS *polarization_results;
extern int nlinear_polarizations, ntotal_polarizations;

extern DATASET *datasets;
extern int d_free;

extern struct gengetopt_args_info args_info;

extern int nbins, first_bin, side_cut, useful_bins;

INT64 spindown_start;

extern SKY_GRID *fine_grid, *patch_grid;
extern SKY_SUPERGRID *super_grid;

extern int do_CutOff;
extern int fake_injection;

extern FILE *LOG;

extern double spindown;
extern double resolution;

extern int subinstance;
extern char *subinstance_name;

SUM_TYPE normalizing_weight;

int stored_fine_bins=0;

int max_shift=0, min_shift=0;


SUM_TYPE  *accum_circ_ul, *accum_circ_ul_freq;  /* lower (over polarizations) accumulation limits */
SUM_TYPE *skymap_circ_ul, *skymap_circ_ul_freq; /* skymaps */
SUM_TYPE *spectral_plot_circ_ul; /* spectral plots */

SUM_TYPE *max_dx=NULL;
short *max_dx_polarization_index=NULL;

float upper_limit_comp;
float lower_limit_comp;

/* Include file with Feldman-Cousins upper limits directly
   so as to benefit from function inlining */
   
#include "fc.c"

/* Power cache - given dataset, SFT number and frequency shift 
	obtain a pointer to array of powers ready for summing */

struct {
	DATASET *d;
	int segment;
	float **power;
	float *shift;
	int free;
	int size;
	int hits;
	long long total_hits;
	long long total_misses;
	} * power_cache;

void init_power_cache(void) 
{
int i;

power_cache=do_alloc(get_max_threads(), sizeof(*power_cache));

for(i=0;i<get_max_threads();i++) {
	power_cache[i].d=NULL;
	power_cache[i].d=NULL;
	power_cache[i].segment=-1;
	power_cache[i].power=NULL;
	power_cache[i].shift=NULL;
	power_cache[i].free=0;
	power_cache[i].size=0;
	power_cache[i].hits=0;
	power_cache[i].total_hits=0;
	power_cache[i].total_misses=0;
	}
}

/* Use matched filter to estimate power */
float *get_matched_power(int thread_id, float shift, DATASET *d, int k)
{
float *bin_re, *bin_im;
int i, m, bin_shift;
float filter[7];
float *power;
float x, y;

thread_id++;
if(thread_id<0) {
	fprintf(stderr, "Aieee ! - thread id is smaller than -1\n");
	exit(-1);
 	}

if( (power_cache[thread_id].d!=d) || power_cache[thread_id].segment!=k) {
	power_cache[thread_id].d=d;
	power_cache[thread_id].segment=k;
	power_cache[thread_id].free=0;
	if(power_cache[thread_id].size==0) {
		power_cache[thread_id].size=100;
		power_cache[thread_id].power=do_alloc(power_cache[thread_id].size, sizeof(*power_cache[thread_id].power));
		power_cache[thread_id].shift=do_alloc(power_cache[thread_id].size, sizeof(*power_cache[thread_id].shift));
		for(i=0;i<power_cache[thread_id].size;i++) {
			power_cache[thread_id].power[i]=do_alloc(nbins-2*side_cut, sizeof(**power_cache[thread_id].power));
			}
		}
	}

for(m=0;m<power_cache[thread_id].free;m++) {
	/* There is no difference in SNR between 0.01 and 0.05, while
           there are substantial savings in run time (at least 14%) */
	if(fabs(power_cache[thread_id].shift[m]-shift)<0.05) {
		power_cache[thread_id].hits++;
		return(power_cache[thread_id].power[m]);
		}
	}
power_cache[thread_id].total_misses++;
power_cache[thread_id].total_hits+=power_cache[thread_id].hits+1;
power_cache[thread_id].hits=0;
m=power_cache[thread_id].free;
power_cache[thread_id].free++;

if(m>=power_cache[thread_id].size) {
	fprintf(stderr, "Aieee ! power cache overflowed\n");
	exit(-1);
	}

power_cache[thread_id].shift[m]=shift;

bin_shift=rintf(shift);

bin_re=&(d->re[k*nbins+side_cut+bin_shift]);
bin_im=&(d->im[k*nbins+side_cut+bin_shift]);

tabulated_fill_hann_filter7(filter, (shift-bin_shift));

power=power_cache[thread_id].power[m];

for(i=nbins-2*side_cut;i>0;i--) {
	
	x=bin_re[-3]*filter[0]+bin_re[-2]*filter[1]+bin_re[-1]*filter[2]+bin_re[0]*filter[3]+bin_re[1]*filter[4]+bin_re[2]*filter[5]+bin_re[3]*filter[6];
	y=bin_im[-3]*filter[0]+bin_im[-2]*filter[1]+bin_im[-1]*filter[2]+bin_im[0]*filter[3]+bin_im[1]*filter[4]+bin_im[2]*filter[5]+bin_im[3]*filter[6];

	*power=x*x+y*y;
	bin_re++;
	bin_im++;
	power++;
	}

return(power_cache[thread_id].power[m]);
}

/* Just return the closest bin - this way we don't have to store precomputed power */
float *get_closest_bin_power(int thread_id, float shift, DATASET *d, int k)
{
float *bin_re, *bin_im;
int i, m, bin_shift;
float *power, *fm;
float x, y;

thread_id++;
if(thread_id<0) {
	fprintf(stderr, "Aieee ! - thread id is smaller than -1\n");
	exit(-1);
 	}

if( (power_cache[thread_id].d!=d) || power_cache[thread_id].segment!=k) {
	power_cache[thread_id].d=d;
	power_cache[thread_id].segment=k;
	power_cache[thread_id].free=0;
	if(power_cache[thread_id].size==0) {
		power_cache[thread_id].size=100;
		power_cache[thread_id].power=do_alloc(power_cache[thread_id].size, sizeof(*power_cache[thread_id].power));
		power_cache[thread_id].shift=do_alloc(power_cache[thread_id].size, sizeof(*power_cache[thread_id].shift));
		for(i=0;i<power_cache[thread_id].size;i++) {
			power_cache[thread_id].power[i]=do_alloc(nbins-2*side_cut+2, sizeof(**power_cache[thread_id].power));
			}
		}
	}

for(m=0;m<power_cache[thread_id].free;m++) {
	/* we want to find closest bin */
	if(fabs(power_cache[thread_id].shift[m]-shift)<0.5) {
		power_cache[thread_id].hits++;
		return(power_cache[thread_id].power[m]+1);
		}
	}
power_cache[thread_id].total_misses++;
power_cache[thread_id].total_hits+=power_cache[thread_id].hits+1;
power_cache[thread_id].hits=0;
m=power_cache[thread_id].free;
power_cache[thread_id].free++;

if(m>=power_cache[thread_id].size) {
	fprintf(stderr, "Aieee ! power cache overflowed\n");
	exit(-1);
	}

bin_shift=rintf(shift);

power_cache[thread_id].shift[m]=bin_shift;

bin_re=&(d->re[k*nbins+side_cut+bin_shift-1]);
bin_im=&(d->im[k*nbins+side_cut+bin_shift-1]);

power=power_cache[thread_id].power[m];

/* We could have computed all nbins bins here - however usually bin_shifts are very close together
and 2*side_cut is comparable (or larger) to nbins */

for(i=nbins-2*side_cut+2;i>0;i--) {
	
	x=bin_re[0];
	y=bin_im[0];

	*power=x*x+y*y;
	bin_re++;
	bin_im++;
	power++;
	}

if(args_info.subtract_background_arg) {
	power=power_cache[thread_id].power[m];
	fm=&(d->expFMedians_plain[side_cut+bin_shift-1]);
	for(i=nbins-2*side_cut+2;i>0;i--) {
		(*power)-=d->expTMedians_plain[k]*(*fm);
		power++;
		fm++;
		}
	}

return(power_cache[thread_id].power[m]+1);
}


void (*process_patch)(int thread_id, DATASET *d, ACCUMULATION_ARRAYS *ar, int pol_index, int pi, int k, float CutOff);



static void process_patch1(int thread_id, DATASET *d, ACCUMULATION_ARRAYS *ar, int pol_index, int pi, int k, float CutOff)
{
int i,kk,b,b0,b1,n,offset;
int bin_shift;
SUM_TYPE mod;
SUM_TYPE a,w,w2;
SUM_TYPE *sum,*sq_sum;
float *p, *power;
float doppler, bin_offset, shift;
float f_plus, f_cross, beta1, beta2;
ACCUMULATION_ARRAYS *aa=&(ar[pol_index]);
POLARIZATION_RESULTS *pr=&(polarization_results[pol_index]);
POLARIZATION *pl=&(d->polarizations[pol_index]);


CutOff=2*CutOff; /* weighted sum can benefit from more SFTs */

bin_offset=d->coherence_time*(args_info.frequency_offset_arg+spindown*(d->gps[k]-spindown_start)+0.5*args_info.fdotdot_arg*(d->gps[k]-spindown_start)*(d->gps[k]-spindown_start));

for(i=0,kk=super_grid->first_map[pi];kk>=0;kk=super_grid->list_map[kk],i++)
		{

		if(fine_grid->band[kk]<0)continue;
		
		/* Get amplitude response */
		f_plus=F_plus(k, fine_grid, kk, pl->AM_coeffs);
		f_cross=F_plus(k, fine_grid, kk, pl->conjugate->AM_coeffs);


		mod=1.0/(pl->plus_factor*f_plus*f_plus+pl->cross_factor*f_cross*f_cross); 

		if(do_CutOff && (mod>CutOff))continue;

		/* this assumes that both relfreq and spindown are small enough that the bin number
		   down not change much within the band - usually true */
		doppler=fine_grid->e[0][kk]*d->detector_velocity[3*k+0]+
			fine_grid->e[1][kk]*d->detector_velocity[3*k+1]+
			fine_grid->e[2][kk]*d->detector_velocity[3*k+2];
		doppler=args_info.doppler_multiplier_arg*doppler;

		shift=((first_bin+nbins*0.5)*doppler+bin_offset);
		bin_shift=-rint(shift);

		if(bin_shift>max_shift)max_shift=bin_shift;
		if(bin_shift<min_shift)min_shift=bin_shift;
		b0=side_cut-bin_shift;
		b1=(nbins-side_cut)-bin_shift;
		if((b0<0)||(b1>nbins)){
			fprintf(stderr,"Working frequency range obscured by bin_shift shift: bin_shift=%d kk=%d i=%d pi=%d\n",
				bin_shift, kk, i, pi);
			exit(-1);
			}
		
		if(args_info.compute_betas_arg){
		        beta1=f_cross*f_plus*mod;	
			beta2=(-pl->cross_factor*f_plus*f_plus+pl->plus_factor*f_cross*f_cross)*mod;
			}

		/* prime array pointers */
		#ifdef WEIGHTED_SUM
		w2=d->expTMedians[k]*d->weight/mod;
		w=w2/mod;
		//pr->skymap.total_weight[kk]+=w;
		aa->total_weight[i]+=w;
	
		if(args_info.compute_betas_arg){
			pr->skymap.beta1[kk]+=w*beta1;
			pr->skymap.beta2[kk]+=w*beta2;
			}
		if(!args_info.no_secondary_skymaps_arg)pr->skymap.total_count[kk]++;
		#else
		pr->skymap.total_count[kk]++;

		if(args_info.compute_betas_arg){
			pr->skymap.beta1[kk]+=beta1;
			pr->skymap.beta2[kk]+=beta2;
			}
		#endif
		
		sum=&(aa->fine_grid_sum[b0-side_cut+bin_shift+useful_bins*i]);
		#ifdef COMPUTE_SIGMA
		sq_sum=&(aa->fine_grid_sq_sum[b0-side_cut+bin_shift+useful_bins*i]);
		#endif

		power=get_closest_bin_power(thread_id, shift, d, k);
		/* p=&(d->power[k*nbins+b0]); */
		
		p=power;
		/* cycle over bins */
		for(b=b0;b<b1;b++){			    
			#ifdef WEIGHTED_SUM
			a=(*p)*w2;
			(*sum)+=a;
			#else
			a=(*p)*mod;
			(*sum)+=a;
			#ifdef COMPUTE_SIGMA
			(*sq_sum)+=a*a;
			sq_sum++;
			#endif			
			#endif
		
			p++;
			sum++;
			}
		/* subtract lines */
		for(n=0;(d->lines_report->lines_list[n]>=0)&&(n<d->lines_report->nlines);n++) {
			b=d->lines_report->lines_list[n];
			if(b<b0)continue;
			if(b>=b1)continue;
			offset=b-side_cut+bin_shift+useful_bins*i;
			/* prime array pointers */
			sum=&(aa->fine_grid_sum[offset]);
			/* p=&(d->power[k*nbins+b]); */
			p=&(power[k*nbins+b-b0]);

			#ifdef WEIGHTED_SUM
			aa->fine_grid_weight[offset]+=w;
			#else
			aa->fine_grid_count[offset]++;
			#endif
			
			#ifdef WEIGHTED_SUM
			a=(*p)*w2;
			(*sum)-=a;
			#else
			a=(*p)*mod;
			(*sum)-=a;
			#ifdef COMPUTE_SIGMA
			pr->fine_grid_sq_sum[offset]-=a*a;
			#endif			
			#endif
			}
		}

}

/* three bin version */

static void process_patch3(int thread_id, DATASET *d, ACCUMULATION_ARRAYS *ar, int pol_index, int pi, int k, float CutOff)
{
int i,kk,b,b0,b1,n,offset;
int bin_shift;
SUM_TYPE mod;
SUM_TYPE a,w,w2;
SUM_TYPE *sum,*sq_sum;
float *p, *power;
float doppler, bin_offset, shift;
float f_plus, f_cross, beta1, beta2;
ACCUMULATION_ARRAYS *aa=&(ar[pol_index]);
POLARIZATION_RESULTS *pr=&(polarization_results[pol_index]);
POLARIZATION *pl=&(d->polarizations[pol_index]);

CutOff=2*CutOff/3.0; /* weighted sum can benefit from more SFTs */

bin_offset=d->coherence_time*(args_info.frequency_offset_arg+spindown*(d->gps[k]-spindown_start)+0.5*args_info.fdotdot_arg*(d->gps[k]-spindown_start)*(d->gps[k]-spindown_start));

for(i=0,kk=super_grid->first_map[pi];kk>=0;kk=super_grid->list_map[kk],i++)
		{
		if(fine_grid->band[kk]<0)continue;

		/* Get amplitude response */
		f_plus=F_plus(k, fine_grid, kk, pl->AM_coeffs);
		f_cross=F_plus(k, fine_grid, kk, pl->conjugate->AM_coeffs);


		mod=1.0/(pl->plus_factor*f_plus*f_plus+pl->cross_factor*f_cross*f_cross); 

		if(do_CutOff && (mod>CutOff))continue;

		
		/* this assumes that both relfreq and spindown are small enough that the bin number
		   down not change much within the band - usually true */
		doppler=fine_grid->e[0][kk]*d->detector_velocity[3*k+0]+
			fine_grid->e[1][kk]*d->detector_velocity[3*k+1]+
			fine_grid->e[2][kk]*d->detector_velocity[3*k+2];
		doppler=args_info.doppler_multiplier_arg*doppler;

		shift=((first_bin+nbins*0.5)*doppler+bin_offset);
		bin_shift=-rint(shift);

		if(bin_shift>max_shift)max_shift=bin_shift;
		if(bin_shift<min_shift)min_shift=bin_shift;
		b0=side_cut-bin_shift;
		b1=(nbins-side_cut)-bin_shift;
		if((b0<1)||(b1>(nbins-1))){
			fprintf(stderr,"Working frequency range obscured by bin_shift shift: bin_shift=%d kk=%d i=%d pi=%d\n",
				bin_shift, kk, i, pi);
			exit(-1);
			}
		
		if(args_info.compute_betas_arg){
		        beta1=f_cross*f_plus*mod;	
			beta2=(-pl->cross_factor*f_plus*f_plus+pl->plus_factor*f_cross*f_cross)*mod;
			}

		/* prime array pointers */
		#ifdef WEIGHTED_SUM
		w2=d->expTMedians[k]*d->weight/mod;
		w=w2/(3.0*mod);
		//pr->skymap.total_weight[kk]+=w;
		aa->total_weight[i]+=w;

		if(args_info.compute_betas_arg){
			pr->skymap.beta1[kk]+=w*beta1;
			pr->skymap.beta2[kk]+=w*beta2;
			}
		if(!args_info.no_secondary_skymaps_arg)pr->skymap.total_count[kk]++;
		#else
		pr->skymap.total_count[kk]++;

		if(args_info.compute_betas_arg){
			pr->skymap.beta1[kk]+=beta1;
			pr->skymap.beta2[kk]+=beta2;
			}
		#endif
		
		sum=&(aa->fine_grid_sum[b0-side_cut+bin_shift+useful_bins*i]);
		#ifdef COMPUTE_SIGMA
		sq_sum=&(pr->fine_grid_sq_sum[b0-side_cut+bin_shift+useful_bins*i]);
		#endif

		power=get_closest_bin_power(thread_id, shift, d, k);
		/* p=&(d->power[k*nbins+b0]); */

		p=power;		
		/* cycle over bins */
		for(b=b0;b<b1;b++){			    
			#ifdef WEIGHTED_SUM
			a=(p[-1]+p[0]+p[1])*w2;
			(*sum)+=a;
			#else
			a=(*p)*mod;
			(*sum)+=a;
			#ifdef COMPUTE_SIGMA
			(*sq_sum)+=a*a;
			sq_sum++;
			#endif			
			#endif
		
			p++;
			sum++;
			}
		/* subtract lines */
		for(n=0;(d->lines_report->lines_list[n]>=0)&&(n<d->lines_report->nlines);n++){
			b=d->lines_report->lines_list[n];
			if(b<b0)continue;
			if(b>=b1)continue;
			offset=b-side_cut+bin_shift+useful_bins*i;
			/* prime array pointers */
			sum=&(aa->fine_grid_sum[offset]);
			p=&(power[k*nbins+b-b0]);

			#ifdef WEIGHTED_SUM
			aa->fine_grid_weight[offset]+=w;
			#else
			aa->fine_grid_count[offset]++;
			#endif
			
			#ifdef WEIGHTED_SUM
			a=(p[-1]+p[0]+p[1])*w2;
			(*sum)-=a;
			#else
			a=(p[-1]+p[0]+p[1])*mod;
			(*sum)-=a;
			#ifdef COMPUTE_SIGMA
			pr->fine_grid_sq_sum[offset]-=a*a;
			#endif			
			#endif
			}
		}

}


static void process_patch_matched(int thread_id, DATASET *d, ACCUMULATION_ARRAYS *ar, int pol_index, int pi, int k, float CutOff)
{
int i,kk,b,b0,b1,n,offset;
int bin_shift;
SUM_TYPE mod;
SUM_TYPE a,w,w2;
SUM_TYPE *sum,*sq_sum;
float *power;
float doppler, bin_offset;
float f_plus, f_cross, beta1, beta2;
ACCUMULATION_ARRAYS *aa=&(ar[pol_index]);
POLARIZATION_RESULTS *pr=&(polarization_results[pol_index]);
POLARIZATION *pl=&(d->polarizations[pol_index]);
float shift;


CutOff=2*CutOff; /* weighted sum can benefit from more SFTs */

bin_offset=d->coherence_time*(args_info.frequency_offset_arg+spindown*(d->gps[k]-spindown_start)+0.5*args_info.fdotdot_arg*(d->gps[k]-spindown_start)*(d->gps[k]-spindown_start));

for(i=0,kk=super_grid->first_map[pi];kk>=0;kk=super_grid->list_map[kk],i++)
		{

		if(fine_grid->band[kk]<0)continue;
		
		/* Get amplitude response */
		f_plus=F_plus(k, fine_grid, kk, pl->AM_coeffs);
		f_cross=F_plus(k, fine_grid, kk, pl->conjugate->AM_coeffs);


		mod=1.0/(pl->plus_factor*f_plus*f_plus+pl->cross_factor*f_cross*f_cross); 

		if(do_CutOff && (mod>CutOff))continue;

		/* this assumes that both relfreq and spindown are small enough that the bin number
		   down not change much within the band - usually true */
		doppler=fine_grid->e[0][kk]*d->detector_velocity[3*k+0]+
			fine_grid->e[1][kk]*d->detector_velocity[3*k+1]+
			fine_grid->e[2][kk]*d->detector_velocity[3*k+2];
		doppler=args_info.doppler_multiplier_arg*doppler;

		shift=((first_bin+nbins*0.5)*doppler+bin_offset);
		bin_shift=rintf(shift);

		if(bin_shift>max_shift)max_shift=bin_shift;
		if(bin_shift<min_shift)min_shift=bin_shift;
		b0=side_cut+bin_shift;
		b1=(nbins-side_cut)+bin_shift;
		if((b0<3)||(b1>nbins-3)){
			fprintf(stderr,"Working frequency range obscured by bin_shift shift: bin_shift=%d kk=%d i=%d pi=%d\n",
				bin_shift, kk, i, pi);
			exit(-1);
			}
		
		if(args_info.compute_betas_arg){
		        beta1=f_cross*f_plus*mod;	
			beta2=(-pl->cross_factor*f_plus*f_plus+pl->plus_factor*f_cross*f_cross)*mod;
			}

		/* prime array pointers */
		#ifdef WEIGHTED_SUM
		w2=d->expTMedians[k]*d->weight/mod;
		w=w2/mod;
		//pr->skymap.total_weight[kk]+=w;
		aa->total_weight[i]+=w;
	
		if(args_info.compute_betas_arg){
			pr->skymap.beta1[kk]+=w*beta1;
			pr->skymap.beta2[kk]+=w*beta2;
			}
		if(!args_info.no_secondary_skymaps_arg)pr->skymap.total_count[kk]++;
		#else
		pr->skymap.total_count[kk]++;

		if(args_info.compute_betas_arg){
			pr->skymap.beta1[kk]+=beta1;
			pr->skymap.beta2[kk]+=beta2;
			}
		#endif
		
		sum=&(aa->fine_grid_sum[useful_bins*i]);
		#ifdef COMPUTE_SIGMA
		sq_sum=&(aa->fine_grid_sq_sum[useful_bins*i]);
		#endif

		power=get_matched_power(thread_id, shift, d, k);

		/* cycle over bins */
		for(b=b0;b<b1;b++) {
			
			#ifdef WEIGHTED_SUM
			a=(*power)*w2;
			(*sum)+=a;
			#else
			a=(*power)*mod;
			(*sum)+=a;
			#ifdef COMPUTE_SIGMA
			(*sq_sum)+=a*a;
			sq_sum++;
			#endif
			#endif

			power++;	
			sum++;
			}

		#if 0  /* filter is 7 bins wide - do not do line subtraction for now */
		/* subtract lines */
		for(n=0;(d->lines_report->lines_list[n]>=0)&&(n<d->lines_report->nlines);n++) {
			b=d->lines_report->lines_list[n];
			if(b<b0)continue;
			if(b>=b1)continue;
			offset=b-side_cut+bin_shift+useful_bins*i;
			/* prime array pointers */
			sum=&(pr->fine_grid_sum[offset]);
			p=&(d->power[k*nbins+b]);

			#ifdef WEIGHTED_SUM
			pr->fine_grid_weight[offset]+=w;
			#else
			pr->fine_grid_count[offset]++;
			#endif
			
			#ifdef WEIGHTED_SUM
			a=(*p)*w2;
			(*sum)-=a;
			#else
			a=(*p)*mod;
			(*sum)-=a;
			#ifdef COMPUTE_SIGMA
			pr->fine_grid_sq_sum[offset]-=a*a;
			#endif			
			#endif
			}
		#endif
		}
}


void dump_pic(char *file, double *z)
{
RGBPic *p;

if(!clear_name_png(file))return;

p=make_RGBPic(fine_grid->max_n_ra+140, fine_grid->max_n_dec);

plot_grid_d(p, fine_grid, z, 1);
RGBPic_dump_png(file, p);

free_RGBPic(p);
}

int float_cmp(float *a, float *b);

void make_limits(int thread_id, POLARIZATION_RESULTS *pol, ACCUMULATION_ARRAYS *ar, int pi)
{
SUM_TYPE M, S,dx;
SUM_TYPE a,b,c;
SUM_TYPE *tmp=NULL;
int i,j,k,offset,band;
NORMAL_STATS nstats;
SUM_TYPE *circ_ul, *circ_ul_freq;

circ_ul=&(accum_circ_ul[stored_fine_bins*useful_bins*(thread_id+1)]);
circ_ul_freq=&(accum_circ_ul[stored_fine_bins*useful_bins*(thread_id+1)]);

/* allocate on stack, for speed */
tmp=alloca(useful_bins*sizeof(*tmp));

/* sort to compute robust estimates */
nstats.flag=STAT_FLAG_INPLACE_SORT_DATA
	| STAT_FLAG_ESTIMATE_MEAN
	| STAT_FLAG_ESTIMATE_SIGMA;

if(args_info.ks_test_arg){
	nstats.flag|=STAT_FLAG_ESTIMATE_KS_LEVEL
		| STAT_FLAG_COMPUTE_KS_TEST;
	}

#define SECONDARY if(!args_info.no_secondary_skymaps_arg)

for(i=0,offset=super_grid->first_map[pi];offset>=0;offset=super_grid->list_map[offset],i++){

        band=fine_grid->band[offset];
	if(band<0)continue;

	/* figure out offset to put results into */

	memcpy(tmp, ar->fine_grid_sum+i*useful_bins,useful_bins*sizeof(*tmp));
	/* compute correlations - for diagnostics only */
	a=0.0;
	b=0.0;
	for(j=0;j<useful_bins-2;j++){
		SECONDARY pol->skymap.cor1[offset]+=(tmp[j]*normalizing_weight)*(tmp[j+1]*normalizing_weight);
		SECONDARY pol->skymap.cor2[offset]+=(tmp[j]*normalizing_weight)*(tmp[j+2]*normalizing_weight);
		b+=(tmp[j]*normalizing_weight);
		a+=(tmp[j]*normalizing_weight)*(tmp[j]*normalizing_weight);
		}
	c=(b-tmp[0]*normalizing_weight-tmp[1]*normalizing_weight+tmp[useful_bins-2]*normalizing_weight+tmp[useful_bins-1]*normalizing_weight);
	SECONDARY pol->skymap.cor2[offset]=(pol->skymap.cor2[offset]-b*c/(useful_bins-2))/
		(sqrt((a-b*b/(useful_bins-2))*
		(a-tmp[0]*normalizing_weight*tmp[0]*normalizing_weight-tmp[1]*normalizing_weight*tmp[1]*normalizing_weight
			+tmp[useful_bins-2]*normalizing_weight*tmp[useful_bins-2]*normalizing_weight+tmp[useful_bins-1]*normalizing_weight*tmp[useful_bins-1]*normalizing_weight
			-c*c/(useful_bins-2))));
			
	b+=tmp[useful_bins-2]*normalizing_weight;
	a+=tmp[useful_bins-2]*normalizing_weight*tmp[useful_bins-2]*normalizing_weight;
	SECONDARY pol->skymap.cor1[offset]+=tmp[useful_bins-2]*normalizing_weight*tmp[useful_bins-1]*normalizing_weight;
	c=b-tmp[0]*normalizing_weight+tmp[useful_bins-1]*normalizing_weight;
	SECONDARY pol->skymap.cor1[offset]=(pol->skymap.cor1[offset]-b*c/(useful_bins-1))/
		(sqrt(
		(a-b*b/(useful_bins-1))*
		(a-tmp[0]*normalizing_weight*tmp[0]*normalizing_weight+tmp[useful_bins-1]*normalizing_weight*tmp[useful_bins-1]*normalizing_weight-c*c/(useful_bins-1))
		));
	
	
	/* Output point data if requested */
	if(args_info.dump_points_arg) {
		char s[20000];
		snprintf(s, 20000, "points/%s%s_%d.png", subinstance_name, pol->name, offset);
		if(clear_name_png(s)) {
			RGBPic *p;
			PLOT *plot;
			float *freq_f;

			freq_f=do_alloc(useful_bins, sizeof(*freq_f));
			for(j=0;j<useful_bins;j++) {
				freq_f[j]=(first_bin+side_cut+j)/1800.0;
				}
		
			if(fine_grid->max_n_dec<800)
				p=make_RGBPic(fine_grid->max_n_ra*(800/fine_grid->max_n_dec)+140, fine_grid->max_n_dec*(800/fine_grid->max_n_dec));
				else 
				p=make_RGBPic(fine_grid->max_n_ra+140, fine_grid->max_n_dec);	

			plot=make_plot(p->width, p->height);

			adjust_plot_limits_f(plot, freq_f, tmp, useful_bins, 1, 1, 1);
			draw_grid(p, plot, 0, 0);
			draw_points_f(p, plot, COLOR(255,0,0), freq_f, tmp, useful_bins, 1, 1);
			RGBPic_dump_png(s, p);


			free_plot(plot);
			free_RGBPic(p);
			free(freq_f);
			}
			
		snprintf(s, 20000, "points/%s%s_%d.dat", subinstance_name, pol->name, offset);
		dump_floats(s, tmp, useful_bins, 1);
		}

	compute_normal_stats(tmp, useful_bins, &nstats);

	pol->skymap.ks_test[offset]=nstats.ks_test;
	SECONDARY pol->skymap.ks_count[offset]=nstats.ks_count;
	
	M=nstats.mean;
	S=nstats.sigma;
	
	SECONDARY pol->skymap.M_map[offset]=M;
	SECONDARY pol->skymap.S_map[offset]=S;
	pol->skymap.max_upper_limit[offset]=0;
	
	for(k=0;k<useful_bins;k++){
		dx=(ar->fine_grid_sum[i*useful_bins+k]-M)/S;		
		a=upper_limit95(dx)*S;
		if(a>pol->skymap.max_upper_limit[offset]) {
			pol->skymap.max_upper_limit[offset]=a;
			pol->skymap.freq_map[offset]=(first_bin+side_cut+k)/1800.0;
			}
		if(a>pol->spectral_plot.max_upper_limit[k+band*useful_bins]){
			pol->spectral_plot.max_upper_limit[k+band*useful_bins]=a;
			pol->spectral_plot.ul_ra[k+band*useful_bins]=fine_grid->longitude[offset];
			pol->spectral_plot.ul_dec[k+band*useful_bins]=fine_grid->latitude[offset];
			}

		if(dx>pol->skymap.max_dx[offset]){
			pol->skymap.max_dx[offset]=dx;
			}
		if(dx>pol->spectral_plot.max_dx[k+band*useful_bins]){
			pol->spectral_plot.max_dx[k+band*useful_bins]=dx;
			pol->spectral_plot.dx_ra[k+band*useful_bins]=fine_grid->longitude[offset];
			pol->spectral_plot.dx_dec[k+band*useful_bins]=fine_grid->latitude[offset];
			}
			
		/* circ_ul describes limit on circularly polarized signals */
		/* this formula should only use linear polarizations */
		if(args_info.compute_betas_arg && (pol->cross_factor==0)) {
			a*=1.0/(1.0+pol->skymap.beta2[offset]);
			if(a<circ_ul[i*useful_bins+k]){
				circ_ul[i*useful_bins+k]=a;
				circ_ul_freq[i*useful_bins+k]=(first_bin+side_cut+k)/1800.0;
				}
			}

			
		a=lower_limit95(dx)*S;
		if(a>pol->skymap.max_lower_limit[offset]){
			pol->skymap.max_lower_limit[offset]=a;
			}
			
		#ifdef WEIGHTED_SUM
		a=ar->fine_grid_weight[i*useful_bins+k]/pol->skymap.total_weight[offset];
		#else
		a=ar->fine_grid_count[i*useful_bins+k]/pol->skymap.total_count[offset];
		#endif

		if(a>pol->spectral_plot.max_mask_ratio[k+band*useful_bins]){
			pol->spectral_plot.max_mask_ratio[k+band*useful_bins]=a;
			}
		}
	}

}

void make_unified_limits(int thread_id, int pi)
{
int i, offset;
int band, k;
SUM_TYPE a,b,c;
SUM_TYPE *circ_ul, *circ_ul_freq;

circ_ul=&(accum_circ_ul[stored_fine_bins*useful_bins*(thread_id+1)]);
circ_ul_freq=&(accum_circ_ul[stored_fine_bins*useful_bins*(thread_id+1)]);

for(i=0,offset=super_grid->first_map[pi];offset>=0;offset=super_grid->list_map[offset],i++) {

        band=fine_grid->band[offset];
	if(band<0)continue;

	for(k=0;k<useful_bins;k++){
		a=circ_ul[i*useful_bins+k];
		if(a>spectral_plot_circ_ul[k+band*useful_bins]){
			spectral_plot_circ_ul[k+band*useful_bins]=a;
			}
		if(a>skymap_circ_ul[offset]){
			skymap_circ_ul[offset]=a;
			skymap_circ_ul_freq[offset]=circ_ul_freq[i*useful_bins+k];
			}
		}
	}
}

void output_limits(POLARIZATION_RESULTS *pol)
{
char s[20000];
RGBPic *p;
PLOT *plot;
int i, max_dx_i, largest_i, masked, k;
SUM_TYPE max_dx, largest;
float *max_band, *masked_max_band;
int *max_band_arg, *masked_max_band_arg;
float *freq_f;
float max_ratio;
HISTOGRAM *hist;

freq_f=do_alloc(useful_bins, sizeof(*freq_f));
for(i=0;i<useful_bins;i++)freq_f[i]=(first_bin+side_cut+i)/1800.0;

max_band=do_alloc(fine_grid->nbands, sizeof(*max_band));
masked_max_band=do_alloc(fine_grid->nbands, sizeof(*max_band));
max_band_arg=do_alloc(fine_grid->nbands, sizeof(*max_band_arg));
masked_max_band_arg=do_alloc(fine_grid->nbands, sizeof(*max_band_arg));
hist=new_histogram(args_info.hist_bins_arg, fine_grid->nbands);

if(fine_grid->max_n_dec<800){
	p=make_RGBPic(fine_grid->max_n_ra*(800/fine_grid->max_n_dec)+140, fine_grid->max_n_dec*(800/fine_grid->max_n_dec));
	} else 
	p=make_RGBPic(fine_grid->max_n_ra+140, fine_grid->max_n_dec);	

plot=make_plot(p->width, p->height);

#define OUTPUT_SKYMAP(format, field)	{ \
	if(pol->skymap.field!=NULL) { \
		snprintf(s,19999, "%s" format ".png", subinstance_name, pol->name); \
		if(clear_name_png(s)){ \
			plot_grid_f(p, fine_grid, pol->skymap.field, 1); \
			RGBPic_dump_png(s, p); \
			} \
		snprintf(s,19999, "%s" format ".dat", subinstance_name, pol->name); \
		dump_floats(s, pol->skymap.field, fine_grid->npoints, 1); \
		} \
	}

OUTPUT_SKYMAP("%s_weight", total_weight);
OUTPUT_SKYMAP("%s_count", total_count);

if(args_info.compute_betas_arg){
	OUTPUT_SKYMAP("%s_beta1", beta1);
	OUTPUT_SKYMAP("%s_beta2", beta2);
	}

OUTPUT_SKYMAP("%s_cor1", cor1);
OUTPUT_SKYMAP("%s_cor2", cor2);

if(args_info.ks_test_arg){
	OUTPUT_SKYMAP("%s_ks_test", ks_test);

	compute_histogram_f(hist, pol->skymap.ks_test, fine_grid->band, fine_grid->npoints);
	snprintf(s,19999,"%shist_%s_ks_test", subinstance_name, pol->name);
	print_histogram(LOG, hist, s);
	
	OUTPUT_SKYMAP("%s_ks_count", ks_count);
	}
	
OUTPUT_SKYMAP("%s_max_upper_limit", max_upper_limit);
OUTPUT_SKYMAP("%s_max_lower_limit", max_lower_limit);
OUTPUT_SKYMAP("%s_arg_freq", freq_map);

snprintf(s,19999,"%s%s_max_dx.dat", subinstance_name, pol->name);
dump_floats(s,pol->skymap.max_dx,fine_grid->npoints,1);

OUTPUT_SKYMAP("%s_S_map", S_map);
OUTPUT_SKYMAP("%s_M_map", M_map);
	
for(i=0;i<fine_grid->npoints;i++){
	if(fine_grid->band[i]<0){
		pol->skymap.max_upper_limit[i]=-1.0;
		pol->skymap.max_lower_limit[i]=-1.0;
		continue;
		}
	pol->skymap.max_upper_limit[i]=sqrt(pol->skymap.max_upper_limit[i])*upper_limit_comp;
	pol->skymap.max_lower_limit[i]=sqrt(pol->skymap.max_lower_limit[i])*lower_limit_comp;
	}

/* output interesting points around fake injection */
largest_i=0;
largest=0.0;
if(fake_injection) {
	double ds, best_ds=10;
	int best_i=-1;
	fprintf(LOG,"Interesting points: index longitude latitude pol max_dx upper_strain lower_strain freq beta1 beta2\n");
	for(i=0;i<fine_grid->npoints;i++){
		/* e[2][i] is just cos of latitude */
		/* Approximate spherical distance */
		#if 1
		ds=sqr_f(fine_grid->latitude[i]-args_info.fake_dec_arg)+
			sqr_f((fine_grid->longitude[i]-args_info.fake_ra_arg)*fine_grid->e[3][i]);
		if(ds<9*resolution*resolution){
		#elif 0
		/* Simplified spherical distance - should be better, different from exact
		   by a monotonic function, but harder to find interesting points.  */
		ds=1.0-(fine_grid->e[2][i]*sin(args_info.fake_dec_arg)+
			fine_grid->e[3][i]*cos(args_info.fake_dec_arg)*
			cos(fine_grid->longitude[i]-args_info.fake_ra_arg));
		if(args_info.fake_dec_arg*fine_grid->latitude[i]<0)ds=1;
		if(ds<9*resolution){
		#else
		/* Exact spherical distance */
		ds=acos(fine_grid->e[2][i]*sin(args_info.fake_dec_arg)+
			fine_grid->e[3][i]*cos(args_info.fake_dec_arg)*
			cos(fine_grid->longitude[i]-args_info.fake_ra_arg));		
		if(args_info.fake_dec_arg*fine_grid->latitude[i]<0)ds=1;
		if(ds<3*resolution){
		#endif
		
		   	fprintf(LOG, "%d %f %f %s %f %g %g %f %f %f\n",
				i,
				fine_grid->longitude[i], fine_grid->latitude[i], 
				pol->name, pol->skymap.max_dx[i], 
				pol->skymap.max_upper_limit[i], pol->skymap.max_lower_limit[i],
				pol->skymap.freq_map[i], 
				args_info.compute_betas_arg?pol->skymap.beta1[i]:NAN,
				args_info.compute_betas_arg?pol->skymap.beta2[i]:NAN);

			if(largest_i<0){
				largest=pol->skymap.max_upper_limit[i];
				largest_i=i;
				} else 
			if(largest<pol->skymap.max_upper_limit[i]){
				largest=pol->skymap.max_upper_limit[i];
				largest_i=i;
				}
		   	}

		if((best_i<0) || (ds<best_ds)){
			best_ds=ds;
			best_i=i;
			}
		
		}
	if(best_i>=0)
	fprintf(LOG, "i_closest: %d %f %f %s %f %g %g %f %f %f\n", best_i, fine_grid->longitude[best_i], fine_grid->latitude[best_i], 
		pol->name, pol->skymap.max_dx[best_i], 
		pol->skymap.max_upper_limit[best_i], pol->skymap.max_lower_limit[best_i], pol->skymap.freq_map[best_i],
		args_info.compute_betas_arg?pol->skymap.beta1[best_i]:NAN, args_info.compute_betas_arg?pol->skymap.beta2[best_i]:NAN);
	if(largest_i>=0)
	fprintf(LOG, "i_largest: %d %f %f %s %f %g %g %f %f %f\n", largest_i, fine_grid->longitude[largest_i], fine_grid->latitude[largest_i], 
		pol->name, pol->skymap.max_dx[largest_i], 
		pol->skymap.max_upper_limit[largest_i], pol->skymap.max_lower_limit[largest_i], pol->skymap.freq_map[largest_i],
		args_info.compute_betas_arg?pol->skymap.beta1[largest_i]:NAN, 
		args_info.compute_betas_arg?pol->skymap.beta2[largest_i]:NAN);
	}

snprintf(s,19999,"%s%s_max_strain.dat", subinstance_name, pol->name);
dump_floats(s,pol->skymap.max_upper_limit,fine_grid->npoints,1);

max_dx=0.0;
max_dx_i=0;
masked=0;
largest_i=0;
largest=0.0;
for(i=0;i<fine_grid->nbands;i++){
	max_band[i]=-1.0;
	masked_max_band[i]=-1.0;
	max_band_arg[i]=-1;
	masked_max_band_arg[i]=-1;
	}
	
for(i=0;i<fine_grid->npoints;i++){
	k=fine_grid->band[i];
	if(k<0)continue;

	if(pol->skymap.max_upper_limit[i]>max_band[k]){
		max_band[k]=pol->skymap.max_upper_limit[i];
		max_band_arg[k]=i;
		}


	if(pol->skymap.max_sub_weight[i]>=pol->skymap.total_weight[i]*(1-args_info.small_weight_ratio_arg)){
		pol->skymap.max_upper_limit[i]=0.0;
		pol->skymap.max_lower_limit[i]=0.0;
		pol->skymap.max_dx[i]=0.0;
		masked++;
		}

	if(pol->skymap.max_upper_limit[i]>largest){
		largest=pol->skymap.max_upper_limit[i];
		largest_i=i;
		}

	if(pol->skymap.max_dx[i]>max_dx){
		max_dx=pol->skymap.max_dx[i];
		max_dx_i=i;
		}

	if(pol->skymap.max_upper_limit[i]>masked_max_band[k]){
		masked_max_band[k]=pol->skymap.max_upper_limit[i];
		masked_max_band_arg[k]=i;
		}
	}
fprintf(LOG, "masked: %s %d\n", pol->name, masked);
fprintf(LOG, "strongest signal: longitude latitude pol max_dx upper_strain lower_strain freq beta1 beta2 spindown\n");	
fprintf(LOG, "max_dx: %f %f %s %f %g %g %f %f %f %g\n",fine_grid->longitude[max_dx_i], fine_grid->latitude[max_dx_i], 
				pol->name, pol->skymap.max_dx[max_dx_i],
				pol->skymap.max_upper_limit[max_dx_i],
				pol->skymap.max_lower_limit[max_dx_i],
				pol->skymap.freq_map[max_dx_i],
				args_info.compute_betas_arg?pol->skymap.beta1[max_dx_i]:NAN,
				args_info.compute_betas_arg?pol->skymap.beta2[max_dx_i]:NAN,
				spindown);

fprintf(LOG, "largest signal: longitude latitude pol max_dx upper_strain lower_strain freq beta1 beta2 spindown\n");	
fprintf(LOG, "largest: %f %f %s %f %g %g %f %f %f %g\n",fine_grid->longitude[largest_i], fine_grid->latitude[largest_i], 
				pol->name, pol->skymap.max_dx[largest_i],
				pol->skymap.max_upper_limit[largest_i],
				pol->skymap.max_lower_limit[largest_i],
				pol->skymap.freq_map[largest_i],
				args_info.compute_betas_arg?pol->skymap.beta1[largest_i]:NAN,
				args_info.compute_betas_arg?pol->skymap.beta2[largest_i]:NAN,
				spindown);

fprintf(LOG, "max/masked band format: band_num longitude latitude pol max_dx upper_strain freq beta1 beta2 spindown\n");
for(i=0;i<fine_grid->nbands;i++){
	if(max_band_arg[i]<0){
		fprintf(LOG, "max_band: %d NAN NAN %s NAN NAN NAN NAN NAN NAN\n", i, pol->name); 
		fprintf(LOG, "masked_max_band: %d NAN NAN %s NAN NAN NAN NAN NAN NAN\n", i, pol->name);
		fprintf(LOG,"max_ratio: %d %s NAN NAN\n", i, pol->name);
		continue;
		}

	fprintf(LOG, "max_band: %d %f %f %s %f %g %f %f %f %g\n", i, fine_grid->longitude[max_band_arg[i]], fine_grid->latitude[max_band_arg[i]], 
				pol->name, pol->skymap.max_dx[max_band_arg[i]], 
				max_band[i], 
				pol->skymap.freq_map[max_band_arg[i]],
				args_info.compute_betas_arg?pol->skymap.beta1[max_band_arg[i]]:NAN, 
				args_info.compute_betas_arg?pol->skymap.beta2[max_band_arg[i]]:NAN,
				spindown);

	fprintf(LOG, "masked_max_band: %d %f %f %s %f %g %f %f %f %g\n", i, fine_grid->longitude[masked_max_band_arg[i]], fine_grid->latitude[masked_max_band_arg[i]],
				pol->name, pol->skymap.max_dx[masked_max_band_arg[i]],
				masked_max_band[i],
				pol->skymap.freq_map[masked_max_band_arg[i]],
				args_info.compute_betas_arg?pol->skymap.beta1[masked_max_band_arg[i]]:NAN,
				args_info.compute_betas_arg?pol->skymap.beta2[masked_max_band_arg[i]]:NAN,
				spindown);

	snprintf(s,19999,"%s%s_max_upper_limit_band_%d.png", subinstance_name, pol->name, i);
	if(clear_name_png(s)){
		adjust_plot_limits_f(plot, freq_f, &(pol->spectral_plot.max_upper_limit[i*useful_bins]), useful_bins, 1, 1, 1);
		draw_grid(p, plot, 0, 0);
		draw_points_f(p, plot, COLOR(255,0,0), freq_f, &(pol->spectral_plot.max_upper_limit[i*useful_bins]), useful_bins, 1, 1);
		RGBPic_dump_png(s, p);
		}
	snprintf(s,19999,"%s%s_max_upper_limit_band_%d.dat", subinstance_name, pol->name, i);
	dump_floats(s, &(pol->spectral_plot.max_upper_limit[i*useful_bins]), useful_bins, 1);
	
	snprintf(s,19999,"%s%s_max_dx_band_%d.png", subinstance_name, pol->name, i);
	if(clear_name_png(s)){
		adjust_plot_limits_f(plot, freq_f, &(pol->spectral_plot.max_dx[i*useful_bins]), useful_bins, 1, 1, 1);
		draw_grid(p, plot, 0, 0);
		draw_points_f(p, plot, COLOR(255,0,0), freq_f, &(pol->spectral_plot.max_dx[i*useful_bins]), useful_bins, 1, 1);
		RGBPic_dump_png(s, p);
		}
	snprintf(s,19999,"%s%s_max_dx_band_%d.dat", subinstance_name, pol->name, i);
	dump_floats(s, &(pol->spectral_plot.max_dx[i*useful_bins]), useful_bins, 1);
	
	snprintf(s,19999,"%s%s_max_mask_ratio_band_%d.png", subinstance_name, pol->name, i);
	if(clear_name_png(s)){
		adjust_plot_limits_f(plot, freq_f, &(pol->spectral_plot.max_mask_ratio[i*useful_bins]), useful_bins, 1, 1, 1);
		draw_grid(p, plot, 0, 0);
		draw_points_f(p, plot, COLOR(255,0,0), freq_f, &(pol->spectral_plot.max_mask_ratio[i*useful_bins]), useful_bins, 1, 1);
		RGBPic_dump_png(s, p);
		}
	snprintf(s,19999,"%s%s_max_mask_ratio_band_%d.dat", subinstance_name, pol->name, i);
	dump_floats(s, &(pol->spectral_plot.max_mask_ratio[i*useful_bins]), useful_bins, 1);

	max_ratio=pol->spectral_plot.max_mask_ratio[i*useful_bins];
	for(k=1;k<useful_bins;k++)
		if(max_ratio<pol->spectral_plot.max_mask_ratio[i*useful_bins+k]){
			max_ratio=pol->spectral_plot.max_mask_ratio[i*useful_bins+k];
			}
	fprintf(LOG,"max_ratio: %d %s %f %g\n", i, pol->name, max_ratio, spindown);
        /* old 
	fprintf(LOG, "max_band: %d %s %g\n", i, pol->name, max_band[i]);
	fprintf(LOG, "masked_max_band: %d %s %g\n", i, pol->name, masked_max_band[i]);
	*/
	}

for(i=0;i<fine_grid->nbands*useful_bins;i++){
	pol->spectral_plot.max_upper_limit[i]=sqrt(pol->spectral_plot.max_upper_limit[i])*upper_limit_comp;
	}

for(i=0;i<fine_grid->nbands;i++){
	snprintf(s,19999,"%s%s_max_upper_strain_band_%d.png", subinstance_name, pol->name, i);
	if(clear_name_png(s)){
		adjust_plot_limits_f(plot, freq_f, &(pol->spectral_plot.max_upper_limit[i*useful_bins]), useful_bins, 1, 1, 1);
		draw_grid(p, plot, 0, 0);
		draw_points_f(p, plot, COLOR(255,0,0), freq_f, &(pol->spectral_plot.max_upper_limit[i*useful_bins]), useful_bins, 1, 1);
		RGBPic_dump_png(s, p);
		}
	snprintf(s,19999,"%s%s_max_upper_strain_band_%d.dat", subinstance_name, pol->name, i);
	dump_floats(s, &(pol->spectral_plot.max_upper_limit[i*useful_bins]), useful_bins, 1);
	}
	
snprintf(s,19999,"%s%s_max_upper_strain.png", subinstance_name, pol->name);
if(clear_name_png(s)){
	plot_grid_f(p, fine_grid, pol->skymap.max_upper_limit, 1);
	RGBPic_dump_png(s, p);
	}
compute_histogram_f(hist, pol->skymap.max_upper_limit, fine_grid->band, fine_grid->npoints);
snprintf(s,19999,"%shist_%s_max_upper_strain", subinstance_name, pol->name);
print_histogram(LOG, hist, s);

snprintf(s,19999,"%s%s_max_lower_strain.png", subinstance_name, pol->name);
if(clear_name_png(s)){
	plot_grid_f(p, fine_grid, pol->skymap.max_lower_limit, 1);
	RGBPic_dump_png(s, p);
	}

snprintf(s,19999,"%s%s_max_dx.png", subinstance_name, pol->name);
if(clear_name_png(s)){
	plot_grid_f(p, fine_grid, pol->skymap.max_dx, 1);
	RGBPic_dump_png(s, p);
	}

fflush(LOG);
free_histogram(hist);
free_plot(plot);
free_RGBPic(p);
free(max_band);
free(masked_max_band);
free(max_band_arg);
free(masked_max_band_arg);
free(freq_f);
}

void output_unified_limits(void)
{
RGBPic *p;
PLOT *plot;
int i,k, m;
char s[20000];
SUM_TYPE *skymap_high_ul, *skymap_high_ul_freq;
SUM_TYPE *spectral_plot_high_ul;
float *freq_f;
SUM_TYPE max_high_ul, max_circ_ul;
int max_high_ul_i, max_circ_ul_i;
HISTOGRAM *hist;
SUM_TYPE a;
SUM_TYPE *max_dx_band;
int *max_dx_band_index;

freq_f=do_alloc(useful_bins, sizeof(*freq_f));
for(i=0;i<useful_bins;i++)freq_f[i]=(first_bin+side_cut+i)/1800.0;

max_dx_band=do_alloc(fine_grid->nbands, sizeof(*max_dx_band));
max_dx_band_index=do_alloc(fine_grid->nbands, sizeof(*max_dx_band_index));
for(i=0;i<fine_grid->nbands;i++) {
	max_dx_band[i]=0;
	max_dx_band_index[i]=-1;
	}

skymap_high_ul=do_alloc(fine_grid->npoints, sizeof(*skymap_high_ul));
skymap_high_ul_freq=do_alloc(fine_grid->npoints, sizeof(*skymap_high_ul_freq));
spectral_plot_high_ul=do_alloc(useful_bins*fine_grid->nbands, sizeof(*spectral_plot_high_ul));
hist=new_histogram(args_info.hist_bins_arg, fine_grid->nbands);

if(fine_grid->max_n_dec<800){
	p=make_RGBPic(fine_grid->max_n_ra*(800/fine_grid->max_n_dec)+140, fine_grid->max_n_dec*(800/fine_grid->max_n_dec));
	} else 
	p=make_RGBPic(fine_grid->max_n_ra+140, fine_grid->max_n_dec);	

plot=make_plot(p->width, p->height);

max_high_ul_i=-1;
max_circ_ul_i=-1;
for(i=0;i<fine_grid->npoints;i++) {
	if(fine_grid->band[i]<0){
		skymap_circ_ul[i]=-1.0;
		skymap_high_ul[i]=-1.0;
		skymap_high_ul_freq[i]=-1.0;
		continue;
		}
	skymap_circ_ul[i]=sqrt(skymap_circ_ul[i])*upper_limit_comp;
	
	skymap_high_ul[i]=polarization_results[0].skymap.max_upper_limit[i];
	skymap_high_ul_freq[i]=polarization_results[0].skymap.freq_map[i];
	for(k=1;k<ntotal_polarizations;k++){	
		if(skymap_high_ul[i]<polarization_results[k].skymap.max_upper_limit[i]){
			skymap_high_ul[i]=polarization_results[k].skymap.max_upper_limit[i];
			skymap_high_ul_freq[i]=polarization_results[k].skymap.freq_map[i];
			}
		}
	if(max_high_ul_i<0){
		max_high_ul_i=i;
		max_high_ul=skymap_high_ul[i];
		} else {
		if(max_high_ul<skymap_high_ul[i]){
			max_high_ul_i=i;
			max_high_ul=skymap_high_ul[i];
			}
		}
	if(max_circ_ul_i<0){
		max_circ_ul_i=i;
		max_circ_ul=skymap_circ_ul[i];
		} else {
		if(max_circ_ul<skymap_circ_ul[i]){
			max_circ_ul_i=i;
			max_circ_ul=skymap_circ_ul[i];
			}
		}
	for(k=0;k<ntotal_polarizations;k++) {
		a=polarization_results[k].skymap.max_dx[i];
		if(a<0)continue;
		if(a>max_dx[i]) {
			max_dx[i]=a;
			max_dx_polarization_index[i]=k;
			}
		}

	if(max_dx_band_index[fine_grid->band[i]]<0 || (max_dx[i]>max_dx_band[fine_grid->band[i]])) {
		max_dx_band_index[fine_grid->band[i]]=i;
		max_dx_band[fine_grid->band[i]]=max_dx[i];
		}
	}

fprintf(LOG,"band SNR: band band_name max_dx pol freq ra dec pt_index spindown\n");
for(i=0;i<fine_grid->nbands;i++) {

	k=max_dx_band_index[i];
	if(k<0) {
		fprintf(LOG, "max_dx_band: %d \"%s\" NaN NaN NaN NaN NaN -1 %g\n",
			i, fine_grid->band_name[i], spindown);
		continue;
		}

	m=max_dx_polarization_index[k];
	if(m<0) {
		fprintf(LOG, "max_dx_band: %d \"%s\" NaN NaN NaN %f %f %d %g\n",
			i, fine_grid->band_name[i],
			fine_grid->longitude[k],
			fine_grid->latitude[k],
			k, spindown);
		continue;
		}

	fprintf(LOG, "max_dx_band: %d \"%s\" %f \"%s\" %f %f %f %d %g\n",
		i, fine_grid->band_name[i], max_dx_band[i], 
		polarization_results[m].name,
		polarization_results[m].skymap.freq_map[k],
		fine_grid->longitude[k],
		fine_grid->latitude[k],
		k,
		spindown);
	}

snprintf(s, 19999, "%smax_dx.png", subinstance_name);
if(clear_name_png(s)){
	plot_grid_f(p, fine_grid, max_dx, 1);
	RGBPic_dump_png(s, p);
	}
snprintf(s, 19999, "%smax_dx.dat", subinstance_name);
dump_floats(s, max_dx, fine_grid->npoints, 1);

if(max_high_ul_i>=0){
	fprintf(LOG, "max_high_ul legend: RA DEC high_ul freq\n");
	fprintf(LOG, "max_high_ul: %f %f %g %f %g\n", 
		fine_grid->longitude[max_high_ul_i],
		fine_grid->latitude[max_high_ul_i],
		max_high_ul,
		skymap_high_ul_freq[max_high_ul_i],
		spindown
		);
	}
if(max_circ_ul_i>=0){
	fprintf(LOG, "max_circ_ul legend: RA DEC circ_ul freq\n");
	fprintf(LOG, "max_circ_ul: %f %f %g %f %g\n", 
		fine_grid->longitude[max_circ_ul_i],
		fine_grid->latitude[max_circ_ul_i],
		max_circ_ul,
		skymap_circ_ul_freq[max_circ_ul_i],
		spindown
		);
	}

if(args_info.compute_betas_arg){
	if(clear_name_png("circ_ul.png")){
		plot_grid_f(p, fine_grid, skymap_circ_ul, 1);
		RGBPic_dump_png("circ_ul.png", p);
		}
	dump_floats("circ_ul.dat", skymap_circ_ul, fine_grid->npoints, 1);
	compute_histogram_f(hist, skymap_circ_ul, fine_grid->band, fine_grid->npoints);
	print_histogram(LOG, hist, "hist_circ_ul");
	}

snprintf(s, 19999, "%shigh_ul.png", subinstance_name);
if(clear_name_png(s)){
	plot_grid_f(p, fine_grid, skymap_high_ul, 1);
	RGBPic_dump_png(s, p);
	}
snprintf(s, 19999, "%shigh_ul.dat", subinstance_name);
dump_floats(s, skymap_high_ul, fine_grid->npoints, 1);
compute_histogram_f(hist, skymap_high_ul, fine_grid->band, fine_grid->npoints);
print_histogram(LOG, hist, "hist_high_ul");


for(i=0;i<useful_bins*fine_grid->nbands;i++){
	spectral_plot_circ_ul[i]=sqrt(spectral_plot_circ_ul[i])*upper_limit_comp;
	
	spectral_plot_high_ul[i]=polarization_results[0].spectral_plot.max_upper_limit[i];
	for(k=1;k<ntotal_polarizations;k++){	
		if(spectral_plot_high_ul[i]<polarization_results[k].spectral_plot.max_upper_limit[i])spectral_plot_high_ul[i]=polarization_results[k].spectral_plot.max_upper_limit[i];
		}
	}

fprintf(LOG,"band upper limits: band UL freq\n");

for(i=0;i<fine_grid->nbands;i++){

	max_high_ul_i=0;
	max_high_ul=spectral_plot_high_ul[i*useful_bins];
	for(k=1;k<useful_bins;k++){
		if(max_high_ul<spectral_plot_high_ul[i*useful_bins+k]){
			max_high_ul_i=k;
			max_high_ul=spectral_plot_high_ul[i*useful_bins+k];
			}
		}
	fprintf(LOG, "max_high_ul_band: %d %g %f %g\n", 
		i, max_high_ul, freq_f[max_high_ul_i], spindown);
	
	max_circ_ul_i=0;
	max_circ_ul=spectral_plot_circ_ul[i*useful_bins];
	for(k=1;k<useful_bins;k++){
		if(max_circ_ul<spectral_plot_circ_ul[i*useful_bins+k]){
			max_circ_ul_i=k;
			max_circ_ul=spectral_plot_circ_ul[i*useful_bins+k];
			}
		}
	fprintf(LOG, "max_circ_ul_band: %d %g %f %g\n", 
		i, max_circ_ul, freq_f[max_circ_ul_i], spindown);
	
	snprintf(s,19999,"%slow_band_%d_ul.png", subinstance_name, i);
	if(clear_name_png(s)){
		adjust_plot_limits_f(plot, freq_f, &(spectral_plot_circ_ul[i*useful_bins]), useful_bins, 1, 1, 1);
		draw_grid(p, plot, 0, 0);
		draw_points_f(p, plot, COLOR(255,0,0), freq_f, &(spectral_plot_circ_ul[i*useful_bins]), useful_bins, 1, 1);
		RGBPic_dump_png(s, p);
		}
	snprintf(s,19999,"%slow_band_%d_ul.dat", subinstance_name, i);
	dump_floats(s, &(spectral_plot_circ_ul[i*useful_bins]), useful_bins, 1);

	snprintf(s,19999,"%shigh_band_%d_ul.png", subinstance_name, i);
	if(clear_name_png(s)){
		adjust_plot_limits_f(plot, freq_f, &(spectral_plot_high_ul[i*useful_bins]), useful_bins, 1, 1, 1);
		draw_grid(p, plot, 0, 0);
		draw_points_f(p, plot, COLOR(255,0,0), freq_f, &(spectral_plot_high_ul[i*useful_bins]), useful_bins, 1, 1);
		RGBPic_dump_png(s, p);
		}
	snprintf(s,19999,"%shigh_band_%d_ul.dat", subinstance_name, i);
	dump_floats(s, &(spectral_plot_high_ul[i*useful_bins]), useful_bins, 1);
	}

free(max_dx_band);
free(max_dx_band_index);
free_histogram(hist);
free_plot(plot);
free_RGBPic(p);
free(skymap_high_ul);
free(freq_f);
}

void compute_mean(int thread_id, ACCUMULATION_ARRAYS *ar, int pi)
{
SUM_TYPE a,c;
int i,k,m;
int offset;
for(k=0,offset=super_grid->first_map[pi];offset>=0;offset=super_grid->list_map[offset],k++){
	if(fine_grid->band[offset]<0)continue;
	
	for(m=0;m<ntotal_polarizations;m++){
		polarization_results[m].skymap.max_sub_weight[offset]=0.0;
		polarization_results[m].skymap.total_weight[offset]=ar[m].total_weight[k];

		#ifdef WEIGHTED_SUM		
		c=(polarization_results[m].skymap.total_weight[offset]);
		#else
		c=(polarization_results[m].skymap.total_count[offset]);
		#endif

		if(c>0){	
			if(args_info.compute_betas_arg){
			    polarization_results[m].skymap.beta1[offset]/=c;
			    polarization_results[m].skymap.beta2[offset]/=c;
			    }

		    for(i=0;i<useful_bins;i++){

			#ifdef WEIGHTED_SUM		
			if(ar[m].fine_grid_weight[i+k*useful_bins]>polarization_results[m].skymap.max_sub_weight[offset]){
				polarization_results[m].skymap.max_sub_weight[offset]=ar[m].fine_grid_weight[i+k*useful_bins];
				}
	
			c=(polarization_results[m].skymap.total_weight[offset]-ar[m].fine_grid_weight[i+k*useful_bins]);
			#else
			c=(polarization_results[m].skymap.total_count[offset]-ar[m].fine_grid_count[i+k*useful_bins]);
			#endif
			if(c>0){
				a=ar[m].fine_grid_sum[i+k*useful_bins];
				#ifdef COMPUTE_SIGMA
				b=ar[m].fine_grid_sq_sum[i+k*useful_bins];
				#endif
		
				ar[m].fine_grid_sum[i+k*useful_bins]=a/c;

				#ifdef COMPUTE_SIGMA
				ar[m].fine_grid_sq_sum[i+k*useful_bins]=sqrt((b*ar[m].fine_grid_count[i+k*useful_bins]-a)/c);
				#endif
				}
			}
		   }
		}
	}
}

void compute_mean_no_lines(int thread_id, ACCUMULATION_ARRAYS *ar, int pi)
{
SUM_TYPE a,c;
int i,k,offset,m;
for(k=0,offset=super_grid->first_map[pi];offset>=0;offset=super_grid->list_map[offset],k++){
	if(fine_grid->band[offset]<0)continue;

	for(m=0;m<ntotal_polarizations;m++){
		polarization_results[m].skymap.max_sub_weight[offset]=0.0;
		polarization_results[m].skymap.total_weight[offset]=ar[m].total_weight[k];
		
		#ifdef WEIGHTED_SUM		
		c=(polarization_results[m].skymap.total_weight[offset]);
		#else
		c=(polarization_results[m].skymap.total_count[offset]);
		#endif


		if(c>0){
			if(args_info.compute_betas_arg){
				polarization_results[m].skymap.beta1[offset]/=c;
				polarization_results[m].skymap.beta2[offset]/=c;
				}

			for(i=0;i<useful_bins;i++){
				a=ar[m].fine_grid_sum[i+k*useful_bins];
				#ifdef COMPUTE_SIGMA
				b=ar[m].fine_grid_sq_sum[i+k*useful_bins];
				#endif
		
				ar[m].fine_grid_sum[i+k*useful_bins]=a/c;

				#ifdef COMPUTE_SIGMA
				ar[m].fine_grid_sq_sum[i+k*useful_bins]=sqrt((b*ar[m].fine_grid_count[i+k*useful_bins]-a)/c);
				#endif
				}
			}
		}
	}
}

void init_fine_grid_stage(void)
{
init_fc_ul();
init_fc_ll();
verify_limits();

normalizing_weight=datasets_normalizing_weight();

if(!strcasecmp(args_info.averaging_mode_arg, "matched")) {
       process_patch=process_patch_matched;
       fprintf(LOG,"mode: matched filter\n");
       } else
if(!strcasecmp(args_info.averaging_mode_arg, "3") || !strcasecmp(args_info.averaging_mode_arg, "three")){
       process_patch=process_patch3;
       fprintf(LOG,"mode: 3 bins\n");
       } else
       {
       process_patch=process_patch1;
       fprintf(LOG,"mode: 1 bin\n");
       }


if(!strcasecmp("Hann", args_info.upper_limit_comp_arg)){
	if(!strcasecmp(args_info.averaging_mode_arg, "matched")) {
		upper_limit_comp=1.0; /* TODO - make sure this number is right */
		} else
	if(!strcasecmp(args_info.averaging_mode_arg, "3") || !strcasecmp(args_info.averaging_mode_arg, "three")){
		/* 3 bins should contain the entire signal, regardless
		   of positioning */
		upper_limit_comp=sqrt(3.0);
		} else 
		{
		/* 0.85 is a ratio between amplitude of 
		   half-bin centered signal and bin centered signal
		   *amplitude*

		   */
		upper_limit_comp=1.0/0.85;
		}
	} else {
	upper_limit_comp=atof(args_info.upper_limit_comp_arg);
	}
fprintf(LOG, "upper limit compensation factor: %8f\n", upper_limit_comp);

	/* Extra factor to convert to amplitude from RMS power */
upper_limit_comp*=sqrt(2.0);
	/* Extra factor to convert to strain from raw SFT units */
upper_limit_comp/=(1800.0*16384.0);
	/* Extra factor to account for the fact that only half of SFT
	   coefficients is stored */
upper_limit_comp*=sqrt(2.0);
	/* Revert strain normalization */
upper_limit_comp*=args_info.strain_norm_factor_arg;

if(!strcasecmp("Hann", args_info.lower_limit_comp_arg)){
	if(!strcasecmp(args_info.averaging_mode_arg, "matched")) {
		lower_limit_comp=1.0; /* TODO - make sure this number is right */
		} else
	if(!strcasecmp(args_info.averaging_mode_arg, "3") || !strcasecmp(args_info.averaging_mode_arg, "three")){
		lower_limit_comp=sqrt(3.0);
		} else 
		{
		lower_limit_comp=1.0;
		}
	} else {
	lower_limit_comp=atof(args_info.lower_limit_comp_arg);
	}
fprintf(LOG, "lower limit compensation factor: %8f\n", lower_limit_comp);

	/* Extra factor to  convert to amplitude from RMS power */
lower_limit_comp*=sqrt(2.0);
	/* Extra factor to convert to strain from raw SFT units */
lower_limit_comp/=(1800.0*16384.0);
	/* Extra factor to account for the fact that only half of SFT
	   coefficients is stored */
lower_limit_comp*=sqrt(2.0);
	/* Revert strain normalization */
lower_limit_comp*=args_info.strain_norm_factor_arg;

}

void fine_grid_allocate_arrays(void)
{
stored_fine_bins=super_grid->max_npatch;

allocate_polarization_arrays();

accum_circ_ul=do_alloc(stored_fine_bins*useful_bins*get_max_threads(), sizeof(*accum_circ_ul));
accum_circ_ul_freq=do_alloc(stored_fine_bins*useful_bins*get_max_threads(), sizeof(*accum_circ_ul_freq));
skymap_circ_ul=do_alloc(fine_grid->npoints, sizeof(*skymap_circ_ul));
skymap_circ_ul_freq=do_alloc(fine_grid->npoints, sizeof(*skymap_circ_ul_freq));
spectral_plot_circ_ul=do_alloc(useful_bins*fine_grid->nbands, sizeof(*spectral_plot_circ_ul));
max_dx=do_alloc(fine_grid->npoints, sizeof(*max_dx));
max_dx_polarization_index=do_alloc(fine_grid->npoints, sizeof(*max_dx_polarization_index));
}

void fine_grid_free_arrays(void)
{
free_polarization_arrays();

#define FREE(var) { free(var); var=NULL; }

FREE(accum_circ_ul);
FREE(accum_circ_ul_freq);
FREE(skymap_circ_ul);
FREE(skymap_circ_ul_freq);
FREE(spectral_plot_circ_ul);
FREE(max_dx);
FREE(max_dx_polarization_index);
}

void fine_grid_cruncher(int thread_id, void *data)
{
ACCUMULATION_ARRAYS *ar;
int i, j, k, m;
double a, b;
int pi=(long)data;
SUM_TYPE *circ_ul, *circ_ul_freq;

circ_ul=&(accum_circ_ul[stored_fine_bins*useful_bins*(thread_id+1)]);
circ_ul_freq=&(accum_circ_ul[stored_fine_bins*useful_bins*(thread_id+1)]);

ar=get_thread_accumulation_arrays(thread_id);

clear_accumulation_arrays(ar);

/* loop over datasets */
for(j=0;j<d_free;j++) {

	/* process single patch */
	for(k=0;k<datasets[j].free;k++) {
		if(datasets[j].sft_veto[k])continue;

		a=datasets[j].expTMedians[k];
		
		for(m=0;m<ntotal_polarizations;m++) {
			b=datasets[j].polarizations[m].patch_CutOff[pi];
			/* process polarization */
			if(!do_CutOff || (b*a*AM_response(k, patch_grid, pi, datasets[j].polarizations[m].AM_coeffs)<4))
				process_patch(thread_id, &(datasets[j]), ar, m, pi, k, b*sqrt(a));
			}

		/* accumulate intermediate results */
		if((k % 100)==0) {
			update_d_accumulation_arrays(ar);
			}
		}
	
	}

finalize_accumulation_arrays(ar);

/* compute means */
if(0)compute_mean_no_lines(thread_id, ar, pi);
	else compute_mean(thread_id, ar, pi);
	
for(i=0;i<stored_fine_bins*useful_bins;i++){
	circ_ul[i]=1.0e23; /* Sufficiently large number, 
				even for SFTs done with 
				make_fake_data */
	}


/* compute upper limits */
for(i=0;i<ntotal_polarizations;i++){
	make_limits(thread_id, &(polarization_results[i]), &(ar[i]), pi);
	}
make_unified_limits(thread_id, pi);
}

void fine_grid_stage(void)
{
int pi,i,j, k,m, last_pi;
double a,b;
struct timeval start_time, end_time, last_time;
double delta, delta_avg;

clear_polarization_arrays();

min_shift=0;
max_shift=0;

for(i=0;i<fine_grid->npoints;i++) {
	skymap_circ_ul[i]=-1.0;
	skymap_circ_ul_freq[i]=-1.0;
	max_dx[i]=0;
	max_dx_polarization_index[i]=-1;
	}
for(i=0;i<useful_bins*fine_grid->nbands;i++){
	spectral_plot_circ_ul[i]=-1.0;
	}
	
fprintf(stderr,"Main loop: %d patches to process.\n", patch_grid->npoints);
gettimeofday(&start_time, NULL);
last_time=start_time;
last_pi=0;
reset_jobs_done_ratio();

for(pi=0;pi<patch_grid->npoints;pi++) {
	if(patch_grid->band[pi]<0)continue;
	if(super_grid->first_map[pi]<0)continue;

	submit_job(fine_grid_cruncher, (void *)((long)pi));

	/* For debugging only: */
	#if 0
	if(pi>300)break;
	#endif
	}

k=0;
while(do_single_job(-1)) {
	if(k % 100 == 0)fprintf(stderr, "% 3.1f ", jobs_done_ratio()*100);
	k++;
	}
wait_for_all_done();
fprintf(stderr, "\n");

gettimeofday(&end_time, NULL);

delta=end_time.tv_sec-last_time.tv_sec+(end_time.tv_usec-last_time.tv_usec)*1e-6;
if(fabs(delta)==0.0)delta=1.0;

fprintf(stderr, "Patch speed: %f\n", (pi-last_pi)/delta);
fprintf(LOG, "Patch speed: %f\n", (pi-last_pi)/delta);

for(i=0;i<get_max_threads();i++) {
	fprintf(stderr, "Power cache %d hits: %lld\n", i, power_cache[i].total_hits);
	fprintf(stderr, "Power cache %d misses: %lld\n", i, power_cache[i].total_misses);
	
	fprintf(LOG, "Power cache %d hits: %lld\n", i, power_cache[i].total_hits);
	fprintf(LOG, "Power cache %d misses: %lld\n", i, power_cache[i].total_misses);
	
	/* reset power cache */
	power_cache[i].free=0;
	}

fprintf(LOG,"Maximum bin shift: %d\n", max_shift);
fprintf(LOG,"Minimum bin shift: %d\n", min_shift);
fflush(LOG);

fprintf(stderr,"Maximum bin shift is %d\n", max_shift);
fprintf(stderr,"Minimum bin shift is %d\n", min_shift);

fprintf(stderr, "Writing polarization specific results\n");
for(i=0;i<ntotal_polarizations;i++){
	output_limits(&(polarization_results[i]));
	}
	
fprintf(stderr, "Writing unified results\n");
output_unified_limits();

fflush(LOG);

if(!args_info.no_candidates_arg)identify_candidates();
}

