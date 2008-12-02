#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/* We need this define to get NAN values */
#define __USE_ISOC99
#include <math.h>

#include "polarization.h"
#include "dataset.h"
#include "power_cache.h"
#include "cmdline.h"

extern DATASET *datasets;
extern int d_free;
extern int nbins, first_bin, side_cut, useful_bins;
extern struct gengetopt_args_info args_info;

extern FILE *LOG;


SEGMENT_INFO *find_segments(double gps_start, double gps_end, int veto_mask, int *count)
{
SEGMENT_INFO *r, *p;
int length=(gps_end-gps_start)/450+1; /* this is enough for two ifos */
int i, k;
DATASET *d;

*count=0;
r=do_alloc(length, sizeof(*r));

for(k=0;k<d_free;k++) {
	d=&(datasets[k]);

	for(i=0;i<d->free;i++) {
		if(d->sft_veto[i] & ~veto_mask)continue;
		if(d->gps[i]<gps_start)continue;
		if(d->gps[i]>=gps_end)continue;

		if(*count>=length) {
			length*=2;
			p=do_alloc(length, sizeof(*p));
			if(*count>0)memcpy(p, r, (*count) * sizeof(*p));
			free(r);
			r=p;
			}

		r[*count].gps=d->gps[i];
		r[*count].detector_velocity[0]=d->detector_velocity[3*i];
		r[*count].detector_velocity[1]=d->detector_velocity[3*i+1];
		r[*count].detector_velocity[2]=d->detector_velocity[3*i+2];
		r[*count].coherence_time=d->coherence_time;
		r[*count].dataset=k;
		r[*count].segment=i;
		r[*count].ra=NAN;
		r[*count].dec=NAN;
		r[*count].freq_shift=NAN;
		
		(*count)++;
		}
	}
return(r);
}

#define REAL double
#define SUFFIX(a) a
#include "partial_power_sum.c"

#undef REAL
#undef SUFFIX

#define REAL float
#define SUFFIX(a) a##_F
#include "partial_power_sum.c"

#undef REAL
#undef SUFFIX

void accumulate_partial_power_sum_F1(PARTIAL_POWER_SUM_F *accum, PARTIAL_POWER_SUM *partial)
{
int i;
if(accum->type!=sizeof(float)) {
	fprintf(stderr, "*** INTERNAL ERROR: %s power sum real type mismatch type=%d sizeof(REAL)=%ld\n",
		__FUNCTION__,
		accum->type, sizeof(float));
	exit(-1);
	}
if(partial->type!=sizeof(double)) {
	fprintf(stderr, "*** INTERNAL ERROR: %s power sum real type mismatch type=%d sizeof(REAL)=%ld\n",
		__FUNCTION__,
		partial->type, sizeof(double));
	exit(-1);
	}

for(i=0;i<useful_bins;i++) {
	accum->power_pp[i]+=partial->power_pp[i];
	accum->power_pc[i]+=partial->power_pc[i];
	accum->power_cc[i]+=partial->power_cc[i];
	}

if(partial->weight_arrays_non_zero) {
	for(i=0;i<useful_bins;i++) {
		accum->weight_pppp[i]+=partial->weight_pppp[i];
		accum->weight_pppc[i]+=partial->weight_pppc[i];
		accum->weight_ppcc[i]+=partial->weight_ppcc[i];
		accum->weight_pccc[i]+=partial->weight_pccc[i];
		accum->weight_cccc[i]+=partial->weight_cccc[i];
		}
	accum->weight_arrays_non_zero=1;
	}

accum->c_weight_pppp+=partial->c_weight_pppp;
accum->c_weight_pppc+=partial->c_weight_pppc;
accum->c_weight_ppcc+=partial->c_weight_ppcc;
accum->c_weight_pccc+=partial->c_weight_pccc;
accum->c_weight_cccc+=partial->c_weight_cccc;
}

void accumulate_partial_power_sum_F2(PARTIAL_POWER_SUM *accum, PARTIAL_POWER_SUM_F *partial)
{
int i;
if(accum->type!=sizeof(double)) {
	fprintf(stderr, "*** INTERNAL ERROR: %s power sum real type mismatch type=%d sizeof(REAL)=%ld\n",
		__FUNCTION__,
		accum->type, sizeof(double));
	exit(-1);
	}
if(partial->type!=sizeof(float)) {
	fprintf(stderr, "*** INTERNAL ERROR: %s power sum real type mismatch type=%d sizeof(REAL)=%ld\n",
		__FUNCTION__,
		partial->type, sizeof(float));
	exit(-1);
	}

for(i=0;i<useful_bins;i++) {
	accum->power_pp[i]+=partial->power_pp[i];
	accum->power_pc[i]+=partial->power_pc[i];
	accum->power_cc[i]+=partial->power_cc[i];
	}

if(partial->weight_arrays_non_zero) {
	for(i=0;i<useful_bins;i++) {
		accum->weight_pppp[i]+=partial->weight_pppp[i];
		accum->weight_pppc[i]+=partial->weight_pppc[i];
		accum->weight_ppcc[i]+=partial->weight_ppcc[i];
		accum->weight_pccc[i]+=partial->weight_pccc[i];
		accum->weight_cccc[i]+=partial->weight_cccc[i];
		}
	accum->weight_arrays_non_zero=1;
	}

accum->c_weight_pppp+=partial->c_weight_pppp;
accum->c_weight_pppc+=partial->c_weight_pppc;
accum->c_weight_ppcc+=partial->c_weight_ppcc;
accum->c_weight_pccc+=partial->c_weight_pccc;
accum->c_weight_cccc+=partial->c_weight_cccc;
}

void get_uncached_single_bin_power_sum(SEGMENT_INFO *si, int count, PARTIAL_POWER_SUM *pps)
{
int i,k,n,m;
int bin_shift;
SEGMENT_INFO *si_local;
DATASET *d;
POLARIZATION *pl;
float *im, *re, *pp, *pc, *cc, *fm;
float *power;
float a;
float TM;
float weight;
float f_plus, f_cross;

float pmax;
float sum, sum_sq;
float useful_bins_inv=1.0/useful_bins;
float pmax_factor=args_info.power_max_median_factor_arg;

double weight_pppp=0;
double weight_pppc=0;
double weight_ppcc=0;
double weight_pccc=0;
double weight_cccc=0;

pp=aligned_alloca(useful_bins*sizeof(*pp));
pc=aligned_alloca(useful_bins*sizeof(*pc));
cc=aligned_alloca(useful_bins*sizeof(*cc));
power=aligned_alloca(useful_bins*sizeof(*power));

for(i=0;i<useful_bins;i++) {
	pp[i]=0.0;
	pc[i]=0.0;
	cc[i]=0.0;
	}

pps->weight_arrays_non_zero=0;

for(i=0;i<useful_bins;i++) {
	pps->weight_pppp[i]=0.0;
	pps->weight_pppc[i]=0.0;
	pps->weight_ppcc[i]=0.0;
	pps->weight_pccc[i]=0.0;
	pps->weight_cccc[i]=0.0;
	}

for(k=0;k<count;k++) {
	si_local=&(si[k]);

	bin_shift=rint(si_local->freq_shift*si_local->coherence_time);
	if((bin_shift+side_cut<0) || (bin_shift>useful_bins+side_cut)) {
		fprintf(stderr, "*** Attempt to sample outside loaded range bin_shift=%d freq_shift=%lg, aborting\n", 
			bin_shift, si_local->freq_shift);
		exit(-1);
		}

	d=&(datasets[si_local->dataset]);
	pl=&(d->polarizations[0]);

	f_plus=F_plus_coeff(k, si_local->e, pl->AM_coeffs);
	f_cross=F_plus_coeff(k, si_local->e, pl->conjugate->AM_coeffs);


	re=&(d->re[si_local->segment*nbins+side_cut+bin_shift]);
	im=&(d->im[si_local->segment*nbins+side_cut+bin_shift]);
	fm=&(d->expFMedians_plain[side_cut+bin_shift]);

	if(args_info.subtract_background_arg) {
		TM=d->expTMedians[si_local->segment];
		} else {
		TM=0;
		}

	pmax=0.0;
	sum=0.0;
	sum_sq=0.0;
	for(i=0;i<useful_bins;i++) {
		a=((*re)*(*re)+(*im)*(*im)-TM*(*fm));
		if(a>pmax)pmax=a;
		power[i]=a;
		sum+=a;
		sum_sq+=a*a;

		im++;
		re++;		
		fm++;
		}

	if(args_info.tmedian_noise_level_arg) {
		weight=d->expTMedians[si_local->segment]*d->weight;
		} else {
		sum*=useful_bins_inv;
		sum_sq*=useful_bins_inv;
		sum_sq-=sum*sum;
	// 	weight=sum_sq;
	// 
	// 	pmax*=pmax_factor; /* scale factor to bring it down to power median */
	// 	a=pmax*pmax;
	// 	//if(a>weight)weight=a;
	// 
		weight=1.0/weight;
		}

	weight_pppp+=weight*f_plus*f_plus*f_plus*f_plus;
	weight_pppc+=weight*f_plus*f_plus*f_plus*f_cross;
	weight_ppcc+=weight*f_plus*f_plus*f_cross*f_cross;
	weight_pccc+=weight*f_plus*f_cross*f_cross*f_cross;
	weight_cccc+=weight*f_cross*f_cross*f_cross*f_cross;

	for(i=0;i<useful_bins;i++) {
		a=power[i]*weight;

		pp[i]+=a*f_plus*f_plus;
		pc[i]+=a*f_plus*f_cross;
		cc[i]+=a*f_cross*f_cross;
		/**/
		}

	for(n=0;(d->lines_report->lines_list[n]>=0)&&(n<d->lines_report->nlines);n++) {
		m=d->lines_report->lines_list[n];
		i=m-side_cut-bin_shift;
		if(i<0)continue;
		if(i>=useful_bins)continue;

		a=power[i]*weight;

		pp[i]-=a*f_plus*f_plus;
		pc[i]-=a*f_plus*f_cross;
		cc[i]-=a*f_cross*f_cross;

		pps->weight_pppp[i]-=weight*f_plus*f_plus*f_plus*f_plus;
		pps->weight_pppc[i]-=weight*f_plus*f_plus*f_plus*f_cross;
		pps->weight_ppcc[i]-=weight*f_plus*f_plus*f_cross*f_cross;
		pps->weight_pccc[i]-=weight*f_plus*f_cross*f_cross*f_cross;
		pps->weight_cccc[i]-=weight*f_cross*f_cross*f_cross*f_cross;		

		pps->weight_arrays_non_zero=1;
		}

	}

pps->c_weight_pppp=weight_pppp;
pps->c_weight_pppc=weight_pppc;
pps->c_weight_ppcc=weight_ppcc;
pps->c_weight_pccc=weight_pccc;
pps->c_weight_cccc=weight_cccc;

for(i=0;i<useful_bins;i++) {
	pps->power_pp[i]=pp[i];
	pps->power_pc[i]=pc[i];
	pps->power_cc[i]=cc[i];
	}

pps->collapsed_weight_arrays=0;
}

struct {
	long hits;
	long misses;
	int count;
	SEGMENT_INFO *si;
	PARTIAL_POWER_SUM *pps;
	} SIMPLE_CACHE={0, 0, 0, NULL, NULL};

void accumulate_single_bin_power_sum_cached1(SEGMENT_INFO *si, int count, PARTIAL_POWER_SUM *pps)
{
int i;
int match=1;
if(count!=SIMPLE_CACHE.count) match=0;
if(SIMPLE_CACHE.pps==NULL)match=0;
if(match) {
	for(i=0;i<count;i++)
		if(abs(si[i].freq_shift-SIMPLE_CACHE.si[i].freq_shift)>(0.5/1800) || si[i].gps!=SIMPLE_CACHE.si[i].gps) {
			match=0;
			break;
			}
	}
if(!match) {
	SIMPLE_CACHE.misses++;

	if(SIMPLE_CACHE.pps==NULL) {
		SIMPLE_CACHE.pps=allocate_partial_power_sum();
		}
	if(SIMPLE_CACHE.count!=count) {
		if(SIMPLE_CACHE.si!=NULL)free(SIMPLE_CACHE.si);
		SIMPLE_CACHE.si=NULL;
		SIMPLE_CACHE.count=count;
		}
	if(SIMPLE_CACHE.si==NULL) {
		SIMPLE_CACHE.si=do_alloc(count, sizeof(*SIMPLE_CACHE.si));
		}
	memcpy(SIMPLE_CACHE.si, si, count*sizeof(*si));
	for(i=0;i<count;i++)SIMPLE_CACHE.si[i].freq_shift=round(SIMPLE_CACHE.si[i].freq_shift*1800)/1800.0;
	get_uncached_single_bin_power_sum(SIMPLE_CACHE.si, SIMPLE_CACHE.count, SIMPLE_CACHE.pps);
	} else {
	SIMPLE_CACHE.hits++;
	}

accumulate_partial_power_sum(pps, SIMPLE_CACHE.pps);
}

void print_cache_stats(void)
{
fprintf(stderr, "SIMPLE_CACHE stats: hits=%ld misses=%ld miss ratio %f\n", SIMPLE_CACHE.hits, SIMPLE_CACHE.misses,
	SIMPLE_CACHE.misses/(SIMPLE_CACHE.hits+SIMPLE_CACHE.misses+0.0));

fprintf(LOG, "SIMPLE_CACHE stats: hits=%ld misses=%ld miss ratio %f\n", SIMPLE_CACHE.hits, SIMPLE_CACHE.misses,
	SIMPLE_CACHE.misses/(SIMPLE_CACHE.hits+SIMPLE_CACHE.misses+0.0));
}
