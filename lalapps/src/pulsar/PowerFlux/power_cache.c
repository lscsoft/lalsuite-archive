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
		r[*count].f_plus=NAN;
		r[*count].f_cross=NAN;
		r[*count].bin_shift=NAN;
		
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
int pps_bins=accum->nbins;
int shift;
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
shift=partial->offset-accum->offset;
if( (shift<0) || (shift+pps_bins>partial->nbins)) {
	fprintf(stderr, "*** INTERNAL ERROR: %s power sums do not overlap shift=%d pps_bins=%d partial->offset=%d accum->offset=%d partial->nbins=%d\n",
		__FUNCTION__,
		shift, pps_bins, partial->offset, accum->offset, partial->nbins);
	exit(-1);
	}

for(i=0;i<pps_bins;i++) {
	accum->power_pp[i]+=partial->power_pp[i+shift];
	accum->power_pc[i]+=partial->power_pc[i+shift];
	accum->power_cc[i]+=partial->power_cc[i+shift];
	}

if(partial->weight_arrays_non_zero) {
	for(i=0;i<pps_bins;i++) {
		accum->weight_pppp[i]+=partial->weight_pppp[i+shift];
		accum->weight_pppc[i]+=partial->weight_pppc[i+shift];
		accum->weight_ppcc[i]+=partial->weight_ppcc[i+shift];
		accum->weight_pccc[i]+=partial->weight_pccc[i+shift];
		accum->weight_cccc[i]+=partial->weight_cccc[i+shift];
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
int pps_bins=accum->nbins;
int shift;
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
shift=partial->offset-accum->offset;
if( (shift<0) || (shift+pps_bins>partial->nbins)) {
	fprintf(stderr, "*** INTERNAL ERROR: %s power sums do not overlap shift=%d pps_bins=%d partial->offset=%d accum->offset=%d partial->nbins=%d\n",
		__FUNCTION__,
		shift, pps_bins, partial->offset, accum->offset, partial->nbins);
	exit(-1);
	}

for(i=0;i<pps_bins;i++) {
	accum->power_pp[i]+=partial->power_pp[i+shift];
	accum->power_pc[i]+=partial->power_pc[i+shift];
	accum->power_cc[i]+=partial->power_cc[i+shift];
	}

if(partial->weight_arrays_non_zero) {
	for(i=0;i<pps_bins;i++) {
		accum->weight_pppp[i]+=partial->weight_pppp[i+shift];
		accum->weight_pppc[i]+=partial->weight_pppc[i+shift];
		accum->weight_ppcc[i]+=partial->weight_ppcc[i+shift];
		accum->weight_pccc[i]+=partial->weight_pccc[i+shift];
		accum->weight_cccc[i]+=partial->weight_cccc[i+shift];
		}
	accum->weight_arrays_non_zero=1;
	}

accum->c_weight_pppp+=partial->c_weight_pppp;
accum->c_weight_pppc+=partial->c_weight_pppc;
accum->c_weight_ppcc+=partial->c_weight_ppcc;
accum->c_weight_pccc+=partial->c_weight_pccc;
accum->c_weight_cccc+=partial->c_weight_cccc;
}

void get_uncached_single_bin_power_sum(SEGMENT_INFO *si, int count, PARTIAL_POWER_SUM_F *pps)
{
int i,k,n,m;
int bin_shift;
SEGMENT_INFO *si_local;
DATASET *d;
//POLARIZATION *pl;
float *im, *re, *pp, *pc, *cc, *fm;
float *power, *p;
float a;
float TM;
float weight;
float f_plus, f_cross, f_pp, f_pc, f_cc;

int pps_bins=pps->nbins;

float pmax;
float sum, sum_sq;
float pps_bins_inv=1.0/pps_bins;
//float pmax_factor=args_info.power_max_median_factor_arg;

float weight_pppp=0;
float weight_pppc=0;
float weight_ppcc=0;
float weight_pccc=0;
float weight_cccc=0;

/*pp=aligned_alloca(pps_bins*sizeof(*pp));
pc=aligned_alloca(pps_bins*sizeof(*pc));
cc=aligned_alloca(pps_bins*sizeof(*cc));*/
power=aligned_alloca(pps_bins*sizeof(*power));

pp=pps->power_pp;
pc=pps->power_pc;
cc=pps->power_cc;
for(i=0;i<pps_bins;i++) {
	(*pp)=0.0;
	(*pc)=0.0;
	(*cc)=0.0;

	pp++;
	pc++;
	cc++;
	}

pps->weight_arrays_non_zero=0;

for(i=0;i<pps_bins;i++) {
	pps->weight_pppp[i]=0.0;
	pps->weight_pppc[i]=0.0;
	pps->weight_ppcc[i]=0.0;
	pps->weight_pccc[i]=0.0;
	pps->weight_cccc[i]=0.0;
	}

//fprintf(stderr, "%d\n", count);

for(k=0;k<count;k++) {
	si_local=&(si[k]);

	bin_shift=rintf(si_local->bin_shift)-pps->offset;
	if((bin_shift+side_cut<0) || (bin_shift>pps_bins+side_cut)) {
		fprintf(stderr, "*** Attempt to sample outside loaded range bin_shift=%d bin_shift=%lg, aborting\n", 
			bin_shift, si_local->bin_shift);
		exit(-1);
		}

 	d=&(datasets[si_local->dataset]);
// 	pl=&(d->polarizations[0]);

	//f_plus=F_plus_coeff(si_local->segment, si_local->e, pl->AM_coeffs);
	//f_cross=F_plus_coeff(si_local->segment, si_local->e, pl->conjugate->AM_coeffs);

	f_plus=si_local->f_plus;
	f_cross=si_local->f_cross;


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
	p=power;
	for(i=0;i<pps_bins;i++) {
		a=((*re)*(*re)+(*im)*(*im)-TM*(*fm));
		//if(a>pmax)pmax=a;
		(*p)=a;
		sum+=a;
		sum_sq+=a*a;

		im++;
		re++;		
		fm++;
		p++;
		}

	if(args_info.tmedian_noise_level_arg) {
		weight=d->expTMedians[si_local->segment]*d->weight;
		} else {
		sum*=pps_bins_inv;
		sum_sq*=pps_bins_inv;
		sum_sq-=sum*sum;
	// 
	// 	pmax*=pmax_factor; /* scale factor to bring it down to power median */
	// 	a=pmax*pmax;
	// 	//if(a>weight)weight=a;
	// 
		weight=1.0/sum_sq;
		}

	f_pp=f_plus*f_plus;
	f_pc=f_plus*f_cross;
	f_cc=f_cross*f_cross;
	

	weight_pppp+=weight*f_pp*f_pp;
	weight_pppc+=weight*f_pp*f_pc;
	weight_ppcc+=weight*f_pp*f_cc;
	weight_pccc+=weight*f_pc*f_cc;
	weight_cccc+=weight*f_cc*f_cc;

	f_pp*=weight;
	f_pc*=weight;
	f_cc*=weight;

	p=power;
	pp=pps->power_pp;
	pc=pps->power_pc;
	cc=pps->power_cc;
	for(i=0;i<pps_bins;i++) {
		a=(*p);

		(*pp)+=a*f_pp;
		(*pc)+=a*f_pc;
		(*cc)+=a*f_cc;

		p++;
		pp++;
		pc++;
		cc++;
		/**/
		}

	pp=pps->power_pp;
	pc=pps->power_pc;
	cc=pps->power_cc;
	for(n=0;(d->lines_report->lines_list[n]>=0)&&(n<d->lines_report->nlines);n++) {
		m=d->lines_report->lines_list[n];
		i=m-side_cut-bin_shift;
		if(i<0)continue;
		if(i>=pps_bins)continue;

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

pps->collapsed_weight_arrays=0;
}

struct {
	long hits;
	long misses;
	long overwrites;
	long large_shifts;
	int segment_count;
	int size;
	int free;
	int *key;
	SEGMENT_INFO **si;
	PARTIAL_POWER_SUM_F **pps;
	} SIMPLE_CACHE={0, 0, 0, 0, 0, 0, 0, NULL, NULL, NULL};

void reset_simple_cache(int segment_count, int template_count)
{
int i;
for(i=0;i<SIMPLE_CACHE.free;i++) {
	free_partial_power_sum_F(SIMPLE_CACHE.pps[i]);
	free(SIMPLE_CACHE.si[i]);
	}
free(SIMPLE_CACHE.pps);
free(SIMPLE_CACHE.si);
free(SIMPLE_CACHE.key);

SIMPLE_CACHE.segment_count=segment_count;
SIMPLE_CACHE.size=template_count;
SIMPLE_CACHE.free=0;
SIMPLE_CACHE.key=do_alloc(SIMPLE_CACHE.size, sizeof(*SIMPLE_CACHE.key));
SIMPLE_CACHE.pps=do_alloc(SIMPLE_CACHE.size, sizeof(*SIMPLE_CACHE.pps));
SIMPLE_CACHE.si=do_alloc(SIMPLE_CACHE.size, sizeof(*SIMPLE_CACHE.si));
}

/* Note: si is modified in place to have bin_shift rounded as appropriate to the power generating algorithm */
void accumulate_single_bin_power_sum_cached1(SEGMENT_INFO *si, int count, PARTIAL_POWER_SUM_F *pps)
{
int i, k;
int match;
int key;
float a;
int first_shift;
int max_shift=args_info.max_first_shift_arg;
SEGMENT_INFO *si_local, *sc_si_local;
if(count!=SIMPLE_CACHE.segment_count) { 
	fprintf(stderr, "Internal error: segment counts don't match.\n");
	exit(-1);
	}

si_local=si;
key=0.0;
k=rintf(si_local->bin_shift);
for(i=0;i<count;i++) {
	//fprintf(stderr, "%0.1f ", si_local->bin_shift);
	a=rintf(si_local->bin_shift);
	key+=((int)a-k)*2*i;
	si_local->bin_shift=a;
	si_local++;
	//fprintf(stderr, "%0.1f ", a);
	}
//fprintf(stderr, "%d ", key);

for(k=0;k<SIMPLE_CACHE.free;k++) {
	/* the reason we use exact equality for floating point numbers is because the number in cache have been placed there by the same function. */
	if(key==SIMPLE_CACHE.key[k]) {
		/* we found the box holding our data, double check it is the right one */
		si_local=si;
		sc_si_local=SIMPLE_CACHE.si[k];
		first_shift=si_local->bin_shift-sc_si_local->bin_shift;
		if( (first_shift>max_shift) || (first_shift< -max_shift)) {
			SIMPLE_CACHE.large_shifts++;
			break;
			}
		match=1;
		for(i=1;i<count;i++) {
			si_local++;
			sc_si_local++;

			if(si_local->bin_shift!=sc_si_local->bin_shift+first_shift) {
				match=0;
				break;
				}
			}
		if(match) {
			SIMPLE_CACHE.hits++;
			/*  align pps with stored data */
			pps->offset-=first_shift;
			accumulate_partial_power_sum_F(pps, SIMPLE_CACHE.pps[k]);
			pps->offset+=first_shift;
			//fprintf(stderr, "hit\n");
			return;
			}
		SIMPLE_CACHE.overwrites++;
		break;
		}
	}

SIMPLE_CACHE.misses++;
//fprintf(stderr, "miss\n");

if(k>=SIMPLE_CACHE.size) {
	fprintf(stderr, "*** INTERNAL ERROR: cache overflow\n");
	exit(-1);
	}

if(k>=SIMPLE_CACHE.free) {
	SIMPLE_CACHE.si[k]=do_alloc(SIMPLE_CACHE.segment_count, sizeof(*si));
	SIMPLE_CACHE.pps[k]=allocate_partial_power_sum_F(useful_bins+2*max_shift);
	SIMPLE_CACHE.free++;
	}

SIMPLE_CACHE.key[k]=key;
memcpy(SIMPLE_CACHE.si[k], si, SIMPLE_CACHE.segment_count*sizeof(*si));

get_uncached_single_bin_power_sum(SIMPLE_CACHE.si[k], SIMPLE_CACHE.segment_count, SIMPLE_CACHE.pps[k]);
accumulate_partial_power_sum_F(pps, SIMPLE_CACHE.pps[k]);
}

void print_cache_stats(void)
{
fprintf(stderr, "SIMPLE_CACHE stats: hits=%ld misses=%ld overwrites=%ld large_shifts=%ld miss ratio %f overwrite ratio %f\n", SIMPLE_CACHE.hits, SIMPLE_CACHE.misses, SIMPLE_CACHE.overwrites, SIMPLE_CACHE.large_shifts,
	SIMPLE_CACHE.misses/(SIMPLE_CACHE.hits+SIMPLE_CACHE.misses+0.0), SIMPLE_CACHE.overwrites/(SIMPLE_CACHE.hits+SIMPLE_CACHE.misses+0.0));

fprintf(LOG, "SIMPLE_CACHE stats: hits=%ld misses=%ld overwrites=%ld  large_shifts=%ld miss ratio %f overwrite ratio %f\n", SIMPLE_CACHE.hits, SIMPLE_CACHE.misses, SIMPLE_CACHE.overwrites, SIMPLE_CACHE.large_shifts,
	SIMPLE_CACHE.misses/(SIMPLE_CACHE.hits+SIMPLE_CACHE.misses+0.0), SIMPLE_CACHE.overwrites/(SIMPLE_CACHE.hits+SIMPLE_CACHE.misses+0.0));
}
