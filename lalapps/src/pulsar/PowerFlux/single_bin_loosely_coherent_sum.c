#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/* We need this define to get NAN values */
#define __USE_ISOC99
#include <math.h>
#include <gsl/gsl_cblas.h>

#include <pmmintrin.h>
#include <xmmintrin.h>

#include "polarization.h"
#include "dataset.h"
#include "util.h"
#include "power_cache.h"
#include "summing_context.h"
#include "cmdline.h"

extern DATASET *datasets;
extern int d_free;
extern int nbins, first_bin, side_cut, useful_bins;
extern struct gengetopt_args_info args_info;

extern FILE *LOG;


#define LOOSE_SEARCH_TOLERANCE 0.3

/* This computes exp(a) for a<=0.22 with 4% precision */
float fast_negexp(float a)
{
if(a>0.22) {
	fprintf(stderr, "INTERNAL ERROR: fast_negexp is not correct for argument values >0.22");
	exit(-1);
	}
if(a< -3.7)return(0.0);
return(1.0/(1.0-a+a*a));
}


/* 
	This function is passed a list of which is split in two parts, A and B, and
        it computes contribution of all terms of the form A_i*B_j
	if A_1!=B_1 it multiplies the result by 2.0

	The first ctx->loose_first_half_count segments got to part A

	This function is meant to work with accumulate_loose_power_sums_sidereal_step
*/

void get_uncached_loose_single_bin_partial_power_sum(SUMMING_CONTEXT *ctx, SEGMENT_INFO *si, int count, PARTIAL_POWER_SUM_F *pps)
{
int i,k,n,m;
int bin_shift, bin_shift2;
SEGMENT_INFO *si_local, *si_local2;
DATASET *d, *d2;
//POLARIZATION *pl;
float *im, *re, *im2, *re2, *pp, *pc, *cc;
float a;
float weight;
float f_plus, f_cross, f_plus2, f_cross2, f_pp, f_pc, f_cc;
float x, y;
float alpha;
double phase_offset, phase_increment;
float f0m_c, f0m_s, inc_c, inc_s;
int same_halfs=(si[0].segment==si[ctx->loose_first_half_count].segment) && (si[0].dataset==si[ctx->loose_first_half_count].dataset);

int pps_bins=pps->nbins;

float weight_pppp=0;
float weight_pppc=0;
float weight_ppcc=0;
float weight_pccc=0;
float weight_cccc=0;

float threshold;

if(ctx->loose_first_half_count<0) {
	fprintf(stderr, "**** INTERNAL ERROR: loose_first_half_count=%d is not set.\n", ctx->loose_first_half_count);
	exit(-1);
	}

alpha= ctx->loose_coherence_alpha;

/*pp=aligned_alloca(pps_bins*sizeof(*pp));
pc=aligned_alloca(pps_bins*sizeof(*pc));
cc=aligned_alloca(pps_bins*sizeof(*cc));
power=aligned_alloca(pps_bins*sizeof(*power)); */

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

threshold=-logf(LOOSE_SEARCH_TOLERANCE/2.0)/ctx->loose_coherence_alpha;

//fprintf(stderr, "%d\n", count);

n=0;
for(k=0;k<ctx->loose_first_half_count;k++)
for(m=(same_halfs?k:0);m<(count-ctx->loose_first_half_count);m++) {
	si_local=&(si[k]);
	si_local2=&(si[m+ctx->loose_first_half_count]);

	/* off diagonal entries are x2 */
	if(same_halfs && (k==m))x=1.0;
		else x=2.0;

// 	y=0*sinf((si_local->bin_shift-si_local2->bin_shift)*M_PI*0.5);
// 	if(fabs(y)<=0.001)y=1.0;
// 		else y=sinf(y)/y;
// 	y=-log(fabs(y))/1800.0;

	if(fabs(si_local->gps-si_local2->gps)>threshold)continue;

	x*=fast_negexp(-(alpha)*fabs(si_local->gps-si_local2->gps));

	if(x<LOOSE_SEARCH_TOLERANCE)continue;

	n++;

 	d=&(datasets[si_local->dataset]);
 	d2=&(datasets[si_local2->dataset]);


	//fprintf(stderr, "%d %d %f %f\n", k, m, si_local->gps-si_local2->gps, x);

	f_plus=si_local->f_plus;
	f_cross=si_local->f_cross;

	f_plus2=si_local2->f_plus;
	f_cross2=si_local2->f_cross;

	f_pp=f_plus*f_plus2;
	f_pc=0.5*(f_plus*f_cross2+f_plus2*f_cross);
	f_cc=f_cross*f_cross2;
	
	weight=x*d->expTMedians[si_local->segment]*d->weight*d2->expTMedians[si_local2->segment]*d2->weight;

	weight_pppp+=weight*f_pp*f_pp;
	weight_pppc+=weight*f_pp*f_pc;
	weight_ppcc+=weight*(0.6666667*f_pc*f_pc+0.3333333*f_pp*f_cc); /* 2/3 and 1/3 */
	weight_pccc+=weight*f_pc*f_cc;
	weight_cccc+=weight*f_cc*f_cc;

	f_pp*=weight;
	f_pc*=weight;
	f_cc*=weight;

	/* contribution from frequency mismatch */
// 	phase_offset=0.5*(si_local->bin_shift+si_local2->bin_shift)*(si_local->gps-si_local2->gps)*2*M_PI/1800.0;
// 	phase_offset+=M_PI*(si_local->bin_shift-si_local2->bin_shift-rintf(si_local->bin_shift)+rintf(si_local2->bin_shift));

	phase_offset=(((first_bin+side_cut) % 1800)+0.5*(si_local->bin_shift+si_local2->bin_shift)-0.5*(0.5*nbins-side_cut)*(si_local->diff_bin_shift+si_local2->diff_bin_shift))*(si_local->gps-si_local2->gps)*2*M_PI/1800.0;
	//phase_offset+=M_PI*(si_local->bin_shift-si_local2->bin_shift-rintf(si_local->bin_shift)+rintf(si_local2->bin_shift));
	phase_offset+=M_PI*(si_local->bin_shift-si_local2->bin_shift-rintf(si_local->bin_shift)+rintf(si_local2->bin_shift));

	f0m_c=cosf(phase_offset);
	f0m_s=sinf(-phase_offset);

	phase_increment=(1.0+0.5*(si_local->diff_bin_shift+si_local2->diff_bin_shift))*(si_local->gps-si_local2->gps)*2*M_PI/1800.0+
			(si_local->diff_bin_shift-si_local2->diff_bin_shift)*M_PI;

	inc_c=cosf(phase_increment);
	inc_s=sinf(-phase_increment);

 	d=&(datasets[si_local->dataset]);
 	d2=&(datasets[si_local2->dataset]);

	bin_shift=rintf(si_local->bin_shift)-pps->offset;
	if((bin_shift+side_cut<0) || (bin_shift>pps_bins+side_cut)) {
		fprintf(stderr, "*** Attempt to sample outside loaded range bin_shift=%d bin_shift=%lg, aborting\n", 
			bin_shift, si_local->bin_shift);
		exit(-1);
		}

	bin_shift2=rintf(si_local2->bin_shift)-pps->offset;
	if((bin_shift2+side_cut<0) || (bin_shift2>pps_bins+side_cut)) {
		fprintf(stderr, "*** Attempt to sample outside loaded range bin_shift=%d bin_shift=%lg, aborting\n", 
			bin_shift2, si_local2->bin_shift);
		exit(-1);
		}

	re=&(d->re[si_local->segment*nbins+side_cut+bin_shift]);
	im=&(d->im[si_local->segment*nbins+side_cut+bin_shift]);

	re2=&(d2->re[si_local2->segment*nbins+side_cut+bin_shift2]);
	im2=&(d2->im[si_local2->segment*nbins+side_cut+bin_shift2]);

	pp=pps->power_pp;
	pc=pps->power_pc;
	cc=pps->power_cc;

	for(i=0;i<pps_bins;i++) {
		x=(*re)*f0m_c-(*im)*f0m_s;
		y=(*re)*f0m_s+(*im)*f0m_c;

		a=(x*(*re2)+y*(*im2));

		x=f0m_c*inc_c-f0m_s*inc_s;
		y=f0m_c*inc_s+f0m_s*inc_c;

		f0m_c=x;
		f0m_s=y;

		(*pp)+=a*f_pp;
		(*pc)+=a*f_pc;
		(*cc)+=a*f_cc;

		pp++;
		pc++;
		cc++;

		re++;
		im++;

		re2++;
		im2++;
		/**/
		}

	pp=pps->power_pp;
	pc=pps->power_pc;
	cc=pps->power_cc;

	/*

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

	*/
	}

//fprintf(stderr, "n=%d out of %d, %f fraction\n", n, count*(count+1)/2, n/(count*(count+1)*0.5));

pps->c_weight_pppp=weight_pppp;
pps->c_weight_pppc=weight_pppc;
pps->c_weight_ppcc=weight_ppcc;
pps->c_weight_pccc=weight_pccc;
pps->c_weight_cccc=weight_cccc;

pps->collapsed_weight_arrays=0;
}


int is_nonzero_loose_partial_power_sum(SUMMING_CONTEXT *ctx, SEGMENT_INFO *si1, int count1, SEGMENT_INFO *si2, int count2)
{
int k,m;
SEGMENT_INFO *si_local1, *si_local2;
//POLARIZATION *pl;
float x;
float alpha;
int same_halfs=(si1[0].segment==si2[0].segment) && (si1[0].dataset==si2[0].dataset);
float threshold;

if(ctx->loose_first_half_count<0) {
	fprintf(stderr, "**** INTERNAL ERROR: loose_first_half_count=%d is not set.\n", ctx->loose_first_half_count);
	exit(-1);
	}

if(same_halfs)return(1);

threshold=-logf(LOOSE_SEARCH_TOLERANCE/2.0)/ctx->loose_coherence_alpha;

si_local1=si1;
for(k=0;k<count1;k++) {
	si_local2=si2;
	for(m=0;m<count2;m++) {
		
		if(fabs(si_local1->gps-si_local2->gps)>threshold) {
			si_local2++;
			continue;
			}
		return(1);
		}
	si_local1++;
	}
return(0);
}

#define KEY_MULT 8388623
#define KEY_DIV ((1<<31)-1)

/* Note: si is modified in place to have bin_shift rounded as appropriate to the power generating algorithm */
void accumulate_power_sum_cached_diff(SUMMING_CONTEXT *ctx, SEGMENT_INFO *si, int count, PARTIAL_POWER_SUM_F *pps)
{
int i, k;
int match;
int key;
float a;
int first_shift;
int max_shift=args_info.max_first_shift_arg;
SEGMENT_INFO *si_local, *sc_si_local;
int *sc_key;

SIMPLE_CACHE *sc=(SIMPLE_CACHE *)ctx->cache;
if((sc==NULL) || (sc->id!=SIMPLE_CACHE_ID)) {
	fprintf(stderr, "INTERNAL ERROR: simple cache id does not match in %s, aborting\n", __FUNCTION__);
	exit(-1);
	}

if(count!=sc->segment_count) { 
	fprintf(stderr, "Internal error: segment counts don't match.\n");
	exit(-1);
	}

/* normalize frequency shifts by rounding off to nearest cache granularity value */
si_local=si;
key=0.0;
k=rintf(si_local->bin_shift)*ctx->cache_granularity;
for(i=0;i<count;i++) {
	//fprintf(stderr, "%0.1f ", si_local->bin_shift);
	a=rintf(si_local->bin_shift*ctx->cache_granularity);

	key=( (key+(int)a-k) * KEY_MULT) & KEY_DIV;

	si_local->bin_shift=a*ctx->inv_cache_granularity;
	si_local++;
	//fprintf(stderr, "%0.1f ", a);
	}
//fprintf(stderr, "k=%d key=%d %f\n", k, key, a);
sc_key=sc->key;
for(k=0;k<sc->free;k++) {
	/* the reason we use exact equality for floating point numbers is because the numbers in cache have been placed there by the same function. */
	if(key==sc_key[k]) {
		/* we found the box holding our data, double check it is the right one */
		si_local=si;
		sc_si_local=sc->si[k];
		first_shift=rintf(si_local->bin_shift-sc_si_local->bin_shift);
		if( (first_shift>max_shift) || (first_shift< -max_shift)) {
			sc->large_shifts++;
			break;
			}
		match=1;
		for(i=0;i<count;i++) {
			if(fabs(si_local->bin_shift-(sc_si_local->bin_shift+first_shift))> ctx->half_inv_cache_granularity) {
				match=0;
				fprintf(stderr, "OVERWRITE: i=%d key=%d count=%d %f %f %d\n", i, key, count, si_local->bin_shift, sc_si_local->bin_shift, first_shift);
				break;
				}

			si_local++;
			sc_si_local++;
			}
		if(match) {
			sc->hits++;
			/*  align pps with stored data */
			pps->offset-=first_shift;
			sse_accumulate_partial_power_sum_F(pps, sc->pps[k]);
			pps->offset+=first_shift;
			//fprintf(stderr, "hit\n");
			return;
			}
		sc->overwrites++;
		break;
		}
	}

sc->misses++;
//fprintf(stderr, "miss\n");

if(k>=sc->size) {
	fprintf(stderr, "*** INTERNAL ERROR: cache overflow\n");
	exit(-1);
	}

if(k>=sc->free) {
	sc->si[k]=do_alloc(sc->segment_count, sizeof(*si));
	sc->pps[k]=allocate_partial_power_sum_F(useful_bins+2*max_shift);
	sc->free++;
	}

sc->key[k]=key;
memcpy(sc->si[k], si, sc->segment_count*sizeof(*si));

ctx->get_uncached_power_sum(ctx, sc->si[k], sc->segment_count, sc->pps[k]);
sse_accumulate_partial_power_sum_F(pps, sc->pps[k]);
}
