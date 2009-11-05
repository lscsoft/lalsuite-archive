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

static void inline sum_F(int count, float *partial, float *accum)
{
int i;
for(i=0;i<count;i++)accum[i]+=partial[i];
}

static void inline sse_sum_unaligned1_F(int count, float *partial, float *accum)
{
int i;
__m128 v4in, v4acc;
float *tmp=aligned_alloca(4*sizeof(*tmp));

for(i=0; ((long)&(accum[i]) & 0x0f) && (i<count);i++) {
	accum[i]+=partial[i];
	}

//fprintf(stderr, "unaligned %d %p %p\n", i, partial, accum);

for(;i<(count-3);i+=4) {
	memcpy(tmp, &(partial[i]), 4*sizeof(*tmp));
	v4in=_mm_load_ps(tmp);
	v4acc=_mm_load_ps(&(accum[i]));
	v4acc=_mm_add_ps(v4acc, v4in);
	_mm_store_ps(&(accum[i]), v4acc);
	}
for(;i<count;i++) {
	accum[i]+=partial[i];
	}
}

static void inline sse_sum_unaligned_F(int count, float *partial, float *accum)
{
int i;
__m128 v4in, v4acc;
float *tmp=aligned_alloca(16*sizeof(*tmp));

for(i=0; ((long)&(accum[i]) & 0x0f) && (i<count);i++) {
	accum[i]+=partial[i];
	}

//fprintf(stderr, "unaligned %d %p %p\n", i, partial, accum);

for(;i<(count-15);i+=16) {
	memcpy(tmp, &(partial[i]), 16*sizeof(*tmp));

	v4in=_mm_load_ps(tmp);
	v4acc=_mm_load_ps(&(accum[i]));
	v4acc=_mm_add_ps(v4acc, v4in);
	_mm_store_ps(&(accum[i]), v4acc);

	v4in=_mm_load_ps(&(tmp[4]));
	v4acc=_mm_load_ps(&(accum[i+4]));
	v4acc=_mm_add_ps(v4acc, v4in);
	_mm_store_ps(&(accum[i+4]), v4acc);

	v4in=_mm_load_ps(&(tmp[8]));
	v4acc=_mm_load_ps(&(accum[i+8]));
	v4acc=_mm_add_ps(v4acc, v4in);
	_mm_store_ps(&(accum[i+8]), v4acc);

	v4in=_mm_load_ps(&(tmp[12]));
	v4acc=_mm_load_ps(&(accum[i+12]));
	v4acc=_mm_add_ps(v4acc, v4in);
	_mm_store_ps(&(accum[i+12]), v4acc);
	}
for(;i<count;i++) {
	accum[i]+=partial[i];
	}
}

static void inline sse_sum_aligned_F(int count, float *partial, float *accum)
{
int i;
__m128 v4in, v4acc;

for(i=0; ((long)&(accum[i]) & 0x0f) && (i<count);i++) {
	accum[i]+=partial[i];
	}

//fprintf(stderr, "aligned %d %p %p\n", i, partial, accum);

for(;i<(count-3);i+=4) {
	v4in=_mm_load_ps(&(partial[i]));
	v4acc=_mm_load_ps(&(accum[i]));
	v4acc=_mm_add_ps(v4acc, v4in);
	_mm_store_ps(&(accum[i]), v4acc);
	}
for(;i<count;i++) {
	accum[i]+=partial[i];
	}
}

static void inline sse_sum_F(int count, float *partial, float *accum)
{
//sum_F(count, partial, accum);

//return;
if( ((long)accum & 0x0f) == ((long)partial & 0x0f)) sse_sum_aligned_F(count, partial, accum);
	else sse_sum_unaligned_F(count, partial, accum);
//	else sum_F(count, partial, accum);
//fprintf(stderr, "done\n");
}

#define REAL double
#define SUFFIX(a) a
#define CBLAS_AXPY cblas_daxpy
#define SSE_SUM(...)
#include "partial_power_sum.c"

#undef SSE_SUM
#undef CBLAS_AXPY
#undef REAL
#undef SUFFIX

#define REAL float
#define SUFFIX(a) a##_F
#define CBLAS_AXPY cblas_saxpy
#define SSE_SUM sse_sum_F
#include "partial_power_sum.c"

#undef SSE_SUM
#undef CBLAS_AXPY
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

void get_uncached_single_bin_power_sum(SUMMING_CONTEXT *ctx, SEGMENT_INFO *si, int count, PARTIAL_POWER_SUM_F *pps)
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

void sse_get_uncached_single_bin_power_sum(SUMMING_CONTEXT *ctx, SEGMENT_INFO *si, int count, PARTIAL_POWER_SUM_F *pps)
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

float sum, sum_sq;
float pps_bins_inv=1.0/pps_bins;
//float pmax_factor=args_info.power_max_median_factor_arg;

float weight_pppp=0;
float weight_pppc=0;
float weight_ppcc=0;
float weight_pccc=0;
float weight_cccc=0;
__m128 v4power, v4tm, v4sum, v4sum_sq, v4pp, v4pc, v4cc, v4a, v4b;
float *tmp1, *tmp2, *tmp3;

/*pp=aligned_alloca(pps_bins*sizeof(*pp));
pc=aligned_alloca(pps_bins*sizeof(*pc));
cc=aligned_alloca(pps_bins*sizeof(*cc));*/
power=aligned_alloca(pps_bins*sizeof(*power));

tmp1=aligned_alloca(16*sizeof(*tmp1));
tmp2=aligned_alloca(16*sizeof(*tmp2));
tmp3=aligned_alloca(16*sizeof(*tmp3));

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
		TM=-d->expTMedians[si_local->segment];
		} else {
		TM=0;
		}

	p=power;

	v4tm=_mm_load1_ps(&TM);
	v4sum=_mm_setzero_ps();
	v4sum_sq=_mm_setzero_ps();
	for(i=0;i<(pps_bins-15);i+=16) {
		memcpy(tmp1, re, 16*sizeof(*tmp1));
		memcpy(tmp2, im, 16*sizeof(*tmp2));
		memcpy(tmp3, fm, 16*sizeof(*tmp3));

		#define S1(tmp1, tmp2, tmp3, p) \
			v4a=_mm_load_ps(tmp1); \
			v4power=_mm_mul_ps(v4a, v4a); \
				\
			v4b=_mm_load_ps(tmp2); \
			v4b=_mm_mul_ps(v4b, v4b); \
			v4power=_mm_add_ps(v4power, v4b); \
				\
			v4a=_mm_load_ps(tmp3); \
			v4a=_mm_mul_ps(v4tm, v4a); \
			v4power=_mm_add_ps(v4power, v4a); \
				\
			_mm_store_ps(p, v4power); \
				\
			v4sum=_mm_add_ps(v4sum, v4power); \
			v4a=_mm_mul_ps(v4power, v4power); \
			v4sum_sq=_mm_add_ps(v4sum_sq, v4a); 

		S1(tmp1, tmp2, tmp3, p)
		S1(&(tmp1[4]), &(tmp2[4]), &(tmp3[4]), &(p[4]))
		S1(&(tmp1[8]), &(tmp2[8]), &(tmp3[8]), &(p[8]))
		S1(&(tmp1[12]), &(tmp2[12]), &(tmp3[12]), &(p[12]))

		im+=16;
		re+=16;
		fm+=16;
		p+=16;
		}

	_mm_store_ps(tmp1, v4sum);
	sum=tmp1[0]+tmp1[1]+tmp1[2]+tmp1[3];

	_mm_store_ps(tmp2, v4sum_sq);
	sum_sq=tmp2[0]+tmp2[1]+tmp2[2]+tmp2[3];

	for(;i<pps_bins;i++) {
		a=((*re)*(*re)+(*im)*(*im)+TM*(*fm));
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

	v4pp=_mm_load1_ps(&f_pp);
	v4pc=_mm_load1_ps(&f_pc);
	v4cc=_mm_load1_ps(&f_cc);

	for(i=0;i<(pps_bins-3);i+=4) {
		v4power=_mm_load_ps(p);

		_mm_store_ps(pp, _mm_add_ps(_mm_load_ps(pp), _mm_mul_ps(v4power, v4pp)));
		_mm_store_ps(pc, _mm_add_ps(_mm_load_ps(pc), _mm_mul_ps(v4power, v4pc)));
		_mm_store_ps(cc, _mm_add_ps(_mm_load_ps(cc), _mm_mul_ps(v4power, v4cc)));

		p+=4;
		pp+=4;
		pc+=4;
		cc+=4;
		/**/
		}

	for(;i<pps_bins;i++) {
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

void get_uncached_matched_power_sum(SUMMING_CONTEXT *ctx, SEGMENT_INFO *si, int count, PARTIAL_POWER_SUM_F *pps)
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
float filter[7];
float x,y;

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

	tabulated_fill_hann_filter7(filter, (si_local->bin_shift-rintf(si_local->bin_shift)));

	pmax=0.0;
	sum=0.0;
	sum_sq=0.0;
	p=power;
	for(i=0;i<pps_bins;i++) {
		x=re[-3]*filter[0]+re[-2]*filter[1]+re[-1]*filter[2]+re[0]*filter[3]+re[1]*filter[4]+re[2]*filter[5]+re[3]*filter[6];
		y=im[-3]*filter[0]+im[-2]*filter[1]+im[-1]*filter[2]+im[0]*filter[3]+im[1]*filter[4]+im[2]*filter[5]+im[3]*filter[6];
	
		a=x*x+y*y;

		//a=((*re)*(*re)+(*im)*(*im)-TM*(*fm));
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

pps->c_weight_pppp=weight_pppp;
pps->c_weight_pppc=weight_pppc;
pps->c_weight_ppcc=weight_ppcc;
pps->c_weight_pccc=weight_pccc;
pps->c_weight_cccc=weight_cccc;

pps->collapsed_weight_arrays=0;
}

void sse_get_uncached_matched_power_sum(SUMMING_CONTEXT *ctx, SEGMENT_INFO *si, int count, PARTIAL_POWER_SUM_F *pps)
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

float sum, sum_sq;
float pps_bins_inv=1.0/pps_bins;
//float pmax_factor=args_info.power_max_median_factor_arg;

float filter[8];

float weight_pppp=0;
float weight_pppc=0;
float weight_ppcc=0;
float weight_pccc=0;
float weight_cccc=0;
__m128 v4power, v4power0, v4power1, v4tm, v4sum, v4sum_sq, v4pp, v4pc, v4cc, v4a, v4b, v4a1, v4b1, v4filt0, v4filt1;
float *tmp1, *tmp2;

/*pp=aligned_alloca(pps_bins*sizeof(*pp));
pc=aligned_alloca(pps_bins*sizeof(*pc));
cc=aligned_alloca(pps_bins*sizeof(*cc));*/
power=aligned_alloca(pps_bins*sizeof(*power));

tmp1=aligned_alloca(8*sizeof(*tmp1));
tmp2=aligned_alloca(8*sizeof(*tmp2));

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
		TM=-d->expTMedians[si_local->segment];
		} else {
		TM=0;
		}

	tabulated_fill_hann_filter7(filter, (si_local->bin_shift-rintf(si_local->bin_shift)));
	filter[7]=0.0;
	v4filt0=_mm_load_ps(filter);
	v4filt1=_mm_load_ps(&(filter[4]));

	p=power;

	sum=0;
	sum_sq=0;
	tmp1[7]=0;
	tmp2[7]=0;
	for(i=0;i<pps_bins;i++) {
		memcpy(tmp1, &(re[-3]), 8*sizeof(*tmp1));
		memcpy(tmp2, &(im[-3]), 8*sizeof(*tmp2));

		v4a=_mm_load_ps(&(tmp1[0]));
		v4a1=_mm_load_ps(&(tmp1[4]));

		v4a=_mm_mul_ps(v4filt0, v4a);
		v4a1=_mm_mul_ps(v4filt1, v4a1);

		v4a=_mm_add_ps(v4a, v4a1);

		v4b=_mm_load_ps(&(tmp2[0]));
		v4b1=_mm_load_ps(&(tmp2[4]));

		v4b=_mm_mul_ps(v4filt0, v4b);
		v4b1=_mm_mul_ps(v4filt1, v4b1);

		v4b=_mm_add_ps(v4b, v4b1);

//		v4a=_mm_add_ps(_mm_mul_ps(v4filt0, _mm_loadu_ps(&(re[-3]))), _mm_mul_ps(v4filt1, _mm_load_ps(&(re[1]))));
//		v4b=_mm_add_ps(_mm_mul_ps(v4filt0, _mm_loadu_ps(&(im[-3]))), _mm_mul_ps(v4filt1, _mm_load_ps(&(im[1]))));
	
		v4power0=_mm_hadd_ps(v4a, v4b);
		v4power0=_mm_hadd_ps(v4power0, v4power0);
		v4power1=_mm_mul_ps(v4power0, v4power0);
		v4power1=_mm_hadd_ps(v4power1, v4power1);

		_mm_store_ss(&a, v4power1);

		*p=a;

		sum+=a;
		sum_sq+=a*a;

		im++;
		re++;
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

	v4pp=_mm_load1_ps(&f_pp);
	v4pc=_mm_load1_ps(&f_pc);
	v4cc=_mm_load1_ps(&f_cc);

	for(i=0;i<(pps_bins-3);i+=4) {
		v4power=_mm_load_ps(p);

		_mm_store_ps(pp, _mm_add_ps(_mm_load_ps(pp), _mm_mul_ps(v4power, v4pp)));
		_mm_store_ps(pc, _mm_add_ps(_mm_load_ps(pc), _mm_mul_ps(v4power, v4pc)));
		_mm_store_ps(cc, _mm_add_ps(_mm_load_ps(cc), _mm_mul_ps(v4power, v4cc)));

		p+=4;
		pp+=4;
		pc+=4;
		cc+=4;
		/**/
		}

	for(;i<pps_bins;i++) {
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

	/*
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
	*/
	}

pps->c_weight_pppp=weight_pppp;
pps->c_weight_pppc=weight_pppc;
pps->c_weight_ppcc=weight_ppcc;
pps->c_weight_pccc=weight_pccc;
pps->c_weight_cccc=weight_cccc;

pps->collapsed_weight_arrays=0;
}

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

#define LOOSE_SEARCH_TOLERANCE 0.3

void get_uncached_loose_power_sum(SUMMING_CONTEXT *ctx, SEGMENT_INFO *si, int count, PARTIAL_POWER_SUM_F *pps)
{
int i,k,n,m;
int bin_shift;
SEGMENT_INFO *si_local, *si_local2;
DATASET *d;
//POLARIZATION *pl;
float *im, *re, *pp, *pc, *cc, *fm;
float *power, *p;
float a;
float TM;
float weight;
float f_plus, f_cross, f_plus2, f_cross2, f_pp, f_pc, f_cc;
float filter[7];
float x, y;
float *m_re, *m_im, *w;
float alpha;
float f0_mismatch;
double phase_offset;
float f0m_c, f0m_s;
double first_gps=si->gps;

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

m_re=do_alloc(count*pps_bins, sizeof(*m_re));
m_im=do_alloc(count*pps_bins, sizeof(*m_im));
w=do_alloc(count, sizeof(*w));

alpha=ctx->loose_coherence_alpha;

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



	re=&(d->re[si_local->segment*nbins+side_cut+bin_shift]);
	im=&(d->im[si_local->segment*nbins+side_cut+bin_shift]);
	fm=&(d->expFMedians_plain[side_cut+bin_shift]);

	if(args_info.subtract_background_arg) {
		TM=d->expTMedians[si_local->segment];
		} else {
		TM=0;
		}

	f0_mismatch=(si_local->bin_shift-rintf(si_local->bin_shift));
	tabulated_fill_hann_filter7(filter, f0_mismatch);

	/* we use the fact that si_local->gps is always integer */
	//phase_offset=f0_mismatch*M_PI+f0_mismatch*(si_local->gps-first_gps)*2*M_PI;
	phase_offset=f0_mismatch*M_PI;

	f0m_c=cos(phase_offset);
	f0m_s=sin(-phase_offset);

	pmax=0.0;
	sum=0.0;
	sum_sq=0.0;
	p=power;
	for(i=0;i<pps_bins;i++) {
		x=re[-3]*filter[0]+re[-2]*filter[1]+re[-1]*filter[2]+re[0]*filter[3]+re[1]*filter[4]+re[2]*filter[5]+re[3]*filter[6];
		y=im[-3]*filter[0]+im[-2]*filter[1]+im[-1]*filter[2]+im[0]*filter[3]+im[1]*filter[4]+im[2]*filter[5]+im[3]*filter[6];

		m_re[i+k*pps_bins]=x*f0m_c-y*f0m_s;	
		m_im[i+k*pps_bins]=x*f0m_s+y*f0m_c;	

		a=x*x+y*y;

		//a=((*re)*(*re)+(*im)*(*im)-TM*(*fm));
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

	w[k]=sqrt(weight);
	}

//fprintf(stderr, "alpha=%.8g\n", alpha);

n=0;
for(k=0;k<count;k++)
for(m=k;m<count;m++) {
	si_local=&(si[k]);
	si_local2=&(si[m]);

	/* off diagonal entries are x2 */
	if(k==m)x=1.0;
		else x=2.0;

	y=0*sin((si_local->bin_shift-si_local2->bin_shift)*M_PI*0.5);
	if(fabs(y)<=0.001)y=1.0;
		else y=sin(y)/y;
	y=-log(fabs(y))/1800.0;

	x*=exp(-(alpha+y)*fabs(si_local->gps-si_local2->gps));

	if(x<LOOSE_SEARCH_TOLERANCE)continue;

	n++;

	//fprintf(stderr, "%d %d %f %f\n", k, m, si_local->gps-si_local2->gps, x);

	f_plus=si_local->f_plus;
	f_cross=si_local->f_cross;

	f_plus2=si_local2->f_plus;
	f_cross2=si_local2->f_cross;

	f_pp=f_plus*f_plus2;
	f_pc=0.5*(f_plus*f_cross2+f_plus2*f_cross);
	f_cc=f_cross*f_cross2;
	
	weight=x*w[k]*w[m];

	weight_pppp+=weight*f_pp*f_pp;
	weight_pppc+=weight*f_pp*f_pc;
	weight_ppcc+=weight*f_pp*f_cc;
	weight_pccc+=weight*f_pc*f_cc;
	weight_cccc+=weight*f_cc*f_cc;

	f_pp*=weight;
	f_pc*=weight;
	f_cc*=weight;

	pp=pps->power_pp;
	pc=pps->power_pc;
	cc=pps->power_cc;

	phase_offset=0.5*(si_local->bin_shift+si_local2->bin_shift)*(si_local->gps-si_local2->gps)*2*M_PI/1800.0;

	f0m_c=cos(phase_offset);
	f0m_s=sin(-phase_offset);

	for(i=0;i<pps_bins;i++) {
		x=m_re[i+k*pps_bins]*f0m_c-m_im[i+k*pps_bins]*f0m_s;
		y=m_re[i+k*pps_bins]*f0m_s+m_im[i+k*pps_bins]*f0m_c;

		a=(x*m_re[i+m*pps_bins]+y*m_im[i+m*pps_bins]);


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

free(m_re);
free(m_im);
free(w);

pps->c_weight_pppp=weight_pppp;
pps->c_weight_pppc=weight_pppc;
pps->c_weight_ppcc=weight_ppcc;
pps->c_weight_pccc=weight_pccc;
pps->c_weight_cccc=weight_cccc;

pps->collapsed_weight_arrays=0;
}

void get_uncached_single_bin_loose_power_sum(SUMMING_CONTEXT *ctx, SEGMENT_INFO *si, int count, PARTIAL_POWER_SUM_F *pps)
{
int i,k,n,m;
int bin_shift;
SEGMENT_INFO *si_local, *si_local2;
DATASET *d;
//POLARIZATION *pl;
float *im, *re, *pp, *pc, *cc, *fm;
float *power, *p;
float a;
float TM;
float weight;
float f_plus, f_cross, f_plus2, f_cross2, f_pp, f_pc, f_cc;
float filter[7];
float x, y;
float *m_re, *m_im, *w;
float alpha;
float f0_mismatch;
double phase_offset;
float f0m_c, f0m_s;

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

m_re=do_alloc(count*pps_bins, sizeof(*m_re));
m_im=do_alloc(count*pps_bins, sizeof(*m_im));
w=do_alloc(count, sizeof(*w));

alpha=ctx->loose_coherence_alpha;

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



	re=&(d->re[si_local->segment*nbins+side_cut+bin_shift]);
	im=&(d->im[si_local->segment*nbins+side_cut+bin_shift]);
	fm=&(d->expFMedians_plain[side_cut+bin_shift]);

	if(args_info.subtract_background_arg) {
		TM=d->expTMedians[si_local->segment];
		} else {
		TM=0;
		}

	f0_mismatch=(si_local->bin_shift-rintf(si_local->bin_shift));
	tabulated_fill_hann_filter7(filter, f0_mismatch);

	/* we use the fact that si_local->gps is always integer */
	//phase_offset=f0_mismatch*M_PI+f0_mismatch*(si_local->gps-first_gps)*2*M_PI;
	phase_offset=f0_mismatch*M_PI;

	f0m_c=cos(phase_offset);
	f0m_s=sin(-phase_offset);

	pmax=0.0;
	sum=0.0;
	sum_sq=0.0;
	p=power;
	for(i=0;i<pps_bins;i++) {
// 		x=re[-3]*filter[0]+re[-2]*filter[1]+re[-1]*filter[2]+re[0]*filter[3]+re[1]*filter[4]+re[2]*filter[5]+re[3]*filter[6];
// 		y=im[-3]*filter[0]+im[-2]*filter[1]+im[-1]*filter[2]+im[0]*filter[3]+im[1]*filter[4]+im[2]*filter[5]+im[3]*filter[6];

		x=re[0];
		y=im[0];

		m_re[i+k*pps_bins]=x*f0m_c-y*f0m_s;	
		m_im[i+k*pps_bins]=x*f0m_s+y*f0m_c;	

		a=x*x+y*y;

		//a=((*re)*(*re)+(*im)*(*im)-TM*(*fm));
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

	w[k]=sqrt(weight);
	}

for(k=0;k<count;k++)
for(m=k;m<count;m++) {
	si_local=&(si[k]);
	si_local2=&(si[m]);

	/* off diagonal entries are x2 */
	if(k==m)x=1.0;
		else x=2.0;

	y=sin((si_local->bin_shift-si_local2->bin_shift)*M_PI);
	if(fabs(y)<=0.001)y=1.0;
		else y=sin(y)/y;

	x*=exp((-alpha+0*log(y)/1800.0)*fabs(si_local->gps-si_local2->gps));

	if(x<LOOSE_SEARCH_TOLERANCE)continue;

	//fprintf(stderr, "%d %d %f %f\n", k, m, si_local->gps-si_local2->gps, x);

	f_plus=si_local->f_plus;
	f_cross=si_local->f_cross;

	f_plus2=si_local2->f_plus;
	f_cross2=si_local2->f_cross;

	f_pp=f_plus*f_plus2;
	f_pc=0.5*(f_plus*f_cross2+f_plus2*f_cross);
	f_cc=f_cross*f_cross2;
	
	weight=x*w[k]*w[m];

	weight_pppp+=weight*f_pp*f_pp;
	weight_pppc+=weight*f_pp*f_pc;
	weight_ppcc+=weight*f_pp*f_cc;
	weight_pccc+=weight*f_pc*f_cc;
	weight_cccc+=weight*f_cc*f_cc;

	f_pp*=weight;
	f_pc*=weight;
	f_cc*=weight;

	pp=pps->power_pp;
	pc=pps->power_pc;
	cc=pps->power_cc;

	phase_offset=0.5*(si_local->bin_shift+si_local2->bin_shift)*(si_local->gps-si_local2->gps)*2*M_PI/1800.0;

	f0m_c=cos(phase_offset);
	f0m_s=sin(-phase_offset);

	for(i=0;i<pps_bins;i++) {
		x=m_re[i+k*pps_bins]*f0m_c-m_im[i+k*pps_bins]*f0m_s;
		y=m_re[i+k*pps_bins]*f0m_s+m_im[i+k*pps_bins]*f0m_c;

		a=(x*m_re[i+m*pps_bins]+y*m_im[i+m*pps_bins]);


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

free(m_re);
free(m_im);
free(w);

pps->c_weight_pppp=weight_pppp;
pps->c_weight_pppc=weight_pppc;
pps->c_weight_ppcc=weight_ppcc;
pps->c_weight_pccc=weight_pccc;
pps->c_weight_cccc=weight_cccc;

pps->collapsed_weight_arrays=0;
}

/* 
	This function is passed a list of which is split in two parts, A and B, and
        it computes contribution of all terms of the form A_i*B_j
	if A_1!=B_1 it multiplies the result by 2.0

	The first ctx->loose_first_half_count segments got to part A

	This function is meant to work with accumulate_loose_power_sums_sidereal_step
*/
void get_uncached_loose_partial_power_sum(SUMMING_CONTEXT *ctx, SEGMENT_INFO *si, int count, PARTIAL_POWER_SUM_F *pps)
{
int i,k,n,m;
int bin_shift;
SEGMENT_INFO *si_local, *si_local2;
DATASET *d;
//POLARIZATION *pl;
float *im, *re, *im2, *re2, *pp, *pc, *cc, *fm;
float *power, *p;
float a;
float TM;
float weight;
float f_plus, f_cross, f_plus2, f_cross2, f_pp, f_pc, f_cc;
float filter[7];
float x, y;
float *m_re, *m_im, *w;
float alpha;
float f0_mismatch;
double phase_offset;
float f0m_c, f0m_s;
int same_halfs=(si[0].segment==si[ctx->loose_first_half_count].segment) && (si[0].dataset==si[ctx->loose_first_half_count].dataset);

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

if(ctx->loose_first_half_count<0) {
	fprintf(stderr, "**** INTERNAL ERROR: loose_first_half_count=%d is not set.\n", ctx->loose_first_half_count);
	exit(-1);
	}

m_re=do_alloc(count*pps_bins, sizeof(*m_re));
m_im=do_alloc(count*pps_bins, sizeof(*m_im));
w=do_alloc(count, sizeof(*w));

alpha=ctx->loose_coherence_alpha;

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



	re=&(d->re[si_local->segment*nbins+side_cut+bin_shift]);
	im=&(d->im[si_local->segment*nbins+side_cut+bin_shift]);
	fm=&(d->expFMedians_plain[side_cut+bin_shift]);

	if(args_info.subtract_background_arg) {
		TM=d->expTMedians[si_local->segment];
		} else {
		TM=0;
		}

	f0_mismatch=(si_local->bin_shift-rintf(si_local->bin_shift));
	tabulated_fill_hann_filter7(filter, f0_mismatch);

	/* we use the fact that si_local->gps is always integer */
	//phase_offset=f0_mismatch*M_PI+f0_mismatch*(si_local->gps-first_gps)*2*M_PI;
	phase_offset=f0_mismatch*M_PI;

	f0m_c=cosf(phase_offset);
	f0m_s=sinf(-phase_offset);

	pmax=0.0;
	sum=0.0;
	sum_sq=0.0;
	p=power;
	for(i=0;i<pps_bins;i++) {
		x=re[-3]*filter[0]+re[-2]*filter[1]+re[-1]*filter[2]+re[0]*filter[3]+re[1]*filter[4]+re[2]*filter[5]+re[3]*filter[6];
		y=im[-3]*filter[0]+im[-2]*filter[1]+im[-1]*filter[2]+im[0]*filter[3]+im[1]*filter[4]+im[2]*filter[5]+im[3]*filter[6];

		m_re[i+k*pps_bins]=x*f0m_c-y*f0m_s;	
		m_im[i+k*pps_bins]=x*f0m_s+y*f0m_c;	

		a=x*x+y*y;

		//a=((*re)*(*re)+(*im)*(*im)-TM*(*fm));
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

	w[k]=sqrt(weight);
	}

//fprintf(stderr, "alpha=%.8g\n", alpha);

n=0;
for(k=0;k<ctx->loose_first_half_count;k++)
for(m=(same_halfs?k:0);m<(count-ctx->loose_first_half_count);m++) {
	si_local=&(si[k]);
	si_local2=&(si[m+ctx->loose_first_half_count]);

	/* off diagonal entries are x2 */
	if(same_halfs && (k==m))x=1.0;
		else x=2.0;

	y=0*sinf((si_local->bin_shift-si_local2->bin_shift)*M_PI*0.5);
	if(fabs(y)<=0.001)y=1.0;
		else y=sinf(y)/y;
	y=-log(fabs(y))/1800.0;

	x*=expf(-(alpha+y)*fabs(si_local->gps-si_local2->gps));

	if(x<LOOSE_SEARCH_TOLERANCE)continue;

	n++;

	//fprintf(stderr, "%d %d %f %f\n", k, m, si_local->gps-si_local2->gps, x);

	f_plus=si_local->f_plus;
	f_cross=si_local->f_cross;

	f_plus2=si_local2->f_plus;
	f_cross2=si_local2->f_cross;

	f_pp=f_plus*f_plus2;
	f_pc=0.5*(f_plus*f_cross2+f_plus2*f_cross);
	f_cc=f_cross*f_cross2;
	
	weight=x*w[k]*w[m];

	weight_pppp+=weight*f_pp*f_pp;
	weight_pppc+=weight*f_pp*f_pc;
	weight_ppcc+=weight*f_pp*f_cc;
	weight_pccc+=weight*f_pc*f_cc;
	weight_cccc+=weight*f_cc*f_cc;

	f_pp*=weight;
	f_pc*=weight;
	f_cc*=weight;

	pp=pps->power_pp;
	pc=pps->power_pc;
	cc=pps->power_cc;

	phase_offset=0.5*(si_local->bin_shift+si_local2->bin_shift)*(si_local->gps-si_local2->gps)*2*M_PI/1800.0;

	f0m_c=cosf(phase_offset);
	f0m_s=sinf(-phase_offset);

	re=&(m_re[i+k*pps_bins]);
	im=&(m_im[i+k*pps_bins]);

	re2=&(m_re[i+(m+ctx->loose_first_half_count)*pps_bins]);
	im2=&(m_im[i+(m+ctx->loose_first_half_count)*pps_bins]);

	for(i=0;i<pps_bins;i++) {
		x=(*re)*f0m_c-(*im)*f0m_s;
		y=(*re)*f0m_s+(*im)*f0m_c;

		a=(x*(*re2)+y*(*im2));


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

free(m_re);
free(m_im);
free(w);

pps->c_weight_pppp=weight_pppp;
pps->c_weight_pppc=weight_pppc;
pps->c_weight_ppcc=weight_ppcc;
pps->c_weight_pccc=weight_pccc;
pps->c_weight_cccc=weight_cccc;

pps->collapsed_weight_arrays=0;
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
double phase_offset;
float f0m_c, f0m_s;
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
	weight_ppcc+=weight*f_pp*f_cc;
	weight_pccc+=weight*f_pc*f_cc;
	weight_cccc+=weight*f_cc*f_cc;

	f_pp*=weight;
	f_pc*=weight;
	f_cc*=weight;

	/* contribution from frequency mismatch */
	phase_offset=0.5*(si_local->bin_shift+si_local2->bin_shift)*(si_local->gps-si_local2->gps)*2*M_PI/1800.0;
	phase_offset+=M_PI*(si_local->bin_shift-si_local2->bin_shift-rintf(si_local->bin_shift)+rintf(si_local2->bin_shift));

	f0m_c=cosf(phase_offset);
	f0m_s=sinf(-phase_offset);

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

/* This helper function populates real and imaginary part of an FFT array, and also computes weight */
void sse_compute_matched_fft(SUMMING_CONTEXT *ctx, SEGMENT_INFO *si, int offset, int pps_bins, float *out_re, float *out_im, float *out_w)
{
int i;
int bin_shift;
DATASET *d;
//POLARIZATION *pl;
float *im, *re, *out_im_p, *out_re_p;
float a;

float sum, sum_sq;
float pps_bins_inv=1.0/pps_bins;
//float pmax_factor=args_info.power_max_median_factor_arg;

float filter[8];

__m128 v4power0, v4power1, v4a, v4b, v4a1, v4b1, v4filt0, v4filt1;
float *tmp1, *tmp2, *tmp3;

tmp1=aligned_alloca(8*sizeof(*tmp1));
tmp2=aligned_alloca(8*sizeof(*tmp2));
tmp3=aligned_alloca(4*sizeof(*tmp3));

//fprintf(stderr, "%d\n", count);

bin_shift=rintf(si->bin_shift)-offset;
if((bin_shift+side_cut<0) || (bin_shift>pps_bins+side_cut)) {
	fprintf(stderr, "*** Attempt to sample outside loaded range bin_shift=%d bin_shift=%lg, aborting\n", 
		bin_shift, si->bin_shift);
	exit(-1);
	}

d=&(datasets[si->dataset]);
// 	pl=&(d->polarizations[0]);

//f_plus=F_plus_coeff(si_local->segment, si_local->e, pl->AM_coeffs);
//f_cross=F_plus_coeff(si_local->segment, si_local->e, pl->conjugate->AM_coeffs);

re=&(d->re[si->segment*nbins+side_cut+bin_shift]);
im=&(d->im[si->segment*nbins+side_cut+bin_shift]);
out_re_p=out_re;
out_im_p=out_im;

tabulated_fill_hann_filter7(filter, (si->bin_shift-rintf(si->bin_shift)));
filter[7]=0.0;
v4filt0=_mm_load_ps(filter);
v4filt1=_mm_load_ps(&(filter[4]));

sum=0;
sum_sq=0;
tmp1[7]=0;
tmp2[7]=0;
for(i=0;i<pps_bins;i++) {
	memcpy(tmp1, &(re[-3]), 8*sizeof(*tmp1));
	memcpy(tmp2, &(im[-3]), 8*sizeof(*tmp2));

	v4a=_mm_load_ps(&(tmp1[0]));
	v4a1=_mm_load_ps(&(tmp1[4]));

	v4a=_mm_mul_ps(v4filt0, v4a);
	v4a1=_mm_mul_ps(v4filt1, v4a1);

	v4a=_mm_add_ps(v4a, v4a1);

	v4b=_mm_load_ps(&(tmp2[0]));
	v4b1=_mm_load_ps(&(tmp2[4]));

	v4b=_mm_mul_ps(v4filt0, v4b);
	v4b1=_mm_mul_ps(v4filt1, v4b1);

	v4b=_mm_add_ps(v4b, v4b1);

//		v4a=_mm_add_ps(_mm_mul_ps(v4filt0, _mm_loadu_ps(&(re[-3]))), _mm_mul_ps(v4filt1, _mm_load_ps(&(re[1]))));
//		v4b=_mm_add_ps(_mm_mul_ps(v4filt0, _mm_loadu_ps(&(im[-3]))), _mm_mul_ps(v4filt1, _mm_load_ps(&(im[1]))));

	v4power0=_mm_hadd_ps(v4a, v4b);
	v4power0=_mm_hadd_ps(v4power0, v4power0);
	_mm_store_ps(tmp3, v4power0);
	*out_re_p=tmp3[0];
	*out_im_p=tmp3[1];

	v4power1=_mm_mul_ps(v4power0, v4power0);
	v4power1=_mm_hadd_ps(v4power1, v4power1);

	_mm_store_ss(&a, v4power1);

	sum+=a;
	sum_sq+=a*a;

	im++;
	re++;
	out_re_p++;
	out_im_p++;
	}

if(args_info.tmedian_noise_level_arg) {
	*out_w=d->expTMedians[si->segment]*d->weight;
	} else {
	sum*=pps_bins_inv;
	sum_sq*=pps_bins_inv;
	sum_sq-=sum*sum;

	*out_w=1.0/sum_sq;
	}

}

/* 
	This function is passed a list of which is split in two parts, A and B, and
        it computes contribution of all terms of the form A_i*B_j
	if A_1!=B_1 it multiplies the result by 2.0

	The first ctx->loose_first_half_count segments got to part A

	This function is meant to work with accumulate_loose_power_sums_sidereal_step
*/

void get_uncached_loose_matched_partial_power_sum(SUMMING_CONTEXT *ctx, SEGMENT_INFO *si, int count, PARTIAL_POWER_SUM_F *pps)
{
int i,k,n,m;
int bin_shift, bin_shift2;
SEGMENT_INFO *si_local, *si_local2;
/*DATASET *d, *d2;*/
//POLARIZATION *pl;
float *im, *re, *im2, *re2, *pp, *pc, *cc;
float *m_re, *m_im, *w;
char *computed;
float a;
float weight;
float f_plus, f_cross, f_plus2, f_cross2, f_pp, f_pc, f_cc;
float x, y;
float alpha;
double phase_offset;
float f0m_c, f0m_s;
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

m_re=do_alloc(count*pps_bins, sizeof(*m_re));
m_im=do_alloc(count*pps_bins, sizeof(*m_im));
w=do_alloc(count, sizeof(*w));
computed=do_alloc(count, sizeof(*computed));
memset(computed, 0, count*sizeof(*computed));

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

//  	d=&(datasets[si_local->dataset]);
//  	d2=&(datasets[si_local2->dataset]);

	if(!computed[k]) {
		sse_compute_matched_fft(ctx, si_local, pps->offset, pps_bins,
			&(m_re[k*pps_bins]),
			&(m_im[k*pps_bins]),
			&(w[k]));
		computed[k]=1; 
		}

	if(!computed[m+ctx->loose_first_half_count]) {
		sse_compute_matched_fft(ctx, si_local2, pps->offset, pps_bins,
			&(m_re[(m+ctx->loose_first_half_count)*pps_bins]),
			&(m_im[(m+ctx->loose_first_half_count)*pps_bins]),
			&(w[m+ctx->loose_first_half_count]));
		computed[m+ctx->loose_first_half_count]=1; 
		}

	//fprintf(stderr, "%d %d %f %f\n", k, m, si_local->gps-si_local2->gps, x);

	f_plus=si_local->f_plus;
	f_cross=si_local->f_cross;

	f_plus2=si_local2->f_plus;
	f_cross2=si_local2->f_cross;

	f_pp=f_plus*f_plus2;
	f_pc=0.5*(f_plus*f_cross2+f_plus2*f_cross);
	f_cc=f_cross*f_cross2;
	
	weight=x*w[k]*w[m+ctx->loose_first_half_count];

	weight_pppp+=weight*f_pp*f_pp;
	weight_pppc+=weight*f_pp*f_pc;
	weight_ppcc+=weight*f_pp*f_cc;
	weight_pccc+=weight*f_pc*f_cc;
	weight_cccc+=weight*f_cc*f_cc;

	f_pp*=weight;
	f_pc*=weight;
	f_cc*=weight;

	/* contribution from frequency mismatch */
	phase_offset=0.5*(si_local->bin_shift+si_local2->bin_shift)*(si_local->gps-si_local2->gps)*2*M_PI/1800.0;
	phase_offset+=M_PI*(si_local->bin_shift-si_local2->bin_shift-rintf(si_local->bin_shift)+rintf(si_local2->bin_shift));

	f0m_c=cosf(phase_offset);
	f0m_s=sinf(-phase_offset);

//  	d=&(datasets[si_local->dataset]);
//  	d2=&(datasets[si_local2->dataset]);
// 
// 	bin_shift=rintf(si_local->bin_shift)-pps->offset;
// 	if((bin_shift+side_cut<0) || (bin_shift>pps_bins+side_cut)) {
// 		fprintf(stderr, "*** Attempt to sample outside loaded range bin_shift=%d bin_shift=%lg, aborting\n", 
// 			bin_shift, si_local->bin_shift);
// 		exit(-1);
// 		}
// 
// 	bin_shift2=rintf(si_local2->bin_shift)-pps->offset;
// 	if((bin_shift2+side_cut<0) || (bin_shift2>pps_bins+side_cut)) {
// 		fprintf(stderr, "*** Attempt to sample outside loaded range bin_shift=%d bin_shift=%lg, aborting\n", 
// 			bin_shift2, si_local2->bin_shift);
// 		exit(-1);
// 		}

	re=&(m_re[k*pps_bins]);
	im=&(m_im[k*pps_bins]);

	re2=&(m_re[(m+ctx->loose_first_half_count)*pps_bins]);
	im2=&(m_im[(m+ctx->loose_first_half_count)*pps_bins]);

	pp=pps->power_pp;
	pc=pps->power_pc;
	cc=pps->power_cc;

	for(i=0;i<pps_bins;i++) {
		x=(*re)*f0m_c-(*im)*f0m_s;
		y=(*re)*f0m_s+(*im)*f0m_c;

		a=(x*(*re2)+y*(*im2));


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

free(m_re);
free(m_im);
free(w);
free(computed);

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

#define SIMPLE_CACHE_ID 1

typedef struct {
	long id;

	/* statistics */
	long hits;
	long misses;
	long overwrites;
	long large_shifts;
	int max_size;

	/* cache contents */
	int segment_count;
	int size;
	int free;
	int *key;
	SEGMENT_INFO **si;
	PARTIAL_POWER_SUM_F **pps;
	} SIMPLE_CACHE;

void free_simple_cache(SUMMING_CONTEXT *ctx)
{
int i;
SIMPLE_CACHE *sc=(SIMPLE_CACHE *)ctx->cache;
if(sc==NULL)return;
if(sc->id!=SIMPLE_CACHE_ID) {
	fprintf(stderr, "INTERNAL ERROR: simple cache id does not match in %s, aborting\n", __FUNCTION__);
	exit(-1);
	}
ctx->cache=NULL;
for(i=0;i<sc->free;i++) {
	free_partial_power_sum_F(sc->pps[i]);
	free(sc->si[i]);
	}
free(sc->pps);
free(sc->si);
free(sc->key);
free(sc);
}

void print_simple_cache_stats(SUMMING_CONTEXT *ctx)
{
SIMPLE_CACHE *sc=(SIMPLE_CACHE *)ctx->cache;
if(sc==NULL) {
	fprintf(stderr, "SIMPLE CACHE has not been initialized yet\n");
	fprintf(LOG, "SIMPLE CACHE has not been initialized yet\n");
	return;
	}
if(sc->id!=SIMPLE_CACHE_ID) {
	fprintf(stderr, "INTERNAL ERROR: simple cache id does not match in %s, aborting\n", __FUNCTION__);
	exit(-1);
	}
fprintf(stderr, "SIMPLE_CACHE stats: max_size=%d hits=%ld misses=%ld overwrites=%ld large_shifts=%ld miss ratio %f overwrite ratio %f\n", sc->max_size, sc->hits, sc->misses, sc->overwrites, sc->large_shifts,
	sc->misses/(sc->hits+sc->misses+0.0), sc->overwrites/(sc->hits+sc->misses+0.0));

fprintf(LOG, "SIMPLE_CACHE stats: max_size=%d hits=%ld misses=%ld overwrites=%ld  large_shifts=%ld miss ratio %f overwrite ratio %f\n", sc->max_size, sc->hits, sc->misses, sc->overwrites, sc->large_shifts,
	sc->misses/(sc->hits+sc->misses+0.0), sc->overwrites/(sc->hits+sc->misses+0.0));
}

void reset_simple_cache(SUMMING_CONTEXT *ctx, int segment_count, int template_count)
{
int i;
SIMPLE_CACHE *sc=(SIMPLE_CACHE *)ctx->cache;
if((sc==NULL) || (sc->id!=SIMPLE_CACHE_ID)) {
	fprintf(stderr, "INTERNAL ERROR: simple cache id does not match in %s, aborting\n", __FUNCTION__);
	exit(-1);
	}
if(sc->free>sc->max_size)sc->max_size=sc->free;
for(i=0;i<sc->free;i++) {
	free_partial_power_sum_F(sc->pps[i]);
	free(sc->si[i]);
	}
free(sc->pps);
free(sc->si);
free(sc->key);

sc->segment_count=segment_count;
sc->size=template_count;
sc->free=0;
sc->key=do_alloc(sc->size, sizeof(*sc->key));
sc->pps=do_alloc(sc->size, sizeof(*sc->pps));
sc->si=do_alloc(sc->size, sizeof(*sc->si));
}

void allocate_simple_cache(SUMMING_CONTEXT *ctx)
{
SIMPLE_CACHE *sc;
if(ctx->cache!=NULL)ctx->free_cache(ctx);
ctx->cache=NULL;
sc=do_alloc(1, sizeof(SIMPLE_CACHE));

memset(sc, 0, sizeof(*sc));
sc->id=SIMPLE_CACHE_ID;

ctx->free_cache=free_simple_cache;
ctx->print_cache_stats=print_simple_cache_stats;
ctx->reset_cache=reset_simple_cache;
ctx->cache=sc;
}


#define KEY_MULT 8388623
#define KEY_DIV ((1<<31)-1)

/* Note: si is modified in place to have bin_shift rounded as appropriate to the power generating algorithm */
void accumulate_power_sum_cached1(SUMMING_CONTEXT *ctx, SEGMENT_INFO *si, int count, PARTIAL_POWER_SUM_F *pps)
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

void power_cache_selftest(void)
{
SEGMENT_INFO *si=NULL;
PARTIAL_POWER_SUM_F *ps1=NULL, *ps2=NULL, *ps3=NULL, *ps4=NULL;
SUMMING_CONTEXT *ctx;
int i;
int count;
int result=0;

ctx=create_summing_context();

ps1=allocate_partial_power_sum_F(useful_bins+10);
ps2=allocate_partial_power_sum_F(useful_bins+10);
ps3=allocate_partial_power_sum_F(useful_bins+10);
ps4=allocate_partial_power_sum_F(useful_bins+10);

randomize_partial_power_sum_F(ps1);
randomize_partial_power_sum_F(ps2);

/* reference implementation */
zero_partial_power_sum_F(ps3);
accumulate_partial_power_sum_F(ps3, ps1);
accumulate_partial_power_sum_F(ps3, ps2);

/* cblas */
zero_partial_power_sum_F(ps4);
cblas_accumulate_partial_power_sum_F(ps4, ps1);
cblas_accumulate_partial_power_sum_F(ps4, ps2);
result+=compare_partial_power_sums_F("cblas test:", ps3, ps4);

/* sse */
zero_partial_power_sum_F(ps4);
sse_accumulate_partial_power_sum_F(ps4, ps1);
sse_accumulate_partial_power_sum_F(ps4, ps2);
result+=compare_partial_power_sums_F("sse test:", ps3, ps4);

si=find_segments(min_gps(), max_gps()+1, ~0, &count);
if(count>100)count=100;
for(i=0;i<count;i++) {
	si[i].bin_shift= (i%7)-3.1;
	si[i].f_plus= (i%11)/11.0;
	si[i].f_cross= (i%17)/17.0;
	}
/* reference implementation */
get_uncached_single_bin_power_sum(ctx, si, count, ps3);

/* sse implementation */
sse_get_uncached_single_bin_power_sum(ctx, si, count, ps4);
result+=compare_partial_power_sums_F("sse_get_uncached_single_bin_power_sum:", ps3, ps4);

/* reference implementation */
get_uncached_matched_power_sum(ctx, si, count, ps3);

/* sse implementation */
sse_get_uncached_matched_power_sum(ctx, si, count, ps4);
result+=compare_partial_power_sums_F("sse_get_uncached_matched_power_sum:", ps3, ps4);

free(si);
free_partial_power_sum_F(ps1);
free_partial_power_sum_F(ps2);
free_partial_power_sum_F(ps3);
free_partial_power_sum_F(ps4);

free_summing_context(ctx);

if(result<0) {
	fprintf(stderr, "*** POWER CACHE selftest failed, exiting\n");
	exit(-1);
	}

fprintf(stderr, "Power cache selftest: passed\n");
fprintf(LOG, "Power cache selftest: passed\n");
}
