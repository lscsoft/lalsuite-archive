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
		sum_sq=a*a;

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
