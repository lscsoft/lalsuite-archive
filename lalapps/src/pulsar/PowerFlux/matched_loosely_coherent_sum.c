#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/* We need this define to get NAN values */
#define __USE_ISOC99
#include <math.h>
#include <gsl/gsl_cblas.h>

#include <pmmintrin.h>
#include <xmmintrin.h>

#include <lal/DetectorSite.h>
#include <lal/LALBarycenter.h>
#include <lal/LALDetectors.h>

#include "polarization.h"
#include "dataset.h"
#include "util.h"
#include "power_cache.h"
#include "summing_context.h"
#include "power_sums.h"
#include "hookup.h"
#include "cmdline.h"
#include "single_bin_loosely_coherent_sum.h"

extern DATASET *datasets;
extern int d_free;
extern int nbins, first_bin, side_cut, useful_bins;
extern struct gengetopt_args_info args_info;

extern FILE *LOG;
extern INT64 spindown_start;

typedef struct {
	long signature;
	
	EmissionTime	*emission_time;
	long emission_time_size;

	/* template parameters */
	double freq_shift;
	double spindown;
	double inv_coherence_length;

	/* these entries hold cache of computed matched power sums */
	long computed_size;
	int pps_bins;
	float *m_re; /* real part */
	float *m_im; /* imaginary */
	float *w;  /* weight */
	unsigned char *computed; /* boolean flag */
	
	} MATCHED_LOOSELY_COHERENT_PATCH_PRIVATE_DATA;


/* Lanczos window actually vanishes */
#define LOOSE_SEARCH_TOLERANCE 0.0


/* 
	This function is passed a list of which is split in two parts, A and B, and
        it computes contribution of all terms of the form A_i*B_j
	if A_1!=B_1 it multiplies the result by 2.0

	The first ctx->loose_first_half_count segments got to part A

	This function is meant to work with accumulate_loose_power_sums_sidereal_step
*/


static double inline exp_kernel(double delta, double gps_delta)
{
return(exp(gps_delta*log(sin(delta)/delta)/1800.0));
}

static double inline sinc_kernel(double delta, double gps_delta)
{
double b;
if(gps_delta<=0)return 1.0;
b=delta*gps_delta/1800.0;
return(sin(b)/b);
}

static double inline lanczos_kernel3(double delta, double gps_delta)
{
double b;
if(gps_delta<=0)return 1.0;
b=delta*gps_delta/1800.0;
if(b>3.0*M_PI)return 0.0;
return(sinc_kernel(delta, gps_delta)*sinc_kernel(delta, gps_delta/3.0));
}

static double inline lanczos_kernel2(double delta, double gps_delta)
{
double b;
if(gps_delta<=0)return 1.0;
b=delta*gps_delta/1800.0;
if(b>2.0*M_PI)return 0.0;
return(sinc_kernel(delta, gps_delta)*sinc_kernel(delta, gps_delta/2.0));
}

#define loose_kernel lanczos_kernel3

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

float *filter;

__m128 v4power0, v4power1, v4a, v4b, v4a1, v4b1, v4filt0, v4filt1;
float *tmp1, *tmp2, *tmp3;

tmp1=aligned_alloca(8*sizeof(*tmp1));
tmp2=aligned_alloca(8*sizeof(*tmp2));
tmp3=aligned_alloca(4*sizeof(*tmp3));
filter=aligned_alloca(8*sizeof(*filter));

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


/* This helper function populates real and imaginary part of an FFT array, and also computes weight
   An additional wrinkle is that uses diff_bin_shift to adjust mismatch constant which can vary
   this only becomes an issue with 1/10 or finer subbin resolution */

void sse_compute_matched_floating_fft(SUMMING_CONTEXT *ctx, SEGMENT_INFO *si, int offset, int pps_bins, float *out_re, float *out_im, float *out_w)
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

float mismatch, mismatch_inc;

float *filter;

__m128 v4power0, v4power1, v4a, v4b, v4a1, v4b1, v4filt0, v4filt1;
float *tmp1, *tmp2, *tmp3;

tmp1=aligned_alloca(8*sizeof(*tmp1));
tmp2=aligned_alloca(8*sizeof(*tmp2));
tmp3=aligned_alloca(4*sizeof(*tmp3));
filter=aligned_alloca(8*sizeof(*filter));

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

/* Do this for every 8 frequency bins */
mismatch=si->bin_shift-rintf(si->bin_shift)-si->diff_bin_shift*(0.5*nbins-side_cut+4.0);
mismatch_inc=si->diff_bin_shift*8.0;

if(fabsf(mismatch_inc)>0.02)
	fprintf(stderr, "mismatch=%g mismatch_inc=%g bin_shift=%g diff_bin_shift=%g\n", mismatch, mismatch_inc, si->bin_shift, si->diff_bin_shift);

sum=0;
sum_sq=0;
tmp1[7]=0.0;
tmp2[7]=0.0;
for(i=0;i<pps_bins;i++) {
	if((i & 0x7)==0) {
		tabulated_fill_hann_filter7(filter, mismatch);
		filter[7]=0.0;
		v4filt0=_mm_load_ps(filter);
		v4filt1=_mm_load_ps(&(filter[4]));
		mismatch+=mismatch_inc;
		}

	memcpy(tmp1, &(re[-3]), 7*sizeof(*tmp1));
	memcpy(tmp2, &(im[-3]), 7*sizeof(*tmp2));

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

void get_uncached_loose_matched_partial_power_sum(SUMMING_CONTEXT *ctx, SEGMENT_INFO *si, int count, PARTIAL_POWER_SUM_F *pps)
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
double phase_offset, phase_increment, gps1, gps2, gps_delta, gps_mid;
float f0m_c, f0m_s, inc_c, inc_s;
int same_halfs=(si[0].segment==si[ctx->loose_first_half_count].segment) && (si[0].dataset==si[ctx->loose_first_half_count].dataset);

int pps_bins=pps->nbins;

float weight_pppp=0;
float weight_pppc=0;
float weight_ppcc=0;
float weight_pccc=0;
float weight_cccc=0;

MATCHED_LOOSELY_COHERENT_PATCH_PRIVATE_DATA *priv=(MATCHED_LOOSELY_COHERENT_PATCH_PRIVATE_DATA *)ctx->patch_private_data;


if(ctx->loose_first_half_count<0) {
	fprintf(stderr, "**** INTERNAL ERROR: loose_first_half_count=%d is not set.\n", ctx->loose_first_half_count);
	exit(-1);
	}

if(priv==NULL || priv->signature!=PATCH_PRIVATE_MATCHED_LOOSELY_COHERENT_SIGNATURE) {
	fprintf(stderr, "**** INTERNAL ERROR: invalid patch structure, priv=%p signature=%ld\n", priv, priv==NULL?-1 : priv->signature);
	exit(-1); 
	}
	
if(priv->pps_bins<pps_bins) {
	fprintf(stderr, "**** INTERNAL ERROR: not enough allocated memory, priv->pps_bins=%d while pps_bins=%d\n", priv->pps_bins, pps_bins);
	exit(-1);
	}

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


n=0;
for(k=0;k<ctx->loose_first_half_count;k++)
for(m=(same_halfs?k:0);m<(count-ctx->loose_first_half_count);m++) {
	si_local=&(si[k]);
	si_local2=&(si[m+ctx->loose_first_half_count]);

	/* off diagonal entries are x2 */
	if(same_halfs && (k==m))x=1.0;
		else x=2.0;

	x*=loose_kernel(args_info.phase_mismatch_arg, fabs(si_local->gps-si_local2->gps));

	if(fabs(x)<=LOOSE_SEARCH_TOLERANCE)continue;

	n++;

 	d=&(datasets[si_local->dataset]);
 	d2=&(datasets[si_local2->dataset]);


	if(!priv->computed[k]) {
		sse_compute_matched_floating_fft(ctx, si_local, pps->offset, pps_bins,
			&(priv->m_re[k*priv->pps_bins]),
			&(priv->m_im[k*priv->pps_bins]),
			&(priv->w[k]));
		priv->computed[k]=1; 
		}

	if(!priv->computed[m+ctx->loose_first_half_count]) {
		sse_compute_matched_floating_fft(ctx, si_local2, pps->offset, pps_bins,
			&(priv->m_re[(m+ctx->loose_first_half_count)*priv->pps_bins]),
			&(priv->m_im[(m+ctx->loose_first_half_count)*priv->pps_bins]),
			&(priv->w[m+ctx->loose_first_half_count]));
		priv->computed[m+ctx->loose_first_half_count]=1; 
		}

	//fprintf(stderr, "%d %d %f %f\n", k, m, si_local->gps-si_local2->gps, x);

	f_plus=si_local->f_plus;
	f_cross=si_local->f_cross;

	f_plus2=si_local2->f_plus;
	f_cross2=si_local2->f_cross;

	f_pp=f_plus*f_plus2;
	f_pc=0.5*(f_plus*f_cross2+f_plus2*f_cross);
	f_cc=f_cross*f_cross2;

	//fprintf(stderr, "pp=%f pc=%f cc=%f\n", f_pp, f_pc, f_cc);	

	weight=x*sqrt(priv->w[k]*priv->w[m+ctx->loose_first_half_count]);

	weight_pppp+=weight*f_pp*f_pp;
	weight_pppc+=weight*f_pp*f_pc;
	weight_ppcc+=weight*(0.6666667*f_pc*f_pc+0.3333333*f_pp*f_cc); /* 2/3 and 1/3 */
	weight_pccc+=weight*f_pc*f_cc;
	weight_cccc+=weight*f_cc*f_cc;

	f_pp*=weight;
	f_pc*=weight;
	f_cc*=weight;

#if 0
	/* The code below uses Doppler shifts alone (same as injection code in dataset.c), and does not
	correctly treat multiple detectors, let alone Shapiro delay, etc. */
	
	/* contribution from frequency mismatch */

	/* This effectively rounds off phase offset to units of pi/900, good enough ! */
	phase_offset=((int)rint((((first_bin+side_cut) % 1800))*(si_local->gps-si_local2->gps)) % 1800 )*2*M_PI/1800.0;

	phase_offset+=((int)rint(((0.5*(si_local->bin_shift+si_local2->bin_shift)-0.5*(0.5*nbins-side_cut)*(si_local->diff_bin_shift+si_local2->diff_bin_shift))*(si_local->gps-si_local2->gps))) %1800)*2*M_PI/1800.0;
	//phase_offset+=M_PI*(si_local->bin_shift-si_local2->bin_shift-rintf(si_local->bin_shift)+rintf(si_local2->bin_shift));
	phase_offset+=M_PI*(si_local->bin_shift-si_local2->bin_shift-rint(si_local->bin_shift)+rint(si_local2->bin_shift));

	phase_increment=(1.0+0.5*(si_local->diff_bin_shift+si_local2->diff_bin_shift))*(si_local->gps-si_local2->gps)*2*M_PI/1800.0+
			(si_local->diff_bin_shift-si_local2->diff_bin_shift)*M_PI;

//	fprintf(stderr, "(%f %f)  %.4f %.4f", si_local->bin_shift, si_local2->bin_shift, phase_offset, phase_increment);

#else
	/* This uses LALBarycenter() timing */
	
	/* First compute phase offsets in units of 2*M_PI */
	gps1=(priv->emission_time[si_local->index].te.gpsSeconds-spindown_start)+1e-9*priv->emission_time[si_local->index].te.gpsNanoSeconds;
	gps2=(priv->emission_time[si_local2->index].te.gpsSeconds-spindown_start)+1e-9*priv->emission_time[si_local2->index].te.gpsNanoSeconds;
	gps_delta=(priv->emission_time[si_local->index].te.gpsSeconds-priv->emission_time[si_local2->index].te.gpsSeconds)+1e-9*(priv->emission_time[si_local->index].te.gpsNanoSeconds-priv->emission_time[si_local2->index].te.gpsNanoSeconds);
	gps_mid=0.5*(gps1+gps2);

	phase_offset=((first_bin+side_cut)*priv->inv_coherence_length+priv->freq_shift+priv->spindown*gps_mid)*gps_delta;
	/* we get an extra M_PI in phase from jumping one bin
	 * This happens because SFT is computed from t=0 but our gps refers to middle of the interval
	 * Every other bin picks one pie of phase.
	 */
	phase_offset+=0.5*(rintf(si_local->bin_shift)-rintf(si_local2->bin_shift));

	phase_increment=gps_delta*priv->inv_coherence_length;

	//fprintf(stderr, "phase_offset=%f phase_increment=%f\n", phase_offset, phase_increment);

	phase_offset=2*M_PI*(phase_offset-floor(phase_offset));
	phase_increment=2*M_PI*(phase_increment-floor(phase_increment));

	//fprintf(stderr, "tdot=(%g,%g)\n", priv->emission_time[si_local->index].tDot, priv->emission_time[si_local2->index].tDot); 

	//fprintf(stderr, " [%d, %d] (%f %f) %.4f %.4f (%f %f %f %g)\n", k, m, gps_mid, gps_delta, phase_offset, phase_increment, priv->emission_time[si_local->index].deltaT, priv->emission_time[si_local2->index].deltaT, (priv->emission_time[si_local->index].deltaT*gps1-priv->emission_time[si_local2->index].deltaT*gps2), priv->spindown);
#endif

	f0m_c=cosf(phase_offset);
	f0m_s=sinf(-phase_offset);

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

	re=&(priv->m_re[k*priv->pps_bins]);
	im=&(priv->m_im[k*priv->pps_bins]);

	re2=&(priv->m_re[(m+ctx->loose_first_half_count)*priv->pps_bins]);
	im2=&(priv->m_im[(m+ctx->loose_first_half_count)*priv->pps_bins]);

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

	}

pps->c_weight_pppp=weight_pppp;
pps->c_weight_pppc=weight_pppc;
pps->c_weight_ppcc=weight_ppcc;
pps->c_weight_pccc=weight_pccc;
pps->c_weight_cccc=weight_cccc;

pps->collapsed_weight_arrays=0;
}


static int is_nonzero_loose_partial_power_sum(SUMMING_CONTEXT *ctx, SEGMENT_INFO *si1, int count1, SEGMENT_INFO *si2, int count2)
{
int k,m;
SEGMENT_INFO *si_local1, *si_local2;
//POLARIZATION *pl;
double x, beta;
int same_halfs=(si1[0].segment==si2[0].segment) && (si1[0].dataset==si2[0].dataset);
MATCHED_LOOSELY_COHERENT_PATCH_PRIVATE_DATA *priv=(MATCHED_LOOSELY_COHERENT_PATCH_PRIVATE_DATA *)ctx->patch_private_data;

if(ctx->loose_first_half_count<0) {
	fprintf(stderr, "**** INTERNAL ERROR: loose_first_half_count=%d is not set.\n", ctx->loose_first_half_count);
	exit(-1);
	}

if(priv==NULL || priv->signature!=PATCH_PRIVATE_MATCHED_LOOSELY_COHERENT_SIGNATURE) {
	fprintf(stderr, "**** INTERNAL ERROR: invalid patch structure, priv=%p signature=%ld\n", priv, priv==NULL?-1 : priv->signature);
	exit(-1); 
	}

if(same_halfs)return(1);


si_local1=si1;
for(k=0;k<count1;k++) {
	si_local2=si2;
	for(m=0;m<count2;m++) {

		beta=fabs(si_local1->gps-si_local2->gps)*args_info.phase_mismatch_arg*priv->inv_coherence_length;

		if(same_halfs && (k==m))x=1.0;
			else x=2.0;
	
		x*=loose_kernel(args_info.phase_mismatch_arg, fabs(si_local1->gps-si_local2->gps));
		
		if(fabs(x)<=LOOSE_SEARCH_TOLERANCE) {
			si_local2++;
			continue;
			}
		return(1);
		}
	si_local1++;
	}
return(0);
}


extern EphemerisData ephemeris;

/* This function is meant to work with get_uncached_loose_partial_power_sum */
void accumulate_matched_loose_power_sums_sidereal_step(SUMMING_CONTEXT *ctx, POWER_SUM *ps, int count, double gps_start, double gps_stop, int veto_mask)
{
int segment_count;
SEGMENT_INFO *si, *si_local;
POWER_SUM *ps_local;
DATASET *d;
POLARIZATION *pl;
int i, j, k, m;
float min_shift, max_shift, a;
double gps_idx, gps_idx_next;
double gps_step=ctx->summing_step;
double center_frequency=(first_bin+nbins*0.5);
int group_count=ctx->sidereal_group_count*ctx->time_group_count;
SEGMENT_INFO **groups;
SEGMENT_INFO *tmp;
int tmp_count;
int max_group_segment_count;
int *group_segment_count;
double avg_spindown=args_info.spindown_start_arg+0.5*args_info.spindown_step_arg*(args_info.spindown_count_arg-1);
MATCHED_LOOSELY_COHERENT_PATCH_PRIVATE_DATA *priv;
LALStatus status={level:0, statusPtr:NULL};
EarthState earth_state;
LIGOTimeGPS tGPS;
BarycenterInput baryinput;

float *patch_e=ps[0].patch_e; /* set of coefficients for this patch, used for amplitude response and bin shift estimation */

if(ctx->patch_private_data!=NULL) {
	fprintf(stderr, "*** INTERNAL ERROR: patch private data is not NULL, (%p)\n", ctx->patch_private_data);
	exit(-1);
	}

//fprintf(stderr, "%p %p %d %lf %lf 0x%08x\n", ctx, ps, count, gps_start, gps_stop, veto_mask);

for(gps_idx=gps_start; gps_idx<gps_stop; gps_idx+=gps_step) {

	gps_idx_next=gps_idx+gps_step;
	si=find_segments(gps_idx, (gps_idx_next<=gps_stop ? gps_idx_next : gps_stop), veto_mask, &segment_count);
	if(segment_count<1) {
		free(si);
		continue;
		}

	/* This assumes that we are patch bound !! *** WARNING ***
	TODO: make this assumption automatic in the data layout.
	 */
	si_local=si;
	min_shift=1000000; /* we should never load this many bins */
	max_shift=-1000000;
	for(j=0;j<segment_count;j++) {
		/* assign index */
		si[j].index=j;

// 		si_local->ra=ps[0].patch_ra;
// 		si_local->dec=ps[0].patch_dec;
// 		memcpy(si_local->e, ps[0].patch_e, GRID_E_COUNT*sizeof(SKY_GRID_TYPE));

		d=&(datasets[si_local->dataset]);
		pl=&(d->polarizations[0]);

// 		si_local->f_plus=F_plus_coeff(si_local->segment,  patch_e, pl->AM_coeffs);
// 		si_local->f_cross=F_plus_coeff(si_local->segment,  patch_e, pl->conjugate->AM_coeffs);

		si_local->f_plus=F_plus_coeff(si_local->segment,  patch_e, d->AM_coeffs_plus);
		si_local->f_cross=F_plus_coeff(si_local->segment,  patch_e, d->AM_coeffs_cross);


		a=center_frequency*args_info.doppler_multiplier_arg*(patch_e[0]*si_local->detector_velocity[0]
						+patch_e[1]*si_local->detector_velocity[1]
						+patch_e[2]*si_local->detector_velocity[2])
			+si_local->coherence_time*avg_spindown*(float)(si_local->gps-spindown_start);
		if(a<min_shift)min_shift=a;
		if(a>max_shift)max_shift=a;
		si_local++;
		}

	if(group_count>200) {
		fprintf(stderr, "Warning group count too large: %d\n", group_count);
		group_count=200;
		}

	group_segment_count=do_alloc(group_count, sizeof(*group_segment_count));
	groups=do_alloc(group_count, sizeof(*groups));

	for(k=0;k<group_count;k++) {
		group_segment_count[k]=0;
		groups[k]=do_alloc(segment_count, sizeof(SEGMENT_INFO));
		}

	/* group segments into bunches with similar shifts - mostly by sidereal time
           this way there is larger correllation of frequency shifts during summing and better use of power cache */
	si_local=si;
	for(j=0;j<segment_count;j++) {
		a=(center_frequency*args_info.doppler_multiplier_arg*(patch_e[0]*si_local->detector_velocity[0]
						+patch_e[1]*si_local->detector_velocity[1]
						+patch_e[2]*si_local->detector_velocity[2])
			+si_local->coherence_time*avg_spindown*(float)(si_local->gps-spindown_start));
		//a*=0.25;
		k=floorf((a-floorf(a))*ctx->sidereal_group_count)+ctx->sidereal_group_count*floorf((si_local->gps-gps_idx)*ctx->time_group_count/gps_step);
//		k=floorf((a-floorf(a))*ctx->sidereal_group_count);
		if(k<0)k=0;
		if(k>=group_count)k=group_count-1;

		memcpy(&(groups[k][group_segment_count[k]]), si_local, sizeof(SEGMENT_INFO));
		group_segment_count[k]++;
		
		si_local++;
		}

	max_group_segment_count=group_segment_count[0];
 	for(k=1;k<group_count;k++) {
		if(group_segment_count[k]>max_group_segment_count)max_group_segment_count=group_segment_count[k];
 		//fprintf(stderr, "group %d has %d segments\n", k, group_segment_count[k]);
 		}

	priv=do_alloc(1, sizeof(*priv));
	priv->signature=PATCH_PRIVATE_MATCHED_LOOSELY_COHERENT_SIGNATURE;
	priv->emission_time_size=count*segment_count;
	priv->emission_time=do_alloc(priv->emission_time_size, sizeof(*priv->emission_time));
	priv->inv_coherence_length=1.0/si->coherence_time;
	
	priv->computed_size=max_group_segment_count*2;
	priv->pps_bins=ps->pps->nbins+args_info.max_first_shift_arg*2;
	priv->m_re=do_alloc(priv->computed_size*priv->pps_bins, sizeof(*priv->m_re));
	priv->m_im=do_alloc(priv->computed_size*priv->pps_bins, sizeof(*priv->m_im));
	priv->w=do_alloc(priv->computed_size, sizeof(*priv->w));
	priv->computed=do_alloc(priv->computed_size, sizeof(*priv->computed));


	for(j=0;j<segment_count;j++) {
		tGPS.gpsSeconds=si[j].gps+si[j].coherence_time*0.5;
		tGPS.gpsNanoSeconds=0;
		LALBarycenterEarth(&status, &earth_state, &tGPS, &ephemeris);
		TESTSTATUS(&status);

		for(i=0;i<count;i++) {
			baryinput.tgps.gpsSeconds=si[j].gps+si[j].coherence_time*0.5;
			baryinput.tgps.gpsNanoSeconds=0;
			baryinput.site=get_detector_struct(datasets[si[j].dataset].detector);
			baryinput.site.location[0]=baryinput.site.location[0]/LAL_C_SI;
			baryinput.site.location[1]=baryinput.site.location[1]/LAL_C_SI;
			baryinput.site.location[2]=baryinput.site.location[2]/LAL_C_SI;
			baryinput.alpha=ps[i].ra;
			baryinput.delta=ps[i].dec;
			baryinput.dInv=args_info.dInv_arg;

			LALBarycenter(&status, &(priv->emission_time[i*segment_count+j]), &baryinput, &earth_state);
			TESTSTATUS(&status);
			}
		}
	
	ctx->patch_private_data=priv;

	tmp=do_alloc(max_group_segment_count*2, sizeof(*tmp));

	/* Do a double loop over groups */

	for(k=0;k<group_count;k++) {
		if(group_segment_count[k]<1)continue;
		ctx->loose_first_half_count=group_segment_count[k];
		memcpy(tmp, groups[k], group_segment_count[k]*sizeof(*groups[k]));

		for(m=k;m<group_count;m++) {
			//fprintf(stderr, "group %d has %d segments\n", k, group_segment_count[k]);
			if(group_segment_count[m]<1)continue;

			/* Loosely coherent code main coefficient depends on GPS time alone, so we can
			   make a quick test and skip all the templates */
			if(!is_nonzero_loose_partial_power_sum(ctx, groups[k], group_segment_count[k], groups[m], group_segment_count[m]))continue;

			memcpy(&(tmp[ctx->loose_first_half_count]), groups[m], group_segment_count[m]*sizeof(*groups[m]));
			tmp_count=ctx->loose_first_half_count+group_segment_count[m];

			/* reset partial sum cache */
			ctx->reset_cache(ctx, tmp_count, count);
		
			/* loop over templates */
			ps_local=ps;
			for(i=0;i<count;i++) {
				/* fill in segment info appropriate to this template */
				si_local=tmp;
				priv->freq_shift=ps_local->freq_shift;
				priv->spindown=ps_local->spindown;
				
				/* reset matched power cache */
				memset(priv->computed, 0, priv->computed_size*sizeof(*priv->computed));

				for(j=0;j<tmp_count;j++) {
					si_local->index= (j<ctx->loose_first_half_count ? groups[k][j].index : groups[m][j-ctx->loose_first_half_count].index)+segment_count*i;
		
					si_local->bin_shift=si_local->coherence_time*(ps_local->freq_shift+ps_local->spindown*(si_local->gps+si_local->coherence_time*0.5-spindown_start))+
						(center_frequency+ps_local->freq_shift)*args_info.doppler_multiplier_arg*(ps_local->e[0]*si_local->detector_velocity[0]
							+ps_local->e[1]*si_local->detector_velocity[1]
							+ps_local->e[2]*si_local->detector_velocity[2]);
					si_local->diff_bin_shift=args_info.doppler_multiplier_arg*(ps_local->e[0]*si_local->detector_velocity[0]
							+ps_local->e[1]*si_local->detector_velocity[1]
							+ps_local->e[2]*si_local->detector_velocity[2]);
					si_local++;
					}
				ctx->accumulate_power_sum_cached(ctx, tmp, tmp_count, ps_local->pps);
				ps_local++;
				}
			}
		}
	//print_simple_cache_stats(ctx);
	free(tmp);
	tmp=NULL;
	for(k=0;k<group_count;k++) {
		free(groups[k]);
		}
	free(groups);
	free(group_segment_count);
	free(si);
	free(priv->emission_time);
	free(priv->m_re);
	free(priv->m_im);
	free(priv->w);
	free(priv->computed);
	memset(priv, 0, sizeof(*priv));
	free(priv);
	ctx->patch_private_data=NULL;
	}

for(i=0;i<count;i++) {
	if(ps[i].min_gps<0 || ps[i].min_gps>gps_start)ps[i].min_gps=gps_start;
	if(ps[i].max_gps<0 || ps[i].max_gps<gps_stop)ps[i].max_gps=gps_stop;
	}
}
