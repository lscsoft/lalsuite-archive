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
#include "cmdline.h"

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

	double freq_shift;
	double spindown;
	double inv_coherence_length;

	} SINGLE_BIN_LOOSELY_COHERENT_PATCH_PRIVATE_DATA;


/* Lanczos window actually vanishes */
#define LOOSE_SEARCH_TOLERANCE 0.0

/* This computes exp(a) for a<=0.22 with 4% precision */
static float fast_negexp(float a)
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
double phase_offset, phase_increment, gps1, gps2, gps_delta, gps_mid;
float f0m_c, f0m_s, inc_c, inc_s;
int same_halfs=(si[0].segment==si[ctx->loose_first_half_count].segment) && (si[0].dataset==si[ctx->loose_first_half_count].dataset);

int pps_bins=pps->nbins;

float weight_pppp=0;
float weight_pppc=0;
float weight_ppcc=0;
float weight_pccc=0;
float weight_cccc=0;

SINGLE_BIN_LOOSELY_COHERENT_PATCH_PRIVATE_DATA *priv=(SINGLE_BIN_LOOSELY_COHERENT_PATCH_PRIVATE_DATA *)ctx->patch_private_data;


if(ctx->loose_first_half_count<0) {
	fprintf(stderr, "**** INTERNAL ERROR: loose_first_half_count=%d is not set.\n", ctx->loose_first_half_count);
	exit(-1);
	}

if(priv==NULL || priv->signature!=PATCH_PRIVATE_SINGLE_BIN_LOOSELY_COHERENT_SIGNATURE) {
	fprintf(stderr, "**** INTERNAL ERROR: invalid patch structure, priv=%p signature=%ld\n", priv, priv==NULL?-1 : priv->signature);
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

//fprintf(stderr, "%d\n", count);

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


	//fprintf(stderr, "%d %d %f %f\n", k, m, si_local->gps-si_local2->gps, x);

	f_plus=si_local->f_plus;
	f_cross=si_local->f_cross;

	f_plus2=si_local2->f_plus;
	f_cross2=si_local2->f_cross;

	f_pp=f_plus*f_plus2;
	f_pc=0.5*(f_plus*f_cross2+f_plus2*f_cross);
	f_cc=f_cross*f_cross2;

	//fprintf(stderr, "pp=%f pc=%f cc=%f\n", f_pp, f_pc, f_cc);	

	#if 0 
	/* The weights below (commented and not) are used in completely ad-hoc manner, a proper way to go is to solve the filtering problem in the presence of non-stationary noise */
	
	weight=d->expTMedians[si_local->segment]*d->weight;
	if(weight> d2->expTMedians[si_local2->segment]*d2->weight)weight=d2->expTMedians[si_local2->segment]*d2->weight;
	weight*=x;
	#endif 
	
	weight=x*sqrt(d->expTMedians[si_local->segment]*d->weight*d2->expTMedians[si_local2->segment]*d2->weight);
	
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
	if((bin_shift+side_cut<0) || (bin_shift+pps_bins>nbins)) {
		fprintf(stderr, "*** Attempt to sample outside loaded range bin_shift=%d bin_shift=%lg, aborting\n", 
			bin_shift, si_local->bin_shift);
		exit(-1);
		}

	bin_shift2=rintf(si_local2->bin_shift)-pps->offset;
	if((bin_shift2+side_cut<0) || (bin_shift2+pps_bins>nbins)) {
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
//exit(-1);

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
SINGLE_BIN_LOOSELY_COHERENT_PATCH_PRIVATE_DATA *priv=(SINGLE_BIN_LOOSELY_COHERENT_PATCH_PRIVATE_DATA *)ctx->patch_private_data;

if(ctx->loose_first_half_count<0) {
	fprintf(stderr, "**** INTERNAL ERROR: loose_first_half_count=%d is not set.\n", ctx->loose_first_half_count);
	exit(-1);
	}

if(priv==NULL || priv->signature!=PATCH_PRIVATE_SINGLE_BIN_LOOSELY_COHERENT_SIGNATURE) {
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
	key+=((int)a-k);
	si_local->bin_shift=a*ctx->inv_cache_granularity;

	a=rintf(si_local->diff_bin_shift*ctx->diff_shift_granularity);
	//fprintf(stderr, "%g %g\n", a, si_local->diff_bin_shift);
	key=((int)a)+(key*KEY_MULT);
	si_local->diff_bin_shift=a*ctx->inv_diff_shift_granularity;

	key=(key*KEY_MULT) & KEY_DIV;

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
			if((fabs(si_local->bin_shift-(sc_si_local->bin_shift+first_shift))> ctx->half_inv_cache_granularity)
			   || (fabs(si_local->diff_bin_shift-sc_si_local->diff_bin_shift)> ctx->half_inv_diff_shift_granularity)) {
				match=0;
				fprintf(stderr, "OVERWRITE: i=%d key=%d count=%d %f %f %g %g %d\n", i, key, count, si_local->bin_shift, sc_si_local->bin_shift, si_local->diff_bin_shift, sc_si_local->diff_bin_shift, first_shift);
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

extern EphemerisData ephemeris;

LALDetector get_detector_struct(char *det) 
{
if(!strcasecmp("LHO", det)){
	return lalCachedDetectors[LALDetectorIndexLHODIFF];
	} else 
if(!strcasecmp("LLO", det)){
	return lalCachedDetectors[LALDetectorIndexLLODIFF];
	} else {
	fprintf(stderr,"Unrecognized detector site: \"%s\"\n", args_info.detector_arg);
	exit(-1);
	}
}


/* This function is meant to work with get_uncached_loose_partial_power_sum */
void accumulate_single_bin_loose_power_sums_sidereal_step(SUMMING_CONTEXT *ctx, POWER_SUM *ps, int count, double gps_start, double gps_stop, int veto_mask)
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
SINGLE_BIN_LOOSELY_COHERENT_PATCH_PRIVATE_DATA *priv;
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

	priv=do_alloc(1, sizeof(*priv));
	priv->signature=PATCH_PRIVATE_SINGLE_BIN_LOOSELY_COHERENT_SIGNATURE;
	priv->emission_time_size=count*segment_count;
	priv->emission_time=do_alloc(priv->emission_time_size, sizeof(*priv->emission_time));
	priv->inv_coherence_length=1.0/si->coherence_time;

	for(j=0;j<segment_count;j++) {
		/* assign index */
		si[j].index=j;

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
			baryinput.dInv=0; /* TODO: pass this from command line */

			LALBarycenter(&status, &(priv->emission_time[i*segment_count+j]), &baryinput, &earth_state);
			TESTSTATUS(&status);
			}
		}
	
	ctx->patch_private_data=priv;

	/* This assumes that we are patch bound !! *** WARNING ***
	TODO: make this assumption automatic in the data layout.
	 */
	si_local=si;
	min_shift=1000000; /* we should never load this many bins */
	max_shift=-1000000;
	for(j=0;j<segment_count;j++) {
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

	//group_count=ceil(max_shift-min_shift);
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

			ctx->reset_cache(ctx, tmp_count, count);
		
			/* loop over templates */
			ps_local=ps;
			for(i=0;i<count;i++) {
				/* fill in segment info appropriate to this template */
				si_local=tmp;
				priv->freq_shift=ps_local->freq_shift;
				priv->spindown=ps_local->spindown;

				for(j=0;j<tmp_count;j++) {
		// 			si[j].ra=ps[i].patch_ra;
		// 			si[j].dec=ps[i].patch_dec;
		// 			memcpy(si[j].e, ps[i].patch_e, GRID_E_COUNT*sizeof(SKY_GRID_TYPE));
		
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
	priv->emission_time=NULL;
	free(priv);
	ctx->patch_private_data=NULL;
	priv=NULL;
	}

for(i=0;i<count;i++) {
	if(ps[i].min_gps<0 || ps[i].min_gps>gps_start)ps[i].min_gps=gps_start;
	if(ps[i].max_gps<0 || ps[i].max_gps<gps_stop)ps[i].max_gps=gps_stop;
	}
if(status.statusPtr)FREESTATUSPTR(&status);
}
