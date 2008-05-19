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
#include <string.h>
#include <time.h>

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

//#define DISTANCE_PRINTF		fprintf(stderr, "rank=%d distance=%f alignment=%f fdist=%f sdist=%f snr=%f strain=%g\n", cand->rank, candidate_distance(cand, args_info.focus_ra_arg, args_info.focus_dec_arg), candidate_alignment(cand, args_info.fake_iota_arg, args_info.fake_psi_arg), (cand->frequency-args_info.fake_freq_arg)*1800, (cand->spindown-args_info.fake_spindown_arg)*(max_gps()-min_gps())*1800.0, cand->snr, cand->strain);

#define DISTANCE_PRINTF

#define WINDOW 240
#define SEARCH_WINDOW 10


extern DATASET *datasets;
extern int d_free;


extern POLARIZATION_RESULTS *polarization_results;
extern int nlinear_polarizations, ntotal_polarizations;

extern SKY_GRID *fine_grid, *patch_grid;
extern SKY_SUPERGRID *super_grid;
extern FILE *LOG;

extern double spindown;
extern double resolution;

extern int subinstance;
extern char *subinstance_name;

extern int nbins, first_bin, side_cut, useful_bins;

extern INT64 spindown_start;

CANDIDATE *candidate=NULL;
int candidate_free=0;
int candidate_size=-1;

extern struct gengetopt_args_info args_info;

char s[20000];

extern SUM_TYPE *max_dx;
extern short *max_dx_polarization_index;

int *max_dx_order=NULL;
int *dec_order=NULL;
int *inverse_dec_order=NULL;
int *ra_order=NULL;
int *inverse_ra_order=NULL;
int *max_dx_local_map=NULL;

float noise_floor;

float candidate_distance(CANDIDATE *cand, float ra, float dec)
{
return(acos(sin(cand->dec)*sin(dec)+cos(cand->dec)*cos(dec)*cos(cand->ra-ra)-1e-14));
}

float candidate_alignment(CANDIDATE *cand, float iota, float psi)
{
double a, b, a2, b2, c, dist;
a=cos(iota);
a2=fabs(a);
a2=(1.0-a2)/(1.0+a2);
a2=a2*a2;

b=cos(cand->iota);
b2=fabs(b);
b2=(1.0-b2)/(1.0+b2);
b2=b2*b2;

c=cand->psi-psi;
while(c>M_PI*0.5)c-=M_PI;
while(c< -M_PI*0.5)c+=M_PI;

dist=fabs(a2-b2)+0.5*fmin((1.0-a)*(1.0-a), (1.0-b)*(1.0-b))*fabs(c);

return(dist);
//return(acos(fabs(sin(cand->iota)*sin(iota)*cos(cand->psi-psi)+cos(cand->iota)*cos(iota))-1e-14));
}

int max_dx_compare(int *a1, int *a2) 
{
if(max_dx[*a1]<max_dx[*a2])return 1;
if(max_dx[*a1]>max_dx[*a2])return -1;
return 0;
}

int dec_compare(int *a1, int *a2) 
{
float ra_dist, dec_dist;
dec_dist=fine_grid->latitude[*a1]-fine_grid->latitude[*a2];
if(dec_dist>0)return 1;
if(dec_dist<0)return -1;

ra_dist=cos(fine_grid->latitude[*a1]+fine_grid->latitude[*a2])*(fine_grid->longitude[*a1]-fine_grid->longitude[*a2]);
if(ra_dist<0)return -1;
if(ra_dist>0)return 1;

return 0;
}

int ra_compare(int *a1, int *a2) 
{
float ra_dist, dec_dist;

ra_dist=cos(0.5*(fine_grid->latitude[*a1]+fine_grid->latitude[*a2]))*(fine_grid->longitude[*a1]-fine_grid->longitude[*a2]);
if(ra_dist<-1.5*resolution)return -1;
if(ra_dist>1.5*resolution)return 1;

dec_dist=fine_grid->latitude[*a1]-fine_grid->latitude[*a2];
if(dec_dist>0)return 1;
if(dec_dist<0)return -1;

return 0;
}

int min_i=0, max_i=0;

static int sweep_points_ra(int index, int mark) 
{
int j,m, count=0;
for(j=inverse_dec_order[index]+1;j<fine_grid->npoints;j++) {
	m=dec_order[j];
	if(max_dx_local_map[m]>=0)break;
	if(max_dx[m]<noise_floor)break;
	if(max_dx[m]>max_dx[dec_order[j-1]]+args_info.snr_precision_arg)break;
	max_dx_local_map[m]=mark;
	if(m<min_i)min_i=m;
	if(m>max_i)max_i=m;
	count++;
	}
for(j=inverse_dec_order[index]-1;j>=0;j--) {
	m=dec_order[j];
	if(max_dx_local_map[m]>=0)break;
	if(max_dx[m]<noise_floor)break;
	if(max_dx[m]>max_dx[dec_order[j+1]]+args_info.snr_precision_arg)break;
	max_dx_local_map[m]=mark;
	if(m<min_i)min_i=m;
	if(m>max_i)max_i=m;
	count++;
	}
return count;
}

static int sweep_points_dec(int index, int mark) 
{
int j,m, count=0;
for(j=inverse_ra_order[index]+1;j<fine_grid->npoints;j++) {
	m=ra_order[j];
	if(max_dx_local_map[m]>=0)break;
	if(max_dx[m]<noise_floor)break;
	if(fabs(fine_grid->latitude[m]-fine_grid->latitude[ra_order[j-1]])>0.1)break;
	if(max_dx[m]>max_dx[ra_order[j-1]]+args_info.snr_precision_arg)break;
	max_dx_local_map[m]=mark;
	if(m<min_i)min_i=m;
	if(m>max_i)max_i=m;
	count++;
	}
for(j=inverse_ra_order[index]-1;j>=0;j--) {
	m=ra_order[j];
	if(max_dx_local_map[m]>=0)break;
	if(max_dx[m]<noise_floor)break;
	if(fabs(fine_grid->latitude[m]-fine_grid->latitude[ra_order[j-1]])>0.1)break;
	if(max_dx[m]>max_dx[ra_order[j+1]]+args_info.snr_precision_arg)break;
	max_dx_local_map[m]=mark;
	if(m<min_i)min_i=m;
	if(m>max_i)max_i=m;
	count++;
	}
return count;
}

/* slow but sure.. */
static int sweep_neighbourhood(int index, int mark)
{
int j,m, count=0;
float fudge=resolution*3.1;
float cosfudge=1.0-cos(fudge);
float frequency;

if(max_dx_polarization_index[index]<0)return 0;

frequency=polarization_results[max_dx_polarization_index[index]].skymap.freq_map[index];

for(j=inverse_dec_order[index]+1;j<fine_grid->npoints;j++) {
	m=dec_order[j];

	if(fabs(fine_grid->latitude[m]-fine_grid->latitude[index])>fudge)break;

	if(max_dx_local_map[m]>=0)continue;
	
	if(max_dx[m]>max_dx[index])continue;
	
	if(max_dx_polarization_index[m]<0)continue;

	if(((max_dx[m]+1.0)>max_dx[index]) && 
		(fabs(frequency-polarization_results[max_dx_polarization_index[m]].skymap.freq_map[m])>3/1800.0) &&
		(fabs(polarization_results[max_dx_polarization_index[m]].skymap.freq_map[index]-polarization_results[max_dx_polarization_index[m]].skymap.freq_map[m])>3/1800.0)
		)continue;

	if(fast_spherical_distance(fine_grid->longitude[m], fine_grid->latitude[m],
				fine_grid->longitude[index], fine_grid->latitude[index])<cosfudge) {
		
		max_dx_local_map[m]=mark;
		count++;
		}
	}
	
for(j=inverse_dec_order[index]-1;j>=0;j--) {
	m=dec_order[j];

	if(fabs(fine_grid->latitude[m]-fine_grid->latitude[index])>fudge)break;

	if(max_dx_local_map[m]>=0)continue;

	if(max_dx[m]>max_dx[index])continue;
	
	if(max_dx_polarization_index[m]<0)continue;

	if(fast_spherical_distance(fine_grid->longitude[m], fine_grid->latitude[m],
				fine_grid->longitude[index], fine_grid->latitude[index])<cosfudge) {
	
		max_dx_local_map[m]=mark;
		count++;
		}
	}
return count;
}

void sweep_points(int index, int mark)
{
int i, count;
unsigned char *p;

max_i=index;
min_i=index;
max_dx_local_map[index]=mark;

p=do_alloc(fine_grid->npoints, sizeof(*p));
for(i=0;i<fine_grid->npoints;i++)*p=0;

while(1) {
	count=0;
	for(i=min_i;i<=max_i;i++) {
		if(p[i])continue;
		if(max_dx_local_map[i]==mark){
			count+=sweep_neighbourhood(i, mark);
			p[i]=1;
			}
		}
	fprintf(stderr, "mark=%d count=%d\n", mark, count);
	if(count==0){ 
		free(p);
		return;
		}
	}
}

void expand_candidates(void)
{
CANDIDATE *p;
candidate_size=candidate_size*2+10;
p=do_alloc(candidate_size, sizeof(*p));
if(candidate_free>0)memcpy(p, candidate, candidate_free*sizeof(*p));
free(candidate);
candidate=p;
}

void dump_candidate(int i)
{
float doppler;
FILE *fout;
char s[20];
DATASET *d;
int index=candidate[i].point_index;
float spindown=candidate[i].spindown;
int b0, b1, bin_shift, j, k, b;
float a, a_plus, a_cross;
POLARIZATION *pl;

snprintf(s, 20, "cand_%d.dat", i);
fprintf(stderr, "%s\n", s);
fout=fopen(s, "w");
fprintf(fout, "dataset\tk\tgps\tdoppler\tspindown\ta_plus\ta_cross\texpTMedian");
for(b=0;b<nbins-2*side_cut;b++)
	fprintf(fout, "\tb%d", b);
fprintf(fout, "\n");
/* loop over datasets */
for(j=0;j<d_free;j++) {
	d=&(datasets[j]);
	pl=&(d->polarizations[candidate[i].max_dx_polarization_index]);	
	
	/* process single SFT */
	for(k=0;k<d->free;k++){

		if(d->sft_veto[k])continue;

		/* Get amplitude response */
		a_plus=F_plus(k, fine_grid, index, pl->AM_coeffs);
		a_cross=F_plus(k, fine_grid, index, pl->conjugate->AM_coeffs);

		doppler=fine_grid->e[0][index]*d->detector_velocity[3*k+0]+
			fine_grid->e[1][index]*d->detector_velocity[3*k+1]+
			fine_grid->e[2][index]*d->detector_velocity[3*k+2];

		bin_shift=-rint((first_bin+nbins*0.5)*doppler+1800.0*spindown*(d->gps[k]-spindown_start+d->coherence_time*0.5));

		b0=side_cut-bin_shift;
		b1=(nbins-side_cut)-bin_shift;

		if((b0<0)||(b1>nbins)){
			fprintf(stderr,"Working frequency range obscured by bin_shift shift: bin_shift=%d index=%d i=%d\n",
				bin_shift, index, i);
			exit(-1);
			}
		fprintf(fout, "%s\t%d\t%lld\t%g\t%g\t%f\t%f\t%g", 
			d->name, 
			k, 
			d->gps[k], 
			doppler, 
			spindown*(d->gps[k]-spindown_start+d->coherence_time*0.5),
			a_plus,
			a_cross,
			d->expTMedians[k]);
		for(b=b0; b<b1;b++) {
			a=d->im[k*nbins+b];
			if(a>0)
				fprintf(fout, "\t%g+%gi", d->re[k*nbins+b], a);
				else 
				fprintf(fout, "\t%g%gi", d->im[k*nbins+b], a);
			}

		fprintf(fout, "\n");
		}
	}
fclose(fout);
}

void output_candidate_header(FILE *fout)
{
fprintf(fout, "candidates: label polarization_index rank opt_rank score point_index domain_size ul S M max_dx frequency psi iota ra dec spindown weight_ratio skyband coherence_score power_cor snr strain strain_err total_weight f_max ifo_freq ifo_freq_sd total\n");
}

void output_candidate(FILE *fout, char * suffix, CANDIDATE *cand)
{
fprintf(fout, "candidate%s: \"%s\" %d %d %d %g %d %d %g %g %g %f %f %f %f %f %f %g %f %d %f %f %f %g %g %f %f %f %f %d\n",
	suffix,
	args_info.label_given ? args_info.label_arg : "",
	cand->max_dx_polarization_index,
	cand->rank,
	cand->opt_rank,
	cand->score,
	cand->point_index,
	cand->domain_size,
	cand->ul,
	cand->S_band,
	cand->M_band,
	cand->max_dx,
	cand->frequency,
	cand->psi,
	cand->iota,
	cand->ra,
	cand->dec,
	cand->spindown,
	cand->weight_ratio,
	cand->skyband,
	cand->coherence_score,
	cand->power_cor,
	cand->snr,
	cand->strain,
	cand->strain_err,
	cand->total_weight,
	cand->f_max,
	cand->ifo_freq,
	cand->ifo_freq_sd,
	candidate_free);
}

#define SNR_FORMULA1 cand->snr=(demod_signal_sum[b_max]-mean_power)/power_sd;
#define SNR_FORMULA2 cand->snr=(demod_signal_sum[b_max]-0.25*fabs(demod_signal_sum[b_max-1]-demod_signal_sum[b_max+1])-mean_power)/power_sd;
#define SNR_FORMULA3 cand->snr=(0.6*demod_signal_sum[b_max]+0.15*(demod_signal_sum[b_max-1]+demod_signal_sum[b_max+1])-mean_power)/power_sd;
#define SNR_FORMULA4 cand->snr=(0.6*demod_signal_sum[b_max]+0.15*(demod_signal_sum[b_max-1]+demod_signal_sum[b_max+1])-mean_power)/power_sd;

#define SNR_FORMULA SNR_FORMULA1

typedef struct {
	/* Precomputed quantities and input parameters */
	CANDIDATE *cand;
	float *e;
	float a_plus;
	float a_cross;
	float a_plus_sq;
	float a_cross_sq;
	float p_proj;
	float c_proj;
	int remainder;
	int divisor;

	/* Accumulation variables */
	double demod_signal_sum[2*WINDOW+1];
	double signal_sum[2*WINDOW+1];
	double signal_sq_sum[2*WINDOW+1];
	double cov_sum[2*WINDOW+1];
	struct {
		double re;
		double im;
		} phase_sum[2*WINDOW+1];


	double total_weight;
	double response_sum;
	double response_sq_sum;
	double response_weight;
	double total_demod_weight;
	double total_demod_weight2;
	double ifo_freq;
	double ifo_freq_sq;

	int error;
	} COMPUTE_SCORES_DATA;

void compute_scores_cruncher(int thread_id, COMPUTE_SCORES_DATA *csd)
{
int b0, b1, j, k, b;
int *signal_bin;
float *bin_re, *bin_im;
double *response;
double *mismatch;
DATASET *d;
POLARIZATION *pl;
float a, f_plus, f_cross, f_plus_sq, f_cross_sq, doppler;
double f, x, y, demod_weight, power, weight;
float filter[7];
float spindown=csd->cand->spindown;
double frequency=csd->cand->frequency;

for(j=0;j<d_free;j++) {
	d=&(datasets[j]);
	//pl=&(d->polarizations[cand->max_dx_polarization_index]);	
	pl=&(d->polarizations[0]);	
	
	response=do_alloc(d->free, sizeof(*response));
	signal_bin=do_alloc(d->free, sizeof(*signal_bin));
	mismatch=do_alloc(d->free, sizeof(*mismatch));

	b0=nbins;
	b1=-1;
		
	/* process single SFTs */
	for(k=0;k<d->free;k++) {

		if(d->sft_veto[k])continue;

		/* skip sfts of other units */
		if((k % csd->divisor)!=csd->remainder)continue;

		/* Get amplitude response */
		#if 0
		a_plus=F_plus(k, fine_grid, index, pl->AM_coeffs);
		a_cross=F_plus(k, fine_grid, index, pl->conjugate->AM_coeffs);
		
		response[k]=pl->plus_factor*a_plus*a_plus+pl->cross_factor*a_cross*a_cross;
		#else 

		f_plus=F_plus_coeff(k, csd->e, pl->AM_coeffs);
		f_cross=F_plus_coeff(k, csd->e, pl->conjugate->AM_coeffs);

// 		f_plus_sq=f_plus*f_plus;
// 		f_cross_sq=f_cross*f_cross;
// 
// 		a=f_plus*a_plus+f_cross*a_cross;
// 		response[k]=0.25*((f_plus_sq+f_cross_sq)*(a_plus_sq+a_cross_sq)+(f_plus_sq-f_cross_sq)*(a_plus_sq-a_cross_sq)*p_proj+2.0*f_plus*f_cross*(a_plus_sq-a_cross_sq))*c_proj;

		a=f_plus*csd->p_proj+f_cross*csd->c_proj;
		f_cross=f_cross*csd->p_proj-f_plus*csd->c_proj;
		f_plus=a;

		f_plus_sq=f_plus*f_plus;
		f_cross_sq=f_cross*f_cross;

		//response[k]=0.25*((f_plus_sq+f_cross_sq)*(a_plus_sq+a_cross_sq)+(f_plus_sq-f_cross_sq)*(a_plus_sq-a_cross_sq));
		response[k]=f_plus_sq*csd->a_plus_sq+f_cross_sq*csd->a_cross_sq;
		#endif

		doppler=csd->e[0]*d->detector_velocity[3*k+0]+
			csd->e[1]*d->detector_velocity[3*k+1]+
			csd->e[2]*d->detector_velocity[3*k+2];

		f=frequency+frequency*doppler+spindown*(d->gps[k]-spindown_start+d->coherence_time*0.5);

		signal_bin[k]=rintf(1800.0*f-first_bin);
		mismatch[k]=1800.0*f-first_bin-signal_bin[k];

		if((signal_bin[k]<b0))b0=signal_bin[k];
		if((signal_bin[k]>b1))b1=signal_bin[k];

		}

	if( (b0<0) || (b1>=nbins)) {
		csd->error=-1;
		free(response);
		free(signal_bin);
		free(mismatch);
		return;
		}
	//fprintf(stderr, "b0=%d b1=%d\n", b0, b1);

	for(k=0;k<d->free;k++) {

		if(d->sft_veto[k])continue;

		/* skip sfts of other units */
		if((k % csd->divisor)!=csd->remainder)continue;

		/* skip SFTs with low weight */
		//if(d->expTMedians[k]*d->weight*response[k]<0.05)continue;

		/* power_cor computation */
		weight=d->expTMedians[k]*d->weight;
		csd->response_sum+=response[k]*weight;
		csd->response_sq_sum+=response[k]*response[k]*weight;
		csd->response_weight+=weight;

		demod_weight=weight*response[k];
		csd->total_demod_weight2+=demod_weight;
		csd->total_demod_weight+=demod_weight*response[k];

		bin_re=&(d->re[k*nbins+signal_bin[k]-WINDOW]);
		bin_im=&(d->im[k*nbins+signal_bin[k]-WINDOW]);

		f=signal_bin[k]+mismatch[k];
		csd->ifo_freq+=f*demod_weight;
		csd->ifo_freq_sq+=f*f*demod_weight;

		fill_hann_filter7(filter, mismatch[k]);

		for(b=0; b< (2*WINDOW+1); b++) {
			x=bin_re[-3]*filter[0]+bin_re[-2]*filter[1]+bin_re[-1]*filter[2]+bin_re[0]*filter[3]+bin_re[1]*filter[4]+bin_re[2]*filter[5]+bin_re[3]*filter[6];
			y=bin_im[-3]*filter[0]+bin_im[-2]*filter[1]+bin_im[-1]*filter[2]+bin_im[0]*filter[3]+bin_im[1]*filter[4]+bin_im[2]*filter[5]+bin_im[3]*filter[6];

			power=x*x+y*y;
			#if 1
			csd->signal_sum[b]+=power*weight;
			csd->signal_sq_sum[b]+=power*power*weight;
			csd->cov_sum[b]+=power*response[k]*weight;
			#endif

			csd->demod_signal_sum[b]+=power*demod_weight;
			bin_re++;
			bin_im++;
			}
		}
		
	/* phase computation */
	#if 0
	for(k=1;k< (d->free-1);k++) {
	
	
		/* skip SFT with asymmetric gaps */
		if(d->gps[k]-d->gps[k-1]!=d->gps[k+1]-d->gps[k])continue;

		/* skip SFTs with low weight */
		if(d->expTMedians[k]<0.01 || d->expTMedians[k-1]<0.01 || d->expTMedians[k+1]<0.01)continue;
	
		weight=((response[k-1]*d->expTMedians[k-1]+2*response[k]*d->expTMedians[k]+response[k+1]*d->expTMedians[k+1])*d->weight);
		csd->total_weight+=weight;
	
		p0=&(d->bin[(k-1)*nbins+signal_bin[k-1]-WINDOW]);
		p1=&(d->bin[k*nbins+signal_bin[k]-WINDOW]);
		p2=&(d->bin[(k+1)*nbins+signal_bin[k+1]-WINDOW]);

		for(b=0; b< (2*WINDOW+1); b++) {
			x=sqrt(p0[b].re*p0[b].re+p0[b].im*p0[b].im);
			y=sqrt(p1[b].re*p1[b].re+p1[b].im*p1[b].im);

			w.re=(p0[b].re*p1[b].re+p0[b].im*p1[b].im)/(x*y);
			w.im=(p0[b].re*p1[b].im-p0[b].im*p1[b].re)/(x*y);

			x=sqrt(p2[b].re*p2[b].re+p2[b].im*p2[b].im);
			z.re=(p2[b].re*p1[b].re+p2[b].im*p1[b].im)/(x*y);
			z.im=(p2[b].re*p1[b].im-p2[b].im*p1[b].re)/(x*y);
			
			csd->phase_sum[b].re+=weight*(w.re*z.re-w.im*z.im);
			csd->phase_sum[b].im+=weight*(w.im*z.re+w.re*z.im);
			}
		}
	#endif
		
	free(response);
	free(signal_bin);
	free(mismatch);
	response=NULL;
	signal_bin=NULL;
	mismatch=NULL;
	}
}

void compute_scores(CANDIDATE *cand, int debug)
{
int b0, b1, i, b, b_max;
float a, coherence_score, power_cor, a_plus, a_cross, a_plus_sq, a_cross_sq;
double total_weight, *demod_signal_sum, *signal_sum, *signal_sq_sum, response_sum, response_sq_sum, response_weight, total_demod_weight, total_demod_weight2, *cov_sum, mean_power, mean_power_sq, power_sd;
double ifo_freq, ifo_freq_sq;
float p_proj, c_proj;
float e[26];
struct {
	double re;
	double im;
	} *phase_sum;
int n_units;
COMPUTE_SCORES_DATA *units;

total_weight=0.0;
response_sq_sum=0.0;
response_sum=0.0;
response_weight=0.0;
total_demod_weight=0.0;
total_demod_weight2=0.0;

phase_sum=do_alloc(2*WINDOW+1, sizeof(*phase_sum));
demod_signal_sum=do_alloc(2*WINDOW+1, sizeof(*signal_sum));
signal_sum=do_alloc(2*WINDOW+1, sizeof(*signal_sum));
signal_sq_sum=do_alloc(2*WINDOW+1, sizeof(*signal_sq_sum));
cov_sum=do_alloc(2*WINDOW+1, sizeof(*cov_sum));

precompute_am_constants(e, cand->ra, cand->dec);
p_proj=cos(2*cand->psi);
c_proj=sin(2*cand->psi);

a=cos(cand->iota);
a_plus=(1.0+a*a)/2.0;
a_cross=a;

/* precompute squares */
a_plus_sq=a_plus*a_plus;
a_cross_sq=a_cross*a_cross;

//fprintf(stderr, "psi=%f %f %f\n", cand->psi, p_proj, c_proj);


/* loop over datasets */
n_units=get_max_threads()*2;
units=do_alloc(n_units, sizeof(*units));

for(i=0;i<n_units;i++) {
	units[i].e=e;
	units[i].a_plus=a_plus;
	units[i].a_cross=a_cross;
	units[i].a_plus_sq=a_plus_sq;
	units[i].a_cross_sq=a_cross_sq;
	units[i].p_proj=p_proj;
	units[i].c_proj=c_proj;

	units[i].cand=cand;

	units[i].remainder=i;
	units[i].divisor=n_units;


	for(b=0;b<2*WINDOW+1;b++) {
		units[i].phase_sum[b].re=0.0;
		units[i].phase_sum[b].im=0.0;
		units[i].signal_sum[b]=0.0;
		units[i].signal_sq_sum[b]=0.0;
		units[i].cov_sum[b]=0.0;
		units[i].demod_signal_sum[b]=0.0;
		}
	units[i].total_weight=0.0;
	units[i].response_sq_sum=0.0;
	units[i].response_sum=0.0;
	units[i].response_weight=0.0;
	units[i].total_demod_weight=0.0;
	units[i].total_demod_weight2=0.0;

	units[i].ifo_freq=0.0;
	units[i].ifo_freq_sq=0.0;

	units[i].error=0;

	submit_job(compute_scores_cruncher, &(units[i]));
	}

cand->ifo_freq=0.0;
cand->ifo_freq_sd=0.0;
ifo_freq_sq=0.0;
ifo_freq=0.0;

for(b=0;b<2*WINDOW+1;b++) {
	phase_sum[b].re=0.0;
	phase_sum[b].im=0.0;
	signal_sum[b]=0.0;
	signal_sq_sum[b]=0.0;
	cov_sum[b]=0.0;
	demod_signal_sum[b]=0.0;
	}

while(do_single_job(-1));

wait_for_all_done();

for(i=0;i<n_units;i++) {
	if(units[i].error) {
		/* we are outside loaded bin range */
		fprintf(stderr, "Outside loaded bin range: b0=%d b1=%d\n", b0, b1);
		cand->coherence_score=-1.0;
		cand->power_cor=-1.0;
		cand->snr=-1.0;
		cand->strain=-1.0;
		cand->f_max=cand->frequency;
		cand->ifo_freq=-1;
		cand->ifo_freq_sd=-1;

		free(units);
		free(phase_sum);
		free(signal_sum);
		free(signal_sq_sum);
		free(cov_sum);
		free(demod_signal_sum);
		return;
		}

	for(b=0;b<2*WINDOW+1;b++) {
		phase_sum[b].re+=units[i].phase_sum[b].re;
		phase_sum[b].im+=units[i].phase_sum[b].im;
		signal_sum[b]+=units[i].signal_sum[b];
		signal_sq_sum[b]+=units[i].signal_sq_sum[b];
		cov_sum[b]+=units[i].cov_sum[b];
		demod_signal_sum[b]+=units[i].demod_signal_sum[b];
		}
	total_weight+=units[i].total_weight;
	response_sq_sum+=units[i].response_sq_sum;
	response_sum+=units[i].response_sum;
	response_weight+=units[i].response_weight;
	total_demod_weight+=units[i].total_demod_weight;
	total_demod_weight2+=units[i].total_demod_weight2;

	ifo_freq+=units[i].ifo_freq;
	ifo_freq_sq+=units[i].ifo_freq_sq;
	}


//fprintf(stderr, "total_weight=%g\n", total_weight);

coherence_score=0;
power_cor=0;
response_sum/=response_weight;
response_sq_sum/=response_weight;

ifo_freq/=total_demod_weight2;
ifo_freq_sq/=total_demod_weight2;

cand->ifo_freq=(ifo_freq+first_bin)/1800.0;
cand->ifo_freq_sd=sqrt(ifo_freq_sq-ifo_freq*ifo_freq);

for(b=0;b<2*WINDOW+1;b++) {
	/* power_cor */
	cov_sum[b]/=response_weight;
	signal_sum[b]/=response_weight;
	signal_sq_sum[b]/=response_weight;
	demod_signal_sum[b]/=total_demod_weight;
	}

b_max=WINDOW;
for(b=WINDOW-SEARCH_WINDOW;b<=WINDOW+SEARCH_WINDOW+1;b++) {
	if(demod_signal_sum[b_max]<demod_signal_sum[b]) b_max=b;
	
	a=(cov_sum[b]-signal_sum[b]*response_sum)/sqrt((signal_sq_sum[b]-signal_sum[b]*signal_sum[b])*(response_sq_sum-response_sum*response_sum));
	if(fabs(a)>fabs(power_cor))power_cor=a;

	/* coherence_score */
	phase_sum[b].re/=total_weight;
	phase_sum[b].im/=total_weight;
	a=sqrt(phase_sum[b].re*phase_sum[b].re+phase_sum[b].im*phase_sum[b].im);
	if(a>coherence_score)coherence_score=a;
	}

cand->f_max=cand->frequency+(b_max-WINDOW)/1800.0;

mean_power=0.0;
mean_power_sq=0.0;
for(b=0;b<WINDOW-SEARCH_WINDOW;b++){
	mean_power+=demod_signal_sum[b];
	mean_power_sq+=demod_signal_sum[b]*demod_signal_sum[b];
	}
for(b=WINDOW+SEARCH_WINDOW+1;b<2*WINDOW+1;b++){
	mean_power+=demod_signal_sum[b];
	mean_power_sq+=demod_signal_sum[b]*demod_signal_sum[b];
	}
	
mean_power/=2.0*(WINDOW-SEARCH_WINDOW);
mean_power_sq/=2.0*(WINDOW-SEARCH_WINDOW);
power_sd=sqrt((mean_power_sq-mean_power*mean_power)*2.0*(WINDOW-SEARCH_WINDOW)/(2.0*(WINDOW-SEARCH_WINDOW)-1));

if(debug) {
	fprintf(stderr, "total_demod_weight=%f mean_power=%g mean_power_sq=%g power_sd=%g\n", total_demod_weight, mean_power, mean_power_sq, power_sd);
	fprintf(stderr, "Power: ");
	for(b=0;b<2*WINDOW+1;b++) {
		fprintf(stderr, "%.2f ", (demod_signal_sum[b]-mean_power)/power_sd);
		}
	fprintf(stderr, "\n");
	}


b_max=WINDOW;

cand->coherence_score=coherence_score;
cand->power_cor=power_cor;
/*cand->snr=(demod_signal_sum[b_max]-mean_power)/power_sd;*/
cand->total_weight=total_demod_weight;
/*cand->snr=(demod_signal_sum[b_max]-0.25*fabs(demod_signal_sum[b_max-1]-demod_signal_sum[b_max+1])-mean_power)/power_sd;*/
//cand->snr=(demod_signal_sum[b_max]-0.25*fabs(demod_signal_sum[b_max-1]-demod_signal_sum[b_max+1])-0.5*mean_power)/power_sd;
SNR_FORMULA

if(demod_signal_sum[b_max]< mean_power)cand->strain=0.0;
	else
	cand->strain=2.0*sqrt(demod_signal_sum[b_max]-mean_power)*args_info.strain_norm_factor_arg/(1800.0*16384.0);
if(cand->strain>0)
	cand->strain_err=0.5*4.0*power_sd*args_info.strain_norm_factor_arg/(cand->strain*1800.0*1800.0*16384.0*16384.0);
	else
	cand->strain_err=1;

free(units);
free(phase_sum);
free(signal_sum);
free(signal_sq_sum);
free(cov_sum);
free(demod_signal_sum);
}

typedef struct {
	long *offset;
	float *response;
	float *doppler;
// 	float *mismatch;
// 	int *signal_bin;

	float e[26];

	/* max/min values needed to shift path by 1 bin */
	float d_spindown[2];
	float d_ra[2];
	float d_dec[2];
	float d_freq[2];

	int valid;
	} SCORE_AUX_DATA;

#define VALID_RESPONSE	1
#define VALID_SHIFT	2
#define VALID_DOPPLER	4
#define VALID_E		8
#define VALID_DIFF_SHIFT 16

#define ALL_VALID (VALID_RESPONSE | VALID_SHIFT | VALID_DOPPLER)

static SCORE_AUX_DATA * allocate_score_aux_data(void)
{
SCORE_AUX_DATA *a;
int i;
a=do_alloc(1, sizeof(*a));

a->offset=do_alloc(d_free+1, sizeof(*a->offset));
a->offset[0]=0;
for(i=0;i<d_free;i++)a->offset[i+1]=a->offset[i]+datasets[i].free;

a->response=do_alloc(a->offset[d_free], sizeof(a->response));
a->doppler=do_alloc(a->offset[d_free], sizeof(a->doppler));
/*a->mismatch=do_alloc(a->offset[d_free], sizeof(a->mismatch));
a->signal_bin=do_alloc(a->offset[d_free], sizeof(a->signal_bin));*/
a->valid=0;

return(a);
}

static void free_score_aux_data(SCORE_AUX_DATA *a)
{
free(a->offset);
free(a->response);
free(a->doppler);
/*free(a->mismatch);
free(a->signal_bin);*/
free(a);
}

static void fill_e(SCORE_AUX_DATA *ad, CANDIDATE *cand)
{
precompute_am_constants(ad->e, cand->ra, cand->dec);
ad->valid |= VALID_E;
}

static void fill_response(SCORE_AUX_DATA *ad, CANDIDATE *cand)
{
int j,k;
float p_proj, c_proj;
float a, a_plus, a_cross, a_plus_sq, a_cross_sq, f_plus, f_cross, f_plus_sq, f_cross_sq;
float *response;
POLARIZATION *pl;
DATASET *d;

if(! (ad->valid & VALID_E) ) {
	fprintf(stderr, "E vector has not been computed !\n");
	exit(-1);
	}

p_proj=cos(2*cand->psi);
c_proj=sin(2*cand->psi);

a=cos(cand->iota);
a_plus=(1.0+a*a)/2.0;
a_cross=a;

/* precompute squares */
a_plus_sq=a_plus*a_plus;
a_cross_sq=a_cross*a_cross;

for(j=0;j<d_free;j++) {
	d=&(datasets[j]);
	//pl=&(d->polarizations[cand->max_dx_polarization_index]);	
	pl=&(d->polarizations[0]);	
	
	response=&(ad->response[ad->offset[j]]);	
	/* process single SFTs */
	for(k=0;k<d->free;k++) {

		if(d->sft_veto[k])continue;

		f_plus=F_plus_coeff(k, ad->e, pl->AM_coeffs);
		f_cross=F_plus_coeff(k, ad->e, pl->conjugate->AM_coeffs);

		a=f_plus*p_proj+f_cross*c_proj;
		f_cross=f_cross*p_proj-f_plus*c_proj;
		f_plus=a;

		f_plus_sq=f_plus*f_plus;
		f_cross_sq=f_cross*f_cross;

		response[k]=f_plus_sq*a_plus_sq+f_cross_sq*a_cross_sq;
		}
	}
ad->valid |= VALID_RESPONSE;
}

static void fill_doppler(SCORE_AUX_DATA *ad)
{
int j,k;
float *doppler;
DATASET *d;

if(! (ad->valid & VALID_E) ) {
	fprintf(stderr, "E vector has not been computed !\n");
	exit(-1);
	}

for(j=0;j<d_free;j++) {
	d=&(datasets[j]);

	doppler=&(ad->doppler[ad->offset[j]]);	

	/* process single SFTs */
	for(k=0;k<d->free;k++) {

		if(d->sft_veto[k])continue;

		doppler[k]=ad->e[0]*d->detector_velocity[3*k+0]+
			ad->e[1]*d->detector_velocity[3*k+1]+
			ad->e[2]*d->detector_velocity[3*k+2];
		}

	}
ad->valid |= VALID_DOPPLER;
}

static void fill_diff_shift(SCORE_AUX_DATA *ad, CANDIDATE *cand)
{
float doppler, ra_doppler, dec_doppler;
float *response;
float demod_weight, total_weight;
float spindown=cand->spindown;
double frequency=cand->frequency;
double f;
float signal_bin, mismatch;
float a;
float one1800;
DATASET *d;
int j, k;

if(! (ad->valid & VALID_E) || !(ad->valid & VALID_RESPONSE) ) {
	fprintf(stderr, "E vector or response has not been computed !\n");
	exit(-1);
	}

ad->d_freq[0]=0;
ad->d_freq[1]=0;

ad->d_spindown[0]=0;
ad->d_spindown[1]=0;

ad->d_ra[0]=0;
ad->d_ra[1]=0;

ad->d_dec[0]=0;
ad->d_dec[1]=0;

total_weight=0;

for(j=0;j<d_free;j++) {
	d=&(datasets[j]);

	response=&(ad->response[ad->offset[j]]);

	one1800=1.0/d->coherence_time;

	/* process single SFTs */
	for(k=0;k<d->free;k++) {

		if(d->sft_veto[k])continue;

		demod_weight=d->expTMedians[k]*d->weight*response[k];
		total_weight+=demod_weight;

		doppler=ad->e[0]*d->detector_velocity[3*k+0]+
			ad->e[1]*d->detector_velocity[3*k+1]+
			ad->e[2]*d->detector_velocity[3*k+2];

		ra_doppler=-ad->e[1]*d->detector_velocity[3*k+0]+
			ad->e[0]*d->detector_velocity[3*k+1];

		dec_doppler=-ad->e[2]*ad->e[4]*d->detector_velocity[3*k+0]
			-ad->e[2]*ad->e[5]*d->detector_velocity[3*k+1]+
			ad->e[3]*d->detector_velocity[3*k+2];


		f=frequency+frequency*doppler+spindown*(d->gps[k]-spindown_start+d->coherence_time*0.5);
		signal_bin=rintf(d->coherence_time*f)*one1800;

		mismatch=(f-signal_bin);

		if(mismatch<0)mismatch=0.5*one1800+mismatch;

		mismatch+=0.25*one1800;

		ad->d_freq[0]+=-demod_weight/mismatch;
		ad->d_freq[1]+=demod_weight/(one1800-mismatch);
		
		a=fabs(d->gps[k]-spindown_start+d->coherence_time*0.5);
		if(a>0) {
			ad->d_spindown[0]+=-demod_weight*a/mismatch;
			ad->d_spindown[1]+=demod_weight*a/(one1800-mismatch);
			}

		a=fabs(frequency*ra_doppler);
		if(fabs(a)>0) {
			ad->d_ra[0]+=-demod_weight*a/mismatch;
			ad->d_ra[1]+=demod_weight*a/(one1800-mismatch);
			}

		a=fabs(frequency*dec_doppler);
		if(fabs(a)>0) {
			ad->d_dec[0]+=-demod_weight*a/mismatch;
			ad->d_dec[1]+=demod_weight*a/(one1800-mismatch);
			}
		}

	}

ad->d_freq[0]=total_weight/ad->d_freq[0];
ad->d_freq[1]=total_weight/ad->d_freq[1];

ad->d_spindown[0]=total_weight/ad->d_spindown[0];
ad->d_spindown[1]=total_weight/ad->d_spindown[1];

ad->d_ra[0]=total_weight/ad->d_ra[0];
ad->d_ra[1]=total_weight/ad->d_ra[1];

ad->d_dec[0]=total_weight/ad->d_dec[0];
ad->d_dec[1]=total_weight/ad->d_dec[1];

ad->valid |= VALID_DIFF_SHIFT;
}

static void fill_matched_diff_shift(SCORE_AUX_DATA *ad, CANDIDATE *cand)
{
float doppler, ra_doppler, dec_doppler;
float *response;
float demod_weight, total_weight;
float spindown=cand->spindown;
double frequency=cand->frequency;
double f;
float mismatch;
float a;
float one1800;
DATASET *d;
int i, j, k, signal_bin;
float filter[7];
float diff_filter[7];
float *bin_re, *bin_im;
float x, y, dx, dy;
double power, dpower;

if(! (ad->valid & VALID_E) || !(ad->valid & VALID_RESPONSE) ) {
	fprintf(stderr, "E vector or response has not been computed !\n");
	exit(-1);
	}

ad->d_freq[0]=0;
ad->d_freq[1]=0;

ad->d_spindown[0]=0;
ad->d_spindown[1]=0;

ad->d_ra[0]=0;
ad->d_ra[1]=0;

ad->d_dec[0]=0;
ad->d_dec[1]=0;

total_weight=0;

for(j=0;j<d_free;j++) {
	d=&(datasets[j]);

	response=&(ad->response[ad->offset[j]]);

	one1800=1.0/d->coherence_time;
	/* TODO */

	/* process single SFTs */
	for(k=0;k<d->free;k++) {

		if(d->sft_veto[k])continue;

		demod_weight=d->expTMedians[k]*d->weight*response[k];

		doppler=ad->e[0]*d->detector_velocity[3*k+0]+
			ad->e[1]*d->detector_velocity[3*k+1]+
			ad->e[2]*d->detector_velocity[3*k+2];

		ra_doppler=-ad->e[1]*d->detector_velocity[3*k+0]+
			ad->e[0]*d->detector_velocity[3*k+1];

		dec_doppler=-ad->e[2]*ad->e[4]*d->detector_velocity[3*k+0]
			-ad->e[2]*ad->e[5]*d->detector_velocity[3*k+1]+
			ad->e[3]*d->detector_velocity[3*k+2];


		f=frequency+frequency*doppler+spindown*(d->gps[k]-spindown_start+d->coherence_time*0.5);
		signal_bin=rintf(d->coherence_time*f-first_bin);

		mismatch=f*d->coherence_time-first_bin-signal_bin;

		fill_hann_filter7(filter, mismatch);
		fill_diff_hann_filter7(diff_filter, mismatch);

		bin_re=&(d->re[k*nbins+signal_bin-3]);
		bin_im=&(d->im[k*nbins+signal_bin-3]);

		x=0;
		y=0;
		dx=0;
		dy=0;

		for(i=0;i<7;i++) {
			x+=bin_re[i]*filter[i];
			y+=bin_im[i]*filter[i];
			dx+=bin_re[i]*diff_filter[i];
			dy+=bin_im[i]*diff_filter[i];
			}

		power=x*x+y*y;
		dpower=2.0*(x*dx+y*dy)*d->coherence_time;

		total_weight+=power*demod_weight;

		//ad->d_freq[0]+=-dpower*demod_weight;
		ad->d_freq[1]+=(1.0+doppler)*dpower*demod_weight;
		
		a=d->gps[k]-spindown_start+d->coherence_time*0.5;
		//ad->d_spindown[0]+=-demod_weight*a*dpower;
		ad->d_spindown[1]+=demod_weight*a*dpower;

		a=frequency*ra_doppler;
		//ad->d_ra[0]+=-demod_weight*a*dpower;
		ad->d_ra[1]+=demod_weight*a*dpower;

		a=frequency*dec_doppler;
		//ad->d_dec[0]+=-demod_weight*a*dpower;
		ad->d_dec[1]+=demod_weight*a*dpower;
		}

	}

ad->d_freq[0]=-ad->d_freq[1]/total_weight;
ad->d_freq[1]=ad->d_freq[1]/total_weight;

ad->d_spindown[0]=-ad->d_spindown[1]/total_weight;
ad->d_spindown[1]=ad->d_spindown[1]/total_weight;

ad->d_ra[0]=-ad->d_ra[1]/total_weight;
ad->d_ra[1]=ad->d_ra[1]/total_weight;

ad->d_dec[0]=-ad->d_dec[1]/total_weight;
ad->d_dec[1]=ad->d_dec[1]/total_weight;

ad->valid |= VALID_DIFF_SHIFT;
}

static void print_diff_shift(FILE *f, SCORE_AUX_DATA *ad)
{
return;
fprintf(f, "ra=[%f,%f] dec=[%f,%f] freq=[%g,%g] spindown=[%g,%g]\n",
	ad->d_ra[0], ad->d_ra[1],
	ad->d_dec[0], ad->d_dec[1],
	ad->d_freq[0], ad->d_freq[1],
	ad->d_spindown[0], ad->d_spindown[1]
	);
}

#if 0
static void fill_shift(SCORE_AUX_DATA *ad, CANDIDATE *cand)
{
int j,k, b0,b1;
int *sb;
float *m;
float f, *doppler;
DATASET *d;
float spindown=cand->spindown;
float frequency=cand->frequency;

if(!(ad->valid & VALID_DOPPLER)) {
	ad->valid &= ~VALID_SHIFT;
	return;
	}

for(j=0;j<d_free;j++) {
	d=&(datasets[j]);
	
	doppler=&(ad->doppler[ad->offset[j]]);
	sb=&(ad->signal_bin[ad->offset[j]]);
	m=&(ad->mismatch[ad->offset[j]]);

	/* process single SFTs */
	for(k=0;k<d->free;k++) {

		if(d->sft_veto[k])continue;

		f=frequency+frequency*doppler[k]+spindown*(d->gps[k]-spindown_start+d->coherence_time*0.5);

		sb[k]=rint(1800.0*f-first_bin);
		m[k]=1800.0*f-first_bin-sb[k];

		if(!k || (sb[k]<b0))b0=sb[k];
		if(!k || (sb[k]>b1))b1=sb[k];

		}

	if( (b0<0) || ((b1+1)>=nbins)) {
		ad->valid &= ~VALID_SHIFT;
		return;
		}
	}
ad->valid |= VALID_SHIFT;
}
#endif

inline static void compute_simple_snr(SCORE_AUX_DATA *ad, CANDIDATE *cand, int debug)
{
DATASET *d;
float spindown=cand->spindown;
double frequency=cand->frequency;
int j, k, b, b_max;
double weight, f, total_demod_weight, mean_power, mean_power_sq, power_sd;
float *response;
float *doppler;
int signal_bin;
float demod_signal_sum[(2*WINDOW+1)*2], *x, *y, *efm, power, *pout, demod_weight;

if((ad->valid & (VALID_RESPONSE | VALID_DOPPLER))!=(VALID_RESPONSE | VALID_DOPPLER)) {
	cand->strain=-1;
	cand->snr=-1;
	cand->f_max=cand->frequency;
	fprintf(stderr, "compute_snr: score_aux_data is not valid\n");
	return;
	}

total_demod_weight=0.0;

for(b=0;b<(2*WINDOW+1);b++) {
	demod_signal_sum[b]=0.0;
	}


cand->ifo_freq=0.0;
cand->ifo_freq_sd=0.0;

/* loop over datasets */
for(j=0;j<d_free;j++) {
	d=&(datasets[j]);

	response=&(ad->response[ad->offset[j]]);
	doppler=&(ad->doppler[ad->offset[j]]);

	for(k=0;k<d->free;k++) {
		if(d->sft_veto[k])continue;

		/* skip SFTs with low weight */
		if(d->expTMedians[k]<0.05)continue;

		/* power_cor computation */
		weight=d->expTMedians[k]*d->weight;

		f=frequency+frequency*doppler[k]+spindown*(d->gps[k]-spindown_start+d->coherence_time*0.5);
		signal_bin=rintf(1800.0*f-first_bin);

		if(signal_bin<WINDOW) {
			cand->strain=-1;
			cand->snr=-1;
			cand->f_max=cand->frequency;
			fprintf(stderr, "compute_snr: outside loaded frequency range\n");
			}

		demod_weight=weight*response[k];
		total_demod_weight+=demod_weight*response[k];

		x=&(d->re[k*nbins+signal_bin-WINDOW]);
		y=&(d->im[k*nbins+signal_bin-WINDOW]);
		efm=&(d->expFMedians_plain[signal_bin-WINDOW]);
		pout=demod_signal_sum;

		for(b=0; b< (2*WINDOW+1); b++) {
			power=(*x)*(*x)+(*y)*(*y);

			if(args_info.subtract_background_arg) {
				power-=d->expTMedians_plain[k]*(*efm);
				}

			(*pout)+=(power)*demod_weight;
			pout++;
			x++;
			y++;
			efm++;
			}
		}
	}

total_demod_weight=1.0/total_demod_weight;
for(b=0;b<2*WINDOW+1;b++) {
	demod_signal_sum[b]*=total_demod_weight;
	}

b_max=WINDOW;
for(b=WINDOW-SEARCH_WINDOW;b<WINDOW+SEARCH_WINDOW+1;b++) {
	if(demod_signal_sum[b_max]<demod_signal_sum[b]) b_max=b;
	}

mean_power=0.0;
mean_power_sq=0.0;
for(b=0;b<WINDOW-SEARCH_WINDOW;b++){
	mean_power+=demod_signal_sum[b];
	mean_power_sq+=demod_signal_sum[b]*demod_signal_sum[b];
	}
for(b=WINDOW+SEARCH_WINDOW+1;b<2*WINDOW+1;b++){
	mean_power+=demod_signal_sum[b];
	mean_power_sq+=demod_signal_sum[b]*demod_signal_sum[b];
	}
	
mean_power*=0.5/(WINDOW-SEARCH_WINDOW);
mean_power_sq*=0.5/(WINDOW-SEARCH_WINDOW);
power_sd=sqrt((mean_power_sq-mean_power*mean_power)*(2.0*(WINDOW-SEARCH_WINDOW)/(2.0*(WINDOW-SEARCH_WINDOW)-1)));

if(debug) {
	fprintf(stderr, "mean_power=%g mean_power_sq=%g power_sd=%g\n", mean_power, mean_power_sq, power_sd);
	fprintf(stderr, "Power: ");
	for(b=0;b<2*WINDOW+1;b++) {
		fprintf(stderr, "%.2f ", (demod_signal_sum[b]-mean_power)/power_sd);
		}
	fprintf(stderr, "\n");
	}

cand->f_max=cand->frequency+(b_max-WINDOW)/1800.0;

b_max=WINDOW;
cand->total_weight=total_demod_weight;
/*cand->snr=(demod_signal_sum[b_max]-mean_power)/power_sd;
cand->snr=(demod_signal_sum[b_max]-0.25*fabs(demod_signal_sum[b_max-1]-demod_signal_sum[b_max+1])-mean_power)/power_sd;*/
SNR_FORMULA

if(demod_signal_sum[b_max]< mean_power)cand->strain=0.0;
	else
	cand->strain=2.0*sqrt(demod_signal_sum[b_max]-mean_power)*args_info.strain_norm_factor_arg/(1800.0*16384.0);

if(cand->strain>0)
	cand->strain_err=0.5*4.0*power_sd*args_info.strain_norm_factor_arg/(cand->strain*(1800.0*1800.0*16384.0*16384.0));
	else
	cand->strain_err=1;

}

void single_pass_compute_simple_snr(CANDIDATE *cand)
{
SCORE_AUX_DATA *ad;

ad=allocate_score_aux_data();

fill_e(ad, cand);
fill_response(ad, cand);
fill_doppler(ad);

compute_simple_snr(ad, cand, 0);

free_score_aux_data(ad);
}

typedef struct {
	double demod_signal_sum_d[(2*WINDOW+1)];
	double total_demod_weight;
	SCORE_AUX_DATA *ad;
	CANDIDATE *cand;
	int remainder;
	int divisor;
	int error;
	} MATCHED_SNR_UNIT;

static void compute_matched_snr_cruncher(int thread_id, MATCHED_SNR_UNIT *msu)
{
float demod_signal_sum[(2*WINDOW+1)], power, *pout, demod_weight, x, y;
float *response;
float *doppler;
int signal_bin;
int j, k, b;
DATASET *d;
float *bin_re, *bin_im;
float mismatch;
int count;
double f, weight;
float filter[7];
float spindown=msu->cand->spindown;
double frequency=msu->cand->frequency;

for(b=0;b<(2*WINDOW+1);b++) {
	demod_signal_sum[b]=0.0;
	}

for(j=0;j<d_free;j++) {
	d=&(datasets[j]);

	response=&(msu->ad->response[msu->ad->offset[j]]);
	doppler=&(msu->ad->doppler[msu->ad->offset[j]]);

	count=0;

	for(k=0;k<d->free;k++) {

		if(d->sft_veto[k])continue;

		/* skip SFTs not summed by this unit */
		if((k % msu->divisor)!=msu->remainder)continue;

		/* skip SFTs with low weight */
		if(d->expTMedians[k]<0.05)continue;


		/* power_cor computation */
		weight=d->expTMedians[k]*d->weight;

		demod_weight=weight*response[k];
		msu->total_demod_weight+=demod_weight*response[k];


		f=frequency+frequency*doppler[k]+spindown*(d->gps[k]-spindown_start+d->coherence_time*0.5);
		signal_bin=rintf(1800.0*f-first_bin);
		mismatch=1800.0*f-(signal_bin+first_bin);

// 		if(fabs(spindown)<2e-12)
// 		fprintf(stderr, "%d %g %g f=%f %d %f\n", k, spindown*(d->gps[k]-spindown_start+d->coherence_time*0.5), spindown*(d->gps[k]-spindown_start+d->coherence_time*0.5)*1800.0, f, signal_bin, mismatch);

		fill_hann_filter7(filter, mismatch);

		if((signal_bin<WINDOW+3) || (signal_bin+WINDOW+3>=nbins)) {
			msu->error=-1;
			return;
			}

		bin_re=&(d->re[k*nbins+signal_bin-WINDOW]);
		bin_im=&(d->im[k*nbins+signal_bin-WINDOW]);

		pout=demod_signal_sum;

		for(b=0; b< (2*WINDOW+1); b++) {

			x=bin_re[-3]*filter[0]+bin_re[-2]*filter[1]+bin_re[-1]*filter[2]+bin_re[0]*filter[3]+bin_re[1]*filter[4]+bin_re[2]*filter[5]+bin_re[3]*filter[6];
			y=bin_im[-3]*filter[0]+bin_im[-2]*filter[1]+bin_im[-1]*filter[2]+bin_im[0]*filter[3]+bin_im[1]*filter[4]+bin_im[2]*filter[5]+bin_im[3]*filter[6];

			power=x*x+y*y;

			(*pout)+=(power)*demod_weight;
			pout++;
			bin_re++;
			bin_im++;
			}
		count++;
		if(count>100) {
			for(b=0; b< (2*WINDOW+1); b++) {
				msu->demod_signal_sum_d[b]+=demod_signal_sum[b];
				demod_signal_sum[b]=0.0;
				}
			count=0;
			}
		}
	for(b=0; b< (2*WINDOW+1); b++) {
		msu->demod_signal_sum_d[b]+=demod_signal_sum[b];
		demod_signal_sum[b]=0.0;
		}
	}
}

static void compute_matched_snr(SCORE_AUX_DATA *ad, CANDIDATE *cand, int debug)
{
int i, j, b, b_max;
double total_demod_weight, mean_power, mean_power_sq, power_sd;
float demod_signal_sum[(2*WINDOW+1)*2];
double demod_signal_sum_d[(2*WINDOW+1)*2];
MATCHED_SNR_UNIT *units;
int n_units;

if((ad->valid & (VALID_RESPONSE | VALID_DOPPLER))!=(VALID_RESPONSE | VALID_DOPPLER)) {
	cand->strain=-1;
	cand->snr=-1;
	cand->f_max=cand->frequency;
	fprintf(stderr, "compute_snr: score_aux_data is not valid\n");
	return;
	}


cand->ifo_freq=0.0;
cand->ifo_freq_sd=0.0;

// fprintf(stderr, "(1) spindown =%g\n", spindown);

/* loop over datasets */
n_units=get_max_threads();
units=alloca(n_units*sizeof(*units));

for(i=0;i<n_units;i++) {
	for(j=0;j<(2*WINDOW+1);j++)units[i].demod_signal_sum_d[j]=0.0;
	units[i].total_demod_weight=0;
	units[i].ad=ad;
	units[i].cand=cand;
	units[i].remainder=i;
	units[i].divisor=n_units;
	units[i].error=0;
	
	submit_job(compute_matched_snr_cruncher, &(units[i]));
	}

total_demod_weight=0.0;

for(b=0;b<(2*WINDOW+1);b++) {
	demod_signal_sum_d[b]=0.0;
	}

while(do_single_job(-1));

wait_for_all_done();

for(i=0;i<n_units;i++) {
	if(units[i].error) {
		fprintf(stderr, "Attempting to sample signal outside loaded band, aborting\n");
		cand->snr=-1;
		cand->strain=-1;
		cand->strain_err=1;
		return;
		}
	total_demod_weight+=units[i].total_demod_weight;
	
	for(b=0;b<2*WINDOW+1;b++) {
		demod_signal_sum_d[b]+=units[i].demod_signal_sum_d[b];
		}
	}

total_demod_weight=1.0/total_demod_weight;
for(b=0;b<2*WINDOW+1;b++) {
	demod_signal_sum[b]=demod_signal_sum_d[b]*total_demod_weight;
	}

b_max=WINDOW;
for(b=WINDOW-SEARCH_WINDOW;b<WINDOW+SEARCH_WINDOW+1;b++) {
	if(demod_signal_sum[b_max]<demod_signal_sum[b]) b_max=b;
	}

mean_power=0.0;
mean_power_sq=0.0;
for(b=0;b<WINDOW-SEARCH_WINDOW;b++){
	mean_power+=demod_signal_sum[b];
	mean_power_sq+=demod_signal_sum[b]*demod_signal_sum[b];
	}
for(b=WINDOW+SEARCH_WINDOW+1;b<2*WINDOW+1;b++){
	mean_power+=demod_signal_sum[b];
	mean_power_sq+=demod_signal_sum[b]*demod_signal_sum[b];
	}
	
mean_power*=0.5/(WINDOW-SEARCH_WINDOW);
mean_power_sq*=0.5/(WINDOW-SEARCH_WINDOW);
power_sd=sqrt((mean_power_sq-mean_power*mean_power)*(2.0*(WINDOW-SEARCH_WINDOW)/(2.0*(WINDOW-SEARCH_WINDOW)-1)));

if(debug) {
	fprintf(stderr, "mean_power=%g mean_power_sq=%g power_sd=%g\n", mean_power, mean_power_sq, power_sd);
	fprintf(stderr, "Power: ");
	for(b=0;b<2*WINDOW+1;b++) {
		fprintf(stderr, "%.2f ", (demod_signal_sum[b]-mean_power)/power_sd);
		}
	fprintf(stderr, "\n");
	}

cand->f_max=cand->frequency+(b_max-WINDOW)/1800.0;

b_max=WINDOW;
cand->total_weight=total_demod_weight;
/*cand->snr=(demod_signal_sum[b_max]-mean_power)/power_sd;
cand->snr=(demod_signal_sum[b_max]-0.25*fabs(demod_signal_sum[b_max-1]-demod_signal_sum[b_max+1])-mean_power)/power_sd;*/
SNR_FORMULA

if(demod_signal_sum[b_max]< mean_power)cand->strain=0.0;
	else
	cand->strain=2.0*sqrt(demod_signal_sum[b_max]-mean_power)*args_info.strain_norm_factor_arg/(1800.0*16384.0);

if(cand->strain>0)
	cand->strain_err=0.5*4.0*power_sd*args_info.strain_norm_factor_arg/(cand->strain*(1800.0*1800.0*16384.0*16384.0));
	else
	cand->strain_err=1;

}

void single_pass_compute_matched_snr(CANDIDATE *cand)
{
SCORE_AUX_DATA *ad;

ad=allocate_score_aux_data();

fill_e(ad, cand);
fill_response(ad, cand);
fill_doppler(ad);

compute_matched_snr(ad, cand, 0);

free_score_aux_data(ad);
}

/* This is used to speedup optimization in iota and psi and changes in the frequency path*/

typedef struct {
	CANDIDATE cand; /* candidate data this was computed for */
	double f_plus_sq[(2*WINDOW+1)*2];
	double f_cross_sq[(2*WINDOW+1)*2];
	double f_plus_cross[(2*WINDOW+1)*2];
	int iter_count;
	} SNR_DATA;


static void precompute_snr_data(SNR_DATA *sd, SCORE_AUX_DATA *ad, CANDIDATE *cand, int debug)
{
DATASET *d;
float spindown=cand->spindown;
double frequency=cand->frequency;
int j, k, b;
double weight, f;
float *doppler;
int signal_bin;
float *x, *y, *efm, power, *pout, *cout, *pcout, pweight, cweight, pcweight, f_plus, f_cross;
POLARIZATION *pl;
float f_plus_sq[(2*WINDOW+1)*2];
float f_cross_sq[(2*WINDOW+1)*2];
float f_plus_cross[(2*WINDOW+1)*2];
int count;


memcpy(&(sd->cand), cand, sizeof(*cand));
sd->iter_count=0;
for(b=0;b<(2*WINDOW+1);b++) {
	sd->f_plus_sq[b]=0.0;
	sd->f_cross_sq[b]=0.0;
	sd->f_plus_cross[b]=0.0;
	}

if((ad->valid & VALID_DOPPLER)!=VALID_DOPPLER) {
	fprintf(stderr, "precompute_snr_data: score_aux_data is not valid\n");
	return;
	}

/* loop over datasets */
for(j=0;j<d_free;j++) {
	d=&(datasets[j]);

	doppler=&(ad->doppler[ad->offset[j]]);
	pl=&(d->polarizations[0]);	

	for(b=0;b<(2*WINDOW+1);b++) {
		f_plus_sq[b]=0.0;
		f_cross_sq[b]=0.0;
		f_plus_cross[b]=0.0;
		}
	count=0;

	for(k=0;k<d->free;k++) {

		if(d->sft_veto[k])continue;

		/* skip SFTs with low weight */
		if(d->expTMedians[k]<0.05)continue;
		sd->iter_count++;

		f_plus=F_plus_coeff(k, ad->e, pl->AM_coeffs);
		f_cross=F_plus_coeff(k, ad->e, pl->conjugate->AM_coeffs);

		/* power_cor computation */
		weight=d->expTMedians[k]*d->weight;

		pweight=weight*f_plus*f_plus;
		cweight=weight*f_cross*f_cross;
		pcweight=weight*f_plus*f_cross;

		f=frequency+frequency*doppler[k]+spindown*(d->gps[k]-spindown_start+d->coherence_time*0.5);
		signal_bin=rintf(1800.0*f-first_bin);

		x=&(d->re[k*nbins+signal_bin-WINDOW]);
		y=&(d->im[k*nbins+signal_bin-WINDOW]);
		efm=&(d->expFMedians_plain[signal_bin-WINDOW]);

		pout=f_plus_sq;
		cout=f_cross_sq;
		pcout=f_plus_cross;

		for(b=0; b< (2*WINDOW+1); b++) {
			power=(*x)*(*x)+(*y)*(*y);

			if(args_info.subtract_background_arg) {
				power-=d->expTMedians_plain[k]*(*efm);
				}

			(*pout)+=power*pweight;
			(*cout)+=power*cweight;
			(*pcout)+=power*pcweight;

			pout++;
			cout++;
			pcout++;
			x++;
			y++;
			efm++;
			}
		count++;
		if(count>100) {
			for(b=0;b<(2*WINDOW+1);b++) {
				sd->f_plus_sq[b]+=f_plus_sq[b];
				sd->f_cross_sq[b]+=f_cross_sq[b];
				sd->f_plus_cross[b]+=f_plus_cross[b];
				f_plus_sq[b]=0;
				f_cross_sq[b]=0;
				f_plus_cross[b]=0;
				}
			count=0;
			}
		}
	for(b=0;b<(2*WINDOW+1);b++) {
		sd->f_plus_sq[b]+=f_plus_sq[b];
		sd->f_cross_sq[b]+=f_cross_sq[b];
		sd->f_plus_cross[b]+=f_plus_cross[b];
		}
	}
}

static void update_snr_data(SNR_DATA *sd, SCORE_AUX_DATA *ad, CANDIDATE *cand, int debug)
{
DATASET *d;
float spindown=cand->spindown;
double frequency=cand->frequency;
int j, k, b;
double weight, f, f0;
float *doppler;
int signal_bin, signal_bin0;
float *x, *y, *efm, power, *x0, *y0, *efm0, power0, p, *pout, *cout, *pcout, pweight, cweight, pcweight, f_plus, f_cross;
POLARIZATION *pl;
float f_plus_sq[(2*WINDOW+1)*2];
float f_cross_sq[(2*WINDOW+1)*2];
float f_plus_cross[(2*WINDOW+1)*2];
int count;


if((ad->valid & VALID_DOPPLER)!=VALID_DOPPLER) {
	fprintf(stderr, "update_snr_data: score_aux_data is not valid\n");
	return;
	}

sd->iter_count=0;

/* loop over datasets */
for(j=0;j<d_free;j++) {
	d=&(datasets[j]);

	doppler=&(ad->doppler[ad->offset[j]]);
	pl=&(d->polarizations[0]);	

	for(b=0;b<(2*WINDOW+1);b++) {
		f_plus_sq[b]=0.0;
		f_cross_sq[b]=0.0;
		f_plus_cross[b]=0.0;
		}
	count=0;

	for(k=0;k<d->free;k++) {

		if(d->sft_veto[k])continue;

		/* skip SFTs with low weight */
		if(d->expTMedians[k]<0.05)continue;

		f=frequency+frequency*doppler[k]+spindown*(d->gps[k]-spindown_start+d->coherence_time*0.5);
		signal_bin=rintf(1800.0*f-first_bin);

		f0=sd->cand.frequency+sd->cand.frequency*doppler[k]+sd->cand.spindown*(d->gps[k]-spindown_start+d->coherence_time*0.5);
		signal_bin0=rintf(1800.0*f0-first_bin);

		if(signal_bin==signal_bin0)continue;

		sd->iter_count++;

		f_plus=F_plus_coeff(k, ad->e, pl->AM_coeffs);
		f_cross=F_plus_coeff(k, ad->e, pl->conjugate->AM_coeffs);

		/* power_cor computation */
		weight=d->expTMedians[k]*d->weight;

		pweight=weight*f_plus*f_plus;
		cweight=weight*f_cross*f_cross;
		pcweight=weight*f_plus*f_cross;


		x=&(d->re[k*nbins+signal_bin-WINDOW]);
		y=&(d->im[k*nbins+signal_bin-WINDOW]);
		efm=&(d->expFMedians_plain[signal_bin-WINDOW]);
		x0=&(d->re[k*nbins+signal_bin0-WINDOW]);
		y0=&(d->im[k*nbins+signal_bin0-WINDOW]);
		efm0=&(d->expFMedians_plain[signal_bin0-WINDOW]);

		pout=f_plus_sq;
		cout=f_cross_sq;
		pcout=f_plus_cross;

		for(b=0; b< (2*WINDOW+1); b++) {
			power=(*x)*(*x)+(*y)*(*y);
			power0=(*x0)*(*x0)+(*y0)*(*y0);

			if(args_info.subtract_background_arg) {
				power-=d->expTMedians_plain[k]*(*efm);
				power0-=d->expTMedians_plain[k]*(*efm0);
				}

			p=(power)-(power0);
			(*pout)+=p*pweight;
			(*cout)+=p*cweight;
			(*pcout)+=p*pcweight;

			pout++;
			cout++;
			pcout++;
			x++;
			y++;
			efm++;
			x0++;
			y0++;
			efm0++;
			}
		count++;
		if(count>100) {
			for(b=0;b<(2*WINDOW+1);b++) {
				sd->f_plus_sq[b]+=f_plus_sq[b];
				sd->f_cross_sq[b]+=f_cross_sq[b];
				sd->f_plus_cross[b]+=f_plus_cross[b];
				f_plus_sq[b]=0;
				f_cross_sq[b]=0;
				f_plus_cross[b]=0;
				}
			count=0;
			}
		}
	for(b=0;b<(2*WINDOW+1);b++) {
		sd->f_plus_sq[b]+=f_plus_sq[b];
		sd->f_cross_sq[b]+=f_cross_sq[b];
		sd->f_plus_cross[b]+=f_plus_cross[b];
		}
	}
memcpy(&(sd->cand), cand, sizeof(*cand));
}

static void matched_precompute_snr_data(SNR_DATA *sd, SCORE_AUX_DATA *ad, CANDIDATE *cand, int debug)
{
DATASET *d;
float spindown=cand->spindown;
double frequency=cand->frequency;
int j, k, b;
double weight, f;
float *doppler;
int signal_bin;
float power, *pout, *cout, *pcout, pweight, cweight, pcweight, f_plus, f_cross, x, y;
POLARIZATION *pl;
float *bin_re, *bin_im;
float mismatch;
float filter[7];
float f_plus_sq[(2*WINDOW+1)*2];
float f_cross_sq[(2*WINDOW+1)*2];
float f_plus_cross[(2*WINDOW+1)*2];
int count;

memcpy(&(sd->cand), cand, sizeof(*cand));
sd->iter_count=0;
for(b=0;b<(2*WINDOW+1);b++) {
	sd->f_plus_sq[b]=0.0;
	sd->f_cross_sq[b]=0.0;
	sd->f_plus_cross[b]=0.0;
	}

if((ad->valid & VALID_DOPPLER)!=VALID_DOPPLER) {
	fprintf(stderr, "matched_precompute_snr_data: score_aux_data is not valid\n");
	return;
	}

/* loop over datasets */
for(j=0;j<d_free;j++) {
	d=&(datasets[j]);

	doppler=&(ad->doppler[ad->offset[j]]);
	pl=&(d->polarizations[0]);	

	for(b=0;b<(2*WINDOW+1);b++) {
		f_plus_sq[b]=0.0;
		f_cross_sq[b]=0.0;
		f_plus_cross[b]=0.0;
		}
	count=0;

	for(k=0;k<d->free;k++) {

		if(d->sft_veto[k])continue;

		/* skip SFTs with low weight */
		if(d->expTMedians[k]<0.05)continue;
		sd->iter_count++;

		f_plus=F_plus_coeff(k, ad->e, pl->AM_coeffs);
		f_cross=F_plus_coeff(k, ad->e, pl->conjugate->AM_coeffs);

		/* power_cor computation */
		weight=d->expTMedians[k]*d->weight;

		pweight=weight*f_plus*f_plus;
		cweight=weight*f_cross*f_cross;
		pcweight=weight*f_plus*f_cross;

		f=frequency+frequency*doppler[k]+spindown*(d->gps[k]-spindown_start+d->coherence_time*0.5);
		signal_bin=rintf(1800.0*f-first_bin);
		mismatch=1800.0*f-signal_bin-first_bin;

		fill_hann_filter7(filter, mismatch);

		bin_re=&(d->re[k*nbins+signal_bin-WINDOW]);
		bin_im=&(d->im[k*nbins+signal_bin-WINDOW]);

		pout=f_plus_sq;
		cout=f_cross_sq;
		pcout=f_plus_cross;

		for(b=0; b< (2*WINDOW+1); b++) {

			/*  Tcl code to generate these lines:
			for { set i -3 } { $i < 4 } { incr i } { puts -nonewline "bin\[$i\].re*filter\[[expr $i+3]\]+" }
			*/

			x=bin_re[-3]*filter[0]+bin_re[-2]*filter[1]+bin_re[-1]*filter[2]+bin_re[0]*filter[3]+bin_re[1]*filter[4]+bin_re[2]*filter[5]+bin_re[3]*filter[6];
			y=bin_im[-3]*filter[0]+bin_im[-2]*filter[1]+bin_im[-1]*filter[2]+bin_im[0]*filter[3]+bin_im[1]*filter[4]+bin_im[2]*filter[5]+bin_im[3]*filter[6];

			power=x*x+y*y;

			(*pout)+=(power)*pweight;
			(*cout)+=(power)*cweight;
			(*pcout)+=(power)*pcweight;

			pout++;
			cout++;
			pcout++;
			bin_re++;
			bin_im++;
			}
		count++;
		if(count>100) {
			for(b=0;b<(2*WINDOW+1);b++) {
				sd->f_plus_sq[b]+=f_plus_sq[b];
				sd->f_cross_sq[b]+=f_cross_sq[b];
				sd->f_plus_cross[b]+=f_plus_cross[b];
				f_plus_sq[b]=0;
				f_cross_sq[b]=0;
				f_plus_cross[b]=0;
				}
			count=0;
			}
		}
	for(b=0;b<(2*WINDOW+1);b++) {
		sd->f_plus_sq[b]+=f_plus_sq[b];
		sd->f_cross_sq[b]+=f_cross_sq[b];
		sd->f_plus_cross[b]+=f_plus_cross[b];
		}
	}
}

void compute_alignment_snr(SNR_DATA *sd, CANDIDATE *cand, int debug)
{
int b, b_max;
float demod_signal_sum[(2*WINDOW+1)];
double mean_power, mean_power_sq, power_sd;
float ac, ap, apc, ci;

ci=cos(cand->iota);
ci=ci*ci;

ap=(1+ci)*0.5;
ap=ap*ap;

ac=ap-ci;
ap=ap+ci;

apc=ac;
ac=ac*cos(4*cand->psi);

apc=2.0*apc*sin(4*cand->psi);

// fprintf(stderr, "ap=%f ac=%f apc=%f\n", ap, ac, apc);

for(b=0;b<(2*WINDOW+1);b++)
	demod_signal_sum[b]=(sd->f_plus_sq[b]+sd->f_cross_sq[b])*ap+(sd->f_plus_sq[b]-sd->f_cross_sq[b])*ac+sd->f_plus_cross[b]*apc;

b_max=WINDOW;
for(b=WINDOW-SEARCH_WINDOW;b<WINDOW+SEARCH_WINDOW+1;b++) {
	if(demod_signal_sum[b_max]<demod_signal_sum[b]) b_max=b;
	}

mean_power=0.0;
mean_power_sq=0.0;

for(b=0;b< WINDOW-SEARCH_WINDOW;b++){
	mean_power+=demod_signal_sum[b];
	mean_power_sq+=demod_signal_sum[b]*demod_signal_sum[b];
	}
for(b=WINDOW+SEARCH_WINDOW+1;b<2*WINDOW+1;b++){
	mean_power+=demod_signal_sum[b];
	mean_power_sq+=demod_signal_sum[b]*demod_signal_sum[b];
	}
	
mean_power*=0.5/((WINDOW-SEARCH_WINDOW));
mean_power_sq*=0.5/((WINDOW-SEARCH_WINDOW));
power_sd=sqrt((mean_power_sq-mean_power*mean_power)*(2.0*(WINDOW-SEARCH_WINDOW)/(2.0*(WINDOW-SEARCH_WINDOW)-1)));

if(debug) {
	fprintf(stderr, "mean_power=%g mean_power_sq=%g power_sd=%g\n", mean_power, mean_power_sq, power_sd);
	fprintf(stderr, "Power: ");
	for(b=0;b<2*WINDOW+1;b++) {
		fprintf(stderr, "%.2f ", (demod_signal_sum[b]-mean_power)/power_sd);
		}
	fprintf(stderr, "\n");
	}

cand->f_max=cand->frequency+(b_max-WINDOW)/1800.0;

b_max=WINDOW;
SNR_FORMULA
}

static void test_alignment_snr(CANDIDATE *cand)
{
CANDIDATE c1, c2;
SCORE_AUX_DATA *ad;
SNR_DATA sd;
int i,j, N=16;
double err_max, err;

memcpy(&c1, cand, sizeof(*cand));
memcpy(&c2, cand, sizeof(*cand));

ad=allocate_score_aux_data();

fill_e(ad, &c1);
fill_response(ad, &c1);
fill_doppler(ad);

precompute_snr_data(&sd, ad, &c1, 0);

err_max=0;

for(i=-N; i<=N; i++) {
	for(j=-N; j<=N; j++) {
		c1.psi=cand->psi+j*M_PI/128.0;
		c1.iota=cand->iota+i*M_PI/128.0;

		//c1.psi=M_PI/2.0;
		//c1.iota=M_PI/2.0;

		c2.psi=c1.psi;
		c2.iota=c1.iota;

		fill_response(ad, &c1);
		compute_simple_snr(ad, &c1, 0);
		compute_alignment_snr(&sd, &c2, 0);
		err=c1.snr-c2.snr;
		if(fabs(err)>(1e-3)*(fabs(c1.snr)+fabs(c2.snr)))
			fprintf(stderr, "iota=%f psi=%f c1.snr=%f c2.snr=%f snr_diff=%f\n", c1.iota, c1.psi, c1.snr, c2.snr, err);
		if(err>err_max)err_max=err;
		}
	}

fprintf(stderr, "candidate_debug rank=%d alignment error %f\n", cand->rank, err_max);
fprintf(LOG, "candidate_debug rank=%d alignment error %f\n", cand->rank, err_max);

free_score_aux_data(ad);
}


#define FREQ_STEPS 1

inline static void compute_snr(SCORE_AUX_DATA *ad, CANDIDATE *cand, int debug)
{
DATASET *d;
float spindown=cand->spindown;
double frequency=cand->frequency;
int i, j, k, b, b_max;
double weight, f, demod_weight, total_demod_weight, mean_power, mean_power_sq, power_sd;
float *response;
float *doppler;
int signal_bin, offset;
float demod_signal_sum[(2*WINDOW+1)*FREQ_STEPS], *x, *y, *efm, power, *pout;

if((ad->valid & (VALID_RESPONSE | VALID_DOPPLER))!=(VALID_RESPONSE | VALID_DOPPLER)) {
	cand->strain=-1;
	cand->snr=-1;
	cand->f_max=cand->frequency;
	fprintf(stderr, "compute_snr: score_aux_data is not valid\n");
	return;
	}

total_demod_weight=0.0;

for(b=0;b<(2*WINDOW+1)*FREQ_STEPS;b++) {
	demod_signal_sum[b]=0.0;
	}


cand->ifo_freq=0.0;
cand->ifo_freq_sd=0.0;

/* loop over datasets */
for(j=0;j<d_free;j++) {
	d=&(datasets[j]);

	response=&(ad->response[ad->offset[j]]);
	doppler=&(ad->doppler[ad->offset[j]]);

	for(k=0;k<d->free;k++) {

		if(d->sft_veto[k])continue;

		/* skip SFTs with low weight */
		if(d->expTMedians[k]<0.05)continue;

		/* power_cor computation */
		weight=d->expTMedians[k]*d->weight;

		demod_weight=weight*response[k];
		total_demod_weight+=demod_weight*response[k];

		for(i=0;i<FREQ_STEPS;i++) {
			f=frequency+frequency*doppler[k]+spindown*(d->gps[k]-spindown_start+d->coherence_time*0.5)+i*(1.0/(1800.0*FREQ_STEPS));
			signal_bin=rintf(1800.0*f-first_bin);
			offset=i*(2*WINDOW+1);
	
			x=&(d->re[k*nbins+signal_bin-WINDOW]);
			y=&(d->im[k*nbins+signal_bin-WINDOW]);
			efm=&(d->expFMedians_plain[signal_bin-WINDOW]);
			pout=&(demod_signal_sum[offset]);
	
			for(b=0; b< (2*WINDOW+1); b++) {
				power=(*x)*(*x)+(*y)*(*y);

				if(args_info.subtract_background_arg) {
					power-=d->expTMedians_plain[k]*(*efm);
					}
	
				(*pout)+=power*demod_weight;
				x++;
				y++;
				efm++;
				pout++;
				}
			}
		}
	}

total_demod_weight=1.0/total_demod_weight;
for(b=0;b< (2*WINDOW+1)*FREQ_STEPS;b++) {
	demod_signal_sum[b]*=total_demod_weight;
	}

b_max=WINDOW;
for(b=WINDOW-SEARCH_WINDOW;b<WINDOW+SEARCH_WINDOW+1;b++) {
	if(demod_signal_sum[b_max]<demod_signal_sum[b]) b_max=b;
	}

mean_power=0.0;
mean_power_sq=0.0;
for(i=0;i<FREQ_STEPS;i++) {
	offset=i*(2*WINDOW+1);
	for(b=offset;b< offset+WINDOW-SEARCH_WINDOW;b++){
		mean_power+=demod_signal_sum[b];
		mean_power_sq+=demod_signal_sum[b]*demod_signal_sum[b];
		}
	for(b=offset+WINDOW+SEARCH_WINDOW+1;b<offset+2*WINDOW+1;b++){
		mean_power+=demod_signal_sum[b];
		mean_power_sq+=demod_signal_sum[b]*demod_signal_sum[b];
		}
	}
	
mean_power*=0.5/((WINDOW-SEARCH_WINDOW)*FREQ_STEPS);
mean_power_sq*=0.5/((WINDOW-SEARCH_WINDOW)*FREQ_STEPS);
power_sd=sqrt((mean_power_sq-mean_power*mean_power)*(2.0*FREQ_STEPS*(WINDOW-SEARCH_WINDOW)/(2.0*FREQ_STEPS*(WINDOW-SEARCH_WINDOW)-1)));

if(debug) {
	fprintf(stderr, "mean_power=%g mean_power_sq=%g power_sd=%g\n", mean_power, mean_power_sq, power_sd);
	fprintf(stderr, "Power: ");
	for(b=0;b<2*WINDOW+1;b++) {
		fprintf(stderr, "%.2f ", (demod_signal_sum[b]-mean_power)/power_sd);
		}
	fprintf(stderr, "\n");
	}

cand->f_max=cand->frequency+(b_max-WINDOW)/1800.0;

b_max=WINDOW;
cand->total_weight=total_demod_weight;
/*cand->snr=(demod_signal_sum[b_max]-mean_power)/power_sd;
cand->snr=(demod_signal_sum[b_max]-0.25*fabs(demod_signal_sum[b_max-1]-demod_signal_sum[b_max+1])-mean_power)/power_sd;*/
SNR_FORMULA
if(demod_signal_sum[b_max]< mean_power)cand->strain=0.0;
	else
	cand->strain=2.0*sqrt(demod_signal_sum[b_max]-mean_power)*args_info.strain_norm_factor_arg/(1800.0*16384.0);

if(cand->strain>0)
	cand->strain_err=0.5*4.0*power_sd*args_info.strain_norm_factor_arg/(cand->strain*(1800.0*1800.0*16384.0*16384.0));
	else
	cand->strain_err=1;


}

typedef struct {
	float frequency_shift;
	float spindown_shift;
	float ra_shift;
	float dec_shift;
	double max_power;
	} HINTS;

#define FREQ_STEPS	32

inline static void compute_hints(SCORE_AUX_DATA *ad, CANDIDATE *cand, HINTS *h)
{
DATASET *d;
float spindown=cand->spindown;
double frequency=cand->frequency;
int j, k, i, b;
float mismatch;
double weight, f, x, y, demod_weight, power[3], freq_hint_right[FREQ_STEPS], freq_hint_left[FREQ_STEPS], max_hint;
float *response;
float *doppler;
int signal_bin;
float *bin_re, *bin_im;

if((ad->valid & (VALID_RESPONSE | VALID_DOPPLER))!=(VALID_RESPONSE | VALID_DOPPLER)) {
	cand->strain=-1;
	cand->snr=-1;
	cand->f_max=cand->frequency;
	fprintf(stderr, "compute_snr: score_aux_data is not valid\n");
	return;
	}

for(k=0;k<FREQ_STEPS;k++) {
	freq_hint_right[k]=0.0;
	freq_hint_left[k]=0.0;
	}

/* loop over datasets */
for(j=0;j<d_free;j++) {
	d=&(datasets[j]);

	response=&(ad->response[ad->offset[j]]);
	doppler=&(ad->doppler[ad->offset[j]]);

	for(k=0;k<d->free;k++) {

		if(d->sft_veto[k])continue;

		/* skip SFTs with low weight */
		if(d->expTMedians[k]<0.05)continue;

		/* power_cor computation */
		weight=d->expTMedians[k]*d->weight;

		demod_weight=weight*response[k];

		f=frequency+frequency*doppler[k]+spindown*(d->gps[k]-spindown_start+d->coherence_time*0.5);
		signal_bin=rint(1800.0*f-first_bin);
		mismatch=1800.0*f-first_bin-signal_bin;

		bin_re=&(d->re[k*nbins+signal_bin-1]);
		bin_im=&(d->im[k*nbins+signal_bin-1]);

		for(b=0;b<=2; b++) {
			x=bin_re[b];
			y=bin_im[b];
			power[b]=(x*x+y*y)*demod_weight;
			}
		for(i=0;i<FREQ_STEPS;i++) {
			if(mismatch*FREQ_STEPS*2.0>(1.0*FREQ_STEPS-i-0.5))
				freq_hint_right[i]+=(power[2]-power[1]);
			if(mismatch*FREQ_STEPS*2.0< (i-1.0*FREQ_STEPS+0.5))
				freq_hint_left[i]+=(power[0]-power[1]);
			}
		}
	}
max_hint=0.0;
h->frequency_shift=0;

for(i=0;i<FREQ_STEPS;i++) {
	if(freq_hint_right[i]>max_hint) {
		h->frequency_shift=((i+0.5)*0.5)/FREQ_STEPS;
		}
	if(freq_hint_left[i]>max_hint) {
		h->frequency_shift=-((i+0.5)*0.5)/FREQ_STEPS;
		}
	}
h->max_power=max_hint;
h->frequency_shift/=1800.0;
}

/* candidate cache */

VARRAY *candidate_cache=NULL;
VARRAY **bin_index=NULL;
double index_queries_total=0.0;
double index_snr_total=0.0;
double index_hits=0.0;
double improvement_queries_total=0.0;
double improvements=0.0;
double improvement_snr_total=0.0;

void init_candidates(void)
{
int i;
candidate_cache=new_varray(sizeof(CANDIDATE));
bin_index=do_alloc(nbins, sizeof(*bin_index));
for(i=0;i<nbins;i++) {
	bin_index[i]=new_varray(sizeof(int));
	}
}

int find_next_candidate(CANDIDATE *cand)
{
int i;
int k;
int best;
double best_snr;
VARRAY *v;
CANDIDATE *c;
double timebase=max_gps()-min_gps();

index_queries_total++;

k=round(cand->frequency*1800.0 -first_bin);
if(k<0)k=0;
if(k>=nbins)k=nbins-1;
v=bin_index[k];

best_snr=cand->snr;
best=-1;
for(i=v->free-1;i>=0;i--) {
	k=VELT(v, int, i);
	c=&(VELT(candidate_cache, CANDIDATE, k));

	if(c->snr <= best_snr)continue;
	if(fabs(cand->spindown-c->spindown)*timebase*1800.0>3)continue;
	if(candidate_distance(c, cand->ra, cand->dec)>5.5*resolution)continue;

	best=k;
	best_snr=c->snr;
	}
if(best>0) {
	index_hits++;
	index_snr_total+=best_snr-cand->snr;
	}
return best;
}

int find_better_candidate(int index)
{
int k;
float start_snr;
improvement_queries_total++;

if(index<0)return index;

start_snr=VELT(candidate_cache, CANDIDATE, index).snr;

while( (k=VELT(candidate_cache, CANDIDATE, index).better_candidate)>=0)index=k;

if(index>0) {
	improvements++;
	improvement_snr_total+=VELT(candidate_cache, CANDIDATE, index).snr - start_snr;
	}
return index;
}

static CANDIDATE *get_candidate(int index)
{
if(index<0)return NULL;
if(index>=candidate_cache->free)return NULL;
return &(VELT(candidate_cache, CANDIDATE, index));
}

int add_candidate(CANDIDATE *cand, int improves)
{
int i, k;
int bin;

k=varray_add(candidate_cache, cand);
VELT(candidate_cache, CANDIDATE, k).better_candidate=-1;

if(improves>=0)VELT(candidate_cache, CANDIDATE, improves).better_candidate=k;

bin=round(cand->frequency*1800.0-first_bin);
for(i=bin-3;i<=bin+3;i++) {
	if(i<0)continue;
	if(i>=nbins)continue;
	varray_add(bin_index[i], &k);
	}
return(k);
}

#define SNR_TOLERANCE 0.01

#define BETTER_SNR(c1, c2)	(((c1).snr>(c2).snr) || (((c1).snr>=(c2).snr)) && ((c1).strain>(c2).strain))
#define BETTER_POWER_COR(c1, c2)	(((c1).snr>=(c2).snr) && ((c1).power_cor>(c2).power_cor))


#define BETTER_SNR_PC(c1, c2)	(((c1).snr+(c1).power_cor)>((c2).snr+(c2).power_cor))

#define BETTER_SNR_PC(c1, c2)	(((c1).snr>(c2).snr+0.1) || (((c1).snr>=(c2).snr) && ((c1).power_cor>(c2).power_cor)))

//#define BETTER_SNR_PC(c1, c2)	(((c1).snr>(c2).snr) || (((c1).snr>=(c2).snr) && ((c1).coherence_score>(c2).coherence_score)))

//#define BETTER_SNR_PC(c1, c2)	(((c1).snr>(c2).snr+0.1) || (((c1).snr>=(c2).snr) && ((c1).strain>(c2).strain)))

#define BETTER_SNR_COH(c1, c2)	(((c1).snr>(c2).snr+0.1) || (((c1).snr>=(c2).snr)) && ((c1).coherence_score>(c2).coherence_score))

#define PLAIN_SNR(c1, c2)	((c1).snr>(c2).snr)

#define CHASE1(opt_cond, opt_var, first_step) \
	int chase_##opt_var(CANDIDATE *cand)  \
	{  \
	CANDIDATE c;  \
	SCORE_AUX_DATA *ad; \
	float fs; \
	float step;  \
	int status=0;  \
	int improv=0; \
	  \
	fs=first_step; \
	step=fs; \
	  \
	ad=allocate_score_aux_data(); \
	\
	memcpy(&c, cand, sizeof(*cand));  \
	\
	fill_e(ad, cand); \
	fill_response(ad, cand); \
	fill_doppler(ad); \
	\
	compute_snr(ad, cand, 0); \
	  \
	while(fabs(step)>fs/128.0) {  \
		c.frequency=cand->f_max; \
		c.opt_var+=step;  \
		\
		fill_e(ad, &c); \
		fill_response(ad, &c); \
		fill_doppler(ad); \
		\
		compute_snr(ad, &c, 0); \
		\
		if(opt_cond(c, *cand)) {  \
			/* c.frequency=c.f_max; */ \
			memcpy(cand, &c, sizeof(CANDIDATE));  \
			status=1;  \
			improv=1; \
			continue;  \
			} \
		c.frequency=cand->f_max+0.25/1800.0; \
		compute_snr(ad, &c, 0); \
		if(opt_cond(c, *cand)) {  \
			/* c.frequency=c.f_max; */ \
			memcpy(cand, &c, sizeof(CANDIDATE));  \
			status=1;  \
			improv=1; \
			continue;  \
			} \
		c.frequency=cand->f_max+0.5/1800.0; \
		compute_snr(ad, &c, 0); \
		if(opt_cond(c, *cand)) {  \
			/* c.frequency=c.f_max; */ \
			memcpy(cand, &c, sizeof(CANDIDATE));  \
			status=1;  \
			improv=1; \
			continue;  \
			} \
		c.frequency=cand->f_max+0.75/1800.0; \
		compute_snr(ad, &c, 0); \
		if(opt_cond(c, *cand)) {  \
			/* c.frequency=c.f_max; */ \
			memcpy(cand, &c, sizeof(CANDIDATE));  \
			status=1;  \
			improv=1; \
			continue;  \
			} \
		c.opt_var-=step;  \
		if(status==0) {  \
			step=-step;  \
			}  \
		status=0;  \
		step=step/2.0;  \
		}  \
	free_score_aux_data(ad); \
	if(improv) {  \
		return 1;  \
		} else {  \
		return 0;  \
		}  \
	}

#define CHASE(opt_cond, opt_var, first_step) \
	int chase_##opt_var(CANDIDATE *cand)  \
	{  \
	CANDIDATE c;  \
	SCORE_AUX_DATA *ad; \
	float fs; \
	float step;  \
	int status=0;  \
	int improv=0; \
	  \
	fs=first_step; \
	step=fs; \
	  \
	ad=allocate_score_aux_data(); \
	\
	memcpy(&c, cand, sizeof(*cand));  \
	\
	fill_e(ad, cand); \
	fill_response(ad, cand); \
	fill_doppler(ad); \
	\
	compute_snr(ad, cand, 0); \
	  \
	while(fabs(step)>fs*1e-10) {  \
		/* c.frequency=cand->f_max; */ \
		c.opt_var+=step;  \
		\
		fill_e(ad, &c); \
		fill_response(ad, &c); \
		fill_doppler(ad); \
		\
		compute_matched_snr(ad, &c, 0); \
		\
		if(opt_cond(c, *cand)) {  \
			/* c.frequency=c.f_max; */ \
			memcpy(cand, &c, sizeof(CANDIDATE));  \
			status=1;  \
			improv=1; \
			continue;  \
			} \
		c.opt_var-=step;  \
		if(status==0) {  \
			step=-step;  \
			}  \
		status=0;  \
		step=step/2.0;  \
		}  \
	free_score_aux_data(ad); \
	if(improv) {  \
		return 1;  \
		} else {  \
		return 0;  \
		}  \
	}

#define FIT(opt_expr, opt_var) \
	int fit_##opt_var(CANDIDATE *cand, float step0) \
	{ \
	CANDIDATE c; \
	int i, max_i, N=128, count;\
	float a[3], alpha, beta, s[2], f, max; \
	float step=step0; \
	\
	start: \
	memcpy(&c, cand, sizeof(*cand));  \
	\
	a[0]=0; \
	a[1]=0; \
	a[2]=0; \
	s[0]=0; \
	s[1]=0; \
	max_i=-N-1; \
	count=0; \
	fprintf(stderr, "" # opt_var "=%f step=%f ", c.opt_var, step); \
	for(i=-N;i<=N;i++) { \
		c.opt_var=cand->opt_var+i*step; \
		compute_scores(&c, 0); \
		f=opt_expr; \
		fprintf(stderr, ",%f", f); \
		if(max_i<-N || (f>max)) { \
			max_i=i; \
			max=f; \
			} \
		a[0]+=f; \
		a[1]+=f*i*step; \
		a[2]+=f*i*i*step*step; \
		f=i*i*step*step; \
		s[0]+=f; \
		s[1]+=f*f; \
		count++; \
		} \
	fprintf(stderr, "\n"); \
	alpha=a[1]/s[0]; \
	beta=(a[2]-a[0]*s[0]/count)/(s[1]-s[0]*s[0]/(count*count)); \
	f=-0.5*alpha/beta; \
	fprintf(stderr, "a=(%f,%f,%f) s=(%f, %f) alpha=%f beta=%f max=%f max_offset=%f offset=%f\n", a[0], a[1], a[2], s[0], s[1], alpha, beta, max, 1.0*max_i*step, f); \
	if(fabs(beta*s[0])<fabs(a[0])*0.05*count) { \
		if(fabs(alpha)<fabs(s[0]*0.05)) return 0; \
		if(alpha*max_i>0) { \
			cand->opt_var+=max_i*step; \
			return 1; \
			} \
		return 0; \
		} \
	if(beta>0) { \
		if(f*max_i<0) { \
			cand->opt_var+=max_i*step; \
			fprintf(stderr, "branch 1\n"); \
			return max_i!=0; \
			} \
		cand->opt_var+= ((f<0)*2-1)*N*step*0.5; \
		fprintf(stderr, "branch 1a\n"); \
		return 1; \
		step=0.5*step; \
		goto start; \
		} \
/*	if(f*max_i>0) { \
		cand->opt_var+=max_i*step; \
		fprintf(stderr, "branch 0\n"); \
		return 1; \
		} \
	if((fabs(N*step)<fabs(f))) { \
		if(max_i) { \
			cand->opt_var+=max_i*step; \
			fprintf(stderr, "branch 2\n"); \
			return 1;\
			} \
		return 0; \
		if(fabs(f)> N*step*0.5) { \
			step=0.5*step; \
			if(step< (step0/256)) return 0; \
			goto start; \
			} \
		} \*/ \
	if(f>N*step)f=N*step; \
	if(f<-N*step)f=-N*step; \
	cand->opt_var+=f; \
	fprintf(stderr, "branch 3\n"); \
	return (fabs(f)>step/4); \
	}

#define FIT2(opt_expr, opt_var) \
	int fit_##opt_var(CANDIDATE *cand, float step0) \
	{ \
	CANDIDATE c; \
	int i, max_i, N=64;\
	float a, alpha, beta, s, f, max; \
	float step=step0; \
	\
	start: \
	memcpy(&c, cand, sizeof(*cand));  \
	\
	a=0; \
	s=0; \
	max_i=0; \
	fprintf(stderr, "" # opt_var "=%f step=%f ", c.opt_var, step); \
	for(i=-N;i<=N;i++) { \
		c.opt_var=cand->opt_var+i*step; \
		compute_scores(&c, 0); \
		f=opt_expr; \
		fprintf(stderr, " %f", f); \
		if(max_i<-N || (f>max)) { \
			max_i=i; \
			max=f; \
			} \
		a+=f*i; \
		s+=f; \
		} \
	fprintf(stderr, "\n"); \
	alpha=a/s; \
	\
	c.opt_var=cand->opt_var+alpha*step; \
	compute_scores(&c, 0); \
	f=opt_expr; \
	if(max> 3*f) {\
		return 0; \
		} \
	if(fabs(alpha)<0.1)return 0; \
	cand->opt_var+=alpha*step; \
	return 1; \
	}

#define SEARCH(opt_expr, opt_var) \
	int search_##opt_var(CANDIDATE *cand, float step0, int N) \
	{ \
	CANDIDATE c; \
	int i, max_i;\
	float f, max; \
	float step=step0; \
	SCORE_AUX_DATA *ad; \
	\
	ad=allocate_score_aux_data(); \
	\
	start: \
	memcpy(&c, cand, sizeof(*cand));  \
	\
	max_i=-N-1; \
	max=0; \
	fprintf(stderr, "" # opt_var "=%f step=%f ", c.opt_var, step);  \
	for(i=-N;i<=N;i++) { \
		c.opt_var=cand->opt_var+i*step; \
		fill_e(ad, &c); \
		fill_response(ad, &c); \
		fill_doppler(ad); \
		\
		compute_matched_snr(ad, &c, (i==-N)); \
		f=opt_expr; \
		fprintf(stderr, ",%f", f); \
		if(max_i<-N || (f>max) || (f>=max && !max_i)) { \
			max_i=i; \
			max=f; \
			} \
		} \
	free_score_aux_data(ad); \
	fprintf(stderr, "\n"); \
	if(max<cand->snr+0.001 && cand->snr>0)return 0; \
	cand->opt_var+=max_i*step; \
	cand->snr=max; \
	return max_i!=0; \
	}

#define ZOOMED_SEARCH(opt_expr, opt_var) \
	int zoomed_search_##opt_var(CANDIDATE *cand, float step0, int N) \
	{ \
	while(step0>0) { \
		if(search_##opt_var(cand, step0, N))return 1; \
		step0=step0/4; \
		} \
	return 0; \
	}


#define RECURSIVE_SEARCH(opt_expr, opt_var) \
	int recursive_search_##opt_var(CANDIDATE *cand, float step0, int N) \
	{ \
	CANDIDATE c; \
	int i, max_i;\
	float f, max; \
	float step=step0, step_max; \
	float *values; \
	\
	values=alloca((2*N+1) *sizeof(*values)); \
	\
	memcpy(&c, cand, sizeof(*cand));  \
	\
	/* first pass - sample */ \
	max_i=-N-1; \
	max=0.0; \
	fprintf(stderr, "" # opt_var "=%f step=%f\n", c.opt_var, step); \
	for(i=-N;i<=N;i++) { \
		c.opt_var=cand->opt_var+i*step; \
		compute_scores(&c, 0); \
		f=opt_expr; \
		values[i+N]=f; \
		if(max_i<-N || (f>max) || (f>=max && !max_i)) { \
			max_i=i; \
			max=f; \
			} \
		} \
	step_max=max_i*step; \
	for(i=2;i<2*N-1;i++) { \
		if((values[i]>values[i-1]) && (values[i]>values[i+1]) && (values[i]>values[i-2]) && (values[i]>values[i+2])) { \
			fprintf(stderr, " %d %f\n", i, values[i]); \
			c.opt_var=cand->opt_var+(i-N)*step; \
			if(!search_##opt_var(&c, step0/16, 32))continue; \
			f=opt_expr; \
			if(f>max) { \
				max=f; \
				 step_max=c.opt_var -cand->opt_var; \
				} \
			} \
		} \
	cand->opt_var+=step_max; \
	cand->snr=max; \
	fprintf(stderr, "step_max=%g\n", step_max); \
	return fabs(step_max)>0; \
	}

/* this is useful for forcing known values for tinkering with the algorithm */
#define ORACLE(opt_char, opt_var, value, force) \
	int chase_##opt_var(CANDIDATE *cand)  \
	{  \
	CANDIDATE c;  \
	  \
	memcpy(&c, cand, sizeof(*cand));  \
	  \
	c.opt_var=(value);  \
	compute_scores(&c, 0);  \
	\
	if((force) || (c.opt_char>cand->opt_char)) {  \
		memcpy(cand, &c, sizeof(CANDIDATE));  \
		return (c.opt_char>cand->opt_char);  \
		} else {  \
		return 0; \
		} \
	}

CHASE(PLAIN_SNR, psi, M_PI/32.0)
CHASE(PLAIN_SNR, iota, M_PI/256.0)
CHASE(PLAIN_SNR, ra, resolution/128.0)
CHASE(PLAIN_SNR, dec, resolution/128.0)
CHASE(PLAIN_SNR, spindown, 0.01/(1800.0*(max_gps()-min_gps())))
//CHASE(PLAIN_SNR, frequency, 1/3600.0)


FIT(c.snr, psi)
FIT(c.snr, iota)
FIT(c.snr, ra)
FIT(c.snr, dec)
// 
SEARCH(c.snr, frequency)
SEARCH(c.snr, psi)
SEARCH(c.snr, iota)
SEARCH(c.snr, ra)
SEARCH(c.snr, dec)
SEARCH(c.snr, spindown)

ZOOMED_SEARCH(c.snr, frequency)
ZOOMED_SEARCH(c.snr, psi)
ZOOMED_SEARCH(c.snr, iota)
ZOOMED_SEARCH(c.snr, ra)
ZOOMED_SEARCH(c.snr, dec)
ZOOMED_SEARCH(c.snr, spindown)

RECURSIVE_SEARCH(c.snr, psi)
RECURSIVE_SEARCH(c.snr, iota)
RECURSIVE_SEARCH(c.snr, ra)
RECURSIVE_SEARCH(c.snr, dec)

// CHASE(PLAIN_SNR, psi, M_PI/16.0)
// CHASE(PLAIN_SNR, iota, M_PI/16.0)
// CHASE(PLAIN_SNR, ra, 2.0*resolution)
// CHASE(PLAIN_SNR, dec, 2.0*resolution)
// CHASE(PLAIN_SNR, spindown, 1.0/(1800.0*(max_gps()-min_gps())))
// CHASE(PLAIN_SNR, frequency, 1/3600.0)

//  i           h0        psi       iota       ra         dec     spindown f0
//  5 6.606819e-24 2.86832137 0.35185479 3.830254 -0.18512499 3.928912e-10 149.1051

// ORACLE(snr, psi, 2.8683, 1)
// ORACLE(snr, iota, 0.3518, 1)
// ORACLE(snr, ra, 3.830254, 1)
// ORACLE(snr, dec, -0.18512, 1)
//ORACLE(snr, spindown, 3.441452e-11, 1)
// ORACLE(snr, frequency, 149.10508763, 1)



int search_vec(CANDIDATE *cand, float freq_step, float psi_step, float iota_step, float ra_step, float dec_step, float spindown_step, int N)
{
CANDIDATE c, best_c;
int i, max_i;
float f, max;

/* skip 0 vectors */
if(fabs(psi_step)+fabs(iota_step)+fabs(ra_step)+fabs(dec_step)+fabs(spindown_step)==0)return 0;

memcpy(&c, cand, sizeof(*cand)); 
compute_scores(&c, 0);
c.frequency=c.f_max;
compute_scores(&c, 0);

max_i=0;
max=cand->snr;
if(c.snr>max)max=c.snr;

/*fprintf(stderr, " vec(freq=%g, psi=%f, iota=%f, ra=%f, dec=%f, spindown=%f) ", freq_step, psi_step, iota_step, ra_step, dec_step, spindown_step);*/
for(i=-N;i<=N;i++) {
	c.frequency=cand->frequency+i*freq_step;
	c.psi=cand->psi+i*psi_step;
	c.iota=cand->iota+i*iota_step;
	c.ra=cand->ra+i*ra_step;
	c.dec=cand->dec+i*dec_step;
	c.spindown=cand->spindown+i*spindown_step;

	compute_scores(&c, 0);
	f=c.snr;
/*	fprintf(stderr, ",%f", f); */
	if(f>max) {
		memcpy(&best_c, &c, sizeof(c));
		max_i=i;
		max=f;
		}
	}
/*fprintf(stderr, "\n");*/
if(max<cand->snr+0.001)return 0;

if(max_i)memcpy(cand, &best_c, sizeof(best_c));

return max_i!=0;
}

int search_monte_carlo_vec(CANDIDATE *cand)
{
int niter, ngap;
float freq_dir=0, psi_dir=0, iota_dir=0, ra_dir=0, dec_dir=0, spindown_dir=0;

niter=0;
ngap=0;
while(1) {

	freq_dir=2.0*rand()/RAND_MAX-1.0;
/*	psi_dir=2.0*rand()/RAND_MAX-1.0;
	iota_dir=2.0*rand()/RAND_MAX-1.0;*/
	ra_dir=2.0*rand()/RAND_MAX-1.0;
	dec_dir=2.0*rand()/RAND_MAX-1.0;
	spindown_dir=2.0*rand()/RAND_MAX-1.0;	

	if(search_vec(cand, 
			freq_dir*0.05/1800.0,
			psi_dir*M_PI/256.0, 
			iota_dir*M_PI/128, 
			ra_dir*resolution*0.5*0.125/(cos(cand->dec)+0.001), 
			dec_dir*resolution*0.5*0.125, 
			spindown_dir*1e-11, 
			128)) {

			compute_scores(cand, 0);
			DISTANCE_PRINTF
			output_candidate(stderr, "", cand);

			ngap=0;
			}

	ngap++;
	niter++;
	if(ngap> niter*0.5+20) {
		return(ngap<niter);
		}
	}
}

int search_four_vec(CANDIDATE *cand)
{
int ra_dir, dec_dir, spindown_dir, frequency_dir;
int found=0;
CANDIDATE c, best_c;
int i, N=16;
float max;
SCORE_AUX_DATA *ad;


ad=allocate_score_aux_data();


memcpy(&c, cand, sizeof(*cand)); 

fill_e(ad, &c);
fill_response(ad, &c);
fill_doppler(ad);
fill_diff_shift(ad, &c);
print_diff_shift(stderr, ad);

compute_simple_snr(ad, &c, 0);

c.frequency=c.f_max;

fprintf(stderr, "%f\n", c.frequency);

max=cand->snr;
if(c.snr>max)max=c.snr;

compute_simple_snr(ad, &c, 0);
if(c.snr>max)max=c.snr;


for(i=1;i<=N;i++) {

for(ra_dir=-1; ra_dir<=1;ra_dir++) {
for(dec_dir=-1; dec_dir<=1;dec_dir++) {

// 	c.ra=cand->ra+i*ra_dir*resolution*0.5*0.125/(cos(cand->dec)+0.001);
// 	c.dec=cand->dec+i*dec_dir*resolution*0.5*0.125;

	c.ra=cand->ra+i*ra_dir*(ad->d_ra[1]-ad->d_ra[0])*0.1;
	c.dec=cand->dec+i*dec_dir*(ad->d_ra[1]-ad->d_ra[0])*0.1;

	fill_e(ad, &c);
	fill_response(ad, &c);
	fill_doppler(ad);

	for(spindown_dir=-1; spindown_dir<=1;spindown_dir++) {
		c.spindown=cand->spindown+i*spindown_dir*(ad->d_spindown[1]-ad->d_spindown[0])*0.1;
		c.frequency=cand->f_max;

	for(frequency_dir=-5; frequency_dir<5; frequency_dir++) {

		c.frequency=cand->frequency+(fabs(spindown_dir)*(i-1)+1.0)*frequency_dir*(ad->d_freq[1]-ad->d_freq[0])*0.1;
	
		compute_simple_snr(ad, &c, 0);

		//compute_scores(&c, 0);
	
	/*	fprintf(stderr, ",%f", f); */
		if(c.snr>max) {
			fprintf(stderr, "%d %d %d %d\n", ra_dir, dec_dir, spindown_dir, frequency_dir);
			memcpy(&best_c, &c, sizeof(c));
			max=c.snr;
			found=1;
			}

		if(fabs(c.f_max-c.frequency)>(0.5/1800.0)) {
			c.frequency=c.f_max;
			compute_simple_snr(ad, &c, 0);
			if(c.snr>max) {
				fprintf(stderr, "%d %d %d %d\n", ra_dir, dec_dir, spindown_dir, frequency_dir);
				memcpy(&best_c, &c, sizeof(c));
				max=c.snr;
				found=1;
				}

			}

		}}
	}}

	if(found) {
		memcpy(cand, &best_c, sizeof(best_c));
		free_score_aux_data(ad);
		return 1;
		}

	}

free_score_aux_data(ad);

return 0;
/*fprintf(stderr, "\n");*/
if(max<cand->snr+0.001)return 0;

memcpy(cand, &best_c, sizeof(best_c));

return 1;
}

int search_gradient_vec(CANDIDATE *cand, int N)
{
CANDIDATE c, best_c;
int i, max_i, dir_ra, dir_dec, dir_f, dir_sp;
float f, max;
float a, step=1.0, norm;
int original_index=-1;

SCORE_AUX_DATA *ad;

ad=allocate_score_aux_data();

step=0.1;

init:

cand->better_candidate=find_better_candidate(find_next_candidate(cand));
original_index=add_candidate(cand, original_index);
if(cand->better_candidate>=0) {
	original_index=cand->better_candidate;
	memcpy(cand, get_candidate(original_index), sizeof(CANDIDATE));
	free_score_aux_data(ad);
	return 1;
	}

search_alignment(cand);

memcpy(&c, cand, sizeof(*cand)); 
memcpy(&best_c, cand, sizeof(c));

fill_e(ad, &c);
fill_response(ad, &c);
fill_doppler(ad);
fill_matched_diff_shift(ad, &c);
print_diff_shift(stderr, ad);
compute_matched_snr(ad, &c, 0);


max_i=0;
max=cand->snr;
if(c.snr>max)max=c.snr;
max*=1.000001; /* this makes sure we don't have infinite loops due to precision issues */

// if(fabs(step*ad->d_freq[1]) > 0.1/1800.0)step=0.1/(1800.0*ad->d_freq[1]);
// if(fabs(step*ad->d_spindown[1]) > 1e-12)step= 1e-14/(ad->d_spindown[1]);

norm=(ad->d_freq[1]*ad->d_freq[1]+ad->d_ra[1]*ad->d_ra[1]+ad->d_dec[1]*ad->d_dec[1]+ad->d_spindown[1]*ad->d_spindown[1]);

max_i=0;
//fprintf(stderr, "rank=%d gradient search init snr=%f step=%g norm=%g\n", c.rank, max, step, norm);

start:

/* R code:

x<-matrix(0, 81, 4)
m<-1
for(i in -1:1)
for(j in -1:1)
for(k in -1:1)
for(l in -1:1)
{
x[m, 1]<-i
x[m, 2]<-j
x[m, 3]<-k
x[m, 4]<-l
m<-m+1
}
y<-x[,1]+x[,2]+x[,3]+x[,4]

Total: 81
  y>=0 50
  y>=0 & y<=1 35

*/

//for(i=-N;i<=N;i++)
for(dir_f=-1;dir_f<=1;dir_f++)
for(dir_sp=-1;dir_sp<=1;dir_sp++)
for(dir_ra=-1;dir_ra<=1;dir_ra++)
for(dir_dec=-1;dir_dec<=1;dir_dec++) {

	/* do not go against the gradient */
	if(dir_f+dir_sp+dir_ra+dir_dec< -1)continue;
	/* I have never seen the search to go along the gradient "too much": */
	if(dir_f+dir_sp+dir_ra+dir_dec>1)continue;

	for(i=1;;i++) {
	//if(!i)continue;

	c.frequency=cand->frequency+dir_f*i*step/ad->d_freq[1];
	c.psi=cand->psi-i*step*0;
	c.iota=cand->iota-i*step*0;
	c.ra=cand->ra+dir_ra*i*step/ad->d_ra[1];
	c.dec=cand->dec+dir_dec*i*step/ad->d_dec[1];
	c.spindown=cand->spindown+dir_sp*i*step/ad->d_spindown[1];

// 	norm=0;
// 	if(dir_f)norm+=ad->d_freq[1]*ad->d_freq[1];
// 	if(dir_sp)norm+=ad->d_spindown[1]*ad->d_spindown[1];
// 	if(dir_ra)norm+=ad->d_ra[1]*ad->d_ra[1];
// 	if(dir_dec)norm+=ad->d_dec[1]*ad->d_dec[1];
// 	if(norm==0)continue;
// 	norm=1.0/norm;

// 	c.frequency=cand->frequency-dir_f*i*step*norm*ad->d_freq[1];
// 	c.psi=cand->psi+i*step*0;
// 	c.iota=cand->iota+i*step*0;
// 	c.ra=cand->ra-dir_ra*i*step*norm*ad->d_ra[1];
// 	c.dec=cand->dec-dir_dec*i*step*norm*ad->d_dec[1];
// 	c.spindown=cand->spindown-dir_sp*i*step*norm*ad->d_spindown[1];

	fill_e(ad, &c);
	fill_response(ad, &c);
	fill_doppler(ad);
	compute_matched_snr(ad, &c, 0);

	f=c.snr;
	if(f>max) {
		//fprintf(stderr, "%f * [%d,%d,%d,%d]\n", f, dir_f, dir_sp, dir_ra, dir_dec);
		memcpy(&best_c, &c, sizeof(c));
		//output_candidate(LOG, "_intermediate", &(candidate[i]));
		max_i=i;
		max=f;
		continue;
		}

	if(best_c.snr>cand->snr && best_c.snr>=max) {
		memcpy(cand, &best_c, sizeof(best_c));
		goto init;
		}
	/* failed at this i - try different combination */
	break;
	}
}

if((max-cand->snr<=0)) {
	step*=0.1;
	if(step>=1e-7)goto start;
	}

free_score_aux_data(ad);

if(!max_i || (max-cand->snr)<=0)return 0;

fprintf(stderr, " max_i=%d\n", max_i);

/*fprintf(stderr, "\n");*/

memcpy(cand, &best_c, sizeof(best_c));

return 1;
}


int search_all(CANDIDATE *cand)
{
int ra_dir, dec_dir, spindown_dir, frequency_dir;
int found=0, lfound;
double frequency;
CANDIDATE c, best_c;
int i, N=5, N_sky=10, sign;
float max, lmax;
SCORE_AUX_DATA *ad;
SNR_DATA sd0, sd;
float timebase=max_gps()-min_gps();


ad=allocate_score_aux_data();


memcpy(&c, cand, sizeof(*cand)); 

fill_e(ad, &c);
fill_response(ad, &c);
fill_doppler(ad);
fill_diff_shift(ad, &c);
print_diff_shift(stderr, ad);

compute_matched_snr(ad, &c, 0);

if(cand->snr<0)cand->snr=c.snr;

c.frequency=c.f_max;

fprintf(stderr, "%f\n", c.frequency);

max=cand->snr;
if(c.snr>max)max=c.snr;

compute_simple_snr(ad, &c, 0);
if(c.snr>max)max=c.snr;

fprintf(stderr, "c.snr=%f\n", c.snr);


for(dec_dir=-N_sky; dec_dir<=N_sky;dec_dir++) {
for(ra_dir=-N_sky; ra_dir<=N_sky;ra_dir++) {
// 	c.ra=cand->ra+i*ra_dir*resolution*0.5*0.125/(cos(cand->dec)+0.001);
// 	c.dec=cand->dec+i*dec_dir*resolution*0.5*0.125;

	c.ra=cand->ra+ra_dir*(ad->d_ra[1]-ad->d_ra[0])*0.025;
	c.dec=cand->dec+dec_dir*(ad->d_ra[1]-ad->d_ra[0])*0.025;
	c.spindown=cand->spindown;
	c.frequency=cand->frequency;

	fill_e(ad, &c);
	fill_response(ad, &c);
	fill_doppler(ad);

	compute_matched_snr(ad, &c, 0);
	frequency=c.f_max;
	//fprintf(stderr,"[%2d]", (int)round(10*c.snr/cand->snr));

	lmax=0;
	lfound=0;

	for(spindown_dir=-N; spindown_dir<=N;spindown_dir++) {
		c.spindown=cand->spindown+spindown_dir*(ad->d_spindown[1]-ad->d_spindown[0])*0.025;
		//c.frequency=cand->f_max;

	//memcpy(&sd, &sd0, sizeof(sd0));
	sign=(((spindown_dir+N) & 1)<<1)-1;
	if(sign!=-1 && sign!=1)fprintf(stderr, "sign=%d\n", sign);
	for(frequency_dir=-5; frequency_dir<=5; frequency_dir++) {

		c.frequency=cand->frequency-(c.spindown-cand->spindown)*timebase*0.5+sign*frequency_dir*(ad->d_freq[1]-ad->d_freq[0])*0.025;

		if(0) {	
		if(spindown_dir==-N && frequency_dir==-5)
			matched_precompute_snr_data(&sd, ad, &c, 0);
			else {
			matched_precompute_snr_data(&sd, ad, &c, 0);
			//update_snr_data(&sd, ad, &c, 0);
			if(0 && sd.iter_count*10>4001) {
				fprintf(stderr, "*%d %d %d|", sign, spindown_dir, frequency_dir);
				matched_precompute_snr_data(&sd, ad, &c, 0);
				//memcpy(&sd, &sd0, sizeof(sd0));
				}
			}
//		fprintf(stderr, "[%d %d] ", sd.iter_count, sd0.iter_count);
		compute_alignment_snr(&sd, &c, 0);
		} else {
		compute_matched_snr(ad, &c, 0);
		}

		//fprintf(stderr, " %f %f %f\n", c.snr, lmax, cand->snr);
		//compute_scores(&c, 0);
	
	/*	fprintf(stderr, ",%f", f); */
		if(c.snr>lmax)lmax=c.snr;
		if(c.snr>max) {
			//fprintf(stderr, "%d %d %d %d %f\n", ra_dir, dec_dir, spindown_dir, frequency_dir, c.snr);
			memcpy(&best_c, &c, sizeof(c));
			max=c.snr;
			found=1;
			lfound=1;
			}

		if(0 && fabs(c.f_max-c.frequency)>(0.5/1800.0)) {
			c.frequency=c.f_max;

			memcpy(&sd, &sd0, sizeof(sd0));
			update_snr_data(&sd, ad, &c, 0);
			if(sd.iter_count*10>sd0.iter_count) {
				fprintf(stderr, "#");
				precompute_snr_data(&sd0, ad, &c, 0);
				memcpy(&sd, &sd0, sizeof(sd0));
				}
			compute_alignment_snr(&sd, &c, 0);

			if(c.snr>lmax)lmax=c.snr;
			if(c.snr>max) {
				fprintf(stderr, "%d %d %d %d %f\n", ra_dir, dec_dir, spindown_dir, frequency_dir, c.snr);
				memcpy(&best_c, &c, sizeof(c));
				max=c.snr;
				found=1;
				}

			}

		}}
	if(lfound)fprintf(stderr, " *");
		else fprintf(stderr,"%2d", (int)round(10*lmax/cand->snr));
	if(lmax<0.8*cand->snr) {
		ra_dir++;
		fprintf(stderr, " .");
		}
	if(lmax<0.7*cand->snr) {
		ra_dir++;
		fprintf(stderr, " .");
		}
	if(lmax<0.6*cand->snr) {
		ra_dir++;
		fprintf(stderr, " .");
		}
	}
	fprintf(stderr, "\n");
	}

if(found) {
	memcpy(cand, &best_c, sizeof(best_c));
	free_score_aux_data(ad);
	return 1;
	}


free_score_aux_data(ad);

return 0;
/*fprintf(stderr, "\n");*/
if(max<cand->snr+0.001)return 0;

memcpy(cand, &best_c, sizeof(best_c));

return 1;
}

int chase_frequency1(CANDIDATE *cand)
{
CANDIDATE c;
float step=1/3600.0;
int status=0;
float var_start;
SCORE_AUX_DATA *ad;

//fprintf(stderr, "%s %d\n", __FUNCTION__, __LINE__);

memcpy(&c, cand, sizeof(*cand));
ad=allocate_score_aux_data();

fill_e(ad, cand);
fill_response(ad, cand);
fill_doppler(ad);

compute_matched_snr(ad, cand, 0);

var_start=cand->snr;

while(fabs(step)>1e-6) {
	//c.frequency=cand->f_max+step;
	c.frequency=cand->frequency+step;

	compute_matched_snr(ad, &c, 0);

	//fprintf(stderr, "%s step=%g  %g\n", __FUNCTION__, step, c.snr);

	if(c.snr>cand->snr) {
		//fprintf(stderr, "better frequency found\n");
		memcpy(cand, &c, sizeof(CANDIDATE));
		status=1;
		step=step/2.0;
		continue;
		} else {
		c.frequency-=step;
		if(status==0) {
			step=-step;
			}
		status=0;
		step=step/2.0;
		}
	}
if(var_start< cand->snr) {
//	fprintf(stderr, "New frequency=%f\n", best_c.frequency);
	return 1;
	} else {
	return 0;
	}
}

int search_sky(CANDIDATE *cand)
{
int i,j, k, N=50;
CANDIDATE c, best_c;

memcpy(&c, cand, sizeof(CANDIDATE));
memcpy(&best_c, cand, sizeof(CANDIDATE));

compute_scores(&best_c, 0);
fprintf(stderr, "start snr=%f\n", best_c.snr);

for(i=-N;i<=N;i++) {
	c.dec=cand->dec+i*resolution*0.125;
	for(j=-N;j<=N;j++) {
		c.ra=cand->ra+j*resolution*0.125/(cos(c.dec)+0.001);
		for(k=0;k<=3;k++) {
			c.frequency=cand->f_max+k*0.25/1800.0;
			compute_scores(&c, 0);
			if(c.snr>best_c.snr) {
				memcpy(&best_c, &c, sizeof(CANDIDATE));
				//fprintf(stderr, "found snr=%f\n", c.snr);
				//return 1;
				}
			}
		}
	}
//return 0;
if(fabs(best_c.ra-cand->ra)>0 || fabs(best_c.dec-cand->dec)>0) {
	memcpy(cand, &best_c, sizeof(best_c));
	return 1;
	}
return 0;
}

int search_alignment1(CANDIDATE *cand)
{
int i,j, k, N=8;
CANDIDATE c, best_c;

memcpy(&c, cand, sizeof(CANDIDATE));
memcpy(&best_c, cand, sizeof(CANDIDATE));

compute_scores(&best_c, 0);
//fprintf(stderr, "search_alignment start snr=%f\n", best_c.snr);

for(i=-N;i<N;i++) {
	c.iota=cand->iota+i*M_PI/128.0;
	for(j=-N;j<N;j++) {
		c.psi=cand->psi+j*M_PI/128.0;
			compute_scores(&c, 0);
/*			fprintf(stderr, "% 2d", (int)floor(10.0*c.snr/best_c.snr));*/
			if(c.snr>best_c.snr) {
				memcpy(&best_c, &c, sizeof(CANDIDATE));
				//fprintf(stderr, "found snr=%f\n", c.snr);
				//return 1;
				}
		}
/*	fprintf(stderr, "\n");*/
	}
//return 0;
if(fabs(best_c.iota-cand->iota)>0 || fabs(best_c.psi-cand->psi)>0) {
	memcpy(cand, &best_c, sizeof(best_c));
	return 1;
	}
return 0;
}



int search_alignment(CANDIDATE *cand)
{
int i,j, k, N=8, found;
CANDIDATE c, best_c;
SCORE_AUX_DATA *ad;
SNR_DATA sd;

memcpy(&c, cand, sizeof(CANDIDATE));
memcpy(&best_c, cand, sizeof(CANDIDATE));

ad=allocate_score_aux_data();

fill_e(ad, &c);
fill_response(ad, &c);
fill_doppler(ad);

matched_precompute_snr_data(&sd, ad, &best_c, 0);

compute_alignment_snr(&sd, &best_c, 0);
//fprintf(stderr, "search_alignment start snr=%f\n", best_c.snr);

while(1) {
	found=0;

	c.psi=best_c.psi;
	for(i=1;i<N;i++) {
		c.iota=best_c.iota+i*M_PI/128.0;

		compute_alignment_snr(&sd, &c, 0);
		if(c.snr>best_c.snr) {
			memcpy(&best_c, &c, sizeof(CANDIDATE));
			//fprintf(stderr, "found snr=%f\n", c.snr);
			found=1;
			break;
			}
		c.iota=best_c.iota-i*M_PI/128.0;
		compute_alignment_snr(&sd, &c, 0);
		if(c.snr>best_c.snr) {
			memcpy(&best_c, &c, sizeof(CANDIDATE));
			//fprintf(stderr, "found snr=%f\n", c.snr);
			found=1;
			break;
			}
		}

	c.iota=best_c.iota;
	for(j=1;j<N;j++) {
		c.psi=best_c.psi+j*M_PI/128.0;
		compute_alignment_snr(&sd, &c, 0);
		if(c.snr>best_c.snr) {
			memcpy(&best_c, &c, sizeof(CANDIDATE));
			//fprintf(stderr, "found snr=%f\n", c.snr);
			found=1;
			break;
			}
		c.psi=best_c.psi-j*M_PI/128.0;
		compute_alignment_snr(&sd, &c, 0);
		if(c.snr>best_c.snr) {
			memcpy(&best_c, &c, sizeof(CANDIDATE));
			//fprintf(stderr, "found snr=%f\n", c.snr);
			found=1;
			break;
			}
		}
	if(!found)break;
	}

free_score_aux_data(ad);

//fprintf(stderr, "search_alignment stop snr=%f\n", best_c.snr);
if(fabs(best_c.iota-cand->iota)>0 || fabs(best_c.psi-cand->psi)>0) {
	memcpy(cand, &best_c, sizeof(best_c));
	return 1;
	}
return 0;
}


int search_spindown1(CANDIDATE *cand, float step, int N)
{
int i,j, k;
CANDIDATE c, best_c;
SCORE_AUX_DATA *ad;
double timebase=max_gps()-min_gps();

memcpy(&c, cand, sizeof(CANDIDATE));
memcpy(&best_c, cand, sizeof(CANDIDATE));

ad=allocate_score_aux_data();

fill_e(ad, &c);
fill_response(ad, &c);
fill_doppler(ad);

compute_matched_snr(ad, &best_c, 0);
fprintf(stderr, "spindown search start snr=%f step=%g initial %g", best_c.snr, step, cand->spindown);

for(i=-N;i<=N;i++) {
	c.spindown=cand->spindown+i*step;
	for(j=-1;j<=1;j++) {
		c.frequency=cand->frequency-i*step*timebase*0.5+j*0.05/1800.0;
		//c.frequency=c.f_max;
		//chase_frequency1(&c);
		compute_matched_snr(ad, &c, 0);
		fprintf(stderr, ",%f", c.snr);
		if(0 && c.snr>best_c.snr-3) {
			for(k=0;k<=7;k++) {
				c.frequency=cand->f_max-i*step*timebase*0.5+(j+k*0.125)*0.25/1800.0;
				compute_matched_snr(ad, &c, 0);
				if(c.snr>best_c.snr) {
					memcpy(&best_c, &c, sizeof(CANDIDATE));
					fprintf(stderr, "found spindown=%g snr=%f\n", c.spindown, c.snr);
					//return 1;
					}
				}
			}
		}
	}

free_score_aux_data(ad);

fprintf(stderr, "\n");
//return 0;
if(fabs(best_c.spindown-cand->spindown)>0) {
	memcpy(cand, &best_c, sizeof(best_c));
	return 1;
	}
return 0;
}

int search_monte_carlo(CANDIDATE *cand, float sky_r, float spindown_r, int limit)
{
CANDIDATE c, best_c;
SCORE_AUX_DATA *ad;
HINTS h;
float timebase=max_gps()-min_gps();
long i, best_i;

i=0;
best_i=-1;

memcpy(&c, cand, sizeof(CANDIDATE));
memcpy(&best_c, cand, sizeof(CANDIDATE));

ad=allocate_score_aux_data();

fill_e(ad, &best_c);
fill_response(ad, &best_c);
fill_doppler(ad);
compute_snr(ad, &best_c, 0);

//fprintf(stderr, "Monte Carlo search start snr=%f\n", best_c.snr);

while(1) {
	i++;

	c.dec=cand->dec+2.0*(sky_r*rand())/RAND_MAX-sky_r;
	c.ra=cand->ra+(2.0*(sky_r*rand())/RAND_MAX-sky_r)/(cos(cand->dec)+0.001);
	c.spindown=cand->spindown+2.0*(spindown_r*rand())/RAND_MAX-spindown_r;
	c.frequency=cand->frequency-0.5*(c.spindown-cand->spindown)*timebase;

	fill_e(ad, &c);
	fill_response(ad, &c);
	fill_doppler(ad);
	compute_snr(ad, &c, 0);

	if(c.snr>best_c.snr) {
		memcpy(&best_c, &c, sizeof(CANDIDATE));
		best_i=i;
		//fprintf(stderr, "found1 (%ld) snr=%f\n", i, c.snr);
		}

	c.frequency=c.f_max;

	compute_hints(ad, &c, &h);

	if(h.frequency_shift==0.0)continue;

	c.frequency=cand->frequency-0.5*(c.spindown-cand->spindown)*timebase+h.frequency_shift;
	compute_snr(ad, &c, 0);

	if(c.snr>best_c.snr) {
		memcpy(&best_c, &c, sizeof(CANDIDATE));
		best_i=i;
		//fprintf(stderr, "found2 (%ld) snr=%f\n", i, c.snr);
		}

	if((best_i<0) && (i> limit)) {
		free_score_aux_data(ad);		
		return 0;
		}

	if((best_i>0) && (i>2*best_i)) {
		memcpy(cand, &best_c, sizeof(CANDIDATE));
		free_score_aux_data(ad);		
		return 1;
		}
	}
}

int search_four1(CANDIDATE *cand)
{
int i,j, k, m, n,  N=16;
CANDIDATE c, best_c;
SCORE_AUX_DATA *ad;
long start;

start=time(NULL);
memcpy(&c, cand, sizeof(CANDIDATE));
memcpy(&best_c, cand, sizeof(CANDIDATE));

compute_scores(&best_c, 0);
fprintf(stderr, "four param search start snr=%f\n", best_c.snr);

ad=allocate_score_aux_data();

for(i=-N;i<=N;i++) {
	c.dec=cand->dec+i*resolution*0.125;
	for(j=-N;j<=N;j++) {
		c.ra=cand->ra+j*resolution*0.125/(cos(c.dec)+0.001);
		fill_e(ad, &c);
		fill_response(ad, &c);
		fill_doppler(ad);

		for(k=0;k<=3;k++) {
			c.frequency=cand->f_max+k*0.25/1800.0;
			c.spindown=cand->spindown;
			compute_snr(ad, &c, 0);
			if(c.snr>best_c.snr*0.70) {
				for(m=-5;m<=5;m++) {
					c.spindown=cand->spindown+m*1e-11;
					for(n=0;n<=7;n++) {
						c.frequency=cand->f_max+(k+n*0.125)*0.25/1800.0;
						compute_snr(ad, &c, 0);

						if(c.snr>best_c.snr) {
							memcpy(&best_c, &c, sizeof(CANDIDATE));
							//fprintf(stderr, "found snr=%f\n", c.snr);
							//return 1;
							}
						}
					}
				}
			}
		}
	}

free_score_aux_data(ad);
fprintf(stderr, "Time elapsed: %ld\n", time(NULL)-start);
//return 0;
if(fabs(best_c.ra-cand->ra)>0 || fabs(best_c.dec-cand->dec)>0 || fabs(best_c.spindown-cand->spindown)>0) {
	memcpy(cand, &best_c, sizeof(best_c));
	return 1;
	}
return 0;
}

int search_four(CANDIDATE *cand)
{
int i,j, k, m, N=16, N_small;
CANDIDATE c, best_c;
SCORE_AUX_DATA *ad;
HINTS h;
float a;
float sky_step=resolution*0.25;
float timebase=max_gps()-min_gps();
long start;
float *map;
int improved=0;
float spindown_step;

spindown_step=0.87/(1800.0*timebase);

fprintf(stderr, "spindown_step=%g\n", spindown_step);

start=time(NULL);
memcpy(&c, cand, sizeof(CANDIDATE));
memcpy(&best_c, cand, sizeof(CANDIDATE));

ad=allocate_score_aux_data();

map=do_alloc(4*N*N+4*N+1, sizeof(*map));

#define MAP(i,j)   map[(i+N)+(j+N)*(2*N+1)]


fill_e(ad, &best_c);
fill_response(ad, &best_c);
fill_doppler(ad);
compute_snr(ad, &best_c, 0);

fprintf(stderr, "four param search start snr=%f\n", best_c.snr);

for(i=-N;i<=N;i++) {
	c.dec=cand->dec+i*sky_step;
	for(j=-N;j<=N;j++) {
		c.ra=cand->ra+j*sky_step/(cos(c.dec)+0.001);

		fill_e(ad, &c);
		fill_response(ad, &c);

		fill_doppler(ad);

		c.frequency=cand->frequency;

		compute_snr(ad, &c, 0);

		MAP(i,j)=c.snr;

		if(c.snr>best_c.snr) {
			memcpy(&best_c, &c, sizeof(CANDIDATE));
			improved= i || j;
			fprintf(stderr, "found0 (%d,%d) snr=%f\n", i, j, c.snr);
			//return 1;
			}

		if(c.snr<0.3*best_c.snr) {
			j++;
			}

		if(c.snr<0.5*best_c.snr) {
			j++;
			}

		#if 0
		c.frequency=c.f_max;

		compute_hints(ad, &c, &h);

		if(h.frequency_shift==0.0)continue;

		c.frequency=cand->frequency+h.frequency_shift;
		compute_snr(ad, &c, 0);

		if(c.snr > MAP(i,j)) MAP(i,j)=c.snr;
		//if(a>c.snr)fprintf(stderr, "snr=%f new snr=%f shift=%f\n", a, c.snr, h.frequency_shift*1800.0);

		if(c.snr>best_c.snr) {
			memcpy(&best_c, &c, sizeof(CANDIDATE));
			improved=1;
			fprintf(stderr, "found0b (%d,%d) snr=%f\n", i, j, c.snr);
			//return 1;
			}
		#endif
		}
	}

if(candidate_distance(cand, best_c.ra, best_c.dec)>2.0*resolution) {
	free_score_aux_data(ad);
	free(map);
	fprintf(stderr, "Time elapsed: %ld\n", time(NULL)-start);
	memcpy(cand, &best_c, sizeof(best_c));
	return 1;
	}

if(improved) {
	N_small=ceil(2.0*resolution/sky_step);
	} else {
	N_small=N;
	}

fprintf(stderr, "N_small=%d\n", N_small);

for(i=-N_small;i<=N_small;i++) {
	c.dec=cand->dec+i*sky_step;
	for(j=-N_small;j<=N_small;j++) {

		#define L(a, b)  (MAP(i, j) < MAP(i+a, j+b))

		#if 0
		if( (i!=N_small) && (i != -N_small) &&
		    (j!=N_small) && (j != -N_small)) {

			if( (L(0,1) && L(-1, 1) && L(1, 1)) || 
			(L(0, -1) && L(-1, -1) && L(1, -1)) ||
			(L(1, 0) && L(1, -1) && L(1, 1)) || 
			(L(-1, 0) && L(-1, -1) && L(-1, 1)) ) {
				fprintf(stderr, " .");
				continue;
				}
			}
		#endif
			
		if(MAP(i,j)<0.5*best_c.snr) {
			fprintf(stderr, " .");
			continue;
			}

	//	fprintf(stderr, "%2d", (int)trunc(10*MAP(i,j)/best_c.snr));

		a=MAP(i,j);

		c.ra=cand->ra+j*sky_step/(cos(c.dec)+0.001);
		fill_e(ad, &c);
		fill_response(ad, &c);
		fill_doppler(ad);

		for(m=-5;m<=5;m++) {
			c.spindown=cand->spindown+m*spindown_step*0.2;

			c.frequency=cand->frequency-m*spindown_step*0.2*timebase*0.5;
			compute_snr(ad, &c, 0);

			if(c.snr>best_c.snr) {
				memcpy(&best_c, &c, sizeof(CANDIDATE));
				improved=1;
				fprintf(stderr, "found1 (%d,%d;%d) snr=%f\n", i, j, m, c.snr);
				//return 1;
				}


			if(c.snr>a)a=c.snr;

			c.frequency=c.f_max;

			compute_hints(ad, &c, &h);

			if(h.frequency_shift==0.0)continue;

			c.frequency=cand->frequency-m*(0.25*0.1/1800.0)+h.frequency_shift;
			compute_snr(ad, &c, 0);

			if(c.snr>best_c.snr) {
				memcpy(&best_c, &c, sizeof(CANDIDATE));
				improved=1;
				fprintf(stderr, "found2 (%d,%d;%d) snr=%f\n", i, j, m, c.snr);
				//return 1;
				}

			if(c.snr>a)a=c.snr;


/*			c.frequency=cand->frequency+k/3600.0+h.frequency_shift+0.25/(1800*FREQ_STEPS);
			compute_snr(ad, &c, 0);

 			//if(a>c.snr)fprintf(stderr, "snr=%f new snr=%f shift=%f\n", a, c.snr, h.frequency_shift*1800.0);

			if(c.snr>best_c.snr) {
				memcpy(&best_c, &c, sizeof(CANDIDATE));
				fprintf(stderr, "found3 snr=%f k=%d\n", c.snr, k);
				//return 1;
				}*/
			}

		if(a>0.80*best_c.snr) {
			c.spindown=cand->spindown;
			c.frequency=cand->frequency;
			search_monte_carlo(&c, sky_step, spindown_step, 2.0/(1.001-a/best_c.snr));

			if(c.snr>best_c.snr) {
				memcpy(&best_c, &c, sizeof(CANDIDATE));
				improved=1;
				fprintf(stderr, "found3 (%d,%d) snr=%f\n", i, j, c.snr);
				//return 1;
				}
			}

		fprintf(stderr, "%2d", (int)trunc(10*a/best_c.snr));
		}
	fprintf(stderr, "\n");
	}

free(map);
free_score_aux_data(ad);
fprintf(stderr, "Time elapsed: %ld\n", time(NULL)-start);

if(improved) {
	memcpy(cand, &best_c, sizeof(best_c));
	return 1;
	}
return 0;
}

int search_four3(CANDIDATE *cand)
{
int i,j, k, m, N=32;
CANDIDATE c, best_c;
SCORE_AUX_DATA *ad;
HINTS h;
float a;
float sky_step=resolution*0.125;
float timebase=max_gps()-min_gps();
long start;
int improved=0;
float spindown_step;

spindown_step=1.0/(1800.0*timebase);

fprintf(stderr, "spindown_step=%g\n", spindown_step);

start=time(NULL);
memcpy(&c, cand, sizeof(CANDIDATE));
memcpy(&best_c, cand, sizeof(CANDIDATE));

ad=allocate_score_aux_data();

fill_e(ad, &best_c);
fill_response(ad, &best_c);
fill_doppler(ad);
compute_snr(ad, &best_c, 0);

fprintf(stderr, "four param search start snr=%f\n", best_c.snr);


for(i=-N;i<=N;i++) {
	c.dec=cand->dec+i*sky_step;
	for(j=-N;j<=N;j++) {

		c.ra=cand->ra+j*sky_step/(cos(c.dec)+0.001);

		fill_e(ad, &c);
		fill_response(ad, &c);
		fill_doppler(ad);

		c.spindown=cand->spindown;
		c.frequency=cand->frequency;
		compute_snr(ad, &c, 0);
		a=c.snr;
		if(a<0.7*best_c.snr) {
			fprintf(stderr, " .");
			continue;
			}
		fprintf(stderr, "%2d", (int)trunc(10*a/best_c.snr));

		for(m=-10;m<=10;m++) {
			c.spindown=cand->spindown+m*spindown_step*0.1;

			for(k=-4;k<=4;k++) {
				c.frequency=cand->frequency-m*spindown_step*0.1*timebase*0.5+k*0.015/1800.0;
				compute_snr(ad, &c, 0);
	
				if(c.snr>a)a=c.snr;
	
				if(c.snr>best_c.snr) {
					memcpy(&best_c, &c, sizeof(CANDIDATE));
					improved=1;
					fprintf(stderr, "found (%d,%d;%d;%d) snr=%f\n", i, j, m, k, c.snr);
					//return 1;
					}
	
			}

			if(a<best_c.snr*0.5)m++;
			}

		}
	fprintf(stderr, "\n");
	}

free_score_aux_data(ad);
fprintf(stderr, "Time elapsed: %ld\n", time(NULL)-start);

if(improved) {
	memcpy(cand, &best_c, sizeof(best_c));
	return 1;
	}
return 0;
}


void tabulate_neighbourhood(char *filename, CANDIDATE *cand)
{
FILE *fout;
int i, j, k, m, N=16;
CANDIDATE c;
fout=fopen(filename, "w");
if(fout==NULL){
	perror("");
	return;
	}
fprintf(fout, "ra_ind\tdec_ind\tfreq_ind\tspindown_ind\tra\tdec\tfrequency\tspindown\tsnr\tsnr_max\tf_max\n");
memcpy(&c, cand, sizeof(c));

for(i=-N;i<=N;i++) {
for(j=-N;j<=N;j++) {
for(k=-N;k<=N;k++) {
for(m=-N;m<=N;m++) {
		c.dec=cand->dec+i*resolution*0.5*0.125;
		c.ra=cand->ra+j*resolution*0.5*0.125/(cos(c.dec)+0.001);
		c.spindown=cand->spindown+k*1e-11;
		c.frequency=cand->frequency+m*0.02/1800.0;
		compute_scores(&c, 0);
		fprintf(fout, "%d\t%d\t%d\t\%d\t%f\t%f\t%f\t%g\t%f", j, i, m, k, c.ra, c.dec, c.frequency, c.spindown, c.snr);
		c.frequency=c.f_max;
		compute_scores(&c, 0);
		fprintf(fout, "\t%f\t%f\n", c.snr, c.f_max);
	}}}}
fclose(fout);
}

void optimize_candidate(CANDIDATE *cand)
{
compute_scores(cand, 0);
output_candidate_header(stderr);
output_candidate(stderr, "", cand);
DISTANCE_PRINTF

if(cand->rank==0)test_alignment_snr(cand);


// tabulate_neighbourhood("cand_0.txt", cand);
// exit(0);

search_gradient_vec(cand, 16);


compute_scores(cand, 0);
//fprintf(LOG, "distance=%f snr=%f\n", candidate_distance(cand, args_info.focus_ra_arg, args_info.focus_dec_arg), cand->snr);
DISTANCE_PRINTF
output_candidate(stderr, "", cand);
}

int compare_opt_ranks(int *a, int *b)
{
if(candidate[*a].opt_rank<0 && candidate[*b].opt_rank<0)return 0;
if(candidate[*a].opt_rank<0)return 1;
if(candidate[*b].opt_rank<0)return -1;
if(candidate[*a].snr<candidate[*b].snr)return 1;
if(candidate[*a].snr>candidate[*b].snr)return -1;
return 0;
}

void assign_opt_ranks(void)
{
int *index;
int i;

index=do_alloc(candidate_free, sizeof(*index));

for(i=0;i<candidate_free;i++)index[i]=i;

qsort(index, candidate_free, sizeof(*index), compare_opt_ranks);

for(i=0;i<candidate_free;i++) {
	if(candidate[index[i]].opt_rank<0)continue;
	candidate[index[i]].opt_rank=i;
	}
free(index);
}

void identify_candidates(void)
{
RGBPic *p;
PLOT *plot;
SUM_TYPE a;
int i,k,m;
float *a_f;
time_t start, now;
int opt_candidates_count;

/* compute power without applying whitening procedure */

fprintf(stderr, "Identifying candidates - pass 1\n");

candidate_free=0;
candidate_size=100;
candidate=do_alloc(candidate_size, sizeof(*candidate));

max_dx_local_map=do_alloc(fine_grid->npoints, sizeof(*max_dx_local_map));
max_dx_order=do_alloc(fine_grid->npoints, sizeof(*max_dx_order));

dec_order=do_alloc(fine_grid->npoints, sizeof(*dec_order));
inverse_dec_order=do_alloc(fine_grid->npoints, sizeof(*inverse_dec_order));

ra_order=do_alloc(fine_grid->npoints, sizeof(*ra_order));
inverse_ra_order=do_alloc(fine_grid->npoints, sizeof(*inverse_ra_order));

a_f=do_alloc(fine_grid->npoints, sizeof(*a_f));

if(fine_grid->max_n_dec<800){
	p=make_RGBPic(fine_grid->max_n_ra*(800/fine_grid->max_n_dec)+140, fine_grid->max_n_dec*(800/fine_grid->max_n_dec));
	} else 
	p=make_RGBPic(fine_grid->max_n_ra+140, fine_grid->max_n_dec);	

plot=make_plot(p->width, p->height);

for(i=0;i<fine_grid->npoints;i++) {
	#if 0
	max_dx[i]=0;
	max_dx_polarization_index[i]=-1;
	#endif
	max_dx_order[i]=i;
	dec_order[i]=i;
	ra_order[i]=i;
	max_dx_local_map[i]=-1;
	}

#if 0
for(i=0;i<fine_grid->npoints;i++) {
	for(k=0;k<ntotal_polarizations;k++){
		a=polarization_results[k].skymap.max_dx[i];
		if(a<0)continue;
		if(a>max_dx[i]) {
			max_dx[i]=a;
			max_dx_polarization_index[i]=k;
			}
		}
	}
#endif

qsort(max_dx_order, fine_grid->npoints, sizeof(*max_dx_order), max_dx_compare);
qsort(dec_order, fine_grid->npoints, sizeof(*dec_order), dec_compare);
qsort(ra_order, fine_grid->npoints, sizeof(*ra_order), ra_compare);

for(i=0;i<fine_grid->npoints;i++) {
	inverse_dec_order[dec_order[i]]=i;
	inverse_ra_order[ra_order[i]]=i;
	}

/* Take noise floor to be bottom 1/3 of max_dx above 2 */
for(i=0;i<fine_grid->npoints;i++)
	if(max_dx[max_dx_order[i]]<2.0)break;
noise_floor=max_dx[max_dx_order[(i*2)/3]];
fprintf(LOG, "noise_floor: %f\n", noise_floor);
fprintf(stderr, "noise_floor: %f\n", noise_floor);
	
/* Find all local maxima, starting with largest.. */
for(i=0;i<fine_grid->npoints;i++) {
	k=max_dx_order[i];
	/* skip masked points.. */
	if(max_dx_polarization_index[k]<0)continue;
	
	/* is the point marked already ? */
	if(max_dx_local_map[k]>=0) {
		a=polarization_results[max_dx_polarization_index[k]].skymap.max_dx[k];
		if(a>candidate[max_dx_local_map[k]].max_dx)candidate[max_dx_local_map[k]].max_dx=a;
		sweep_neighbourhood(k, max_dx_local_map[k]);
		continue;
		}
	
	/* record found maximum.. */
	candidate[candidate_free].point_index=k;
	
	/* mark and sweep nearby points */
	max_dx_local_map[k]=candidate_free;
	sweep_neighbourhood(k, candidate_free);

	candidate_free++;	
	if(candidate_free>=candidate_size)expand_candidates();
	}

for(i=0;i<candidate_free;i++) {
	candidate[i].domain_size=0;
	}
	
for(i=0;i<fine_grid->npoints;i++) {
	k=max_dx_local_map[i];
	if(k>=0)candidate[k].domain_size++;
	}

fprintf(stderr, "Writing skymaps\n");

for(i=0;i<fine_grid->npoints;i++) {
	a_f[max_dx_order[i]]=fine_grid->npoints-i;
	}
snprintf(s, 19999, "%smax_dx_order.png", subinstance_name);
if(clear_name_png(s)){
	plot_grid_f(p, fine_grid, a_f, 1);
	RGBPic_dump_png(s, p);
	}
snprintf(s, 19999, "%smax_dx_order.dat", subinstance_name);
dump_ints(s, max_dx_order, fine_grid->npoints, 1);

for(i=0;i<fine_grid->npoints;i++) {
	a_f[dec_order[i]]=i;
	}
snprintf(s, 19999, "%sdec_order.png", subinstance_name);
if(clear_name_png(s)){
	plot_grid_f(p, fine_grid, a_f, 1);
	RGBPic_dump_png(s, p);
	}

for(i=0;i<fine_grid->npoints;i++) {
	a_f[ra_order[i]]=i;
	}
snprintf(s, 19999, "%sra_order.png", subinstance_name);
if(clear_name_png(s)){
	plot_grid_f(p, fine_grid, a_f, 1);
	RGBPic_dump_png(s, p);
	}

for(i=0;i<fine_grid->npoints;i++) {
	a_f[i]=max_dx_local_map[i]<0 ? -1 : candidate_free-max_dx_local_map[i];
	}
snprintf(s, 19999, "%smax_dx_local_map.png", subinstance_name);
if(clear_name_png(s)){
	plot_grid_f(p, fine_grid, a_f, 1);
	RGBPic_dump_png(s, p);
	}
snprintf(s, 19999, "%smax_dx_local_map.dat", subinstance_name);
dump_ints(s, max_dx_local_map, fine_grid->npoints, 1);

for(i=0;i<fine_grid->npoints;i++) {
	a_f[i]=max_dx_polarization_index[i];
	}
snprintf(s, 19999, "%smax_dx_max_dx_polarization_index.png", subinstance_name);
if(clear_name_png(s)){
	plot_grid_f(p, fine_grid, a_f, 1);
	RGBPic_dump_png(s, p);
	}
snprintf(s, 19999, "%smax_dx_max_dx_polarization_index.dat", subinstance_name);
dump_shorts(s, max_dx_polarization_index, fine_grid->npoints, 1);

/* Now do detailed processing of the candidates */
fprintf(stderr, "Optimizing candidates (pass 2)\n");

output_candidate_header(LOG);

time(&start);
m=0;	
for(i=0;i<candidate_free;i++) {
	k=candidate[i].point_index;
	//fprintf(stderr, "k=%d\n", k);
	candidate[i].max_dx_polarization_index=max_dx_polarization_index[k];
	//fprintf(stderr, "max_dx_polarization_index=%d\n", max_dx_polarization_index[k]);
	candidate[i].rank=i;
	candidate[i].score=max_dx[k];
	if(candidate[i].score>noise_floor)m++;
	candidate[i].ra=fine_grid->longitude[k];
	candidate[i].dec=fine_grid->latitude[k];
	candidate[i].spindown=spindown;
	/* candidate[i].max_dx=polarization_results[max_dx_polarization_index[k]].skymap.max_dx[k]; */
	candidate[i].psi=polarization_results[max_dx_polarization_index[k]].orientation;
	/* currently we have either circular or linear polarizations */
	if(polarization_results[max_dx_polarization_index[k]].cross_factor>0)candidate[i].iota=0.0;
		else	candidate[i].iota=M_PI/2.0;
	//candidate[i].a_plus=polarization_results[max_dx_polarization_index[k]].plus_proj;
	//candidate[i].a_cross=polarization_results[max_dx_polarization_index[k]].cross_proj;
	candidate[i].ul=polarization_results[max_dx_polarization_index[k]].skymap.max_upper_limit[k];
	candidate[i].S_band=args_info.no_secondary_skymaps_arg ? -1 : polarization_results[max_dx_polarization_index[k]].skymap.S_map[k];
	candidate[i].M_band=args_info.no_secondary_skymaps_arg ? -1 : polarization_results[max_dx_polarization_index[k]].skymap.M_map[k];
	candidate[i].frequency=polarization_results[max_dx_polarization_index[k]].skymap.freq_map[k];
	candidate[i].weight_ratio=1.0-polarization_results[max_dx_polarization_index[k]].skymap.max_sub_weight[k]/polarization_results[max_dx_polarization_index[k]].skymap.total_weight[k];
	candidate[i].skyband=fine_grid->band[k];
	
	//fprintf(stderr, "Computing scores\n");
	compute_scores(&(candidate[i]), 0);
	
	if(args_info.output_initial_arg)output_candidate(LOG, "_initial", &(candidate[i]));

	if( ((args_info.max_candidates_arg<0) || (i<args_info.max_candidates_arg)) && (max_dx[k]>args_info.min_candidate_snr_arg)) {
		//time(&now);
		optimize_candidate(&(candidate[i]));
		//fprintf(stderr, "Time per candidate %f\n", ((now-start)*1.0)/(i+1));
		candidate[i].opt_rank=0;
		} else {
		candidate[i].opt_rank=-1;
		}

	if(i<args_info.dump_candidates_arg)dump_candidate(i);
	}

assign_opt_ranks();

opt_candidates_count=0;
for(i=0;i<candidate_free;i++) {
	if(candidate[i].opt_rank<0)continue;
	if(args_info.output_optimized_arg)output_candidate(LOG, "_optimized", &(candidate[i]));
	opt_candidates_count++;
	}

fprintf(LOG, "optimized_candidates_count: %d\n", opt_candidates_count);
fprintf(stderr, "optimized_candidates_count: %d\n", opt_candidates_count);

fprintf(LOG, "high_candidates_count: %d\n", m);
fprintf(stderr, "high_candidates_count: %d\n", m);

fprintf(LOG, "candidates_count: %d\n", candidate_free);
fprintf(stderr, "candidates_count: %d\n", candidate_free);

time(&now);
fprintf(stderr, "Second pass processing time: %d\n", now-start);
fprintf(LOG, "Second pass processing time: %d\n", now-start);


free(a_f);
a_f=NULL;
free_plot(plot);
free_RGBPic(p);

free(max_dx_local_map);
free(max_dx_order);
max_dx_local_map=NULL;
max_dx_order=NULL;

free(dec_order);
free(inverse_dec_order);
dec_order=NULL;
inverse_dec_order=NULL;

free(ra_order);
free(inverse_ra_order);
ra_order=NULL;
inverse_ra_order=NULL;

free(candidate);
candidate=NULL;
}

void output_candidates(FILE *fout)
{
int i;
int non_zero_count=0;
int leaf_count=0;
float sum=0.0;
int max=0;
for(i=0;i<nbins;i++) {
	if(bin_index[i]->free>0) {
		non_zero_count++;
		sum+=bin_index[i]->free;
		if(bin_index[i]->free>max)max=bin_index[i]->free;
		}
	}
output_candidate_header(fout);
for(i=0;i<candidate_cache->free;i++) {
	if(VELT(candidate_cache, CANDIDATE, i).better_candidate>=0)continue;
	if(find_next_candidate(&VELT(candidate_cache, CANDIDATE, i))>=0)continue;
	output_candidate(fout, "_final", &VELT(candidate_cache, CANDIDATE, i));
	leaf_count++;
	}
fprintf(fout, "candidates cache length: %d\n", candidate_cache->free);
fprintf(fout, "candidates cache leaves: %d\n", leaf_count);
fprintf(fout, "candidates index non zero count: %d\n", non_zero_count);
fprintf(fout, "candidates index max length: %d\n", max);
fprintf(fout, "candidates index average length: %f\n", non_zero_count > 0 ? sum/non_zero_count : 0.0 );
fprintf(fout, "index queries total: %lf\n", index_queries_total);
fprintf(fout, "index hits: %lf\n", index_hits);
fprintf(fout, "index hit ratio: %lf\n", index_hits/index_queries_total);
fprintf(fout, "index average snr change: %lf\n", index_snr_total/index_queries_total);
fprintf(fout, "improvement queries total: %lf\n", improvement_queries_total);
fprintf(fout, "improvements: %lf\n", improvements);
fprintf(fout, "improvement ratio: %lf\n", improvements/improvement_queries_total);
fprintf(fout, "improvement average snr change: %lf\n", improvement_snr_total/improvement_queries_total);
}
