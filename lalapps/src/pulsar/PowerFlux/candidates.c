#include <stdio.h>
#include <stdlib.h>
/* We need this define to get NAN values */
#define __USE_ISOC99
#include <math.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
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

INT64 spindown_start;

CANDIDATE *candidate=NULL;
int candidate_free=0;
int candidate_size=-1;

extern struct gengetopt_args_info args_info;

char s[20000];

SUM_TYPE *max_dx=NULL;
short *polarization_index=NULL;
int *max_dx_order=NULL;
int *dec_order=NULL;
int *inverse_dec_order=NULL;
int *ra_order=NULL;
int *inverse_ra_order=NULL;
int *max_dx_local_map=NULL;

float noise_floor;

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

for(j=inverse_dec_order[index]+1;j<fine_grid->npoints;j++) {
	m=dec_order[j];

	if(fabs(fine_grid->latitude[m]-fine_grid->latitude[index])>fudge)break;

	if(max_dx_local_map[m]>=0)continue;
	
	if(max_dx[m]>max_dx[index])continue;
	
	if(polarization_index[m]<0)continue;

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
	
	if(polarization_index[m]<0)continue;

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
	pl=&(d->polarizations[candidate[i].polarization_index]);	
	
	/* process single SFT */
	for(k=0;k<d->free;k++){
		/* Get amplitude response */
		a_plus=F_plus(k, fine_grid, index, pl->AM_coeffs);
		a_cross=F_plus(k, fine_grid, index, pl->conjugate->AM_coeffs);

		doppler=fine_grid->e[0][index]*d->detector_velocity[3*k+0]+
			fine_grid->e[1][index]*d->detector_velocity[3*k+1]+
			fine_grid->e[2][index]*d->detector_velocity[3*k+2];

		bin_shift=-rint((first_bin+nbins*0.5)*doppler+1800.0*spindown*(d->gps[k]-spindown_start));

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
			spindown*(d->gps[k]-spindown_start),
			a_plus,
			a_cross,
			d->expTMedians[k]);
		for(b=b0; b<b1;b++) {
			a=d->bin[k*nbins+b].im;
			if(a>0)
				fprintf(fout, "\t%g+%gi", d->bin[k*nbins+b].re, a);
				else 
				fprintf(fout, "\t%g%gi", d->bin[k*nbins+b].re, a);
			}

		fprintf(fout, "\n");
		}
	}
fclose(fout);
}

void output_candidate_header(FILE *fout)
{
fprintf(fout, "candidates: label polarization_index rank score point_index domain_size ul S M max_dx frequency psi iota ra dec spindown  weight_ratio skyband coherence_score power_cor snr strain f_max ifo_freq ifo_freq_sd total\n");
}

void output_candidate(FILE *fout, char * suffix, CANDIDATE *cand)
{
fprintf(fout, "candidate%s: \"%s\" %d %d %g %d %d %g %g %g %f %f %f %f %f %f %g %f %d %f %f %f %g %f %f %f %d\n",
	suffix,
	args_info.label_given ? args_info.label_arg : "",
	cand->polarization_index,
	cand->rank,
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
	cand->f_max,
	cand->ifo_freq,
	cand->ifo_freq_sd,
	candidate_free);
}

void compute_scores(CANDIDATE *cand, int debug)
{
float doppler;
DATASET *d;
int index=cand->point_index;
int patch_index=super_grid->reverse_map[index];
float spindown=cand->spindown;
float frequency=cand->frequency;
int b0, b1, j, k, b, b_max;
float a, coherence_score, power_cor, a_plus, a_cross, a_plus_sq, a_cross_sq, f_plus, f_cross, f_plus_sq, f_cross_sq;
POLARIZATION *pl;
double weight, total_weight, *response, *mismatch, f, x, y, *demod_signal_sum, *signal_sum, *signal_sq_sum, response_sum, response_sq_sum, response_weight, demod_weight, total_demod_weight, total_demod_weight2, *cov_sum, power, mean_power, mean_power_sq, power_sd;
int window=20;
int search_window=5;
int *signal_bin;
COMPLEX *p0, *p1, *p2, w, z;
float e[26];
float p_proj, c_proj;
struct {
	double re;
	double im;
	} *phase_sum;

total_weight=0;
response_sq_sum=0.0;
response_sum=0.0;
response_weight=0.0;
total_demod_weight=0.0;
total_demod_weight2=0.0;

phase_sum=do_alloc(2*window+1, sizeof(*phase_sum));
demod_signal_sum=do_alloc(2*window+1, sizeof(*signal_sum));
signal_sum=do_alloc(2*window+1, sizeof(*signal_sum));
signal_sq_sum=do_alloc(2*window+1, sizeof(*signal_sq_sum));
cov_sum=do_alloc(2*window+1, sizeof(*cov_sum));

for(b=0;b<2*window+1;b++) {
	phase_sum[b].re=0.0;
	phase_sum[b].im=0.0;
	signal_sum[b]=0.0;
	signal_sq_sum[b]=0.0;
	cov_sum[b]=0.0;
	demod_signal_sum[b]=0.0;
	}

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

cand->ifo_freq=0.0;
cand->ifo_freq_sd=0.0;

/* loop over datasets */
for(j=0;j<d_free;j++) {
	d=&(datasets[j]);
	//pl=&(d->polarizations[cand->polarization_index]);	
	pl=&(d->polarizations[0]);	
	
	response=do_alloc(d->free, sizeof(*response));
	signal_bin=do_alloc(d->free, sizeof(*signal_bin));
	mismatch=do_alloc(d->free, sizeof(*mismatch));
		
	/* process single SFTs */
	for(k=0;k<d->free;k++) {
		/* Get amplitude response */
		#if 0
		a_plus=F_plus(k, fine_grid, index, pl->AM_coeffs);
		a_cross=F_plus(k, fine_grid, index, pl->conjugate->AM_coeffs);
		
		response[k]=pl->plus_factor*a_plus*a_plus+pl->cross_factor*a_cross*a_cross;
		#else 

		f_plus=F_plus_coeff(k, e, pl->AM_coeffs);
		f_cross=F_plus_coeff(k, e, pl->conjugate->AM_coeffs);

// 		f_plus_sq=f_plus*f_plus;
// 		f_cross_sq=f_cross*f_cross;
// 
// 		a=f_plus*a_plus+f_cross*a_cross;
// 		response[k]=0.25*((f_plus_sq+f_cross_sq)*(a_plus_sq+a_cross_sq)+(f_plus_sq-f_cross_sq)*(a_plus_sq-a_cross_sq)*p_proj+2.0*f_plus*f_cross*(a_plus_sq-a_cross_sq))*c_proj;

		a=f_plus*p_proj+f_cross*c_proj;
		f_cross=f_cross*p_proj-f_plus*c_proj;
		f_plus=a;

		f_plus_sq=f_plus*f_plus;
		f_cross_sq=f_cross*f_cross;

		//response[k]=0.25*((f_plus_sq+f_cross_sq)*(a_plus_sq+a_cross_sq)+(f_plus_sq-f_cross_sq)*(a_plus_sq-a_cross_sq));
		response[k]=f_plus_sq*a_plus_sq+f_cross_sq*a_cross_sq;
		#endif

		doppler=e[0]*d->detector_velocity[3*k+0]+
			e[1]*d->detector_velocity[3*k+1]+
			e[2]*d->detector_velocity[3*k+2];

		f=frequency+frequency*doppler+spindown*(d->gps[k]-spindown_start);

		signal_bin[k]=rint(1800.0*f-first_bin);
		mismatch[k]=1800.0*f-first_bin-signal_bin[k];

		if(!k || (signal_bin[k]<b0))b0=signal_bin[k];
		if(!k || (signal_bin[k]>b1))b1=signal_bin[k];

		}

	if( (b0<0) || (b1>=nbins)) {
		/* we are outside loaded bin range */
		fprintf(stderr, "b0=%d b1=%d\n", b0, b1);
		cand->coherence_score=-1.0;
		cand->power_cor=-1.0;
		cand->snr=-1.0;
		cand->strain=-1.0;
		cand->f_max=cand->frequency;
		cand->ifo_freq=-1;
		cand->ifo_freq_sd=-1;
		return;
		}
	//fprintf(stderr, "b0=%d b1=%d\n", b0, b1);

	for(k=0;k<d->free;k++) {
		/* skip SFTs with low weight */
		if(d->expTMedians[k]<0.05)continue;

		/* power_cor computation */
		weight=d->expTMedians[k]*d->weight;
		response_sum+=response[k]*weight;
		response_sq_sum+=response[k]*response[k]*weight;
		response_weight+=weight;

		demod_weight=weight*response[k];
		total_demod_weight2+=demod_weight;
		total_demod_weight+=demod_weight*response[k];

		p0=&(d->bin[k*nbins+signal_bin[k]-window]);

		f=signal_bin[k]+mismatch[k];
		cand->ifo_freq+=f*demod_weight;		
		cand->ifo_freq_sd+=f*f*demod_weight;		

		for(b=0; b< (2*window+1); b++) {
			x=p0[b].re;
			y=p0[b].im;
			power=x*x+y*y;
			signal_sum[b]+=power*weight;
			signal_sq_sum[b]+=power*power*weight;
			cov_sum[b]+=power*response[k]*weight;

			demod_signal_sum[b]+=power*demod_weight;
			}
		}
		
	/* phase computation */
	for(k=1;k< (d->free-1);k++) {
	
	
		/* skip SFT with asymmetric gaps */
		if(d->gps[k]-d->gps[k-1]!=d->gps[k+1]-d->gps[k])continue;

		/* skip SFTs with low weight */
		if(d->expTMedians[k]<0.01 || d->expTMedians[k-1]<0.01 || d->expTMedians[k+1]<0.01)continue;
	
		weight=((response[k-1]*d->expTMedians[k-1]+2*response[k]*d->expTMedians[k]+response[k+1]*d->expTMedians[k+1])*d->weight);
		total_weight+=weight;
	
		p0=&(d->bin[(k-1)*nbins+signal_bin[k-1]-window]);
		p1=&(d->bin[k*nbins+signal_bin[k]-window]);
		p2=&(d->bin[(k+1)*nbins+signal_bin[k+1]-window]);

		for(b=0; b< (2*window+1); b++) {
			x=sqrt(p0[b].re*p0[b].re+p0[b].im*p0[b].im);
			y=sqrt(p1[b].re*p1[b].re+p1[b].im*p1[b].im);

			w.re=(p0[b].re*p1[b].re+p0[b].im*p1[b].im)/(x*y);
			w.im=(p0[b].re*p1[b].im-p0[b].im*p1[b].re)/(x*y);

			x=sqrt(p2[b].re*p2[b].re+p2[b].im*p2[b].im);
			z.re=(p2[b].re*p1[b].re+p2[b].im*p1[b].im)/(x*y);
			z.im=(p2[b].re*p1[b].im-p2[b].im*p1[b].re)/(x*y);
			
			phase_sum[b].re+=weight*(w.re*z.re-w.im*z.im);
			phase_sum[b].im+=weight*(w.im*z.re+w.re*z.im);
			}
		}
		
	free(response);
	free(signal_bin);
	free(mismatch);
	response=NULL;
	signal_bin=NULL;
	mismatch=NULL;
	}

//fprintf(stderr, "total_weight=%g\n", total_weight);

coherence_score=0;
power_cor=0;
response_sum/=response_weight;
response_sq_sum/=response_weight;

cand->ifo_freq/=total_demod_weight2;
cand->ifo_freq_sd/=total_demod_weight2;

cand->ifo_freq_sd=sqrt(cand->ifo_freq_sd-cand->ifo_freq*cand->ifo_freq);

cand->ifo_freq=(cand->ifo_freq+first_bin)/1800.0;

for(b=0;b<2*window+1;b++) {
	/* power_cor */
	cov_sum[b]/=response_weight;
	signal_sum[b]/=response_weight;
	signal_sq_sum[b]/=response_weight;
	demod_signal_sum[b]/=total_demod_weight;
	}

b_max=window;
for(b=window-search_window;b<window+search_window+1;b++) {
	if(demod_signal_sum[b_max]<demod_signal_sum[b]) b_max=b;
	
	a=(cov_sum[b]-signal_sum[b]*response_sum)/sqrt((signal_sq_sum[b]-signal_sum[b]*signal_sum[b])*(response_sq_sum-response_sum*response_sum));
	if(fabs(a)>fabs(power_cor))power_cor=a;

	/* coherence_score */
	phase_sum[b].re/=total_weight;
	phase_sum[b].im/=total_weight;
	a=sqrt(phase_sum[b].re*phase_sum[b].re+phase_sum[b].im*phase_sum[b].im);
	if(a>coherence_score)coherence_score=a;		
	}

mean_power=0.0;
mean_power_sq=0.0;
for(b=0;b<window-search_window;b++){
	mean_power+=demod_signal_sum[b];
	mean_power_sq+=demod_signal_sum[b]*demod_signal_sum[b];
	}
for(b=window+search_window+1;b<2*window+1;b++){
	mean_power+=demod_signal_sum[b];
	mean_power_sq+=demod_signal_sum[b]*demod_signal_sum[b];
	}
	
mean_power/=2.0*(window-search_window);
mean_power_sq/=2.0*(window-search_window);
power_sd=sqrt((mean_power_sq-mean_power*mean_power)*2.0*(window-search_window)/(2.0*(window-search_window)-1));

if(debug) {
	fprintf(stderr, "mean_power=%g mean_power_sq=%g power_sd=%g\n", mean_power, mean_power_sq, power_sd);
	fprintf(stderr, "Power: ");
	for(b=0;b<2*window+1;b++) {
		fprintf(stderr, "%.2f ", (demod_signal_sum[b]-mean_power)/power_sd);
		}
	fprintf(stderr, "\n");
	}

cand->coherence_score=coherence_score;
cand->power_cor=power_cor;
cand->snr=(demod_signal_sum[b_max]-mean_power)/power_sd;
cand->snr=(demod_signal_sum[b_max]-fabs(demod_signal_sum[b_max-1]-demod_signal_sum[b_max+1])-mean_power)/power_sd;
if(demod_signal_sum[b_max]< mean_power)cand->strain=0.0;
	else
	cand->strain=2.0*sqrt(demod_signal_sum[b_max]-mean_power)/(1800.0*16384.0);

cand->f_max=cand->frequency+(b_max-window)/1800.0;

free(phase_sum);
free(signal_sum);
free(signal_sq_sum);
free(cov_sum);
free(demod_signal_sum);
}

//		fprintf(stderr, "%s step=%g  %g\n", __FUNCTION__, step, c.opt_char-cand->); \

//		fprintf(stderr, "New " #opt_var "=%f\n", cand->opt_var); \

#define BETTER_SNR(c1, c2)	(((c1).snr>(c2).snr) || (((c1).snr>=(c2).snr)) && ((c1).strain>(c2).strain))
#define BETTER_POWER_COR(c1, c2)	(((c1).snr>=(c2).snr) && ((c1).power_cor>(c2).power_cor))


#define BETTER_SNR_PC(c1, c2)	(((c1).snr+(c1).power_cor)>((c2).snr+(c2).power_cor))

#define BETTER_SNR_PC(c1, c2)	(((c1).snr>(c2).snr) || (((c1).snr>=(c2).snr) && ((c1).power_cor>(c2).power_cor)))

#define BETTER_SNR_PC(c1, c2)	(((c1).snr>(c2).snr+0.1) || (((c1).snr>=(c2).snr) && ((c1).strain>(c2).strain)))

#define BETTER_SNR_COH(c1, c2)	(((c1).snr>(c2).snr+0.1) || (((c1).snr>=(c2).snr)) && ((c1).coherence_score>(c2).coherence_score))

#define PLAIN_SNR(c1, c2)	((c1).snr>(c2).snr)

#define CHASE(opt_cond, opt_var, first_step) \
	int chase_##opt_var(CANDIDATE *cand)  \
	{  \
	CANDIDATE c;  \
	float fs; \
	float step;  \
	int status=0;  \
	int improv=0; \
	  \
	fs=first_step; \
	step=fs; \
	  \
	memcpy(&c, cand, sizeof(*cand));  \
	  \
	while(fabs(step)>fs/128.0) {  \
		c.opt_var+=step;  \
		compute_scores(&c, 0);  \
	  \
		if(opt_cond(c, *cand)) {  \
			c.frequency=c.f_max; \
			memcpy(cand, &c, sizeof(CANDIDATE));  \
			status=1;  \
			improv=1; \
			continue;  \
			} else {  \
			c.opt_var-=step;  \
			if(status==0) {  \
				step=-step;  \
				}  \
			status=0;  \
			step=step/2.0;  \
			}  \
		}  \
	if(improv) {  \
		return 1;  \
		} else {  \
		return 0;  \
		}  \
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

CHASE(BETTER_SNR_PC, psi, M_PI/16.0)
CHASE(BETTER_SNR_PC, iota, M_PI/16.0)
CHASE(BETTER_SNR_PC, ra, 2.0*resolution)
CHASE(BETTER_SNR_PC, dec, 2.0*resolution)
CHASE(BETTER_SNR_PC, spindown, 1.0/(1800.0*(max_gps()-min_gps())))
CHASE(BETTER_SNR_PC, frequency, 1/3600.0)

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
// ORACLE(snr, spindown, 3.9289e-10, 1)
// ORACLE(snr, frequency, 149.10508763, 1)

int chase_frequency1(CANDIDATE *cand)
{
CANDIDATE c;
float step=1/3600.0;
int status=0;
float var_start;

//fprintf(stderr, "%s %d\n", __FUNCTION__, __LINE__);

memcpy(&c, cand, sizeof(*cand));
var_start=cand->snr;

while(fabs(step)>1/720000.0) {
	c.frequency=c.f_max+step;
	compute_scores(&c, 0);

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


void optimize_candidate(CANDIDATE *cand)
{
CANDIDATE *tries;
int cont, i;

output_candidate_header(stderr);
/*compute_scores(cand, 1);*/
output_candidate(stderr, "", cand);
for(i=0;i<1000;i++) {
	cont=0;
	cont|=chase_psi(cand);
	cont|=chase_ra(cand);
	cont|=chase_dec(cand);
	cont|=chase_iota(cand);
	cont|=chase_frequency(cand);
	cont|=chase_spindown(cand);

	if(!cont)break;
/*	compute_scores(cand, 1);*/
	output_candidate(stderr, "", cand);
	}
compute_scores(cand, 1);
output_candidate(stderr, "", cand);
}

void identify_candidates(void)
{
RGBPic *p;
PLOT *plot;
SUM_TYPE a;
int i,k,m;
float *a_f;


candidate_free=0;
candidate_size=args_info.max_candidates_arg;
candidate=do_alloc(candidate_size, sizeof(*candidate));

max_dx=do_alloc(fine_grid->npoints, sizeof(*max_dx));
max_dx_local_map=do_alloc(fine_grid->npoints, sizeof(*max_dx_local_map));
max_dx_order=do_alloc(fine_grid->npoints, sizeof(*max_dx_order));

dec_order=do_alloc(fine_grid->npoints, sizeof(*dec_order));
inverse_dec_order=do_alloc(fine_grid->npoints, sizeof(*inverse_dec_order));

ra_order=do_alloc(fine_grid->npoints, sizeof(*ra_order));
inverse_ra_order=do_alloc(fine_grid->npoints, sizeof(*inverse_ra_order));

a_f=do_alloc(fine_grid->npoints, sizeof(*a_f));
polarization_index=do_alloc(fine_grid->npoints, sizeof(*polarization_index));

if(fine_grid->max_n_dec<800){
	p=make_RGBPic(fine_grid->max_n_ra*(800/fine_grid->max_n_dec)+140, fine_grid->max_n_dec*(800/fine_grid->max_n_dec));
	} else 
	p=make_RGBPic(fine_grid->max_n_ra+140, fine_grid->max_n_dec);	

plot=make_plot(p->width, p->height);

for(i=0;i<fine_grid->npoints;i++) {
	max_dx[i]=0;
	polarization_index[i]=-1;
	max_dx_order[i]=i;
	dec_order[i]=i;
	ra_order[i]=i;
	max_dx_local_map[i]=-1;
	}

for(i=0;i<fine_grid->npoints;i++) {
	for(k=0;k<ntotal_polarizations;k++){
		a=polarization_results[k].skymap.max_dx[i];
		if(a<0)continue;
		if(a>max_dx[i]) {
			max_dx[i]=a;
			polarization_index[i]=k;
			}
		}
	}
qsort(max_dx_order, fine_grid->npoints, sizeof(*max_dx_order), max_dx_compare);
qsort(dec_order, fine_grid->npoints, sizeof(*dec_order), dec_compare);
qsort(ra_order, fine_grid->npoints, sizeof(*ra_order), ra_compare);

for(i=0;i<fine_grid->npoints;i++) {
	inverse_dec_order[dec_order[i]]=i;
	inverse_ra_order[ra_order[i]]=i;
	}

/* Take noise floor to be bottom 1/3 of max_dx */
noise_floor=max_dx[max_dx_order[(fine_grid->npoints*2)/3]];
fprintf(LOG, "noise_floor: %f\n", noise_floor);
fprintf(stderr, "noise_floor: %f\n", noise_floor);
	
/* Find all local maxima, starting with largest.. */
for(i=0;i<fine_grid->npoints;i++) {
	k=max_dx_order[i];
	/* skip masked points.. */
	if(polarization_index[k]<0)continue;
	
	/* is the point marked already ? */
	if(max_dx_local_map[k]>=0){
		a=polarization_results[polarization_index[k]].skymap.max_dx[k];
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

output_candidate_header(LOG);

m=0;	
for(i=0;i<candidate_free;i++) {
	k=candidate[i].point_index;
	//fprintf(stderr, "k=%d\n", k);
	candidate[i].polarization_index=polarization_index[k];
	//fprintf(stderr, "polarization_index=%d\n", polarization_index[k]);
	candidate[i].rank=i;
	candidate[i].score=max_dx[k];
	if(candidate[i].score>noise_floor)m++;
	candidate[i].ra=fine_grid->longitude[k];
	candidate[i].dec=fine_grid->latitude[k];
	candidate[i].spindown=spindown;
	/* candidate[i].max_dx=polarization_results[polarization_index[k]].skymap.max_dx[k]; */
	candidate[i].psi=polarization_results[polarization_index[k]].orientation;
	/* currently we have either circular or linear polarizations */
	if(polarization_results[polarization_index[k]].cross_proj>0)candidate[i].iota=0.0;
		else	candidate[i].iota=M_PI/2.0;
	//candidate[i].a_plus=polarization_results[polarization_index[k]].plus_proj;
	//candidate[i].a_cross=polarization_results[polarization_index[k]].cross_proj;
	candidate[i].ul=polarization_results[polarization_index[k]].skymap.max_upper_limit[k];
	candidate[i].S_band=polarization_results[polarization_index[k]].skymap.S_map[k];
	candidate[i].M_band=polarization_results[polarization_index[k]].skymap.M_map[k];
	candidate[i].frequency=polarization_results[polarization_index[k]].skymap.freq_map[k];
	candidate[i].weight_ratio=1.0-polarization_results[polarization_index[k]].skymap.max_sub_weight[k]/polarization_results[polarization_index[k]].skymap.total_weight[k];
	candidate[i].skyband=fine_grid->band[k];
	
	//fprintf(stderr, "Computing scores\n");
	compute_scores(&(candidate[i]), 0);
	
	output_candidate(LOG, "_initial", &(candidate[i]));

	optimize_candidate(&(candidate[i]));

	output_candidate(LOG, "_optimized", &(candidate[i]));

	if(i<args_info.dump_candidates_arg)dump_candidate(i);
	}

fprintf(LOG, "high_candidates_count: %d\n", m);
fprintf(stderr, "high_candidates_count: %d\n", m);

fprintf(LOG, "candidates_count: %d\n", candidate_free);
fprintf(stderr, "candidates_count: %d\n", candidate_free);

snprintf(s, 19999, "%smax_dx.png", subinstance_name);
if(clear_name_png(s)){
	plot_grid_f(p, fine_grid, max_dx, 1);
	RGBPic_dump_png(s, p);
	}
snprintf(s, 19999, "%smax_dx.dat", subinstance_name);
dump_floats(s, max_dx, fine_grid->npoints, 1);

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
	a_f[i]=polarization_index[i];
	}
snprintf(s, 19999, "%smax_dx_polarization_index.png", subinstance_name);
if(clear_name_png(s)){
	plot_grid_f(p, fine_grid, a_f, 1);
	RGBPic_dump_png(s, p);
	}
snprintf(s, 19999, "%smax_dx_polarization_index.dat", subinstance_name);
dump_shorts(s, polarization_index, fine_grid->npoints, 1);


free(a_f);
a_f=NULL;
free_plot(plot);
free_RGBPic(p);

free(max_dx);
free(max_dx_local_map);
free(max_dx_order);
max_dx=NULL;
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
free(polarization_index);
candidate=NULL;
polarization_index=NULL;
}
