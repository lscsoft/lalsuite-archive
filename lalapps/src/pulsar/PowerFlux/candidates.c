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

#define DISTANCE_PRINTF		fprintf(stderr, "rank=%d distance=%f alignment=%f fdist=%f sdist=%f snr=%f strain=%g\n", cand->rank, candidate_distance(cand, args_info.focus_ra_arg, args_info.focus_dec_arg), candidate_alignment(cand, args_info.fake_iota_arg, args_info.fake_psi_arg), (cand->frequency-args_info.fake_freq_arg)*1800, (cand->spindown-args_info.fake_spindown_arg)*(max_gps()-min_gps())*1800.0, cand->snr, cand->strain);

#define WINDOW 240
#define SEARCH_WINDOW 5


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
float frequency;

if(polarization_index[index]<0)return 0;

frequency=polarization_results[polarization_index[index]].skymap.freq_map[index];

for(j=inverse_dec_order[index]+1;j<fine_grid->npoints;j++) {
	m=dec_order[j];

	if(fabs(fine_grid->latitude[m]-fine_grid->latitude[index])>fudge)break;

	if(max_dx_local_map[m]>=0)continue;
	
	if(max_dx[m]>max_dx[index])continue;
	
	if(polarization_index[m]<0)continue;

	if(((max_dx[m]+1.0)>max_dx[index]) && (fabs(frequency-polarization_results[polarization_index[m]].skymap.freq_map[m])>3/1800.0))continue;

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
fprintf(fout, "candidates: label polarization_index rank score point_index domain_size ul S M max_dx frequency psi iota ra dec spindown weight_ratio skyband coherence_score power_cor snr strain strain_err total_weight f_max ifo_freq ifo_freq_sd total\n");
}

void output_candidate(FILE *fout, char * suffix, CANDIDATE *cand)
{
fprintf(fout, "candidate%s: \"%s\" %d %d %g %d %d %g %g %g %f %f %f %f %f %f %g %f %d %f %f %f %g %g %f %f %f %f %d\n",
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

#define SNR_FORMULA SNR_FORMULA1

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
int window=WINDOW;
int search_window=SEARCH_WINDOW;
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

		signal_bin[k]=rintf(1800.0*f-first_bin);
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
		//if(d->expTMedians[k]*d->weight*response[k]<0.05)continue;

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
			#if 1
			signal_sum[b]+=power*weight;
			signal_sq_sum[b]+=power*power*weight;
			cov_sum[b]+=power*response[k]*weight;
			#endif

			demod_signal_sum[b]+=power*demod_weight;
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
	#endif
		
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

cand->f_max=cand->frequency+(b_max-window)/1800.0;

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
	fprintf(stderr, "total_demod_weight=%f mean_power=%g mean_power_sq=%g power_sd=%g\n", total_demod_weight, mean_power, mean_power_sq, power_sd);
	fprintf(stderr, "Power: ");
	for(b=0;b<2*window+1;b++) {
		fprintf(stderr, "%.2f ", (demod_signal_sum[b]-mean_power)/power_sd);
		}
	fprintf(stderr, "\n");
	}


b_max=window;

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

	int valid;
	} SCORE_AUX_DATA;

#define VALID_RESPONSE	1
#define VALID_SHIFT	2
#define VALID_DOPPLER	4
#define VALID_E		8

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
	//pl=&(d->polarizations[cand->polarization_index]);	
	pl=&(d->polarizations[0]);	
	
	response=&(ad->response[ad->offset[j]]);	
	/* process single SFTs */
	for(k=0;k<d->free;k++) {

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

		doppler[k]=ad->e[0]*d->detector_velocity[3*k+0]+
			ad->e[1]*d->detector_velocity[3*k+1]+
			ad->e[2]*d->detector_velocity[3*k+2];
		}

	}
ad->valid |= VALID_DOPPLER;
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

		f=frequency+frequency*doppler[k]+spindown*(d->gps[k]-spindown_start);

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
float frequency=cand->frequency;
int j, k, b, b_max;
double weight, f, total_demod_weight, mean_power, mean_power_sq, power_sd;
float *response;
float *doppler;
int signal_bin;
float demod_signal_sum[(2*WINDOW+1)*2], *power, *pout, demod_weight;

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
		/* skip SFTs with low weight */
		if(d->expTMedians[k]<0.05)continue;

		/* power_cor computation */
		weight=d->expTMedians[k]*d->weight;

		demod_weight=weight*response[k];
		total_demod_weight+=demod_weight*response[k];

		f=frequency+frequency*doppler[k]+spindown*(d->gps[k]-spindown_start);
		signal_bin=rintf(1800.0*f-first_bin);

		power=&(d->power[k*nbins+signal_bin-WINDOW]);
		pout=demod_signal_sum;

		for(b=0; b< (2*WINDOW+1); b++) {

			(*pout)+=(*power)*demod_weight;
			pout++;
			power++;
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

/* This is used to speedup optimization in iota and psi */

typedef struct {
	float f_plus_sq[(2*WINDOW+1)*2];
	float f_cross_sq[(2*WINDOW+1)*2];
	float f_plus_cross[(2*WINDOW+1)*2];
	} SNR_DATA;


static void precompute_snr_data(SNR_DATA *sd, SCORE_AUX_DATA *ad, CANDIDATE *cand, int debug)
{
DATASET *d;
float spindown=cand->spindown;
float frequency=cand->frequency;
int j, k, b, b_max;
double weight, f, mean_power, mean_power_sq, power_sd;
float *response;
float *doppler;
int signal_bin;
float *power, *pout, *cout, *pcout, pweight, cweight, pcweight, f_plus, f_cross;
POLARIZATION *pl;


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

	for(k=0;k<d->free;k++) {
		/* skip SFTs with low weight */
		if(d->expTMedians[k]<0.05)continue;

		f_plus=F_plus_coeff(k, ad->e, pl->AM_coeffs);
		f_cross=F_plus_coeff(k, ad->e, pl->conjugate->AM_coeffs);

		/* power_cor computation */
		weight=d->expTMedians[k]*d->weight;

		pweight=weight*f_plus*f_plus;
		cweight=weight*f_cross*f_cross;
		pcweight=weight*f_plus*f_cross;

		f=frequency+frequency*doppler[k]+spindown*(d->gps[k]-spindown_start);
		signal_bin=rintf(1800.0*f-first_bin);

		power=&(d->power[k*nbins+signal_bin-WINDOW]);

		pout=sd->f_plus_sq;
		cout=sd->f_cross_sq;
		pcout=sd->f_plus_cross;

		for(b=0; b< (2*WINDOW+1); b++) {

			(*pout)+=(*power)*pweight;
			(*cout)+=(*power)*cweight;
			(*pcout)+=(*power)*pcweight;

			pout++;
			cout++;
			pcout++;
			power++;
			}
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

		c1.psi=M_PI/2.0;
		c1.iota=M_PI/2.0;

		c2.psi=c1.psi;
		c2.iota=c1.iota;

		fill_response(ad, &c1);
		compute_simple_snr(ad, &c1, 0);
		compute_alignment_snr(&sd, &c2, 0);
		err=c1.snr-c2.snr;
		if(fabs(err)>1e-5)
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
float frequency=cand->frequency;
int i, j, k, b, b_max;
double weight, f, demod_weight, total_demod_weight, mean_power, mean_power_sq, power_sd;
float *response;
float *doppler;
int signal_bin, offset;
float demod_signal_sum[(2*WINDOW+1)*FREQ_STEPS], *power, *pout;

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
		/* skip SFTs with low weight */
		if(d->expTMedians[k]<0.05)continue;

		/* power_cor computation */
		weight=d->expTMedians[k]*d->weight;

		demod_weight=weight*response[k];
		total_demod_weight+=demod_weight*response[k];

		for(i=0;i<FREQ_STEPS;i++) {
			f=frequency+frequency*doppler[k]+spindown*(d->gps[k]-spindown_start)+i*(1.0/(1800.0*FREQ_STEPS));
			signal_bin=rintf(1800.0*f-first_bin);
			offset=i*(2*WINDOW+1);
	
			power=&(d->power[k*nbins+signal_bin-WINDOW]);
			pout=&(demod_signal_sum[offset]);
	
			for(b=0; b< (2*WINDOW+1); b++) {
				(*pout)+=*power*demod_weight;
				power++;
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
float frequency=cand->frequency;
int j, k, i, b;
float mismatch;
double weight, f, x, y, demod_weight, power[3], freq_hint_right[FREQ_STEPS], freq_hint_left[FREQ_STEPS], max_hint;
float *response;
float *doppler;
int signal_bin;
COMPLEX *p0;

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
		/* skip SFTs with low weight */
		if(d->expTMedians[k]<0.05)continue;

		/* power_cor computation */
		weight=d->expTMedians[k]*d->weight;

		demod_weight=weight*response[k];

		f=frequency+frequency*doppler[k]+spindown*(d->gps[k]-spindown_start);
		signal_bin=rint(1800.0*f-first_bin);
		mismatch=1800.0*f-first_bin-signal_bin;

		p0=&(d->bin[k*nbins+signal_bin-1]);

		for(b=0;b<=2; b++) {
			x=p0[b].re;
			y=p0[b].im;
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

#define SNR_TOLERANCE 0.01

//		fprintf(stderr, "%s step=%g  %g\n", __FUNCTION__, step, c.opt_char-cand->); \

//		fprintf(stderr, "New " #opt_var "=%f\n", cand->opt_var); \

#define BETTER_SNR(c1, c2)	(((c1).snr>(c2).snr) || (((c1).snr>=(c2).snr)) && ((c1).strain>(c2).strain))
#define BETTER_POWER_COR(c1, c2)	(((c1).snr>=(c2).snr) && ((c1).power_cor>(c2).power_cor))


#define BETTER_SNR_PC(c1, c2)	(((c1).snr+(c1).power_cor)>((c2).snr+(c2).power_cor))

#define BETTER_SNR_PC(c1, c2)	(((c1).snr>(c2).snr+0.1) || (((c1).snr>=(c2).snr) && ((c1).power_cor>(c2).power_cor)))

//#define BETTER_SNR_PC(c1, c2)	(((c1).snr>(c2).snr) || (((c1).snr>=(c2).snr) && ((c1).coherence_score>(c2).coherence_score)))

//#define BETTER_SNR_PC(c1, c2)	(((c1).snr>(c2).snr+0.1) || (((c1).snr>=(c2).snr) && ((c1).strain>(c2).strain)))

#define BETTER_SNR_COH(c1, c2)	(((c1).snr>(c2).snr+0.1) || (((c1).snr>=(c2).snr)) && ((c1).coherence_score>(c2).coherence_score))

#define PLAIN_SNR(c1, c2)	((c1).snr>(c2).snr)

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
	/* fprintf(stderr, "" # opt_var "=%f step=%f ", c.opt_var, step); */ \
	for(i=-N;i<=N;i++) { \
		c.opt_var=cand->opt_var+i*step; \
		fill_e(ad, &c); \
		fill_response(ad, &c); \
		fill_doppler(ad); \
		\
		compute_simple_snr(ad, &c, 0); \
		f=opt_expr; \
		/* fprintf(stderr, ",%f", f); */ \
		if(max_i<-N || (f>max) || (f>=max && !max_i)) { \
			max_i=i; \
			max=f; \
			} \
		} \
	free_score_aux_data(ad); \
	/* fprintf(stderr, "\n"); */ \
	if(max<cand->snr+0.001)return 0; \
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
CHASE(PLAIN_SNR, ra, 0.5*resolution)
CHASE(PLAIN_SNR, dec, 0.5*resolution)
CHASE(PLAIN_SNR, spindown, 0.5/(1800.0*(max_gps()-min_gps())))
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
while(c>M_PI)c-=M_PI;
while(c<0)c+=M_PI;

dist=fabs(a2-b2)+0.5*fmin((1.0-a)*(1.0-a), (1.0-b)*(1.0-b))*c;

return(dist);
//return(acos(fabs(sin(cand->iota)*sin(iota)*cos(cand->psi-psi)+cos(cand->iota)*cos(iota))-1e-14));
}

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
float freq_dir, psi_dir, iota_dir, ra_dir, dec_dir, spindown_dir;

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
int i, N=8;
float max;
SCORE_AUX_DATA *ad;


ad=allocate_score_aux_data();


memcpy(&c, cand, sizeof(*cand)); 

fill_e(ad, &c);
fill_response(ad, &c);
fill_doppler(ad);

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

	c.ra=cand->ra+i*ra_dir*resolution*0.5*0.125/(cos(cand->dec)+0.001);
	c.dec=cand->dec+i*dec_dir*resolution*0.5*0.125;

	fill_e(ad, &c);
	fill_response(ad, &c);
	fill_doppler(ad);

	for(spindown_dir=-1; spindown_dir<=1;spindown_dir++) {
		c.spindown=cand->spindown+i*spindown_dir*1e-11;
		c.frequency=cand->f_max;

	for(frequency_dir=-5; frequency_dir<5; frequency_dir++) {

		c.frequency=cand->frequency+frequency_dir*0.05/1800.0;
	
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

compute_snr(ad, cand, 0);

var_start=cand->snr;

while(fabs(step)>1/720000.0) {
	c.frequency=cand->f_max+step;

	compute_snr(ad, &c, 0);

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
				fprintf(stderr, "found snr=%f\n", c.snr);
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
fprintf(stderr, "search_alignment start snr=%f\n", best_c.snr);

for(i=-N;i<N;i++) {
	c.iota=cand->iota+i*M_PI/128.0;
	for(j=-N;j<N;j++) {
		c.psi=cand->psi+j*M_PI/128.0;
			compute_scores(&c, 0);
/*			fprintf(stderr, "% 2d", (int)floor(10.0*c.snr/best_c.snr));*/
			if(c.snr>best_c.snr) {
				memcpy(&best_c, &c, sizeof(CANDIDATE));
				fprintf(stderr, "found snr=%f\n", c.snr);
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

precompute_snr_data(&sd, ad, &best_c, 0);

compute_alignment_snr(&sd, &best_c, 0);
fprintf(stderr, "search_alignment start snr=%f\n", best_c.snr);

while(1) {
	found=0;

	c.psi=best_c.psi;
	for(i=1;i<N;i++) {
		c.iota=best_c.iota+i*M_PI/128.0;

		compute_alignment_snr(&sd, &c, 0);
		if(c.snr>best_c.snr) {
			memcpy(&best_c, &c, sizeof(CANDIDATE));
			fprintf(stderr, "found snr=%f\n", c.snr);
			found=1;
			break;
			}
		c.iota=best_c.iota-i*M_PI/128.0;
		compute_alignment_snr(&sd, &c, 0);
		if(c.snr>best_c.snr) {
			memcpy(&best_c, &c, sizeof(CANDIDATE));
			fprintf(stderr, "found snr=%f\n", c.snr);
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
			fprintf(stderr, "found snr=%f\n", c.snr);
			found=1;
			break;
			}
		c.psi=best_c.psi-j*M_PI/128.0;
		compute_alignment_snr(&sd, &c, 0);
		if(c.snr>best_c.snr) {
			memcpy(&best_c, &c, sizeof(CANDIDATE));
			fprintf(stderr, "found snr=%f\n", c.snr);
			found=1;
			break;
			}
		}
	if(!found)break;
	}

free_score_aux_data(ad);

fprintf(stderr, "search_alignment stop snr=%f\n", best_c.snr);
if(fabs(best_c.iota-cand->iota)>0 || fabs(best_c.psi-cand->psi)>0) {
	memcpy(cand, &best_c, sizeof(best_c));
	return 1;
	}
return 0;
}


int search_spindown1(CANDIDATE *cand)
{
int i,j, k, N=40;
CANDIDATE c, best_c;
SCORE_AUX_DATA *ad;

memcpy(&c, cand, sizeof(CANDIDATE));
memcpy(&best_c, cand, sizeof(CANDIDATE));

ad=allocate_score_aux_data();

fill_e(ad, &c);
fill_response(ad, &c);
fill_doppler(ad);

compute_snr(ad, &best_c, 0);
fprintf(stderr, "spindown search start snr=%f ", best_c.snr);

for(i=-N;i<=N;i++) {
	c.spindown=cand->spindown+i*1e-11;
	for(j=0;j<=3;j++) {
		c.frequency=cand->f_max-i*1e-11*(max_gps()-min_gps())*0.5+j*0.25/1800.0;
		//c.frequency=c.f_max;
		//chase_frequency1(&c);
		compute_snr(ad, &c, 0);
		fprintf(stderr, ",%f", c.snr);
		if(c.snr>best_c.snr-3) {
			for(k=0;k<=7;k++) {
				c.frequency=cand->f_max-i*1e-11*(max_gps()-min_gps())*0.5+(j+k*0.125)*0.25/1800.0;
				compute_snr(ad, &c, 0);
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
							fprintf(stderr, "found snr=%f\n", c.snr);
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
CANDIDATE *tries;
int cont, i;
int frequency_dir, psi_dir, iota_dir, ra_dir, dec_dir, spindown_dir;
float factor;

//#define DISTANCE_PRINTF		fprintf(stderr, "distance=%f fdist=%f sdist=%f snr=%f\n", candidate_distance(cand, args_info.focus_ra_arg, args_info.focus_dec_arg), (cand->frequency-140.5972)*1800, (cand->spindown-3.794117e-10)*(max_gps()-min_gps())*1800.0, cand->snr);



//cand->ra=6.112886;
//cand->dec=-0.6992048;


output_candidate_header(stderr);
compute_scores(cand, 1);
DISTANCE_PRINTF
output_candidate(stderr, "", cand);

if(cand->rank==0)test_alignment_snr(cand);

/*cand->ra=args_info.fake_ra_arg;
cand->dec=args_info.fake_dec_arg;*/
/*cand->iota=0.0;
cand->psi=0.0;
cand->spindown=0.0;*/
// cand->frequency=200.1;

/*for(i=0;i<1000;i++) {
	cont=0;
	cont|=chase_iota(cand);
	fprintf(stderr, "distance=%f snr=%f\n", candidate_distance(cand, 6.112886, -0.6992048), cand->snr);
	cont|=chase_psi(cand);
	fprintf(stderr, "distance=%f snr=%f\n", candidate_distance(cand, 6.112886, -0.6992048), cand->snr);
	cont|=chase_ra(cand);
	fprintf(stderr, "distance=%f snr=%f\n", candidate_distance(cand, 6.112886, -0.6992048), cand->snr);
	cont|=chase_dec(cand);
	fprintf(stderr, "distance=%f snr=%f\n", candidate_distance(cand, 6.112886, -0.6992048), cand->snr);
	cont|=chase_frequency(cand);
	fprintf(stderr, "distance=%f snr=%f\n", candidate_distance(cand, 6.112886, -0.6992048), cand->snr);
	//cont|=chase_spindown(cand);

	if(!cont)break;
	compute_scores(cand, 1);
	output_candidate(stderr, "", cand);
	}*/

// tabulate_neighbourhood("cand_0.txt", cand);
// exit(0);

factor=1.0;

for(i=0;i<1000;i++) {
	cont=0;
/*	cont=search_ra(cand, resolution*0.5/(cos(cand->dec)+0.001));
	compute_scores(cand, 0);
	fprintf(stderr, "distance=%f fdist=%f snr=%f\n", candidate_distance(cand, 6.112886, -0.6992048), (cand->frequency-140.5285)*1800, cand->snr);

	cont|=search_dec(cand, resolution*0.5);*/
/*	compute_scores(cand, 0);*/
// 	compute_scores(cand, 0);
// 	DISTANCE_PRINTF

// 	cont=chase_ra(cand);
// /*	compute_scores(cand, 0);*/
// 	DISTANCE_PRINTF
// 
// 	cont|=chase_dec(cand);
// /*	compute_scores(cand, 0);*/
// 	DISTANCE_PRINTF

// /*	while(search_iota(cand, M_PI/256.0))cont|=1;
// /*	cont|=chase_iota(cand);*/
// /*	compute_scores(cand, 0);*/
// 	compute_scores(cand, 0);
// 	DISTANCE_PRINTF
// 	output_candidate(stderr, "", cand);*/
// 	
// /*	cont|=search_psi(cand, M_PI/256.0);*/
// 	cont|=chase_psi(cand);
// /*	compute_scores(cand, 0);*/
// 	DISTANCE_PRINTF
// 
/*	cont|=chase_frequency1(cand);
	compute_scores(cand, 0);
	DISTANCE_PRINTF
	output_candidate(stderr, "", cand);*/
// 
//  	cont|=chase_spindown(cand);
// /* 	compute_scores(cand, 0);*/
//  	DISTANCE_PRINTF

// 	cont|=fit_psi(cand, M_PI/64.0);
//  	DISTANCE_PRINTF
// 
// 	cont+=fit_ra(cand, resolution*0.125/(cos(cand->dec)+0.001));
// 	compute_scores(cand, 0);
//  	DISTANCE_PRINTF
// 	output_candidate(stderr, "", cand);
// 
// 	cont+=fit_dec(cand, resolution*0.125);
// 	compute_scores(cand, 0);
//  	DISTANCE_PRINTF
// 	output_candidate(stderr, "", cand);

// 	cont|=search_monte_carlo(cand, resolution*5, 2e-10, 1000);

	cont+=search_frequency(cand, 0.02/1800.0, 256);
/*	compute_scores(cand, 0);*/
 	DISTANCE_PRINTF
	output_candidate(stderr, "", cand);

	cont+=search_alignment(cand);

/*	cont+=search_psi(cand, factor*M_PI/512.0, 128);
	compute_scores(cand, 0);
 	DISTANCE_PRINTF
	output_candidate(stderr, "", cand);

	cont+=search_iota(cand, factor*M_PI/128.0, 64);
	compute_scores(cand, 0);
	DISTANCE_PRINTF
	output_candidate(stderr, "", cand);*/
/*
	cont+=search_ra(cand, factor*0.5*resolution*0.125/(cos(cand->dec)+0.001), 256);
	compute_scores(cand, 0);
 	DISTANCE_PRINTF
	output_candidate(stderr, "", cand);

	cont+=search_dec(cand, factor*0.5*resolution*0.125, 256);
	compute_scores(cand, 0);
 	DISTANCE_PRINTF
	output_candidate(stderr, "", cand);
*/
// 	cont+=search_spindown(cand, factor*1e-11, 128);
// 	compute_scores(cand, 0);
//  	DISTANCE_PRINTF
// 	output_candidate(stderr, "", cand);

// 	for(ra_dir=0; ra_dir<=1;ra_dir++) {
// 	for(dec_dir=-1; dec_dir<=1;dec_dir++) {
// 	for(spindown_dir=-1; spindown_dir<=1;spindown_dir++) {
// 	for(frequency_dir=-1; frequency_dir<=1; frequency_dir++) {
// 		fprintf(stderr, "%d %d %d %d\n", ra_dir, dec_dir, spindown_dir, frequency_dir);
// 
// 		if(search_vec(cand, 
// 				frequency_dir*0.02/1800.0,
// 				0.0, 
// 				0.0, 
// 				ra_dir*resolution*0.5*0.125/(cos(cand->dec)+0.001), 
// 				dec_dir*resolution*0.5*0.125, 
// 				spindown_dir*1e-11, 
// 				128)) {
// 			cont+=1;
// 			compute_scores(cand, 0);
// 			DISTANCE_PRINTF
// 			output_candidate(stderr, "", cand);
// 			}
// 
// 		}}}}

	cont+=search_four_vec(cand);

//    	cont|=search_monte_carlo_vec(cand);

// 	if(!cont && (cand->f_max!=cand->frequency)) {
// 		cand->frequency=cand->f_max;
// 		cont=1;
// 		}

//  	if(!cont)cont|=search_spindown(cand);

//	if(!cont)cont|=search_monte_carlo(cand, resolution*5, 2e-10);

/* 	if(!cont)cont|=search_four(cand);*/
// 

// 	if(!cont) { 
// 		cont|=search_sky(cand);
// 		compute_scores(cand, 0);
// 		DISTANCE_PRINTF
// 		}

	if(!cont) {
		factor=factor*0.5;
		//cont= factor > 1e-6;
		} else 
	if( (cont>2) && (factor<1.0) ) {
		//factor=factor*2.0;
		}

	if(!cont)break;
/*	compute_scores(cand, 1);*/
/*	compute_scores(cand, 0);*/
/*	output_candidate(stderr, "", cand);*/
	//exit(0);
	}

compute_scores(cand, 1);
fprintf(LOG, "distance=%f snr=%f\n", candidate_distance(cand, args_info.focus_ra_arg, args_info.focus_dec_arg), cand->snr);
DISTANCE_PRINTF
output_candidate(stderr, "", cand);

// tabulate_neighbourhood(cand);
// exit(0);

}

void identify_candidates(void)
{
RGBPic *p;
PLOT *plot;
SUM_TYPE a;
int i,k,m;
float *a_f;

/* compute power without applying whitening procedure */
recompute_power();

candidate_free=0;
candidate_size=100;
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
	if(polarization_index[k]<0)continue;
	
	/* is the point marked already ? */
	if(max_dx_local_map[k]>=0) {
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

/* Now do detailed processing of the candidates */

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
	if(polarization_results[polarization_index[k]].cross_factor>0)candidate[i].iota=0.0;
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

	if(i<args_info.max_candidates_arg) {
		optimize_candidate(&(candidate[i]));
		output_candidate(LOG, "_optimized", &(candidate[i]));
		}

	if(i<args_info.dump_candidates_arg)dump_candidate(i);
	}

fprintf(LOG, "high_candidates_count: %d\n", m);
fprintf(stderr, "high_candidates_count: %d\n", m);

fprintf(LOG, "candidates_count: %d\n", candidate_free);
fprintf(stderr, "candidates_count: %d\n", candidate_free);



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
