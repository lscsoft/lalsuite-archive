#define _GNU_SOURCE
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <alloca.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_integration.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <fcntl.h>

#include "dataset.h"
#include "statistics.h"
#include "cmdline.h"
#include "intervals.h"
#include "hookup.h"
#include "grid.h"
#include "rastermagic.h"

#include <lal/GeneratePulsarSignal.h>

// typedef double REAL8;
// typedef int INT4;
// typedef float REAL4;

extern FILE *LOG;
extern struct gengetopt_args_info args_info;
extern int first_bin;
extern int nbins;

extern int ntotal_polarizations;
extern SKY_GRID *fine_grid;
extern SKY_GRID *patch_grid;

extern INT64 spindown_start;

extern int fake_injection;

char s[20000];

DATASET * datasets=NULL;
int d_free=0;
int d_size=0;

int lock_file=-1;

/* this is a busy wait.. as condor does not allow sleep() */
static void spin_wait(int seconds)
{
time_t start_time, end_time;
int i;
time(&start_time);
i=0;
while(1) {
	i++;
	if(i>1000000) {
		time(&end_time);
		if(end_time>start_time+seconds)return;
		i=0;
		}
	}
}

static void acquire_lock(char *filename)
{
struct flock fl;
int i;

lock_file=open(filename, O_CREAT | O_RDWR, 0666);
if(lock_file<0) {
	fprintf(stderr, "Could not open file \"%s\" for locking\n", filename);
	return;
	}
fl.l_type=F_WRLCK;
fl.l_whence=SEEK_SET;
fl.l_start=0;
fl.l_len=1;
errno=0;
i=0;
while(fcntl(lock_file, F_SETLK, &fl)<0){
	if(i>100) {
		fprintf(stderr, "Waiting for lock: %s\n", strerror(errno));
		i=0;
		}
	spin_wait(10);
	i++;
	}
}

static void release_lock(void)
{
if(lock_file>=0)close(lock_file);
lock_file=-1;
}

typedef struct {
	float e[26];
	int bin;

	double a_plus;
	double a_cross;
	double cos_e;
	double sin_e;

	double segment_start;

	double ref_time;

	double freq;
	double spindown;

	double ra;
	double dec;

	} SIGNAL_PARAMS;

static void compute_signal(double *re, double *im, double *f, double t, SIGNAL_PARAMS *p)
{
float det_vel[3];
float f_plus, f_cross;
double doppler, freq, omega_t;
double response;

get_detector_vel(round(t), det_vel);

doppler=p->e[0]*det_vel[0]+p->e[1]*det_vel[1]+p->e[2]*det_vel[2];

get_AM_response(round(t), p->dec, p->ra, 0.0, &f_plus, &f_cross);

*f=p->freq*(1.0+doppler)+p->spindown*(t-p->ref_time);

*f=p->freq*(1.0+doppler);

/* this loses the phase coherence between segments, but avoids the need to accumulate
  phase from the start of the run */
omega_t=2.0*M_PI*(*f-p->bin/1800.0)*(t-p->segment_start);

// response=f_plus*(p->a_plus*cos(omega_t)*p->cos_e-p->a_cross*sin(omega_t)*p->sin_e)+
// 	f_cross*(p->a_plus*cos(omega_t)*p->sin_e-p->a_cross*sin(omega_t)*p->cos_e);

*re=f_plus*(p->a_plus*cos(omega_t)*p->cos_e-p->a_cross*sin(omega_t)*p->sin_e)+
	f_cross*(p->a_plus*cos(omega_t)*p->sin_e+p->a_cross*sin(omega_t)*p->cos_e);

*im=f_plus*(p->a_plus*sin(omega_t)*p->cos_e+p->a_cross*cos(omega_t)*p->sin_e)+
	f_cross*(p->a_plus*sin(omega_t)*p->sin_e-p->a_cross*cos(omega_t)*p->cos_e);

//fprintf(stderr, "response=%g a_plus=%g a_cross=%g \n", response, p->a_plus, p->a_cross);

//fprintf(stderr, "f=%.4f (%.3f) re=%.2f im=%.2f omega_t=%.2f t=%.1f\n", (*f), (*f)*1800.0-p->bin*1.0, *re, *im, omega_t, t-p->ref_time);

}

static double signal_re(double t, void *params)
{
double re, im, f;
compute_signal(&re, &im, &f, t, params);
return(re);
}

static double signal_im(double t, void *params)
{
double re, im, f;
compute_signal(&re, &im, &f, t, params);
return(im);
}

static void inject_fake_signal(DATASET *d, int segment)
{
double result, abserr;
size_t neval;
int err;
int window=5, bin;
int i;
double re, im, f, cos_i;
double a,b;
SIGNAL_PARAMS p;
gsl_function F; 

F.params=&p;
p.bin=0;

p.ra=args_info.fake_ra_arg;
p.dec=args_info.fake_dec_arg;
p.freq=args_info.fake_freq_arg;
p.spindown=args_info.fake_spindown_arg;
p.ref_time=args_info.fake_ref_time_arg;
p.segment_start=d->gps[segment];

cos_i=cos(args_info.fake_iota_arg);

p.a_plus=(1+cos_i*cos_i)*0.5;
p.a_cross=cos_i;

p.cos_e=cos(2.0*(args_info.fake_psi_arg-args_info.fake_phi_arg));
p.sin_e=sin(2.0*(args_info.fake_psi_arg-args_info.fake_phi_arg));

precompute_am_constants(p.e, args_info.fake_ra_arg, args_info.fake_dec_arg);
compute_signal(&re, &im, &f, d->gps[segment]+(int)(d->coherence_time/2), &p);

bin=round(f*1800.0-d->first_bin);

if(bin+window>nbins)bin=nbins-window;
if(bin<window)bin=window;

for(i=bin-window; i<=bin+window; i++) {
	F.function=signal_re;
	p.bin=i+d->first_bin;

	err=gsl_integration_qng(&F, d->gps[segment], d->gps[segment]+d->coherence_time,
		1, 1e-2,
		&result, &abserr, &neval);
/*	fprintf(stderr, "re %d %d result=%g abserr=%g %d %s\n", segment, i, result, abserr, neval, gsl_strerror(err)); */
	d->bin[segment*d->nbins+i].re+=args_info.fake_strain_arg*result*16384.0;


	F.function=signal_im;
	err=gsl_integration_qng(&F, d->gps[segment], d->gps[segment]+d->coherence_time,
		1, 1e-2,
		&result, &abserr, &neval);
/*	fprintf(stderr, "im %d %d result=%g abserr=%g %d %s\n", segment, i, result, abserr, neval, gsl_strerror(err)); */
	d->bin[segment*d->nbins+i].im+=args_info.fake_strain_arg*result*16384.0;


	}

}

static void inject_fake_signal2(DATASET *d, int segment)
{
SFTandSignalParams sft_params;
PulsarSignalParams pulsar_params;

sft_params.resTrig=0;
sft_params.Dterms=5;
sft_params.nSamples=5;

}

static void init_dataset(DATASET *d)
{
memset(d, 0, sizeof(*d));
d->name="unknown";
d->lock_file=NULL;
d->validated=0;
d->detector="unknown";
d->gps_start=0;
d->gps_stop=-1;

d->coherence_time=1800.0;
d->nbins=nbins;
d->first_bin=first_bin;

d->segment_list=NULL;
d->veto_segment_list=NULL;

d->size=1000;
d->free=0;
d->gps=do_alloc(d->size, sizeof(*d->gps));
d->bin=do_alloc(d->size*d->nbins, sizeof(*d->bin));
d->power=NULL;

d->apply_hanning_filter=0;
d->weight=1.0;

d->polarizations=NULL;
d->polarizations=do_alloc(ntotal_polarizations, sizeof(*(d->polarizations)));
memset(d->polarizations, 0, (ntotal_polarizations)*sizeof(*(d->polarizations)));
}

static void compute_detector_speed(DATASET *d)
{
int i;
float det_vel_ra, det_vel_dec;
double orbital_axis[3];
float band_axis_ra, band_axis_dec;
fprintf(stderr,"Computing detector speed for dataset %s\n", d->name);
d->detector_velocity=do_alloc(3*d->free, sizeof(*d->detector_velocity));
d->average_detector_velocity[0]=0.0;
d->average_detector_velocity[1]=0.0;
d->average_detector_velocity[2]=0.0;
for(i=0;i<d->free;i++){
	/* middle of the 30min interval */
	get_detector_vel(d->gps[i]+(int)(d->coherence_time/2),&(d->detector_velocity[3*i]));
		
	d->average_detector_velocity[0]+=d->detector_velocity[3*i];
	d->average_detector_velocity[1]+=d->detector_velocity[3*i+1];
	d->average_detector_velocity[2]+=d->detector_velocity[3*i+2];
	}
d->average_detector_velocity[0]/=d->free;
d->average_detector_velocity[1]/=d->free;
d->average_detector_velocity[2]/=d->free;
fprintf(LOG,"dataset %s average detector velocity: %g %g %g\n", d->name, 
	d->average_detector_velocity[0],
	d->average_detector_velocity[1],
	d->average_detector_velocity[2]);
det_vel_dec=atan2f(d->average_detector_velocity[2], 
	sqrt(d->average_detector_velocity[0]*d->average_detector_velocity[0]+
		d->average_detector_velocity[1]*d->average_detector_velocity[1]));
det_vel_ra=atan2f(d->average_detector_velocity[1], d->average_detector_velocity[0]);
if(det_vel_ra<0)det_vel_ra+=2.0*M_PI;
fprintf(LOG,"dataset %s average detector velocity RA (degrees) : %f\n", d->name, det_vel_ra*180.0/M_PI);
fprintf(LOG,"dataset %s average detector velocity DEC (degrees): %f\n", d->name, det_vel_dec*180.0/M_PI);
fprintf(stderr,"dataset %s average detector velocity RA (degrees) : %f\n", d->name, det_vel_ra*180.0/M_PI);
fprintf(stderr,"dataset %s average detector velocity DEC (degrees): %f\n", d->name, det_vel_dec*180.0/M_PI);

orbital_axis[0]=0.0;
orbital_axis[1]=-sin(M_PI*23.44/180.0);
orbital_axis[2]=cos(M_PI*23.44/180.0);

/* crossproduct gives the vector perpedicular to both the average doppler shift and
  orbital axis */
d->band_axis[0]=orbital_axis[1]*d->average_detector_velocity[2]-orbital_axis[2]*d->average_detector_velocity[1];
d->band_axis[1]=orbital_axis[2]*d->average_detector_velocity[0]-orbital_axis[0]*d->average_detector_velocity[2];
d->band_axis[2]=orbital_axis[0]*d->average_detector_velocity[1]-orbital_axis[1]*d->average_detector_velocity[0];

/* Normalize */
d->band_axis_norm=sqrt(d->band_axis[0]*d->band_axis[0]+
	d->band_axis[1]*d->band_axis[1]+
	d->band_axis[2]*d->band_axis[2]);
/* replace 0.0 with something more reasonable later */
if(d->band_axis_norm<=0.0){
	d->band_axis[0]=0.0;
	d->band_axis[1]=0.0;
	d->band_axis[2]=1.0;
	} else {
	d->band_axis[0]/=d->band_axis_norm;
	d->band_axis[1]/=d->band_axis_norm;
	d->band_axis[2]/=d->band_axis_norm;
	}
/* 
  Normalize band_axis_norm so it matches the definition of \vec{u}
*/
d->band_axis_norm*=2.0*M_PI/(365.0*24.0*3600.0);

fprintf(LOG, "dataset %s auto band axis norm: %g\n", d->name, d->band_axis_norm);
fprintf(LOG, "dataset %s maximum S contribution from Doppler shifts: %g\n", d->name, d->band_axis_norm*(d->first_bin+d->nbins*0.5)/d->coherence_time);

if(args_info.band_axis_norm_given) {
	d->band_axis_norm=args_info.band_axis_norm_arg;
	}

fprintf(LOG, "dataset %s actual band axis norm: %g\n", d->name, d->band_axis_norm);

d->large_S=6.0/(d->coherence_time*(d->gps[d->free-1]-d->gps[0]+d->coherence_time));

fprintf(LOG, "dataset %s auto large S: %g\n", d->name, d->large_S);

if(args_info.large_S_given){
	d->large_S=args_info.large_S_arg;
	}
fprintf(LOG, "dataset %s large S: %g\n", d->name, d->large_S);
	

fprintf(LOG,"dataset %s auto band axis: %g %g %g\n", d->name, 
	d->band_axis[0],
	d->band_axis[1],
	d->band_axis[2]);
fprintf(stderr,"dataset %s auto band axis: %g %g %g\n", d->name, 
	d->band_axis[0],
	d->band_axis[1],
	d->band_axis[2]);

fprintf(LOG, "dataset %s band_axis: %s\n", d->name, args_info.band_axis_arg);
fprintf(stderr, "dataset %s band_axis: %s\n", d->name, args_info.band_axis_arg);
if(!strcasecmp(args_info.band_axis_arg, "equatorial")){
	d->band_axis[0]=0.0;
	d->band_axis[1]=0.0;
	d->band_axis[2]=1.0;	
	} else
if(!strncasecmp(args_info.band_axis_arg, "explicit", 8)){
	int q;
	q=sscanf(args_info.band_axis_arg+8, "(%lf,%lf,%lf)", 
		&(d->band_axis[0]),
		&(d->band_axis[1]),
		&(d->band_axis[2]));
	if(q!=3){
		fprintf(stderr,"Warning ! In dataset %s not all explicit band axis values were assigned. Format error ?\n", d->name);
		fprintf(LOG,"Warning ! In dataset %s not all explicit band axis values were assigned. Format error ?\n", d->name);
		}
	}
	
fprintf(LOG,"dataset %s actual band axis: %g %g %g\n", d->name, 
	d->band_axis[0],
	d->band_axis[1],
	d->band_axis[2]);
fprintf(stderr,"dataset %s actual band axis: %g %g %g\n", d->name, 
	d->band_axis[0],
	d->band_axis[1],
	d->band_axis[2]);

band_axis_dec=atan2f(d->band_axis[2], 
	sqrt(d->band_axis[0]*d->band_axis[0]+d->band_axis[1]*d->band_axis[1]));
band_axis_ra=atan2f(d->band_axis[1], d->band_axis[0]);

if(band_axis_ra<0)band_axis_ra+=2.0*M_PI;
fprintf(stderr,"dataset %s band axis RA (degrees) : %f\n", d->name, band_axis_ra*180.0/M_PI);
fprintf(stderr,"dataset %s band axis DEC (degrees): %f\n", d->name, band_axis_dec*180.0/M_PI);
fprintf(LOG,"dataset %s band axis RA (degrees) : %f\n", d->name, band_axis_ra*180.0/M_PI);
fprintf(LOG,"dataset %s band axis DEC (degrees): %f\n", d->name, band_axis_dec*180.0/M_PI);
}

/* Note: this messes with tm variable,
  tm[i]=exp(log(10)*TMedians[i])/Modulation[i] */
float FindCutOff(float *tm, int nsegments)
{
int i;
double sum,sum_squared,mean,sigma;
double best_cutoff,smallest_sigma;
double a;
int best_i;
sort_floats(tm, nsegments);
sum=0;
sum_squared=0;
best_i=0;
smallest_sigma=10;
best_cutoff=tm[0];
for(i=0;i<nsegments;i++) {
	a=tm[i]/tm[0];
	sum+=a;
	sum_squared+=a*a;
	mean=sum/(i+1);
	sigma=10*sqrt((sum_squared))/(i+1);
	
	if(sigma<smallest_sigma) {
		smallest_sigma=sigma;
		best_i=i;
		best_cutoff=tm[i];
		}
	}
//fprintf(stderr,"Cutoff: i=%d sigma=%g cutoff=%g (%g,%g,..)\n",best_i, smallest_sigma, best_cutoff,tm[0],tm[1]);
return best_cutoff;
}

void apply_hanning_filter(DATASET *d)
{
int i,k;
COMPLEX *tmp, *p;
tmp=do_alloc(d->nbins, sizeof(*tmp));
for(i=0;i<d->free;i++) {
	tmp[0].re=0;
	tmp[0].im=0;
	tmp[d->nbins-1].re=0;
	tmp[d->nbins-1].im=0;
	p=&(d->bin[i*d->nbins]);
	for(k=1;k<(d->nbins-1);k++) {
		tmp[k].re=p[k].re-(p[k-1].re+p[k+1].re)*0.5;
		tmp[k].im=p[k].im-(p[k-1].im+p[k+1].im)*0.5;
		}
	memcpy(p, tmp, nbins*sizeof(*p));	
	}
free(tmp);
}

typedef struct {
	int index;
	INT64 gps;
	} ELEMENT;

int element_cmp(ELEMENT *e1, ELEMENT *e2)
{
if(e1->gps<e2->gps)return(-1);
if(e1->gps>e2->gps)return(1);
return 0;
}

void sort_dataset(DATASET *d)
{
ELEMENT *p;
COMPLEX *bin;
INT64 *gps;
int i, k;

p=do_alloc(d->free, sizeof(*p));
for(i=0;i<d->free;i++) {
	p[i].index=i;
	p[i].gps=d->gps[i];
	}

qsort(p, d->free, sizeof(*p), element_cmp);

bin=do_alloc(d->free*d->nbins, sizeof(*bin));
gps=do_alloc(d->free, sizeof(*gps));
for(i=0;i<d->free;i++) {
	k=p[i].index;
	gps[i]=d->gps[k];
	memcpy(&(bin[i*d->nbins]), &(d->bin[k*d->nbins]), d->nbins*sizeof(*bin));
	}

d->size=d->free;
free(d->bin);
free(d->gps);
d->bin=bin;
d->gps=gps;

free(p);
}

static int validate_dataset(DATASET *d)
{
int i,j, m;
float x,y, *tm;
if(d->validated)return 1;
fprintf(stderr, "Validating dataset \"%s\"\n", d->name);
fprintf(stderr, "Found %d segments\n", d->free);

fprintf(LOG, "dataset %s nsegments : %d\n", d->name, d->free);
fprintf(LOG, "dataset %s detector : %s\n", d->name, d->detector);

if(d->free<1){
	/* free large chunks of ram that might have been reserved */
	if(d->size>0) {
		free(d->gps);
		free(d->bin);
		d->gps=NULL;
		d->bin=NULL;
		}
	if(d->power!=NULL){
		free(d->power);
		d->power=NULL;
		}
	return 0;
	}

/* TODO: fixme 
fprintf(LOG, "dataset %s first gps : %lld\n", d->name, d->gps[0]);
fprintf(LOG, "dataset %s last gps  : %lld\n", d->name, d->gps[d->free-1]);
*/

if(!strcmp(d->detector, "unknown")) {
	fprintf(stderr, "Each dataset must specify detector being used\n");
	return 0;
	}

get_detector(d->detector);

sort_dataset(d);

if(d->apply_hanning_filter)apply_hanning_filter(d);

if(fake_injection) {
	fprintf(stderr, "Injecting fake signal.\n");
	fprintf(LOG, "Injecting fake signal.\n");
	for(i=0;i<d->free;i++) {
		inject_fake_signal(d, i);
		}
	}

d->power=do_alloc(d->free*d->nbins, sizeof(*d->power));
for(i=0;i<d->free*d->nbins;i++){
	x=d->bin[i].re;
	y=d->bin[i].im;
	d->power[i]=x*x+y*y;
	}

/* compute time from the start of the run */
d->hours=do_alloc(d->free, sizeof(*(d->hours)));
d->hours_d=do_alloc(d->free, sizeof(*(d->hours_d)));
for(i=0;i<d->free;i++){
	d->hours[i]=(1.0*(d->gps[i]-d->gps[0]))/3600.0;
	d->hours_d[i]=(1.0*(d->gps[i]-d->gps[0]))/3600.0;
	}
	
/* compute frequency array */
d->frequencies=do_alloc(d->nbins, sizeof(*(d->frequencies)));
d->freq_d=do_alloc(d->nbins, sizeof(*(d->freq_d)));
for(i=0;i<d->nbins;i++){
	d->freq_d[i]=(1.0*(d->first_bin+i))/d->coherence_time;
	d->frequencies[i]=(1.0*(d->first_bin+i))/d->coherence_time;
	}

compute_detector_speed(d);

d->TMedians=do_alloc(d->free, sizeof(*(d->TMedians)));
d->expTMedians=do_alloc(d->free, sizeof(*(d->expTMedians)));
d->FMedians=do_alloc(d->nbins, sizeof(*(d->TMedians)));

compute_noise_curves(d);

for(i=0;i<d->free;i++){
	d->expTMedians[i]=exp(-M_LN10*2.0*(d->TMedians[i]-d->TMedian));
	}
d->expTMedian=exp(-M_LN10*2.0*d->TMedian);

fprintf(LOG, "dataset %s TMedian : %f\n", d->name, d->TMedian);

if(args_info.subtract_background_arg){
	fprintf(LOG, "dataset %s subtract background: yes\n", d->name);
	for(i=0;i<d->free;i++){
		for(j=0;j<d->nbins;j++){
			d->power[i*d->nbins+j]-=exp(M_LN10*(d->TMedians[i]+d->FMedians[j]));
			}
		}
	}

fprintf(stderr, "Obtaining whole sky AM response for dataset %s\n", d->name);
get_whole_sky_AM_response(d->gps, d->free, args_info.orientation_arg, &(d->AM_coeffs_plus), &(d->AM_coeffs_cross), &(d->AM_coeffs_size));

fprintf(stderr, "Initializing polarizations for dataset %s\n", d->name);
init_polarizations1(d->polarizations, d->AM_coeffs_plus, d->AM_coeffs_cross, d->AM_coeffs_size);

/* Check AM_response for correctness */
fprintf(stderr, "Verifying AM response computation for dataset %s\n", d->name);
for(i=0;i<ntotal_polarizations;i++){
	verify_whole_sky_AM_response(d->gps, d->free, d->polarizations[i].orientation, 
		fine_grid, 
		d->polarizations[i].AM_coeffs, d->polarizations[i].name);	
	}

fprintf(stderr,"Computing cutoff values for dataset %s\n", d->name);
/* compute CutOff values for each patch */

tm=do_alloc(d->free, sizeof(*tm));

for(i=0;i<patch_grid->npoints;i++){
	
	for(m=0;m<ntotal_polarizations;m++){
		if(patch_grid->band[i]<0){
			d->polarizations[m].patch_CutOff[i]=0.0;
			continue;
			}
		for(j=0;j<d->free;j++)tm[j]=1.0/(sqrt(d->expTMedians[j])*AM_response(j, patch_grid, i, d->polarizations[m].AM_coeffs));
		d->polarizations[m].patch_CutOff[i]=FindCutOff(tm, d->free);
		}
	}

free(tm);

characterize_dataset(d);
d->validated=1;
return 1;
}

void output_dataset_info(DATASET *d)
{
RGBPic *p;
PLOT *plot;
int m;

fprintf(stderr, "dataset %s adjusted weight: %g\n", d->name, d->weight);
fprintf(LOG, "dataset %s adjusted weight: %g\n", d->name, d->weight);

if(fine_grid->max_n_dec<800){
	p=make_RGBPic(fine_grid->max_n_ra*(800/fine_grid->max_n_dec)+140, fine_grid->max_n_dec*(800/fine_grid->max_n_dec));
	} else 
	p=make_RGBPic(fine_grid->max_n_ra+140, fine_grid->max_n_dec);

plot=make_plot(p->width, p->height);

adjust_plot_limits_f(plot, d->hours, d->TMedians, d->free, 1, 1, 1);
draw_grid(p, plot, 0, 0);
draw_points_f(p, plot, COLOR(255,0,0), d->hours, d->TMedians, d->free, 1, 1);
snprintf(s, 2000, "%s_TMedians.png", d->name);
RGBPic_dump_png(s, p);
snprintf(s, 2000, "%s_TMedians.dat", d->name);
dump_floats(s, d->TMedians, d->free, 1);

adjust_plot_limits_f(plot, d->frequencies, d->FMedians, d->nbins, 1, 1, 1);
draw_grid(p, plot, 0, 0);
draw_points_f(p, plot, COLOR(255,0,0), d->frequencies, d->FMedians, d->nbins, 1, 1);
snprintf(s, 2000, "%s_FMedians.png", d->name);
RGBPic_dump_png(s, p);
snprintf(s, 2000, "%s_FMedians.dat", d->name);
dump_floats(s, d->FMedians, d->nbins, 1);

adjust_plot_limits_d(plot, d->freq_d, d->mean, d->nbins, 1, 1, 1);
draw_grid(p, plot, 0, 0);
draw_points_d(p, plot, COLOR(255,0,0), d->freq_d, d->mean, d->nbins, 1, 1);
snprintf(s, 2000, "%s_mean.png", d->name);
RGBPic_dump_png(s, p);
snprintf(s, 2000, "%s_mean.dat", d->name);
dump_doubles(s, d->mean, d->nbins, 1);

adjust_plot_limits_d(plot, d->freq_d, d->weighted_mean, d->nbins, 1, 1, 1);
draw_grid(p, plot, 0, 0);
draw_points_d(p, plot, COLOR(255,0,0), d->freq_d, d->weighted_mean, d->nbins, 1, 1);
snprintf(s, 2000, "%s_weighted_mean.png", d->name);
RGBPic_dump_png(s, p);
snprintf(s, 2000, "%s_weighted_mean.dat", d->name);
dump_doubles(s, d->weighted_mean, d->nbins, 1);

adjust_plot_limits_d(plot, d->freq_d, d->sigma, d->nbins, 1, 1, 1);
draw_grid(p, plot, 0, 0);
draw_points_d(p, plot, COLOR(255,0,0), d->freq_d, d->sigma, d->nbins, 1, 1);
snprintf(s, 2000, "%s_sigma.png", d->name);
RGBPic_dump_png(s, p);

adjust_plot_limits_d(plot, d->freq_d, d->new_mean, d->nbins, 1, 1, 1);
draw_grid(p, plot, 0, 0);
draw_points_d(p, plot, COLOR(255,0,0), d->freq_d, d->new_mean, d->nbins, 1, 1);
snprintf(s, 2000, "%s_new_mean.png", d->name);
RGBPic_dump_png(s, p);
snprintf(s, 2000, "%s_new_mean.dat", d->name);
dump_doubles(s, d->new_mean, d->nbins, 1);

adjust_plot_limits_d(plot, d->freq_d, d->new_weighted_mean, d->nbins, 1, 1, 1);
draw_grid(p, plot, 0, 0);
draw_points_d(p, plot, COLOR(255,0,0), d->freq_d, d->new_weighted_mean, d->nbins, 1, 1);
snprintf(s, 2000, "%s_new_weighted_mean.png", d->name);
RGBPic_dump_png(s, p);
snprintf(s, 2000, "%s_new_weighted_mean.dat", d->name);
dump_doubles(s, d->new_weighted_mean, d->nbins, 1);

adjust_plot_limits_d(plot, d->freq_d, d->new_sigma, d->nbins, 1, 1, 1);
draw_grid(p, plot, 0, 0);
draw_points_d(p, plot, COLOR(255,0,0), d->freq_d, d->new_sigma, d->nbins, 1, 1);
snprintf(s, 2000, "%s_new_sigma.png", d->name);
RGBPic_dump_png(s, p);

for(m=0;m<ntotal_polarizations;m++){
	plot_grid_f(p, patch_grid, d->polarizations[m].patch_CutOff, 1);
	snprintf(s,20000,"%s_patch_CutOff_%s.png", d->name, d->polarizations[m].name);
	RGBPic_dump_png(s, p);
	}

free_plot(plot);
free_RGBPic(p);
}

void output_datasets_info(void) 
{
int i;
fprintf(stderr, "Output dataset summary data\n");
for(i=0;i<d_free;i++)output_dataset_info(&(datasets[i]));
}

static int get_geo_range(DATASET *d, char *filename, long startbin, long count, float *data, INT64 *gps)
{
FILE *fin;
REAL8 a, timebase;
INT4 b, bin_start, nbins;
REAL4 *tmp;
float factor;
long i;
long retries;
errno=0;
retries=0;
while((fin=fopen(filename,"r"))==NULL) {
	//if((fin==NULL) && (errno==ENOENT))return -1; /* no such file */
	fprintf(stderr,"Error opening file \"%s\":", filename);
	perror("");
	retries++;
	sleep(1);
	}
if(retries>0) { 
	fprintf(stderr, "Successfully opened file \"%s\"\n", filename);
	}
/* read header */	
/* Key */
a=0.0;
fread(&a, sizeof(a), 1, fin);
if(a!=1.0){
	fprintf(stderr,"Cannot read file \"%s\": wrong endianness\n", filename);
	return -1;
	}
/* gps */
fread(&b, sizeof(b), 1, fin);
*gps=b;

if(!check_intervals(d->segment_list, *gps)){
	fclose(fin);
	return -1;
	}
if(check_intervals(d->veto_segment_list, *gps)>0){
	fclose(fin);
	return -1;
	}


/* skip nsec */
fread(&b, sizeof(b), 1, fin);
/* timebase */
fread(&timebase, sizeof(a), 1, fin);

fread(&bin_start, sizeof(bin_start), 1, fin);
fread(&nbins, sizeof(nbins), 1, fin);

tmp=do_alloc(count*2, sizeof(*tmp));

fseek(fin, (startbin-bin_start)*8,SEEK_CUR);
if(fread(tmp,4,count*2,fin)<count*2){
	fprintf(stderr,"Not enough data in file \"%s\" gps=%lld.\n",filename,*gps);
	free(tmp);
	fclose(fin);
	return -1;
	}
/* reverse normalization applied to geo format files */
if(timebase < 0) {
	factor=1.0; /* make_sft_op did not apply normalization .. */
	fprintf(stderr,"** Timebase is negative, assuming unnormalized data\n");
	fprintf(LOG,"** Timebase is negative, assuming unnormalized data\n");
	} else {
	factor=(0.5*1800.0*16384.0)/nbins; /* use fixed normalization for 1800 sec SFTs .. */
	}
for(i=0;i<2*count;i++){
	data[i]=tmp[i]*factor;
	if(!isfinite(data[i])) {
		free(tmp);
		fprintf(stderr, "Infinite value encountered in file \"%s\"\n", filename);
		return -2;
		}
	}
free(tmp);
fclose(fin);
return 0;
}

static void expand_sft_array(DATASET *d, int count)
{
void *p;

d->size=2*d->size+count;
fprintf(stderr, "Growing %s SFT array to %f MB count=%d\n", d->name, d->size*d->nbins*sizeof(*d->bin)/(1024.0*1024.0), d->free);

p=do_alloc(d->size*d->nbins, sizeof(*(d->bin)));
if(d->free>0)memcpy(p, d->bin, d->free*d->nbins*sizeof(*(d->bin)));
free(d->bin);
d->bin=p;

p=do_alloc(d->size, sizeof(*(d->gps)));
if(d->free>0)memcpy(p, d->gps, d->free*sizeof(*(d->gps)));
free(d->gps);
d->gps=p;
}

static void add_file(DATASET *d, char *filename)
{
if(d->free>=d->size) {
	expand_sft_array(d, 1);
	}
d->gps[d->free]=0;
if(!get_geo_range(d, filename, d->first_bin, d->nbins, &(d->bin[d->free*d->nbins]), &(d->gps[d->free]))) {
	/* fprintf(stderr, "Loaded file %s (%lld)\n", filename, d->gps[d->free]); */
	d->free++;
	} else {
	fprintf(stderr, "Skipped file %s (%lld)\n", filename, d->gps[d->free]);
	}
}

static void d_read_directory(DATASET *dst, char *dir, int length)
{
DIR * d;
struct dirent *de;
char *s;
struct timeval start_time, end_time, last_time;
int i, last_i, limit=100;
double delta, delta_avg;
long retries;
s=do_alloc(length+20001, sizeof(*s));
memcpy(s, dir, length);
s[length]=0;
fprintf(stderr, "Reading directory %s\n", s);

retries=0;
errno=0;
while((d=opendir(s))==NULL) {
	int errsv=errno;
	fprintf(stderr, "Error reading directory %s: %s\n", s, strerror(errsv));
	retries++;
	sleep(1);
	}
if(retries>0) {
	fprintf(stderr, "Successfully opened directory %s\n", s);
	}
if(args_info.enable_dataset_locking_arg && (dst->lock_file != NULL)) {
	acquire_lock(dst->lock_file);
	} else 
if(args_info.lock_file_given) {
	acquire_lock(args_info.lock_file_arg);
	}
gettimeofday(&start_time, NULL);
last_time=start_time;
i=0;
last_i=0;
while((de=readdir(d))!=NULL) {
	if(de->d_name[0]!='.') {
		/* fprintf(stderr, "Found file \"%s\"\n", de->d_name); */
		snprintf(s+length, 20000, "/%s", de->d_name);
		add_file(dst, s);
		i++;
		}
	if(i > limit+last_i) {
		gettimeofday(&end_time, NULL);

		delta=end_time.tv_sec-last_time.tv_sec+(end_time.tv_usec-last_time.tv_usec)*1e-6;
		if(fabs(delta)==0.0)delta=1.0;

		delta_avg=end_time.tv_sec-start_time.tv_sec+(end_time.tv_usec-start_time.tv_usec)*1e-6;
		if(fabs(delta_avg)==0.0)delta_avg=1.0;

		fprintf(stderr, "Read rate %f SFTs/sec current, %f SFTs/sec avg\n", (i-last_i)/(delta), i/delta_avg);
		last_time=end_time;
		limit=ceil(10*limit/delta);
		if(limit <1 )limit=100;
		last_i=i;
		}
	}
closedir(d);
gettimeofday(&end_time, NULL);
delta=end_time.tv_sec-start_time.tv_sec+(end_time.tv_usec-start_time.tv_usec)*1e-6;
if(fabs(delta)==0.0)delta=1.0;
fprintf(stderr, "Read rate %f SFTs/sec\n", i/(delta));
free(s);
release_lock();
}

static void locate_arg(char *line, int length, int arg, int *arg_start, int *arg_stop)
{
int i;
int k;
*arg_start=0;
*arg_stop=0;
for(i=0,k=0;i<length;i++) {
	if((line[i]=='\n') || (line[i]=='\r'))break;
	if((line[i]==' ') || (line[i]=='\t')) {
		k++;
		while((i<length) && ((line[i]==' ') || (line[i]=='\t')))i++;
		}
	if(k==arg)break;
	}
if(k!=arg)return;
*arg_start=i;
while((i<length) && !((line[i]==' ') || (line[i]=='\t') || (line[i]=='\n')))i++;
*arg_stop=i;
if(line[*arg_start]=='"')(*arg_start)++;
if(((*arg_stop)>(*arg_start)) && (line[(*arg_stop)-1]=='"'))(*arg_stop)--;
}

static void gaussian_fill(DATASET *d, INT64 gps_start, int step, int count, double amp)
{
int i,k;
gsl_rng *rng=NULL;

fprintf(stderr, "Generating %d gaussian SFTs starting with gps %lld step %d amplitude %g for dataset %s\n",
	count, gps_start, step, amp, d->name);

fprintf(LOG, "Generating %d gaussian SFTs starting with gps %lld step %d amplitude %g for dataset %s\n",
	count, gps_start, step, amp, d->name);

if((d->free+count)>=d->size)expand_sft_array(d, count);

rng=gsl_rng_alloc(gsl_rng_default);

for(i=d->free;i<(d->free+count);i++) {
	d->gps[i]=gps_start+(i-d->free)*step;
	
	for(k=0;k<d->nbins;k++) {
		d->bin[i*d->nbins+k].re=amp*gsl_ran_gaussian(rng, 1.0);
		d->bin[i*d->nbins+k].im=amp*gsl_ran_gaussian(rng, 1.0);
		}
	}
d->free+=count;

gsl_rng_free(rng);
}

static void process_dataset_definition_line(char *line, int length)
{
int ai,aj;
DATASET *p;
/* skip whitespace in the beginning */
while(((*line)==' ') || ((*line)=='\t'))line++;
/* skip comments */
if((*line)=='#')return;
/* skip empty lines */
if((*line)=='\n')return;
if((*line)=='\r')return;
if((*line)==0)return;
if(!strncasecmp(line, "new_dataset", 11)) {
	if(d_free>=d_size) {
		d_size=2*d_size+5;
		p=do_alloc(d_size, sizeof(*datasets));
		if(d_free>0)memcpy(p, datasets, d_free*sizeof(*datasets));
		if(datasets!=NULL)free(datasets);
		datasets=p;
		}
	if( (d_free==0) || validate_dataset(&(datasets[d_free-1])))d_free++;
	memset(&(datasets[d_free-1]), 0, sizeof(*datasets));
	init_dataset(&(datasets[d_free-1]));
	locate_arg(line, length, 1, &ai, &aj);
	datasets[d_free-1].name=strndup(&(line[ai]), aj-ai);
	/* Automatically set correct offset */
	if(d_free==1)datasets[d_free-1].offset=0;
		else datasets[d_free-1].offset=datasets[d_free-2].offset+datasets[d_free-2].free;
	return;
	}
if(d_free<1) {
	fprintf(stderr, "One must use new_dataset to start a new dataset first\n");
	exit(-1);
	}
if(!strncasecmp(line, "detector", 8)) {
	locate_arg(line, length, 1, &ai, &aj);
	datasets[d_free-1].detector=strndup(&(line[ai]), aj-ai);
	} else 
if(!strncasecmp(line, "lock_file", 9)) {
	locate_arg(line, length, 1, &ai, &aj);
	if(datasets[d_free-1].lock_file!=NULL)free(datasets[d_free-1].lock_file);
	datasets[d_free-1].lock_file=strndup(&(line[ai]), aj-ai);
	} else 
if(!strncasecmp(line, "gps_start", 9)) {
	locate_arg(line, length, 1, &ai, &aj);
	datasets[d_free-1].gps_start=atoll(&(line[ai]));
	} else 
if(!strncasecmp(line, "gps_stop", 8)) {
	locate_arg(line, length, 1, &ai, &aj);
	datasets[d_free-1].gps_stop=atoll(&(line[ai]));
	} else 
if(!strncasecmp(line, "directory", 9)) {
	locate_arg(line, length, 1, &ai, &aj);
	d_read_directory(&(datasets[d_free-1]), &(line[ai]), aj-ai);
	} else 
if(!strncasecmp(line, "weight", 8)) {
	locate_arg(line, length, 1, &ai, &aj);
	datasets[d_free-1].weight=atof(&(line[ai]));
	} else 
if(!strncasecmp(line, "segments_file", 13)) {
	char *s2;
	locate_arg(line, length, 1, &ai, &aj);
	if(datasets[d_free-1].segment_list!=NULL)free_interval_set(datasets[d_free-1].segment_list);
	datasets[d_free-1].segment_list=new_interval_set();
	s2=strndup(&(line[ai]), aj-ai);
	add_intervals_from_file(datasets[d_free-1].segment_list, s2);
	free(s2);
	} else
if(!strncasecmp(line, "veto_segments_file", 18)) {
	char *s2;
	locate_arg(line, length, 1, &ai, &aj);
	if(datasets[d_free-1].veto_segment_list!=NULL)free_interval_set(datasets[d_free-1].veto_segment_list);
	datasets[d_free-1].veto_segment_list=new_interval_set();
	s2=strndup(&(line[ai]), aj-ai);
	add_intervals_from_file(datasets[d_free-1].veto_segment_list, s2);
	free(s2);
	} else
if(!strncasecmp(line, "apply_hanning_filter", 20)) {
	datasets[d_free-1].apply_hanning_filter=1;	
	} else
if(!strncasecmp(line, "gaussian_fill", 13)) {
	INT64 gps_start;
	int step, count;
	double amp;

	locate_arg(line, length, 1, &ai, &aj);
	sscanf(&(line[ai]), "%lld", &gps_start);

	locate_arg(line, length, 2, &ai, &aj);
	sscanf(&(line[ai]), "%d", &step);

	locate_arg(line, length, 3, &ai, &aj);
	sscanf(&(line[ai]), "%d", &count);

	locate_arg(line, length, 4, &ai, &aj);
	sscanf(&(line[ai]), "%lg", &amp);

	gaussian_fill(&(datasets[d_free-1]), gps_start, step, count, amp*datasets[d_free-1].coherence_time*16384.0);
	} else {
	fprintf(stderr, "*** Could not parse line: \n%*s\n *** Exiting.\n", length, line);
	exit(-1);
	}	
}

void load_datasets(char *definition)
{
char *p=definition;
int k;

if(d_size<1) {
	d_size=10;
	d_free=-1;
	datasets=do_alloc(d_size, sizeof(*datasets));
	}

while(1){
	/* skip whitespace in the beginning */
	while((*p==' ') | (*p=='\t'))p++;

	for(k=0;p[k] && p[k]!='\n' && p[k]!='\r';k++);

	/* skip comments */
	if(*p=='#'){
		p+=k;
		if(p[k])p++;
		continue;
		}
	/* skip empty lines */
	if(k==0){
		if(!p[k])break;
		p++;
		continue;
		}
	process_dataset_definition_line(p, k);
	}
}

void load_dataset_from_file(char *file)
{
FILE *fin;
char *s;
int size=60000;
int k;
fin=fopen(file, "r");
if(fin==NULL){
	fprintf(stderr, "*** could not read dataset file \"%s\":", file);
	perror("");
	exit(-1);
	}
s=do_alloc(size, sizeof(*s));
while(1) {
	if(fgets(s, size-1, fin)==NULL) {
		if(feof(fin))break;
		fprintf(stderr, "Error reading file \"%s\": ", file);
		perror("");
		exit(-1);
		}
	s[size-1]=0;
	k=strlen(s);
	/* fprintf(stderr, "%d %s\n", k, s); */
	process_dataset_definition_line(s, k);
	}
if(d_free>0)
	if(!validate_dataset(&(datasets[d_free-1])))d_free--;
fclose(fin);
}

long total_segments(void)
{
int i;
long total=0;
for(i=0;i<d_free;i++){
	total+=datasets[i].free;
	}
return(total);
}

void compute_powers(DATASET *dataset)
{
COMPLEX *tmp;
int i,k;
float *p;
float x,y;
COMPLEX *c;
int nbins=dataset->nbins;
int nsegments=dataset->free;
tmp=do_alloc(dataset->nbins, sizeof(*tmp));

for(i=0;i<nsegments;i++) {
	memcpy(tmp, &(dataset->bin[i*nbins]), nbins*sizeof(*tmp));
	/* inject_fake_signal(dataset, i, tmp); */
	c=tmp;
	p=&(dataset->power[i*nbins]);
	for(k=0;i<nbins;k++) {
		x=(*c).re;
		y=(*c).im;
		*p=x*x+y*y;
		p++;
		c++;
		}
	}
}

static float compute_median(float *firstbin, int step, int count)
{
float *tmp;
int i;
tmp=alloca(count*sizeof(float));
for(i=0;i<count;i++)tmp[i]=firstbin[i*step];
sort_floats(tmp, count);
if(!(count & 1))return (tmp[count>>1]+tmp[(count>>1)-1])/2.0;
return tmp[count>>1];
}

void compute_noise_curves(DATASET *dataset)
{
float *tmp;
float *p,*t;
float a;
int i,j;
float b, b_initial;
int nsegments=dataset->free;
int nbins=dataset->nbins;
HISTOGRAM *hist;

tmp=do_alloc(nsegments*nbins,sizeof(float));
for(i=0;i<nsegments;i++){
	t=&(tmp[i*nbins]);
	p=&(dataset->power[i*nbins]);
	for(j=0;j<nbins;j++){
		t[j]=log10(p[j]);
		}
	}
b=0;
for(i=0;i<nsegments;i++){
	a=compute_median(tmp+i*nbins,1,nbins);
	dataset->TMedians[i]=a;
	b+=a*a;
	t=&(tmp[i*nbins]);
	for(j=0;j<nbins;j++){
		t[j]-=a;
		}
	}
for(j=0;j<nbins;j++){
	a=compute_median(tmp+j,nbins,nsegments);
	dataset->FMedians[j]=a;
	b+=a*a;
	t=&(tmp[j]);
	for(i=0;i<nsegments;i++){
		t[i*nbins]-=a;
		}
	}
b_initial=b;
while(b>0){
	b=0;
	for(i=0;i<nsegments;i++){
		a=compute_median(tmp+i*nbins,1,nbins);
		dataset->TMedians[i]+=a;
		b+=a*a;
		t=&(tmp[i*nbins]);
		for(j=0;j<nbins;j++){
			t[j]-=a;
			}
		}
	for(j=0;j<nbins;j++){
		a=compute_median(tmp+j,nbins,nsegments);
		dataset->FMedians[j]+=a;
		b+=a*a;
		t=&(tmp[j]);
		for(i=0;i<nsegments;i++){
			t[i*nbins]-=a;
			}
		}
	fprintf(stderr,"%g\n",b);
	/* if one of r nbins or nsegments are even break because of order of
	   magnitude, do not try to solve exactly */
	if(!((nsegments &1)&&(nbins&1)) && (b<(b_initial*(1E-16))))break;
	}

hist=new_histogram(args_info.hist_bins_arg, 1);
compute_histogram_f(hist, tmp, NULL, nsegments*nbins);
snprintf(s, 20000, "hist_residuals: %s", dataset->name);
print_histogram(LOG, hist, s);
print_histogram(stderr, hist, s);
free_histogram(hist);

free(tmp);
dataset->TMedian=compute_median(dataset->TMedians,1,nsegments);
fprintf(LOG,"Median noise level (TMedian): %s %g\n",dataset->name, dataset->TMedian);
fflush(LOG);
}

void characterize_dataset(DATASET *d)
{
double a,w, CutOff;
int i, j, count;
float *tmp;
/*
for(i=0;i<d->free;i++){
	d->expTMedians[i]=exp(-M_LN10*2.0*(d->TMedians[i]-d->TMedian));
	}
d->expTMedian=exp(-M_LN10*2.0*d->TMedian);
*/


d->mean=do_alloc(d->nbins, sizeof(*(d->mean)));
d->weighted_mean=do_alloc(d->nbins, sizeof(*(d->weighted_mean)));
d->sigma=do_alloc(d->nbins, sizeof(*(d->sigma)));
d->bin_weight=do_alloc(d->nbins, sizeof(*(d->bin_weight)));

d->new_mean=do_alloc(d->nbins, sizeof(*(d->new_mean)));
d->new_weighted_mean=do_alloc(d->nbins, sizeof(*(d->new_weighted_mean)));
d->new_sigma=do_alloc(d->nbins, sizeof(*(d->new_sigma)));
d->new_bin_weight=do_alloc(d->nbins, sizeof(*(d->new_bin_weight)));

/* compute CutOff for non-demodulated analysis */
tmp=do_alloc(d->free, sizeof(*tmp));
for(i=0;i<d->free;i++)tmp[i]=exp(M_LN10*d->TMedians[i]);
CutOff=log10(FindCutOff(tmp, d->free));
fprintf(stderr,"%s CutOff=%g\n", d->name, CutOff);
free(tmp);
tmp=NULL;

for(i=0;i<d->nbins;i++){
	d->weighted_mean[i]=0.0;
	d->bin_weight[i]=0.0;
	}
for(i=0;i<d->free;i++){
	w=d->expTMedians[i];
	for(j=0;j<d->nbins;j++){
		a=d->power[i*d->nbins+j];
		d->bin_weight[j]+=w;
		d->weighted_mean[j]+=a*w;
		}
	}
for(i=0;i<d->nbins;i++){
	d->weighted_mean[i]/=d->bin_weight[i];
	}

#if 0
adjust_plot_limits_d(plot, freq_d, weighted_mean, d->nbins, 1, 1, 1);
draw_grid(p, plot, 0, 0);
draw_points_d(p, plot, COLOR(255,0,0), freq_d, weighted_mean, d->nbins, 1, 1);
RGBPic_dump_png("weighted_mean.png", p);
dump_doubles("weighted_mean.dat", weighted_mean, d->nbins, 1);
#endif

/* second pass - apply CutOff */

for(i=0;i<d->nbins;i++){
	d->new_mean[i]=0;
	d->new_sigma[i]=0;
	d->new_weighted_mean[i]=0;
	d->new_bin_weight[i]=0;
	}
count=0;
for(i=0;i<d->free;i++) {
	if(d->TMedians[i]>=CutOff)continue;
	w=d->expTMedians[i];
	count++;
	for(j=0;j<d->nbins;j++){
		a=d->power[i*d->nbins+j];
		d->new_mean[j]+=a;
		d->new_sigma[j]+=a*a;
		d->new_bin_weight[j]+=w;
		d->new_weighted_mean[j]+=a*w;
		}
	}
for(i=0;i<d->nbins;i++){
	d->new_mean[i]/=count;
	d->new_sigma[i]=sqrt((d->new_sigma[i]-count*d->new_mean[i]*d->new_mean[i])/(count-1));
	d->new_weighted_mean[i]/=d->new_bin_weight[i];
	}

fprintf(stderr,"Checking background lines\n");

d->lines_report=make_lines_report(0, d->nbins-1, 5);
sprintf(s, "%s background", d->name);
detect_lines_d(d->new_weighted_mean, d->lines_report, s);
print_lines_report(LOG, d->lines_report, s);
print_lines_report(stderr, d->lines_report, s);

if(!args_info.filter_lines_arg) {
	fprintf(stderr,"Background lines removal will not be performed\n");
	fprintf(LOG, "Background lines removal will not be performed\n");
	for(i=0;i<d->lines_report->nlines;i++)d->lines_report->lines_list[i]=-1; /* no line */
	}
fflush(LOG);
}

float datasets_normalizing_weight(void)
{
int i;
float max_tm=0;
for(i=0;i<d_free;i++){
	if(!i || datasets[i].TMedian<max_tm)max_tm=datasets[i].TMedian;
	}
return(exp(-M_LN10*max_tm));
}

INT64 min_gps(void)
{
int i,j;
INT64 min_gps=0;

for(i=0;i<d_free;i++){
	for(j=0;j<datasets[i].free;j++)
		if(!min_gps || datasets[i].gps[j]<min_gps)min_gps=datasets[i].gps[j];
	}
return(min_gps);
}

INT64 max_gps(void)
{
int i,j;
INT64 max_gps=0;

for(i=0;i<d_free;i++){
	for(j=0;j<datasets[i].free;j++)
		if(datasets[i].gps[j]>max_gps)max_gps=datasets[i].gps[j];
	}
return(max_gps);
}

void post_init_datasets(void)
{
float a, b;
int i;
DATASET *d;

/* adjust d->weight in proportion to exp(TMedian*C)  */
a=datasets_normalizing_weight();
for(i=0;i<d_free;i++){
	d=&(datasets[i]);
	b=a*exp(M_LN10*d->TMedian); 
	d->weight*=b*b;
	}
}

void datasets_average_detector_speed(double *average_det_velocity)
{
int i, j;
double total_weight=0.0, w;
DATASET *d;
for(i=0;i<3;i++)average_det_velocity[i]=0.0;

for(i=0;i<d_free;i++){
	d=&(datasets[i]);
	for(j=0;j<d->free;j++) {
		w=d->expTMedians[j]*d->weight;
		total_weight+=w;
		if(isnan(total_weight))exit(-1);
		average_det_velocity[0]+=d->detector_velocity[3*j]*w;
		average_det_velocity[1]+=d->detector_velocity[3*j+1]*w;
		average_det_velocity[2]+=d->detector_velocity[3*j+2]*w;
		}
	}
for(i=0;i<3;i++)average_det_velocity[i]/=total_weight;
}

void inject_signal(void) 
{
int i,j;
DATASET *d;
for(i=0;i<d_free;i++){
	d=&(datasets[i]);
	for(j=0;j<d->free;j++) {
		inject_fake_signal(d, j);
		}
	}
}


void dump_datasets(char *filename) 
{
FILE *fout;
int i,j,k;
DATASET *d;
double x,y;

fout=fopen(filename, "w");
if(fout==NULL) {
	fprintf(stderr, "Error opening files %s:", filename);
	perror("");
	return;
	}
fprintf(fout, "Dataset\tDetector\tGPS");
for(k=0;k<nbins;k++) fprintf(fout, "\tR%d", k+first_bin);
for(k=0;k<nbins;k++) fprintf(fout, "\tI%d", k+first_bin);
for(k=0;k<nbins;k++) fprintf(fout, "\tP%d", k+first_bin);
fprintf(fout, "\n");

for(i=0;i<d_free;i++) {
	d=&(datasets[i]);
	for(j=0;j<d->free;j++) {
		fprintf(fout, "%s\t%s\t%lld", d->name, d->detector, d->gps[j]);
		for(k=0;k<d->nbins;k++) {
			fprintf(fout, "\t%8g", d->bin[j*d->nbins+k].re);
			}
		for(k=0;k<d->nbins;k++) {
			fprintf(fout, "\t%8g",d->bin[j*d->nbins+k].im);
			}
		for(k=0;k<d->nbins;k++) {
			x=d->bin[j*d->nbins+k].re;
			y=d->bin[j*d->nbins+k].im;
			fprintf(fout, "\t%8g", x*x+y*y);
			}
		fprintf(fout, "\n");
		}
	}
fclose(fout);
}
