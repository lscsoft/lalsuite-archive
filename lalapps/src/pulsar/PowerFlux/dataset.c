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

#define _GNU_SOURCE
#define _FILE_OFFSET_BITS 64
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <alloca.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <regex.h>

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
#include "hookup.h"
#include "util.h"
#include "jobs.h"

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
extern double spindown;

char s[20000];

DATASET * datasets=NULL;
int d_free=0;
int d_size=0;

extern int do_CutOff;

int lock_file=-1;

extern EphemerisData ephemeris;

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

/* we need the structure to be global as a workaround for a bug in condor libraries:
   this way the pointer is 32bit even on 64bit machines */

static int acquire_lock(char *filename, int wait)
{
int i;
struct flock fl;

lock_file=open(filename, O_CREAT | O_RDWR, 0666);

memset(&fl, 0, sizeof(fl));
fl.l_type=F_WRLCK;
fl.l_whence=SEEK_SET;
fl.l_start=0;
fl.l_len=1;
errno=0;
i=0;
while((lock_file<0) || (direct_fcntl(lock_file, F_SETLK, &fl)<0)){
	if(lock_file<0) {
		fprintf(stderr, "Could not open file \"%s\" for locking\n", filename);
		}
	if(i>10) {
		fprintf(stderr, "Waiting for lock: %s\n", strerror(errno));
		i=0;
		}
	if(lock_file>=0)close(lock_file);
	if(!wait)return(0);
	condor_safe_sleep(args_info.lock_retry_delay_arg);
	lock_file=open(filename, O_CREAT | O_RDWR, 0666);
	i++;
	}
return(1);
}

static void release_lock(void)
{
if(lock_file>=0)close(lock_file);
lock_file=-1;
}

typedef struct {
	double a_plus;
	double a_cross;
	double cos_e;
	double sin_e;

	double segment_start;

	double ref_time;

	double freq;
	double spindown;
	double dInv; /* Inverse distance to source in seconds */

	double ra;
	double dec;
	double iota;
	double psi;
	double phi;
	
	/* Parameters for additional sinusoidal phase modulation */
	double phase_modulation_depth;
	double phase_modulation_freq;
	double phase_modulation_phase;

	/* Parameters for additional sinusoidal frequency modulation */
	double freq_modulation_depth;
	double freq_modulation_freq;
	double freq_modulation_phase;

	double strain;

	double coherence_time;

	double extra_phase; /* phase accumulated from the start of the run */
	double phase_integration_factor;

	float e[26];
	int bin;
	
	char *detector;
	} SIGNAL_PARAMS;

	
static void compute_signal(double *re, double *im, double *f, double t, SIGNAL_PARAMS *p)
{
double doppler, omega_t, c_omega_t, s_omega_t, modomega_t, fmodomega_t, a;
double hann;
float f_plus, f_cross;
EmissionTime emission_time;
LIGOTimeGPS tGPS;
EarthState earth_state;
double te, phase_spindown, phase_freq;
LALStatus status={level:0, statusPtr:NULL};

tGPS.gpsSeconds=floor(t);
tGPS.gpsNanoSeconds=(t-floor(t))*1e9;

LALBarycenterEarth(&status, &(earth_state), &tGPS, &ephemeris);
TESTSTATUS(&status);
if(status.statusPtr)FREESTATUSPTR(&status);

get_emission_time(&emission_time, &(earth_state), p->ra, p->dec, p->dInv, p->detector, tGPS);

get_AM_response_d(t, p->dec, p->ra, 0.0, p->detector, &f_plus, &f_cross);
// fprintf(stderr, "%d f_plus=%f f_cross=%f\n", (int)round(t), f_plus, f_cross);

doppler=p->e[0]*emission_time.vDetector[0]+p->e[1]*emission_time.vDetector[1]+p->e[2]*emission_time.vDetector[2];

/* Compute SSB time since ref time */
te=(emission_time.te.gpsSeconds-p->ref_time)+((double)(1e-9))*emission_time.te.gpsNanoSeconds;

fmodomega_t=2.0*M_PI*(te*p->freq_modulation_freq-floor(te*p->freq_modulation_freq));

*f=(p->freq+p->spindown*te+p->freq_modulation_depth*cos(fmodomega_t+p->freq_modulation_phase))*(1.0+doppler);

phase_spindown=0.5*te*te*p->spindown;


phase_freq=(p->freq-(double)(p->bin)/(double)(p->coherence_time))*te
	+(double)p->bin*(te-(t-p->segment_start))/(double)(p->coherence_time);
	
if(fabs(p->freq_modulation_freq)>1e-28) {
	/* Split up computation according to the formula below.
	   Note that for small freq_modulation_depth the second component will be zero due to precision loss and we need random phi to properly simulate global phase.

	   sin(wt+a)    sin(wt)cos(a)     cos(wt)sin(a) 
	   --------- =  ------------- +   ------------- 
	       w             w                 w        
	*/
	phase_freq+=p->freq_modulation_depth*sin(fmodomega_t)*cos(p->freq_modulation_phase)/(2*M_PI*p->freq_modulation_freq);
	a=p->freq_modulation_depth*cos(fmodomega_t)*sin(p->freq_modulation_phase)/(2*M_PI*p->freq_modulation_freq);
	phase_freq+=a-floor(a);
	} else {
	/* we just have a constant frequency offset for practical purposes which corresponds to freezing 
	   omega at the level of definition of *f, not phase_freq which would diverge for anything but sin(), but that is an extra random constant offset in phase which should be taken care of by using random phi :
	   
	   sin(wt+a)    sin(wt)cos(a)     cos(wt)sin(a)               sin(a)     wt^2 sin(a) 
	   --------- =  ------------- +   ------------- = t cos(a) +  ------ -   ----------- + o(w^2)
	       w             w                 w                        w             2
	*/
	phase_freq+=p->freq_modulation_depth*cos(p->freq_modulation_phase)*te;
	}

omega_t=2.0*M_PI*((phase_freq-floor(phase_freq))+(phase_spindown-floor(phase_spindown)))+p->phi;

/* add contribution from sinusoidal phase modulation */
modomega_t=2.0*M_PI*(te*p->phase_modulation_freq-floor(te*p->phase_modulation_freq))+p->phase_modulation_phase;
omega_t+=p->phase_modulation_depth*sin(modomega_t);

hann=0.5*(1.0-cos(2.0*M_PI*(t-p->segment_start)/p->coherence_time));

sincos(omega_t, &s_omega_t, &c_omega_t);

/* 
	cos(a)*cos(-b)=0.5*(cos(a-b)+cos(a+b))
	sin(a)*cos(-b)=0.5*(sin(a-b)+sin(a+b))
	cos(a)*sin(-b)=0.5*(sin(a-b)-sin(a+b))
	sin(a)*sin(-b)=0.5*(cos(a-b)-cos(a+b))

	|f-w|<<1
	cos(ft)*exp(-iwt)~= 0.5*(cos((f-w)t)+i*sin((f-w)t))
	sin(ft)*exp(-iwt)~= 0.5*(sin((f-w)t)-i*cos((f-w)t))
*/

// fprintf(stderr, "a_plus=%f a_cross=%f cos_e=%f sin_e=%f\n", p->a_plus, p->a_cross, p->cos_e, p->sin_e);

/* Note: we do not have an extra factor of 0.5 because of normalization convention used, 
  there is an extra factor of 2.0 to compensate for Hann windowing */

*re=hann*(f_plus*(p->a_plus*c_omega_t*p->cos_e-p->a_cross*s_omega_t*p->sin_e)+
	f_cross*(p->a_plus*c_omega_t*p->sin_e+p->a_cross*s_omega_t*p->cos_e));

*im=hann*(f_plus*(p->a_plus*s_omega_t*p->cos_e+p->a_cross*c_omega_t*p->sin_e)+
	f_cross*(p->a_plus*s_omega_t*p->sin_e-p->a_cross*c_omega_t*p->cos_e));

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

static double signal_omega(double t, void *params)
{
SIGNAL_PARAMS *p=params;
float det_vel[3];
double doppler, omega;

get_detector_vel(round(t), det_vel);

doppler=p->e[0]*det_vel[0]+p->e[1]*det_vel[1]+p->e[2]*det_vel[2];
omega=2.0*M_PI*p->freq*doppler-p->phase_integration_factor;

return(omega);
}


static void precompute_signal_params(SIGNAL_PARAMS *p)
{
double cos_i;

cos_i=cos(p->iota);

p->a_plus=(1+cos_i*cos_i)*0.5;
p->a_cross=cos_i;

p->cos_e=cos(2.0*p->psi);
p->sin_e=sin(2.0*p->psi);
	
precompute_am_constants(p->e, p->ra, p->dec);
}

static void fill_signal_params_with_defaults(SIGNAL_PARAMS *p)
{
memset(p, 0, sizeof(SIGNAL_PARAMS));
p->freq_modulation_freq=1.0;

precompute_signal_params(p);
}

static void fill_signal_params_from_args_info(SIGNAL_PARAMS *p)
{
p->bin=0;

p->ra=args_info.fake_ra_arg;
p->dec=args_info.fake_dec_arg;
p->freq=args_info.fake_freq_arg;
p->spindown=args_info.fake_spindown_arg;
p->ref_time=args_info.fake_ref_time_arg;
p->strain=args_info.fake_strain_arg;
p->iota=args_info.fake_iota_arg;
p->psi=args_info.fake_psi_arg;
p->phi=args_info.fake_phi_arg;
p->dInv=args_info.fake_dInv_arg;
p->freq_modulation_depth=args_info.fake_freq_modulation_depth_arg;
p->freq_modulation_freq=args_info.fake_freq_modulation_freq_arg;
p->freq_modulation_phase=args_info.fake_freq_modulation_phase_arg;
p->phase_modulation_depth=args_info.fake_phase_modulation_depth_arg;
p->phase_modulation_freq=args_info.fake_phase_modulation_freq_arg;
p->phase_modulation_phase=args_info.fake_phase_modulation_phase_arg;

precompute_signal_params(p);
}

static void test_inject_fake_signal05(void)
{
SIGNAL_PARAMS p;
double re, im, f, re2, im2, f2;
int status=0;

memset(&p, 0, sizeof(p));

p.bin=0;

p.ra=2.56;
p.dec=1.43;
p.freq=500.12;
p.spindown=1.23e-9;
p.ref_time=793154935;
p.segment_start=793161250-900.0;
p.coherence_time=1800.0;
p.extra_phase=0;
p.iota=1.23;
p.psi=0.45;
p.phi=0.123;
p.dInv=0;
p.detector="LHO";

precompute_signal_params(&p);

get_detector("LHO");

compute_signal(&re, &im, &f, 793161250.0, &p);
compute_signal(&re2, &im2, &f2, 793161250.0+43*3600.0+123, &p);
fprintf(stderr, "compute_signal_test05a: %g %g %f %g %g %f %g\n", re, im, f, re2, im2, f2, re*re2+im*im2);
fprintf(LOG, "compute_signal_test05a: %g %g %f %g %g %f %g\n", re, im, f, re2, im2, f2, re*re2+im*im2);

if(fabs(re*re2+im*im2+0.0140269)>1e-5 ||
   fabs(f-500.101774)>1e-5 ||
   fabs(f2-500.101466)>1e-5
	) status|=1;

fprintf(stderr, "compute_signal_test05b: %g %g %f\n", re, im, f);
fprintf(LOG, "compute_signal_test05b: %g %g %f\n", re, im, f);

p.cos_e=1.0;
p.sin_e=0.0;

compute_signal(&re, &im, &f, 793161250.0, &p);
fprintf(stderr, "compute_signal_test05b: %g %g %f\n", re, im, f);
fprintf(LOG, "compute_signal_test05b: %g %g %f\n", re, im, f);

p.cos_e=1.0;
p.sin_e=0.0;
p.phi+=M_PI*0.5;

compute_signal(&re2, &im2, &f2, 793161250.0, &p);
fprintf(stderr, "compute_signal_test05c: %g %g %f\n", re2, im2, f2);
fprintf(LOG, "compute_signal_test05c: %g %g %f\n", re2, im2, f2);

if(fabs(f2-f)>1e-5 ||
   fabs(re-im2)>1e-11 ||
   fabs(im+re2)>1e-11) status|=2;

if(status) {
	fprintf(stderr, "compute_signal_test05: failed %d\n", status);
	fprintf(LOG, "compute_signal_test05: failed %d\n", status);
	exit(-1);
	}
}

static void test_inject_fake_signal09(void)
{
SIGNAL_PARAMS p;
double re, im, f, re2, im2, f2;
int status=0;

memset(&p, 0, sizeof(p));

p.bin=0;

p.ra=2.56;
p.dec=1.43;
p.freq=500.12;
p.spindown=1.23e-9;
p.ref_time=910000000.0-12345.0;
p.segment_start=910000000.0-900.0;
p.coherence_time=1800.0;
p.extra_phase=0;
p.iota=1.23;
p.psi=0.45;
p.phi=0.123;
p.dInv=0;
p.detector="LHO";

precompute_signal_params(&p);

get_detector("LHO");

compute_signal(&re, &im, &f, 930000000.0, &p);
compute_signal(&re2, &im2, &f2, 930000000.0+43*3600.0+123, &p);
fprintf(stderr, "compute_signal_test09a: %g %g %f %g %g %f %g\n", re, im, f, re2, im2, f2, re*re2+im*im2);
fprintf(LOG, "compute_signal_test09a: %g %g %f %g %g %f %g\n", re, im, f, re2, im2, f2, re*re2+im*im2);

if(fabs(re*re2+im*im2-0.0470977)>1e-5 ||
   fabs(f-500.140589)>1e-5 ||
   fabs(f2-500.141517)>1e-5) status|=1;

p.cos_e=1.0;
p.sin_e=0.0;

compute_signal(&re, &im, &f, 930000000.0, &p);
fprintf(stderr, "compute_signal_test09b: %g %g %f\n", re, im, f);
fprintf(LOG, "compute_signal_test09b: %g %g %f\n", re, im, f);

p.cos_e=1.0;
p.sin_e=0.0;
p.phi+=M_PI*0.5;

compute_signal(&re2, &im2, &f2, 930000000.0, &p);
fprintf(stderr, "compute_signal_test09c: %g %g %f\n", re2, im2, f2);
fprintf(LOG, "compute_signal_test09c: %g %g %f\n", re2, im2, f2);

if(fabs(f2-f)>1e-5 ||
   fabs(re-im2)>1e-11 ||
   fabs(im+re2)>1e-11) status|=2;

if(status) {
	fprintf(stderr, "compute_signal_test09: failed %d\n", status);
	fprintf(LOG, "compute_signal_test09: failed %d\n", status);
	exit(-1);
	}
}

static void inject_fake_signal(SIGNAL_PARAMS *p, DATASET *d, int segment)
{
double result, abserr;
int err;
int window=args_info.fake_injection_window_arg;
int bin;
int i;
double re, im, f;
gsl_function F; 
gsl_integration_workspace *w;
gsl_integration_qawo_table *t_sine, *t_cosine;
int w_size=1024*32;

F.params=p;
p->bin=0;

p->segment_start=d->gps[segment];
p->coherence_time=d->coherence_time;
p->detector=d->detector;


compute_signal(&re, &im, &f, d->gps[segment]+(int)(d->coherence_time/2), p);


bin=round(f*1800.0-d->first_bin);
if((bin+window>nbins) || (bin<window))fprintf(stderr, "Injected signal outside loaded band: bin=%d, segment=%d\n", bin, segment);
if(bin+window>nbins)bin=nbins-window;
if(bin<window)bin=window;

p->extra_phase=0.0;

w=gsl_integration_workspace_alloc(w_size);

for(i=bin-window; i<=bin+window; i++) {
	p->bin=i+d->first_bin;

	F.function=signal_re;
	err=gsl_integration_qag(&F, d->gps[segment], d->gps[segment]+d->coherence_time,
		1e-4, 1e-3, w_size, GSL_INTEG_GAUSS61, w,
		&result, &abserr);
	//fprintf(stderr, "re %d %d result=%g abserr=%g %s\n", segment, i, result, abserr, gsl_strerror(err)); 
	d->re[segment*d->nbins+i]+=p->strain*result*16384.0/args_info.strain_norm_factor_arg;


	F.function=signal_im;
	err=gsl_integration_qag(&F, d->gps[segment], d->gps[segment]+d->coherence_time,
		1e-4, 1e-3, w_size, GSL_INTEG_GAUSS61, w,
		&result, &abserr);
	//fprintf(stderr, "im %d %d result=%g abserr=%g %s\n", segment, i, result, abserr, gsl_strerror(err)); 
	d->im[segment*d->nbins+i]+=p->strain*result*16384.0/args_info.strain_norm_factor_arg;

	#if 0
	/* this code computes fast contribution from 2f it is not run because the values returned are many orders of magnitude smaller than injection strength (1e-15 for 1000 Hz) */
	
	t_sine=gsl_integration_qawo_table_alloc(-2*2*M_PI*p->bin*1.0/d->coherence_time, d->coherence_time, GSL_INTEG_SINE, 16);
	t_cosine=gsl_integration_qawo_table_alloc(-2*2*M_PI*p->bin*1.0/d->coherence_time, d->coherence_time, GSL_INTEG_COSINE, 16);
	
	F.function=signal_re;
	err=gsl_integration_qawo(&F, d->gps[segment], 
		1e-4, 1e-3, w_size, w, t_cosine,
		&result, &abserr);
	fprintf(stderr, "re2 %d %d result=%g abserr=%g %s\n", segment, i, result, abserr, gsl_strerror(err)); 
	d->re[segment*d->nbins+i]+=p->strain*result*16384.0/args_info.strain_norm_factor_arg;

	F.function=signal_im;
	err=gsl_integration_qawo(&F, d->gps[segment], 
		1e-4, 1e-3, w_size, w, t_sine,
		&result, &abserr);
	fprintf(stderr, "re3 %d %d result=%g abserr=%g %s\n", segment, i, result, abserr, gsl_strerror(err)); 
	d->re[segment*d->nbins+i]+=p->strain*result*16384.0/args_info.strain_norm_factor_arg;

	F.function=signal_re;
	err=gsl_integration_qawo(&F, d->gps[segment], 
		1e-4, 1e-3, w_size, w, t_sine,
		&result, &abserr);
	fprintf(stderr, "im2 %d %d result=%g abserr=%g %s\n", segment, i, result, abserr, gsl_strerror(err)); 
	d->im[segment*d->nbins+i]+=p->strain*result*16384.0/args_info.strain_norm_factor_arg;

	F.function=signal_im;
	err=gsl_integration_qawo(&F, d->gps[segment], 
		1e-4, 1e-3, w_size, w, t_cosine,
		&result, &abserr);
	fprintf(stderr, "im3 %d %d result=%g abserr=%g %s\n", segment, i, result, abserr, gsl_strerror(err)); 
	d->im[segment*d->nbins+i]-=p->strain*result*16384.0/args_info.strain_norm_factor_arg;

	gsl_integration_qawo_table_free(t_sine);
	gsl_integration_qawo_table_free(t_cosine);
	#endif
	}
	
gsl_integration_workspace_free(w);

}

static void init_dataset(DATASET *d)
{
int i;
memset(d, 0, sizeof(*d));
d->name="unknown";
for(i=0;i<MAX_LOCKS;i++)
	d->lock_file[i]=NULL;
d->validated=0;
d->detector="unknown";
d->gps_start=0;
d->gps_stop=-1;

d->coherence_time=1800.0;
d->nbins=nbins;
d->first_bin=first_bin;

d->segment_list=NULL;
d->veto_segment_list=NULL;
d->no_duplicate_gps=1;

d->size=1000;
d->free=0;
d->gps=do_alloc(d->size, sizeof(*d->gps));
d->dc_factor=1.0;
d->dc_factor_touched=0;
d->dc_factor_blocked=0;
d->re=do_alloc(d->size*d->nbins, sizeof(*d->re));
d->im=do_alloc(d->size*d->nbins, sizeof(*d->im));
d->sft_veto=NULL;
d->veto_level=1e-2; /* this takes care of SFTs that have 1/100 weight... */
d->veto_spike_level=1.7; /* this takes care of SFTs with spikes 50 times median level */

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
	get_detector_vel(d->gps[i]+(int)(d->coherence_time*0.5),&(d->detector_velocity[3*i]));
		
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
float *tmp, *p;
tmp=do_alloc(d->nbins, sizeof(*tmp));
for(i=0;i<d->free;i++) {
	p=&(d->re[i*d->nbins]);

	/* wrap around, just so that we don't have 0 exactly */
	tmp[0]=p[0]-(p[d->nbins-1]+p[1])*0.5;
	tmp[d->nbins-1]=p[d->nbins-1]-(p[d->nbins-2]+p[0])*0.5;

	for(k=1;k<(d->nbins-1);k++) {
		tmp[k]=p[k]-(p[k-1]+p[k+1])*0.5;
		}
	memcpy(p, tmp, nbins*sizeof(*p));	
	}
for(i=0;i<d->free;i++) {
	p=&(d->im[i*d->nbins]);

	/* wrap around, just so that we don't have 0 exactly */
	tmp[0]=p[0]-(p[d->nbins-1]+p[1])*0.5;
	tmp[d->nbins-1]=p[d->nbins-1]-(p[d->nbins-2]+p[0])*0.5;

	for(k=1;k<(d->nbins-1);k++) {
		tmp[k]=p[k]-(p[k-1]+p[k+1])*0.5;
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

/* check whether dataset is sorted in increasing order by GPS time 
   This is useful to reduce influence of phase errors when doing software injections
   and to increase effectiveness of caches.

   Loading SFTs in increasing order by GPS reduces memory requirements as we don't need to reserve space for sorting.
*/
int is_dataset_sorted(DATASET *d)
{
int i;
for(i=1;i<d->free;i++) {
	if(d->gps[i]<=d->gps[i-1])return 0;
	}
return 1;
}

void sort_dataset(DATASET *d)
{
ELEMENT *p;
float *bin_re, *bin_im;
INT64 *gps;
int i, k;

if(d->free==d->size && is_dataset_sorted(d))return;

p=do_alloc(d->free, sizeof(*p));
for(i=0;i<d->free;i++) {
	p[i].index=i;
	p[i].gps=d->gps[i];
	}

qsort(p, d->free, sizeof(*p), element_cmp);

bin_re=do_alloc(d->free*d->nbins, sizeof(*bin_re));
for(i=0;i<d->free;i++) {
	k=p[i].index;
	memcpy(&(bin_re[i*d->nbins]), &(d->re[k*d->nbins]), d->nbins*sizeof(*bin_re));
	}
free(d->re);
d->re=bin_re;

bin_im=do_alloc(d->free*d->nbins, sizeof(*bin_im));
for(i=0;i<d->free;i++) {
	k=p[i].index;
	memcpy(&(bin_im[i*d->nbins]), &(d->im[k*d->nbins]), d->nbins*sizeof(*bin_im));
	}
free(d->im);
d->im=bin_im;

gps=do_alloc(d->free, sizeof(*gps));
for(i=0;i<d->free;i++) {
	k=p[i].index;
	gps[i]=d->gps[k];
	}
free(d->gps);
d->gps=gps;

d->size=d->free;

free(p);
}

void veto_sfts(DATASET *d)
{
int i, j;
float power, *x, *y;
d->sft_veto=do_alloc(d->free, sizeof(*d->sft_veto));
for(i=0;i<d->free;i++) {
	d->sft_veto[i]=0;

	if(d->expTMedians[i]<d->veto_level) {
		d->sft_veto[i]=1;
		continue;
		}

	x=&(d->re[i*d->nbins]);
	y=&(d->im[i*d->nbins]);
	for(j=0;j<d->nbins;j++) {
		power=(*x)*(*x)+(*y)*(*y);
		if(log10(power)>d->TMedians[i]+d->veto_spike_level) {
			d->sft_veto[i]=2;
			break;
			}
		x++;
		y++;
		}
	}
}

typedef struct {
	int point;
	DATASET *d;
	
	} CUTOFF_DATA;

void cutoff_cruncher(int thread_id, CUTOFF_DATA *data)
{
int j, m;
float *tm;
tm=do_alloc(data->d->free, sizeof(*tm));
for(m=0;m<ntotal_polarizations;m++) {
/*	if(patch_grid->band[i]<0){
		d->polarizations[m].patch_CutOff[i]=0.0;
		continue;
		}*/
	if(!do_CutOff) {
		data->d->polarizations[m].patch_CutOff[data->point]=1e30;
		continue;
		}
	for(j=0;j<data->d->free;j++)tm[j]=1.0/(sqrt(data->d->expTMedians[j])*AM_response(j, patch_grid, data->point, data->d->polarizations[m].AM_coeffs));
	data->d->polarizations[m].patch_CutOff[data->point]=FindCutOff(tm, data->d->free);
	}
free(tm);
}

static int validate_dataset(DATASET *d)
{
int i, k;
CUTOFF_DATA *cd;

if(d->validated)return 1;
fprintf(stderr, "Validating dataset \"%s\"\n", d->name);
fprintf(stderr, "Found %d segments\n", d->free);

fprintf(LOG, "dataset %s nsegments : %d\n", d->name, d->free);
fprintf(LOG, "dataset %s detector : %s\n", d->name, d->detector);

if(d->free<1){
	/* free large chunks of ram that might have been reserved */
	if(d->size>0) {
		free(d->gps);
		free(d->re);
		free(d->im);
		d->gps=NULL;
		d->re=NULL;
		d->im=NULL;
		}
	return 0;
	}

if(!strcmp(d->detector, "unknown")) {
	fprintf(stderr, "Each dataset must specify detector being used\n");
	return 0;
	}

get_detector(d->detector);

sort_dataset(d);

fprintf(stderr, "Dataset \"%s\" sorted memory: %g MB\n", d->name, (MEMUSAGE*10.0/(1024.0*1024.0))/10.0);
fprintf(LOG, "Dataset \"%s\" sorted  memory: %g MB\n", d->name, (MEMUSAGE*10.0/(1024.0*1024.0))/10.0);

if(fake_injection) {
	SIGNAL_PARAMS sp;

	fprintf(stderr, "Injecting fake signal.\n");
	fprintf(LOG, "Injecting fake signal.\n");

	fill_signal_params_from_args_info(&sp);

	for(i=0;i<d->free;i++) {
		inject_fake_signal(&sp, d, i);
		}

	fprintf(stderr, "Dataset \"%s\" injected memory: %g MB\n", d->name, (MEMUSAGE*10.0/(1024.0*1024.0))/10.0);
	fprintf(LOG, "Dataset \"%s\" injected  memory: %g MB\n", d->name, (MEMUSAGE*10.0/(1024.0*1024.0))/10.0);
	}

/* compute_power(d); */

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
d->expTMedians_plain=do_alloc(d->free, sizeof(*(d->expTMedians_plain)));

d->FMedians=do_alloc(d->nbins, sizeof(*(d->TMedians)));
d->expFMedians_plain=do_alloc(d->nbins, sizeof(*(d->expFMedians_plain)));

compute_noise_curves(d);

for(i=0;i<d->free;i++){
	d->expTMedians[i]=exp(-M_LN10*2.0*(d->TMedians[i]-d->TMedian));
	d->expTMedians_plain[i]=exp(M_LN10*d->TMedians[i]);
	}
d->expTMedian=exp(-M_LN10*2.0*d->TMedian);

for(i=0;i<d->nbins;i++){
	d->expFMedians_plain[i]=exp(M_LN10*d->FMedians[i]);
	}

fprintf(LOG, "dataset %s TMedian : %f\n", d->name, d->TMedian);

veto_sfts(d);

fprintf(stderr, "Obtaining whole sky AM response for dataset %s\n", d->name);
get_whole_sky_AM_response(d->gps, d->free, args_info.orientation_arg, &(d->AM_coeffs_plus), &(d->AM_coeffs_cross), &(d->AM_coeffs_size));

fprintf(stderr, "Initializing polarizations for dataset %s\n", d->name);
init_polarizations1(d->polarizations, d->AM_coeffs_plus, d->AM_coeffs_cross, d->AM_coeffs_size);

if(args_info.extended_test_arg) {
	precompute_values(fine_grid);
	/* Check AM_response for correctness */
	fprintf(stderr, "Verifying AM response computation for dataset %s\n", d->name);
	for(i=0;i<ntotal_polarizations;i++) {
		verify_whole_sky_AM_response(d->gps, d->free, d->polarizations[i].orientation, 
			fine_grid, 
			d->polarizations[i].AM_coeffs, d->polarizations[i].name);	
		}
	
	free_values(fine_grid);
	}

fprintf(stderr,"Computing cutoff values for dataset %s\n", d->name);
/* compute CutOff values for each patch */

//tm=do_alloc(d->free, sizeof(*tm));

cd=do_alloc(patch_grid->npoints, sizeof(*cd));
reset_jobs_done_ratio();

precompute_values(patch_grid);

for(i=0;i<patch_grid->npoints;i++) {
	cd[i].point=i;
	cd[i].d=d;

	submit_job(cutoff_cruncher, &(cd[i]));
	
	#if 0
	for(m=0;m<ntotal_polarizations;m++) {
/*		if(patch_grid->band[i]<0){
			d->polarizations[m].patch_CutOff[i]=0.0;
			continue;
			}*/
		for(j=0;j<d->free;j++)tm[j]=1.0/(sqrt(d->expTMedians[j])*AM_response(j, patch_grid, i, d->polarizations[m].AM_coeffs));
		d->polarizations[m].patch_CutOff[i]=FindCutOff(tm, d->free);
		}
	#endif
	}

k=0;
while(do_single_job(-1)) {
	if(k % 100 == 0)fprintf(stderr, "% 3.1f ", jobs_done_ratio()*100);
	k++;
	}
wait_for_all_done();

fprintf(stderr, " 100\n");
free_values(patch_grid);
free(cd);

//free(tm);

characterize_dataset(d);
d->validated=1;
fprintf(stderr, "Dataset \"%s\" validated memory: %g MB\n", d->name, (MEMUSAGE*10.0/(1024.0*1024.0))/10.0);
fprintf(LOG, "Dataset \"%s\" validated  memory: %g MB\n", d->name, (MEMUSAGE*10.0/(1024.0*1024.0))/10.0);
return 1;
}

void output_dataset_info(DATASET *d)
{
RGBPic *p;
PLOT *plot;
int m;


fprintf(stderr, "dataset %s allocated size: %f MB\n", d->name, (12.0*d->free*nbins)/(1024.0*1024.0));
fprintf(LOG, "dataset %s allocated size: %f MB\n", d->name, (12.0*d->free*nbins)/(1024.0*1024.0));

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

static int gps_exists(DATASET *d, INT64 gps)
{
int i;
for(i=0;i<d->free;i++) {
	if(d->gps[i]==gps)return 1;
	}
return 0;
}

static int get_geo_range(DATASET *d, char *filename, FILE *fin, long startbin, long count, float *re, float *im, INT64 *gps)
{
REAL8 timebase;
INT4 b, bin_start, nbins;
REAL4 *tmp;
float factor;
long i;

/* read header */	
/* Key */

/* gps */
fread(&b, sizeof(b), 1, fin);
*gps=b;

if(!check_intervals(d->segment_list, *gps)){
	return -1;
	}
if(check_intervals(d->veto_segment_list, *gps)>0){
	return -1;
	}

if(d->no_duplicate_gps && gps_exists(d, *gps)) {
	return -1;
	}


/* skip nsec */
fread(&b, sizeof(b), 1, fin);
/* timebase */
fread(&timebase, sizeof(timebase), 1, fin);

fread(&bin_start, sizeof(bin_start), 1, fin);
fread(&nbins, sizeof(nbins), 1, fin);

tmp=do_alloc(count*2, sizeof(*tmp));

fseek(fin, (startbin-bin_start)*8,SEEK_CUR);
if(fread(tmp,4,count*2,fin)<count*2){
	fprintf(stderr,"Not enough data in file \"%s\" gps=%lld.\n",filename,*gps);
	free(tmp);
	fclose(fin);
	return -3;
	}
/* reverse normalization applied to geo format files */
if(timebase < 0) {
	factor=1.0/args_info.strain_norm_factor_arg; /* make_sft_op did not apply normalization .. */
	fprintf(stderr,"** Timebase is negative, assuming unnormalized data\n");
	fprintf(LOG,"** Timebase is negative, assuming unnormalized data\n");
	} else {
	factor=(0.5*1800.0*16384.0)/(args_info.strain_norm_factor_arg*nbins); /* use fixed normalization for 1800 sec SFTs .. */
	}
factor*=d->dc_factor;
for(i=0;i<count;i++){
	re[i]=tmp[2*i]*factor;
	im[i]=tmp[2*i+1]*factor;
	if(!isfinite(re[i]) || !isfinite(im[i])) {
		free(tmp);
		fprintf(stderr, "Infinite value encountered in file \"%s\"\n", filename);
		return -2;
		}
	}
free(tmp);
return 0;
}

typedef struct {
	INT4 gps_sec;
	INT4 gps_nsec;
	REAL8 tbase;
	INT4 first_frequency_index;
	INT4 nsamples;
	UINT8 crc64;
	CHAR detector[2];
	CHAR padding[2];
	INT4 comment_length;
	} SFTv2_header2;

static int get_sftv2_range(DATASET *d, char *filename, FILE *fin, long startbin, long count, float *re, float *im, INT64 *gps)
{
REAL8 timebase;
INT4 bin_start, nbins;
REAL4 *tmp;
float factor;
long i;
SFTv2_header2 ht;
/* read header */	

if(fread(&ht, sizeof(ht), 1, fin)<1) {
	fprintf(stderr,"Not enough data in file \"%s\" gps=%lld.\n",filename,*gps);
	fclose(fin);
	return -1;
	}

/* fprintf(stderr, "%s gps=%d nsamples=%d\n", filename, ht.gps_sec, ht.nsamples); */

*gps=ht.gps_sec;

if(!check_intervals(d->segment_list, *gps)){
	return -(48+ht.nsamples*8+ht.comment_length);
	}
if(check_intervals(d->veto_segment_list, *gps)>0){
	return -(48+ht.nsamples*8+ht.comment_length);
	}

if(d->no_duplicate_gps && gps_exists(d, *gps)) {
	return -(48+ht.nsamples*8+ht.comment_length);
	}


/* timebase */
timebase=ht.tbase;
bin_start=ht.first_frequency_index;
nbins=ht.nsamples;

if(startbin<bin_start) {
	fprintf(stderr, "Requested frequency range is below region covered by SFT %s gps=%lld bin_start=%d requested=%ld\n", filename, *gps, bin_start, startbin);
	return -(48+ht.nsamples*8+ht.comment_length);
	}

if(startbin+count>bin_start+ht.nsamples) {
	fprintf(stderr, "Requested frequency range is above region covered by SFT %s gps=%lld SFT end=%d requested=%ld\n", filename, *gps, bin_start+ht.nsamples, startbin+count);
	return -(48+ht.nsamples*8+ht.comment_length);
	}

if(fseek(fin, ht.comment_length+(startbin-bin_start)*8,SEEK_CUR)<0) {
	fprintf(stderr,"Not enough data in file \"%s\" gps=%lld.\n",filename,*gps);
	fclose(fin);
	return -1;
	}

tmp=do_alloc(count*2, sizeof(*tmp));

if(fread(tmp,4,count*2,fin)<count*2){
	fprintf(stderr,"Not enough data in file \"%s\" gps=%lld.\n",filename,*gps);
	free(tmp);
	fclose(fin);
	return -1;
	}
/* reverse normalization applied to geo format files */
if(timebase < 0) {
	fprintf(stderr,"** Timebase is negative, ERROR !\n");
	fprintf(LOG,"** Timebase is negative, ERROR !\n");
	exit(-1);
	} else {
	factor=16384.0/args_info.strain_norm_factor_arg; /* fixed normalization for v2 SFTs ? */
	}
factor*=d->dc_factor;
for(i=0;i<count;i++){
	re[i]=tmp[2*i]*factor;
	im[i]=tmp[2*i+1]*factor;
	if(!isfinite(re[i]) || !isfinite(im[i])) {
		free(tmp);
		fprintf(stderr, "Infinite value encountered in file \"%s\"\n", filename);
		return -(48+ht.nsamples*8+ht.comment_length);
		}
	}
free(tmp);
/* return length of this record - so we can seek to the next header */
return (48+ht.nsamples*8+ht.comment_length);
}


static void expand_sft_array(DATASET *d, int count)
{
void *p;

if(count<0)return;
d->size=d->size+count;
fprintf(stderr, "Growing %s SFT array to %f MB count=%d\n", d->name, d->size*d->nbins*2*sizeof(*d->re)/(1024.0*1024.0), d->free);

p=do_alloc(d->size*d->nbins, sizeof(*(d->re)));
if(d->free>0)memcpy(p, d->re, d->free*d->nbins*sizeof(*(d->re)));
free(d->re);
d->re=p;

p=do_alloc(d->size*d->nbins, sizeof(*(d->im)));
if(d->free>0)memcpy(p, d->im, d->free*d->nbins*sizeof(*(d->im)));
free(d->im);
d->im=p;

p=do_alloc(d->size, sizeof(*(d->gps)));
if(d->free>0)memcpy(p, d->gps, d->free*sizeof(*(d->gps)));
free(d->gps);
d->gps=p;
}

static void add_file(DATASET *d, char *filename)
{
FILE *fin;
REAL8 a;
long i;
long retries;
long header_offset;
errno=0;
retries=0;
while((fin=fopen(filename,"r"))==NULL) {
	//if((fin==NULL) && (errno==ENOENT))return -1; /* no such file */
	fprintf(stderr,"Error opening file \"%s\":", filename);
	perror("");
	retries++;
	condor_safe_sleep(args_info.retry_delay_arg);
	}
if(retries>0) {
	fprintf(stderr, "Successfully opened file \"%s\" after %ld attempts.\n", filename, retries);
	}
/* read header */
header_offset=0;
while(1) {
	/* find header */
	if(fseek(fin, header_offset, SEEK_SET)<0) {
		/* no header to find */
		if(!feof(fin)) {
			fprintf(stderr, "*** Error seeking to offset %ld in file %s !!\n", header_offset, filename);
			if(header_offset>=(1<<31)) fprintf(stderr, "*** Possible 2GB filelimit condor issue\n");
			exit(-1);
			}
		return;
		}

	if(d->free>=d->size) {
		expand_sft_array(d, d->size+1);
		}
	d->gps[d->free]=0;

	/* Key */
	a=0.0;
	fread(&a, sizeof(a), 1, fin);
	if(a==1.0) {
		i=get_geo_range(d, filename, fin,  d->first_bin, d->nbins, &(d->re[d->free*d->nbins]), &(d->im[d->free*d->nbins]), &(d->gps[d->free]));
		fclose(fin);
		if(!i)d->free++;
			else if(i< -1)fprintf(stderr, "Skipped file %s (%lld)\n", filename, d->gps[d->free]);
		return;
		} else
	if(a==2.0) {
		i=get_sftv2_range(d, filename, fin,  d->first_bin, d->nbins, &(d->re[d->free*d->nbins]), &(d->im[d->free*d->nbins]), &(d->gps[d->free]));
		if((i<0) && (i> -48)) {
			fclose(fin);
			fprintf(stderr,"Cannot read file \"%s\": wrong endianness or invalid data\n", filename);
			return;
			}
		if(i>0) {
			header_offset+=i;
			d->free++;
			} else 
			header_offset+=-i;
		} else {
		if(!header_offset)fprintf(stderr,"Cannot read file \"%s\": wrong endianness or invalid data\n", filename);
		fclose(fin);
		return;
		}
	}
}

static void d_read_directory(DATASET *dst, char *line, int length)
{
DIR * d;
struct dirent *de;
char *s;
struct timeval start_time, end_time, last_time;
int i, last_i, limit=100;
double delta, delta_avg;
long retries;
int alternative, ai, aj, lock_ok;
s=do_alloc(length+20001, sizeof(*s));


retries=0;
alternative=1;
while(1) {

	locate_arg(line, length, alternative, &ai, &aj);

	if(ai==aj) {
		if(alternative==1) {
			fprintf(stderr, "Could not parse directory line \"%s\"\n", line);
			exit(-1);
			}
		alternative=1;
		continue;
		}

	memcpy(s, &(line[ai]), aj-ai);
	s[aj-ai]=0;

	lock_ok=1;
	if(args_info.enable_dataset_locking_arg && (alternative<=MAX_LOCKS) && (dst->lock_file[alternative-1] != NULL)) {
		lock_ok=acquire_lock(dst->lock_file[alternative-1], 0);
		} else 
	if(args_info.lock_file_given) {
		lock_ok=acquire_lock(args_info.lock_file_arg, 0);
		}

	if(lock_ok) {
		fprintf(stderr, "Reading directory %s\n", s);

		errno=0;
		if((d=opendir(s))!=NULL)break;
		int errsv=errno;

		fprintf(stderr, "Error reading directory %s: %s\n", s, strerror(errsv));
		release_lock();
		}

	alternative++;
	retries++;
	condor_safe_sleep(args_info.retry_delay_arg);
	}
if(retries>0) {
	fprintf(stderr, "Successfully opened directory %s\n", s);
	}

gettimeofday(&start_time, NULL);
last_time=start_time;
i=0;
last_i=0;
while((de=readdir(d))!=NULL) {
	if(de->d_name[0]!='.') {
		/* fprintf(stderr, "Found file \"%s\"\n", de->d_name); */
		snprintf(s+aj-ai, 20000, "/%s", de->d_name);
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

void d_read_file(DATASET *dst, char *line, int length)
{
char *s;
long retries;
int alternative, ai, aj;
struct stat stat_buf;
int lock_ok;

s=do_alloc(length+20001, sizeof(*s));

retries=0;
alternative=1;
while(1) {
		
	locate_arg(line, length, alternative, &ai, &aj);

	if(ai==aj) {
		if(alternative==1) {
			fprintf(stderr, "Could not parse file line \"%s\"\n", line);
			exit(-1);
			}
		alternative=1;
		continue;
		}

	memcpy(s, &(line[ai]), aj-ai);
	s[aj-ai]=0;

	lock_ok=1;
	if(args_info.enable_dataset_locking_arg && (alternative<=MAX_LOCKS) && (dst->lock_file[alternative-1] != NULL)) {
		lock_ok=acquire_lock(dst->lock_file[alternative-1], 0);
		} else 
	if(args_info.lock_file_given) {
		lock_ok=acquire_lock(args_info.lock_file_arg, 0);
		}

	if(lock_ok) {
		fprintf(stderr, "Accessing file %s\n", s);
		errno=0;
		if(stat(s, &stat_buf)==0)break;
		int errsv=errno;

		fprintf(stderr, "Error accessing file %s: %s\n", s, strerror(errsv));
		release_lock();
		}

	alternative++;
	retries++;
	condor_safe_sleep(args_info.retry_delay_arg);
	}
if(retries>0) {
	fprintf(stderr, "Successfully accessed file %s\n", s);
	}

add_file(dst, s);
free(s);
release_lock();
}

long int fill_seed=0;

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
gsl_rng_set(rng, fill_seed);

amp/=args_info.strain_norm_factor_arg;
amp*=d->dc_factor;

for(i=d->free;i<(d->free+count);i++) {
	d->gps[i]=gps_start+(i-d->free)*step;
	
	for(k=0;k<d->nbins;k++) {
		d->re[i*d->nbins+k]=amp*gsl_ran_gaussian(rng, 1.0);
		d->im[i*d->nbins+k]=amp*gsl_ran_gaussian(rng, 1.0);
		}
	}
d->free+=count;

fill_seed=gsl_rng_get(rng);
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
	get_detector(datasets[d_free-1].detector);
	} else 
if(!strncasecmp(line, "expand_sft_array", 16)) {
	int count;
	locate_arg(line, length, 1, &ai, &aj);
	count=atoll(&(line[ai]));
	if(count> 0)expand_sft_array(&(datasets[d_free-1]), count);
	} else 
if(!strncasecmp(line, "expected_sft_count", 18)) {
	int count;
	locate_arg(line, length, 1, &ai, &aj);
	count=atoll(&(line[ai]));
	if(count>= 0) {
		if(count!=datasets[d_free-1].free) {
			fprintf(stderr, "**** ERROR: SFT count mismatch have %d versus expected %d\n", datasets[d_free-1].free, count);
			exit(-1);
			}
		} else {
		if(-count >datasets[d_free-1].free) {
			fprintf(stderr, "**** ERROR: SFT count mismatch have only %d but expected at least %d\n", datasets[d_free-1].free, -count);
			exit(-1);
			}
		}
	} else 
if(!strncasecmp(line, "lock_file", 9)) {
	int i;
	for(i=0;i<MAX_LOCKS;i++) {
		locate_arg(line, length, i+1, &ai, &aj);
		if(datasets[d_free-1].lock_file[i]!=NULL)free(datasets[d_free-1].lock_file[i]);
		if(aj>ai)datasets[d_free-1].lock_file[i]=strndup(&(line[ai]), aj-ai);
			else datasets[d_free-1].lock_file[i]=NULL;
		}
	} else 
if(!strncasecmp(line, "sleep", 5)) {
	locate_arg(line, length, 1, &ai, &aj);
	condor_safe_sleep(atoi(&(line[ai])));
	} else 
if(!strncasecmp(line, "gps_start", 9)) {
	locate_arg(line, length, 1, &ai, &aj);
	datasets[d_free-1].gps_start=atoll(&(line[ai]));
	} else 
if(!strncasecmp(line, "gps_stop", 8)) {
	locate_arg(line, length, 1, &ai, &aj);
	datasets[d_free-1].gps_stop=atoll(&(line[ai]));
	} else 
if(!strncasecmp(line, "block_dc_factor", 15)) {
	datasets[d_free-1].dc_factor_blocked=1;
	if(datasets[d_free-1].dc_factor_touched) {
		fprintf(stderr, "Dataset \"%s\" has dc_factor blocked, but attempted to apply dc correction, exiting\n", datasets[d_free].name);
		fprintf(LOG, "Dataset \"%s\" has dc_factor blocked, but attempted to apply dc correction, exiting\n", datasets[d_free].name);
		exit(-1);
		}
	} else 
if(!strncasecmp(line, "dc_factor", 9)) {
	locate_arg(line, length, 1, &ai, &aj);
	datasets[d_free-1].dc_factor=atof(&(line[ai]));
	datasets[d_free-1].dc_factor_touched=1;
	fprintf(LOG, "Set dc_factor=%f for dataset \"%s\" sft count=%d\n", datasets[d_free-1].dc_factor, datasets[d_free-1].name, datasets[d_free-1].free);
	fprintf(stderr, "Set dc_factor=%f for dataset \"%s\" sft count=%d\n", datasets[d_free-1].dc_factor, datasets[d_free-1].name, datasets[d_free-1].free);
	if(datasets[d_free-1].dc_factor_blocked) {
		fprintf(stderr, "Dataset \"%s\" has dc_factor blocked, but attempted to apply dc correction, exiting\n", datasets[d_free].name);
		fprintf(LOG, "Dataset \"%s\" has dc_factor blocked, but attempted to apply dc correction, exiting\n", datasets[d_free].name);
		exit(-1);
		}
	} else 
if(!strncasecmp(line, "directory", 9)) {
	//locate_arg(line, length, 1, &ai, &aj);
	//d_read_directory(&(datasets[d_free-1]), &(line[ai]), aj-ai);
	d_read_directory(&(datasets[d_free-1]), line, length);
	} else 
if(!strncasecmp(line, "file", 4)) {
	d_read_file(&(datasets[d_free-1]), line, length);
	} else 
if(!strncasecmp(line, "weight", 6)) {
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
if(!strncasecmp(line, "no_duplicate_gps", 15)) {
	locate_arg(line, length, 1, &ai, &aj);
	datasets[d_free-1].no_duplicate_gps=atoi(&(line[ai]));	
	} else 
if(!strncasecmp(line, "veto_level", 10)) {
	locate_arg(line, length, 1, &ai, &aj);
	datasets[d_free-1].veto_level=atof(&(line[ai]));
	} else
if(!strncasecmp(line, "veto_spike_level", 16)) {
	locate_arg(line, length, 1, &ai, &aj);
	datasets[d_free-1].veto_spike_level=atof(&(line[ai]));
	} else
if(!strncasecmp(line, "apply_hanning_filter", 20)) {
	apply_hanning_filter(&(datasets[d_free-1]));
	} else
if(!strncasecmp(line, "gaussian_fill_seed", 18)) {
	locate_arg(line, length, 1, &ai, &aj);
	sscanf(&(line[ai]), "%ld", &fill_seed);
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
	} else
if(!strncasecmp(line, "inject_wandering_line1", 22)) {
	INT64 gps_start, gps_stop;
	double strain, frequency, ref_time, spindown, spread, omega;
	double f, phi;
	gsl_rng *rng=NULL;
	int bin;
	int i;
	int once=1;
	DATASET *d=&(datasets[d_free-1]);

	locate_arg(line, length, 1, &ai, &aj);
	sscanf(&(line[ai]), "%lld", &gps_start);

	locate_arg(line, length, 2, &ai, &aj);
	sscanf(&(line[ai]), "%lld", &gps_stop);

	locate_arg(line, length, 3, &ai, &aj);
	sscanf(&(line[ai]), "%lg", &ref_time);

	locate_arg(line, length, 4, &ai, &aj);
	sscanf(&(line[ai]), "%lg", &strain);

	locate_arg(line, length, 5, &ai, &aj);
	sscanf(&(line[ai]), "%lg", &frequency);

	locate_arg(line, length, 6, &ai, &aj);
	sscanf(&(line[ai]), "%lg", &spindown);

	locate_arg(line, length, 7, &ai, &aj);
	sscanf(&(line[ai]), "%lg", &spread);

	locate_arg(line, length, 8, &ai, &aj);
	sscanf(&(line[ai]), "%lg", &omega);

	fprintf(stderr, "Injecting fake wandering line artifact 1 for period %lld-%lld ref_time=%lf strain=%lg frequency=%lf spindown=%lg spread=%lg omega=%lg\n", gps_start, gps_stop,
		ref_time,
		strain,
		frequency,
		spindown,
		spread,
		omega);

	fprintf(LOG, "Injecting fake wandering line artifact 1 for period %lld-%lld ref_time=%lf strain=%lg frequency=%lf spindown=%lg spread=%lg omega=%lg\n", gps_start, gps_stop,
		ref_time,
		strain,
		frequency,
		spindown,
		spread,
		omega);

	/* Inject wandering artifact into a single bin in each SFT */

	rng=gsl_rng_alloc(gsl_rng_default);
	gsl_rng_set(rng, fill_seed);

	for(i=0;i<d->free;i++) {
		if(gps_start>0 && d->gps[i]<gps_start)continue;
		if(gps_stop>0 && d->gps[i]>=gps_stop)continue;

		f=frequency+spindown*(d->gps[i]-ref_time)+spread*sin(omega*(d->gps[i]-ref_time));
		phi=gsl_ran_flat(rng, 0, 2*M_PI);
		bin=round(f*d->coherence_time-d->first_bin);

		if(bin<0 || bin>=d->nbins) {
			if(once) {
				fprintf(stderr, "Requested to inject power outside loaded range, skipping affected SFTs\n");
				fprintf(LOG, "Requested to inject power outside loaded range, skipping affected SFTs\n");
				once=0;
				}
			continue;
			}

		d->re[i*d->nbins+bin]+=0.5*strain*cos(phi)*d->coherence_time*16384.0/args_info.strain_norm_factor_arg;
		d->im[i*d->nbins+bin]+=0.5*strain*sin(phi)*d->coherence_time*16384.0/args_info.strain_norm_factor_arg;
		}

	fill_seed=gsl_rng_get(rng);
	gsl_rng_free(rng);
	} else
if(!strncasecmp(line, "inject_wandering_line2", 22)) {
	INT64 gps_start, gps_stop;
	double strain, frequency, ref_time, spindown, spread;
	int count;
	double f, phi;
	gsl_rng *rng=NULL;
	int bin;
	int i, k;
	int once=1;
	DATASET *d=&(datasets[d_free-1]);

	locate_arg(line, length, 1, &ai, &aj);
	sscanf(&(line[ai]), "%lld", &gps_start);

	locate_arg(line, length, 2, &ai, &aj);
	sscanf(&(line[ai]), "%lld", &gps_stop);

	locate_arg(line, length, 3, &ai, &aj);
	sscanf(&(line[ai]), "%lg", &ref_time);

	locate_arg(line, length, 4, &ai, &aj);
	sscanf(&(line[ai]), "%lg", &strain);

	locate_arg(line, length, 5, &ai, &aj);
	sscanf(&(line[ai]), "%lg", &frequency);

	locate_arg(line, length, 6, &ai, &aj);
	sscanf(&(line[ai]), "%lg", &spindown);

	locate_arg(line, length, 7, &ai, &aj);
	sscanf(&(line[ai]), "%lg", &spread);

	locate_arg(line, length, 8, &ai, &aj);
	sscanf(&(line[ai]), "%d", &count);

	fprintf(stderr, "Injecting fake wandering line artifact 2 for period %lld-%lld ref_time=%lf strain=%lg frequency=%lf spindown=%lg spread=%lg count=%d\n", gps_start, gps_stop,
		ref_time,
		strain,
		frequency,
		spindown,
		spread,
		count);

	fprintf(LOG, "Injecting fake wandering line artifact 2 for period %lld-%lld ref_time=%lf strain=%lg frequency=%lf spindown=%lg spread=%lg count=%d\n", gps_start, gps_stop,
		ref_time,
		strain,
		frequency,
		spindown,
		spread,
		count);

	rng=gsl_rng_alloc(gsl_rng_default);
	gsl_rng_set(rng, fill_seed);

	/* Inject wandering artifact that affects count bins (with replacement) in each SFT - this simulates artifact that randomly changes frequency multiple times during SFT period. */

	for(i=0;i<d->free;i++) {
		if(gps_start>0 && d->gps[i]<gps_start)continue;
		if(gps_stop>0 && d->gps[i]>=gps_stop)continue;

		for(k=0;k<count;k++) {
			f=frequency+spindown*(d->gps[i]-ref_time)+gsl_ran_gaussian(rng, spread);
			phi=gsl_ran_flat(rng, 0, 2*M_PI);
			bin=round(f*d->coherence_time-d->first_bin);
	
			if(bin<0 || bin>=d->nbins) {
				if(once) {
					fprintf(stderr, "Requested to inject power outside loaded range, skipping affected SFTs\n");
					fprintf(LOG, "Requested to inject power outside loaded range, skipping affected SFTs\n");
					once=0;
					}
				continue;
				}
	
			d->re[i*d->nbins+bin]+=0.5*strain*cos(phi)*d->coherence_time*16384.0/args_info.strain_norm_factor_arg;
			d->im[i*d->nbins+bin]+=0.5*strain*sin(phi)*d->coherence_time*16384.0/args_info.strain_norm_factor_arg;
			}
		}

	fill_seed=gsl_rng_get(rng);
	gsl_rng_free(rng);
	} else
if(!strncasecmp(line, "inject_cw_signal", 16)) {
	INT64 gps_start, gps_stop;
	SIGNAL_PARAMS sp;
	int i;
	DATASET *d=&(datasets[d_free-1]);

	/* Fill with defaults from command line */
	fill_signal_params_with_defaults(&sp);
	
	
	locate_arg(line, length, 1, &ai, &aj);
	sscanf(&(line[ai]), "%lld", &gps_start);

	locate_arg(line, length, 2, &ai, &aj);
	sscanf(&(line[ai]), "%lld", &gps_stop);

	locate_arg(line, length, 3, &ai, &aj);
	sscanf(&(line[ai]), "%lg", &sp.ref_time);

	locate_arg(line, length, 4, &ai, &aj);
	sscanf(&(line[ai]), "%lg", &sp.strain);

	locate_arg(line, length, 5, &ai, &aj);
	sscanf(&(line[ai]), "%lg", &sp.freq);

	locate_arg(line, length, 6, &ai, &aj);
	sscanf(&(line[ai]), "%lg", &sp.spindown);

	locate_arg(line, length, 7, &ai, &aj);
	sscanf(&(line[ai]), "%lg", &sp.ra);

	locate_arg(line, length, 8, &ai, &aj);
	sscanf(&(line[ai]), "%lg", &sp.dec);

	locate_arg(line, length, 9, &ai, &aj);
	sscanf(&(line[ai]), "%lg", &sp.iota);

	locate_arg(line, length, 10, &ai, &aj);
	sscanf(&(line[ai]), "%lg", &sp.psi);

	locate_arg(line, length, 11, &ai, &aj);
	sscanf(&(line[ai]), "%lg", &sp.phi);

	fprintf(stderr, "Injecting fake signal for period %lld-%lld ref_time=%lf strain=%lg frequency=%lf spindown=%lg ra=%lf dec=%lf iota=%lf psi=%lf phi=%lf dInv=%lg freq_mod_depth=%lg freq_mod_freq=%lg freq_mod_phase=%lf phase_mod_depth=%lg phase_mod_freq=%lg phase_mod_phase=%lf\n", gps_start, gps_stop,
			sp.ref_time,
			sp.strain,
			sp.freq,
			sp.spindown,
			sp.ra,
			sp.dec,
			sp.iota,
			sp.psi,
			sp.phi,
			sp.dInv,
			sp.freq_modulation_depth,
			sp.freq_modulation_freq,
			sp.freq_modulation_phase,
			sp.phase_modulation_depth,
			sp.phase_modulation_freq,
			sp.phase_modulation_phase);
	fprintf(LOG, "Injecting fake signal for period %lld-%lld ref_time=%lf strain=%lg frequency=%lf spindown=%lg ra=%lf dec=%lf iota=%lf psi=%lf phi=%lf dInv=%lg freq_mod_depth=%lg freq_mod_freq=%lg freq_mod_phase=%lf phase_mod_depth=%lg phase_mod_freq=%lg phase_mod_phase=%lf\n", gps_start, gps_stop,
			sp.ref_time,
			sp.strain,
			sp.freq,
			sp.spindown,
			sp.ra,
			sp.dec,
			sp.iota,
			sp.psi,
			sp.phi,
			sp.dInv,
			sp.freq_modulation_depth,
			sp.freq_modulation_freq,
			sp.freq_modulation_phase,
			sp.phase_modulation_depth,
			sp.phase_modulation_freq,
			sp.phase_modulation_phase);

	precompute_signal_params(&sp);

	sort_dataset(d);

	for(i=0;i<d->free;i++) {
		if(gps_start>0 && d->gps[i]<gps_start)continue;
		if(gps_stop>0 && d->gps[i]>=gps_stop)continue;

		inject_fake_signal(&sp, d, i);
		}

	} else
if(!strncasecmp(line, "inject_cw_signal2", 16)) {
	INT64 gps_start, gps_stop;
	SIGNAL_PARAMS sp;
	int i;
	DATASET *d=&(datasets[d_free-1]);

	/* Fill with defaults from command line */
	fill_signal_params_with_defaults(&sp);
	
	
	locate_arg(line, length, 1, &ai, &aj);
	sscanf(&(line[ai]), "%lld", &gps_start);

	locate_arg(line, length, 2, &ai, &aj);
	sscanf(&(line[ai]), "%lld", &gps_stop);

	locate_arg(line, length, 3, &ai, &aj);
	sscanf(&(line[ai]), "%lg", &sp.ref_time);

	locate_arg(line, length, 4, &ai, &aj);
	sscanf(&(line[ai]), "%lg", &sp.strain);

	locate_arg(line, length, 5, &ai, &aj);
	sscanf(&(line[ai]), "%lg", &sp.freq);

	locate_arg(line, length, 6, &ai, &aj);
	sscanf(&(line[ai]), "%lg", &sp.spindown);

	locate_arg(line, length, 7, &ai, &aj);
	sscanf(&(line[ai]), "%lg", &sp.ra);

	locate_arg(line, length, 8, &ai, &aj);
	sscanf(&(line[ai]), "%lg", &sp.dec);

	locate_arg(line, length, 9, &ai, &aj);
	sscanf(&(line[ai]), "%lg", &sp.iota);

	locate_arg(line, length, 10, &ai, &aj);
	sscanf(&(line[ai]), "%lg", &sp.psi);

	locate_arg(line, length, 11, &ai, &aj);
	sscanf(&(line[ai]), "%lg", &sp.phi);

	locate_arg(line, length, 12, &ai, &aj);
	sscanf(&(line[ai]), "%lg", &sp.dInv);

	locate_arg(line, length, 13, &ai, &aj);
	sscanf(&(line[ai]), "%lg", &sp.freq_modulation_depth);

	locate_arg(line, length, 14, &ai, &aj);
	sscanf(&(line[ai]), "%lg", &sp.freq_modulation_freq);

	locate_arg(line, length, 15, &ai, &aj);
	sscanf(&(line[ai]), "%lg", &sp.freq_modulation_phase);

	locate_arg(line, length, 16, &ai, &aj);
	sscanf(&(line[ai]), "%lg", &sp.phase_modulation_depth);

	locate_arg(line, length, 17, &ai, &aj);
	sscanf(&(line[ai]), "%lg", &sp.phase_modulation_freq);

	locate_arg(line, length, 18, &ai, &aj);
	sscanf(&(line[ai]), "%lg", &sp.phase_modulation_phase);

	fprintf(stderr, "Injecting fake signal (mode 2) for period %lld-%lld ref_time=%lf strain=%lg frequency=%lf spindown=%lg ra=%lf dec=%lf iota=%lf psi=%lf phi=%lf dInv=%lg freq_mod_depth=%lg freq_mod_freq=%lg freq_mod_phase=%lf phase_mod_depth=%lg phase_mod_freq=%lg phase_mod_phase=%lf\n", gps_start, gps_stop,
			sp.ref_time,
			sp.strain,
			sp.freq,
			sp.spindown,
			sp.ra,
			sp.dec,
			sp.iota,
			sp.psi,
			sp.phi,
			sp.dInv,
			sp.freq_modulation_depth,
			sp.freq_modulation_freq,
			sp.freq_modulation_phase,
			sp.phase_modulation_depth,
			sp.phase_modulation_freq,
			sp.phase_modulation_phase);
	fprintf(LOG, "Injecting fake signal for period %lld-%lld ref_time=%lf strain=%lg frequency=%lf spindown=%lg ra=%lf dec=%lf iota=%lf psi=%lf phi=%lf dInv=%lg freq_mod_depth=%lg freq_mod_freq=%lg freq_mod_phase=%lf phase_mod_depth=%lg phase_mod_freq=%lg phase_mod_phase=%lf\n", gps_start, gps_stop,
			sp.ref_time,
			sp.strain,
			sp.freq,
			sp.spindown,
			sp.ra,
			sp.dec,
			sp.iota,
			sp.psi,
			sp.phi,
			sp.dInv,
			sp.freq_modulation_depth,
			sp.freq_modulation_freq,
			sp.freq_modulation_phase,
			sp.phase_modulation_depth,
			sp.phase_modulation_freq,
			sp.phase_modulation_phase);

	precompute_signal_params(&sp);

	sort_dataset(d);

	for(i=0;i<d->free;i++) {
		if(gps_start>0 && d->gps[i]<gps_start)continue;
		if(gps_stop>0 && d->gps[i]>=gps_stop)continue;

		inject_fake_signal(&sp, d, i);
		}

	} else {
	fprintf(stderr, "*** Could not parse line: \n%*s\n *** Exiting.\n", length, line);
	exit(-1);
	}	
}

void load_datasets(char *definition)
{
char *p=definition;
int k;

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
free(s);
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

long vetoed_segments(void)
{
int i, j;
long total=0;
for(i=0;i<d_free;i++) {
	for(j=0;j<datasets[i].free;j++)
		if(datasets[i].sft_veto[j]!=0)total++;
	}
return(total);
}

static float compute_median(float *firstbin, int step, int count)
{
float *tmp;
int i;
tmp=aligned_alloca(count*sizeof(float));
for(i=0;i<count;i++)tmp[i]=firstbin[i*step];
sort_floats(tmp, count);
if(!(count & 1))return (tmp[count>>1]+tmp[(count>>1)-1])*0.5;
return tmp[count>>1];
}

static void test_compute_median(void)
{
float test1[11]={-0.256, -0.096, 1.357, 0.442, -0.728, 1.084, 0.178, -1.527, 0.333, 0.651, 0.809};
float test2[12]={-0.308, -2.124, 0.153, 0.314, 1.91, -1.646, 3.299, -0.226, 0.201, 1.22, 0.122, 0.751};
float err;

err=compute_median(test1, 1, 11)-0.333;
fprintf(stderr, "compute_median_test1: %g\n", err); 
fprintf(LOG, "compute_median_test1: %g\n", err);
if(fabs(err)>1e-6)exit(-1); 

err=compute_median(test2, 1, 12)-0.177;
fprintf(stderr, "compute_median_test2: %g\n", err); 
fprintf(LOG, "compute_median_test2: %g\n", err);
if(fabs(err)>1e-6)exit(-1); 
}

void compute_noise_curves(DATASET *dataset)
{
float *tmp;
float *x, *y, *t;
float a;
int i,j;
double b, b_initial;
int nsegments=dataset->free;
int nbins=dataset->nbins;
HISTOGRAM *hist;

tmp=do_alloc(nsegments*nbins,sizeof(float));
for(i=0;i<nsegments;i++) {
	t=&(tmp[i*nbins]);
	x=&(dataset->re[i*nbins]);
	y=&(dataset->im[i*nbins]);
	for(j=0;j<nbins;j++) {
		t[j]=log10(x[j]*x[j]+y[j]*y[j]);
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
fprintf(stderr,"%g\n",b);
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
	//if(b<(b_initial*(1E-16)))break;
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
float x, y;
float *tmp;
/*
for(i=0;i<d->free;i++){
	d->expTMedians[i]=exp(-M_LN10*2.0*(d->TMedians[i]-d->TMedian));
	}
d->expTMedian=exp(-M_LN10*2.0*d->TMedian);
*/

fprintf(stderr, "%s dc_factor_blocked: %d\n", d->name, d->dc_factor_blocked);
fprintf(LOG, "%s dc_factor_blocked: %d\n", d->name, d->dc_factor_blocked);
fprintf(stderr, "%s dc_factor_touched: %d\n", d->name, d->dc_factor_touched);
fprintf(LOG, "%s dc_factor_touched: %d\n", d->name, d->dc_factor_touched);

fprintf(stderr, "%s SFTs veto level: %f\n", d->name, d->veto_level);
fprintf(LOG, "%s SFTs veto level: %f\n", d->name, d->veto_level);
fprintf(stderr, "%s SFTs veto spike level: %f\n", d->name, d->veto_spike_level);
fprintf(LOG, "%s SFTs veto spike level: %f\n", d->name, d->veto_spike_level);



count=0;
for(i=0;i<d->free;i++)
	if(d->sft_veto[i]) {
		if(count<args_info.max_sft_report_arg) {
			fprintf(stderr, "%s vetoed SFT: %lld %d\n", d->name, d->gps[i], d->sft_veto[i]+0);
			fprintf(LOG, "%s vetoed SFT: %lld %d\n", d->name, d->gps[i], d->sft_veto[i]+0);
			}
		count++;
		}

fprintf(stderr, "%s SFTs veto count: %d\n", d->name, count);
fprintf(LOG, "%s SFTs veto count: %d\n", d->name, count);


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
	if(d->sft_veto[i])continue;

	w=d->expTMedians[i];
	for(j=0;j<d->nbins;j++){
		x=d->re[i*d->nbins+j];
		y=d->im[i*d->nbins+j];

		a=x*x+y*y;

		if(args_info.subtract_background_arg) {
			a-=d->expTMedians_plain[i]*d->expFMedians_plain[j];
			}

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
	if(d->sft_veto[i])continue;
	if(d->TMedians[i]>=CutOff)continue;
	w=d->expTMedians[i];
	count++;
	for(j=0;j<d->nbins;j++){
		x=d->re[i*d->nbins+j];
		y=d->im[i*d->nbins+j];
		a=x*x+y*y;

		if(args_info.subtract_background_arg) {
			/* the line below is faster a-=exp(M_LN10*(d->TMedians[i]+d->FMedians[j])); */
			a-=d->expTMedians_plain[i]*d->expFMedians_plain[j];
			}

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
float min_tm=0;
for(i=0;i<d_free;i++){
	if(!i || datasets[i].TMedian<min_tm)min_tm=datasets[i].TMedian;
	}
return(exp(-M_LN10*min_tm));
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
	d->weight/=b*b;
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

#define TRACE_EFFECTIVE 0

float effective_weight_ratio(float target_ra, float target_dec, float source_ra, float source_dec, float source_spindown, float bin_tolerance, float spindown_tolerance)
{
int i, k;
double total_weight=0.0, w, fdiff, f1, f2, c1, c2;
DATASET *d;
float e1[26], e2[26], ed[26];
double offset, fdot;
double timebase=max_gps()-min_gps()+1800.0;
double inv_timebase=1.0/timebase;
double max_fdot;
double delta_spindown=spindown-source_spindown;
#if TRACE_EFFECTIVE
static once=1;
#endif

precompute_am_constants(e1, source_ra, source_dec);
precompute_am_constants(e2, target_ra, target_dec);

for(i=0;i<26;i++)ed[i]=e2[i]-e1[i];

/* fit first */
f1=0;
f2=0;
c1=0;
c2=0;
for(i=0;i<d_free;i++) {
	d=&(datasets[i]);
	for(k=0;k<d->free;k++) {
		w=d->expTMedians[k]*d->weight;
		total_weight+=w;

		fdiff=(first_bin+nbins*0.5)*(
				ed[0]*d->detector_velocity[3*k+0]+
				ed[1]*d->detector_velocity[3*k+1]+
				ed[2]*d->detector_velocity[3*k+2]
				)*args_info.doppler_multiplier_arg+(d->gps[k]-spindown_start+d->coherence_time*0.5)*d->coherence_time*delta_spindown;
		f1+=fdiff*w;
		#if TRACE_EFFECTIVE
		if(once)fprintf(LOG, "%8f %8f %8f\n",w,  (d->gps[k]-spindown_start+d->coherence_time*0.5)*inv_timebase, fdiff);
		#endif
		w*=(d->gps[k]-spindown_start+d->coherence_time*0.5)*inv_timebase;
		f2+=fdiff*w;
		c1+=w;
		c2+=w*(d->gps[k]-spindown_start+d->coherence_time*0.5)*inv_timebase;
		}
	}
f1/=total_weight;
f2/=total_weight;
c1/=total_weight;
c2/=total_weight;

/* Note: the denominator will only be zero if we have only one SFT - in which case this has issues working anyway */
offset=(c2*f1-c1*f2)/(c2-c1*c1);
fdot=(-c1*f1+f2)/(c2-c1*c1);

/* make it a width parameter */
bin_tolerance*=0.5;

max_fdot=0.5*fabs(spindown_tolerance)*timebase*datasets[0].coherence_time;
if(fdot>max_fdot) {
	offset=offset+(fdot-max_fdot)*c1;
	fdot=max_fdot;
	}
if(fdot< -max_fdot) {
	offset=offset+(fdot+max_fdot)*c1;
	fdot=-max_fdot;
	}

#if TRACE_EFFECTIVE
if(once)fprintf(LOG, "%g %g %g %g %g\n", f1, f2, total_weight, offset, fdot);
once=0;
#endif
/* now compare */
f1=0;
for(i=0;i<d_free;i++) {
	d=&(datasets[i]);
	for(k=0;k<d->free;k++) {
		w=d->expTMedians[k]*d->weight;
		total_weight+=w;

		fdiff=(first_bin+nbins*0.5)*(
				ed[0]*d->detector_velocity[3*k+0]+
				ed[1]*d->detector_velocity[3*k+1]+
				ed[2]*d->detector_velocity[3*k+2]
				)*args_info.doppler_multiplier_arg
			+(d->gps[k]-spindown_start+d->coherence_time*0.5)*d->coherence_time*delta_spindown
			-offset-fdot*(d->gps[k]-spindown_start+d->coherence_time*0.5)*inv_timebase;
		if(fabs(fdiff)<bin_tolerance)f1+=w;
		}
	}
/* fprintf(stderr, "%g %g\n", f1, total_weight); */
return f1/total_weight;
}

float stationary_effective_weight_ratio(float target_ra, float target_dec, float bin_tolerance)
{
int i, k;
double total_weight=0.0, w, fdiff, f1;
DATASET *d;
float e[26];
double offset;

precompute_am_constants(e, target_ra, target_dec);

/* fit first */
f1=0;
for(i=0;i<d_free;i++) {
	d=&(datasets[i]);
	for(k=0;k<d->free;k++) {
		w=d->expTMedians[k]*d->weight;
		total_weight+=w;

		fdiff=(first_bin+nbins*0.5)*(
				e[0]*d->detector_velocity[3*k+0]+
				e[1]*d->detector_velocity[3*k+1]+
				e[2]*d->detector_velocity[3*k+2]
				)*args_info.doppler_multiplier_arg
			+d->coherence_time*spindown*(d->gps[k]-spindown_start+d->coherence_time*0.5);
		f1+=fdiff*w;
		}
	}
offset=f1/total_weight;

bin_tolerance*=0.5; /* make it a width */

/* fprintf(stderr, "%g %g %g %g %g\n", f1, f2, total_weight, offset, fdot); */
/* now compare */
f1=0;
for(i=0;i<d_free;i++) {
	d=&(datasets[i]);
	for(k=0;k<d->free;k++) {
		w=d->expTMedians[k]*d->weight;
		total_weight+=w;

		fdiff=(first_bin+nbins*0.5)*(
				e[0]*d->detector_velocity[3*k+0]+
				e[1]*d->detector_velocity[3*k+1]+
				e[2]*d->detector_velocity[3*k+2]
				)*args_info.doppler_multiplier_arg
			+d->coherence_time*spindown*(d->gps[k]-spindown_start+d->coherence_time*0.5)
			-offset;

		if(fabs(fdiff)<bin_tolerance)f1+=w;
		}
	}
/* fprintf(stderr, "%g %g\n", f1, total_weight); */
return f1/total_weight;
}

void inject_signal(SIGNAL_PARAMS *sp)
{
int i,j;
DATASET *d;

for(i=0;i<d_free;i++){

	d=&(datasets[i]);
	
	for(j=0;j<d->free;j++) {
		inject_fake_signal(sp, d, i);
		}

	}
}

void apply_phase(char *phase_string)
{
int i,j;
DATASET *d;
float x, y, a, b, phase;
regex_t regexp;
char *match;

match=strdup(phase_string);
for(i=strlen(match)-1; (i>0) && (match[i]==' ' || match[i]==0 || match[i]=='\n' || match[i]=='\t');i--);
for(; (i>0) && (match[i]!=' ' && match[i]!='\t');i--);
match[i]=0;

fprintf(stderr, "Applying phase of %s degrees to datasets matching \"%s\"\n", &(match[i+1]), match);

if(regcomp(&regexp, match, REG_EXTENDED | REG_NOSUB)){
	fprintf(stderr,"Cannot compile \"%s\" in extra phase specifier \"%s\"\n", match, phase_string);
	exit(-1);
	}
	
phase=atof(&(match[i+1]));

a=cos(phase*M_PI/180.0);
b=sin(phase*M_PI/180.0);

for(i=0;i<d_free;i++){

	d=&(datasets[i]);
	
	if(regexec(&regexp, d->name, 0, NULL, 0))continue;
	
	fprintf(stderr, "Applying phase %f to dataset %s\n", phase, d->name);
	fprintf(LOG, "Applying phase %f to dataset %s\n", phase, d->name);
	
	for(j=0;j<d->free;j++) {
		x=d->re[j]*a-d->im[j]*b;
		y=d->re[j]*b+d->im[j]*a;
		
		#if 0
		fprintf(stderr, "%f+%fi --> %f+%fi\n", d->re[j], d->im[j], x, y);
		#endif

		d->re[j]=x;
		d->im[j]=y;
		
		}

	}

regfree(&regexp);
free(match);
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
fprintf(fout, "Dataset\tDetector\tGPS\tdopplerX\tdopplerY\tdopplerZ");
for(k=0;k<nbins;k++) fprintf(fout, "\tR%d", k+first_bin);
for(k=0;k<nbins;k++) fprintf(fout, "\tI%d", k+first_bin);
for(k=0;k<nbins;k++) fprintf(fout, "\tP%d", k+first_bin);
fprintf(fout, "\n");

for(i=0;i<d_free;i++) {
	d=&(datasets[i]);
	for(j=0;j<d->free;j++) {
		fprintf(fout, "%s\t%s\t%lld\t%8g\t%8g\t%8g", d->name, d->detector, d->gps[j],
			d->detector_velocity[3*j],
			d->detector_velocity[3*j+1],
			d->detector_velocity[3*j+2]
			);
		for(k=0;k<d->nbins;k++) {
			fprintf(fout, "\t%8g", d->re[j*d->nbins+k]*args_info.strain_norm_factor_arg/(0.5*1800.0*16384.0));
			}
		for(k=0;k<d->nbins;k++) {
			fprintf(fout, "\t%8g",d->im[j*d->nbins+k]*args_info.strain_norm_factor_arg/(0.5*1800.0*16384.0));
			}
		for(k=0;k<d->nbins;k++) {
			x=d->re[j*d->nbins+k]*args_info.strain_norm_factor_arg/(0.5*1800.0*16384.0);
			y=d->im[j*d->nbins+k]*args_info.strain_norm_factor_arg/(0.5*1800.0*16384.0);
			fprintf(fout, "\t%8g", x*x+y*y);
			}
		fprintf(fout, "\n");
		}
	}
fclose(fout);
}

void sftv2_dump_datasets(char *directory) 
{
FILE *fout;
FILE *fdst;
int i,j,k;
DATASET *d;
float x,y;
char filename[1000];
char *path;
double key=2.0;
SFTv2_header2 sfth;

fprintf(stderr, "Dumping datasets into %s\n", directory);
mkdir(directory, ~0);

path=do_alloc(strlen(directory)+30, 1);
sprintf(path, "%s/dataset.dst", directory);
fdst=fopen(path, "w");

if(fdst==NULL) {
	perror(path);
	free(path);
	return;
	}
free(path);

for(i=0;i<d_free;i++) {
	d=&(datasets[i]);

	if(d->free<1)continue;

	snprintf(filename+1, 999, "-%s-%lld-%lld.sft", d->name, d->gps[0], (d->gps[d->free-1]-d->gps[0]+1800));
	if(!strcasecmp(d->detector, "LHO"))filename[0]='H';
		else
	if(!strcasecmp(d->detector, "LLO"))filename[0]='L';
		else filename[0]='u';

	fprintf(stderr, "\twriting %s\n", filename);

	fprintf(fdst, "new_dataset \"%s\"\n", d->name);
	fprintf(fdst, "detector \"%s\"\n\n", d->detector);
	if(d->dc_factor_touched || d->dc_factor_blocked)fprintf(fdst, "block_dc_factor\n\n");

	fprintf(fdst, "# weight %f\n", d->weight);
	fprintf(fdst, "# firstbin %d\n", d->first_bin);
	fprintf(fdst, "# nbins %d\n", d->nbins);
	fprintf(fdst, "# gps_start %lld\n", d->gps[0]);
	fprintf(fdst, "# gps_stop %lld\n\n", d->gps[d->free-1]);
	fprintf(fdst, "# dc_factor_touched %d\n\n", d->dc_factor_touched);
	fprintf(fdst, "# dc_factor_blocked %d\n\n", d->dc_factor_blocked);

	path=do_alloc(strlen(directory)+strlen(filename)+2, 1);
	sprintf(path, "%s/%s", directory, filename);

	fprintf(fdst, "file \"%s\"\n\n", path);

	fout=fopen(path, "w");
	if(fout==NULL) {
		perror(path);
		fclose(fdst);
		free(path);
		return;
		}
	free(path);

	for(j=0;j<d->free;j++) {

		key=2.0;
		sfth.gps_sec=d->gps[j];
		sfth.gps_nsec=0;
		sfth.tbase=d->coherence_time;
		sfth.first_frequency_index=d->first_bin;
		sfth.nsamples=d->nbins;
		sfth.crc64=0;
		sfth.detector[0]=filename[0];
		sfth.detector[1]='1'; /* TODO: we can't tell H2 from H1 here */
		sfth.padding[0]=0;
		sfth.padding[1]=0;
		sfth.comment_length=0;

		fwrite(&key, sizeof(key), 1, fout);
		fwrite(&sfth, sizeof(sfth), 1, fout);

		for(k=0;k<d->nbins;k++) {
			x=d->re[j*d->nbins+k]*args_info.strain_norm_factor_arg/(16384.0);
			y=d->im[j*d->nbins+k]*args_info.strain_norm_factor_arg/(16384.0);
			fwrite(&x, sizeof(x), 1, fout);
			fwrite(&y, sizeof(y), 1, fout);
			}
		}
	fclose(fout);
	}
fclose(fdst);
}

void verify_dataset_whole_sky_AM_response(void)
{
DATASET *d;
int i,j;
/* Check AM_response for correctness */
for(j=0;j<d_free;j++) {
	d=&(datasets[j]);
	fprintf(stderr, "Verifying AM response computation for dataset %s\n", d->name);
	get_detector(d->detector);
	for(i=0;i<ntotal_polarizations;i++){
		verify_whole_sky_AM_response(d->gps, d->free, d->polarizations[i].orientation, 
			fine_grid, 
			d->polarizations[i].AM_coeffs, d->polarizations[i].name);	
		}
	}
}

void fake_dataset_test(void)
{
DATASET ds;
long int old_fill_seed=fill_seed;
int old_fake_injection=fake_injection;
int pass=1;

memset(&ds, 0, sizeof(ds));
init_dataset(&ds);
ds.name="internal_test";
ds.detector="LHO";
ds.coherence_time=1800.0;

fill_seed=0;
gaussian_fill(&ds, 793161250, 900, 1000, 1e-24*ds.coherence_time*16384.0);
fill_seed=old_fill_seed;
fake_injection=0;
validate_dataset(&ds);
fake_injection=old_fake_injection;
if(fabs(ds.TMedian-7.08164)>1e-2) {
	fprintf(stderr, "ERROR test TMedian has unexpected value: %g\n", ds.TMedian); 
	fprintf(LOG, "ERROR test TMedian has unexpected value: %g\n", ds.TMedian); 
	pass=0;
	}
if(fabs(ds.average_detector_velocity[0]+3.59906e-05)>1e-10 ||
   fabs(ds.average_detector_velocity[1]+8.57193e-05)>1e-10 ||
   fabs(ds.average_detector_velocity[2]+3.71474e-05)>1e-10 ) {
	fprintf(stderr, "ERROR test average detector velocity has unexpected value: (%g, %g, %g)\n", 
		ds.average_detector_velocity[0],
		ds.average_detector_velocity[1],
		ds.average_detector_velocity[2]
		); 
	fprintf(LOG, "ERROR test average detector velocity has unexpected value: (%g, %g, %g)\n", 
		ds.average_detector_velocity[0],
		ds.average_detector_velocity[1],
		ds.average_detector_velocity[2]
		); 
	pass=0;
	}
if(!pass)exit(-1);
}

void test_datasets(void)
{
test_compute_median();
test_inject_fake_signal09();

if(args_info.extended_test_arg) {
	fake_dataset_test();
	test_inject_fake_signal05();
	}
}

void init_datasets(void)
{
fill_seed=args_info.initial_dataset_seed_arg;
}
