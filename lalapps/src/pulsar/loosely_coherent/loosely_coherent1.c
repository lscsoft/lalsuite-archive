/**
 * \file
 * \ingroup pulsarApps
 * \author Vladimir Dergachev
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fenv.h>
#include <unistd.h>
#include <time.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>

/* This likely has something to do with CUDA */
#define restrict

#include <lal/LALDatatypes.h>
#include <lal/AVFactories.h>
#include <lal/ComplexFFT.h>

FILE *LOG=NULL, *FILE_LOG=NULL, *DATA_LOG=NULL;

#include "global.h"
#include "dataset.h"
#include "cmdline.h"
#include "hookup.h"
#include "jobs.h"
#include "statistics.h"
#include "util.h"
#include "skymarks.h"
#include "bessel.h"
#include "context.h"

struct gengetopt_args_info args_info;

extern DATASET *datasets;
extern int d_free;


char *output_dir=NULL;

char *earth_ephemeris=NULL, *sun_ephemeris=NULL;

int first_bin=-1, useful_bins=-1, nbins=-1, side_cut=-1;

int nsegments=-1;

int no_am_response=0, fake_injection=0;

char *sky_marks=NULL;
int sky_marks_free=0;
int sky_marks_size=0;
SKYMARK *compiled_skymarks=NULL;

double spindown;
double spindown_start;

double average_det_velocity[3];
double det_vel_ra, det_vel_dec;
double orbital_axis[3];
double band_axis_norm;
double band_axis[3];

double resolution;

SKY_GRID *main_grid=NULL;

void *do_alloc(long a, long b)
{
void *r;
int i=0;
if(a<1)a=1;
if(b<1)b=1;
r=calloc(a,b);
while(r==NULL){
	fprintf(stderr,"Could not allocate %ld chunks of %ld bytes each (%ld bytes total), current memory usage %ld\n",a,b,a*b, MEMUSAGE);
	if(i>10)exit(-1);
	condor_safe_sleep(10);
	r=calloc(a,b);
	i++;
	}
return r;
}

char s[20000];

typedef struct {
	double dc;
	double sidereal_day[4];
	double sidereal_month[4];
	double sidereal_month2[4];
	double julian_year[4];
	} EPICYCLIC_PHASE_COEFFICIENTS;

typedef struct {
	double mean;
	double sd;
	double max_snr;
	int max_snr_index;
	} POWER_STATS;
	
	
void compute_power_stats(POWER_STATS *ps, double *power, int N)
{
double sum, sum_sq;
double max=power[0];
double min=power[0];
int idx=0, k;
double mean, mean2;
double sd, sd2;
int i, m;
int window=10;
double *p;

sum=0;
for(i=0;i<N;i++) {
	if(power[i]>max) {
		max=power[i];
		idx=i;
		} else
	if(power[i]<min) min=power[i];
	sum+=power[i];
	}

mean=sum/N;
sd=(mean-min)*0.25;

m=0;
if(max>mean+20*sd)while(1) {
	k=0;
	sum=0;
	sum_sq=0;
	m++;
	
	p=power;
	for(i=0;i<N;i++,p++) {
		if((abs(i-idx)<window)||(fabs(power[i]-mean)>20*sd)) {
			continue;
			}
		sum+=*p;
		sum_sq+=*p * *p;
		k++;
		}
	mean2=sum/k;
	sd2=sqrt((sum_sq-k*mean2*mean2)/(k-1.0));
	if(fabs(mean-mean2)<0.01*sd2 && fabs(sd-sd2)<0.01*sd2)break;
	mean=mean2;
	sd=sd2;
	}
//fprintf(stderr, "m=%d\n", m);
ps->mean=mean;
ps->sd=sd;
ps->max_snr=(max-mean)/sd;
ps->max_snr_index=idx;
}

struct {
	double ra;
	double dec;
	double dInv;
	LIGOTimeGPS tGPS1;
	LIGOTimeGPS tGPS2;
	char *detector;
	EmissionTime et1;
	EmissionTime et2;
	double range;
	}  ETC = {-10, -10, -10, -1, -1, -1, -1, NULL};

extern EphemerisData ephemeris;

void fast_get_emission_time(EmissionTime *emission_time, EarthState *earth_state, double ra, double dec, double dInv, char *detector, LIGOTimeGPS tGPS)
{
double dx, dt;
double step=90;

if(ETC.ra!=ra || ETC.dec!=dec || ETC.dInv!=dInv || strcmp(ETC.detector, detector)) {

	ETC.ra=ra;
	ETC.dec=dec;
	ETC.dInv=dInv;
	ETC.tGPS1=tGPS;
	ETC.tGPS2.gpsSeconds=-1;
	free(ETC.detector);
	ETC.detector=strdup(detector);
	
	if(earth_state==NULL) {
		EarthState earth_state1;
		LALStatus status={level:0, statusPtr:NULL};
		LALBarycenterEarth(&status, &(earth_state1), &ETC.tGPS1, &ephemeris);
		TESTSTATUS(&status);
		if(status.statusPtr)FREESTATUSPTR(&status);

		get_emission_time(&(ETC.et1), &earth_state1, ra, dec, dInv, detector, tGPS);
		} else {
		get_emission_time(&(ETC.et1), earth_state, ra, dec, dInv, detector, tGPS);
		}
	memcpy(emission_time, &(ETC.et1), sizeof(ETC.et1));
	return;	
	}
	
if(tGPS.gpsSeconds+step*0.5<ETC.tGPS1.gpsSeconds || tGPS.gpsSeconds>ETC.tGPS1.gpsSeconds+step*1.5) {
	ETC.ra=-10;
	fast_get_emission_time(emission_time, earth_state, ra, dec, dInv, detector, tGPS);
	return;
	}
	
if(ETC.tGPS2.gpsSeconds<0) {
	EarthState earth_state2;
	LALStatus status={level:0, statusPtr:NULL};

	ETC.tGPS2=ETC.tGPS1;
	ETC.tGPS2.gpsSeconds+=step;
	LALBarycenterEarth(&status, &(earth_state2), &ETC.tGPS2, &ephemeris);
	TESTSTATUS(&status);
	if(status.statusPtr)FREESTATUSPTR(&status);
	get_emission_time(&(ETC.et2), &earth_state2, ra, dec, dInv, detector, ETC.tGPS2);
	ETC.range=(double)(ETC.et2.te.gpsSeconds-ETC.et1.te.gpsSeconds)+1e-9*(double)(ETC.et2.te.gpsNanoSeconds-ETC.et1.te.gpsNanoSeconds);
	}
	
dx=((tGPS.gpsSeconds-ETC.tGPS1.gpsSeconds)+1e-9*(tGPS.gpsNanoSeconds-ETC.tGPS1.gpsNanoSeconds))/step;

dt=ETC.range*dx+ETC.et1.te.gpsNanoSeconds*1e-9;

memcpy(emission_time, &(ETC.et1), sizeof(ETC.et1));
emission_time->te.gpsSeconds+=floor(dt);
emission_time->te.gpsNanoSeconds=(dt-floor(dt))*1e9;

if(0){
/* Double check computation */
EmissionTime emission_time2;
double err;
get_emission_time(&emission_time2, earth_state, ra, dec, dInv, detector, tGPS);

err=(emission_time2.te.gpsSeconds-emission_time->te.gpsSeconds)+1e-9*(emission_time2.te.gpsNanoSeconds-emission_time->te.gpsNanoSeconds);

if(fabs(err)>1e-6) fprintf(stderr, "err=%g err2=%g dx=%g dt=%g range=%g gps=%g gps_delta2=%g\n", err, 
	(ETC.et2.te.gpsSeconds-emission_time->te.gpsSeconds)+1e-9*(ETC.et2.te.gpsNanoSeconds-emission_time->te.gpsNanoSeconds),
	dx, dt, ETC.range, tGPS.gpsSeconds+1e-9*tGPS.gpsNanoSeconds, (ETC.tGPS2.gpsSeconds-tGPS.gpsSeconds)+1e-9*(ETC.tGPS2.gpsNanoSeconds-tGPS.gpsNanoSeconds));
}

	
}


typedef struct {
	int free;
	int size;
	
	int *bin;
	COMPLEX16 *data;
	COMPLEX16 first9[9];
	
	double slope;
	} SPARSE_CONV;

void free_sparse_conv(SPARSE_CONV *sc)
{
free(sc->bin);
free(sc->data);
sc->bin=NULL;
sc->data=NULL;
free(sc);
}

SPARSE_CONV * compute_te_offset_structure(char *detector, double ra, double dec, double dInv, double gps_start, double step, int count)
{
COMPLEX16Vector *te_offsets;
COMPLEX16Vector *fft;
COMPLEX16FFTPlan *fft_plan=NULL;
LIGOTimeGPS tGPS;
EmissionTime emission_time;
double t;
double slope;
double a, max;
int i;

SPARSE_CONV *sc;
	
te_offsets=XLALCreateCOMPLEX16Vector(count);
if(!te_offsets) {
	fprintf(stderr, "**** ERROR: could not allocate te_offsets\n");
	exit(-1);
	}

fft=XLALCreateCOMPLEX16Vector(count);
if(!fft) {
	fprintf(stderr, "**** ERROR: could not allocate te_offsets\n");
	exit(-1);
	}
	
fft_plan=XLALCreateForwardCOMPLEX16FFTPlan(count, 1);
if(!fft_plan) {
	fprintf(stderr, "**** ERROR: could not create fft plan\n");
	exit(-1);
	}

for(i=0;i<count;i++) {
	t=gps_start+step*i;
	tGPS.gpsSeconds=floor(t);
	tGPS.gpsNanoSeconds=t-tGPS.gpsSeconds;

	fast_get_emission_time(&emission_time, NULL, ra, dec, dInv, detector, tGPS);
	
	te_offsets->data[i].re=(emission_time.te.gpsSeconds-t)+1e-9*emission_time.te.gpsNanoSeconds;
	te_offsets->data[i].im=0;
	}
	
//dump_doubles("te_offsets.dat", (double *)te_offsets->data, count, 2);

//slope=(te_offsets->data[count-1].re-te_offsets->data[0].re)/(count-1);
slope=(te_offsets->data[count-1].re-te_offsets->data[0].re)/count;
for(i=0;i<count;i++) {
	te_offsets->data[i].re-=slope*i;
	}
/* Normalize in units of time */
//fprintf(stderr, "slope/count=%g slope=%g\n", slope, slope/step);
slope=slope/step;

// Operator has 6 elements
// 0 1 -21814.6 -3210.39
// 1 2 -5393.65 -442.205
// 2 3 -2388.69 -140.633
// 3 -3 -2388.69 140.633
// 4 -2 -5393.65 442.205
// 5 -1 -21814.6 3210.39

XLALCOMPLEX16VectorFFT(fft, te_offsets, fft_plan);
//dump_doubles("te_offsets_fft.dat", (double *)fft->data, 2*count, 1);

sc=do_alloc(1, sizeof(*sc));
sc->free=0;
sc->size=20;
sc->bin=do_alloc(sc->size, sizeof(*sc->bin));
sc->data=do_alloc(sc->size, sizeof(*sc->data));

sc->slope=slope;

/* ignore constant term */
max=fft->data[1].re*fft->data[1].re+fft->data[1].im*fft->data[1].im;
for(i=2;i<count;i++) {
	a=fft->data[i].re*fft->data[i].re+fft->data[i].im*fft->data[i].im;
	if(a>max)max=a;
	}

for(i=1;i<count;i++) {
	a=fft->data[i].re*fft->data[i].re+fft->data[i].im*fft->data[i].im;
	if(a>(max*0.01)) {
		if(sc->free>=sc->size) {
			fprintf(stderr, "*** Aeiii ! Operator is too large\n");
			exit(-1);
			}
		sc->bin[sc->free]=2*i<count ? i : i-count;
		sc->data[sc->free].re=fft->data[i].re/count;
		sc->data[sc->free].im=fft->data[i].im/count;
		sc->free++;
		}
	}
	
// fprintf(stderr, "Operator has %d elements, fft[0]=%g\n", sc->free, fft->data[0].re/(count));
// for(i=0;i<sc->free;i++) {
// 	fprintf(stderr, "%d %d %g %g %g\n",i, sc->bin[i], sc->data[i].re, sc->data[i].im, sqrt(sc->data[i].re*sc->data[i].re+sc->data[i].im*sc->data[i].im));
// 	}

for(i=0;i<9;i++) {
	sc->first9[i].re=fft->data[i+1].re/count;
	sc->first9[i].im=fft->data[i+1].im/count;
	}

// for(i=0;i<3;i++) {
//  	fprintf(stderr, "%d %g %g %g\n",i, sc->first3[i].re, sc->first3[i].im, sqrt(sc->first3[i].re*sc->first3[i].re+sc->first3[i].im*sc->first3[i].im));
//  	}

XLALDestroyCOMPLEX16FFTPlan(fft_plan);
XLALDestroyCOMPLEX16Vector(te_offsets);
XLALDestroyCOMPLEX16Vector(fft);

return sc;
}

SPARSE_CONV * compute_spindown_offset_structure(char *detector, double ra, double dec, double dInv, double gps_start, double step, int count)
{
COMPLEX16Vector *te_offsets;
COMPLEX16Vector *fft;
COMPLEX16FFTPlan *fft_plan=NULL;
LIGOTimeGPS tGPS;
EmissionTime emission_time;
double t;
double slope;
double a, max;
int i;

SPARSE_CONV *sc;
	
te_offsets=XLALCreateCOMPLEX16Vector(count);
if(!te_offsets) {
	fprintf(stderr, "**** ERROR: could not allocate te_offsets\n");
	exit(-1);
	}

fft=XLALCreateCOMPLEX16Vector(count);
if(!fft) {
	fprintf(stderr, "**** ERROR: could not allocate te_offsets\n");
	exit(-1);
	}
	
fft_plan=XLALCreateForwardCOMPLEX16FFTPlan(count, 1);
if(!fft_plan) {
	fprintf(stderr, "**** ERROR: could not create fft plan\n");
	exit(-1);
	}

for(i=0;i<count;i++) {
	t=gps_start+step*i;
	tGPS.gpsSeconds=floor(t);
	tGPS.gpsNanoSeconds=t-tGPS.gpsSeconds;

	fast_get_emission_time(&emission_time, NULL, ra, dec, dInv, detector, tGPS);
	
	a=(emission_time.te.gpsSeconds-spindown_start)+1e-9*emission_time.te.gpsNanoSeconds;
	
	/* Count spindown in units of 1/timebase^2 */
	te_offsets->data[i].re=a*a*0.5/(step*step*count*count);
	te_offsets->data[i].im=0;
	}
	
//dump_doubles("te_offsets.dat", (double *)te_offsets->data, count, 2);

//slope=(te_offsets->data[count-1].re-te_offsets->data[0].re)/(count-1);
slope=(te_offsets->data[count-1].re-te_offsets->data[0].re)/count;
for(i=0;i<count;i++) {
	te_offsets->data[i].re-=slope*i;
	}
/* Normalize in units of time */
//fprintf(stderr, "slope/count=%g slope=%g\n", slope, slope/step);
slope=slope/step;

// Operator has 6 elements
// 0 1 -21814.6 -3210.39
// 1 2 -5393.65 -442.205
// 2 3 -2388.69 -140.633
// 3 -3 -2388.69 140.633
// 4 -2 -5393.65 442.205
// 5 -1 -21814.6 3210.39

XLALCOMPLEX16VectorFFT(fft, te_offsets, fft_plan);
//dump_doubles("te_offsets_fft.dat", (double *)fft->data, 2*count, 1);

sc=do_alloc(1, sizeof(*sc));
sc->free=0;
sc->size=20;
sc->bin=do_alloc(sc->size, sizeof(*sc->bin));
sc->data=do_alloc(sc->size, sizeof(*sc->data));

sc->slope=slope;

/* ignore constant term */
max=fft->data[1].re*fft->data[1].re+fft->data[1].im*fft->data[1].im;
for(i=2;i<count;i++) {
	a=fft->data[i].re*fft->data[i].re+fft->data[i].im*fft->data[i].im;
	if(a>max)max=a;
	}

for(i=1;i<count;i++) {
	a=fft->data[i].re*fft->data[i].re+fft->data[i].im*fft->data[i].im;
	if(a>(max*0.01)) {
		if(sc->free>=sc->size) {
			fprintf(stderr, "*** Aeiii ! Operator is too large\n");
			exit(-1);
			}
		sc->bin[sc->free]=2*i<count ? i : i-count;
		sc->data[sc->free].re=fft->data[i].re/count;
		sc->data[sc->free].im=fft->data[i].im/count;
		sc->free++;
		}
	}
	
// fprintf(stderr, "Operator has %d elements, fft[0]=%g\n", sc->free, fft->data[0].re/(count));
// for(i=0;i<sc->free;i++) {
// 	fprintf(stderr, "%d %d %g %g %g\n",i, sc->bin[i], sc->data[i].re, sc->data[i].im, sqrt(sc->data[i].re*sc->data[i].re+sc->data[i].im*sc->data[i].im));
// 	}

for(i=0;i<9;i++) {
	sc->first9[i].re=fft->data[i+1].re/count;
	sc->first9[i].im=fft->data[i+1].im/count;
	}

// for(i=0;i<3;i++) {
//  	fprintf(stderr, "%d %g %g %g\n",i, sc->first3[i].re, sc->first3[i].im, sqrt(sc->first3[i].re*sc->first3[i].re+sc->first3[i].im*sc->first3[i].im));
//  	}

XLALDestroyCOMPLEX16FFTPlan(fft_plan);
XLALDestroyCOMPLEX16Vector(te_offsets);
XLALDestroyCOMPLEX16Vector(fft);

return sc;
}

void compute_sky_basis(double ra, double dec, double step, double *ra_out, double *dec_out)
{
double v1[3], vp1[3], vp2[3];
double a,b;
int i;

v1[0]=cos(ra)*cos(dec);
v1[1]=sin(ra)*cos(dec);
v1[2]=sin(dec);

vp1[0]=sin(ra);
vp1[1]=-cos(ra);
vp1[2]=0;

vp2[0]=cos(ra)*sin(dec);
vp2[1]=sin(ra)*sin(dec);
vp2[2]=-cos(dec);

for(i=0;i<3;i++) {
	vp1[i]=-vp1[i]*step-vp2[i]*step+v1[i];
	vp2[i]=-vp2[i]*step+v1[i];
	}

a=1.0/sqrt(vp1[0]*vp1[0]+vp1[1]*vp1[1]+vp1[2]*vp1[2]);
b=1.0/sqrt(vp2[0]*vp2[0]+vp2[1]*vp2[1]+vp2[2]*vp2[2]);
for(i=0;i<3;i++) {
	vp1[i]*=a;
	vp2[i]*=b;
	}
ra_out[0]=atan2(vp1[1], vp1[0]);
ra_out[1]=atan2(vp2[1], vp2[0]);
dec_out[0]=asin(vp1[2]);
dec_out[1]=asin(vp2[2]);
}

SPARSE_CONV * compute_sky_offset_structure(char *detector, double ra, double dec, double ra2, double dec2, double dInv, double gps_start, double step, int count)
{
COMPLEX16Vector *te_offsets;
COMPLEX16Vector *fft;
COMPLEX16FFTPlan *fft_plan=NULL;
LIGOTimeGPS tGPS;
EmissionTime emission_time;
double t;
double slope;
double a, max;
int i;

SPARSE_CONV *sc;
	
te_offsets=XLALCreateCOMPLEX16Vector(count);
if(!te_offsets) {
	fprintf(stderr, "**** ERROR: could not allocate te_offsets\n");
	exit(-1);
	}

fft=XLALCreateCOMPLEX16Vector(count);
if(!fft) {
	fprintf(stderr, "**** ERROR: could not allocate te_offsets\n");
	exit(-1);
	}
	
fft_plan=XLALCreateForwardCOMPLEX16FFTPlan(count, 1);
if(!fft_plan) {
	fprintf(stderr, "**** ERROR: could not create fft plan\n");
	exit(-1);
	}

/* scan twice to make most use of fast_get_emission_time */
for(i=0;i<count;i++) {
	t=gps_start+step*i;
	tGPS.gpsSeconds=floor(t);
	tGPS.gpsNanoSeconds=t-tGPS.gpsSeconds;

	fast_get_emission_time(&emission_time, NULL, ra, dec, dInv, detector, tGPS);
	
	a=(emission_time.te.gpsSeconds-spindown_start)+1e-9*emission_time.te.gpsNanoSeconds;
	
	/* Count offset in units of 1/f */
	te_offsets->data[i].re=a*args_info.focus_f0_arg;
	te_offsets->data[i].im=0;
	}

for(i=0;i<count;i++) {
	t=gps_start+step*i;
	tGPS.gpsSeconds=floor(t);
	tGPS.gpsNanoSeconds=t-tGPS.gpsSeconds;

	fast_get_emission_time(&emission_time, NULL, ra2, dec2, dInv, detector, tGPS);
	
	a=(emission_time.te.gpsSeconds-spindown_start)+1e-9*emission_time.te.gpsNanoSeconds;
	
	/* Count offset in units of 1/f */
	te_offsets->data[i].re-=a*args_info.focus_f0_arg;
	}

//dump_doubles("te_offsets.dat", (double *)te_offsets->data, count, 2);

//slope=(te_offsets->data[count-1].re-te_offsets->data[0].re)/(count-1);
slope=(te_offsets->data[count-1].re-te_offsets->data[0].re)/count;
for(i=0;i<count;i++) {
	te_offsets->data[i].re-=slope*i;
	}
/* Normalize in units of time */
//fprintf(stderr, "slope/count=%g slope=%g\n", slope, slope/step);
slope=slope/step;

// Operator has 6 elements
// 0 1 -21814.6 -3210.39
// 1 2 -5393.65 -442.205
// 2 3 -2388.69 -140.633
// 3 -3 -2388.69 140.633
// 4 -2 -5393.65 442.205
// 5 -1 -21814.6 3210.39

XLALCOMPLEX16VectorFFT(fft, te_offsets, fft_plan);
//dump_doubles("te_offsets_fft.dat", (double *)fft->data, 2*count, 1);

sc=do_alloc(1, sizeof(*sc));
sc->free=0;
sc->size=20;
sc->bin=do_alloc(sc->size, sizeof(*sc->bin));
sc->data=do_alloc(sc->size, sizeof(*sc->data));

sc->slope=slope;

/* ignore constant term */
max=fft->data[1].re*fft->data[1].re+fft->data[1].im*fft->data[1].im;
for(i=2;i<count;i++) {
	a=fft->data[i].re*fft->data[i].re+fft->data[i].im*fft->data[i].im;
	if(a>max)max=a;
	}

for(i=1;i<count;i++) {
	a=fft->data[i].re*fft->data[i].re+fft->data[i].im*fft->data[i].im;
	if(a>(max*0.001)) {
		if(sc->free>=sc->size) {
			fprintf(stderr, "*** Aeiii ! Operator is too large\n");
			exit(-1);
			}
		sc->bin[sc->free]=2*i<count ? i : i-count;
		sc->data[sc->free].re=fft->data[i].re/count;
		sc->data[sc->free].im=fft->data[i].im/count;
		sc->free++;
		}
	}
	
fprintf(stderr, "Operator has %d elements, fft[0]=%g\n", sc->free, fft->data[0].re/(count));
for(i=0;i<sc->free;i++) {
	fprintf(stderr, "%d %d %g %g %g\n",i, sc->bin[i], sc->data[i].re, sc->data[i].im, sqrt(sc->data[i].re*sc->data[i].re+sc->data[i].im*sc->data[i].im));
	}

for(i=0;i<9;i++) {
	sc->first9[i].re=fft->data[i+1].re/count;
	sc->first9[i].im=fft->data[i+1].im/count;
	}

// for(i=0;i<3;i++) {
//  	fprintf(stderr, "%d %g %g %g\n",i, sc->first3[i].re, sc->first3[i].im, sqrt(sc->first3[i].re*sc->first3[i].re+sc->first3[i].im*sc->first3[i].im));
//  	}

XLALDestroyCOMPLEX16FFTPlan(fft_plan);
XLALDestroyCOMPLEX16Vector(te_offsets);
XLALDestroyCOMPLEX16Vector(fft);

return sc;
}


/* TODO: optimize */


void compute_fft_stats(COMPLEX16Vector *te_fft, double *power, int nsamples, POWER_STATS *ps)
{
int i, i_left, i_right;
double a, b, x, y;
double spread=0.05*0;

for(i=0;i<nsamples;i++) {
	i_left=i-1;
	i_right=i+1;
	if(i_left<0)i_left+=nsamples;
	if(i_right>=nsamples)i_right-=nsamples;
	
	a=te_fft->data[i].re;
	b=te_fft->data[i].im;
	
	x=a*(te_fft->data[i_left].re-te_fft->data[i_right].re)-b*(te_fft->data[i_left].im-te_fft->data[i_right].im);
	y=b*(te_fft->data[i_left].re+te_fft->data[i_right].re)+a*(te_fft->data[i_left].im+te_fft->data[i_right].im);
	
	power[i]=a*a+b*b+spread*sqrt(x*x+y*y);
	}

// sum=0.0;
// for(i=0;i<nsamples;i++) {
// 	a=0.0;
// 	b=0.0;
// 	c=0.0;
// 	for(j=i-half_window; j<=i+half_window;j++) {
// 		k=j;
// 		if(k<0)k+=nsamples;
// 		if(k>=nsamples)k-=nsamples;
// 		x=fft->data[k].re;
// 		y=fft->data[k].im;
// 		a+=x*x+y*y;
// 
// 		k=j+wing_step;
// 		if(k<0)k+=nsamples;
// 		if(k>=nsamples)k-=nsamples;
// 		x1=fft->data[k].re;
// 		y1=fft->data[k].im;
// 		c+=(x1*x1+y1*y1);
// 		b+=(x1*x+y1*y);
// 
// 		k=j-wing_step;
// 		if(k<0)k+=nsamples;
// 		if(k>=nsamples)k-=nsamples;
// 		x1=fft->data[k].re;
// 		y1=fft->data[k].im;
// 		c+=(x1*x1+y1*y1);
// 		b+=(x1*x+y1*y);
// 		}
// 	power[i]=a;
// 	cross[i]=b;
// 	power_aux[i]=c;
// 
// 	sum+=a+c;
// 	}
// 	
// sum/=nsamples;


// for(i=0;i<nsamples;i++) {
// 	if(power[i]+power_aux[i]>3*sum)
// 	fprintf(DATA_LOG, "fft: %d %d %.12f %.12f %.12f %.12f %.12g %.12g %.12g %.12g %.12g\n", 
// 		point, fstep, (i*2>nsamples ? i-nsamples : i)+fstep*1.0/nfsteps, ((i*2>nsamples ? i-nsamples : i)*1.0 +fstep*1.0/nfsteps)/(nsamples*args_info.coherence_length_arg), ra, dec, power[i], cross[i], power_aux[i], fft->data[i].re, fft->data[i].im);
// 	}

compute_power_stats(ps, power, nsamples);
}

void compute_fft_statsf(COMPLEX16Vector *te_fft, double *power, int nsamples, POWER_STATS *ps)
{
int i, i_left, i_right;
float a, b, x, y;
float spread=0.05;

for(i=0;i<nsamples;i++) {
	i_left=i-1;
	i_right=i+1;
	if(i_left<0)i_left+=nsamples;
	if(i_right>=nsamples)i_right-=nsamples;
	
	a=te_fft->data[i].re;
	b=te_fft->data[i].re;
	
	x=a*(te_fft->data[i_left].re-te_fft->data[i_right].re)-b*(te_fft->data[i_left].im-te_fft->data[i_right].im);
	y=b*(te_fft->data[i_left].re+te_fft->data[i_right].re)+a*(te_fft->data[i_left].im+te_fft->data[i_right].im);
	
	power[i]=a*a+b*b+spread*sqrtf(x*x+y*y);
	}

// sum=0.0;
// for(i=0;i<nsamples;i++) {
// 	a=0.0;
// 	b=0.0;
// 	c=0.0;
// 	for(j=i-half_window; j<=i+half_window;j++) {
// 		k=j;
// 		if(k<0)k+=nsamples;
// 		if(k>=nsamples)k-=nsamples;
// 		x=fft->data[k].re;
// 		y=fft->data[k].im;
// 		a+=x*x+y*y;
// 
// 		k=j+wing_step;
// 		if(k<0)k+=nsamples;
// 		if(k>=nsamples)k-=nsamples;
// 		x1=fft->data[k].re;
// 		y1=fft->data[k].im;
// 		c+=(x1*x1+y1*y1);
// 		b+=(x1*x+y1*y);
// 
// 		k=j-wing_step;
// 		if(k<0)k+=nsamples;
// 		if(k>=nsamples)k-=nsamples;
// 		x1=fft->data[k].re;
// 		y1=fft->data[k].im;
// 		c+=(x1*x1+y1*y1);
// 		b+=(x1*x+y1*y);
// 		}
// 	power[i]=a;
// 	cross[i]=b;
// 	power_aux[i]=c;
// 
// 	sum+=a+c;
// 	}
// 	
// sum/=nsamples;


// for(i=0;i<nsamples;i++) {
// 	if(power[i]+power_aux[i]>3*sum)
// 	fprintf(DATA_LOG, "fft: %d %d %.12f %.12f %.12f %.12f %.12g %.12g %.12g %.12g %.12g\n", 
// 		point, fstep, (i*2>nsamples ? i-nsamples : i)+fstep*1.0/nfsteps, ((i*2>nsamples ? i-nsamples : i)*1.0 +fstep*1.0/nfsteps)/(nsamples*args_info.coherence_length_arg), ra, dec, power[i], cross[i], power_aux[i], fft->data[i].re, fft->data[i].im);
// 	}

compute_power_stats(ps, power, nsamples);
}

void fast_shift_fft_R(COMPLEX16Vector *fft_out, COMPLEX16Vector *fft_in, int nsamples, double v0, double v1)
{
int i, i_left, i_right;
double a, b, x, y;
for(i=0;i<nsamples;i++) {
	i_left=i-1;
	i_right=i+1;
	if(i_left<0)i_left+=nsamples;
	if(i_right>=nsamples)i_right-=nsamples;
	
	a=fft_in->data[i_left].re;
	b=fft_in->data[i_left].im;
	
	x=-v0*a;
	y=-v0*b;

	a=fft_in->data[i].re;
	b=fft_in->data[i].im;
	
	x+=v1*a;
	y+=v1*b;

	a=fft_in->data[i_right].re;
	b=fft_in->data[i_right].im;
	
	x+=v0*a;
	y+=v0*b;
	
	fft_out->data[i].re=x;
	fft_out->data[i].im=y;
	}	
}

void fast_shift_fft_I(COMPLEX16Vector *fft_out, COMPLEX16Vector *fft_in, int nsamples, double v0, double v1)
{
int i, i_left, i_right;
double a, b, x, y;
for(i=0;i<nsamples;i++) {
	i_left=i-1;
	i_right=i+1;
	if(i_left<0)i_left+=nsamples;
	if(i_right>=nsamples)i_right-=nsamples;
	
	a=fft_in->data[i_left].re;
	b=fft_in->data[i_left].im;
	
	x=-v0*b;
	y=v0*a;

	a=fft_in->data[i].re;
	b=fft_in->data[i].im;
	
	x+=v1*a;
	y+=v1*b;

	a=fft_in->data[i_right].re;
	b=fft_in->data[i_right].im;
	
	x+=-v0*b;
	y+=v0*a;
	
	fft_out->data[i].re=x;
	fft_out->data[i].im=y;
	}	
}

void scan_fft_stats(LOOSE_CONTEXT *ctx, double *power, int nsamples, COMPLEX16 *v1, COMPLEX16 *v2, POWER_STATS *ps)
{
COMPLEX16Vector *fft2, *fft3, *fft4, *fft5, *fft_tmp;	
#define N_SCAN_FFT_FILTER 9
COMPLEX16 filter1[N_SCAN_FFT_FILTER], filter2[N_SCAN_FFT_FILTER];
POWER_STATS ps2;
double ms;
int i, j;
int nscan=2;

fft2=XLALCreateCOMPLEX16Vector(nsamples);
fft3=XLALCreateCOMPLEX16Vector(nsamples);
fft4=XLALCreateCOMPLEX16Vector(nsamples);
fft5=XLALCreateCOMPLEX16Vector(nsamples);

compute_fft_stats(ctx->plus_te_fft, power, nsamples, ps);
ms=ps->max_snr;

make_bessel_filter(filter1, N_SCAN_FFT_FILTER, v1, 3, 2*M_PI/(nscan));
make_bessel_filter(filter2, N_SCAN_FFT_FILTER, v2, 3, 2*M_PI/(nscan));

if(0) {
	for(i=0;i<3;i++)
		fprintf(stderr, "%d %f %f %f %f\n", i, v1[i].re, v1[i].im, v2[i].re, v2[i].im);
	fprintf(stderr, "\n");

	for(i=0;i<N_SCAN_FFT_FILTER;i++)
		fprintf(stderr, "%d %f %f %f %f\n", i, filter1[i].re, filter1[i].im, filter2[i].re, filter2[i].im);
	fprintf(stderr, "\n");

	for(i=0;i<ctx->nsamples;i++)
		fprintf(DATA_LOG, "te_fft: %d %g %g\n", i, ctx->plus_te_fft->data[i].re, ctx->plus_te_fft->data[i].im);

	fflush(DATA_LOG);
	exit(-1);
	}

for(i=0;i<=nscan;i++) {
	if(i==0)memcpy(fft2->data, ctx->plus_te_fft->data, nsamples*sizeof(COMPLEX16));
		else {
		shift_fft(fft3, fft2, filter1, N_SCAN_FFT_FILTER);
		fft_tmp=fft2;
		fft2=fft3;
		fft3=fft_tmp;
		}
	for(j=0;j<=nscan;j++) {

		if(j==0)memcpy(fft4->data, fft2->data, nsamples*sizeof(COMPLEX16));
			else {
			shift_fft(fft5, fft4, filter2, N_SCAN_FFT_FILTER);
			fft_tmp=fft4;
			fft4=fft5;
			fft5=fft_tmp;
			}

		if(i==0 && j==0)continue;

		compute_fft_stats(fft4, power, nsamples, &ps2);
/*		fprintf(stderr, "+ %.2f ", (ps2.max_snr-ps->max_snr)/ps->max_snr);*/
		if(ps2.max_snr>ps->max_snr) {
			fprintf(stderr, "* %f", ps2.max_snr-ms);
			memcpy(ps, &ps2, sizeof(*ps));
			}
		}
/*	fprintf(stderr, "\n");*/
	}

make_bessel_filter(filter1, N_SCAN_FFT_FILTER, v1, 3, -2*M_PI);
make_bessel_filter(filter2, N_SCAN_FFT_FILTER, v2, 3, 2*M_PI);

for(i=0;i<=nscan;i++) {
	if(i==0) {
		memcpy(fft2->data, ctx->plus_te_fft->data, nsamples*sizeof(COMPLEX16));
		continue;
		}
		else {
		shift_fft(fft3, fft2, filter1, N_SCAN_FFT_FILTER);
		fft_tmp=fft2;
		fft2=fft3;
		fft3=fft_tmp;
		}
	for(j=0;j<=nscan;j++) {

		if(j==0)memcpy(fft4->data, fft2->data, nsamples*sizeof(COMPLEX16));
			else {
			shift_fft(fft5, fft4, filter2, N_SCAN_FFT_FILTER);
			fft_tmp=fft4;
			fft4=fft5;
			fft5=fft_tmp;
			}

		if(i==0 && j==0)continue;

		compute_fft_stats(fft4, power, nsamples, &ps2);
		if(ps2.max_snr>ps->max_snr) {
			fprintf(stderr, "* %f", ps2.max_snr-ms);
			memcpy(ps, &ps2, sizeof(*ps));
			}
		}
	}

make_bessel_filter(filter1, N_SCAN_FFT_FILTER, v1, 3, 2*M_PI);
make_bessel_filter(filter2, N_SCAN_FFT_FILTER, v2, 3, -2*M_PI);

for(i=0;i<=nscan;i++) {
	if(i==0) {
		memcpy(fft2->data, ctx->plus_te_fft->data, nsamples*sizeof(COMPLEX16));
		}
		else {
		shift_fft(fft3, fft2, filter1, N_SCAN_FFT_FILTER);
		fft_tmp=fft2;
		fft2=fft3;
		fft3=fft_tmp;
		}
	for(j=0;j<=nscan;j++) {

		if(j==0) {
			memcpy(fft4->data, fft2->data, nsamples*sizeof(COMPLEX16));
			continue;
			} else {
			shift_fft(fft5, fft4, filter2, N_SCAN_FFT_FILTER);
			fft_tmp=fft4;
			fft4=fft5;
			fft5=fft_tmp;
			}

		if(i==0 && j==0)continue;

		compute_fft_stats(fft4, power, nsamples, &ps2);
		if(ps2.max_snr>ps->max_snr) {
			fprintf(stderr, "* %f", ps2.max_snr-ms);
			memcpy(ps, &ps2, sizeof(*ps));
			}
		}
	}

make_bessel_filter(filter1, N_SCAN_FFT_FILTER, v1, 3, -2*M_PI);
make_bessel_filter(filter2, N_SCAN_FFT_FILTER, v2, 3, -2*M_PI);

for(i=0;i<=nscan;i++) {
	if(i==0) {
		memcpy(fft2->data, ctx->plus_te_fft->data, nsamples*sizeof(COMPLEX16));
		continue;
		}
		else {
		shift_fft(fft3, fft2, filter1, N_SCAN_FFT_FILTER);
		fft_tmp=fft2;
		fft2=fft3;
		fft3=fft_tmp;
		}
	for(j=0;j<=nscan;j++) {

		if(j==0) {
			memcpy(fft4->data, fft2->data, nsamples*sizeof(COMPLEX16));
			continue;
			} else {
			shift_fft(fft5, fft4, filter2, N_SCAN_FFT_FILTER);
			fft_tmp=fft4;
			fft4=fft5;
			fft5=fft_tmp;
			}

		if(i==0 && j==0)continue;

		compute_fft_stats(fft4, power, nsamples, &ps2);
		if(ps2.max_snr>ps->max_snr) {
			fprintf(stderr, "* %f", ps2.max_snr-ms);
			memcpy(ps, &ps2, sizeof(*ps));
			}
		}
	}

XLALDestroyCOMPLEX16Vector(fft2);
XLALDestroyCOMPLEX16Vector(fft3);
XLALDestroyCOMPLEX16Vector(fft4);
XLALDestroyCOMPLEX16Vector(fft5);
}

int main(int argc, char *argv[]) 
{
time_t start_time, input_done_time, end_time;
struct rlimit rl;
int i, j;

time(&start_time);

test_bessel();

fedisableexcept(FE_ALL_EXCEPT);
fprintf(stderr, "Initial memory: %g MB\n", (MEMUSAGE*10.0/(1024.0*1024.0))/10.0);
if(getrlimit(RLIMIT_AS, &rl)<0) {
	perror("Could not obtain virtual memory limit:");
	} else {
	fprintf(stderr, "Virtual memory limits soft=%ld hard=%ld\n", rl.rlim_cur, rl.rlim_max);
	}

if(getrlimit(RLIMIT_DATA, &rl)<0) {
	perror("Could not obtain data segment limit:");
	} else {
	fprintf(stderr, "Data segment limits soft=%ld hard=%ld\n", rl.rlim_cur, rl.rlim_max);
	}

if(cmdline_parser(argc, argv, &args_info))exit(-1);

for(i=0;i<args_info.config_given;i++)
	if(cmdline_parser_configfile(args_info.config_arg[i], &args_info, 0, 0, 0))exit(-1);

if(!args_info.dataset_given){
	fprintf(stderr,"** You must specify dataset description file (--dataset)\n");
	exit(-1);
	}
if(!args_info.ephemeris_path_given &&
	!(args_info.earth_ephemeris_given && args_info.sun_ephemeris_given)){
	fprintf(stderr,"** You must specify path to ephemeris files (--detreponse-path or --XXX-ephemeris)\n");
	exit(-1);
	}
if(!args_info.first_bin_given){
	fprintf(stderr,"** You must specify first bin to analyze (--first-bin)\n");
	exit(-1);
	}

gsl_rng_env_setup();
gsl_set_error_handler_off();

/* create output directories if not present */
if(args_info.output_given){
	mkdir(args_info.output_arg, 0777);
	output_dir=do_alloc(strlen(args_info.output_arg)+30, sizeof(*output_dir));
	sprintf(output_dir, "%s/%d-%f/", args_info.output_arg, args_info.first_bin_arg,args_info.first_bin_arg/args_info.coherence_length_arg);
	} else {
	output_dir=do_alloc(30, sizeof(*output_dir));
	sprintf(output_dir, "%d-%f/", args_info.first_bin_arg,args_info.first_bin_arg/args_info.coherence_length_arg);
	}
mkdir(output_dir, 0777);


snprintf(s,20000,"%s/powerflux.log", output_dir);
LOG=fopen(s,"w");

while(LOG==NULL) {
	i=errno;
	fprintf(stderr, "Could not open log file \"%s\": %s\n", s, strerror(errno));
	condor_safe_sleep(60);
	LOG=fopen(s, "w");
	}

snprintf(s,20000,"%s/file.log", output_dir);
FILE_LOG=fopen(s,"w");

snprintf(s,20000,"%s/data.log", output_dir);
DATA_LOG=fopen(s,"w");
setbuffer(DATA_LOG, do_alloc(1024*1024*32, 1), 1024*1024*32);

if(args_info.label_given){
	fprintf(LOG, "label: \"%s\"\n", args_info.label_arg);
	fprintf(stderr, "label: \"%s\"\n", args_info.label_arg);
	}
	
if(gethostname(s, 19999)>=0){
	fprintf(stderr, "Running on %s\n", s);
	fprintf(LOG, "node: %s\n", s);
	} else {
	fprintf(stderr, "Could not obtain hostname\n");
	fprintf(LOG, "node: unknown\n");
	}

for(i=0;i<args_info.config_given;i++) {
	fprintf(LOG, "using config file: %s\n", args_info.config_arg[i]);
	}

init_threads(args_info.num_threads_arg);
init_jobs();
init_hookup();
init_statistics();
tabulate_hann_filter7();

if(args_info.earth_ephemeris_given){
	earth_ephemeris=args_info.earth_ephemeris_arg;
	} else {
	earth_ephemeris=do_alloc(strlen(args_info.ephemeris_path_arg)+20,1);
	sprintf(earth_ephemeris,"%s/earth00-04.dat",args_info.ephemeris_path_arg);
	}
	
if(args_info.sun_ephemeris_given){
	sun_ephemeris=args_info.sun_ephemeris_arg;
	} else {
	sun_ephemeris=do_alloc(strlen(args_info.ephemeris_path_arg)+20,1);
	sprintf(sun_ephemeris,"%s/sun00-04.dat",args_info.ephemeris_path_arg);
	}
	
init_ephemeris();



no_am_response=args_info.no_am_response_arg;

fprintf(LOG,"loose1 : %s\n",VERSION);
if(no_am_response){
	fprintf(LOG,"no_am_response : true\n");
	fprintf(stderr,"NO_AM_RESPONSE flag passed\n");
	}
fprintf(LOG,"algorithm: %s\n", args_info.algorithm_arg);

if(args_info.dataset_given)fprintf(LOG,"dataset: %s\n",args_info.dataset_arg);

fflush(LOG);


/* do we need to inject fake signal ? */
if(args_info.fake_freq_given) {
	fake_injection=1;

   	fprintf(LOG,"fake signal injection: yes\n");
	fprintf(LOG,"fake psi : %f\n", args_info.fake_psi_arg);
	fprintf(LOG,"fake phi : %f\n", args_info.fake_phi_arg);
	fprintf(LOG,"fake iota : %f\n", args_info.fake_iota_arg);
	fprintf(LOG,"fake ra : %f\n", args_info.fake_ra_arg);
	fprintf(LOG,"fake dec: %f\n", args_info.fake_dec_arg);
	fprintf(LOG,"fake spindown: %g\n", args_info.fake_spindown_arg);
	fprintf(LOG,"fake strain: %g\n", args_info.fake_strain_arg);
	fprintf(LOG,"fake frequency: %f\n", args_info.fake_freq_arg);
	fprintf(LOG,"fake reference time: %f\n", args_info.fake_ref_time_arg);
	

	fprintf(LOG,"fake dInv: %g\n", args_info.fake_dInv_arg);
	fprintf(LOG,"fake frequency modulation depth: %f\n", args_info.fake_freq_modulation_depth_arg);
	fprintf(LOG,"fake frequency modulation frequency: %f\n", args_info.fake_freq_modulation_freq_arg);
	fprintf(LOG,"fake frequency modulation phase: %f\n", args_info.fake_freq_modulation_phase_arg);
	fprintf(LOG,"fake phase modulation depth: %f\n", args_info.fake_phase_modulation_depth_arg);
	fprintf(LOG,"fake phase modulation frequency: %f\n", args_info.fake_phase_modulation_freq_arg);
	fprintf(LOG,"fake phase modulation phase: %f\n", args_info.fake_phase_modulation_phase_arg);
	fprintf(LOG,"fake injection window: %d\n", args_info.fake_injection_window_arg);
   	} else {
   	fprintf(LOG,"fake signal injection: none\n");
	}

side_cut=2;
first_bin=args_info.first_bin_arg-side_cut;
useful_bins=args_info.nbins_arg;
nbins=args_info.nbins_arg+2*side_cut;

fprintf(LOG,"firstbin  : %d\n",first_bin);
fprintf(LOG,"band start: %g Hz\n",first_bin/args_info.coherence_length_arg);
fprintf(LOG,"nbins     : %d\n",nbins);
fprintf(LOG,"side_cut  : %d\n",side_cut);
fprintf(LOG,"useful bins : %d\n",useful_bins);
fprintf(LOG,"useful band start: %g Hz\n",(first_bin+side_cut)/args_info.coherence_length_arg);

init_datasets();
test_datasets();

if(args_info.dataset_given) {
	fprintf(stderr, "Loading data from dataset %s\n", args_info.dataset_arg);
	load_dataset_from_file(args_info.dataset_arg);
	nsegments=total_segments();
	if(nsegments==vetoed_segments()) {
		fprintf(LOG, "All SFTs vetoed, aborting!\n");
		fprintf(stderr, "All SFTs vetoed, aborting!\n");
		exit(-1);
		}
	}

/* This diagnostics should be moved into dataset.c when tested */
for(i=0;i<d_free;i++) {
	fprintf(LOG, "FMedians: \"%s\" \"%s\"", args_info.label_arg, datasets[i].name);
	for(j=0;j<datasets[i].nbins;j++)fprintf(LOG, " %g", datasets[i].TMedians[j]);
	fprintf(LOG, "\n");
	}

for(i=0;i<d_free;i++) {
	fprintf(LOG, "new_weighted_mean: \"%s\" \"%s\"", args_info.label_arg, datasets[i].name);
	for(j=0;j<datasets[i].nbins;j++)fprintf(LOG, " %g", datasets[i].new_weighted_mean[j]);
	fprintf(LOG, "\n");
	}

TODO("compute sky resolution automatically from loaded data and mode")
resolution=args_info.skymap_resolution_arg;
resolution*=args_info.skymap_resolution_ratio_arg;

fprintf(LOG,"resolution : %f\n", resolution);
TODO("warn (and error out in extreme cases) if we get too many templates")
if(!strcasecmp("sin_theta", args_info.sky_grid_arg)){
	main_grid=make_sin_theta_grid(resolution);
	} else
if(!strcasecmp("plain_rectangular", args_info.sky_grid_arg)){
	main_grid=make_rect_grid(ceil(2.0*M_PI/(resolution)), ceil(M_PI/(resolution)));
	} else
if(!strcasecmp("targeted_rectangular", args_info.sky_grid_arg)){
	if(!(args_info.focus_ra_given && args_info.focus_dec_given && args_info.focus_radius_given)) {
		fprintf(stderr, "*** ERROR: focus* options are required for targeted rectangular grid\n"); 
		}
	main_grid=make_targeted_rect_grid(args_info.focus_ra_arg, args_info.focus_dec_arg, args_info.focus_radius_arg, ceil(2*args_info.focus_radius_arg/resolution)+2);
	} else
if(!strcasecmp("arcsin", args_info.sky_grid_arg)){
	main_grid=make_arcsin_grid(ceil(2.0*M_PI/(resolution)), ceil(M_PI/(resolution)));
	} else {
	fprintf(stderr,"*** Unknown sky grid type: \"%s\"\n", args_info.sky_grid_arg);
	exit(-1);
	}

fprintf(stderr, "Main sky grid has %d points.\n", main_grid->npoints);

if(args_info.sky_marks_file_given) {
	FILE *f;
	int a;
	f=fopen(args_info.sky_marks_file_arg, "r");
	if(f==NULL) {
		perror("Could not open skymarks file:");
		exit(-1);
		}
	sky_marks_size=10000;
	sky_marks=do_alloc(sky_marks_size, sizeof(char));
	sky_marks_free=0;

	while(1) {
		a=fread(&(sky_marks[sky_marks_free]), 1, 10000, f);
		if(a>0)sky_marks_free+=a;
		if(a<10000)break;

		if(sky_marks_free+10000>=sky_marks_size) {
			char *ptr;
			sky_marks_size*=2;
			ptr=do_alloc(sky_marks_size, sizeof(char));
			memcpy(ptr, sky_marks, sky_marks_free);
			free(sky_marks);
			sky_marks=ptr;
			}
		}
	fclose(f);
	fprintf(LOG, "sky marks file: \"%s\"\n", args_info.sky_marks_file_arg);
	}
	
if(args_info.dump_data_given) {
	dump_datasets(args_info.dump_data_arg);
	}

if(nsegments==0){
	fprintf(stderr,"ERROR: no input data found !\n");
	return -1;
	}
post_init_datasets();

if(args_info.dump_sftv2_given) {
	sftv2_dump_datasets(args_info.dump_sftv2_arg);
	}

fprintf(LOG,"nsegments : %d\n", nsegments);

fprintf(LOG,"first gps : %lld\n", min_gps());
fprintf(LOG,"last gps  : %lld\n", max_gps());

if(!args_info.spindown_start_time_given){
	spindown_start=min_gps();
	} else {
	spindown_start=args_info.spindown_start_time_arg;
	}
fprintf(LOG, "spindown start time: %f\n", spindown_start);

fflush(LOG);

fprintf(stderr,"Computing detector speed\n");
datasets_average_detector_speed(average_det_velocity);
fprintf(LOG,"average detector velocity: %g %g %g\n", 
	average_det_velocity[0],
	average_det_velocity[1],
	average_det_velocity[2]);
det_vel_dec=atan2f(average_det_velocity[2], 
	sqrt(average_det_velocity[0]*average_det_velocity[0]+average_det_velocity[1]*average_det_velocity[1]));
det_vel_ra=atan2f(average_det_velocity[1], average_det_velocity[0]);
if(det_vel_ra<0)det_vel_ra+=2.0*M_PI;
fprintf(LOG,"average detector velocity RA (degrees) : %f\n", det_vel_ra*180.0/M_PI);
fprintf(LOG,"average detector velocity DEC (degrees): %f\n", det_vel_dec*180.0/M_PI);
fprintf(stderr,"average detector velocity RA (degrees) : %f\n", det_vel_ra*180.0/M_PI);
fprintf(stderr,"average detector velocity DEC (degrees): %f\n", det_vel_dec*180.0/M_PI);

orbital_axis[0]=0.0;
orbital_axis[1]=-sin(M_PI*23.44/180.0);
orbital_axis[2]=cos(M_PI*23.44/180.0);

/* crossproduct gives the vector perpedicular to both the average doppler shift and
  orbital axis */
band_axis[0]=orbital_axis[1]*average_det_velocity[2]-orbital_axis[2]*average_det_velocity[1];
band_axis[1]=orbital_axis[2]*average_det_velocity[0]-orbital_axis[0]*average_det_velocity[2];
band_axis[2]=orbital_axis[0]*average_det_velocity[1]-orbital_axis[1]*average_det_velocity[0];

/* Normalize */
band_axis_norm=sqrt(band_axis[0]*band_axis[0]+
	band_axis[1]*band_axis[1]+
	band_axis[2]*band_axis[2]);
/* replace 0.0 with something more reasonable later */
if(band_axis_norm<=0.0){
	band_axis[0]=0.0;
	band_axis[1]=0.0;
	band_axis[2]=1.0;
	} else {
	band_axis[0]/=band_axis_norm;
	band_axis[1]/=band_axis_norm;
	band_axis[2]/=band_axis_norm;
	}
/* 
  Normalize band_axis_norm so it matches the definition of \vec{u}
*/
band_axis_norm*=2.0*M_PI/(365.0*24.0*3600.0);

fprintf(LOG, "auto band axis norm: %g\n", band_axis_norm);
fprintf(LOG, "maximum S contribution from Doppler shifts: %g\n", band_axis_norm*(first_bin+nbins*0.5)/args_info.coherence_length_arg);

if(args_info.band_axis_norm_given) {
	band_axis_norm=args_info.band_axis_norm_arg;
	}

fprintf(LOG, "actual band axis norm: %g\n", band_axis_norm);


time(&input_done_time);
fprintf(LOG, "input complete: %d\n", (int)(input_done_time-start_time));
fprintf(stderr, "input complete: %d\n", (int)(input_done_time-start_time));

/* For testing - dump out emission time grid for SFT 0 */
if(0) {
EmissionTime emission_time;
int k=0, n=0;
LIGOTimeGPS tGPS;

TODO("analyze emission time dependence on sky location")

for(k=0;k<datasets[n].free;k+=10000) {
	for(i=0;i<main_grid->npoints; i++) {
		tGPS.gpsSeconds=datasets[n].gps[k]+datasets[n].coherence_time*0.5;
		tGPS.gpsNanoSeconds=0;
		get_emission_time(&emission_time, &(datasets[n].earth_state[k]), main_grid->longitude[i], main_grid->latitude[i], 0.0, datasets[n].detector, tGPS);
		fprintf(DATA_LOG, "emission_time: %d %.12g %.12g %lld %d %d\n", k, main_grid->longitude[i], main_grid->latitude[i], datasets[n].gps[k], emission_time.te.gpsSeconds, emission_time.te.gpsNanoSeconds);
		}
	}
fflush(DATA_LOG);
}


/* For testing - pick a bin and do phase correction */
{
double *re;
double *im;
double *power, *cross, *power_aux;
int *index;
COMPLEX16Vector *te_offsets;
double first_gps=min_gps();
int nsamples=1+ceil(2.0*(max_gps()-min_gps())/args_info.coherence_length_arg);
double total_weight;
double f0=args_info.focus_f0_arg;
//double ra=args_info.focus_ra_arg;
//double dec=args_info.focus_dec_arg;
double dInv=args_info.focus_dInv_arg;
int point;
int half_window=0;
int wing_step=round(nsamples*args_info.coherence_length_arg/SIDEREAL_DAY);
int day_samples=round(2.0*SIDEREAL_DAY/args_info.coherence_length_arg);
int fstep, nfsteps=4;
LOOSE_CONTEXT *ctx;

/* round of nsamples to contain integer number of sideral days */

ctx=create_context();

nsamples=day_samples*floor(nsamples/day_samples);

nsamples=ctx->nsamples;

fprintf(stderr, "nsamples=%d day_samples=%d wing_step=%d half_window=%d\n", nsamples, day_samples, wing_step, half_window);

te_offsets=XLALCreateCOMPLEX16Vector(nsamples);
if(!te_offsets) {
	fprintf(stderr, "**** ERROR: could not allocate te_offsets\n");
	exit(-1);
	}

power=do_alloc(nsamples, sizeof(*power));
cross=do_alloc(nsamples, sizeof(*cross));
power_aux=do_alloc(nsamples, sizeof(*power_aux));
index=do_alloc(nsamples, sizeof(*index));

for(fstep=0;fstep<nfsteps; fstep++) {
for(point=0;point< main_grid->npoints;point++) {
double ra=main_grid->longitude[point];
double dec=main_grid->latitude[point];

if(point % 10==0)fprintf(stderr, "point=%d out of %d\n", point, main_grid->npoints);

{
double x,y, x2, y2;
double c,s;
double dt, te;
double phase_spindown;
double phase_barycenter;
double phase_heterodyne;
double phase_bin;
double total_phase;
double f;
float f_plus, f_cross;
float e[GRID_E_COUNT];
int bin;
int i,j,k,n;
EmissionTime emission_time;
EmissionTime emission_time2;
LIGOTimeGPS tGPS;

//fprintf(stderr,"f0=%g bin=%d\n", f0, bin);
	
spindown=args_info.spindown_start_arg;

total_weight=0.0;

for(i=0;i<nsamples;i++) {
	ctx->plus_samples->data[i].re=0.0;
	ctx->plus_samples->data[i].im=0.0;
	ctx->cross_samples->data[i].re=0.0;
	ctx->cross_samples->data[i].im=0.0;
	}

for(i=0;i<nsamples;i++) {
	te_offsets->data[i].re=0.0;
	te_offsets->data[i].im=0.0;
	}

precompute_am_constants(e, ra, dec);

for(n=0;n<d_free;n++) {
	for(j=0;j<datasets[n].free;j++) {
		
	//	tGPS.gpsSeconds=datasets[n].gps[j]+datasets[n].coherence_time*0.5;
		tGPS.gpsSeconds=datasets[n].gps[j]; /* SFTs count phase from 0 */
		tGPS.gpsNanoSeconds=0;
		fast_get_emission_time(&emission_time, &(datasets[n].earth_state[j]), ra, dec, dInv, datasets[n].detector, tGPS);


		i=round(2.0*(datasets[n].gps[j]-first_gps)/args_info.coherence_length_arg);
/*		i=round(2.0*((emission_time.te.gpsSeconds-first_gps)+1e-9*emission_time.te.gpsNanoSeconds)/args_info.coherence_length_arg);*/
		if(i<0)continue;
		if(i>=nsamples)continue;
		
		dt=(emission_time.te.gpsSeconds-datasets[n].gps[j]-datasets[n].coherence_time*0.5)+1e-9*emission_time.te.gpsNanoSeconds;
		
		f=f0+((emission_time.te.gpsSeconds-spindown_start)+1e-9*emission_time.te.gpsNanoSeconds)*spindown+2.0*fstep/(nfsteps*nsamples*datasets[0].coherence_time);
		bin=round(datasets[0].coherence_time*f-first_bin);
				
		te=(emission_time.te.gpsSeconds-spindown_start)+1e-9*emission_time.te.gpsNanoSeconds;
/*		te_offsets->data[i].re=te;*/
				
		phase_spindown=0.5*te*te*spindown;
		
		phase_barycenter=(f0+2.0*fstep/(nfsteps*nsamples*datasets[0].coherence_time))*te;
		
		phase_heterodyne=(f0-(first_bin+bin+0.0)/datasets[0].coherence_time)*(datasets[n].gps[j]-spindown_start);

		phase_bin=0.5*(f0*datasets[0].coherence_time-(first_bin+bin));
		
		total_phase=2.0*M_PI*((phase_spindown -floor(phase_spindown))+(phase_barycenter-floor(phase_barycenter)));

		f_plus=F_plus_coeff(j, e, datasets[n].AM_coeffs_plus);
		f_cross=F_plus_coeff(j, e, datasets[n].AM_coeffs_cross);
		
		x=datasets[n].re[j*datasets[n].nbins+bin];
		y=datasets[n].im[j*datasets[n].nbins+bin];


		c=cos(total_phase);
		s=sin(total_phase);
		
// 		samples->data[i].re=x*c+y*s;
// 		samples->data[i].im=-x*s+y*c;

		x2=x*c+y*s;
		y2=-x*s+y*c;

		/* circular */
// 		samples->data[i].re=x2*(f_plus+f_cross)-y2*(f_plus-f_cross);
// 		samples->data[i].im=x2*(f_plus+f_cross)+y2*(f_plus-f_cross);
		//total_weight+=f_plus*f_plus+f_cross*f_cross;

		/* plus */
		/* Hann window */
		//c=0.5*(1.0-cos(2*M_PI*i/nsamples));
		
		ctx->plus_samples->data[i].re=x2*f_plus;
		ctx->plus_samples->data[i].im=y2*f_plus;

		ctx->cross_samples->data[i].re=x2*f_plus;
		ctx->cross_samples->data[i].im=y2*f_plus;
		//total_weight+=f_plus*f_plus+f_cross*f_cross;
		total_weight+=f_plus*f_plus;

		if(args_info.dump_stream_data_arg) {
			get_emission_time(&emission_time2, &(datasets[n].earth_state[j]), args_info.focus_ra_arg, args_info.focus_dec_arg, dInv, datasets[n].detector, tGPS);
			fprintf(DATA_LOG, "stream: %d %d %d %d %lld %d %d %d %d %.12g %.12g %.12g %.12g %.12f %.12f %.12f %.12f %d %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f\n",
				i, n, j, fstep, datasets[n].gps[j], 
				emission_time.te.gpsSeconds, emission_time.te.gpsNanoSeconds,
				emission_time2.te.gpsSeconds, emission_time2.te.gpsNanoSeconds,
				x, y, ctx->plus_samples->data[i].re, ctx->plus_samples->data[i].im, phase_spindown, phase_barycenter, phase_heterodyne, f, bin, f_plus, f_cross, total_phase,datasets[n].earth_state[j].posNow[0], datasets[n].earth_state[j].posNow[1], datasets[n].earth_state[j].posNow[2], ra, dec, args_info.focus_ra_arg, args_info.focus_dec_arg);
			}
		
		}
	}
//fprintf(stderr, "total_weight=%g\n", total_weight);
}
{
double norm, a, b, a2, b2, c, x, y, x1, y1;
double sum;
int i,j,k,m, i_filter, i_left, i_right;
SPARSE_CONV *sc, *sp_sc, *ra_sc, *dec_sc;
#define N_FREQ_ADJ_FILTER 9
COMPLEX16 filter[N_FREQ_ADJ_FILTER];
double timebase=max_gps()-min_gps();
POWER_STATS ps;
double sb_ra[2], sb_dec[2];

//fprintf(stderr, "Running FFT\n");
XLALCOMPLEX16VectorFFT(ctx->plus_fft, ctx->plus_samples, ctx->fft_plan);
XLALCOMPLEX16VectorFFT(ctx->cross_fft, ctx->cross_samples, ctx->fft_plan);
// XLALCOMPLEX16VectorFFT(te_fft, te_offsets, fft_plan);

sc=compute_te_offset_structure(datasets[0].detector, ra, dec, args_info.focus_dInv_arg, min_gps(), (max_gps()-min_gps())/4095, 4096);
sp_sc=compute_spindown_offset_structure(datasets[0].detector, ra, dec, args_info.focus_dInv_arg, min_gps(), (max_gps()-min_gps())/4095, 4096);

compute_sky_basis(ra, dec, resolution*0.5, sb_ra, sb_dec);
fprintf(stderr, "Sky basis: in=(%f, %f) out=(%f, %f) (%f, %f)\n", ra, dec, sb_ra[0], sb_dec[0], sb_ra[1], sb_dec[1]);
ra_sc=compute_sky_offset_structure(datasets[0].detector, ra, dec, sb_ra[0], sb_dec[0], args_info.focus_dInv_arg, min_gps(), (max_gps()-min_gps())/4095, 4096);
dec_sc=compute_sky_offset_structure(datasets[0].detector, ra, dec, sb_ra[1], sb_dec[1], args_info.focus_dInv_arg, min_gps(), (max_gps()-min_gps())/4095, 4096);

fprintf(stderr, "point %d step %d slope %g\n", point, fstep, sc->slope);
for(i=0;i<9;i++)
	fprintf(stderr, "  %g %g\n", sc->first9[i].re, sc->first9[i].im);

fprintf(stderr, "sp_sc point %d slope %g\n", point, sp_sc->slope);
for(i=0;i<9;i++)
	fprintf(stderr, "  %g %g\n", sp_sc->first9[i].re, sp_sc->first9[i].im);

fprintf(stderr, "ra_sc point %d slope %g\n", point, ra_sc->slope);
for(i=0;i<9;i++)
	fprintf(stderr, "  %g %g\n", ra_sc->first9[i].re, ra_sc->first9[i].im);

fprintf(stderr, "dec_sc point %d slope %g\n", point, dec_sc->slope);
for(i=0;i<9;i++)
	fprintf(stderr, "  %g %g\n", dec_sc->first9[i].re, dec_sc->first9[i].im);

norm=1.0/total_weight;
norm/=args_info.coherence_length_arg*16384.0;
norm*=2.0*sqrt(2.0); /* note - I do not understand where this extra factor of 2 comes from .. second fft ?? */
norm*=args_info.strain_norm_factor_arg;

for(i=0;i<nsamples;i++) {
	ctx->plus_fft->data[i].re*=norm;
	ctx->plus_fft->data[i].im*=norm;
	ctx->cross_fft->data[i].re*=norm;
	ctx->cross_fft->data[i].im*=norm;
	}

for(i=0;i<7;i++) {
	filter[i].re=0.0;
	filter[i].im=0.0;
	}
filter[3].re=1.0;
filter[3].im=0.0;
i_filter=0;	

for(i=0;i<nsamples;i++) {
	if(i>i_filter+200) {
		if((2*i)<nsamples)
			make_bessel_filter(filter, N_FREQ_ADJ_FILTER, sc->first9, 6, -i*2*M_PI/timebase);
			else 
			make_bessel_filter(filter, N_FREQ_ADJ_FILTER, sc->first9, 6, (nsamples-i)*2*M_PI/timebase);
/*		if(i<3000)for(j=0;j<7;j++) 
			fprintf(stderr, "i=%d (%f) j=%d %g %g\n", i, i/timebase, j, filter[j].re, filter[j].im);*/
		i_filter=i;
		}
/*	a=filter[6].re*fft->data[i-3].re-filter[6].im*fft->data[i-3].im;
	b=filter[6].re*fft->data[i-3].im+filter[6].im*fft->data[i-3].re;*/
	a=0.0;
	b=0.0;
	a2=0.0;
	b2=0.0;
	for(j=0;j<7;j++) {
		k=i-3+j;
		if(k<0)k=nsamples+k;
		if(k>=nsamples)k=k-nsamples;
		a+=filter[6-j].re*ctx->plus_fft->data[k].re-filter[6-j].im*ctx->plus_fft->data[k].im;
		b+=filter[6-j].re*ctx->plus_fft->data[k].im+filter[6-j].im*ctx->plus_fft->data[k].re;
		a2+=filter[6-j].re*ctx->cross_fft->data[k].re-filter[6-j].im*ctx->cross_fft->data[k].im;
		b2+=filter[6-j].re*ctx->cross_fft->data[k].im+filter[6-j].im*ctx->cross_fft->data[k].re;
		}
	
	ctx->plus_te_fft->data[i].re=a;
	ctx->plus_te_fft->data[i].im=b;
	ctx->cross_te_fft->data[i].re=a2;
	ctx->cross_te_fft->data[i].im=b2;
	}

scan_fft_stats(ctx, power, nsamples, ra_sc->first9, dec_sc->first9, &ps);

// -241
if((point== -21) || args_info.dump_fft_data_arg) {
	for(i=0;i<nsamples;i++) {
		fprintf(DATA_LOG, "fft: %d %d %.12f %.12f %.12f %.12g %.12g %.12g %.12g %.12g\n", 
			i, fstep, (i*2.0)/(nsamples*args_info.coherence_length_arg), ra, dec, ctx->plus_fft->data[i].re, ctx->plus_fft->data[i].im, ctx->plus_te_fft->data[i].re, ctx->plus_te_fft->data[i].im, power[i]);
		}
	}

if(fabs(ra-args_info.focus_ra_arg)<0.5*resolution && fabs(dec-args_info.focus_dec_arg)<0.5*resolution) {
	fprintf(stderr, "Close point: %d\n", point);
	}


fprintf(DATA_LOG, "max: %d %d %.12f %.12f %.12f %.12f %.12g %.12g %.12f %d %.12f\n",
	point, fstep, (i*2>nsamples ? i-nsamples : i)+fstep*1.0/nfsteps, ((i*2>nsamples ? i-nsamples : i)*1.0 +fstep*1.0/nfsteps)/(nsamples*args_info.coherence_length_arg), ra, dec, ps.mean, ps.sd, ps.max_snr, ps.max_snr_index, (power[0]-ps.mean)/ps.sd);

free_sparse_conv(sc);
free_sparse_conv(sp_sc);

}
}
}

}


time(&end_time);
fprintf(stderr, "exit memory: %g MB\n", (MEMUSAGE*10.0/(1024.0*1024.0))/10.0);
fprintf(LOG,"seconds elapsed: %ld\n",end_time-start_time);
fprintf(stderr,"seconds elapsed: %ld\n",end_time-start_time);
fclose(LOG);
fclose(FILE_LOG);
fclose(DATA_LOG);
return(0);
}
