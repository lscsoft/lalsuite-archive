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
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <regex.h>

#include <lal/LALRCSID.h>
#include <lal/LALInitBarycenter.h>
#include <lal/DetResponse.h>
#include <lal/Velocity.h>
#include <lal/DetectorSite.h>

#include <gsl/gsl_multifit.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>

#include "global.h"
#include "cmdline.h"
#include "grid.h"
#include "intervals.h"
#include "polarization.h"
#include "hookup.h"

#ifndef PATH_MAX
/* just in case it is not defined */
#define PATH_MAX PATH_MAX
#endif

extern FILE *LOG, *FILE_LOG;
extern char *earth_ephemeris;
extern char *sun_ephemeris;
extern struct gengetopt_args_info args_info;
extern char *output_dir;

regex_t write_dat, write_png;

EphemerisData ephemeris;
LALDetector detector;

void init_hookup(void)
{
if(regcomp(&write_dat, args_info.write_dat_arg, REG_EXTENDED | REG_NOSUB)){
	fprintf(stderr,"Cannot compile \"--write-dat=%s\"\n", args_info.write_dat_arg);
	exit(-1);
	}
if(regcomp(&write_png, args_info.write_png_arg, REG_EXTENDED | REG_NOSUB)){
	fprintf(stderr,"Cannot compile \"--write-dat=%s\"\n", args_info.write_dat_arg);
	exit(-1);
	}
}

int clear_name_dat(char *name)
{
return !regexec(&write_dat, name, 0, NULL, 0);
}

int clear_name_png(char *name)
{
return !regexec(&write_png, name, 0, NULL, 0);
}

void dump_shorts(char *name, short *x, long count, long step)
{
FILE *fout;
long i;
char s[PATH_MAX];

if(x==NULL) {
	fprintf(stderr, "Skipping %s\n", name);
	return;
	}
if(!clear_name_dat(name))return;

snprintf(s,PATH_MAX,"%s%s", output_dir, name);
fout=fopen(s, "w");
if(fout==NULL){
	fprintf(FILE_LOG, "Could not open file \"%s\" for writing.\n",
		name);
	return;
	}
/* important special case */
if(step==1)fwrite(x, sizeof(*x), count,fout);
	else
	for(i=0;i<count;i++)fwrite(x+i*step, sizeof(*x), 1,fout);
fclose(fout);
fprintf(FILE_LOG, "shorts: %ld %s\n", count, name);
}

void dump_ints(char *name, int *x, long count, long step)
{
FILE *fout;
long i;
char s[PATH_MAX];

if(x==NULL) {
	fprintf(stderr, "Skipping %s\n", name);
	return;
	}
if(!clear_name_dat(name))return;

snprintf(s,PATH_MAX,"%s%s", output_dir, name);
fout=fopen(s, "w");
if(fout==NULL){
	fprintf(FILE_LOG, "Could not open file \"%s\" for writing.\n",
		name);
	return;
	}
/* important special case */
if(step==1)fwrite(x, sizeof(*x), count,fout);
	else
	for(i=0;i<count;i++)fwrite(x+i*step, sizeof(*x), 1,fout);
fclose(fout);
fprintf(FILE_LOG, "ints: %ld %s\n", count, name);
}

void dump_floats(char *name, float *x, long count, long step)
{
FILE *fout;
long i;
char s[PATH_MAX];

if(x==NULL) {
	fprintf(stderr, "Skipping %s\n", name);
	return;
	}
if(!clear_name_dat(name))return;

snprintf(s,PATH_MAX,"%s%s", output_dir, name);
fout=fopen(s, "w");
if(fout==NULL){
	fprintf(FILE_LOG, "Could not open file \"%s\" for writing.\n",
		name);
	return;
	}
/* important special case */
if(step==1)fwrite(x, sizeof(*x), count,fout);
	else
	for(i=0;i<count;i++)fwrite(x+i*step, sizeof(*x), 1,fout);
fclose(fout);
fprintf(FILE_LOG, "floats: %ld %s\n", count, name);
}

void dump_doubles(char *name, double *x, long count, long step)
{
FILE *fout;
long i;
char s[PATH_MAX];

if(x==NULL) {
	fprintf(stderr, "Skipping %s\n", name);
	return;
	}
if(!clear_name_dat(name))return;

snprintf(s,PATH_MAX,"%s%s", output_dir, name);
fout=fopen(s, "w");
if(fout==NULL){
	fprintf(FILE_LOG, "Could not open file \"%s\" for writing.\n",
		name);
	return;
	}
/* important special case */
if(step==1)fwrite(x, sizeof(*x), count,fout);
	else
	for(i=0;i<count;i++)fwrite(x+i*step, sizeof(*x), 1,fout);
fclose(fout);
fprintf(FILE_LOG, "doubles: %ld %s\n", count, name);
}


void get_detector(char *det) 
{
if(!strcasecmp("LHO", det)){
	detector=lalCachedDetectors[LALDetectorIndexLHODIFF];
	} else 
if(!strcasecmp("LLO", det)){
	detector=lalCachedDetectors[LALDetectorIndexLLODIFF];
	} else
if(!strcasecmp("VIRGO", det)){
	detector=lalCachedDetectors[LALDetectorIndexVIRGODIFF];
	} else {
	fprintf(stderr,"Unrecognized detector site: \"%s\"\n", args_info.detector_arg);
	exit(-1);
	}
}

LALDetector get_detector_struct(char *det) 
{
if(!strcasecmp("LHO", det)){
	return lalCachedDetectors[LALDetectorIndexLHODIFF];
	} else 
if(!strcasecmp("LLO", det)){
	return lalCachedDetectors[LALDetectorIndexLLODIFF];
	} else
if(!strcasecmp("VIRGO", det)){
	return lalCachedDetectors[LALDetectorIndexVIRGODIFF];
	} else {
	fprintf(stderr,"Unrecognized detector site: \"%s\"\n", args_info.detector_arg);
	exit(-1);
	}
}

void init_ephemeris(void)
{
LALStatus status={level:0, statusPtr:NULL};

#if 0
get_detector(args_info.detector_arg);
fprintf(LOG,"detector  : %s (%s)\n", args_info.detector_arg, detector.frDetector.name);
#endif

memset(&ephemeris, 0, sizeof(ephemeris));
ephemeris.ephiles.earthEphemeris=earth_ephemeris;
ephemeris.ephiles.sunEphemeris=sun_ephemeris;
LALInitBarycenter(&status, &ephemeris);
TESTSTATUS(&status);

#if 0
detectorvel_inputs.detector=detector;
detectorvel_inputs.edat=&ephemeris;
#endif
fprintf(stderr,"Successfully initialized ephemeris data\n");
if(status.statusPtr)FREESTATUSPTR(&status);
}

void get_AM_response(INT64 gps, float latitude, float longitude, float orientation,
	float *plus, float *cross)
{
LALStatus status={level:0, statusPtr:NULL};
LALSource source;
LALDetAndSource det_and_source={NULL, NULL};
LALDetAMResponse response;
LIGOTimeGPS ligo_gps;

memset(&ligo_gps, 0, sizeof(ligo_gps));
ligo_gps.gpsSeconds=gps; 
ligo_gps.gpsNanoSeconds=0;

memset(&source, 0, sizeof(source));
source.equatorialCoords.system=COORDINATESYSTEM_EQUATORIAL;
source.orientation=orientation;
source.equatorialCoords.longitude=longitude;
source.equatorialCoords.latitude=latitude;

det_and_source.pDetector=&detector;
det_and_source.pSource=&source;

LALComputeDetAMResponse(&status, &response, &det_and_source, &ligo_gps);
TESTSTATUS(&status);

*cross=response.cross;
*plus=response.plus;
if(status.statusPtr)FREESTATUSPTR(&status);
}

void get_AM_response_d(double gps, float latitude, float longitude, float orientation, char *detector,
	float *plus, float *cross)
{
LALStatus status={level:0, statusPtr:NULL};
LALSource source;
LALDetAndSource det_and_source={NULL, NULL};
LALDetAMResponse response;
LIGOTimeGPS ligo_gps;
LALDetector det;

memset(&ligo_gps, 0, sizeof(ligo_gps));
ligo_gps.gpsSeconds=floor(gps); 
ligo_gps.gpsNanoSeconds=(gps-floor(gps))*1e9;

memset(&source, 0, sizeof(source));
source.equatorialCoords.system=COORDINATESYSTEM_EQUATORIAL;
source.orientation=orientation;
source.equatorialCoords.longitude=longitude;
source.equatorialCoords.latitude=latitude;

det=get_detector_struct(detector);

det_and_source.pDetector=&det;
det_and_source.pSource=&source;

LALComputeDetAMResponse(&status, &response, &det_and_source, &ligo_gps);
TESTSTATUS(&status);

*cross=response.cross;
*plus=response.plus;
if(status.statusPtr)FREESTATUSPTR(&status);
}

void get_detector_vel(INT64 gps, float *velocity)
{
LALStatus status={level:0, statusPtr:NULL};
REAL8 det_velocity[3];
int i;
LIGOTimeGPS ligo_gps;

memset(&ligo_gps, 0, sizeof(ligo_gps));

ligo_gps.gpsSeconds=gps; 
ligo_gps.gpsNanoSeconds=0;
LALDetectorVel(&status, det_velocity, &ligo_gps, detector, &ephemeris);
TESTSTATUS(&status);

#if 0
fprintf(stderr,"powerflux: det_velocity=(%g,%g,%g)\n", 
	det_velocity[0],
	det_velocity[1],
	det_velocity[2]
	);
fprintf(stderr,"gps=%d (nano=%d)\n",ligo_gps.gpsSeconds, ligo_gps.gpsNanoSeconds);
fprintf(stderr,"detector=%s\n", detector.frDetector.name);
fprintf(stderr,"powerflux nE=%d nS=%d dE=%g dS=%g\n",
	ephemeris.nentriesE,
	ephemeris.nentriesS,
	ephemeris.dtEtable,
	ephemeris.dtStable);
#endif

for(i=0;i<3;i++)velocity[i]=det_velocity[i];
if(status.statusPtr)FREESTATUSPTR(&status);
}

void get_emission_time(EmissionTime *emission_time, EarthState *earth_state, double ra, double dec, double dInv, char *detector, LIGOTimeGPS tGPS)
{
BarycenterInput baryinput;
LALStatus status={level:0, statusPtr:NULL};

memset(&baryinput, 0, sizeof(baryinput));
baryinput.tgps=tGPS;
baryinput.site=get_detector_struct(detector);
baryinput.site.location[0]=baryinput.site.location[0]/LAL_C_SI;
baryinput.site.location[1]=baryinput.site.location[1]/LAL_C_SI;
baryinput.site.location[2]=baryinput.site.location[2]/LAL_C_SI;
baryinput.alpha=ra;
baryinput.delta=dec;
baryinput.dInv=dInv;

LALBarycenter(&status, emission_time, &baryinput, earth_state);
TESTSTATUS(&status);
if(status.statusPtr)FREESTATUSPTR(&status);
}


/* there are count*GRID_FIT_COUNT coefficients */
void get_whole_sky_AM_response(INT64 *gps, double coherence_time, char *detector, long count, float orientation, float **coeffs_plus, float **coeffs_cross, long *size)
{
int i, j, k, last_k;
SKY_GRID *sample_grid=NULL;
float plus, cross;

gsl_multifit_linear_workspace *workspace=NULL;
gsl_matrix *X=NULL, *cov=NULL;
gsl_vector *y_plus=NULL, *y_cross=NULL, *c=NULL;
double chisq;

fprintf(stderr,"Computing whole sky AM response for %ld SFTs\n", count);
fprintf(stderr, "AM coeffs size: %f KB\n", 2*count*GRID_FIT_COUNT*sizeof(**coeffs_plus)/1024.0);
fprintf(LOG, "AM coeffs size: %f KB\n", 2*count*GRID_FIT_COUNT*sizeof(**coeffs_plus)/1024.0);
*size=count*GRID_FIT_COUNT;
*coeffs_plus=do_alloc(*size, sizeof(**coeffs_plus));
*coeffs_cross=do_alloc(*size, sizeof(**coeffs_cross));

sample_grid=make_rect_grid(8,5);
precompute_values(sample_grid);

workspace=gsl_multifit_linear_alloc(sample_grid->npoints, GRID_FIT_COUNT);

y_plus=gsl_vector_alloc(sample_grid->npoints);
y_cross=gsl_vector_alloc(sample_grid->npoints);
c=gsl_vector_alloc(GRID_FIT_COUNT);

cov=gsl_matrix_alloc(GRID_FIT_COUNT, GRID_FIT_COUNT);
X=gsl_matrix_alloc(sample_grid->npoints, GRID_FIT_COUNT);

last_k=-1000;
for(k=0;k<count;k++){
	TODO("1 minute is good enough, but it would be better to do this with rounding");
	if(last_k>=0 && (gps[k]<gps[last_k]+60) && (gps[k]+60>gps[last_k])) {
		for(j=0;j<GRID_FIT_COUNT;j++){
			(*coeffs_plus)[k*GRID_FIT_COUNT+j]=(*coeffs_plus)[last_k*GRID_FIT_COUNT+j];
			(*coeffs_cross)[k*GRID_FIT_COUNT+j]=(*coeffs_cross)[last_k*GRID_FIT_COUNT+j];
			}
		continue;
		}
	
	for(i=0;i<sample_grid->npoints;i++){
		get_AM_response_d(gps[k]+coherence_time*0.5, 
			sample_grid->latitude[i], sample_grid->longitude[i], 
			orientation, detector,
			&plus, &cross);
		gsl_vector_set(y_plus, i, plus);
		gsl_vector_set(y_cross, i, cross);
		
		for(j=0;j<GRID_FIT_COUNT;j++){
			gsl_matrix_set(X, i, j, sample_grid->e[j+GRID_FIT_START][i]);
			}
		}
	gsl_multifit_linear(X, y_plus, c, cov, &chisq, workspace);
	if(chisq>1e-12){
		fprintf(stderr,"** Sky grid approximation fault: non-zero chisq %g when computing gps[%d]=%lld, plus polarization - aborting !\n", 
			chisq, k, gps[k]);
		fprintf(LOG,"** Sky grid approximation fault: non-zero chisq %g when computing gps[%d]=%lld, plus polarization - aborting !\n", 
			chisq, k, gps[k]);
		exit(-1);
		}
	for(j=0;j<GRID_FIT_COUNT;j++){
		(*coeffs_plus)[k*GRID_FIT_COUNT+j]=gsl_vector_get(c, j);
		}

	gsl_multifit_linear(X, y_cross, c, cov, &chisq, workspace);
	if(chisq>1e-12){
		fprintf(stderr,"** Sky grid approximation fault: non-zero chisq %g when computing gps[%d]=%lld, cross polarization - aborting !\n", 
			chisq, k, gps[k]);
		fprintf(LOG,"** Sky grid approximation fault: non-zero chisq %g when computing gps[%d]=%lld, cross polarization - aborting !\n", 
			chisq, k, gps[k]);
		exit(-1);
		}
	for(j=0;j<GRID_FIT_COUNT;j++){
		(*coeffs_cross)[k*GRID_FIT_COUNT+j]=gsl_vector_get(c, j);
		}
	last_k=k;
	}

gsl_matrix_free(X);
gsl_matrix_free(cov);
gsl_vector_free(y_plus);
gsl_vector_free(y_cross);
gsl_vector_free(c);

gsl_multifit_linear_free(workspace);
free_grid(sample_grid);
}

/* there are count*GRID_FIT_COUNT coefficients */
void verify_whole_sky_AM_response(INT64 *gps, double coherence_time, char *detector, long count, float orientation,  SKY_GRID *sample_grid, float *coeffs_plus, char *name)
{
int i,j;
long offset;
float plus, cross;
float max_err,err;
gsl_rng *rng=NULL;
rng=gsl_rng_alloc(gsl_rng_default);

max_err=0;
for(i=0;i<count;i++){
	/* 20 points per segment ought to be enough */
	for(j=0;j<20;j++){
		/* test in random grid points */
		offset=floor(sample_grid->npoints*gsl_rng_uniform(rng));
		get_AM_response_d(gps[i]+coherence_time*0.5, 
			sample_grid->latitude[offset], sample_grid->longitude[offset], 
			orientation, detector,
			&plus, &cross);
		err=fabs(plus*plus-AM_response(i, sample_grid, offset, coeffs_plus));
		if(err>max_err)max_err=err;
		}
	}
fprintf(stderr, "%s AM coeffs error: %g\n", name, max_err);
fprintf(LOG, "%s AM coeffs error: %g\n", name, max_err);
gsl_rng_free(rng);
}
