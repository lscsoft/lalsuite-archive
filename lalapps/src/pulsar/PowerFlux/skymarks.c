#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "global.h"
#include "skymarks.h"
#include "util.h"
#include "grid.h"
#include "dataset.h"
#include "cmdline.h"

extern struct gengetopt_args_info args_info;

extern FILE *LOG;
extern double spindown;

extern DATASET *datasets;
extern int d_free;
extern int nbins, first_bin, side_cut, useful_bins;

INT64 spindown_start;

int find_closest(SKY_GRID *sky_grid, SKY_GRID_TYPE ra, SKY_GRID_TYPE dec)
{
SKY_GRID_TYPE ds, ds_min;
int i,k;

k=-1;
ds_min=100;

for(i=0;i<sky_grid->npoints;i++) {
	ds=acos(sin(sky_grid->latitude[i])*sin(dec)+
		cos(sky_grid->latitude[i])*cos(dec)*
		cos(sky_grid->longitude[i]-ra));

	if((ds<ds_min)) {
		ds_min=ds;
		k=i;
		}
	}
return k;
}

void parse_skymarks_line(char *line, int length, SKY_GRID *sky_grid, SKYMARK *sm)
{
int ai, aj;

sm->type=SKYMARK_NOP;

/* skip whitespace in the beginning */
while(((*line)==' ') || ((*line)=='\t'))line++;
/* skip comments */
if((*line)=='#')return;
/* skip empty lines */
if((*line)=='\n')return;
if((*line)=='\r')return;
if((*line)==0)return;

/* General format of the command:
	command band_to band_from [other_args]
*/

locate_arg(line, length, 1, &ai, &aj);
sm->band_to=add_band(sky_grid, &(line[ai]), aj-ai);

locate_arg(line, length, 2, &ai, &aj);
sm->band_from=add_band(sky_grid, &(line[ai]), aj-ai);

if(!strncasecmp(line, "disk", 4)) {
	sm->type=SKYMARK_DISK;
	
	locate_arg(line, length, 3, &ai, &aj);
	sscanf(&(line[ai]), "%g", &sm->p.disk.ra);

	locate_arg(line, length, 4, &ai, &aj);
	sscanf(&(line[ai]), "%g", &sm->p.disk.dec);

	locate_arg(line, length, 5, &ai, &aj);
	sscanf(&(line[ai]), "%g", &sm->p.disk.radius);

	sm->p.disk.closest_point=find_closest(sky_grid, sm->p.disk.ra, sm->p.disk.dec);
	if(sm->p.disk.radius<0)sm->p.disk.cos_radius=100;
		else
	if(sm->p.disk.radius>M_PI)sm->p.disk.cos_radius=-100.0;
		else
		sm->p.disk.cos_radius=cos(sm->p.disk.radius);

	fprintf(stderr, "Compiled disk (%d <- %d) around (%g, %g) with radius %g\n", sm->band_to, sm->band_from, sm->p.disk.ra, sm->p.disk.dec, sm->p.disk.radius);
	fprintf(LOG, "Compiled disk (%d <- %d) around (%g, %g) with radius %g\n", sm->band_to, sm->band_from, sm->p.disk.ra, sm->p.disk.dec, sm->p.disk.radius);

	// TODO: find closest sky position explicitly
	} else
if(!strncasecmp(line, "band", 4)) {
	sm->type=SKYMARK_BAND;

	locate_arg(line, length, 3, &ai, &aj);
	sscanf(&(line[ai]), "%g", &sm->p.band.ra);

	locate_arg(line, length, 4, &ai, &aj);
	sscanf(&(line[ai]), "%g", &sm->p.band.dec);

	locate_arg(line, length, 5, &ai, &aj);
	sscanf(&(line[ai]), "%g", &sm->p.band.level1);

	locate_arg(line, length, 6, &ai, &aj);
	sscanf(&(line[ai]), "%g", &sm->p.band.level2);

	sm->p.band.x0=cos(sm->p.band.ra)*sin(M_PI_2-sm->p.band.dec);
	sm->p.band.y0=sin(sm->p.band.ra)*sin(M_PI_2-sm->p.band.dec);
	sm->p.band.z0=cos(M_PI_2-sm->p.band.dec);

	fprintf(stderr, "Compiled band (%d <- %d) around (%g, %g) with cos in [%g, %g]\n", sm->band_to, sm->band_from, sm->p.band.ra, sm->p.band.dec, sm->p.band.level1, sm->p.band.level2);
	fprintf(LOG, "Compiled band (%d <- %d) around (%g, %g) with cos in [%g, %g]\n", sm->band_to, sm->band_from, sm->p.band.ra, sm->p.band.dec, sm->p.band.level1, sm->p.band.level2);
	} else 
if(!strncasecmp(line, "closest", 7)) {
	sm->type=SKYMARK_CLOSEST;

	locate_arg(line, length, 3, &ai, &aj);
	sscanf(&(line[ai]), "%g", &sm->p.closest.ra);

	locate_arg(line, length, 4, &ai, &aj);
	sscanf(&(line[ai]), "%g", &sm->p.closest.dec);

	sm->p.closest.closest_point=find_closest(sky_grid, sm->p.disk.ra, sm->p.disk.dec);

	fprintf(stderr, "Compiled point (%d <- %d) closest to (%g, %g)\n", sm->band_to, sm->band_from, sm->p.closest.ra, sm->p.closest.dec);
	fprintf(LOG, "Compiled point (%d <- %d) closest to (%g, %g)\n", sm->band_to, sm->band_from, sm->p.closest.ra, sm->p.closest.dec);
	
	// TODO: find closest sky position explicitly

	} else
if(!strncasecmp(line, "response", 8)) {
	sm->type=SKYMARK_RESPONSE;

	locate_arg(line, length, 3, &ai, &aj);
	sscanf(&(line[ai]), "%g", &sm->p.response.ra);

	locate_arg(line, length, 4, &ai, &aj);
	sscanf(&(line[ai]), "%g", &sm->p.response.dec);

	locate_arg(line, length, 5, &ai, &aj);
	sscanf(&(line[ai]), "%g", &sm->p.response.weight_ratio_level);

	locate_arg(line, length, 6, &ai, &aj);
	sscanf(&(line[ai]), "%g", &sm->p.response.bin_tolerance);

	locate_arg(line, length, 7, &ai, &aj);
	sscanf(&(line[ai]), "%g", &sm->p.response.spindown_tolerance);

	fprintf(stderr, "Compiled points (%d <- %d) swept by (%g, %g) weight_ratio=%g bin_width=%g spindown_width=%g\n", sm->band_to, sm->band_from, sm->p.response.ra, sm->p.response.dec, sm->p.response.weight_ratio_level, sm->p.response.bin_tolerance, sm->p.response.spindown_tolerance);
	fprintf(LOG, "Compiled points (%d <- %d) swept by (%g, %g) weight_ratio=%g bin_width=%g spindown_width=%g\n", sm->band_to, sm->band_from, sm->p.response.ra, sm->p.response.dec, sm->p.response.weight_ratio_level, sm->p.response.bin_tolerance, sm->p.response.spindown_tolerance);

	} else
if(!strncasecmp(line, "echo_response", 13)) {
	sm->type=SKYMARK_ECHO_RESPONSE;

	locate_arg(line, length, 3, &ai, &aj);
	sscanf(&(line[ai]), "%g", &sm->p.echo_response.ra);

	locate_arg(line, length, 4, &ai, &aj);
	sscanf(&(line[ai]), "%g", &sm->p.echo_response.dec);

	locate_arg(line, length, 5, &ai, &aj);
	sscanf(&(line[ai]), "%g", &sm->p.echo_response.ref_spindown);

	locate_arg(line, length, 6, &ai, &aj);
	sscanf(&(line[ai]), "%g", &sm->p.echo_response.weight_ratio_level);

	locate_arg(line, length, 7, &ai, &aj);
	sscanf(&(line[ai]), "%g", &sm->p.echo_response.bin_tolerance);

	locate_arg(line, length, 8, &ai, &aj);
	sscanf(&(line[ai]), "%g", &sm->p.echo_response.spindown_tolerance);

	fprintf(stderr, "Compiled points (%d <- %d) swept by (%g, %g) spindown=%g weight_ratio=%g bin_width=%g spindown_width=%g\n", sm->band_to, sm->band_from, sm->p.echo_response.ra, sm->p.echo_response.dec, sm->p.echo_response.ref_spindown, sm->p.echo_response.weight_ratio_level, sm->p.echo_response.bin_tolerance, sm->p.echo_response.spindown_tolerance);
	fprintf(LOG, "Compiled points (%d <- %d) swept by (%g, %g) spindown=%g weight_ratio=%g bin_width=%g spindown_width=%g\n", sm->band_to, sm->band_from, sm->p.echo_response.ra, sm->p.echo_response.dec, sm->p.echo_response.ref_spindown, sm->p.echo_response.weight_ratio_level, sm->p.echo_response.bin_tolerance, sm->p.echo_response.spindown_tolerance);

	} else
if(!strncasecmp(line, "line_response", 13)) {
	sm->type=SKYMARK_LINE_RESPONSE;

	locate_arg(line, length, 3, &ai, &aj);
	sscanf(&(line[ai]), "%g", &sm->p.line_response.weight_ratio_level);

	locate_arg(line, length, 4, &ai, &aj);
	sscanf(&(line[ai]), "%g", &sm->p.line_response.bin_tolerance);

	fprintf(stderr, "Compiled points (%d <- %d) swept by lines weight_ratio=%g bin_width=%g\n", sm->band_to, sm->band_from, sm->p.line_response.weight_ratio_level, sm->p.line_response.bin_tolerance);
	fprintf(LOG, "Compiled points (%d <- %d) swept by lines weight_ratio=%g bin_width=%g\n", sm->band_to, sm->band_from, sm->p.line_response.weight_ratio_level, sm->p.line_response.bin_tolerance);
	} else
	{
	fprintf(stderr, "*** UNKNOWN masking command \"%s\"\n", line);
	sm->type=SKYMARK_ERROR;
	exit(-1);
	}
}

SKYMARK * compile_marks(SKY_GRID *sky_grid, char *s, int length)
{
int ai, aj;
ai=0;
aj=0;
SKYMARK *sm;
int sm_size;
int sm_free=0;

sm_size=10;
sm=do_alloc(sm_size, sizeof(*sm));

while(aj<length) {
	ai=aj;
	while(s[aj] && s[aj]!='\n' && (aj<length))aj++;

	if(sm_free+1>=sm_size) {
		SKYMARK *p;
		sm_size*=2;
		p=do_alloc(sm_size, sizeof(*sm));
		memcpy(p, sm, sm_free*sizeof(*sm));
		free(sm);
		sm=p;
		}

	parse_skymarks_line(&(s[ai]), aj-ai, sky_grid, &(sm[sm_free]));
	if(sm[sm_free].type!=SKYMARK_NOP)sm_free++;
	aj++;
	}
sm[sm_free].type=SKYMARK_END;
return sm;
}

float fast_stationary_effective_weight_ratio(float *e, float spindown, float bin_tolerance)
{
int i, k;
float w, fdiff;
float total_weight=0.0;
float f1=0.0;
DATASET *d;
float offset;
float center_frequency=first_bin+nbins*0.5;
float doppler_mult=args_info.doppler_multiplier_arg;

/* fit first */
f1=0;
for(i=0;i<d_free;i++) {
	d=&(datasets[i]);
	for(k=0;k<d->free;k++) {
		w=d->expTMedians[k]*d->weight;
		total_weight+=w;

		fdiff=center_frequency*(
				e[0]*d->detector_velocity[3*k+0]+
				e[1]*d->detector_velocity[3*k+1]+
				e[2]*d->detector_velocity[3*k+2]
				)*doppler_mult
			+d->coherence_time*spindown*(d->gps[k]-spindown_start);
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

		fdiff=center_frequency*(
				e[0]*d->detector_velocity[3*k+0]+
				e[1]*d->detector_velocity[3*k+1]+
				e[2]*d->detector_velocity[3*k+2]
				)*doppler_mult
			+d->coherence_time*spindown*(d->gps[k]-spindown_start)
			-offset;

		if(fabsf(fdiff)<bin_tolerance)f1+=w;
		}
	}
/* fprintf(stderr, "%g %g\n", f1, total_weight); */

/* TODO: 0.5 is here to fixup a code issue in earlier revision - we were accumulating total_weight twice, 
   this does not affect results since weight_ratio was determined from runs on actual data, but this does make the return value more confusing */
return 0.5*f1/total_weight;
}

int mark_sky_point(SKYMARK *sm, int point, float ra, float dec, float *e, float spindown1)
{
int skyband=-1;
float ds;

if(args_info.focus_radius_given) {
	if(args_info.focus_dec_given && 
	args_info.focus_ra_given) {

		if((sin(dec)*sin((float)args_info.focus_dec_arg)+
			cos(dec)*cos((float)args_info.focus_dec_arg)*
			cos(ra-(float)args_info.focus_ra_arg)) < cos((float)args_info.focus_radius_arg)) return -1;
		}
	}

/* TODO:
if(args_info.only_large_cos_given) {
	fprintf(LOG, "only large cos level: %f\n", args_info.only_large_cos_arg);
	mask_small_cos(fine_grid, band_axis[0], band_axis[1], band_axis[3], args_info.only_large_cos_arg);
	propagate_far_points_from_super_grid(patch_grid, proto_super_grid);
	}
*/

while(sm->type!=SKYMARK_END) {
	if((sm->band_from>=0) && skyband!=sm->band_from)continue;

	switch(sm->type) {
		case SKYMARK_CLOSEST:
			if(point==sm->p.closest.closest_point)skyband=sm->band_to;
			break;
		case SKYMARK_DISK:
			if(point==sm->p.disk.closest_point) {
				skyband=sm->band_to;
				break;
				}
			if((sin(dec)*sin(sm->p.disk.dec)+
				cos(dec)*cos(sm->p.disk.dec)*
				cos(ra-sm->p.disk.ra))>sm->p.disk.cos_radius) {
				skyband=sm->band_to;
				}
			break;
		case SKYMARK_BAND:
			ds=cos(dec)*cos(ra)*sm->p.band.x0+
				cos(dec)*sin(ra)*sm->p.band.y0+
				sin(dec)*sm->p.band.z0;
	
			if((ds>sm->p.band.level1) && (ds<=sm->p.band.level2)) {
				skyband=sm->band_to;
				}
			break;
		case SKYMARK_RESPONSE:
			spindown=spindown1; /* TODO: change code in dataset.c to pass spindown value explicitly instead of relying on built-in version */
			if(effective_weight_ratio(ra, dec, sm->p.response.ra, sm->p.response.dec, spindown1, sm->p.response.bin_tolerance, sm->p.response.spindown_tolerance)> sm->p.response.weight_ratio_level) {
				skyband=sm->band_to;
				}
			break;
		case SKYMARK_ECHO_RESPONSE:
			spindown=spindown1;
			if(effective_weight_ratio(ra, dec, sm->p.echo_response.ra, sm->p.echo_response.dec, sm->p.echo_response.ref_spindown, sm->p.echo_response.bin_tolerance, sm->p.echo_response.spindown_tolerance)> sm->p.echo_response.weight_ratio_level) {
				skyband=sm->band_to;
				}
			break;
		case SKYMARK_LINE_RESPONSE:
			if(fast_stationary_effective_weight_ratio(e, spindown1, sm->p.line_response.bin_tolerance)>sm->p.line_response.weight_ratio_level) {
				skyband=sm->band_to;
				}
			break;
		}
	sm++;
	}
return skyband;
}
