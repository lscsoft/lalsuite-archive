#ifndef __HOOKUP_H__
#define __HOOKUP_H__

#include "global.h"

void init_hookup(void);
int clear_name_dat(char *name);
int clear_name_png(char *name);

void read_directory(char *prefix, long first,long last, 
	long first_bin,long bin_count,
	long *nsegments, float **power, INT64 **gps);

void dump_shorts(char *name, short *x, long count, long step);
void dump_ints(char *name, int *x, long count, long step);
void dump_floats(char *name, float *x, long count, long step);
void dump_doubles(char *name, double *x, long count, long step);

void init_ephemeris(void);
void get_AM_response(INT64 gps, float latitude, float longitude, float orientation,
	float *plus, float *cross);
void get_detector_vel(INT64 gps, float *velocity);
/* there are count*GRID_FIT_COUNT coefficients */
void get_whole_sky_AM_response(INT64 *gps, long count, float orientation, float **coeffs_plus, float **coeffs_cross, long *size);
/* Accepts one set of coefficients for a fixed polarization */
void verify_whole_sky_AM_response(INT64 *gps, long count, float orientation,  SKY_GRID *grid, float *coeffs_plus, char *name);

#endif
