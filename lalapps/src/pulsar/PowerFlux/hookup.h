#ifndef __HOOKUP_H__
#define __HOOKUP_H__

#include "global.h"

void init_hookup(void);
int clear_name_dat(char *name);
int clear_name_png(char *name);

void read_directory(char *prefix, long first,long last, 
	long first_bin,long bin_count,
	long *nsegments, float **power, INT64 **gps);

void get_patch_modulations(char *detresponse_path, char *ifo,
	int Nx, int Ny, int Nsegments, 
	INT64 *gps, double **patch_cross, double **patch_plus);

void get_modulations(char *detresponse_path, char *ifo,
	int Nx, int Ny, INT64 gps,
	 double *cross, double *plus, double *relfreq);

void dump_shorts(char *name, short *x, long count, long step);
void dump_ints(char *name, int *x, long count, long step);
void dump_floats(char *name, float *x, long count, long step);
void dump_doubles(char *name, double *x, long count, long step);

void init_ephemeris(void);
void get_AM_response(INT64 gps, float latitude, float longitude,
	float *plus, float *cross);
void get_detector_vel(INT64 gps, float *velocity);
/* there are count*GRID_FIT_COUNT coefficients */
void get_whole_sky_AM_response(INT64 *gps, long count, float **coeffs_plus, float **coeffs_cross, long *size);

#endif
