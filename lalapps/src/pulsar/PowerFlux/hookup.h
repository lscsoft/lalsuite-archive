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

#ifndef __HOOKUP_H__
#define __HOOKUP_H__

#include "global.h"

#include <lal/DetectorSite.h>
#include <lal/LALBarycenter.h>
#include <lal/LALDetectors.h>

void init_hookup(void);
int clear_name_dat(char *name);
int clear_name_png(char *name);

void read_directory(char *prefix, int first,int last, 
	int first_bin,int bin_count,
	int *nsegments, float **power, INT64 **gps);

void dump_shorts(char *name, short *x, long count, long step);
void dump_ints(char *name, int *x, long count, long step);
void dump_floats(char *name, float *x, long count, long step);
void dump_doubles(char *name, double *x, long count, long step);

void init_ephemeris(void);
void get_detector(char *det);
LALDetector get_detector_struct(char *det);
void get_AM_response(INT64 gps, float latitude, float longitude, float orientation,
	float *plus, float *cross);
void get_AM_response_d(double gps, float latitude, float longitude, float orientation, char *detector,
	float *plus, float *cross);
void get_detector_vel(INT64 gps, float *velocity);
/* there are count*GRID_FIT_COUNT coefficients */
void get_whole_sky_AM_response(INT64 *gps, long count, float orientation, float **coeffs_plus, float **coeffs_cross, long *size);
/* Accepts one set of coefficients for a fixed polarization */
void verify_whole_sky_AM_response(INT64 *gps, long count, float orientation,  SKY_GRID *grid, float *coeffs_plus, char *name);
void get_emission_time(EmissionTime *emission_time, EarthState *earth_state, double ra, double dec, double dInv, char *detector, LIGOTimeGPS tGPS);

#endif
