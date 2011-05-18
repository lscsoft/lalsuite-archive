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

#ifndef __DATASET_H__
#define __DATASET_H__

#include "global.h"
#include "polarization.h"
#include "lines.h"
#include "intervals.h"

#include <lal/DetectorSite.h>
#include <lal/LALBarycenter.h>
#include <lal/LALDetectors.h>

#define MAX_LOCKS 10

typedef struct {
	/* name of the data set - for logging purposes */
	char *name; 
	char *lock_file[MAX_LOCKS];

	int validated; /* has this data set been validated ? see corresponding functions */

	/* these describe which time periods we are interested in */
	INTERVAL_SET *segment_list;
	INTERVAL_SET *veto_segment_list;

	/* when on, ignore segments with duplicate GPS times - useful for implementing poor mans RAID of the nodes */
	int no_duplicate_gps;

	/* detector this dataset was observed on */
	char *detector;

	INT64 gps_start;
	INT64 gps_stop;
	
	INT64 max_gps;

	INT64 *gps;

	float veto_level;
	float veto_spike_level;
	unsigned char *sft_veto;

	/* dc factor to apply to further loaded SFTs */
	float dc_factor;
	int dc_factor_touched;
	int dc_factor_blocked;

	/* real and imaginary parts - they are separate to facilitate use of vectorized operations */
	float *re; 
	float *im;

	int size;
	int free;	/* this used to be called nsegments */

	int offset;	/* this points to our place in global coefficient tables */
	int first_bin;  /* first bin of each segment */
	int nbins;	/* number of bins in each segment */
	double coherence_time; /* how long the SFTs are */

	/* statistics and data derived from this data set */
	double *mean;
	double *weighted_mean;
	double *sigma;
	double *bin_weight;

	double *new_mean;
	double *new_weighted_mean;
	double *new_sigma;
	double *new_bin_weight;

	LINES_REPORT *lines_report;

	EarthState *earth_state;
	float *detector_velocity;
	double average_detector_velocity[3];
	double band_axis[3];
	double band_axis_norm;
	double large_S;

	int input_flags;
	#define FLAG_BAND_AXIS		1
	#define FLAG_BAND_AXIS_NORM	2
	#define FLAG_LARGE_S		4
	double user_band_axis[3];
	double user_band_axis_norm;
	double user_large_S;

	double weight; /* additional, user-specified weight factor */

	/* results of noise decomposition - TMedian and LogMedians are log10(power) */
	float *TMedians;
	float *FMedians;

	float *expTMedians; /* squared relative exponentiated TMedians */
	float *expTMedians_plain; /* plain exponentiated TMedians */

	float *expFMedians_plain; /* plain exponentiated TMedians */

	float TMedian;
	float expTMedian;

	/* helper arrays, for plotting */
	float *hours;
	float *frequencies;
	double *hours_d;
	double *freq_d;

	SKY_GRID_TYPE *AM_coeffs_plus;
	SKY_GRID_TYPE *AM_coeffs_cross;
	long AM_coeffs_size;

	} DATASET;

void output_dataset_info(DATASET *d);
void characterize_dataset(DATASET *d);
void compute_noise_curves(DATASET *dataset);

void load_dataset_from_file(char *file);
long total_segments(void);
long vetoed_segments(void);
float datasets_normalizing_weight(void);
INT64 min_gps(void);
INT64 max_gps(void);
void post_init_datasets(void);
void datasets_average_detector_speed(double *average_det_velocity);
float effective_weight_ratio(float target_ra, float target_dec, float source_ra, float source_dec, float source_spindown, float bin_tolerance, float spindown_tolerance);
float stationary_effective_weight_ratio(float target_ra, float target_dec, float bin_tolerance);
void dump_datasets(char *filename);
void sftv2_dump_datasets(char *directory);
void output_datasets_info(void);
void verify_dataset_whole_sky_AM_response(void);
void test_datasets(void);
void init_datasets(void);

#define PHASE_ACCUMULATION_TIMEOUT 50000

#endif 
