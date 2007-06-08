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

#ifndef __STATISTICS_H__
#define __STATISTICS_H__

#include <stdio.h>

typedef double STAT_TYPE;

#define STAT_FLAG_SORT_DATA		(1<<0)
#define STAT_FLAG_INPLACE_SORT_DATA	(1<<1)
#define STAT_FLAG_ESTIMATE_MEAN		(1<<2)
#define STAT_FLAG_ESTIMATE_SIGMA	(1<<3)
#define STAT_FLAG_ESTIMATE_KS_LEVEL	(1<<4)
#define STAT_FLAG_COMPUTE_KS_TEST	(1<<5)
#define STAT_FLAG_AUTOFIT		(1<<6)

typedef struct {
	int flag;
	STAT_TYPE ks_level;
	/* estimated distribution parameters */
	STAT_TYPE  mean;
	STAT_TYPE  sigma;
	/* statistics results */
	STAT_TYPE  ks_test;
	int ks_count;
	} NORMAL_STATS;



void compute_normal_sorted_stats(float *data, int count, NORMAL_STATS *stats);

void compute_normal_stats(float *data, int count, NORMAL_STATS *stats);

void sort_floats(float *data, int count);

typedef struct {
	int nbands;
	int nbins;  /* these are histogram bins.. */
	STAT_TYPE *max;
	STAT_TYPE *min;
	int *hist; /* counts */
	} HISTOGRAM;

void init_statistics(void);

HISTOGRAM * new_histogram(int nbins, int nbands);
void free_histogram(HISTOGRAM *h);
void compute_histogram_f(HISTOGRAM *h, float *data, int *bands, int count);
void print_histogram(FILE *f, HISTOGRAM *h, char *prefix);

#endif
