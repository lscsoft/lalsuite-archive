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
	long ks_count;
	} NORMAL_STATS;



void compute_normal_sorted_stats(float *data, long count, NORMAL_STATS *stats);

void compute_normal_stats(float *data, long count, NORMAL_STATS *stats);

void sort_floats(float *data, long count);

typedef struct {
	int nbands;
	int nbins;  /* these are histogram bins.. */
	STAT_TYPE *max;
	STAT_TYPE *min;
	long *hist; /* counts */
	} HISTOGRAM;

void init_statistics(void);

HISTOGRAM * new_histogram(int nbins, int nbands);
void free_histogram(HISTOGRAM *h);
void compute_histogram_f(HISTOGRAM *h, float *data, int *bands, long count);
void print_histogram(FILE *f, HISTOGRAM *h, char *prefix);

#endif
