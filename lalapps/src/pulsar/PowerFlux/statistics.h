#ifndef __STATISTICS_H__
#define __STATISTICS_H__

typedef double STAT_TYPE;

typedef struct {
	/* estimated distribution parameters */
	STAT_TYPE  mean;
	STAT_TYPE  sigma;
	/* statistics results */
	STAT_TYPE  ks_test;
	long ks_count;
	} NORMAL_STATS;



void compute_normal_sorted_ks_test(float *data, long count, NORMAL_STATS *stats);


void compute_normal_ks_test(float *data, long count, NORMAL_STATS *stats);


void sort_floats(float *data, long count);

#endif
