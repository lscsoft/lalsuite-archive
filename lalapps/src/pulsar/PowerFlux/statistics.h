#ifndef __STATISTICS_H__
#define __STATISTICS_H__

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

#endif
