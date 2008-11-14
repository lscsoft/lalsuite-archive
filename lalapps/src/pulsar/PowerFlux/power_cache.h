#ifndef __POWER_CACHE_H__
#define __POWER_CACHE_H__

typedef struct {
	/* fields below are filled in when computing power */

	/* sky coordinates for amplitude response */
	double ra;
	double dec;
	/* for convenience and to reduce number of times we call sin/cos */
	float e[26];

	/* frequency shift to apply */
	double freq_shift;

	/* fields below are filled in when locating segments */

	/* for convenience, same as datasets[dataset].gps[segment] */
	double gps;
	double detector_velocity[3]; /* also from datasets[dataset] */
	double coherence_time;


	/* segment coordinates */
	int dataset;
	int segment;
	} SEGMENT_INFO;

#define REAL double
#define SUFFIX(a) a
#include "partial_power_sum.h"

#undef REAL
#undef SUFFIX

#define REAL float
#define SUFFIX(a) a##_F
#include "partial_power_sum.h"

#undef REAL
#undef SUFFIX

void accumulate_partial_power_sum_F1(PARTIAL_POWER_SUM_F *accum, PARTIAL_POWER_SUM *partial);
void accumulate_partial_power_sum_F2(PARTIAL_POWER_SUM *accum, PARTIAL_POWER_SUM_F *partial);

SEGMENT_INFO *find_segments(double gps_start, double gps_end, int veto_mask, int *count);

void get_uncached_single_bin_power_sum(SEGMENT_INFO *si, int count, PARTIAL_POWER_SUM *pps);
void accumulate_single_bin_power_sum_cached1(SEGMENT_INFO *si, int count, PARTIAL_POWER_SUM *pps);

void print_cache_stats(void);

#endif
