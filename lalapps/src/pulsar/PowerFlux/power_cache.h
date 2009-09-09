#ifndef __POWER_CACHE_H__
#define __POWER_CACHE_H__

typedef struct {
	/* fields below are filled in when computing power */

	/* amplitude response factors to use in power sum - note these are kept constant throughout the patch */
	float f_plus;
	float f_cross;

	/* bin shift to apply, this is in units of 1/coherence_time - as opposed to power sums */
	double bin_shift;

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

#include "summing_context.h"

void accumulate_partial_power_sum_F1(PARTIAL_POWER_SUM_F *accum, PARTIAL_POWER_SUM *partial);
void accumulate_partial_power_sum_F2(PARTIAL_POWER_SUM *accum, PARTIAL_POWER_SUM_F *partial);

SEGMENT_INFO *find_segments(double gps_start, double gps_end, int veto_mask, int *count);

void allocate_simple_cache(SUMMING_CONTEXT *ctx);

void get_uncached_single_bin_power_sum(SUMMING_CONTEXT *ctx, SEGMENT_INFO *si, int count, PARTIAL_POWER_SUM_F *pps);
void sse_get_uncached_single_bin_power_sum(SUMMING_CONTEXT *ctx, SEGMENT_INFO *si, int count, PARTIAL_POWER_SUM_F *pps);

void get_uncached_matched_power_sum(SUMMING_CONTEXT *ctx, SEGMENT_INFO *si, int count, PARTIAL_POWER_SUM_F *pps);
void sse_get_uncached_matched_power_sum(SUMMING_CONTEXT *ctx, SEGMENT_INFO *si, int count, PARTIAL_POWER_SUM_F *pps);
void accumulate_power_sum_cached1(SUMMING_CONTEXT *ctx, SEGMENT_INFO *si, int count, PARTIAL_POWER_SUM_F *pps);

void power_cache_selftest(void);

#endif
