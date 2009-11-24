#ifndef __SUMMING_CONTEXT_H__
#define __SUMMING_CONTEXT_H__

#include "power_cache.h"

typedef struct S_SUMMING_CONTEXT {
	void (*get_uncached_power_sum)(struct S_SUMMING_CONTEXT  *ctx, SEGMENT_INFO *si, int count, PARTIAL_POWER_SUM_F *pps);
	void (*accumulate_power_sum_cached)(struct S_SUMMING_CONTEXT  *ctx, SEGMENT_INFO *si, int count, PARTIAL_POWER_SUM_F *pps);
	void (*accumulate_power_sums)(struct S_SUMMING_CONTEXT *ctx, struct S_POWER_SUM *ps, int count, double gps_start, double gps_stop, int veto_mask);

	int cache_granularity;
	float inv_cache_granularity;
	float half_inv_cache_granularity;
	int sidereal_group_count; /* group sfts falling on similar times of the day in this many groups */
	double summing_step; /* process SFTs in blocks of this many seconds each */
	int time_group_count; /* group SFTs by their GPS time within a block into this many groups - used by loosely coherent code */
	float loose_coherence_alpha;

	void *cache;
	void (*free_cache)(struct S_SUMMING_CONTEXT *ctx);
	void (*print_cache_stats)(struct S_SUMMING_CONTEXT *ctx);
	void (*reset_cache)(struct S_SUMMING_CONTEXT *ctx, int segment_count, int template_count);

	/* dynamic parameters */
	int loose_first_half_count;
	} SUMMING_CONTEXT;

SUMMING_CONTEXT *create_summing_context(void);
void free_summing_context(SUMMING_CONTEXT *ctx);

#endif