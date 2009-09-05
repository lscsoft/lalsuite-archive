#ifndef __SUMMING_CONTEXT_H__
#define __SUMMING_CONTEXT_H__

#include "power_cache.h"

typedef struct S_SUMMING_CONTEXT {
	void (*get_uncached_power_sum)(struct S_SUMMING_CONTEXT  *ctx, SEGMENT_INFO *si, int count, PARTIAL_POWER_SUM_F *pps);
	void (*accumulate_power_sum_cached)(struct S_SUMMING_CONTEXT  *ctx, SEGMENT_INFO *si, int count, PARTIAL_POWER_SUM_F *pps);

	int cache_granularity;
	float inv_cache_granularity;
	float half_inv_cache_granularity;

	} SUMMING_CONTEXT;

SUMMING_CONTEXT *create_summing_context(void);
void free_summing_context(SUMMING_CONTEXT *ctx);

#endif