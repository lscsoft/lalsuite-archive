#ifndef __SINGLE_BIN_LOOSELY_COHERENT_SUM_H__
#define __SINGLE_BIN_LOOSELY_COHERENT_SUM_H__


#include "summing_context.h"
#include "power_cache.h"

void get_uncached_loose_partial_power_sum(SUMMING_CONTEXT *ctx, SEGMENT_INFO *si, int count, PARTIAL_POWER_SUM_F *pps);
void get_uncached_loose_single_bin_partial_power_sum(SUMMING_CONTEXT *ctx, SEGMENT_INFO *si, int count, PARTIAL_POWER_SUM_F *pps);
void accumulate_power_sum_cached_diff(SUMMING_CONTEXT *ctx, SEGMENT_INFO *si, int count, PARTIAL_POWER_SUM_F *pps);

#endif