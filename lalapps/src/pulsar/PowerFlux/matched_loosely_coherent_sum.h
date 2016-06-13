#ifndef __MATCHED_LOOSELY_COHERENT_SUM_H__
#define __MATCHED_LOOSELY_COHERENT_SUM_H__


#include "summing_context.h"
#include "power_cache.h"

void get_uncached_loose_matched_partial_power_sum(SUMMING_CONTEXT *ctx, SEGMENT_INFO *si, int count, PARTIAL_POWER_SUM_F *pps);
/* This function is meant to work with get_uncached_loose_matched_partial_power_sum */
void accumulate_matched_loose_power_sums_sidereal_step(SUMMING_CONTEXT *ctx, POWER_SUM *ps, int count, double gps_start, double gps_stop, int veto_mask);

#endif