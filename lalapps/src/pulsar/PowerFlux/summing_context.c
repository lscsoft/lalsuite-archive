#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "global.h"
#include "power_cache.h"
#include "power_sums.h"
#include "summing_context.h"
#include "single_bin_loosely_coherent_sum.h"
#include "matched_loosely_coherent_sum.h"
#include "cmdline.h"

extern struct gengetopt_args_info args_info;

extern FILE * DATA_LOG, * LOG, *FILE_LOG;


SUMMING_CONTEXT *create_summing_context(void)
{
SUMMING_CONTEXT *ctx;

ctx=do_alloc(1, sizeof(*ctx));
memset(ctx, 0, sizeof(*ctx));

fprintf(stderr, "Averaging mode: %s\n", args_info.averaging_mode_arg);
fprintf(stderr, "SSE: %d\n", args_info.sse_arg);
fprintf(LOG, "Averaging mode: %s\n", args_info.averaging_mode_arg);
fprintf(LOG, "SSE: %d\n", args_info.sse_arg);

ctx->diff_shift_granularity=0; 

#define MODE(a)	(args_info.sse_arg ? (sse_ ## a) : (a) )

/* default values appropriate for particular averaging mode */
if(!strcasecmp(args_info.averaging_mode_arg, "matched")) {
	ctx->get_uncached_power_sum=MODE(get_uncached_matched_power_sum);
	ctx->accumulate_power_sum_cached=accumulate_power_sum_cached1;
	ctx->accumulate_power_sums=accumulate_power_sums_sidereal_step;

	ctx->cache_granularity=8;
	ctx->sidereal_group_count=24;
	ctx->summing_step=864000; /* ten days */
	} else
if(!strcasecmp(args_info.averaging_mode_arg, "single_bin_loose")) {
	fprintf(LOG, "single_bin_loose: only tested with delta of pi/2 and pi/5\n");
	fprintf(stderr, "single_bin_loose: only tested with delta of pi/2 and pi/5\n");
	ctx->get_uncached_power_sum=MODE(get_uncached_loose_single_bin_partial_power_sum);
	ctx->accumulate_power_sum_cached=accumulate_power_sum_cached_diff;
	ctx->accumulate_power_sums=accumulate_single_bin_loose_power_sums_sidereal_step;

	ctx->cache_granularity=8; /* TODO: find actual value from experiment */
	ctx->diff_shift_granularity=8192*8; 
	ctx->sidereal_group_count=12;
	ctx->summing_step=86400*3; /* three days */
	ctx->time_group_count=3;
	ctx->cross_terms_present=args_info.compute_cross_terms_arg;
	} else
if(!strcasecmp(args_info.averaging_mode_arg, "matched_loose")) {
	fprintf(LOG, "**** WARNING matched_loose: this experimental code has not been reviewed yet.\n");
	fprintf(stderr, "**** WARNING matched_loose: this experimental code has not been reviewed yet.\n");
	ctx->get_uncached_power_sum=get_uncached_loose_matched_partial_power_sum;
	ctx->accumulate_power_sum_cached=accumulate_power_sum_cached_diff;
	ctx->accumulate_power_sums=accumulate_matched_loose_power_sums_sidereal_step;

	ctx->cache_granularity=8; /* TODO: find actual value from experiment */
	ctx->diff_shift_granularity=8192*8; 
	ctx->sidereal_group_count=12;
	ctx->summing_step=86400*3; /* three days */
	ctx->time_group_count=3;
	ctx->cross_terms_present=args_info.compute_cross_terms_arg;
	} else
if(!strcasecmp(args_info.averaging_mode_arg, "3") || !strcasecmp(args_info.averaging_mode_arg, "three")) {
	fprintf(stderr, "PowerFlux2 does not support 3-bin mode\n");
	exit(-1);
	} else
if(!strcasecmp(args_info.averaging_mode_arg, "1") || !strcasecmp(args_info.averaging_mode_arg, "one")) {
	ctx->get_uncached_power_sum=MODE(get_uncached_single_bin_power_sum);
	ctx->accumulate_power_sum_cached=accumulate_power_sum_cached1;
	ctx->accumulate_power_sums=accumulate_power_sums_sidereal_step;

	ctx->cache_granularity=1;
	ctx->sidereal_group_count=24;
	ctx->summing_step=864000; /* ten days */
	} else {
	fprintf(stderr, "Unrecognized averaging mode: %s\n", args_info.averaging_mode_arg);
	exit(-1);
	}

/* apply command line overrides */
if(args_info.cache_granularity_given) {
	ctx->cache_granularity=args_info.cache_granularity_arg;
	if(ctx->cache_granularity<1) {
		fprintf(stderr, "*** ERROR: cache granularity must be positive\n");
		fprintf(LOG, "*** ERROR: cache granularity must be positive\n");
		exit(-1);
		}
	}

if(args_info.summing_step_given) {
	ctx->summing_step=args_info.summing_step_arg;
	}

if(args_info.sidereal_group_count_given) {
	ctx->sidereal_group_count=args_info.sidereal_group_count_arg;
	if(ctx->sidereal_group_count<1) {
		fprintf(stderr, "*** ERROR: sidereal_group_count must be positive\n");
		fprintf(LOG, "*** ERROR: sidereal_group_count must be positive\n");
		exit(-1);
		}
	}

if(args_info.time_group_count_given) {
	ctx->time_group_count=args_info.time_group_count_arg;
	if(ctx->time_group_count<1) {
		fprintf(stderr, "*** ERROR: time_group_count must be positive\n");
		fprintf(LOG, "*** ERROR: time_group_count must be positive\n");
		exit(-1);
		}
	}
	
ctx->inv_cache_granularity=1.0/ctx->cache_granularity;
ctx->half_inv_cache_granularity=0.5/ctx->cache_granularity;
ctx->inv_diff_shift_granularity=1.0/ctx->diff_shift_granularity;
ctx->half_inv_diff_shift_granularity=0.5/ctx->diff_shift_granularity;

fprintf(LOG, "summing_step: %g\n", ctx->summing_step);
fprintf(LOG, "cache_granularity: %d\n", ctx->cache_granularity);
fprintf(LOG, "diff_shift_granularity: %d\n", ctx->diff_shift_granularity);
fprintf(LOG, "sidereal_group_count: %d\n", ctx->sidereal_group_count);
fprintf(LOG, "time_group_count: %d\n", ctx->time_group_count);
fprintf(LOG, "phase_mismatch: %g\n", args_info.phase_mismatch_arg);
fprintf(LOG, "cross_terms_present: %d\n", ctx->cross_terms_present);

allocate_simple_cache(ctx);

ctx->loose_first_half_count=-1;


return(ctx);
}

void free_summing_context(SUMMING_CONTEXT *ctx)
{
if(ctx->free_cache!=NULL)ctx->free_cache(ctx);
free(ctx);
}
