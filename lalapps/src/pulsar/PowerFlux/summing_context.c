#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "global.h"
#include "power_cache.h"
#include "power_sums.h"
#include "summing_context.h"
#include "cmdline.h"

extern struct gengetopt_args_info args_info;

extern FILE * DATA_LOG, * LOG, *FILE_LOG;


SUMMING_CONTEXT *create_summing_context(void)
{
SUMMING_CONTEXT *ctx;

ctx=do_alloc(1, sizeof(*ctx));

fprintf(stderr, "Averaging mode: %s\n", args_info.averaging_mode_arg);

/* default values appropriate for particular averaging mode */
if(!strcasecmp(args_info.averaging_mode_arg, "matched")) {
	ctx->get_uncached_power_sum=get_uncached_matched_power_sum;
	ctx->accumulate_power_sum_cached=accumulate_power_sum_cached1;

	ctx->cache_granularity=10;
	} else
if(!strcasecmp(args_info.averaging_mode_arg, "3") || !strcasecmp(args_info.averaging_mode_arg, "three")) {
	fprintf(stderr, "PowerFlux2 does not support 3-bin mode\n");
	exit(-1);
	} else
if(!strcasecmp(args_info.averaging_mode_arg, "1") || !strcasecmp(args_info.averaging_mode_arg, "one")) {
	ctx->get_uncached_power_sum=sse_get_uncached_single_bin_power_sum;
	ctx->accumulate_power_sum_cached=accumulate_single_bin_power_sum_cached1;

	ctx->cache_granularity=1;
	} else {
	fprintf(stderr, "Unrecognized averaging mode: %s\n", args_info.averaging_mode_arg);
	exit(-1);
	}

/* apply command line overrides */
if(args_info.cache_granularity_given) {
	ctx->cache_granularity=args_info.cache_granularity_arg;
	}

ctx->inv_cache_granularity=1.0/ctx->cache_granularity;
ctx->half_inv_cache_granularity=0.5/ctx->cache_granularity;
fprintf(LOG, "cache_granularity: %f\n", ctx->cache_granularity);

return(ctx);
}

void free_summing_context(SUMMING_CONTEXT *ctx)
{
free(ctx);
}