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
extern int useful_bins;

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

#if MANUAL_SSE
#define MODE(a)	(args_info.sse_arg ? (sse_ ## a) : (a) )
#else
#define MODE(a)	(a)
#endif

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
	
	if(args_info.freq_modulation_depth_arg!=0 || args_info.freq_modulation_depth_count_arg>1) {
		fprintf(stderr, "*** WARNING: single-bin mode uses approximation to binary evolution equation that usually requires scanning over freq_modulation_phase\n");
		}
	
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

if(args_info.diff_shift_granularity_given) {
	ctx->diff_shift_granularity=args_info.diff_shift_granularity_arg;
	if(ctx->diff_shift_granularity<1) {
		fprintf(stderr, "*** ERROR: diff shift granularity must be positive\n");
		fprintf(LOG, "*** ERROR: diff shift granularity must be positive\n");
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

ctx->pps_pool_size=20000;
ctx->pps_pool_free=0;
ctx->partial_power_sum_pool=do_alloc(ctx->pps_pool_size, sizeof(*ctx->partial_power_sum_pool));
ctx->pps_hits=0;
ctx->pps_misses=0;
ctx->pps_rollbacks=0;

ctx->log_extremes_pps=allocate_partial_power_sum_F(useful_bins, 1);
ctx->log_extremes_pstats_scratch_size=10;
ctx->log_extremes_pstats_scratch=do_alloc(ctx->log_extremes_pstats_scratch_size, sizeof(*ctx->log_extremes_pstats_scratch));

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
int i;

fprintf(stderr, "pps_hits=%ld pps_misses=%ld pps_rollbacks=%ld pool_size=%ld\n", ctx->pps_hits, ctx->pps_misses, ctx->pps_rollbacks, ctx->pps_pool_size);

if(ctx->free_cache!=NULL)ctx->free_cache(ctx);

for(i=0;i<ctx->pps_pool_free;i++)free_partial_power_sum_F(ctx->partial_power_sum_pool[i]);
free(ctx->partial_power_sum_pool);

free_partial_power_sum_F(ctx->log_extremes_pps);

free(ctx);
}

PARTIAL_POWER_SUM_F * get_partial_power_sum_F(SUMMING_CONTEXT *ctx, int pps_bins, int cross_terms_present)
{
PARTIAL_POWER_SUM_F * p;
while(ctx->pps_pool_free>0) {
	ctx->pps_pool_free--;
	
	p=ctx->partial_power_sum_pool[ctx->pps_pool_free];
	ctx->partial_power_sum_pool[ctx->pps_pool_free]=NULL;
	
	if(p->nbins!=pps_bins) {
		free_partial_power_sum_F(p);
		ctx->pps_rollbacks++;
		continue;
		}
		
	if(p->power_im_pc==NULL && cross_terms_present) {
		free_partial_power_sum_F(p);
		ctx->pps_rollbacks++;
		continue;
		}
	ctx->pps_hits++;
	return(p);
	}
ctx->pps_misses++;
return(allocate_partial_power_sum_F(pps_bins, cross_terms_present));
}

void put_partial_power_sum_F(SUMMING_CONTEXT *ctx, PARTIAL_POWER_SUM_F *pps)
{
PARTIAL_POWER_SUM_F ** p;

if(ctx->pps_pool_free>=ctx->pps_pool_size) {
	ctx->pps_pool_size*=2;
	p=do_alloc(ctx->pps_pool_size, sizeof(*ctx->partial_power_sum_pool));
	memcpy(p, ctx->partial_power_sum_pool, ctx->pps_pool_free*sizeof(*ctx->partial_power_sum_pool));
	free(ctx->partial_power_sum_pool);
	ctx->partial_power_sum_pool=p;
	}
ctx->partial_power_sum_pool[ctx->pps_pool_free]=pps;
ctx->pps_pool_free++;
}

