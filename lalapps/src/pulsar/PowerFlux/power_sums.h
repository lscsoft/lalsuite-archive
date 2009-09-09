#ifndef __POWER_SUMS_H__
#define __POWER_SUMS_H__

#include "power_cache.h"
#include "summing_context.h"

typedef struct {
	float freq_shift; /* additional shift e.g. for half-bin sampling */
	float spindown;
	float ra;
	float dec;
	float patch_ra;
	float patch_dec;

	double min_gps;
	double max_gps;

	float e[26];
	float patch_e[26];

	int skyband;

	PARTIAL_POWER_SUM_F *pps;
	} POWER_SUM;

void generate_patch_templates(int pi, POWER_SUM **ps, int *count);
void clone_templates(POWER_SUM *ps, int count, POWER_SUM **ps_out);
void free_templates(POWER_SUM *ps, int count);

void accumulate_power_sums(SUMMING_CONTEXT *ctx, POWER_SUM *ps, int count, double gps_start, double gps_stop, int veto_mask);

#endif
