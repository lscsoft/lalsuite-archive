#ifndef __POWER_SUMS_H__
#define __POWER_SUMS_H__

#include "power_cache.h"

typedef struct {
	double freq_shift; /* additional shift e.g. for half-bin sampling */
	double spindown;
	double ra;
	double dec;
	double patch_ra;
	double patch_dec;

	double min_gps;
	double max_gps;

	float e[26];
	float patch_e[26];

	int skyband;

	PARTIAL_POWER_SUM *pps;
	} POWER_SUM;

void generate_patch_templates(int pi, POWER_SUM **ps, int *count);
void clone_templates(POWER_SUM *ps, int count, POWER_SUM **ps_out);
void free_templates(POWER_SUM *ps, int count);

void accumulate_power_sums(POWER_SUM *ps, int count, double gps_start, double gps_stop, int veto_mask);

#endif
