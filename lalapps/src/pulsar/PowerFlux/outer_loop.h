#ifndef __OUTER_LOOP_H__
#define __OUTER_LOOP_H__

#include "power_sum_stats.h"
#include "jobs.h"

typedef ALIGN_DECLSPEC struct {
	MUTEX mutex;
	char *name;

	float *ul_skymap;
	float *circ_ul_skymap;
	float *snr_skymap;
	float *ul_freq_skymap;
	float *circ_ul_freq_skymap;
	float *snr_freq_skymap;
	float *snr_ul_skymap;
	float *max_weight_skymap;
	float *min_weight_skymap;
	float *weight_loss_fraction_skymap;
	float *ks_skymap;

	POWER_SUM_STATS *band_info;
	int *band_valid_count;
	int *band_masked_count;

	/* convenience info for keeping track of which ei is which */
	int first_chunk;
	int last_chunk;
	int veto_num;
	} EXTREME_INFO;


void outer_loop(void);

#endif
