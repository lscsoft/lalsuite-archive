#ifndef __POWER_SUMS_H__
#define __POWER_SUMS_H__

#include "power_cache.h"

typedef struct {
	/* These are double so we can follow parameters of observed source, e.g. Sco-X1
	   It is possible that depth could be a float */
	double freq_modulation_freq;
	double freq_modulation_phase;
	double freq_modulation_depth;
	
	float spindown;
	float fdotdot;
	float ra;
	float dec;
	
	int skyband;
	int snr_subbin;
	
	/* These are auxiliary values for diverted templates */
	float snr;
	float ul;
	float circ_ul;
	
	int first_chunk;
	int last_chunk;
	int veto_num;
	int pi;   /* this is patch index, but could also be arbitrary grouping parameter. It is useful to group templates to increase cache efficiency in followup stages */
	} TEMPLATE_INFO;

typedef struct S_POWER_SUM {
	float freq_shift; /* additional shift e.g. for half-bin sampling */
	float spindown;
	
	/* These are double so we can follow parameters of observed source, e.g. Sco-X1
	   It is possible that depth could be a float */
	double freq_modulation_freq;
	double freq_modulation_phase;
	double freq_modulation_depth;
	
	float fdotdot;
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

#include "summing_context.h"
	
void generate_patch_templates(SUMMING_CONTEXT *ctx, int pi, POWER_SUM **ps, int *count);
void generate_followup_templates(SUMMING_CONTEXT *ctx, TEMPLATE_INFO *template_info, int ti_count, POWER_SUM **ps, int *count);
void clone_templates(SUMMING_CONTEXT *ctx, POWER_SUM *ps, int count, POWER_SUM **ps_out);
//void free_templates(POWER_SUM *ps, int count);
void free_templates_ctx(SUMMING_CONTEXT *ctx, POWER_SUM *ps, int count);

void accumulate_power_sums_sidereal_step(SUMMING_CONTEXT *ctx, POWER_SUM *ps, int count, double gps_start, double gps_stop, int veto_mask);
void accumulate_power_sums_plain(SUMMING_CONTEXT *ctx, POWER_SUM *ps, int count, double gps_start, double gps_stop, int veto_mask);

#endif
