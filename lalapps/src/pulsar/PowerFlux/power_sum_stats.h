#ifndef __POWER_SUM_STATS_H__
#define __POWER_SUM_STATS_H__

#include "power_sums.h"

typedef struct {
	float iota;
	float psi;
	
	float pp;
	float pc;
	float cc;
	float im_pc;

	float pppp;
	float pppc;
	float ppcc;
	float pccc;
	float cccc;
	float im_ppcc;
	} ALIGNMENT_COEFFS;

typedef struct {
	double iota;
	double psi;

	double ul;
	double ll;
	double centroid;
	double snr;

	double M;
	double S;
	double ks_value;
	double m1_neg;
	double m3_neg;
	double m4;
	double max_weight;
	double weight_loss_fraction;

	int ks_count;

	int bin;

	/* the following fields are for convenience and are filled in by outside code based on value of bin */
	double frequency;
	double spindown;
	double ra;
	double dec;

	} POINT_STATS;

typedef struct {
	POINT_STATS highest_ul;
	POINT_STATS highest_snr;
	POINT_STATS highest_ks;
	POINT_STATS highest_M;
	POINT_STATS highest_S;
	POINT_STATS highest_circ_ul;
	double max_weight_loss_fraction;
	double max_weight;
	double min_weight;
	double max_m1_neg;
	double min_m1_neg;
	double max_m3_neg;
	double min_m3_neg;
	double max_m4;
	double min_m4;
	int ntemplates;
	} POWER_SUM_STATS;

void compute_alignment_coeffs(ALIGNMENT_COEFFS *ac);
void generate_alignment_grid(void);

//void point_power_sum_stats(PARTIAL_POWER_SUM_F *pps, ALIGNMENT_COEFFS *ag, POINT_STATS *pst);
void power_sum_stats(PARTIAL_POWER_SUM_F *pps, POWER_SUM_STATS *stats);

void init_power_sum_stats(void);
void power_sum_stats_selftest(void);

#endif
