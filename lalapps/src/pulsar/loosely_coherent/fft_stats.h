#ifndef __FFT_STATS_H__
#define __FFT_STATS_H__

#include <lal/LALDatatypes.h>

typedef struct {
	double iota;
	double psi;
	double Ap;
	double Ax;
	double w1_re;
	double w1_im;
	double w2_re;
	double w2_im;
	double w11;
	double w12;
	double w22;
	} ALIGNMENT_COEFFS;

typedef struct {
	int free;
	int size;
	
	ALIGNMENT_COEFFS *coeffs;
	
	} ALIGNMENT_DATA;

	
typedef struct {
	double value;
	int fft_bin;
	int alignment_bin;
	double fft_offset;
	double frequency;
	double spindown;
	double ra;
	double dec;
	double iota;
	double psi;
	double phi;
	COMPLEX16 z;
	} STAT_INFO;
	
typedef struct {
	STAT_INFO snr;
	STAT_INFO ul;
	STAT_INFO circ_ul;
	STAT_INFO B_stat;
	STAT_INFO F_stat;
	STAT_INFO power_stat;
	
	double template_count;
	double stat_hit_count;
	double min_noise_ratio;
	double max_noise_ratio;
	} FFT_STATS;

	
#include "context.h"	
	
void init_fft_stats(void);
void init_stats(FFT_STATS *st);
void log_stats(LOOSE_CONTEXT *ctx, FILE *f, char *tag, FFT_STATS *st, double ul_adjust);
void update_stats(FFT_STATS *st_accum, FFT_STATS *st);
void compute_stats_variance(LOOSE_CONTEXT *ctx);
 void compute_fft_stats(LOOSE_CONTEXT *ctx, FFT_STATS *stats, COMPLEX8Vector *fft1, COMPLEX8Vector *fft2, double fft_offset);


#endif