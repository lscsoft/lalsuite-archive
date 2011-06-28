#ifndef __CONTEXT_H__
#define __CONTEXT_H__

#define restrict

#include <lal/LALDatatypes.h>
#include <lal/AVFactories.h>
#include <lal/ComplexFFT.h>
#include <lal/DetectorSite.h>
#include <lal/LALBarycenter.h>

#include "fft_stats.h"

typedef struct {
	int free;
	int size;
	
	int *bin;
	COMPLEX8 *data;
	COMPLEX8 first9[9];
	
	double slope;
	} SPARSE_CONV;

typedef struct {
	double ra;
	double dec;
	double dInv;
	LIGOTimeGPS tGPS1;
	LIGOTimeGPS tGPS2;
	char *detector;
	EmissionTime et1;
	EmissionTime et2;
	double range;
	}  ETC;

typedef struct {
	int nsamples; /* total number of samples in adjusted fft */
	double timebase;
	double first_gps;
	
	/* signal parameters */
	double frequency;
	double spindown;
	double ra;
	double dec;
	double dInv;
	int fstep;
	
	float *power;
	float *cum_power;
	float *variance;
	
	COMPLEX16FFTPlan *fft_plan;

	/* Demodulation code data */
	COMPLEX16Vector *plus_samples;
	COMPLEX16Vector *cross_samples;

	double weight_pp;
	double weight_pc;
	double weight_cc;

	COMPLEX16Vector *plus_fft;
	COMPLEX16Vector *cross_fft;	

	/* Bessel coeffs */
	
	SPARSE_CONV *te_sc;
	SPARSE_CONV *spindown_sc;
	SPARSE_CONV *ra_sc;
	SPARSE_CONV *dec_sc;
	
	/* frequency adjustment */
	int n_freq_adj_filter;
	int n_scan_fft_filter;
	
	/* Number of sub-bin frequency steps */
	int n_fsteps;
	int n_sky_scan;

	int half_window;
	
	int variance_half_window;
	
	int total_segments;
	
	COMPLEX8Vector *plus_te_fft;
	COMPLEX8Vector *cross_te_fft;	
	
	/* compute_*_offset */
	int offset_count;
	COMPLEX16Vector *offset_in;
	COMPLEX16Vector *offset_fft;	
	COMPLEX16FFTPlan *offset_fft_plan;
	
	/* fast_get_emission_time */
	ETC etc;
	
	COMPLEX8Vector *scan_tmp[8];
	
	FFT_STATS stats;
	
	} LOOSE_CONTEXT;

LOOSE_CONTEXT * create_context(void);

#endif