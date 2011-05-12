#ifndef __CONTEXT_H__
#define __CONTEXT_H__

#define restrict

#include <lal/LALDatatypes.h>
#include <lal/AVFactories.h>
#include <lal/ComplexFFT.h>
#include <lal/DetectorSite.h>
#include <lal/LALBarycenter.h>

typedef struct {
	int free;
	int size;
	
	int *bin;
	COMPLEX16 *data;
	COMPLEX16 first9[9];
	
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
	
	/* Number of sub-bin frequency steps */
	int n_fsteps;

	int half_window;
	
	COMPLEX16Vector *plus_te_fft;
	COMPLEX16Vector *cross_te_fft;	
	
	/* scan_ffts */
	double *power;
	
	/* compute_*_offset */
	int offset_count;
	COMPLEX16Vector *offset_in;
	COMPLEX16Vector *offset_fft;	
	COMPLEX16FFTPlan *offset_fft_plan;
	
	/* fast_get_emission_time */
	ETC etc;
	
	COMPLEX16Vector *scan_tmp[4];
	
	} LOOSE_CONTEXT;

LOOSE_CONTEXT * create_context(void);

#endif