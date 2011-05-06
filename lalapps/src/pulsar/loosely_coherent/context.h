#ifndef __CONTEXT_H__
#define __CONTEXT_H__

#define restrict

#include <lal/LALDatatypes.h>
#include <lal/AVFactories.h>
#include <lal/ComplexFFT.h>

typedef struct {
	int free;
	int size;
	
	int *bin;
	COMPLEX16 *data;
	COMPLEX16 first9[9];
	
	double slope;
	} SPARSE_CONV;

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
	
	} LOOSE_CONTEXT;

LOOSE_CONTEXT * create_context(void);

#endif