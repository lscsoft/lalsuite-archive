#ifndef __CONTEXT_H__
#define __CONTEXT_H__

#define restrict

#include <lal/LALDatatypes.h>
#include <lal/AVFactories.h>
#include <lal/ComplexFFT.h>

typedef struct {
	int nsamples; /* total number of samples in adjusted fft */
	
	COMPLEX16FFTPlan *fft_plan;

	/* Demodulation code data */
	COMPLEX16Vector *plus_samples;
	COMPLEX16Vector *cross_samples;

	COMPLEX16Vector *plus_fft;
	COMPLEX16Vector *cross_fft;	

	COMPLEX16Vector *plus_te_fft;
	COMPLEX16Vector *cross_te_fft;	
	} LOOSE_CONTEXT;

LOOSE_CONTEXT * create_context(void);

#endif