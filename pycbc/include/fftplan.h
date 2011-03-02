#include <pycbc.h>
#include <cuda.h>
#include <cufft.h>

#ifndef PYCBC_FFTPLAN_H
#define PYCBC_FFTPLAN_H

real_fft_plan_t* new_real_fft_plan_t( int size, int fwdflag, int level );
void delete_real_fft_plan_t( real_fft_plan_t *plan );

complex_fft_plan_t* new_complex_fft_plan_t( int size, int fwdflag, int level );
void delete_complex_fft_plan_t( complex_fft_plan_t *plan );

#endif /* PYCBC_FFTPLAN_H */
