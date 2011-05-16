#ifndef __BESSEL_H__
#define __BESSEL_H__

#include <lal/LALDatatypes.h>

void make_bessel_filter(COMPLEX8 *filter, int filter_size, COMPLEX8 *coeffs, int coeffs_size, double scale);

void shift_fft(COMPLEX8Vector *fft_out, COMPLEX8Vector *fft_in, COMPLEX8 *filter, int filter_size);

void test_bessel(void);

#endif