#ifndef __BESSEL_H__
#define __BESSEL_H__

#include <lal/LALDatatypes.h>

void make_bessel_filter(COMPLEX16 *filter, int filter_size, COMPLEX16 *coeffs, int coeffs_size, double scale);

#define make_bessel_filter7(scale, coeff, filter)  make_bessel_filter(filter, 7, coeff, 3, scale)

void shift_fft(COMPLEX16Vector *fft_out, COMPLEX16Vector *fft_in, COMPLEX16 *filter, int filter_size);

void test_bessel(void);

#endif