#ifndef __OP_METHOD_H__
#define __OP_METHOD_H__

/* these data structures are only meant for quadratic interpolation
of R_n(t), approx=3 */

typedef struct {
	double a[5];
	double det;
	double inverse[3][3];
	} PHI_DATA3;

typedef struct {
	double phi_r_re[3];
	double phi_r_im[3];
	} PHI_RESPONSE3;
	
void compute_test_fft(REAL4Vector *data_in, COMPLEX8Vector **phi, int approx, long freq_start, long freq_stop);
void compute_calibrated_periodogram3(COMPLEX8Vector **phi, long freq_start, long freq_stop, COMPLEX8Vector *pgram, double window_sum);

#endif
