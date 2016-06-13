/*
*  Copyright (C) 2007 Vladimir Dergachev
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

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
