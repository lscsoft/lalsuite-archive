/*
 * Copyright (C) 2011 J. Clark
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

#include <complex.h>

#include <lal/LALSimRingdown.h>
#include <lal/LALConstants.h>
#include <lal/LALStdlib.h>
#include <lal/FrequencySeries.h>
#include <lal/Units.h>
#include <lal/SphericalHarmonics.h>
#include <lal/Date.h>

#define EPS LAL_REAL4_EPS

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif


/**
 * Computes the frequency domain representation for a generic ringdown.  NOTE:
 * amplitude is not (currently) handled in this function!
 *
 */
int XLALSimRingdownFD(
	COMPLEX16FrequencySeries **htilde,
    double f_min,
    double f_max,
	double deltaF,			/**< sampling interval (s) */
	double omega0,			/* f-mode oscillation frequency (rad) */
	double quality			/* quality factor = pi*decay*frequency*/
)
{
    static LIGOTimeGPS ligotimegps_zero = {0, 0};
    size_t i;

    /* allocate htilde */
    size_t n =ceil(f_max / deltaF) + 1;
    *htilde = XLALCreateCOMPLEX16FrequencySeries("htilde: FD waveform",
            &ligotimegps_zero, 0.0, deltaF, &lalStrainUnit, n);

    memset((*htilde)->data->data, 0, n * sizeof(COMPLEX16));

    /* Get tau from quality factor and frequency */
	/* From T040055-00-Z: Q = sqrt(2)*PI*f0*tau */
    double tau = sqrt(2.0)*quality/omega0;

    if (!(*htilde)) XLAL_ERROR(XLAL_EFUNC);

    /* Generate the waveform */
    size_t ind_max = (size_t) (f_max / deltaF);
    for (i = (size_t) (f_min / deltaF); i < ind_max; i++) {

        /* Fourier frequency corresponding to this bin */
        double f = i * deltaF;
        double w = LAL_TWOPI*f;

        ((*htilde)->data->data)[i] = 1.0/(2.0*(w+omega0)+2.0*I/tau);
        ((*htilde)->data->data)[i] -= 1.0/(2.0*(w-omega0)+2.0*I/tau);
    }

	return XLAL_SUCCESS;
}


