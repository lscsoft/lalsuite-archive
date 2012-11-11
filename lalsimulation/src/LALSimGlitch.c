/*
 * Copyright (C) 2012 V. Lockett, A. Lundgren, J. Smith
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with with program; see the file COPYING. If not, write to the Free
 * Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 */


/*
 * ============================================================================
 *
 *                                  Preamble
 *
 * ============================================================================
 */


#define LAL_USE_OLD_COMPLEX_STRUCTS
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <lal/LALSimBurst.h>
#include <lal/LALSimGlitch.h>
#include <lal/LALComplex.h>
#include <lal/LALConstants.h>
#include <lal/LIGOLwXMLBurstRead.h>
#include <lal/FrequencySeries.h>
#include <lal/Sequence.h>
#include <lal/TimeFreqFFT.h>
#include <lal/TimeSeries.h>
#include <lal/RealFFT.h>
#include <lal/Units.h>
#include <lal/Date.h>
#include "check_series_macros.h"

#define UNUSED(expr) do { (void)(expr); } while (0)

/*
 * ============================================================================
 *
 *          Inject simulated data into a time stream
 *
 * ============================================================================
 */

int XLALSimGlitchInject(
    REAL8TimeSeries *time_series,
    SimBurst *glitchInjections
)
{
    UNUSED(glitchInjections);
    fprintf(stdout, "Data length is %u\n", time_series->data->length);
    fprintf(stdout, "First entry is %f\n", time_series->data->data[0]);
    return 0;
}

void XLALSimGlitchNewTimeSeries(
    REAL8TimeSeries **time_series,
    const CHAR * name,
    INT4 start_time,
    REAL8 deltaT,
    size_t length)
{
    LIGOTimeGPS tStart;
    tStart.gpsSeconds = start_time;
    tStart.gpsNanoSeconds = 0;
    *time_series = XLALCreateREAL8TimeSeries( name, &tStart, 0.0, deltaT, &lalStrainUnit, length );
    memset((*time_series)->data->data, 0, (*time_series)->data->length*sizeof(*(*time_series)->data->data));

    return;
}
/*
	*hplus = XLALCreateREAL8TimeSeries("Impulse +", &epoch, 0.0, delta_t, &lalStrainUnit, length);

	(*hplus)->data->data[(length - 1) / 2] = hpeak;
*/
