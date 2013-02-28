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
#include <lal/LIGOMetadataTables.h>
#include <lal/FrequencySeries.h>
#include <lal/Sequence.h>
#include <lal/TimeFreqFFT.h>
#include <lal/TimeSeries.h>
#include <lal/RealFFT.h>
#include <lal/Units.h>
#include <lal/Date.h>
#include "check_series_macros.h"

#define UNUSED(expr) do { (void)(expr); } while (0)

unsigned int i, j;

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
  REAL8 timeInj, timeData, lag, duration_sampled;
  SimBurst *Dummy;
  REAL8TimeSeries *h_plus, *h_cross;
  gsl_rng * rng;
  /*    fprintf(stdout, "Data length is %u\n", time_series->data->length);
      fprintf(stdout, "First entry is %e\n", time_series->data->data[0]);
    XLALSimGlitchPrintSimBurstTable(glitchInjections);
  */  timeData = (REAL8)time_series->epoch.gpsSeconds + 1.0e-9 * ((REAL8)time_series->epoch.gpsNanoSeconds);
    //  fprintf(stdout, "Data GPS time is %f\n", timeData);
    for(Dummy = glitchInjections; Dummy; Dummy = Dummy->next)
      {	
	timeInj = Dummy->time_geocent_gmst;
	lag = (timeInj - timeData) / time_series->deltaT;
	duration_sampled = Dummy->duration / time_series->deltaT;
	//	fprintf(stdout, "timeInj = %d\t", (int)timeInj);
	if(!strcmp(Dummy->waveform, "SineGaussian"))
	    {
	      //    fprintf(stdout, "This is a Sine Gaussian injection\n");
	      XLALSimBurstSineGaussian(
					&h_plus,
					&h_cross,
					Dummy->q, // Q parameter
					Dummy->frequency, // centre_frequency
					Dummy->hrss, // hrss
					1, // eccentricity
					0, // polarization, not relevant when eccent = 0
					time_series->deltaT // delta_t
					);
	      // fprintf(stdout, "Frequency = %f\t hrss = %e\n\n", Dummy->frequency, Dummy->hrss);
	      for(i=0; i<h_plus->data->length; i++)
		{
		  j = lag - h_plus->data->length / 2 + i;
		  if(j >= time_series->data->length) continue;
		  time_series->data->data[j] += h_plus->data->data[i];
		  // fprintf(stdout, "%d\t%e\n", i, time_series->data->data[j]);
		} 
	    }
	else if(!strcmp(Dummy->waveform, "BTLWNB"))
	  {
	    // fprintf(stdout, "This is a White Noise injection\n");
	    rng = gsl_rng_alloc (gsl_rng_taus);
	    XLALGenerateBandAndTimeLimitedWhiteNoiseBurst(
							  &h_plus,
							  &h_cross,
							  Dummy->duration,
							  Dummy->frequency,
							  Dummy->bandwidth,
							  Dummy->egw_over_rsquared, // int_hdtot_squared
							  time_series->deltaT,
							  rng
							  );
	    // fprintf(stdout, "frequency = %f\t amplitude = %e\n\n", Dummy->frequency, Dummy->amplitude);
	    for(i=0; i<h_plus->data->length; i++)
	      {
		j = lag - h_plus->data->length / 2 + i;
		if(j >= time_series->data->length) continue;
		time_series->data->data[j] += h_plus->data->data[i];
		// fprintf(stdout, "%d\t%e\n", i, time_series->data->data[j]);                  
	      }
	  }
      }
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

void XLALSimGlitchNewSimBurst(
    SimBurst **ret_sim,
    SimBurst *next_sim,
    INT4 gpsSeconds,
    INT4 gpsNanoSeconds,
    const char *waveform,
    REAL8 duration,
    REAL8 frequency,
    REAL8 bandwidth,
    REAL8 q,
    REAL8 amplitude,
    REAL8 hrss,
    long process_id)
{
    SimBurst *new_sim;
    *ret_sim = XLALMalloc(sizeof(*new_sim));
    new_sim = *ret_sim;

    new_sim->next = next_sim;

    new_sim->time_geocent_gps.gpsSeconds = gpsSeconds;
    new_sim->time_geocent_gps.gpsNanoSeconds = gpsNanoSeconds;
    new_sim->time_geocent_gmst = (REAL8) gpsSeconds + 1.0e-9 * ((REAL8) gpsNanoSeconds);

    strcpy(new_sim->waveform, waveform);

    new_sim->process_id = process_id;
    new_sim->duration = duration;
    new_sim->frequency = frequency;
    new_sim->bandwidth = bandwidth;
    new_sim->q = q;
    new_sim->amplitude = amplitude;
    new_sim->hrss = hrss;

    ret_sim = &new_sim;
}

void XLALSimGlitchPrintSimBurstTable(
    SimBurst *head)
{
    SimBurst *thisSim;

    for(thisSim = head; thisSim; thisSim = thisSim->next)
    {
        fprintf(stdout, "\n=== SimBurst ===\n");
        fprintf(stdout, "Waveform = %s\n", thisSim->waveform);
        fprintf(stdout, "Time = %u.%.9u\n", thisSim->time_geocent_gps.gpsSeconds, thisSim->time_geocent_gps.gpsNanoSeconds);
        fprintf(stdout, "Frequency = %f\n", thisSim->frequency);
        fprintf(stdout, "Duration = %f\n", thisSim->duration);
        fprintf(stdout, "Bandwidth = %f\n", thisSim->bandwidth);
        fprintf(stdout, "Q = %f\n", thisSim->q);
        fprintf(stdout, "Amplitude = %g\n", thisSim->amplitude);
        fprintf(stdout, "hrss = %g\n", thisSim->hrss);
    }
}
/*
	*hplus = XLALCreateREAL8TimeSeries("Impulse +", &epoch, 0.0, delta_t, &lalStrainUnit, length);

	(*hplus)->data->data[(length - 1) / 2] = hpeak;
*/
