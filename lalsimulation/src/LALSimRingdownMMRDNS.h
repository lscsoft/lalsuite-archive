#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#ifndef _LALSIM_RINGDOWN_MMRDNS_H
#define _LALSIM_RINGDOWN_MMRDNS_H


/* ************************************************************  */
/*
 * Copyright (C) 2016 Lionel London
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

 #if defined(__cplusplus)
 extern "C" {
 #elif 0
 } /* so that editors will match preceding brace */
 #endif

/* Include the desired Libs */
#include <stdbool.h>
#include <math.h>
#include <complex.h>
/* LAL specific libs  */

#include <lal/Date.h>
#include <lal/FrequencySeries.h>
#include <lal/Units.h>

#include <lal/LALDatatypes.h>
#include <lal/LALStdlib.h>
#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>
#include <lal/LALConfig.h>
#include <lal/LALConstants.h>

#include <lal/LALConstants.h>
#include <lal/LALStdio.h>
#include <lal/LALSimSphHarmSeries.h>
#include <lal/LALStdlib.h>
#include <lal/LALSimInspiral.h>

#include <lal/TimeSeries.h>
#include <lal/FrequencySeries.h>

/* ---------------------------------------- */
/* General model methods and parameters     */
/* ---------------------------------------- */

/* Gnerate Time domain ringdown waveform  */
UNUSED int XLALSimRingdownMMRDNSTD(
        UNUSED REAL8TimeSeries **hplus,        /**< OUTPUT h_+ vector */
        UNUSED REAL8TimeSeries **hcross,       /**< OUTPUT h_x vector */
        UNUSED REAL8 phiRef,                   /**< orbital phase at reference pt. */
        UNUSED REAL8 inclination,              /**< inclination angle */
        UNUSED REAL8 deltaT,                   /**< sampling interval (s) */
        UNUSED REAL8 m1,                       /**< mass of companion 1 (kg) */
        UNUSED REAL8 m2,                       /**< mass of companion 2 (kg) */
        UNUSED REAL8 r,                        /**< distance of source (m) */
        UNUSED const LALSimInspiralTestGRParam *extraParams /**< linked list containing the extra testing GR parameters */
        );

/* XLALSimRingdownSingleModeTD: Time domain waveformgenerator for single QNM without angular dependence (i.e. this function generates multipole moments only ). In  */
UNUSED int XLALSimRingdownGenerateSingleModeTD(
  UNUSED COMPLEX16TimeSeries **hk,        /**< OUTPUT vector for QNM time series */
  UNUSED const LIGOTimeGPS T0,            /**< Start time  */
  UNUSED REAL8 deltaT,                    /**< sampling interval (s) */
  UNUSED REAL8 Nsamples,                  /**< Number of time domain samples */
  UNUSED complex double Ak,               /**< COMPLEX QNM Strain Amplitude */
  UNUSED complex double CWk               /**< COMPLEX QNM Frequency */
);

/* XLALSimRingdownMMRDNSFD: Frequency domain waveformgenerator for all QNM with angular dependence */
int XLALSimRingdownMMRDNSFD(
        COMPLEX16FrequencySeries **hptilde,          /**< OUTPUT FD h_+ polarization */
        COMPLEX16FrequencySeries **hctilde,          /**< OUTPUT FD h_x polarization */
        const REAL8 deltaF,                          /**< Frequency resolution (Hz) */
        const REAL8 fStart,                          /**< Start GW frequency (Hz) */
        const REAL8 fEnd,                            /**< Highest GW frequency (Hz) */
        REAL8 Mf,                                    /**< Final BH Mass (kg) */
        REAL8 jf,                                    /**< Final BH dimensionaless spin */
        REAL8 eta,                                   /**< Symmetric mass ratio of two companions */
        REAL8 iota,                                  /**< inclination angle (in rad) */
        REAL8 phi_offset,                            /**< intrinsic phase offset (in rad) */
        REAL8 r,                                     /**< distance of source (m) */
        LALSimInspiralTestGRParam *nonGRparams       /**< testing GR parameters */
);

/* XLALSimRingdownGenerateSingleModeFD: Frequency domain waveformgenerator for single QNM with angular dependence */
int XLALSimRingdownGenerateSingleModeFD(
      COMPLEX16FrequencySeries **hptilde_lmn,      /**< OUTPUT FD h_+ polarization */
      COMPLEX16FrequencySeries **hctilde_lmn,      /**< OUTPUT FD h_x polarization */
      const REAL8 deltaF,                          /**< Frequency resolution (Hz) */
      const REAL8 fStart,                          /**< Start GW frequency (Hz) */
      const REAL8 fEnd,                            /**< Highest GW frequency (Hz) */
      REAL8 Mf,                                    /**< Final BH Mass (solar mass) */
      REAL8 jf,                                    /**< Final BH dimensionaless spin */
      REAL8 eta,                                   /**< Symmetric mass ratio of two companions */
      REAL8 iota,                                  /**< inclination angle */
      REAL8 phi_offset,                            /**< intrinsic phase offset */
      UINT4 l,                                     /**< Polar eigenvalue */
      UINT4 m,                                     /**< Azimuthal eigenvalue */
      UINT4 n,                                     /**< Overtone Number */
      REAL8 r,                                     /**< distance of source (m) */
      REAL8 dfreq,                                 /**< relative shift in the real frequency parameter */
      REAL8 dtau                                   /**< relative shift in the damping time parameter */
);

/* XLALSimRingdownMMRDNS_time: Time domain waveformgenerator for all QNM with angular dependence */
int XLALSimRingdownMMRDNS_time(
        REAL8TimeSeries **hplus,                     /**< OUTPUT TD waveform */
        REAL8TimeSeries **hcross,                     /**< OUTPUT TD waveform */
        const LIGOTimeGPS *t0,                       /**< start time of ringdown => NEEDS TO BE CHECKED! */
        REAL8 deltaT,                                /**< sampling interval (s) */
        REAL8 Mf,                                    /**< Final BH Mass (kg) */
        REAL8 jf,                                    /**< Final BH dimensionaless spin */
        REAL8 eta,                                   /**< Symmetric mass ratio of two companions */
        REAL8 iota,                                  /**< inclination angle (in rad) */
        REAL8 phi_offset,                            /**< intrinsic phase offset */
        REAL8 r,                                     /**< distance of source (m) */
        LALSimInspiralTestGRParam *nonGRparams       /**< testing GR parameters */
);

/* XLALSimRingdownGenerateSingleModeMMRDNS_time: Time domain waveformgenerator for single QNM with angular dependence */
int XLALSimRingdownGenerateSingleModeMMRDNS_time(
        COMPLEX16TimeSeries **htilde_lmn,            /**< OUTPUT TD waveform mode lmn */
        const LIGOTimeGPS *t0,                       /**< start time of ringdown */
        REAL8 deltaT,                                /**< sampling interval (s) */
        REAL8 Mf,                                    /**< Final BH Mass (kg) */
        REAL8 jf,                                    /**< Final BH dimensionaless spin */
        REAL8 eta,                                   /**< Symmetric mass ratio of two companions */
        REAL8 iota,                                  /**< inclination angle (in rad) */
        REAL8 phi_offset,                            /**< intrinsic phase offset (in rad) */
        UINT4 l,                                     /**< Polar eigenvalue */
        UINT4 m,                                     /**< Azimuthal eigenvalue */
        UINT4 n,                                     /**< Overtone Number */
        REAL8 r,                                     /**< distance of source (m) */
        REAL8 dfreq,                                 /**< relative shift in the real frequency parameter */
        REAL8 dtau,                                  /**< relative shift in the damping time parameter */
        UINT4 Nsamples,                               /**< waveform length */
        REAL8 Tstart                                 /**< starting time of waveform (10M at zero) */

);

/*
* Spheroical Harmonic Functions (Leaver's Formulation circa 1986/85)
*/
COMPLEX16 XLALSpinWeightedSpheroidalHarmonic( double jf, int l, int m, int n, double theta, double phi);

/*
* Domain mapping for dimnesionless BH spin
*/
REAL8 XLALKAPPA( double jf, int l, int m );

/*
* -------------------------------------------------------------------------------- *
* Low level models: QNM Frequencies, Separation Constants and Spheroidal Harmonics
* -------------------------------------------------------------------------------- *
*/

/*
* QNM Separation Constants: Note that name encodes date of writing
*/
COMPLEX16 XLALseparationConstant( double kappa, int l, int input_m, int n );

/*
* Dimensionless QNM Frequencies: Note that name encodes date of writing
*/
COMPLEX16 XLALcomplexOmega( double kappa, int l, int input_m, int n );




/* Convert NR Code Time to Physical Units */
//UNUSED static double XLALNRtoPhysTime( UNUSED double NRtime  );

/* ************************************************************  */
#endif	/* of #ifndef _LALSIM_RINGDOWN_MMRDNS_H */
