/*
 * Copyright (C) 2012 John Whelan
 * Based on XLALComputeAMTest.c Copyright (C) 2010 Reinhard Prix
 *                              Copyright (C) 2006 John Whelan, Reinhard Prix
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

/*********************************************************************************/
/** \author John Whelan
 * \file
 * \brief Test for XLALGetSSBtimes() and XLALGetMultiSSBtimes() by
 * comparison with the equivalent LAL functions LALGetSSBtimes() and
 * LALGetMultiSSBtimes()
 *
 * Note, we run a comparison only for the 2-IFO multiAM functions
 * XLALGetMultiSSBtimes() comparing it to LALGetMultiSSBtimes(), as
 * this excercises the 1-IFO functions as well.
 *
 * Sky-location and time are picked at random each time, which allows
 * a minimal Monte-Carlo validation by simply running this script
 * repeatedly.
 *
 *********************************************************************************/
#include <config.h>

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#ifdef HAVE_GETOPT_H
#include <getopt.h>
#endif

#include <math.h>
#include <sys/times.h>

#include <gsl/gsl_math.h>

#include <lal/ComputeFstat.h>
#include <lal/LALBarycenter.h>
#include <lal/LALInitBarycenter.h>
#include <lal/AVFactories.h>
extern char *optarg;

static const LALStatus empty_status;

extern int lalDebugLevel;

/* ----- internal prototypes ---------- */
int XLALCompareMultiSSBtimes ( MultiSSBtimes *multiSSB1, MultiSSBtimes *multiSSB2, REAL8 tol_dt, REAL8 tol_td );


/* ----- function definitions ---------- */
int main(int argc, char *argv[])
{
  int opt;             /* Command-line option. */

  UINT4 numIFOs = 2;
  const CHAR *sites[2] = {"H1", "V1"};

  UINT4 startTime = 888580733; /* should choose this randomly */
  LIGOTimeGPS refTime = {788600000, 0};
  UINT4 duration = 180000;	/* 50 hours */
  UINT4 Tsft = 1800;		/* assume 30min SFTs */

  REAL8 tol_dt = 2e-6;	/* same algorithm, should be basically identical results */
  REAL8 tol_td = 2e-6;	/* same algorithm, should be basically identical results */

  char earthEphem[] = "earth05-09.dat";
  char sunEphem[] = "sun05-09.dat";
  UINT4 numChecks = 1; /* Number of times to check */

  /* read user input */
  lalDebugLevel = 0;

  while ((opt = getopt( argc, argv, "n:qv:" )) != -1) {
    switch (opt) {
    case 'v': /* set lalDebugLevel */
      lalDebugLevel = atoi( optarg );
      break;
    case 'n': /* number of times to check */
      numChecks = atoi( optarg );
      break;
    }
  }

  /* init random-generator */
  struct tms buf;
  UINT4 seed = times(&buf);
  srand ( seed );
  XLALPrintInfo ("%s: seed used = %d\n", __func__, seed );

  /* ----- init ephemeris ----- */
  EphemerisData *edat;
  if ( (edat = XLALInitBarycenter ( earthEphem, sunEphem )) == NULL ) {
    XLALPrintError ("%s: XLALInitBarycenter() failed with xlalErrno = %d\n", __func__, xlalErrno );
    return XLAL_EFAILED;
  }

  /* ----- init detector info ---------- */
  UINT4 X;
  MultiLALDetector *multiDet;
  if ( (multiDet = XLALCreateMultiLALDetector ( numIFOs )) == NULL ) {
    XLALPrintError ("%s: XLALCreateMultiLALDetector(%d) failed with errno=%d\n", __func__, numIFOs, xlalErrno );
    return XLAL_EFAILED;
  }

  for (X=0; X < numIFOs; X ++ )
    {
      LALDetector *site;
      if ( (site = XLALGetSiteInfo ( sites[X] )) == NULL ) {
        XLALPrintError ("%s: Failed to get site-info for detector '%s'\n", __func__, sites[X] );
        return XLAL_EFAILED;
      }
      multiDet->data[X] = (*site); 	/* copy! */
      XLALFree ( site );
    }

  /* ----- init multi-IFO timestamps vector ---------- */
  UINT4 numSteps = (UINT4) ceil ( duration / Tsft );
  MultiLIGOTimeGPSVector * multiTS;
  if ( (multiTS = XLALCalloc ( 1, sizeof(*multiTS))) == NULL ) {
    XLAL_ERROR ( XLAL_ENOMEM );
  }
  multiTS->length = numIFOs;
  if ( (multiTS->data = XLALCalloc (numIFOs, sizeof(*multiTS->data))) == NULL ) {
    XLAL_ERROR ( XLAL_ENOMEM );
  }
  for ( X=0; X < numIFOs; X ++ )
    {
      if ( (multiTS->data[X] = XLALCreateTimestampVector (numSteps)) == NULL ) {
        XLALPrintError ("%s: XLALCreateTimestampVector(%d) failed.\n", __func__, numSteps );
        return XLAL_EFAILED;
      }
      multiTS->data[X]->deltaT = Tsft;

      UINT4 i;
      for ( i=0; i < numSteps; i ++ )
        {
          UINT4 ti = startTime + i * Tsft;
          multiTS->data[X]->data[i].gpsSeconds = ti;
          multiTS->data[X]->data[i].gpsNanoSeconds = 0;
        } /* for i < numSteps */

    } /* for X < numIFOs */

  /* ---------- compute multi-detector states -------------------- */
  MultiDetectorStateSeries *multiDetStates;
  if ( (multiDetStates = XLALGetMultiDetectorStates ( multiTS, multiDet, edat, 0.5 * Tsft )) == NULL ) {
    XLALPrintError ( "%s: XLALGetMultiDetectorStates() failed.\n", __func__ );
    return XLAL_EFAILED;
  }
  XLALDestroyMultiLALDetector ( multiDet );
  XLALDestroyMultiTimestamps ( multiTS );
  XLALFree(edat->ephemE);
  XLALFree(edat->ephemS);
  XLALFree ( edat );

  /* ========== MAIN LOOP: N-trials of comparisons XLAL <--> LAL SSBtimes functions ========== */
  while ( numChecks-- )
    {
      LALStatus status = empty_status;

      /* ----- pick skyposition at random ----- */
      SkyPosition skypos;
      skypos.longitude = LAL_TWOPI * (1.0 * rand() / ( RAND_MAX + 1.0 ) );  /* uniform in [0, 2pi) */
      skypos.latitude = LAL_PI_2 - acos ( 1 - 2.0 * rand()/RAND_MAX );	/* sin(delta) uniform in [-1,1] */
      skypos.system = COORDINATESYSTEM_EQUATORIAL;

      /* should change this to loop over precision options */
      SSBprecision prec = SSBPREC_NEWTONIAN;

      /* ----- compute multi SSB times using LAL function ----- */
      MultiSSBtimes *multiSSB_LAL  = NULL;
      LALGetMultiSSBtimes ( &status,
			    &multiSSB_LAL,
			    multiDetStates,
			    skypos,
			    refTime,
			    prec );
      if ( status.statusCode ) {
        XLALPrintError ("%s: LALGetMultiSSBtimes() failed with statusCode = %d : %s\n", __func__, status.statusCode, status.statusDescription );
        return XLAL_EFAILED;
      }

      /* ----- compute multi SSB times using XLAL function ----- */
      MultiSSBtimes *multiSSB_XLAL;
      if ( ( multiSSB_XLAL = XLALGetMultiSSBtimes ( multiDetStates, skypos, refTime, prec )) == NULL ) {
        XLALPrintError ("%s: XLALGetMultiSSBtimes() failed with xlalErrno = %d\n", __func__, xlalErrno );
        return XLAL_EFAILED;
      }

      /* now run comparison */
      if ( XLALCompareMultiSSBtimes ( multiSSB_XLAL, multiSSB_LAL, tol_dt, tol_td ) != XLAL_SUCCESS ) {
        XLALPrintError ("%s: comparison between multiAM_XLAL and multiAM_LAL failed.\n", __func__ );
        return XLAL_EFAILED;
      }

      /* free memory created inside this loop */
      XLALDestroyMultiSSBtimes ( multiSSB_LAL );
      XLALDestroyMultiSSBtimes ( multiSSB_XLAL );

    } /* for numChecks */

  /* we're done: free memory */
  XLALDestroyMultiDetectorStateSeries ( multiDetStates );

  LALCheckMemoryLeaks();

  printf ("OK. All tests passed successfully\n\n");

  return 0;	/* OK */

} /* main() */

/** Comparison function for two multiSSB vectors, return success or failure for given tolerances.
 *
 * we compare avg and max of |DeltaT1-DeltaT2| and |Tdot1-Tdot2| to tol_dt and tol_td, respectively; the caller should set reasonable tolerances for these
 *
 */
int
XLALCompareMultiSSBtimes ( MultiSSBtimes *multiSSB1, MultiSSBtimes *multiSSB2, REAL8 tol_dt, REAL8 tol_td )
{
  /* check input */
  if ( !multiSSB1 || !multiSSB2 || tol_dt <= 0 || tol_td <= 0 ) {
    XLALPrintError ("%s: invalid NULL input or non-positive tolerance\n", __func__ );
    XLAL_ERROR ( XLAL_EINVAL );
  }

  UINT4 numDet = multiSSB1->length;
  if ( numDet != multiSSB2->length ) {
    XLALPrintError ("%s: number of detectors differ multiSSB1 = %d, multiSSB2 = %d\n", __func__, multiSSB1->length, multiSSB2->length );
    XLAL_ERROR ( XLAL_EFAILED );
  }

  UINT4 X;
  REAL8 maxerr_dt = 0, maxerr_td = 0, avgerr_dt = 0, avgerr_td = 0;
  UINT4 numTerms = 0;
  for ( X=0; X < numDet; X ++ )
    {
      UINT4 numSteps = multiSSB1->data[X]->DeltaT->length;
      if ( numSteps != multiSSB2->data[X]->DeltaT->length ) {
        XLALPrintError ("%s: number of timesteps differ multiSSB1[%d]->DeltaT = %d, multiSSB2[%d]->DeltaT = %d\n",__func__, X, multiSSB1->data[X]->DeltaT->length, X, multiSSB2->data[X]->DeltaT->length );
        XLAL_ERROR ( XLAL_EFAILED );
      }
      if ( numSteps != multiSSB1->data[X]->Tdot->length || numSteps != multiSSB2->data[X]->Tdot->length) {
        XLALPrintError ("%s: number of timesteps differ multiSSB1[%d]->Tdot = %d, multiSSB2[%d]->Tdot = %d\n",__func__, X, multiSSB1->data[X]->Tdot->length, X, multiSSB2->data[X]->Tdot->length );
        XLAL_ERROR ( XLAL_EFAILED );
      }

      UINT4 i;
      for ( i=0; i < numSteps; i ++ )
        {

          REAL8 err_dt = fabs ( multiSSB1->data[X]->DeltaT->data[i] - multiSSB2->data[X]->DeltaT->data[i] );
          REAL8 err_td = fabs ( multiSSB1->data[X]->Tdot->data[i] - multiSSB2->data[X]->Tdot->data[i] );

          if ( err_dt > maxerr_dt )
            maxerr_dt = err_dt;
          if ( err_td > maxerr_td )
            maxerr_td = err_td;

          avgerr_dt += err_dt;
          avgerr_td += err_td;
          numTerms += 1;

	  XLALPrintInfo("DeltaT: LAL: %.16f XLAL: %.16f diff: %.8g\n",
			multiSSB1->data[X]->DeltaT->data[i],
			multiSSB2->data[X]->DeltaT->data[i],
			multiSSB1->data[X]->DeltaT->data[i]
			- multiSSB2->data[X]->DeltaT->data[i]);

	  XLALPrintInfo("Tdot: LAL: %.16f XLAL: %.16f diff: %.8g\n",
			multiSSB1->data[X]->Tdot->data[i],
			multiSSB2->data[X]->Tdot->data[i],
			multiSSB1->data[X]->Tdot->data[i]
			- multiSSB2->data[X]->Tdot->data[i]);

        } /* for i < numSteps */

    } /* for X < numDet */

  avgerr_dt /= numTerms;
  avgerr_td /= numTerms;

  UINT4 failed = 0;

  /* ----- compare individual DeltaT errors -------------------- */
  if ( maxerr_dt > tol_dt ) {
    XLALPrintError ("%s: maximal difference in DeltaT is %g, which exceeds the tolerance %g\n", __func__, maxerr_dt, tol_dt );
    failed ++;
  }
  else
    XLALPrintInfo ("%s: maxerr_dt = %g (< %g)\n", __func__, maxerr_dt, tol_dt);

  if ( avgerr_dt > tol_dt ) {
    XLALPrintError ("%s: average difference in DeltaT is %g, which exceeds the tolerance %g\n", __func__, avgerr_dt, tol_dt );
    failed ++;
  }
  else
    XLALPrintInfo ("%s: avgerr_dt= %g (< %g)\n", __func__, avgerr_dt, tol_dt);

  /* ----- compare individual Tdot errors -------------------- */
  if ( maxerr_td > tol_td ) {
    XLALPrintError ("%s: maximal difference in Tdot is %g, which exceeds the tolerance %g\n", __func__, maxerr_td, tol_td );
    failed ++;
  }
  else
    XLALPrintInfo ("%s: maxerr_td = %g (< %g)\n", __func__, maxerr_td, tol_td);

  if ( avgerr_td > tol_td ) {
    XLALPrintError ("%s: average difference in Tdot is %g, which exceeds the tolerance %g\n", __func__, avgerr_td, tol_td );
    failed ++;
  }
  else
    XLALPrintInfo ("%s: avgerr_td= %g (< %g)\n", __func__, avgerr_td, tol_td);

  return failed;

} /* XLALCompareMultiSSBtimes() */
