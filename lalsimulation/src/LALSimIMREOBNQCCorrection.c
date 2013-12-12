/*
*  Copyright (C) 2010 Craig Robinson, Yi Pan
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


/**
 * \author Craig Robinson, Yi Pan
 *
 * \brief More recent versions of the EOB models, such as EOBNRv2 and SEOBNRv1, utilise
 * a non-quasicircular correction (NQC) to bring the peak of the EOB frequency
 * into agreement with that of NR simulations. This file contains the functions
 * used to calculate these NQC corrections, described in DCC document T1100433.
 * The fits to NR peak amplitude, frequency, and their derivatives, are taken
 * from Pan et al. PRD 84 124052 (2011), for EOBNRv2, and
 * from Taracchini et al. PRD 86, 024011 (2012), for SEOBNRv1.
 */

#include <complex.h>
#include <math.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_linalg.h>

#include "LALSimIMREOBNRv2.h"

#ifndef _LALSIMIMRNQCCORRECTION_C 
#define _LALSIMIMRNQCCORRECTION_C 

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/* ------------------------------------------------
 *          Non-spin (EOBNRv2)
 * ------------------------------------------------*/

/**
 * Compute the time offset which should be used in computing the
 * non-quasicircular correction and performing the ringdown attachment.
 * These numbers were tuned to numerical relativity simulations, and
 * are taken from Pan et al, PRD84, 124052(2011), lines 1-5 of Table II.
 */
 static REAL8 XLALSimIMREOBGetNRPeakDeltaT( 
                         INT4 l,    /**<< Mode l */ 
                         INT4 m,    /**<< Mode m */
                         REAL8 eta  /**<< Symmetric mass ratio */
                         )
{
  switch ( l )
  {
    case 2:
      if ( m == 2 )
      {
        return 0.0;
      }
      else if ( m == 1 )
      {
        return 10.67 - 2.5 + 9.0*eta - 41.41 * eta + 76.1 * eta*eta;
      }
      else
      {
        XLAL_ERROR_REAL8( XLAL_EINVAL );
      }
      break;
    case 3:
      if ( m == 3 )
      {
        return 3.383 + 3.847 * eta + 8.979 * eta * eta;
      }
      else
      {
        XLAL_ERROR_REAL8( XLAL_EINVAL );
      }
      break;
    case 4:
      if ( m == 4 )
      {
        return 5.57 - 49.86 * eta + 154.3 * eta * eta;
      }
      else
      {
        XLAL_ERROR_REAL8( XLAL_EINVAL );
      }
      break;
    case 5:
      if ( m == 5 )
      {
        return 6.693 - 34.47 * eta + 102.7 * eta*eta;
      }
      else
      {
        XLAL_ERROR_REAL8( XLAL_EINVAL );
      }
      break;
    default:
      XLAL_ERROR_REAL8( XLAL_EINVAL );
      break;
  }
}

/**
 * Function which returns a value of the expected peak amplitude
 * taken from a fit to numerical relativity simulations. The functions
 * are taken from Pan et al, PRD84, 124052(2011), lines 1-5 of Table II.
 */
static inline
REAL8 GetNRPeakAmplitude( 
                        INT4 l,   /**<< Mode l */ 
                        INT4 m,   /**<< Mode m */
                        REAL8 eta /**<< Symmetric mass ratio */
                        )
{
  switch ( l )
  {
    case 2:
      if ( m == 2 )
      {
        return eta * ( 1.422 + 0.3013 * eta + 1.246 * eta * eta );
      }
      else if ( m == 1 )
      {
        return eta * sqrt( 1.0 - 4. * eta ) * (0.4832 - 0.01032 * eta);
      }
      else
      {
        XLAL_ERROR_REAL8( XLAL_EINVAL );
      }
      break;
    case 3:
      if ( m == 3 )
      {
        return eta * sqrt(1.-4.*eta) * ( 0.5761 - 0.09638 * eta + 2.715*eta*eta );
      }
      else
      {
        XLAL_ERROR_REAL8( XLAL_EINVAL );
      }
      break;
    case 4:
      if ( m == 4 )
      {
        return eta * (0.354 - 1.779 * eta + 2.834 * eta*eta );
      }
      else
      {
        XLAL_ERROR_REAL8( XLAL_EINVAL );
      }
      break;
    case 5:
      if ( m == 5 )
      {
        return eta * sqrt(1.-4.*eta) * ( 0.1353 - 0.1485 * eta );
      }
      else
      {
        XLAL_ERROR_REAL8( XLAL_EINVAL );
      }
      break;
    default:
      XLAL_ERROR_REAL8( XLAL_EINVAL );
      break;
  }
}

/**
 * Function which returns second derivative of the amplitude at the peak
 * taken from a fit to numerical relativity simulations. The functions
 * are taken from Pan et al, PRD84, 124052(2011), lines 1-5 of Table II.
 */
static inline
REAL8 GetNRPeakADDot( 
                    INT4 l,   /**<< Mode l */ 
                    INT4 m,   /**<< Mode m */
                    REAL8 eta /**<< Symmetric mass ratio */
                    )
{
  switch ( l )
  {
    case 2:
      if ( m == 2 )
      {
        return -0.01 * eta * ( 0.1679 + 1.44 * eta - 2.001 * eta * eta );
      }
      else if ( m == 1 )
      {
        return -0.01 * eta * sqrt(1.-4.*eta) * (0.1867 + 0.6094 * eta );
      }
      else
      {
        XLAL_ERROR_REAL8( XLAL_EINVAL );
      }
      break;
    case 3:
      if ( m == 3 )
      {
        return -0.01 * eta * sqrt(1.-4.*eta) * (0.2518 - 0.8145*eta + 5.731*eta*eta);
      }
      else
      {
        XLAL_ERROR_REAL8( XLAL_EINVAL );
      }
      break;
    case 4:
      if ( m == 4 )
      {
        return -0.01 * eta * (0.1813 - 0.9935 * eta + 1.858 * eta * eta );
      }
      else
      {
        XLAL_ERROR_REAL8( XLAL_EINVAL );
      }
      break;
    case 5:
      if ( m == 5 )
      {
        return -0.01 * eta * sqrt(1.-4.*eta) * ( 0.09051 - 0.1604 * eta );
      }
      else
      {
        XLAL_ERROR_REAL8( XLAL_EINVAL );
      }
      break;
    default:
      XLAL_ERROR_REAL8( XLAL_EINVAL );
      break;
  }
}


/**
 * Function which returns a value of the expected peak frequency
 * taken from a fit to numerical relativity simulations. The functions
 * are taken from Pan et al, PRD84, 124052(2011), lines 1-5 of Table II.
 */
static inline 
REAL8 GetNRPeakOmega( 
                    INT4 l,   /**<< Mode l */
                    INT4 m,   /**<< Mode m */
                    REAL8 eta /**<< Symmetric mass ratio */
                    )
{
  switch ( l )
  {
    case 2:
      if ( m == 2 )
      {
        return 0.2733 + 0.2316 * eta + 0.4463 * eta * eta;
      }
      else if ( m == 1 )
      {
        return 0.2907 - 0.08338 * eta + 0.587 * eta*eta;
      }
      else
      {
        XLAL_ERROR_REAL8( XLAL_EINVAL );
      }
      break;
    case 3:
      if ( m == 3 )
      {
        return 0.4539 + 0.5376 * eta + 1.042 * eta*eta;
      }
      else
      {
        XLAL_ERROR_REAL8( XLAL_EINVAL );
      }
      break;
    case 4:
      if ( m == 4 )
      {
        return 0.6435 - 0.05103 * eta + 2.216 * eta*eta;
      }
      else
      {
        XLAL_ERROR_REAL8( XLAL_EINVAL );
      }
      break;
    case 5:
      if ( m == 5 )
      {
        return 0.8217 + 0.2346 * eta + 2.599 * eta*eta;
      }
      else
      {
        XLAL_ERROR_REAL8( XLAL_EINVAL );
      }
      break;
    default:
      XLAL_ERROR_REAL8( XLAL_EINVAL );
      break;
  }
}

/**
 * Function which returns the derivative of the expected peak frequency
 * taken from a fit to numerical relativity simulations. The functions
 * are taken from Pan et al, PRD84, 124052(2011), lines 1-5 of Table II.
 */
static inline 
REAL8 GetNRPeakOmegaDot( 
                       INT4 l,   /**<< Mode l */ 
                       INT4 m,   /**<< Mode m */
                       REAL8 eta /**<< Symmetric mass ratio */
                       )
{
  switch ( l )
  {
    case 2:
      if ( m == 2 )
      {
        return 0.005862 + 0.01506 * eta + 0.02625 * eta * eta;
      }
      else if ( m == 1 )
      {
        return 0.00149 + 0.09197 * eta - 0.1909 * eta*eta;
      }
      else
      {
        XLAL_ERROR_REAL8( XLAL_EINVAL );
      }
      break;
    case 3:
      if ( m == 3 )
      {
        return 0.01074 + 0.0293 * eta + 0.02066 * eta*eta;
      }
      else
      {
        XLAL_ERROR_REAL8( XLAL_EINVAL );
      }
      break;
    case 4:
      if ( m == 4 )
      {
        return 0.01486 + 0.08529 * eta - 0.2174 * eta * eta;
      }
      else
      {
        XLAL_ERROR_REAL8( XLAL_EINVAL );
      }
      break;
    case 5:
      if ( m == 5 )
      {
        return 0.01775 + 0.09801 * eta - 0.1686 * eta*eta;
      }
      else
      {
        XLAL_ERROR_REAL8( XLAL_EINVAL );
      }
      break;
    default:
      XLAL_ERROR_REAL8( XLAL_EINVAL );
      break;
  }
}


/**
 * For the 2,2 mode, there are fits available for the NQC coefficients,
 * given in Eqs.(40a)-(40c) of Pan et al, PRD84, 124052(2011).
 * This function provides the values of these coefficients, so the
 * correction can be used in the dynamics prior to finding the more
 * accurate NQC values later on.
 */
UNUSED static int XLALSimIMREOBGetCalibratedNQCCoeffs( 
                                EOBNonQCCoeffs *coeffs, /**<< OUTPUT, Structure for NQC coeffs */
                                INT4            l,      /**<< Mode l */
                                INT4            m,      /**<< Mode m */
                                REAL8           eta     /**<< Symmetric mass ratio */
                                )
{

#ifndef LAL_NDEBUG
  if ( !coeffs )
  {
    XLAL_ERROR( XLAL_EINVAL );
  }
#endif

  if ( l != 2 || m != 2 )
  {
    XLALPrintError( "Mode %d,%d is not supported by this function.\n", l, m );
    XLAL_ERROR( XLAL_EINVAL );
  }

  /* All NQC coefficients are set to zero here */
  /* including coeffs->a3S, coeffs->a4 and coeffs->a5 that are not used in EOBNRv2 */
  memset( coeffs, 0, sizeof( *coeffs ) );

  coeffs->a1 = -4.55919 + 18.761 * eta - 24.226 * eta*eta;
  coeffs->a2 = 37.683 - 201.468 * eta + 324.591 * eta*eta;
  coeffs->a3 = - 39.6024 + 228.899 * eta - 387.222 * eta * eta;

  return XLAL_SUCCESS;
}

/**
 * This function calculates the non-quasicircular correction to apply to
 * the waveform. The form of this correction can be found in Pan et al,
 * PRD84, 124052(2011), Eq.(22), and also in the DCC document T1100433. Note
 * that when calling this function, the NQC coefficients should already
 * have been pre-computed.
 */
UNUSED static int  XLALSimIMREOBNonQCCorrection(
                      COMPLEX16      * restrict nqc,    /**<< OUTPUT, The NQC correction */
                      REAL8Vector    * restrict values, /**<< Dynamics r, phi, pr, pphi */
                      const REAL8               omega,  /**<< Angular frequency */
                      EOBNonQCCoeffs * restrict coeffs  /**<< NQC coefficients */
                     )

{

  REAL8 rOmega, rOmegaSq;
  REAL8 r, p, sqrtR;

  REAL8 mag, phase;


  r = values->data[0];
  p = values->data[2];

  sqrtR = sqrt(r);

  rOmega = r * omega;
  rOmegaSq = rOmega*rOmega;

  /* In EOBNRv2, coeffs->a3S, coeffs->a4 and coeffs->a5 are set to zero */
  /* through XLALSimIMREOBGetCalibratedNQCCoeffs() */
  /* and XLALSimIMREOBCalculateNQCCoefficients() */
  mag = 1. + (p*p / rOmegaSq) * ( coeffs->a1
     + coeffs->a2 / r + ( coeffs->a3 + coeffs->a3S) / (r*sqrtR)
     + coeffs->a4 / (r*r) + coeffs->a5 / (r*r*sqrtR));

  phase = coeffs->b1 * p / rOmega + p*p*p/rOmega * ( coeffs->b2
     + coeffs->b3 / sqrtR + coeffs->b4 / r );

  *nqc = mag * cos(phase);
  *nqc += I * mag * sin(phase);

  return XLAL_SUCCESS;

}

/**
 * This function computes the coefficients a1, a2, etc. used in the
 * non-quasicircular correction. The details of the calculation of these
 * coefficients are found in the DCC document T1100433.
 */
UNUSED static int XLALSimIMREOBCalculateNQCCoefficients(
                 EOBNonQCCoeffs * restrict coeffs,    /**<< OUTPUT, NQC coefficients */
                 REAL8Vector    * restrict amplitude, /**<< Waveform amplitude, func of time */
                 REAL8Vector    * restrict phase,     /**<< Waveform phase(rad), func of time */
                 REAL8Vector    * restrict q1,        /**<< Function of dynamics (see DCC doc) */
                 REAL8Vector    * restrict q2,        /**<< Function of dynamics (see DCC doc) */
                 REAL8Vector    * restrict q3,        /**<< Function of dynamics (see DCC doc) */
                 REAL8Vector    * restrict p1,        /**<< Function of dynamics (see DCC doc) */
                 REAL8Vector    * restrict p2,        /**<< Function of dynamics (see DCC doc) */
                 INT4                      l,         /**<< Mode l */
                 INT4                      m,         /**<< Mode m */
                 REAL8                     timePeak,  /**<< Time of peak orbital frequency */
                 REAL8                     deltaT,    /**<< Sampling interval */
                 REAL8                     eta        /**<< Symmetric mass ratio */
                 )
{

  UINT4 i;

  int signum;

  REAL8Vector * restrict timeVec = NULL;

  /* Since the vectors we actually want are q etc * A, we will have to generate them here */
  REAL8Vector *q1LM = NULL;
  REAL8Vector *q2LM = NULL;
  REAL8Vector *q3LM = NULL; 

  REAL8 a, aDot, aDDot;
  REAL8 omega, omegaDot;

  REAL8 nra, nraDDot;
  REAL8 nromega, nromegaDot;

  REAL8 nrDeltaT, nrTimePeak;

  /* Stuff for finding numerical derivatives */
  gsl_spline    *spline = NULL;
  gsl_interp_accel *acc = NULL;

  /* Matrix stuff for calculating coefficients */
  gsl_matrix *qMatrix = NULL, *pMatrix = NULL;
  gsl_vector *aCoeff  = NULL, *bCoeff  = NULL;

  gsl_vector *amps = NULL, *omegaVec = NULL;

  gsl_permutation *perm1 = NULL, *perm2 = NULL;

  /* All NQC coefficients are set to zero here */ 
  /* including coeffs->a4 that is not used in EOBNRv2 */
  memset( coeffs, 0, sizeof( EOBNonQCCoeffs ) );

  /* Populate the time vector */
  /* It is okay to assume initial t = 0 */
  timeVec = XLALCreateREAL8Vector( q1->length );
  q1LM    = XLALCreateREAL8Vector( q1->length );
  q2LM    = XLALCreateREAL8Vector( q2->length );
  q3LM    = XLALCreateREAL8Vector( q3->length );

  /* Populate vectors as necessary */
  for ( i = 0; i < timeVec->length; i++ )
  {
    timeVec->data[i] = i * deltaT;
    q1LM->data[i]    = q1->data[i] * amplitude->data[i];
    q2LM->data[i]    = q2->data[i] * amplitude->data[i];
    q3LM->data[i]    = q3->data[i] * amplitude->data[i];
  }

  /* Allocate all the memory we need */
  XLAL_CALLGSL(
    /* a stuff */
    qMatrix = gsl_matrix_alloc( 3, 3 );
    aCoeff  = gsl_vector_alloc( 3 );
    amps    = gsl_vector_alloc( 3 );
    perm1   = gsl_permutation_alloc( 3 );

    /* b stuff */
    pMatrix  = gsl_matrix_alloc( 2, 2 );
    bCoeff   = gsl_vector_alloc( 2 );
    omegaVec = gsl_vector_alloc( 2 );
    perm2    = gsl_permutation_alloc( 2 );
  );

  if ( !qMatrix || !aCoeff || !amps || !pMatrix || !bCoeff || !omegaVec )
  {
    gsl_matrix_free( qMatrix );
    gsl_vector_free( amps );
    gsl_vector_free( aCoeff );
    gsl_permutation_free( perm1 );
    gsl_matrix_free( pMatrix );
    gsl_vector_free( omegaVec );
    gsl_vector_free( bCoeff );
    gsl_permutation_free( perm2 );
    XLALDestroyREAL8Vector( q1LM );
    XLALDestroyREAL8Vector( q2LM );
    XLALDestroyREAL8Vector( q3LM );
    XLALDestroyREAL8Vector( timeVec );
    XLAL_ERROR( XLAL_ENOMEM );
  }

  /* The time we want to take as the peak time depends on l and m */
  /* Calculate the adjustment we need to make here */
  nrDeltaT = XLALSimIMREOBGetNRPeakDeltaT( l, m, eta );
  if ( XLAL_IS_REAL8_FAIL_NAN( nrDeltaT ) )
  {
    XLALDestroyREAL8Vector( q1LM );
    XLALDestroyREAL8Vector( q2LM );
    XLALDestroyREAL8Vector( q3LM );
    XLALDestroyREAL8Vector( timeVec );
    XLAL_ERROR( XLAL_EFUNC );
  }

  nrTimePeak = timePeak + nrDeltaT;

  /* We are now in a position to use the interp stuff to calculate the derivatives we need */
  /* We will start with the quantities used in the calculation of the a coefficients */
  spline = gsl_spline_alloc( gsl_interp_cspline, amplitude->length );
  acc    = gsl_interp_accel_alloc();

  /* Q1 */
  gsl_spline_init( spline, timeVec->data, q1LM->data, q1LM->length );
  gsl_matrix_set( qMatrix, 0, 0, gsl_spline_eval( spline, nrTimePeak, acc ) );
  gsl_matrix_set( qMatrix, 1, 0, gsl_spline_eval_deriv( spline, nrTimePeak, acc ) );
  gsl_matrix_set( qMatrix, 2, 0, gsl_spline_eval_deriv2( spline, nrTimePeak, acc ) );

  /* Q2 */
  gsl_spline_init( spline, timeVec->data, q2LM->data, q2LM->length );
  gsl_interp_accel_reset( acc );
  gsl_matrix_set( qMatrix, 0, 1, gsl_spline_eval( spline, nrTimePeak, acc ) );
  gsl_matrix_set( qMatrix, 1, 1, gsl_spline_eval_deriv( spline, nrTimePeak, acc ) );
  gsl_matrix_set( qMatrix, 2, 1, gsl_spline_eval_deriv2( spline, nrTimePeak, acc ) );

  /* Q3 */
  gsl_spline_init( spline, timeVec->data, q3LM->data, q3LM->length );
  gsl_interp_accel_reset( acc );
  gsl_matrix_set( qMatrix, 0, 2, gsl_spline_eval( spline, nrTimePeak, acc ) );
  gsl_matrix_set( qMatrix, 1, 2, gsl_spline_eval_deriv( spline, nrTimePeak, acc ) );
  gsl_matrix_set( qMatrix, 2, 2, gsl_spline_eval_deriv2( spline, nrTimePeak, acc ) );

  /* Amplitude */
  gsl_spline_init( spline, timeVec->data, amplitude->data, amplitude->length );
  gsl_interp_accel_reset( acc );
  a     = gsl_spline_eval( spline, nrTimePeak, acc );
  aDot  = gsl_spline_eval_deriv( spline, nrTimePeak, acc );
  aDDot = gsl_spline_eval_deriv2( spline, nrTimePeak, acc );

  nra = GetNRPeakAmplitude( l, m, eta );
  nraDDot = GetNRPeakADDot( l, m, eta );

  if ( XLAL_IS_REAL8_FAIL_NAN( nra ) || XLAL_IS_REAL8_FAIL_NAN( nraDDot ) )
  {
    XLALDestroyREAL8Vector( q1LM );
    XLALDestroyREAL8Vector( q2LM );
    XLALDestroyREAL8Vector( q3LM );
    XLALDestroyREAL8Vector( timeVec );
    XLAL_ERROR( XLAL_EFUNC );
  }

  gsl_vector_set( amps, 0, nra - a );
  gsl_vector_set( amps, 1, - aDot );
  gsl_vector_set( amps, 2, nraDDot - aDDot );

  /* We have now set up all the stuff to calculate the a coefficients */
  /* So let us do it! */
  gsl_linalg_LU_decomp( qMatrix, perm1, &signum );
  gsl_linalg_LU_solve( qMatrix, perm1, amps, aCoeff );

  /* Now we (should) have calculated the a values. Now we can do the b values */

  /* P1 */
  gsl_spline_init( spline, timeVec->data, p1->data, p1->length );
  gsl_interp_accel_reset( acc );
  gsl_matrix_set( pMatrix, 0, 0, - gsl_spline_eval_deriv( spline, nrTimePeak, acc ) );
  gsl_matrix_set( pMatrix, 1, 0, - gsl_spline_eval_deriv2( spline, nrTimePeak, acc ) );

  /* P2 */
  gsl_spline_init( spline, timeVec->data, p2->data, p2->length );
  gsl_interp_accel_reset( acc );
  gsl_matrix_set( pMatrix, 0, 1, - gsl_spline_eval_deriv( spline, nrTimePeak, acc ) );
  gsl_matrix_set( pMatrix, 1, 1, - gsl_spline_eval_deriv2( spline, nrTimePeak, acc ) );

  /* Phase */
  gsl_spline_init( spline, timeVec->data, phase->data, phase->length );
  gsl_interp_accel_reset( acc );
  omega    = gsl_spline_eval_deriv( spline, nrTimePeak, acc );
  omegaDot = gsl_spline_eval_deriv2( spline, nrTimePeak, acc );

  /* Since the phase can be decreasing, we need to take care not to have a -ve frequency */
  if ( omega * omegaDot > 0.0 )
  {
    omega    = fabs( omega );
    omegaDot = fabs( omegaDot );
  }
  else
  {
    omega    = fabs( omega );
    omegaDot = - fabs( omegaDot );
  }

  nromega = GetNRPeakOmega( l, m, eta );
  nromegaDot = GetNRPeakOmegaDot( l, m, eta );

  if ( XLAL_IS_REAL8_FAIL_NAN( nromega ) || XLAL_IS_REAL8_FAIL_NAN( nromegaDot ) )
  {
    XLALDestroyREAL8Vector( q1LM );
    XLALDestroyREAL8Vector( q2LM );
    XLALDestroyREAL8Vector( q3LM );
    XLALDestroyREAL8Vector( timeVec );
    XLAL_ERROR( XLAL_EFUNC );
  }

  gsl_vector_set( omegaVec, 0, nromega - omega );
  gsl_vector_set( omegaVec, 1, nromegaDot - omegaDot );

  /* And now solve for the b coefficients */
  gsl_linalg_LU_decomp( pMatrix, perm2, &signum );
  gsl_linalg_LU_solve( pMatrix, perm2, omegaVec, bCoeff );

  /* We can now populate the coefficients structure */
  coeffs->a1 = gsl_vector_get( aCoeff, 0 );
  coeffs->a2 = gsl_vector_get( aCoeff, 1 );
  coeffs->a3 = gsl_vector_get( aCoeff, 2 );
  coeffs->b1 = gsl_vector_get( bCoeff, 0 );
  coeffs->b2 = gsl_vector_get( bCoeff, 1 );

  /* Free memory and exit */
  gsl_matrix_free( qMatrix );
  gsl_vector_free( amps );
  gsl_vector_free( aCoeff );
  gsl_permutation_free( perm1 );

  gsl_matrix_free( pMatrix );
  gsl_vector_free( omegaVec );
  gsl_vector_free( bCoeff );
  gsl_permutation_free( perm2 );

  gsl_spline_free( spline );
  gsl_interp_accel_free( acc );

  XLALDestroyREAL8Vector( q1LM );
  XLALDestroyREAL8Vector( q2LM );
  XLALDestroyREAL8Vector( q3LM );
  XLALDestroyREAL8Vector( timeVec );

  return XLAL_SUCCESS;
}

/* ------------------------------------------------
 *          Spin (SEOBNRv1 and SEOBNRv2)
 * ------------------------------------------------*/

/**
 * The time difference between the orbital peak and the peak amplitude
 * of the mode in question (currently only 2,2 implemented ).
 * Eq. 33 of Taracchini et al. PRD 86, 024011 (2012).
 */
UNUSED static inline REAL8 XLALSimIMREOBGetNRSpinPeakDeltaT( 
                 INT4 l,           /**<< Mode l */
                 INT4 m,           /**<< Mode m */
                 REAL8 UNUSED eta, /**<< Symmetric mass ratio */
                 REAL8 a           /**<< Dimensionless spin */
                 )
{

  switch ( l )
  {
    case 2:
      switch ( m )
      {
        case 2:
          /* DeltaT22 defined here is a minus sign different from Eq. (33) of Taracchini et al. */
          if ( a <= 0.0 )
          {
            return 2.5;
          }
          else
          {
            return (2.5 + 1.77*a*a*a*a/(0.43655*0.43655*0.43655*0.43655)/(1.0-2.0*eta)/(1.0-2.0*eta)/(1.0-2.0*eta)/(1.0-2.0*eta));
          }
          break;
        default:
          XLAL_ERROR_REAL8( XLAL_EINVAL );
      }
      break;
    default:
      XLAL_ERROR_REAL8( XLAL_EINVAL );
  }

  /* We should never get here, but I expect a compiler whinge without it... */
  XLALPrintError( "XLAL Error %s - We should never get here!!\n", __func__ );
  XLAL_ERROR_REAL8( XLAL_EINVAL );
}

/* FIXME: Add XLALSimIMREOB to these function names */

/**
 * Peak amplitude predicted by fitting NR results (currently only 2,2 available).
 * Tables IV and V and Eq. 42 of Taracchini et al. PRD 86, 024011 (2012).
 */
UNUSED static inline REAL8 GetNRSpinPeakAmplitude( INT4 UNUSED l, INT4 UNUSED m, REAL8 UNUSED eta, REAL8 UNUSED a )
{
  /* Fit for HOMs missing */
  return 1.3547468629743946*eta + 0.9187885481024214*eta*eta;
}

/**
 * Peak amplitude curvature predicted by fitting NR results (currently only 2,2 available).
 * Tables IV and V and Eq. 42 of Taracchini et al. PRD 86, 024011 (2012).
 */
UNUSED static inline REAL8 GetNRSpinPeakADDot( INT4 UNUSED l, INT4 UNUSED m, REAL8 UNUSED eta, REAL8 UNUSED a )
{
  /* Fit for HOMs missing */
  return eta*(-0.0024971911410897156 + (-0.006128515435641139 + 0.01732656*a/(2.0-4.0*eta))*eta);
}

/**
 * Peak frequency predicted by fitting NR results (currently only 2,2 available).
 * Tables IV and V and Eq. 42 of Taracchini et al. PRD 86, 024011 (2012).
 */
UNUSED static inline REAL8 GetNRSpinPeakOmega( INT4 UNUSED l, INT4 UNUSED m, REAL8 UNUSED eta, REAL8 a )
{
  /* Fit for HOMs missing */
  return 0.27581190323955274 + 0.19347381066059993*eta
       - 0.08898338208573725*log(1.0 - a/(1.0-2.0*eta))
       + eta*eta*(1.78832*(0.2690779744133912 + a/(2.0-4.0*eta))*(1.2056469070395925
       + a/(2.0-4.0*eta)) + 1.423734113371796*log(1.0 - a/(1.0-2.0*eta)));
}

/**
 * Peak frequency slope predicted by fitting NR results (currently only 2,2 available).
 * Tables IV and V and Eq. 42 of Taracchini et al. PRD 86, 024011 (2012).
 */
UNUSED static inline REAL8 GetNRSpinPeakOmegaDot( INT4 UNUSED l, INT4 UNUSED m, REAL8 UNUSED eta, REAL8 UNUSED a )
{
  /* Fit for HOMs missing */
  return 0.006075014646800278 + 0.012040017219351778*eta
       + (0.0007353536801336875 + 0.0015592659912461832*a/(1.0-2.0*eta))*log(1.0-a/(1.0-2.0*eta))
       + eta*eta*(0.03575969677378844 + (-0.011765658882139 - 0.02494825585993893*a/(1.0-2.0*eta))
       * log(1.0 - a/(1.0-2.0*eta)));
}

/**
 * Peak frequency predicted by fitting NR results (currently only 2,2 available).
 * Take from unpublished SEOBNRv2 results.
 */
UNUSED static inline REAL8 GetNRSpinPeakOmegav2( INT4 UNUSED l, INT4 UNUSED m, REAL8 UNUSED eta, REAL8 a )
{
  /* Fit for HOMs missing */
  return 0.43747541927878864 + (-0.10933208665273314 - 0.007325831113333813 * a/(1.0-2.0*eta)) 
       * log( 2.0803657591886267 - a * (1.376497141999324 - 11.513558950322647 
       * (eta - 0.25) ) - 9.681916048928946 * (eta - 0.25));
}

/**
 * Peak frequency slope predicted by fitting NR results (currently only 2,2 available).
 * Take from unpublished SEOBNRv2 results.
 */
UNUSED static inline REAL8 GetNRSpinPeakOmegaDotv2( INT4 UNUSED l, INT4 UNUSED m, REAL8 UNUSED eta, REAL8 UNUSED a )
{
  REAL8 an = a/(1.0-2.0*eta);
  /* Fit for HOMs missing */
  return -0.000016000256004096065 * (0.011334382589915553 + 0.0002412377539787519 * an) 
       - 0.00099601593625498 * (0.01790260950222435 + 0.0187545997359713 * an) 
       + (0.01790260950222435 + 0.0187545997359713 * an) * eta 
       + 1.0000160002560041 * (-0.011209791668429111 + (0.004086795897856442 
       + 0.0006333925136134383 * an) * log( 68.47466578101876 - 58.30148755701496 * an)) 
       + eta*eta * (16.000256004096066 * (0.011334382589915553 + 0.0002412377539787519 * an) 
       - 3.9840637450199203 * (0.01790260950222435 + 0.0187545997359713 * an) - 16.000256004096066
       * (-0.011209791668429111 + (0.004086795897856442 + 0.0006333925136134383 * an) 
       * log( 68.47466578101876 - 58.30148755701496 * an )));
}

/**
 * Function to interpolate known amplitude NQC coeffcients of spin terms,
 * namely a3s, a4 and a5.
 * The a3s, a4 and a5 values were calculated for
 * 11 mass ratios q=1,1.5,2,3,4,5,6,10,20,50 and 100, and
 * 19 spin (\f$\chi\f$ defined in Taracchini et al. PRD 86, 024011 (2012)) values
 * chi = -1, -0.9, -0.8, ......, 0.3, 0.4, 0.5, 0.55, 0.6, 0.65.
 * The calculation was done by Andrea Taracchini using a C++ code of the UMaryland group.
 * In principle, these numbers can be automatically calculated iteratively by the LAL code.
 * However, since such calcualtion increase the cost of each waveform generation by
 * about an order of magnitude, we prepare these numbers in advance reduce cost.
 * These number can be verified by confirming that
 * the peak amplitude and frequency agree well with the NR-fits predicted values,
 * and to get exact NR-fits predicted values, corrections on these numbers are ~1%.
 */
UNUSED static int XLALSimIMRGetEOBCalibratedSpinNQC( EOBNonQCCoeffs * restrict coeffs, 
                                    INT4 UNUSED l, 
                                    INT4 UNUSED m, 
                                    REAL8 eta, 
                                    REAL8 a )
{
  const unsigned int nsqdim = 100;
  const unsigned int   qdim = 50;
  const unsigned int   adim = 69;
  UINT4 i;
  /* REAL8 eta2 = eta*eta;*/
  REAL8 a3slist[qdim], a4list[qdim], a5list[qdim];

  memset( coeffs, 0, sizeof( *coeffs ) );
  const double nsetalist[100] = {0.0025, 0.005, 0.0075, 0.01, 0.0125, 0.015, 0.0175, 0.02, 0.0225, 0.025, 0.0275, 0.03, 0.0325, 0.035, 0.0375, 0.04, 0.0425, 0.045, 0.0475, 0.05, 0.0525, 0.055, 0.0575, 0.06, 0.0625, 0.065, 0.0675, 0.07, 0.0725, 0.075, 0.0775, 0.08, 0.0825, 0.085, 0.0875, 0.09, 0.0925, 0.095, 0.0975, 0.1, 0.1025, 0.105, 0.1075, 0.11, 0.1125, 0.115, 0.1175, 0.12, 0.1225, 0.125, 0.1275, 0.13, 0.1325, 0.135, 0.1375, 0.14, 0.1425, 0.145, 0.1475, 0.15, 0.1525, 0.155, 0.1575, 0.16, 0.1625, 0.165, 0.1675, 0.17, 0.1725, 0.175, 0.1775, 0.18, 0.1825, 0.185, 0.1875, 0.19, 0.1925, 0.195, 0.1975, 0.2, 0.2025, 0.205, 0.2075, 0.21, 0.2125, 0.215, 0.2175, 0.22, 0.2225, 0.225, 0.2275, 0.23, 0.2325, 0.235, 0.2375, 0.24, 0.2425, 0.245, 0.2475, 0.25};
  const double etalist[50] = {0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.065, 0.07, 0.075, 0.08, 0.085, 0.09, 0.095, 0.1, 0.105, 0.11, 0.115, 0.12, 0.125, 0.13, 0.135, 0.14, 0.145, 0.15, 0.155, 0.16, 0.165, 0.17, 0.175, 0.18, 0.185, 0.19, 0.195, 0.2, 0.205, 0.21, 0.215, 0.22, 0.225, 0.23, 0.235, 0.24, 0.245, 0.25};
  const double alist[69]   = {-1., -0.975, -0.95, -0.925, -0.9, -0.875, -0.85, -0.825, -0.8, -0.775, -0.75, -0.725, -0.7, -0.675, -0.65, -0.625, -0.6, -0.575, -0.55, -0.525, -0.5, -0.475, -0.45, -0.425, -0.4, -0.375, -0.35, -0.325, -0.3, -0.275, -0.25, -0.225, -0.2, -0.175, -0.15, -0.125, -0.1, -0.075, -0.05, -0.025, 0., 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475, 0.5, 0.525, 0.55, 0.575, 0.6, 0.625, 0.65, 0.675, 0.7};

  const double a1list[100] = {-21.648106754277453, -20.47585498319721, -19.530321799953853, -18.760004570030972, -18.079213572570588, -17.419521042795658, -16.859312290615286, -16.28756533894047, -15.772560753656803, -15.303367029222859, -14.847798623544204, -14.404050939670224, -13.976897526988164, -13.578815297784823, -13.200763216702892, -12.841976980076495, -12.501692284239752, -12.179144825526786, -11.865421415437213, -11.561721630442458, -11.267369305715642, -10.980678895252732, -10.699964853049693, -10.42829543384451, -10.167524207712544, -9.914674799750443, -9.669372681081281, -9.43124332282813, -9.199912196114067, -8.97500477206216, -8.756146521795488, -8.542962916437121, -8.335079427110134, -8.130867290476582, -7.9313984345769715, -7.736518055359128, -7.545988962963793, -7.35957396753171, -7.177035879203623, -6.998068318971551, -6.8224288840577945, -6.649961333309703, -6.480441962265538, -6.313647066463564, -6.149352941442041, -5.9873358827392344, -5.827372185893404, -5.669238146442815, -5.512513869477759, -5.347444298469859, -5.183197664196798, -5.019955297249104, -4.857898528217308, -4.697208687691938, -4.538067106263523, -4.380655114522592, -4.225154043059675, -4.071745222465302, -3.9206099833299994, -3.771929656244297, -3.6258855717987246, -3.4826590605838113, -3.342431453190086, -3.2053840802080775, -3.0798495456033685, -2.957143501994997, -2.837141996613968, -2.7197210766912865, -2.6047567894579564, -2.4921251821449832, -2.381702301983373, -2.273364196204128, -2.1669869120382543, -2.0624464967167566, -1.9596189974706397, -1.8495581080328225, -1.7409749049276515, -1.6340876818580403, -1.5291147325269023, -1.4262743506371505, -1.3257848298916983, -1.2278644639934593, -1.1327315466453465, -1.0406043715502744, -0.9517012324111543, -0.8662404229309005, -0.7844402368124264, -0.7065189677586451, -0.6326949094724702, -0.5631863556568149, -0.4982116000145925, -0.4379889362487165, -0.3827366580621, -0.3326730591576565, -0.28801643323829923, -0.24898507400694198, -0.21579727516649713, -0.18867133041987855, -0.1678255334699995, -0.15347817801977337};

  const double a2list[100] = {178.73204288078207, 167.7345170263427, 158.85457776976878, 151.63702661020528, 145.25886103554546, 139.08071349079492, 133.85300186061994, 128.49833153024582, 123.69581397508206, 119.34304449990263, 115.1225036726121, 111.02138715464824, 107.08437843884623, 103.42908442981195, 99.96941370234742, 96.69762600347904, 93.6059810802332, 90.68673867963636, 87.86042720725196, 85.13659365364732, 82.50827361644501, 79.95961751850369, 77.47477578268207, 75.08044751075573, 72.79299640785891, 70.58542621544011, 68.45392194344085, 66.39466860180258, 64.40385120046682, 62.47765474937501, 60.61226425846867, 58.80386473768926, 57.04864119697828, 55.32723154023535, 53.65375061411054, 52.0271070701423, 50.44518833065014, 48.905881817953585, 47.407074954372156, 45.94482445612376, 44.51687981462722, 43.123270309615144, 41.76224017946817, 40.43203366256695, 39.13089499729211, 37.85706842202429, 36.60879817514412, 35.384328495032264, 34.180608559988464, 32.932965816541305, 31.701640999519693, 30.48755199799916, 29.291616701055233, 28.11475299776345, 26.95787877719934, 25.821911928438432, 24.70777034055626, 23.616371902628366, 22.54863450373026, 21.50547603293748, 20.487814379325563, 19.496567431970035, 18.53265307994643, 17.596989212330282, 16.733055464585682, 15.895477726977264, 15.083579856295431, 14.296685709330593, 13.53411914287315, 12.795204013713512, 12.079264178642092, 11.385623494449275, 10.71360581792548, 10.06253500586111, 9.431734915046569, 8.755551062926262, 8.098379005319565, 7.462063310457031, 6.848448546569211, 6.259379281886656, 5.69670008463992, 5.162255523059555, 4.657890165376111, 4.1854485798201475, 3.7467753346222055, 3.3437149980128416, 2.978112138222609, 2.651811323482059, 2.366657122021744, 2.124494102072216, 1.9271668318640272, 1.77651987962773, 1.674397813593876, 1.6226452019930173, 1.6231066130557064, 1.6776266150124943, 1.7880497760939342, 1.956220664530578, 2.1839838485529777, 2.4731838963916855};

  const double a3list[100] = {-198.6684740510964, -185.71560983427335, -175.26102024642407, -166.7798147654891, -159.28986302859136, -152.03810090459336, -145.91734205946045, -139.6355596843493, -134.01653939524905, -128.93945143613453, -124.02034829337536, -119.24645127356649, -114.67037356959001, -110.43118685036107, -106.4264567068743, -102.6466423946065, -99.08220316903453, -95.72359828563523, -92.4798228915885, -89.36096908139594, -86.35837723133453, -83.45329697039512, -80.62697792756852, -77.90968987612793, -75.32016588735331, -72.82709457524146, -70.42579806564999, -68.11159848443651, -65.87981795745866, -63.72577861057402, -61.64480256964021, -59.63221196051487, -57.683328909055604, -55.7731074201788, -53.92035933426217, -52.12397499502319, -50.38150684315194, -48.69050731933857, -47.04852886427319, -45.45070004394427, -43.89439067703587, -42.37998904664536, -40.90552002179141, -39.46900847149271, -38.068479264767944, -36.70195727063581, -35.36746735811498, -34.063034396224154, -32.78536874375959, -31.468633548362497, -30.173781283368523, -28.901550623295833, -27.652680242662605, -26.427908815987003, -25.227975017787198, -24.053617522581362, -22.905575004887663, -21.78458613922428, -20.691389600109368, -19.6267240620611, -18.59132819959765, -17.585940687237187, -16.61130019949788, -15.668145410897901, -14.794309398041404, -13.950185357101319, -13.135122694543167, -12.348470816832473, -11.589579130434757, -10.857797041815543, -10.152473957440362, -9.47295928377472, -8.818602427284151, -8.188752794434174, -7.582759791690314, -6.932160232642741, -6.304213548518645, -5.700899806714059, -5.124199074625014, -4.576091419647542, -4.0585569091776765, -3.573575610611449, -3.123127591344892, -2.709192918774041, -2.3337516602949204, -1.9987838833035663, -1.706269655196011, -1.458189043368287, -1.2565221152164259, -1.1032489381364605, -1.0003495795244224, -0.9498041067763443, -0.953592587288258, -1.0136950884561957, -1.1320916776761898, -1.31076242234427, -1.5516873898564725, -1.8568466476088277, -2.2282202629973678, -2.667788303418125};

  const double b1list[100] = {-0.5693500504085347, -0.576434151837257, -0.5827940588889807, -0.5875005969333106, -0.5915255507494274, -0.5970658827548452, -0.6057016775611604, -0.6053160995270499, -0.6070988602490128, -0.6110941958474475, -0.6140262971912503, -0.6172989788502661, -0.6206089494421659, -0.6234488672149939, -0.6258813681301192, -0.6277498501118401, -0.628897711084455, -0.6291683489722621, -0.6308371388180571, -0.6331613723912005, -0.6359581432020665, -0.6393457874435003, -0.6434426413083475, -0.6478966252652467, -0.6524003236766366, -0.6571425380351175, -0.6620695856777987, -0.6671277839417898, -0.6722634501642004, -0.6774229016821397, -0.6825524558327174, -0.6875984299530429, -0.6925071413802255, -0.6961049747777461, -0.6996298196832709, -0.7032241904115, -0.7069570372346382, -0.7108973104248908, -0.7151139602544623, -0.720064149835206, -0.7258460785616149, -0.7320745022572427, -0.7387427060855165, -0.7458439752098636, -0.7533715947937112, -0.7613188500004865, -0.7696790259936167, -0.778445407936529, -0.7876035918776259, -0.7967733093076456, -0.8063159350249182, -0.8162406283658405, -0.8265565486668099, -0.8372728552642233, -0.8483987074944775, -0.8599432646939698, -0.871915686199097, -0.884325131346256, -0.8971807594718443, -0.9104917299122585, -0.9242672020038958, -0.938516335083153, -0.9532482884864274, -0.9684722215501159, -0.9837886712042127, -0.9996512322126136, -1.0160843677317775, -1.0331125409181625, -1.0507602149282274, -1.0690518529184303, -1.08801191804523, -1.1076648734650847, -1.1280351823344532, -1.1491473078097938, -1.171025713047565, -1.1939017528688474, -1.2175927015073384, -1.2421149961497486, -1.2674850739827876, -1.2937193721931661, -1.320834327967594, -1.3488463784927816, -1.377771960955439, -1.4076275125422764, -1.4384294704400042, -1.4701942718353325, -1.5029383539149717, -1.5366781538656316, -1.5714301088740226, -1.6072106561268549, -1.6440362328108387, -1.6819232761126843, -1.7208882232191016, -1.7609475113168012, -1.802117577592493, -1.8444148592328868, -1.8878557934246936, -1.9324568173546235, -1.9782343682093864, -2.0252048831756926};

  const double b2list[100] = {1.6514745488753086, 1.6733593678301482, 1.687838328174986, 1.7031979370575185, 1.712831020475929, 1.7266279186283089, 1.7581796869631672, 1.7499867318965905, 1.7518398412177276, 1.7634468469918447, 1.770740685014047, 1.779639998248617, 1.788893228893931, 1.7964389585725973, 1.8024779675983216, 1.8063408969988246, 1.8073583878018264, 1.8048610810350476, 1.8077536017385523, 1.8130592620946404, 1.8200051503317694, 1.8290042591854243, 1.8404695813910905, 1.8531136802761718, 1.8658884601302625, 1.8795206835759333, 1.8938451458175296, 1.908696642059398, 1.923909967505884, 1.9393199173613336, 1.9547612868300928, 1.9700688711165077, 1.9850774654249244, 1.9960816492286722, 2.0069990026439246, 2.0182845185376195, 2.0301606467081585, 2.042849836953944, 2.056574539073377, 2.0726094682881886, 2.0912563552368337, 2.111506538675078, 2.1333773752936596, 2.1568862217833127, 2.1820504348347747, 2.208887371138781, 2.2374143873860683, 2.2676488402673725, 2.299578947658098, 2.3318064402049257, 2.3657435940219895, 2.4014679228653333, 2.439056940491002, 2.478588160655038, 2.520139097113486, 2.56378726362239, 2.609610173937794, 2.657685341815741, 2.708090281012276, 2.760902505283443, 2.8161995283852854, 2.874058864073847, 2.9345580261051722, 2.9977745282353045, 3.0615596376827443, 3.128412230361185, 3.1984931979224718, 3.271963432018449, 3.3489838243009618, 3.429715266421855, 3.5143186500329735, 3.602954866786163, 3.6957848083332676, 3.792969366326133, 3.894669432416603, 4.002044979667201, 4.114256382857298, 4.23142577617022, 4.353675293789295, 4.48112706989785, 4.613903238679211, 4.752125934316707, 4.895917290993665, 5.04539944289341, 5.200694524199272, 5.361924669094577, 5.529212011762653, 5.702678686386827, 5.882446827150425, 6.068638568236777, 6.261376043829206, 6.460781388111044, 6.666976735265615, 6.880084219476247, 7.100225974926268, 7.327524135799001, 7.5621008362777795, 7.804078210545929, 8.053578392786775, 8.310723517183646};

  const double a3stab[50][69] = {
{1298.87913, 1261.29637, 1227.72908, 1197.25947, 1168.96973, 1141.94207, 1115.2587, 1088.00181, 1059.25362, 1022.90618, 985.307903, 947.617046, 910.991876, 881.346609, 853.181173, 825.751447, 798.313312, 767.87765, 736.843339, 705.36426, 673.594293, 642.125047, 610.497581, 578.690684, 546.683144, 514.497622, 482.051486, 449.305976, 416.222332, 382.495195, 348.459044, 314.181759, 279.731222, 245.411071, 210.959123, 176.348956, 141.554145, 106.628243, 71.43286, 35.9095836, 0., -36.1519021, -72.887899, -110.347365, -148.669675, -186.733245, -226.44279, -268.442067, -313.374835, -358.91977, -409.871742, -468.06054, -535.315953, -602.768787, -687.227409, -794.8012, -931.599543, -1073.78457, -1269.39181, -1536.50955, -1893.22609, -2326.65745, -2947.80866, -3741.68954, -4974.29064, -6365.10282, -9538.63496, -15643.1414, -25826.8766}, 
{1196.32002, 1167.06104, 1137.27475, 1107.12301, 1076.76769, 1046.37066, 1016.09379, 986.098937, 956.547978, 928.8469, 901.415801, 873.918898, 846.020409, 814.568078, 783.16919, 752.614552, 723.694973, 702.877684, 683.006503, 662.601671, 640.183429, 607.718815, 572.902552, 536.876163, 500.78117, 469.351999, 438.700105, 408.529845, 378.545579, 347.703945, 316.75611, 285.705522, 254.555628, 223.39562, 192.108901, 160.664621, 129.031928, 97.2453174, 65.1824526, 32.7863434, 0., -33.0547348, -66.6852171, -101.01997, -136.187517, -170.971669, -207.383546, -246.089556, -287.756108, -330.313665, -378.258958, -433.352771, -497.35589, -562.191692, -643.393334, -746.656564, -877.67713, -1013.79966, -1200.41146, -1454.54875, -1793.24772, -2204.40576, -2792.4755, -3542.23306, -4707.36116, -6024.93464, -9029.55323, -14806.4354, -24440.7998}, 
{1111.34089, 1088.04208, 1061.93086, 1033.65502, 1003.86233, 973.200551, 942.317471, 911.86086, 882.478492, 858.393797, 835.248631, 812.260507, 788.646935, 759.48951, 729.796029, 700.438372, 672.288418, 651.387883, 631.370875, 611.041341, 589.203225, 559.502697, 527.964588, 495.455951, 462.843844, 433.627925, 404.989603, 376.74289, 348.701797, 320.308767, 291.898012, 263.432172, 234.873889, 206.130982, 177.242843, 148.194039, 118.969142, 89.6533759, 60.0903922, 30.2244979, 0., -30.4005164, -61.3663335, -93.0484555, -125.597886, -157.893074, -191.866601, -228.178495, -267.488782, -307.871972, -353.607815, -406.390546, -467.914398, -530.2754, -608.605274, -708.437534, -835.305696, -967.651586, -1148.93709, -1395.53239, -1723.80769, -2122.58454, -2690.87904, -3414.09589, -4533.73512, -5812.62405, -8672.18555, -14119.1936, -23160.4223}, 
{1039.53251, 1020.2822, 997.603313, 972.198968, 944.772298, 916.02643, 886.664491, 857.389611, 828.904917, 804.925294, 781.93741, 759.439688, 736.930555, 712.636669, 687.836926, 662.538457, 636.748394, 609.788544, 582.625489, 555.540492, 528.814811, 504.169828, 479.870634, 455.622441, 431.13046, 404.996091, 378.469884, 351.698574, 324.8289, 298.520437, 272.201947, 245.815032, 219.301294, 192.43116, 165.385874, 138.17551, 110.810137, 83.4601202, 55.9111218, 28.1090968, 0., -28.1232207, -56.8004007, -86.2243824, -116.588008, -146.981092, -179.140716, -213.700934, -251.2958, -290.060701, -334.127824, -385.130693, -444.702827, -504.716717, -580.471327, -677.504592, -801.354447, -931.385611, -1109.77852, -1352.54038, -1675.67843, -2068.96254, -2627.11195, -3337.00139, -4427.97475, -5693.23265, -8428.0796, -13552.2941, -21985.6547}, 
{976.485655, 959.824144, 940.197989, 918.098323, 894.016279, 868.442993, 841.869596, 814.787222, 787.687004, 761.819814, 736.613152, 712.254257, 688.930369, 670.53532, 652.067116, 632.23036, 609.729651, 574.459963, 537.459377, 500.956344, 467.179318, 448.413342, 432.807639, 418.568025, 403.900314, 381.127188, 356.690849, 331.150363, 305.064801, 280.549729, 255.985115, 231.307428, 206.453135, 181.110154, 155.562924, 129.847331, 103.999264, 78.2732522, 52.3990841, 26.3251899, 0., -26.1568213, -52.856571, -80.3393122, -108.845108, -137.719354, -168.454647, -201.64892, -237.900101, -275.345863, -318.028498, -367.530041, -425.432527, -483.212445, -556.599593, -651.218222, -772.692584, -901.046982, -1077.7456, -1318.65266, -1639.63241, -2031.30849, -2585.26688, -3290.67293, -4364.64226, -5631.82205, -8258.783, -13076.6149, -20916.4077}, 
{917.791081, 902.710621, 885.620781, 866.696557, 846.112943, 824.044936, 800.66753, 776.155722, 750.684506, 722.455778, 694.406874, 667.502029, 642.705478, 629.711227, 617.261835, 602.829629, 583.886937, 541.782439, 496.561572, 452.146121, 412.457876, 398.926373, 390.962551, 385.485097, 379.412698, 359.691907, 337.202401, 312.851722, 287.547412, 264.607418, 241.564715, 218.362682, 194.944699, 170.981969, 146.79892, 122.447804, 97.9808711, 73.7004736, 49.308722, 24.7578271, 0., -24.435292, -49.4039968, -75.1848062, -102.056412, -129.591492, -159.057156, -191.014502, -226.024627, -262.193469, -303.519346, -351.545419, -407.814847, -463.459392, -534.598171, -626.938905, -746.189313, -872.680952, -1047.64817, -1286.94915, -1606.44207, -1997.39115, -2549.4365, -3254.83387, -4318.29988, -5593.45384, -8125.84343, -12663.0342, -19952.5918}, 
{859.039559, 844.98436, 829.777581, 813.337142, 795.580959, 776.426953, 755.793041, 733.597141, 709.757173, 680.211609, 650.449591, 621.980817, 596.314984, 586.690157, 578.196318, 567.651816, 551.875002, 508.136269, 460.621105, 413.967044, 372.811617, 362.402057, 358.522317, 357.566049, 355.926907, 338.360939, 317.554443, 294.556114, 270.414647, 248.904279, 227.257944, 205.434118, 183.391278, 160.860605, 138.11879, 115.215226, 92.1993096, 69.3494866, 46.3944783, 23.2920582, 0., -22.8926065, -46.3118306, -70.5524255, -95.9091441, -122.081137, -150.197001, -180.789729, -214.392317, -249.069531, -288.809882, -335.133657, -389.561139, -443.154361, -512.075161, -602.027124, -718.713836, -842.332775, -1014.29608, -1250.50977, -1566.87987, -1954.97924, -2503.71346, -3209.20757, -4263.50982, -5543.18962, -7990.80852, -12282.43, -19094.1177}, 
{795.821854, 782.688084, 768.57428, 753.36355, 736.939001, 719.18374, 699.980874, 679.213512, 656.764759, 628.465728, 599.872319, 572.488437, 547.817986, 537.997872, 529.6458, 520.012472, 506.348594, 469.901748, 430.32701, 391.276333, 356.401671, 345.533527, 339.673884, 336.003275, 331.702235, 314.804975, 295.296878, 274.017005, 251.804419, 231.651088, 211.382003, 190.97506, 170.40816, 149.560065, 128.547461, 107.387899, 86.098929, 64.8279929, 43.4107958, 21.8129331, 0., -21.4627385, -43.4492248, -66.2337315, -90.0905311, -114.671921, -141.12294, -169.966649, -201.726112, -234.440059, -272.109617, -316.251583, -368.382751, -419.994158, -486.638661, -573.843361, -687.135356, -806.047703, -972.49916, -1202.41444, -1511.71827, -1891.84153, -2432.19041, -3133.5174, -4174.83432, -5446.09099, -7815.22592, -11905.6805, -18340.896}, 
{723.728735, 711.864516, 697.916767, 682.119252, 664.705738, 645.909991, 625.965777, 605.106863, 583.567013, 560.596558, 537.806075, 515.822704, 495.273586, 480.160139, 466.385516, 453.227147, 439.962463, 423.459176, 406.36832, 388.931212, 371.389169, 355.013918, 338.6042, 321.989168, 304.997976, 286.694707, 267.979608, 248.987859, 229.854639, 211.058618, 192.254089, 173.438833, 154.610635, 135.894352, 117.109862, 98.2041202, 79.1240795, 59.7436944, 40.1121174, 20.2055015, 0., -20.0796619, -40.6853319, -62.0202854, -84.2877981, -106.847476, -131.083733, -157.53731, -186.748953, -216.771064, -251.628063, -292.856026, -341.991034, -391.675587, -455.896772, -539.748099, -648.323078, -759.870989, -917.067264, -1135.7431, -1431.72971, -1795.74674, -2318.96003, -3007.48672, -4026.8356, -5267.21955, -7560.6433, -11503.6636, -17692.8373}, 
{670.702392, 658.608735, 645.231211, 630.627242, 614.85425, 597.969656, 580.030884, 561.095353, 541.220486, 518.602027, 495.903746, 473.927738, 453.476095, 438.689549, 425.696101, 413.962387, 402.955046, 391.824136, 380.479507, 368.51443, 355.522175, 339.17332, 321.752905, 303.623275, 285.146778, 267.626973, 250.10851, 232.577251, 215.019058, 197.328079, 179.618576, 161.913096, 144.234187, 126.856941, 109.450342, 91.9359207, 74.2352057, 56.1190544, 37.7199375, 19.0196532, 0., -18.9413868, -38.4048807, -58.5750179, -79.636335, -100.928376, -123.818666, -148.82974, -176.48413, -204.795257, -237.798413, -277.019773, -323.985518, -371.643004, -433.528758, -514.600485, -619.815891, -726.601121, -878.458066, -1091.35706, -1381.26842, -1738.68566, -2256.0096, -2944.21016, -3956.54775, -5195.43966, -7424.38522, -11191.6133, -17045.3526}, 
{623.046755, 610.336599, 597.375254, 583.988785, 570.003256, 555.244732, 539.539276, 522.712953, 504.591828, 481.903784, 458.812338, 436.382826, 415.680585, 401.790876, 390.151139, 380.218743, 371.451052, 363.730969, 355.920111, 347.30563, 337.174678, 321.458633, 304.142732, 285.856435, 267.229206, 250.700686, 234.366085, 218.130789, 201.900189, 185.228258, 168.512366, 151.798468, 135.13252, 118.913825, 102.69365, 86.3766111, 69.8673234, 52.8685906, 35.5675654, 17.9495883, 0., -17.9060924, -36.3290217, -55.439354, -75.4076554, -95.5838756, -117.287444, -141.017173, -167.271875, -194.031478, -225.321234, -262.64751, -307.516673, -353.072705, -412.529312, -490.737815, -592.549535, -694.603535, -841.248297, -1048.62004, -1332.855, -1684.25015, -2196.45943, -2885.12536, -3891.89364, -5131.56784, -7300.31693, -10898.9685, -16428.3503}, 
{578.777874, 565.343763, 552.617764, 540.222136, 527.779136, 514.911024, 501.240056, 486.388491, 469.978587, 447.525103, 424.400795, 401.870923, 381.200742, 368.162933, 357.712364, 349.311322, 342.422098, 337.556737, 332.707869, 326.917882, 319.229161, 304.187697, 287.130832, 268.899514, 250.334688, 234.780471, 219.573371, 204.553063, 189.559225, 173.868813, 158.109315, 142.345496, 126.64212, 111.492196, 96.360951, 81.1418536, 65.7283731, 49.7810319, 33.5194242, 16.9301977, 0., -16.9275463, -34.3655585, -52.4701786, -71.3975492, -90.5102836, -111.075464, -133.566646, -158.457381, -183.715399, -213.322409, -248.754293, -291.486936, -334.842038, -391.711338, -466.832391, -564.942754, -661.929092, -802.920208, -1004.19402, -1282.02842, -1626.52282, -2132.4907, -2819.95411, -3819.61549, -5057.44347, -7167.66153, -10607.4294, -15833.9069}, 
{537.705661, 523.592259, 510.97538, 499.323946, 488.10688, 476.793102, 464.851536, 451.751104, 436.960726, 415.185168, 392.563172, 370.469324, 350.278208, 338.111834, 328.698393, 321.5135, 316.032768, 313.155404, 310.363995, 306.564723, 300.663765, 286.44217, 269.9813, 252.237386, 234.16666, 219.614044, 205.491601, 191.600088, 177.74026, 163.024266, 148.21691, 133.394391, 118.632905, 104.472791, 90.3404513, 76.1264269, 61.7212598, 46.7790037, 31.5212837, 15.9332368, 0., -15.974157, -32.4508828, -49.5726937, -67.4821057, -85.5612143, -105.017124, -126.29652, -149.846085, -173.651974, -201.605614, -235.1379, -275.67973, -316.701541, -370.778872, -442.526802, -536.56041, -628.112094, -762.932689, -957.39035, -1227.85323, -1564.25632, -2062.26726, -2745.81353, -3735.46447, -4967.31563, -7018.71577, -10306.7003, -15248.3048}, 
{499.640027, 485.044119, 472.464742, 461.290869, 450.911477, 440.715539, 430.092032, 418.42993, 405.118208, 384.603164, 363.193521, 342.255325, 323.154622, 311.94369, 303.427849, 297.098651, 292.447647, 290.380935, 288.409705, 285.45969, 280.456626, 267.30371, 251.958228, 235.354929, 218.428561, 204.949119, 191.88201, 179.027883, 166.187391, 152.469136, 138.642637, 124.785364, 110.97479, 97.7363464, 84.5203574, 71.2251094, 57.7488886, 43.7851316, 29.5189138, 14.9304613, 0., -15.0143326, -30.5213867, -46.652101, -63.5374142, -80.5902821, -98.9468193, -119.025158, -141.243429, -163.646157, -189.974524, -221.596107, -259.878479, -298.401753, -349.435951, -417.463636, -506.967366, -592.686844, -720.744632, -907.520433, -1169.39395, -1496.20329, -1983.95294, -2659.82073, -3635.19172, -4855.43337, -6845.77638, -9986.48562, -14657.826}, 
{464.390881, 449.661375, 437.102489, 426.119557, 416.117917, 406.502905, 396.679857, 386.054108, 374.030995, 355.498276, 336.185896, 317.306223, 300.071624, 289.964612, 282.219351, 276.34015, 271.831319, 269.087297, 266.366212, 262.816322, 257.585886, 245.853978, 232.325711, 217.737017, 202.823827, 190.533415, 178.505829, 166.592464, 154.644712, 141.977941, 129.193979, 116.358629, 103.537694, 91.1635985, 78.7888757, 66.3326796, 53.7141645, 40.722041, 27.4580845, 13.8936269, 0., -14.0164816, -28.5134619, -43.6136023, -59.439564, -75.4511011, -92.6989449, -111.57092, -132.45485, -153.502899, -178.232816, -207.926688, -243.866605, -279.69321, -327.386613, -391.285482, -475.728482, -555.187646, -675.814927, -853.895645, -1105.71512, -1421.11637, -1895.71161, -2559.09281, -3514.54842, -4716.04577, -6641.14011, -9636.48968, -14048.7528}, 
{431.768135, 417.406059, 404.905261, 393.806663, 383.651192, 373.979771, 364.333324, 354.252775, 343.279049, 327.589687, 311.434349, 295.699313, 281.270856, 272.480712, 265.391519, 259.511373, 254.348369, 249.128454, 243.754732, 237.848156, 231.029682, 221.17463, 210.347843, 198.868528, 187.055894, 176.114645, 165.124294, 154.049849, 142.856319, 131.325202, 119.678423, 107.954395, 96.1915321, 84.6352843, 73.0342128, 61.3439161, 49.5199927, 37.5123573, 25.2845654, 12.794489, 0., -12.9490123, -26.3635003, -40.3623994, -55.0646445, -69.9972857, -86.1078972, -103.752168, -123.285788, -143.027154, -166.184163, -193.927421, -227.427532, -260.326451, -304.334894, -363.634929, -442.40862, -515.148803, -627.602467, -795.827371, -1035.88127, -1337.74821, -1795.70711, -2440.7469, -3369.28571, -4543.40188, -6397.10368, -9246.41689, -13407.3673}, 
{397.382476, 384.161063, 372.406617, 361.789178, 351.978785, 342.645477, 333.459295, 324.090277, 314.208462, 301.130105, 287.820545, 274.891334, 262.954027, 255.327872, 248.833648, 242.999831, 237.354894, 230.264221, 222.884616, 215.20979, 207.233455, 198.951905, 190.355236, 181.436127, 172.187257, 162.452347, 152.432615, 142.180322, 131.747728, 121.278447, 110.696847, 100.018645, 89.2595619, 78.5197505, 67.6967231, 56.7724257, 45.7288044, 34.6183092, 23.3241805, 11.8001626, 0., -11.943738, -24.3276619, -37.2695573, -50.8872099, -64.7777681, -79.7879094, -96.2436743, -114.471103, -132.961668, -154.609806, -180.475383, -211.618268, -241.71635, -282.164266, -336.974677, -410.160241, -476.397825, -580.770198, -739.024339, -966.907224, -1254.16741, -1694.54713, -2321.20376, -3222.44556, -4362.79135, -6146.7668, -8868.89588, -12823.7025}, 
{359.974686, 348.352, 338.125236, 328.975588, 320.58425, 312.632416, 304.801281, 296.772038, 288.225882, 276.655174, 264.805473, 253.233508, 242.496004, 235.582587, 229.643928, 224.263593, 219.02515, 212.45994, 205.624648, 198.523735, 191.16166, 183.591405, 175.749496, 167.620982, 159.190914, 150.290229, 141.119734, 131.726121, 122.156087, 112.539351, 102.806371, 92.9706295, 83.0456113, 73.124859, 63.1097732, 52.9818134, 42.7224396, 32.3549175, 21.8021786, 11.0289603, 0., -11.1220613, -22.6495511, -34.6968934, -47.3785121, -60.3411594, -74.354, -89.7185264, -106.736232, -124.027317, -144.247083, -168.369539, -197.368694, -225.064687, -262.446946, -313.351026, -381.612486, -442.034343, -539.097709, -688.251154, -904.943253, -1178.0799, -1602.7377, -2215.93586, -3095.8159, -4202.28901, -5927.06542, -8553.90562, -12366.5701}, 
{324.284553, 314.363487, 305.717204, 298.046437, 291.051921, 284.434388, 277.894571, 271.133205, 263.851022, 253.758508, 243.342743, 233.100562, 223.528797, 217.260142, 211.801225, 206.794534, 201.882558, 195.779693, 189.427757, 182.840476, 176.031575, 169.10102, 161.9418, 154.533147, 146.854291, 138.731414, 130.358013, 121.774536, 113.021431, 104.2128, 95.2859787, 86.2519542, 77.1217145, 67.9818174, 58.7374528, 49.3693806, 39.8583605, 30.1989457, 20.3525856, 10.2945228, 0., -10.342512, -21.0562901, -32.2513829, -44.0378392, -56.1058801, -69.1533129, -83.4581171, -99.2982726, -115.426321, -134.255855, -156.675029, -183.572, -208.924082, -243.294606, -290.336061, -353.700939, -408.371479, -498.138522, -638.12266, -843.444484, -1102.14389, -1510.58355, -2109.2871, -2967.15339, -4038.77567, -5707.59445, -8249.83842, -11941.7363}, 
{290.624535, 282.418333, 275.336854, 269.105352, 263.449083, 258.093303, 252.763265, 247.184226, 241.08144, 232.406317, 223.367496, 214.39977, 205.937932, 200.246338, 195.198394, 190.497067, 185.845325, 180.150535, 174.229506, 168.103447, 161.793567, 155.43567, 148.890531, 142.133519, 135.140006, 127.739783, 120.11203, 112.290346, 104.308333, 96.2624774, 88.09834, 79.8243672, 71.4490055, 63.0513927, 54.5410077, 45.8980203, 37.1026006, 28.1249184, 18.959144, 9.58944784, 0., -9.60482061, -19.549136, -29.9368594, -40.8719042, -52.0789138, -64.1927792, -77.4691216, -92.1635624, -107.160693, -124.635576, -145.392245, -170.234732, -193.315202, -224.750304, -268.004818, -326.543524, -375.558951, -458.097032, -588.931449, -782.835883, -1027.0032, -1418.94953, -2002.20341, -2837.66181, -3874.62585, -5491.41154, -7959.08117, -11548.697}, 
{259.307093, 252.739349, 247.138517, 242.255957, 237.843023, 233.651073, 229.431464, 224.935552, 219.914693, 212.564812, 204.814872, 197.038402, 189.608933, 184.426978, 179.72829, 175.275605, 170.831659, 165.499518, 159.965458, 154.256084, 148.398001, 142.550276, 136.554067, 130.382995, 124.010678, 117.279218, 110.346361, 103.238334, 95.9813658, 88.6520651, 81.2061239, 73.6496167, 65.9886178, 58.2943521, 50.4816833, 42.5306256, 34.4211932, 26.1073599, 17.6055963, 8.90633292, 0., -8.90871776, -18.1293461, -27.7571564, -37.8874201, -48.2672439, -59.4793298, -71.758215, -85.3384368, -99.2324482, -115.385705, -134.521577, -157.363436, -178.258714, -206.857097, -246.432332, -300.258165, -343.746479, -419.177631, -540.970111, -723.542412, -953.301677, -1328.70045, -1895.63073, -2708.54494, -3712.21403, -5281.5743, -7684.02071, -11186.9482}, 
{230.644686, 225.549343, 221.276528, 217.601877, 214.301026, 211.14961, 207.923267, 204.397631, 200.348339, 194.200203, 187.620012, 180.923731, 174.427326, 169.687863, 165.283766, 161.034559, 156.759767, 151.753696, 146.571178, 141.241823, 135.795243, 130.399758, 124.890791, 119.24247, 113.428926, 107.313599, 101.025583, 94.5832827, 88.0051028, 81.3452456, 74.5719997, 67.6894509, 60.7016853, 53.6714628, 46.5207255, 39.2300897, 31.7801716, 24.1207947, 16.2756849, 8.23777562, 0., -8.25393396, -16.7981776, -25.7161074, -35.0911002, -44.6778543, -55.0198957, -66.3320723, -78.829232, -91.6435998, -106.505696, -124.063416, -144.964658, -163.775287, -189.658045, -225.69364, -274.962784, -313.083784, -381.584712, -494.531239, -665.989037, -881.683134, -1240.70113, -1790.51499, -2581.00657, -3553.91473, -5081.14038, -7427.04393, -10855.9858}, 
{204.949773, 201.071126, 197.905217, 195.246737, 192.890376, 190.630825, 188.262774, 185.580914, 182.379935, 177.278699, 171.718057, 165.963028, 160.278636, 155.914796, 151.757677, 147.678343, 143.547856, 138.840123, 133.982228, 129.004099, 123.935661, 118.939037, 113.859081, 108.672842, 103.357369, 97.8068086, 92.1142731, 86.2899726, 80.3441169, 74.3057015, 68.1586367, 61.9056181, 55.5493414, 49.1434918, 42.6193799, 35.9593056, 29.1455691, 22.1397473, 14.9531524, 7.57637353, 0., -7.6401998, -15.5568876, -23.817546, -32.4896574, -41.3177284, -50.8214078, -61.1973687, -72.6422837, -84.3961618, -97.9950056, -114.018154, -133.044944, -149.885588, -173.196203, -205.863779, -250.775306, -283.720583, -345.522668, -449.907425, -610.600723, -812.791405, -1155.8164, -1687.80212, -2456.25049, -3402.10245, -4893.16741, -7190.53768, -10555.3056}, 
{182.534813, 179.527507, 177.178917, 175.294163, 173.678361, 172.136629, 170.474086, 168.49585, 166.007037, 161.766511, 157.044147, 152.063566, 147.048389, 142.99358, 139.042878, 135.111368, 131.114134, 126.685851, 122.134173, 117.486347, 112.76962, 108.123033, 103.417319, 98.6350069, 93.7586259, 88.7227275, 83.5770077, 78.3231858, 72.9629811, 67.4971152, 61.9287041, 56.2598663, 50.4927201, 44.6712063, 38.738892, 32.6811664, 26.4834191, 20.1387421, 13.621741, 6.91472423, 0., -7.06724582, -14.4067335, -22.0653056, -30.0898049, -38.19385, -46.8907973, -56.3607792, -66.783928, -77.4921481, -89.8530911, -104.386181, -121.610841, -136.610285, -157.514631, -187.017785, -227.813654, -255.806599, -311.195892, -407.391261, -557.802432, -747.270316, -1074.91108, -1588.43805, -2335.48049, -3259.15168, -4720.71301, -6976.88882, -10284.4035}, 
{165.290507, 162.286376, 160.078255, 158.450848, 157.188853, 156.076972, 154.899907, 153.442357, 151.489025, 147.77082, 143.54775, 139.026033, 134.411887, 130.567225, 126.78029, 122.995023, 119.155364, 115.021951, 110.795345, 106.492805, 102.131591, 97.7874608, 93.395776, 88.9503958, 84.4451802, 79.8803286, 75.2408251, 70.5179937, 65.7031583, 60.7565616, 55.7130412, 50.5763533, 45.3502544, 40.0687128, 34.693188, 29.2153516, 23.6268752, 17.9944365, 12.2046983, 6.21932983, 0., -6.54855977, -13.3849676, -20.5247794, -27.983551, -35.3990376, -43.315716, -51.9002622, -61.3193524, -70.9677412, -82.0927948, -95.1699579, -110.674675, -124.050902, -142.81816800000001, -169.464514, -206.47798, -229.664784, -278.867519, -367.246955, -507.963863, -685.99481, -999.053177, -1492.11906, -2218.17943, -3129.51904, -4573.22466, -6797.68259, -10051.2791}, 
{152.137356, 148.560618, 146.028799, 144.288529, 143.086442, 142.169167, 141.283337, 140.175584, 138.592538, 135.133503, 131.151371, 126.851704, 122.440065, 118.735667, 115.084965, 111.448061, 107.78506, 103.948328, 100.048801, 96.0896769, 92.0741556, 87.9923766, 83.8658206, 79.7029095, 75.5120655, 71.371298, 67.1916066, 62.9535786, 58.6378009, 54.159387, 49.5905871, 44.9381778, 40.2089354, 35.4261284, 30.5734448, 25.6510642, 20.6591665, 15.7676707, 10.7391215, 5.50580328, 0., -6.0717388, -12.4611006, -19.1455073, -26.1023811, -32.8645712, -40.0319015, -47.7596232, -56.2029876, -64.7948193, -74.7017665, -86.3680509, -100.237894, -112.156793, -129.007185, -153.07278, -186.637291, -205.292071, -248.690132, -329.792132, -461.558725, -629.404881, -928.92831, -1400.89323, -2107.06889, -3013.74267, -4448.12437, -6647.19318, -9847.92831}, 
{140.781753, 136.682234, 133.820464, 131.915018, 130.684473, 129.847403, 129.122385, 128.227993, 126.882804, 123.587991, 119.766493, 115.623843, 111.365579, 107.785616, 104.265758, 100.776187, 97.2870864, 93.7156472, 90.1062419, 86.4502507, 82.7390537, 78.8895359, 74.9973708, 71.0837365, 67.1698115, 63.3987972, 59.6210392, 55.8089063, 51.9347676, 47.8783341, 43.7416956, 39.5342838, 35.2655307, 30.9463936, 26.5841688, 22.1876782, 17.7657435, 13.5945717, 9.30864537, 4.80983222, 0., -5.61379701, -11.5769522, -17.8296726, -24.3121648, -30.4556994, -36.9129944, -43.8278312, -51.3439913, -58.9189753, -67.6573579, -77.9774331, -90.2974945, -100.815114, -115.857596, -137.531524, -167.94348, -182.552124, -220.771131, -295.366252, -419.10324, -577.762732, -865.065819, -1317.76997, -2006.19433, -2910.75903, -4337.91853, -6512.61869, -9659.80539}, 
{130.494647, 126.16847, 123.140568, 121.124242, 119.832792, 118.979516, 118.277715, 117.440689, 116.181737, 112.994264, 109.299423, 105.29847, 101.192663, 97.752079, 94.3816266, 91.0550345, 87.7460315, 84.3994753, 81.0295143, 77.6214257, 74.1604867, 70.521654, 66.8446537, 63.1588911, 59.4937719, 56.0305044, 52.5859699, 49.1288527, 45.6278369, 41.946012, 38.1998949, 34.400408, 30.5584734, 26.6696049, 22.7662965, 18.8656338, 14.9847023, 11.5026577, 7.92968738, 4.13804898, 0., -5.16541025, -10.7114327, -16.5445265, -22.5711508, -28.1368435, -33.9333535, -40.0915086, -46.7421363, -53.3523053, -60.9821059, -70.0278696, -80.8859276, -90.0081075, -103.313046, -122.774876, -150.36773, -161.586559, -195.476351, -264.602911, -381.532045, -532.136331, -809.061255, -1246.06057, -1920.26372, -2824.73561, -4245.42361, -6395.05448, -9486.35496}, 
{120.218374, 116.363283, 113.606148, 111.699338, 110.395219, 109.44616, 108.604528, 107.622691, 106.253017, 103.165789, 99.6282923, 95.8257276, 91.9432956, 88.7062903, 85.5437823, 82.4249355, 79.3189135, 76.1411207, 72.935984, 69.6941707, 66.4063484, 62.9681364, 59.5032696, 56.0404345, 52.6083179, 49.374891, 46.1738422, 42.9781445, 39.7607707, 36.4023776, 33.0051808, 29.5790797, 26.1339737, 22.6437028, 19.1686496, 15.7331373, 12.361489, 9.52475412, 6.62183918, 3.4983769, 0., -4.71442103, -9.83778635, -15.2497584, -20.8299995, -25.8728603, -31.0774398, -36.5575249, -42.4269028, -48.1422156, -54.7372534, -62.5886611, -72.0730838, -79.7205999, -91.2930477, -108.705699, -133.873825, -142.639015, -173.419696, -238.560614, -350.406513, -494.295262, -763.592234, -1191.33966, -1857.42283, -2763.52395, -4177.2945, -6299.89362, -9332.48046}, 
{110.75343, 107.590506, 105.231651, 103.484465, 102.156549, 101.055506, 99.9889362, 98.7644409, 97.1896215, 94.1790465, 90.790563, 87.1889854, 83.539128, 80.5146295, 77.5679503, 74.660375, 71.7531881, 68.7101427, 65.6290678, 62.5102607, 59.3540188, 56.1009217, 52.8348715, 49.5800526, 46.3606494, 43.309762, 40.2990927, 37.3092595, 34.3208804, 31.2458379, 28.1609798, 25.0744183, 21.9942657, 18.872276, 15.7954629, 12.7944818, 9.89998836, 7.66367431, 5.38674428, 2.89143919, 0., -4.26469321, -8.96129429, -13.9478179, -19.0822789, -23.6261418, -28.2725914, -33.1182624, -38.2597896, -43.1364487, -48.7651769, -55.5055526, -63.7171541, -69.9228172, -79.8535597, -95.4036569, -118.467384, -125.376221, -153.812357, -215.895186, -323.744101, -462.041787, -725.217765, -1146.38373, -1806.37838, -2714.14396, -4118.48155, -6208.22748, -9172.21812}, 
{102.109544, 99.7569725, 97.8623984, 96.2978131, 94.9352071, 93.6465714, 92.3038964, 90.7791731, 88.9443922, 85.9985892, 82.7558921, 79.3574736, 75.9445064, 73.1306553, 70.3956043, 67.6915293, 64.9706062, 62.0343038, 59.0457882, 56.0175184, 52.9619535, 49.8802388, 46.8006721, 43.7402377, 40.7159199, 37.8110855, 34.9497834, 32.1224447, 29.3195008, 26.4924822, 23.6862811, 20.9068892, 18.1602981, 15.3777876, 12.6699457, 10.0726484, 7.62177209, 5.93500366, 4.23368424, 2.32096572, 0., -3.81327259, -8.07360397, -12.6229576, -17.3032971, -21.3601287, -25.4704559, -29.714825, -34.1737824, -38.2677164, -43.0013945, -48.7194263, -55.7664212, -60.5911985, -68.9924741, -82.8731736, -104.136223, -109.710419, -136.462466, -196.284943, -301.070425, -434.866515, -693.100714, -1109.41568, -1764.10962, -2672.48073, -4063.17178, -6111.48937, -8992.74014}, 
{94.2964493, 92.7695158, 91.3437153, 89.9575739, 88.5496178, 87.0583731, 85.422366, 83.5801226, 81.4701691, 78.5889691, 75.4939364, 72.300422, 69.1237773, 66.5079261, 63.9682178, 61.4485745, 58.8929184, 56.0413661, 53.1231674, 50.1637667, 47.1886082, 44.2663165, 41.3618836, 38.4834818, 35.6392833, 32.8548298, 30.1139763, 27.4179471, 24.7679666, 22.1583995, 19.6000738, 17.0969577, 14.6530196, 12.1827008, 9.81530731, 7.59061801, 5.54841194, 4.35432756, 3.17194063, 1.79068637, 0., -3.35720491, -7.16636284, -11.2594299, -15.4683622, -19.038262, -22.6226808, -26.2883164, -30.1018666, -33.4687308, -37.3814243, -42.1711641, -48.1691674, -51.7021831, -58.7076835, -71.1186729, -90.8681555, -95.5538487, -121.178158, -179.408204, -281.911104, -412.260057, -666.403945, -1078.65842, -1727.59583, -2634.41935, -4005.5522, -6001.11258, -8781.21869}, 
{87.7713906, 86.4960091, 85.1811683, 83.7968101, 82.3128763, 80.6993088, 78.9260494, 76.9630399, 74.7802222, 72.0511526, 69.1607128, 66.1973987, 63.2497064, 60.7494344, 58.3044555, 55.8659448, 53.3850774, 50.6152646, 47.784551, 44.9232173, 42.0615442, 39.301138, 36.5724235, 33.8771513, 31.2170719, 28.5794753, 25.9863566, 23.4452507, 20.9636922, 18.5573668, 16.2223978, 13.9630598, 11.783627, 9.60000198, 7.54017972, 5.64378337, 3.95043606, 3.11669301, 2.31847245, 1.34862468, 0., -2.86170759, -6.14562894, -9.68805125, -13.3252617, -16.3257928, -19.3207882, -22.3736369, -25.547728, -28.2713039, -31.4969588, -35.5421401, -40.7242953, -43.2310953, -49.1616752, -60.4853932, -79.1716077, -83.2430703, -108.194389, -165.573564, -266.928597, -396.429249, -647.75824, -1054.36803, -1693.5088, -2596.01584, -3934.58632, -5846.45193, -8468.84437}, 
{82.0309902, 80.8879567, 79.6096357, 78.1899457, 76.622805, 74.9021323, 73.0218458, 70.9758641, 68.7581056, 66.1932829, 63.5122026, 60.7764656, 58.0476726, 55.637402, 53.257286, 50.8689344, 48.4339568, 45.7340129, 42.9826428, 40.2134361, 37.4599828, 34.8511951, 32.2872114, 29.7634926, 27.2754996, 24.774658, 22.3180779, 19.918834, 17.5900013, 15.3641527, 13.2270657, 11.1840158, 9.24027873, 7.31904842, 5.54051468, 3.94278571, 2.56396973, 2.03191979, 1.55910133, 0.947724641, 0., -2.3603199, -5.10118416, -8.06899945, -11.1101725, -13.5280581, -15.9293355, -18.3776321, -20.9365754, -23.0678086, -25.6777369, -29.0707811, -33.5513625, -35.1556526, -40.1636214, -50.5869898, -68.4374784, -72.203814, -96.8299089, -153.736681, -254.34505, -383.750348, -632.350248, -1032.43624, -1659.18033, -2552.95931, -3850.45028, -5664.25527, -8106.97631}, 
{76.9707195, 75.8621394, 74.5611871, 73.0790334, 71.4268495, 69.6158063, 67.6570748, 65.5618261, 63.3412311, 60.9450827, 58.4704813, 55.9531493, 53.4288089, 51.087234, 48.7484747, 46.3866325, 43.9758091, 41.3365635, 38.6579569, 35.975508, 33.3247354, 30.8558433, 28.4437906, 26.0782215, 23.7487807, 21.3748538, 19.0444472, 16.775309, 14.585187, 12.5188412, 10.5562024, 8.70421353, 6.96981753, 5.28812718, 3.76664768, 2.44105406, 1.34702134, 1.06977331, 0.875616728, 0.580407118, 0., -1.85821402, -4.04652853, -6.42570206, -8.85649313, -10.6846322, -12.4919171, -14.3451174, -16.311003, -17.8933143, -19.9470619, -22.7642269, -26.6367909, -27.4525786, -31.6693904, -41.3408701, -58.5206616, -62.2399847, -86.7838771, -143.414952, -243.395824, -373.020461, -618.457409, -1010.56686, -1621.66818, -2501.0912, -3748.66476, -5451.59296, -7697.07988}, 
{72.4860498, 71.3353381, 69.967892, 68.4061262, 66.6724552, 64.7892935, 62.7790558, 60.6641567, 58.4670106, 56.2362749, 53.9576244, 51.6429764, 49.3042485, 47.0143357, 44.6997869, 42.3481283, 39.9468865, 37.3618686, 34.7510078, 32.1505183, 29.5966139, 27.2544382, 24.9797041, 22.7570538, 20.5711298, 18.3145381, 16.1007718, 13.9512876, 11.8875423, 9.96151604, 8.15193297, 6.46804054, 4.91908621, 3.45552525, 2.16891421, 1.09201745, 0.257599308, 0.200019021, 0.249808109, 0.239092969, 0., -1.36056209, -2.99516214, -4.78158666, -6.59762218, -7.83508972, -9.05212751, -10.3209083, -11.7136048, -12.7828905, -14.3282371, -16.629617, -19.9670027, -20.0985972, -23.6348501, -32.6644414, -49.2760512, -53.1554874, -77.7554513, -134.125772, -233.316278, -363.036691, -604.357163, -986.463671, -1578.03011, -2436.25295, -3624.75043, -5205.53533, -7240.62045}, 
{68.4724527, 67.2243335, 65.7618201, 64.1132768, 62.3070676, 60.3715567, 58.3351084, 56.2260867, 54.0728559, 51.9965824, 49.8957072, 47.7614734, 45.5851245, 43.3341123, 41.0329879, 38.682511, 36.2834413, 33.7488806, 31.2023098, 28.6795519, 26.21643, 23.9863352, 21.8324951, 19.7357053, 17.6767614, 15.5281868, 13.4223587, 11.3833816, 9.43536035, 7.63226104, 5.95638238, 4.41988452, 3.0349276, 1.76952966, 0.697649771, -0.150895095, -0.746287965, -0.607577646, -0.336535071, -0.0837969507, 0., -0.872536286, -1.96058503, -3.16008079, -4.36695815, -5.01900493, -5.65356114, -6.34982007, -7.18697498, -7.77160669, -8.84456592, -10.6740909, -13.5284199, -13.0704323, -16.0158687, -24.475111, -40.558541, -44.7542269, -69.4437895, -125.386536, -223.341773, -352.596143, -588.326948, -957.830483, -1525.32388, -2354.28602, -3474.22797, -4923.15275, -6739.06338}, 
{68.2547617, 64.4030038, 61.0624053, 58.1628557, 55.6342446, 53.4064616, 51.4093962, 49.572938, 47.8269765, 45.957039, 44.0951222, 42.2288605, 40.3458883, 38.4667234, 36.5329635, 34.5190896, 32.3995828, 29.9385344, 27.4049711, 24.8575296, 22.354847, 20.1021213, 17.9528037, 15.9069066, 13.9644425, 12.1438868, 10.4194037, 8.78362046, 7.22916423, 5.71567078, 4.28195543, 2.93384198, 1.67715426, 0.346381608, -0.812783878, -1.72598459, -2.3188629, -1.86839962, -1.20836336, -0.523861146, 0., -0.627400983, -1.46345147, -2.37105288, -3.21310661, -3.25085515, -3.18952241, -3.13267337, -3.18387301, -2.80547233, -2.9987359, -4.1237143, -6.54045811, -6.83371701, -10.6489628, -19.8563666, -36.3260991, -39.4859101, -62.6253602, -116.591589, -212.231736, -342.977808, -571.922343, -920.991834, -1451.388, -2230.85247, -3261.692, -4562.57208, -6152.15819}, 
{68.4539243, 61.4483189, 55.8628302, 51.4992593, 48.1594069, 45.6450741, 43.7580616, 42.3001706, 41.0732018, 39.3851077, 37.7290771, 36.1044501, 34.5105672, 33.1173068, 31.6852556, 30.1455386, 28.4292807, 26.0790926, 23.5700193, 20.9885914, 18.4213396, 16.0987041, 13.9057423, 11.871421, 10.0247072, 8.59280239, 7.32714494, 6.17740787, 5.09326418, 3.90357687, 2.727153, 1.56198962, 0.406083746, -1.08429602, -2.43273282, -3.50453824, -4.16502385, -3.33725149, -2.20568241, -1.01252808, 0., -0.497549335, -1.22525202, -1.99042368, -2.60037993, -2.04190726, -1.27106206, -0.423371592, 0.365636893, 1.81885925, 2.59897586, 2.22709019, 0.224305732, -1.26648705, -6.39468643, -16.6879037, -33.6737502, -35.6379855, -56.6468132, -107.524585, -199.095654, -329.609444, -547.615087, -868.664624, -1349.75899, -2057.6127, -2978.63703, -4119.31905, -5486.14585}, 
{65.4719479, 57.4208681, 51.1223649, 46.3288876, 42.792886, 40.2668095, 38.5031077, 37.2542302, 36.2726266, 34.6627386, 33.0842266, 31.5487434, 30.0679416, 28.9024083, 27.7152878, 26.4186592, 24.9246012, 22.6910859, 20.2659417, 17.7428902, 15.2156527, 12.9020086, 10.7219986, 8.71972156, 6.93927587, 5.70709501, 4.6720087, 3.76518156, 2.9177782, 1.89110459, 0.854127432, -0.194045212, -1.25430528, -2.72905117, -4.05706578, -5.07863845, -5.63405853, -4.49776766, -3.002242, -1.41410996, 0., -0.232086191, -0.689233061, -1.15585079, -1.41634954, -0.326326355, 1.02947019, 2.49510466, 3.91464163, 6.02483426, 7.41998307, 7.58707717, 6.01310568, 4.25526237, -1.09775015, -11.3870246, -27.9536538, -30.037432, -49.9212697, -97.7867791, -183.815572, -309.038656, -511.089458, -801.715921, -1230.32357, -1857.50504, -2661.44699, -3640.93761, -4794.76507}, 
{60.4371914, 52.8944439, 47.0061111, 42.5344545, 39.2417357, 36.8902165, 35.2421584, 34.0598231, 33.1054723, 31.4888045, 29.8856696, 28.3193545, 26.8131459, 25.6686144, 24.5194497, 23.2776252, 21.8551143, 19.7300162, 17.4217281, 15.0157733, 12.5976748, 10.34573, 8.215578, 6.25563199, 4.51430532, 3.33212144, 2.34853948, 1.49512871, 0.703458391, -0.282898639, -1.26917813, -2.24861224, -3.21443316, -4.5363881, -5.68058815, -6.48965946, -6.80622818, -5.41661253, -3.64226975, -1.74834914, 0., 0.144274795, 0.0900209157, 0.0494304766, 0.23469559, 1.80159795, 3.64130425, 5.58857079, 7.47815384, 9.93910759, 11.6941713, 12.260382, 11.1547768, 9.8495151, 5.12446293, -4.2853914, -19.6450596, -22.6973878, -42.0384195, -86.7420329, -165.882106, -281.689148, -463.767142, -721.985092, -1095.87966, -1636.34431, -2319.17767, -3138.9521, -4090.23996}, 
{54.4780137, 48.4428385, 43.6791707, 39.9986734, 37.21301, 35.1338437, 33.5728378, 32.3416556, 31.2519604, 29.5621783, 27.8585046, 26.1738973, 24.5413145, 23.2625116, 21.9941304, 20.6616102, 19.1903901, 17.1513851, 14.9663681, 12.7025882, 10.427294, 8.26356395, 6.20048559, 4.282976, 2.55595225, 1.31323844, 0.251281779, -0.684563488, -1.54894311, -2.57958539, -3.57482051, -4.51606121, -5.38472022, -6.45481104, -7.29810532, -7.77897552, -7.76179404, -6.16045051, -4.16999329, -2.03498794, 0., 0.606819966, 1.05792541, 1.54218478, 2.24846652, 4.24757601, 6.49367, 8.82284217, 11.0711862, 13.6883341, 15.6514259, 16.5511399, 15.9781543, 15.634255, 12.1545697, 4.28533388, -9.22721724, -13.6309911, -32.5879525, -73.754209, -144.785868, -247.984626, -407.069825, -631.311504, -949.22515, -1399.94534, -1960.88484, -2624.88686, -3384.79461}, 
{48.7227738, 44.6398442, 41.3066454, 38.604258, 36.4137625, 34.6162396, 33.0927698, 31.7244338, 30.3923122, 28.5817329, 26.7278303, 24.8699861, 23.0475818, 21.5306864, 20.0357196, 18.509788, 16.8999986, 14.9106943, 12.8288514, 10.6986821, 8.56439891, 6.48920575, 4.49072666, 2.60557732, 0.870373416, -0.504197236, -1.72521994, -2.82570787, -3.83867419, -4.96010825, -5.99485657, -6.91074186, -7.67558681, -8.43282424, -8.90442271, -8.98796086, -8.58101735, -6.79594598, -4.62964023, -2.2937687, 0., 1.13083567, 2.15989592, 3.23917679, 4.5206743, 6.91731822, 9.51579733, 12.1637342, 14.7087513, 17.3991688, 19.5216327, 20.7634863, 20.8120732, 21.7274658, 19.8751872, 13.9934892, 2.82062372, -2.85138016, -21.1595584, -58.1871699, -120.017474, -208.348795, -342.419193, -531.534524, -793.157935, -1154.12298, -1595.62429, -2110.26622, -2690.65314}, 
{44.2998304, 42.0592533, 40.0536373, 38.233922, 36.5510472, 34.9559528, 33.3995785, 31.832864, 30.2067492, 28.2463412, 26.2187454, 24.1652348, 22.1270823, 20.3197254, 18.5406064, 16.7613322, 14.9535099, 12.9634457, 10.9381676, 8.89940251, 6.86887778, 4.85635082, 2.90030642, 1.02725964, -0.73627446, -2.27482884, -3.6864212, -4.98011725, -6.16498272, -7.38561978, -8.46134313, -9.34700391, -9.99745327, -10.418932, -10.4943457, -10.1579897, -9.34415937, -7.38976337, -5.06543819, -2.54443376, 0., 1.69160825, 3.34134796, 5.05717119, 6.94702998, 9.71653493, 12.6369161, 15.5770622, 18.4058619, 21.1982666, 23.5346769, 25.2015565, 25.9853686, 28.2471315, 28.1689322, 24.5074123, 16.0192136, 9.6283068, -7.34292715, -39.4047781, -91.0675361, -163.20536, -271.236932, -424.493521, -630.475917, -904.692043, -1232.45178, -1606.61452, -2020.03966}, 
{31.3611692, 33.4138344, 34.6123875, 35.0543707, 34.837326, 34.0587958, 32.8163219, 31.2074467, 29.3297123, 27.404283, 25.3556298, 23.2318459, 21.0810247, 19.0533799, 17.0540359, 15.0902378, 13.1692306, 11.3600238, 9.58339213, 7.82187479, 6.05801097, 4.17746025, 2.29839333, 0.442101287, -1.37012481, -3.08825048, -4.73122545, -6.289256, -7.75254846, -9.29802506, -10.6544898, -11.7374626, -12.4624633, -12.5589911, -12.202995, -11.3844033, -10.0931443, -7.98083331, -5.51103668, -2.80900778, 0., 2.33855583, 4.71445469, 7.18331415, 9.80075179, 13.1193625, 16.4989955, 19.7964776, 22.8686353, 25.1976694, 27.1648828, 28.7769525, 30.0405554, 34.3745852, 37.0086155, 36.5844368, 31.7438393, 24.3660043, 8.56037482, -18.326215, -58.9469312, -111.659095, -192.003407, -312.956989, -466.821193, -649.230038, -863.162643, -1109.76185, -1390.17051}, 
{37.107078, 37.9494898, 38.1475584, 37.7748672, 36.9049999, 35.61154, 33.9680713, 32.0481773, 29.9254415, 27.7748251, 25.5279832, 23.2179484, 20.8777534, 18.5825868, 16.3064628, 14.0655515, 11.8760232, 9.83055734, 7.83821112, 5.88455095, 3.9551433, 1.96138975, -0.00731245436, -1.93573091, -3.80863323, -5.6288923, -7.35592834, -8.96726683, -10.4404333, -11.8918334, -13.1045604, -14.0005876, -14.5018883, -14.3327386, -13.6918882, -12.5803894, -10.9992944, -8.70169555, -6.03578926, -3.10181195, 0., 2.87088123, 5.82853593, 8.89213931, 12.0808666, 15.681342, 19.3383122, 22.9639726, 26.4705188, 29.5587418, 32.4368037, 35.1014621, 37.5494744, 42.4526156, 46.0626187, 47.3062344, 45.1102132, 39.9921181, 28.651562, 9.37897069, -19.5352305, -57.2665312, -113.12676, -196.55736, -298.235907, -413.117734, -537.811019, -668.510728, -801.411827}, 
{54.2496337, 50.3266662, 46.8385324, 43.7179324, 40.8975662, 38.3101336, 35.8883348, 33.5648695, 31.2724379, 28.7494751, 26.2006517, 23.6363737, 21.0670471, 18.525831, 15.9912767, 13.4646886, 10.9473713, 8.38916945, 5.86343143, 3.39204573, 0.99690084, -1.26416265, -3.41958917, -5.46187108, -7.3835007, -9.22601578, -10.9132451, -12.4180629, -13.7133433, -14.8183369, -15.6409909, -16.1356292, -16.2565754, -15.8913591, -15.0878157, -13.8269866, -12.089913, -9.575889, -6.66040207, -3.43719234, 0., 3.32808289, 6.77440878, 10.3369783, 14.0137919, 17.7915291, 21.6840401, 25.6938539, 29.8234995, 34.4683192, 39.0809036, 43.5066565, 47.5909816, 52.0892308, 55.5728806, 57.5233557, 57.4220807, 57.0603078, 52.6857029, 42.8557596, 26.1279713, 0.306807571, -33.7911656, -77.5102898, -128.97305, -194.828091, -250.841278, -280.499173, -267.288336}, 
{45.9632263, 43.5948047, 41.2944117, 39.0497927, 36.8486933, 34.6788586, 32.5280343, 30.3839656, 28.2343979, 26.0514214, 23.8446989, 21.6082379, 19.3360459, 16.9792869, 14.5919492, 12.1851779, 9.7701179, 7.3618403, 4.9659933, 2.59215129, 0.249888704, -2.11170034, -4.39836902, -6.57635081, -8.61187918, -10.3777734, -11.9710468, -13.3952985, -14.6541278, -15.978402, -17.0535448, -17.792248, -18.1072036, -17.711688, -16.7975747, -15.3573219, -13.3833878, -10.6397615, -7.43875759, -3.86422188, 0., 3.91966746, 8.02148764, 12.2817728, 16.6768352, 21.1786397, 25.7695849, 30.427722, 35.1311022, 39.788373, 44.4747508, 49.1960481, 53.9580778, 59.5390344, 64.8633958, 69.6280218, 73.5297724, 77.145323, 78.9397916, 78.2581119, 74.4452177, 65.9571768, 54.8055203, 38.1237308, 26.1576508, 25.2296282, 44.2436292, 91.4582154, 175.131949}, 
{51.3870855, 48.6794894, 46.1831565, 43.8403542, 41.59335, 39.3844114, 37.155806, 34.8498012, 32.4086645, 29.4555478, 26.3794804, 23.2503763, 20.1381492, 17.3421612, 14.6110986, 11.9230958, 9.25628726, 6.47079096, 3.70996456, 0.999149184, -1.63631406, -4.15046897, -6.54683554, -8.80831871, -10.9178234, -12.8574803, -14.6112782, -16.1624317, -17.4941555, -18.7345028, -19.6639141, -20.2076685, -20.2910453, -19.6603376, -18.491405, -16.7811208, -14.5263586, -11.5733769, -8.12991022, -4.25307808, 0., 4.55743405, 9.3827832, 14.4248361, 19.6323815, 24.7770149, 30.0555953, 35.4877886, 41.0932604, 46.9053614, 52.9245989, 59.1651649, 65.6412515, 72.5154609, 79.5942109, 86.8323292, 94.184644, 100.942288, 107.989262, 115.545873, 123.832425, 131.321797, 143.476579, 160.747982, 191.051251, 237.852606, 314.950336, 433.559716, 604.896017}, 
{52.5791609, 52.189412, 50.9809292, 49.0755733, 46.5952055, 43.6616866, 40.3968776, 36.9226396, 33.3608334, 30.1281937, 26.9337584, 23.7814391, 20.6751474, 17.664573, 14.689538, 11.7356427, 8.78848729, 5.69032192, 2.6274371, -0.357226835, -3.22072955, -5.77288536, -8.1768974, -10.4487235, -12.6043213, -14.8262476, -16.8972218, -18.766562, -20.3835865, -21.7939441, -22.8120902, -23.348811, -23.3148923, -22.3006067, -20.6654593, -18.4484415, -15.6885449, -12.4879459, -8.79717705, -4.62995592, 0., 5.28008011, 10.9150172, 16.810651, 22.8728212, 28.6611796, 34.566229, 40.6322844, 46.9036604, 53.1753873, 59.8407787, 67.0438633, 74.92867, 83.7004277, 93.417485, 104.199391, 116.165694, 127.53876, 141.094195, 157.710418, 178.265852, 202.213911, 234.708038, 276.259409, 335.688216, 421.713219, 535.1499, 678.789558, 855.423494}
};

  const double a4tab[50][69] = {
{-4795.25483, -4657.11547, -4533.57386, -4420.99082, -4315.72713, -4214.14359, -4112.601, -4007.46016, -3895.08186, -3749.86323, -3598.9142, -3447.38103, -3300.40998, -3183.58631, -3073.44169, -2966.94678, -2861.07226, -2743.59299, -2624.35375, -2504.00351, -2383.19127, -2264.35938, -2145.64608, -2026.98299, -1908.30176, -1789.65691, -1670.80798, -1551.63744, -1432.02772, -1310.85768, -1189.41483, -1067.98309, -946.846341, -827.165104, -707.996037, -589.272406, -470.927472, -353.153181, -235.520643, -117.859652, 0., 117.602012, 235.993705, 355.595891, 476.829381, 595.994952, 719.281464, 848.757742, 986.492612, 1125.29678, 1280.20043, 1456.97565, 1661.39451, 1865.75367, 2122.69078, 2451.36808, 2870.9478, 3306.89371, 3909.54591, 4735.54605, 5841.53576, 7186.8075, 9120.05042, 11591.8611, 15444.2924, 19758.4141, 29767.7712, 49167.7902, 81653.8973}, 
{-4386.2486, -4276.4998, -4166.87176, -4057.27842, -3947.63375, -3837.85169, -3727.8462, -3617.53124, -3506.82077, -3396.24985, -3284.86287, -3172.32536, -3058.30283, -2934.12142, -2811.12175, -2692.30507, -2580.67265, -2501.29049, -2426.26919, -2349.78411, -2266.0106, -2143.26267, -2011.92156, -1876.50719, -1741.53944, -1625.69861, -1513.68005, -1404.33951, -1296.53273, -1186.17111, -1076.23248, -966.750311, -857.758075, -749.581249, -641.844515, -534.46455, -427.358033, -320.645434, -213.958122, -107.131257, 0., 107.056887, 214.965285, 324.107473, 434.865731, 543.216635, 655.710448, 774.491733, 901.705051, 1031.01986, 1176.44587, 1343.51767, 1537.76987, 1733.97589, 1980.73599, 2295.88922, 2697.27467, 3113.77504, 3687.76832, 4472.67611, 5521.92004, 6796.91537, 8627.10278, 10961.4315, 14608.3709, 18698.1916, 28214.4395, 46677.342, 77607.1262}, 
{-4089.70007, -3999.83272, -3900.40845, -3793.52861, -3681.29454, -3565.8076, -3449.16915, -3333.48054, -3220.84311, -3124.77905, -3031.40055, -2938.24064, -2842.83236, -2728.08787, -2612.00939, -2497.97828, -2389.37592, -2309.8317, -2234.37973, -2158.30216, -2076.88115, -1964.99854, -1846.4969, -1724.8185, -1603.40559, -1496.11382, -1391.80673, -1289.76124, -1189.25426, -1088.05749, -987.55515, -887.626246, -788.149784, -688.774549, -589.701849, -490.902776, -392.348421, -294.336138, -196.38025, -98.321342, 0., 97.9874303, 196.858427, 297.074709, 399.097994, 499.259521, 603.803679, 714.844378, 834.49553, 956.8933, 1095.32044, 1255.08195, 1441.48285, 1629.7946, 1867.36916, 2171.52496, 2559.58043, 2963.77126, 3520.53168, 4281.21321, 5297.16736, 6532.52802, 8300.29955, 10550.9045, 14057.1258, 18028.6165, 27105.1878, 44569.1112, 73702.6582}, 
{-3869.29607, -3792.59815, -3700.71666, -3596.83571, -3484.13943, -3365.81194, -3245.03736, -3124.99981, -3008.88343, -2913.32449, -2822.67409, -2734.73548, -2647.31191, -2552.51261, -2456.11249, -2358.19241, -2258.83326, -2155.5432, -2052.00489, -1949.32827, -1848.62329, -1756.61105, -1666.5459, -1577.29332, -1487.71881, -1392.40141, -1296.20764, -1199.71758, -1103.51131, -1010.11081, -917.377477, -825.114634, -733.125587, -640.567837, -548.148831, -455.930204, -363.973591, -272.876208, -181.949879, -91.0420086, 0., 90.1866618, 181.269486, 273.857901, 368.561335, 462.482307, 561.139919, 666.546362, 780.71383, 897.948208, 1031.05052, 1185.11547, 1365.23778, 1545.93007, 1775.10198, 2070.08108, 2448.19491, 2844.61617, 3391.68922, 4139.60356, 5138.54867, 6355.55344, 8090.28923, 10297.6099, 13712.8828, 17642.831, 26328.2418, 42776.3866, 69994.5366}, 
{-3688.72344, -3620.28003, -3534.3291, -3434.29407, -3323.59833, -3205.66528, -3083.91832, -2961.78084, -2842.67626, -2739.75984, -2642.83036, -2551.41848, -2465.05485, -2394.42261, -2323.43893, -2247.17347, -2160.69588, -2024.75159, -1882.46418, -1742.633, -1614.05739, -1545.14423, -1489.24234, -1439.30806, -1388.29774, -1306.06025, -1217.90239, -1126.02348, -1032.62284, -945.925062, -859.694111, -773.719216, -687.789602, -600.783951, -513.766253, -426.889953, -340.308496, -254.916558, -169.829861, -84.9053578, 0., 83.447602, 167.794816, 253.817352, 342.290922, 431.243692, 525.297932, 626.328371, 736.209734, 849.215705, 977.862472, 1127.06518, 1301.73898, 1475.10257, 1696.44611, 1983.36334, 2353.44796, 2744.04363, 3285.09418, 4026.29337, 5017.33496, 6227.8996, 7947.72035, 10138.8771, 13497.9681, 17433.977, 25771.8272, 41232.4571, 66536.8049}, 
{-3511.669, -3448.36228, -3367.77851, -3272.99804, -3167.10117, -3053.16822, -2934.27953, -2813.5154, -2693.95617, -2581.95877, -2476.01625, -2377.89827, -2289.37451, -2240.84482, -2193.9966, -2139.14746, -2066.61499, -1903.78345, -1729.07712, -1557.98694, -1406.00384, -1357.64213, -1331.76003, -1316.23912, -1298.96098, -1228.58923, -1147.91061, -1060.49389, -969.90785, -889.094259, -808.499694, -727.943732, -647.245948, -565.245724, -483.134905, -401.125143, -319.42809, -239.108105, -159.183049, -79.523491, 0., 77.5632713, 156.030772, 236.313367, 319.321921, 403.902372, 493.856483, 590.921089, 696.833024, 805.726906, 929.982676, 1074.37806, 1243.69077, 1410.03238, 1623.91323, 1903.17752, 2265.66943, 2649.78746, 3184.59978, 3919.72886, 4904.79718, 6111.47446, 7823.24142, 10012.0358, 13334.7074, 17295.1967, 25324.1697, 39870.6114, 63383.5064}, 
{-3301.81958, -3242.32883, -3167.59761, -3080.04195, -2982.07786, -2876.12137, -2764.5885, -2649.89526, -2534.45768, -2417.79493, -2306.37861, -2203.78347, -2113.58426, -2078.80618, -2047.79337, -2008.34038, -1948.2418, -1778.96537, -1595.16323, -1415.16067, -1257.28299, -1221.1488, -1211.27277, -1213.4629, -1213.52717, -1151.48722, -1077.25191, -994.943767, -908.685307, -833.212394, -757.788871, -682.291925, -606.598741, -529.775993, -452.835581, -375.978895, -299.407326, -224.101762, -149.172293, -74.5085093, 0., 72.3266902, 145.573709, 220.706247, 298.689499, 378.817046, 464.394334, 557.055201, 658.433482, 762.51293, 881.637499, 1020.50106, 1183.79748, 1343.43977, 1550.015, 1821.32939, 2175.18916, 2549.58149, 3074.05924, 3798.35623, 4772.2063, 5968.186, 7667.50096, 9854.41554, 13145.4268, 17119.632, 24873.4951, 38624.1385, 60588.6845}, 
{-3022.862, -2967.66361, -2900.31911, -2822.52015, -2735.95834, -2642.32532, -2543.31273, -2440.61219, -2335.91534, -2225.14199, -2118.06432, -2018.68268, -1930.99743, -1895.33367, -1864.83711, -1828.97824, -1777.2275, -1636.62394, -1484.04202, -1333.92478, -1200.71524, -1162.70827, -1144.95436, -1136.35581, -1125.81492, -1066.2531, -996.945919, -921.188061, -842.274195, -771.873467, -701.556286, -631.267537, -560.952101, -490.197592, -419.449071, -348.794329, -278.321158, -208.548445, -138.960447, -69.472514, 0., 67.5308789, 136.01998, 206.356297, 279.428824, 354.34641, 434.49025, 521.461393, 616.860889, 714.604893, 827.05331, 958.88115, 1114.76342, 1268.04501, 1467.2631, 1729.62474, 2072.337, 2431.15956, 2937.32579, 3640.6217, 4590.83328, 5759.94219, 7431.14749, 9603.34594, 12852.4522, 16800.4251, 24308.029, 37426.3271, 58206.3824}, 
{-2638.48309, -2589.85055, -2532.47575, -2467.52697, -2396.17251, -2319.58066, -2238.91972, -2155.35797, -2070.06371, -1981.87362, -1895.22025, -1812.20453, -1734.92739, -1677.45422, -1625.13573, -1575.28704, -1525.22331, -1463.08574, -1399.03301, -1334.04983, -1269.12095, -1209.3646, -1149.97859, -1090.29423, -1029.64285, -964.385738, -898.012251, -831.041724, -763.993491, -698.671476, -633.796583, -569.374308, -505.410146, -442.333357, -379.556167, -316.914565, -254.24454, -191.09907, -127.710361, -64.0276063, 0., 62.9688578, 126.965942, 192.62382, 260.575062, 328.849163, 401.722993, 480.870349, 567.965028, 657.033913, 760.456481, 882.965295, 1029.29292, 1176.56839, 1368.1692, 1619.86934, 1947.44278, 2282.2555, 2758.25267, 3424.97148, 4331.9491, 5448.65099, 7064.82952, 9196.15657, 12378.1095, 16230.7179, 23515.9972, 36210.4661, 56290.6436}, 
{-2405.71268, -2357.37371, -2305.09333, -2248.92728, -2188.93131, -2125.16116, -2057.67258, -1986.52131, -1911.7631, -1825.78366, -1739.37678, -1655.66624, -1577.77578, -1521.53728, -1472.28315, -1428.05392, -1386.8901, -1345.84384, -1304.33938, -1260.8126, -1213.69937, -1153.95915, -1090.49479, -1024.73273, -958.099404, -895.75918, -833.905408, -772.469359, -711.382305, -650.167044, -589.326712, -528.955972, -469.149485, -410.875751, -353.00606, -295.28554, -237.459321, -178.716175, -119.580127, -60.018845, 0., 59.1806414, 119.414746, 181.265883, 245.297621, 309.449323, 377.958449, 452.438249, 534.501976, 618.058134, 715.506621, 831.542589, 970.861189, 1111.62937, 1295.68177, 1538.3248, 1854.86492, 2173.99376, 2632.27845, 3279.67134, 4166.12475, 5260.24742, 6856.02255, 8984.53783, 12142.789, 15991.1783, 23070.5741, 35201.7283, 54205.3928}, 
{-2204.0302, -2154.47867, -2105.90828, -2057.23597, -2007.37868, -1955.25335, -1899.77692, -1839.86634, -1774.43854, -1688.83889, -1600.98454, -1515.22105, -1435.89401, -1383.22182, -1339.32807, -1302.20921, -1269.86166, -1242.18385, -1214.50944, -1184.07407, -1148.11338, -1090.76646, -1027.60414, -961.10068, -893.730354, -835.115058, -777.722379, -721.167533, -665.065735, -607.620119, -550.422817, -493.653877, -437.493349, -383.360258, -329.700085, -276.19729, -222.536329, -167.664731, -112.298656, -56.4173358, 0., 55.7448171, 112.55891, 170.95477, 231.444888, 291.989405, 356.674102, 427.03241, 504.59776, 583.170483, 675.11035, 785.044031, 917.598198, 1051.62411, 1227.83402, 1461.16475, 1766.55314, 2070.16998, 2511.22455, 3140.16011, 4007.4199, 5080.84751, 6658.68514, 8786.8724, 11925.8417, 15777.5862, 22663.3632, 34249.6139, 52202.7796}, 
{-2019.11349, -1967.86498, -1922.03119, -1879.54631, -1838.34453, -1796.36004, -1751.52704, -1701.77973, -1645.05228, -1560.9017, -1472.99026, -1386.60303, -1307.02508, -1257.73185, -1218.54189, -1187.46414, -1162.50754, -1146.0165, -1129.93027, -1110.52355, -1084.07109, -1029.27314, -967.008674, -900.58221, -833.298258, -778.328451, -725.163331, -673.160561, -621.677803, -567.845774, -514.139863, -460.808512, -408.100162, -357.77398, -307.963395, -258.312561, -208.46563, -157.218748, -105.403279, -53.0025778, 0., 52.510405, 106.100252, 161.230464, 218.361962, 275.485215, 336.519778, 402.914754, 476.119245, 549.893005, 636.450223, 740.315742, 866.014403, 993.013679, 1160.91873, 1384.27734, 1677.6373, 1964.76951, 2387.3094, 2996.11552, 3842.0464, 4892.15254, 6448.71667, 8571.55319, 11686.9318, 15534.6254, 22233.9388, 33309.2935, 50285.1113}, 
{-1852.906, -1799.92433, -1755.7931, -1717.7438, -1683.00798, -1648.81714, -1612.4028, -1570.99649, -1521.82971, -1440.70191, -1354.84952, -1270.0769, -1192.18842, -1146.35433, -1111.26673, -1084.98361, -1065.56296, -1056.86211, -1048.81999, -1037.17485, -1017.66497, -965.979445, -905.925352, -841.260613, -775.743151, -724.511292, -675.390392, -627.586209, -580.304504, -530.041543, -479.796376, -429.858557, -380.517643, -333.704815, -287.411353, -241.270163, -194.914152, -147.113112, -98.7083065, -49.6778856, 0., 49.3678455, 99.8185964, 151.765844, 205.623179, 259.438588, 316.93751, 379.479777, 448.425223, 517.590152, 598.895338, 696.718027, 815.435464, 935.025994, 1094.02532, 1306.57025, 1586.79759, 1856.40165, 2258.93873, 2845.52263, 3667.26715, 4690.50381, 6220.69328, 8329.94908, 11413.3409, 15245.1024, 21759.8548, 32352.1221, 48416.4288}, 
{-1707.35118, -1653.0484, -1609.52501, -1573.71394, -1542.54814, -1512.96053, -1481.88405, -1446.25164, -1402.99624, -1326.96932, -1246.01786, -1165.90739, -1092.40343, -1050.37625, -1018.84473, -995.932496, -979.763174, -974.240986, -969.396721, -961.041762, -944.987495, -897.385648, -841.571127, -781.219183, -720.005065, -672.775514, -627.56569, -583.582244, -540.031828, -493.404962, -446.710883, -400.242696, -354.293506, -310.740658, -267.659322, -224.708907, -181.548824, -137.082706, -92.0280493, -46.3465739, 0., 46.2075789, 93.4937656, 142.233789, 192.802879, 243.351361, 297.36933, 356.121975, 420.874489, 485.626378, 561.814791, 653.61119, 765.18704, 876.888988, 1026.24324, 1226.95119, 1492.71423, 1743.67573, 2124.5183, 2686.36654, 3480.34503, 4472.24266, 5969.19114, 8053.42895, 11092.3504, 14891.8235, 21218.6652, 31349.4551, 46560.7729}, 
{-1584.39249, -1529.62886, -1485.55795, -1449.34221, -1418.14411, -1389.1261, -1359.45065, -1326.28022, -1286.77728, -1218.43376, -1145.95085, -1074.35921, -1008.68952, -971.084549, -942.618, -921.475684, -905.843417, -897.67345, -889.878586, -879.138066, -862.131131, -819.992009, -771.162956, -718.541216, -665.024034, -622.233049, -580.851352, -540.286431, -499.94577, -457.133564, -414.20191, -371.399611, -328.975471, -288.469407, -248.322665, -208.267601, -168.036574, -126.862414, -85.1768186, -42.9119574, 0., 42.9200457, 86.9055837, 132.30718, 179.4754, 226.72537, 277.257271, 332.235845, 392.825836, 453.366138, 524.577678, 610.355535, 714.594789, 817.830593, 956.661925, 1144.32784, 1394.06738, 1625.20104, 1982.45386, 2516.63232, 3278.54291, 4233.71041, 5688.78638, 7733.36167, 10711.2418, 14457.5951, 20587.924, 30272.6475, 44682.1844}, 
{-1485.97338, -1432.0574, -1386.22294, -1346.5141, -1310.975, -1277.64975, -1244.58246, -1209.81724, -1171.39821, -1113.82504, -1054.10405, -995.697115, -942.066106, -909.766207, -883.928673, -862.778061, -844.53893, -826.679817, -808.483705, -788.47756, -765.188346, -730.298793, -691.917793, -651.310007, -609.740092, -571.99583, -534.409508, -496.836536, -459.132322, -420.424886, -381.587985, -342.767985, -304.111253, -266.47896, -229.016746, -191.585055, -154.044331, -116.18712, -77.9689255, -39.2773507, 0., 39.3956861, 79.8338744, 121.658895, 165.21508, 209.062449, 256.043366, 307.215884, 363.638057, 420.173885, 486.553096, 566.311367, 662.984373, 757.078737, 884.370809, 1057.60788, 1289.53726, 1499.58688, 1831.15114, 2334.30506, 3059.12368, 3971.24838, 5374.05517, 7361.11611, 10257.2967, 13925.2237, 19845.1853, 29093.0545, 42744.7041}, 
{-1376.76382, -1326.68783, -1282.07861, -1241.75639, -1204.5414, -1169.25385, -1134.71398, -1099.742, -1063.15815, -1015.4217, -967.058202, -920.232267, -877.1085, -849.472809, -826.019977, -805.066087, -784.927224, -759.752074, -733.691079, -706.727285, -678.843735, -650.016673, -620.238666, -589.495479, -557.772878, -524.522938, -490.478592, -455.839083, -420.803651, -385.914211, -350.890265, -315.793989, -280.687555, -245.903895, -211.126125, -176.308115, -141.403738, -106.570871, -71.4777782, -35.9967311, 0., 36.0990347, 73.185657, 111.604041, 151.698359, 192.269829, 235.822765, 283.318523, 335.718461, 388.430658, 450.191058, 524.182331, 613.587143, 698.997506, 815.223007, 974.482576, 1188.99514, 1378.73642, 1685.06585, 2157.09963, 2843.95399, 3710.28866, 5058.5893, 6988.85511, 9802.94449, 13373.0307, 19088.1684, 27958.192, 40992.9357}, 
{-1237.58523, -1193.72216, -1155.05249, -1120.43609, -1088.73286, -1058.80268, -1029.50544, -999.701022, -968.249314, -926.239712, -883.41079, -841.730634, -803.167327, -778.318849, -757.071429, -737.941191, -719.444261, -696.334781, -672.39565, -647.647786, -622.112106, -595.968306, -569.015013, -541.209634, -512.509574, -482.320899, -451.372892, -419.843492, -387.910642, -356.060285, -324.03916, -291.902006, -259.703565, -227.751645, -195.74669, -163.642216, -131.391734, -99.0699117, -66.4606474, -33.4689929, 0., 33.4513572, 67.8059178, 103.394599, 140.548317, 178.214416, 218.660816, 262.771863, 311.431906, 360.448557, 417.813589, 486.442042, 569.248956, 647.275805, 754.06062, 901.267868, 1100.56201, 1272.16848, 1555.76639, 1999.59584, 2651.8969, 3474.58878, 4774.87422, 6664.71547, 9415.52636, 12887.0161, 18426.2608, 27008.4947, 39608.9522}, 
{-1105.58134, -1068.28477, -1035.7158, -1006.80631, -980.488224, -955.693435, -931.35385, -906.401371, -879.767902, -843.331603, -805.899619, -769.225353, -735.062207, -712.723629, -693.378956, -675.757571, -658.588856, -637.292406, -615.231308, -592.458861, -569.028363, -545.287037, -520.876688, -495.733044, -469.791833, -442.442738, -414.385952, -385.775624, -356.765902, -327.780257, -298.595783, -269.258899, -239.816024, -210.549023, -181.174687, -151.645255, -121.912965, -91.969331, -61.7116068, -31.0763203, 0., 30.9521313, 62.7237173, 95.6297075, 129.985051, 164.863675, 202.31796, 243.159265, 288.198947, 333.654277, 386.768336, 450.190117, 526.568614, 597.431846, 694.998168, 830.364963, 1014.62961, 1168.42203, 1429.49406, 1845.13006, 2462.6144, 3241.06015, 4492.26561, 6339.07397, 9025.31546, 12395.8862, 17770.4317, 26099.91, 38335.2793}, 
{-981.843359, -951.135893, -924.577433, -901.1915, -880.001617, -860.031308, -840.304093, -819.843497, -797.673041, -766.547528, -734.266689, -702.361533, -672.363072, -652.260502, -634.543371, -618.159416, -602.056372, -582.353316, -561.958105, -540.949937, -519.408011, -497.806608, -475.671811, -452.924788, -429.486703, -404.761057, -379.393749, -353.513012, -327.24708, -300.95003, -274.433913, -247.736624, -220.896058, -194.167337, -167.28424, -140.19777, -112.858933, -85.1877245, -57.178561, -28.7948502, 0., 28.5987857, 57.9394136, 88.3159941, 120.022638, 152.2321, 186.808389, 224.494156, 266.032056, 308.046801, 357.046158, 415.419955, 485.558019, 549.519157, 638.156626, 761.992664, 931.549505, 1067.95242, 1306.87941, 1694.61147, 2277.42964, 3011.70205, 4213.44837, 6014.88302, 8636.06863, 11907.0428, 17130.0951, 25239.4523, 37169.3408}, 
{-867.462487, -843.035749, -822.146291, -803.916104, -787.467177, -771.921499, -756.401059, -740.027848, -721.923855, -695.737645, -668.254002, -640.784284, -614.639853, -596.502818, -580.165486, -564.790917, -549.542169, -531.245875, -512.336089, -492.910439, -473.066553, -453.360756, -433.248503, -412.643943, -391.461228, -369.148461, -346.272259, -322.93319, -299.231825, -275.445509, -251.427325, -227.207134, -202.814794, -178.477895, -153.949475, -129.180301, -104.121142, -78.6436884, -52.8094149, -26.6007194, 0., 26.3887492, 53.4533645, 81.4600863, 110.675155, 140.334184, 172.146291, 206.789968, 244.943705, 283.625111, 328.637912, 382.124952, 446.229074, 503.591266, 583.656969, 696.36977, 851.673255, 971.215005, 1188.55301, 1548.94926, 2097.66575, 2788.51376, 3941.10738, 5695.09497, 8251.54269, 11427.8877, 16514.6653, 24434.1361, 36108.5609}, 
{-763.529923, -744.744569, -728.931269, -715.304576, -703.079038, -691.469208, -679.689636, -666.954874, -652.479471, -630.752109, -607.603558, -584.138715, -561.462481, -545.02393, -529.846115, -515.296264, -500.741608, -483.698449, -466.12531, -448.129792, -429.819492, -411.783222, -393.454883, -374.749589, -355.582454, -335.477555, -314.897458, -293.913693, -272.597788, -251.142599, -229.449798, -207.542384, -185.443357, -163.352002, -141.044517, -118.473386, -95.5910942, -72.2558187, -48.5520731, -24.4700646, 0., 24.3194502, 49.2659279, 75.0686116, 101.95668, 129.18442, 158.345858, 190.060128, 224.946366, 260.38819, 301.534457, 350.298506, 408.59368, 459.701698, 531.62017, 633.715083, 775.352427, 878.665118, 1075.14545, 1409.05263, 1924.64589, 2573.49455, 3677.92749, 5382.66222, 7875.49447, 10965.8228, 15933.5562, 23690.976, 35150.3634}, 
{-671.136869, -657.022582, -645.441263, -635.681364, -627.031337, -618.779636, -610.214713, -600.625022, -589.299015, -571.441079, -552.057359, -532.069936, -512.400888, -497.39719, -483.186069, -469.319648, -455.35005, -439.439403, -423.08582, -406.397423, -389.482331, -372.907744, -356.139071, -339.100802, -321.717425, -303.620942, -285.145324, -266.332054, -247.222617, -227.917203, -208.375107, -188.61433, -168.652875, -148.660967, -128.443493, -107.957565, -87.1602927, -65.9427114, -44.3544406, -22.3790227, 0., 22.3883175, 45.3774618, 69.1481977, 93.8812897, 118.797302, 145.421279, 174.318067, 206.052511, 238.33502, 275.726649, 319.934016, 372.66374, 417.903983, 482.167203, 574.247404, 702.938588, 790.758109, 967.287275, 1275.83075, 1759.69318, 2368.6437, 3426.5936, 5080.53715, 7511.68081, 10528.25, 15396.182, 23016.9864, 34292.1725}, 
{-591.374525, -580.630018, -572.185168, -565.370918, -559.518208, -553.957981, -548.021179, -541.038742, -532.341613, -517.654709, -501.357407, -484.223055, -467.025006, -453.195949, -439.786162, -426.505259, -413.062857, -398.197102, -382.977668, -367.502757, -351.870573, -336.568062, -321.149189, -305.556662, -289.733187, -273.451226, -256.891831, -240.065808, -222.983963, -205.645227, -188.07703, -170.294927, -152.314472, -134.276096, -116.020528, -97.5133744, -78.7202402, -59.6229627, -40.164422, -20.3037303, 0., 20.5927796, 41.7883242, 63.7054721, 86.463062, 109.187321, 133.386744, 159.577213, 188.274612, 217.464584, 251.205346, 291.024877, 338.451155, 378.251645, 435.419044, 518.185534, 634.783301, 707.949324, 865.609071, 1150.19281, 1604.13079, 2175.96048, 3189.79057, 4791.67214, 7163.85853, 10122.571, 14911.9569, 22419.182, 33531.4121}, 
{-531.224895, -520.626258, -512.797608, -506.999762, -502.493538, -498.539755, -494.399231, -489.332783, -482.601229, -469.841541, -455.387921, -439.950727, -424.240315, -411.213282, -398.435247, -385.718072, -372.873616, -359.085644, -345.04535, -330.815832, -316.46019, -302.247325, -287.952211, -273.555623, -259.038339, -244.396642, -229.589599, -214.591784, -199.377772, -183.808107, -168.017004, -152.02465, -135.851231, -119.609428, -103.189935, -86.5759385, -69.7506277, -52.91163, -35.741918, -18.1389039, 0., 18.9692432, 38.6028073, 58.9265352, 79.9662703, 100.621877, 122.495568, 146.063579, 171.802143, 197.878892, 228.002105, 263.571461, 305.986636, 341.038107, 391.978429, 466.450961, 572.099059, 631.266927, 771.016734, 1033.1115, 1459.31424, 1998.35806, 2971.09571, 4515.64123, 6831.24812, 9763.79368, 14511.5205, 21934.5798, 32893.1229}, 
{-486.583679, -473.979814, -465.030685, -458.861989, -454.599418, -451.368669, -448.295435, -444.505412, -439.124294, -427.323975, -413.765471, -399.155998, -384.202771, -371.710721, -359.450261, -347.28952, -335.096629, -322.375158, -309.503618, -296.495961, -283.36614, -270.088482, -256.732415, -243.32774, -229.904261, -216.721713, -203.487992, -190.140926, -176.618344, -162.629655, -148.432476, -134.056005, -119.529439, -104.937331, -90.2313794, -75.4186386, -60.5061629, -45.9950208, -31.2006471, -15.9324906, 0., 17.4807597, 35.7296022, 54.6597253, 74.1843266, 92.8924074, 112.55104, 133.603101, 156.491467, 179.486059, 206.071889, 237.561016, 275.265497, 306.105316, 351.541434, 418.642738, 514.478116, 560.661054, 683.898001, 925.440002, 1326.53811, 1837.0297, 2772.40681, 4258.3764, 6521.60732, 9452.61168, 14185.5666, 21543.733, 32350.3718}, 
{-448.814905, -434.351947, -424.23195, -417.480755, -413.1242, -410.188124, -407.698365, -404.680764, -400.161159, -388.964091, -375.997217, -361.966893, -347.579478, -335.550049, -323.772757, -312.150471, -300.586061, -288.808322, -276.963829, -265.025081, -252.964578, -240.504183, -227.967286, -215.426643, -202.955006, -191.03281, -179.162057, -167.252429, -155.213609, -142.637141, -129.878101, -116.973426, -103.960054, -90.8894741, -77.7782523, -64.6575056, -51.5583513, -39.2929482, -26.798955, -13.7950721, 0., 16.0603954, 32.9974454, 50.6153159, 68.7181726, 85.5958004, 103.172498, 121.858185, 142.062778, 162.114943, 185.338352, 212.975426, 246.268582, 273.110677, 313.433519, 373.819356, 460.850433, 495.64119, 604.4288, 827.982629, 1207.07205, 1692.46558, 2594.93511, 4028.40848, 6246.18361, 9183.8677, 13908.4607, 21204.1145, 31854.9808}, 
{-415.107526, -399.828275, -389.110135, -381.960008, -377.384796, -374.391399, -371.986718, -369.177655, -364.971111, -354.165601, -341.659766, -328.143861, -314.308141, -302.780952, -291.539221, -280.497968, -269.57221, -258.587577, -247.584233, -236.512952, -225.324507, -213.595855, -201.801116, -190.040593, -178.414589, -167.532131, -156.781306, -146.058926, -135.261805, -123.926816, -112.454684, -100.886198, -89.2621454, -77.5906511, -65.9582303, -54.4187358, -43.0260202, -32.8911881, -22.5879395, -11.7472262, 0., 14.6788678, 30.3399124, 46.6897497, 63.4349954, 78.6160763, 94.2722731, 110.776677, 128.50238, 145.787851, 165.854652, 189.889726, 219.080012, 241.990545, 277.478935, 331.780884, 411.132096, 436.62174, 533.690667, 742.633193, 1103.74364, 1567.96956, 2443.64553, 3835.9396, 6019.38928, 8970.19375, 13688.5361, 20918.6033, 31404.5825}, 
{-381.404703, -367.810516, -358.064516, -351.312694, -346.701042, -343.375552, -340.482213, -337.167018, -332.575957, -322.140584, -310.207104, -297.407282, -284.372885, -273.56541, -263.055001, -252.741533, -242.524879, -232.132597, -221.705805, -211.213306, -200.6239, -189.585152, -178.515596, -167.512531, -156.673252, -146.561658, -136.621805, -126.76435, -116.899948, -106.625098, -96.2902787, -85.9318111, -75.5860158, -65.1968036, -54.9298686, -44.858495, -35.0559672, -26.8956431, -18.630704, -9.81440464, 0., 13.2974927, 27.6712859, 42.7528299, 58.1735747, 71.8319533, 85.7856399, 100.359291, 115.877565, 130.624358, 147.781391, 168.489624, 193.890017, 212.684812, 243.429178, 292.239562, 365.232414, 384.336086, 473.530361, 672.606928, 1021.35747, 1469.13545, 2327.04725, 3698.40244, 5866.54264, 8835.96394, 13546.5125, 20704.1969, 31015.0258}, 
{-350.693994, -339.567532, -331.244199, -325.066325, -320.376238, -316.516268, -312.828743, -308.655992, -303.340345, -293.184706, -281.786597, -269.704117, -257.495364, -247.430419, -237.670605, -228.089227, -218.55959, -208.634708, -198.636295, -188.565772, -178.424563, -168.013895, -157.615462, -147.310764, -137.1813, -127.67297, -118.357114, -109.169472, -100.045784, -90.6878242, -81.3588831, -72.0882854, -62.9053565, -53.6877426, -44.6771193, -35.9634834, -27.6368316, -21.2984096, -14.9224656, -7.99449665, 0., 11.9307478, 25.0149408, 38.8249936, 52.933321, 65.1521334, 77.5181313, 90.3078109, 103.797669, 116.194025, 130.671622, 148.335027, 170.288808, 185.113374, 211.447112, 255.404249, 323.099016, 337.713944, 421.467639, 613.647007, 953.538957, 1388.60832, 2233.60823, 3592.34408, 5751.29326, 8739.15611, 13433.0691, 20498.0598, 30599.1557}, 
{-322.965118, -314.722469, -308.053627, -302.529432, -297.72072, -293.198329, -288.533096, -283.29586, -277.057457, -267.132496, -256.250535, -244.884904, -233.50893, -224.173201, -215.142885, -206.260408, -197.368196, -187.807977, -178.123154, -168.356433, -158.550519, -148.713216, -138.936092, -129.275811, -119.789041, -110.754011, -101.917197, -93.246636, -84.7103687, -76.143424, -67.7000552, -59.4015055, -51.2690183, -43.1182337, -35.2582393, -27.7925197, -20.8245593, -16.1400159, -11.4873314, -6.29712101, 0., 10.5704871, 22.3478263, 34.8625739, 47.6452865, 58.4725586, 69.3304936, 80.4512328, 92.0669177, 102.300478, 114.336951, 129.252163, 148.121941, 159.204636, 181.520537, 221.27246, 284.663221, 296.452606, 376.86367, 564.67644, 898.670942, 1324.54377, 2160.32525, 3511.71103, 5663.53469, 8666.0039, 13328.6919, 20271.3306, 30113.6518}, 
{-298.207795, -292.898473, -287.897243, -283.010547, -278.044824, -272.806515, -267.102061, -260.737901, -253.520476, -243.818486, -233.451208, -222.800179, -212.246936, -203.590978, -195.228694, -186.974434, -178.64255, -169.366469, -159.913835, -150.371369, -140.825792, -131.514249, -122.312866, -113.248194, -104.346783, -95.6927269, -87.2320163, -78.9681837, -70.9047623, -63.0203273, -55.353353, -47.9173557, -40.7258518, -33.5430427, -26.7314855, -20.4044225, -14.675096, -11.4609901, -8.34940844, -4.73189658, 0., 9.2085648, 19.6468917, 30.8219035, 42.2405231, 51.6891709, 61.0834733, 70.6185542, 80.4895374, 88.7473442, 98.5889822, 111.067257, 127.234972, 134.887001, 153.637255, 189.841712, 249.85635, 260.249361, 339.079625, 524.618236, 855.136285, 1275.09739, 2104.19508, 3450.4498, 5593.16045, 8602.74098, 13213.8667, 19995.1477, 29515.1939}, 
{-277.692073, -273.266665, -268.661563, -263.791764, -258.572268, -252.918071, -246.744171, -239.965567, -232.497256, -223.318657, -213.654579, -203.794253, -194.026909, -185.767148, -177.728681, -169.75059, -161.671955, -152.674565, -143.517713, -134.303398, -125.133621, -116.351293, -107.721136, -99.2487811, -90.9398631, -82.7556789, -74.7639327, -66.9879926, -59.4512265, -52.2008259, -45.2268059, -38.543005, -32.163262, -25.856573, -19.9795559, -14.6439865, -9.96164035, -7.81861146, -5.84262963, -3.43574307, 0., 7.74471995, 16.6413325, 26.2149221, 35.990573, 43.8290239, 51.5854432, 59.4506533, 67.6154769, 74.1745522, 82.25336, 92.8811966, 107.087358, 112.157665, 128.36228, 162.227891, 220.281186, 230.228588, 308.945153, 494.485671, 824.904937, 1246.65773, 2072.59887, 3408.55568, 5528.94543, 8536.44835, 13052.1028, 19567.3577, 28573.6617}, 
{-259.939683, -255.968472, -251.491441, -246.505072, -241.005847, -234.990248, -228.454756, -221.395854, -213.810024, -205.1907, -196.23863, -187.151517, -178.127059, -170.170813, -162.349484, -154.537632, -146.609818, -137.847271, -128.955213, -120.045537, -111.230135, -102.937193, -94.8357926, -86.9113071, -79.1491114, -71.3967003, -63.8324798, -56.4969763, -49.4307162, -42.7343354, -36.3642068, -30.3368131, -24.6686367, -19.1508851, -14.1154261, -9.66885215, -5.91775591, -4.66091588, -3.63586439, -2.27231968, 0., 6.27642281, 13.5921656, 21.5074915, 29.5826637, 35.7843028, 41.9037718, 48.1387909, 54.6870805, 59.7199925, 66.272163, 75.3518594, 87.9673494, 90.8661607, 105.021597, 137.146221, 193.952597, 204.451318, 284.137707, 470.805114, 802.246893, 1227.94069, 2050.62697, 3373.0398, 5462.06578, 8452.38814, 12843.3511, 19043.1335, 27459.9144}, 
{-244.613452, -240.742451, -236.178876, -230.976389, -225.188652, -218.869326, -212.072073, -204.850553, -197.258429, -189.207498, -180.950031, -172.596435, -164.257118, -156.526127, -148.836772, -141.106004, -133.250774, -124.687935, -116.034573, -107.407678, -98.9242389, -91.0767404, -83.4564783, -76.0302437, -68.7648278, -61.4070679, -54.2316906, -47.2934683, -40.6471739, -34.4312784, -28.5833766, -23.1247619, -18.0767276, -13.2663693, -8.98685713, -5.33716324, -2.41625999, -1.89649508, -1.67411532, -1.2187429, 0., 4.81904605, 10.5395963, 16.7694066, 23.1162328, 27.6727934, 32.1678966, 36.8153135, 41.8288147, 45.4844063, 50.7087301, 58.4906628, 69.8190813, 70.9233303, 83.4556316, 114.308675, 170.37515, 182.25048, 263.643527, 451.965888, 784.629157, 1214.99982, 2032.62481, 3336.30428, 5382.70068, 8336.57764, 12572.4318, 18412.3265, 26178.3252}, 
{-231.376208, -227.327158, -222.515866, -217.031636, -210.963774, -204.401586, -197.434378, -190.151455, -182.642124, -175.141937, -167.535452, -159.853477, -152.126817, -144.557247, -136.936217, -129.226146, -121.389454, -112.999904, -104.564031, -96.1997155, -88.0248379, -80.5747267, -73.3828348, -66.400063, -59.5773124, -52.5780585, -45.755598, -39.1758021, -32.9045422, -27.1020775, -21.7021359, -16.7328335, -12.2222863, -8.04341573, -4.44161024, -1.50706368, 0.670030087, 0.566059177, 0.0976149603, -0.252129215, 0., 3.38796234, 7.52383021, 12.0704624, 16.6907177, 19.6122814, 22.5072553, 25.6125674, 29.1651461, 31.5685348, 35.6264003, 42.3090248, 52.5866904, 52.2400158, 63.5048119, 93.4272263, 149.053406, 162.959, 246.448854, 436.357316, 769.518734, 1203.88868, 2012.93783, 3290.75124, 5281.02929, 8175.03416, 12224.1654, 17664.7881, 24733.2674}, 
{-219.890779, -215.461152, -210.294408, -204.496731, -198.174303, -191.433308, -184.37993, -177.120352, -169.760758, -162.766899, -155.741564, -148.647109, -141.445889, -133.988329, -126.393489, -118.6685, -110.820491, -102.586526, -94.3518247, -86.2315432, -78.3408367, -71.2359444, -64.414504, -57.8152372, -51.3768655, -44.7009486, -38.1982351, -31.9423113, -26.0067637, -20.5571553, -15.5383054, -10.9870099, -6.94006466, -3.32241444, -0.327446578, 1.96330266, 3.46829701, 2.81815514, 1.73432386, 0.650404896, 0., 1.99854432, 4.58507283, 7.48045392, 10.405556, 11.7205528, 13.0512852, 14.6628993, 16.8205414, 18.073119, 21.0885124, 26.8183632, 36.2143128, 34.727059, 45.009565, 74.2138497, 129.491933, 145.909806, 231.539927, 422.368724, 754.382627, 1190.66081, 1985.91147, 3228.78282, 5147.23075, 7953.77503, 11783.3722, 16790.3694, 23129.1142}, 
{-226.724899, -212.764992, -200.613432, -190.033165, -180.787138, -172.638298, -165.34959, -158.683961, -152.404358, -145.81283, -139.317579, -132.865911, -126.40513, -119.900384, -113.273999, -106.466146, -99.4169927, -91.5262098, -83.4906648, -75.4667261, -67.6107623, -60.4573783, -53.6334118, -47.1439369, -40.9940276, -35.2433287, -29.8205152, -24.7088327, -19.8915271, -15.2693682, -10.9410679, -6.92286233, -3.23098758, 0.61081133, 3.89681083, 6.41377843, 7.9484816, 6.40498557, 4.20584103, 1.89089637, 0., 1.36847131, 3.32249934, 5.4837441, 7.4738656, 7.10299269, 6.52892902, 6.09794707, 6.15631933, 4.78894311, 5.50801613, 9.56436095, 18.2088001, 19.6826667, 33.4500685, 65.965624, 123.683951, 138.560516, 221.348751, 408.302935, 735.677349, 1179.71831, 1956.70398, 3138.69873, 4942.56472, 7595.34243, 11138.5612, 15660.1283, 21247.9508}, 
{-236.84171, -211.295632, -190.740811, -174.495373, -161.877439, -152.205134, -144.796581, -138.969903, -134.043225, -127.670992, -121.500477, -115.515274, -109.698978, -104.443078, -99.1601186, -93.6705368, -87.7947709, -80.3881138, -72.6222059, -64.7035429, -56.8386203, -49.4919717, -42.5088397, -35.9925048, -30.0462478, -25.3821555, -21.2511794, -17.5130774, -14.0276074, -10.327812, -6.73085042, -3.2271669, 0.192794552, 4.54191817, 8.42510028, 11.4505655, 13.2265385, 10.5952775, 7.03736001, 3.26739713, 0., 1.10351971, 2.87733143, 4.77455029, 6.24829142, 4.25773493, 1.74750499, -0.831709239, -3.02921862, -7.34619724, -9.19934742, -6.95723472, 1.01157529, 6.82986723, 25.4411853, 62.2804239, 122.782477, 135.652615, 213.747207, 393.192996, 710.116727, 1156.14026, 1900.90499, 2996.52076, 4644.78844, 7069.78957, 10256.7919, 14250.7667, 19096.6851}, 
{-232.280967, -202.805662, -179.493136, -161.488734, -147.937801, -137.985681, -130.77772, -125.459263, -121.175653, -114.86966, -108.770235, -102.903752, -97.2965875, -92.6305548, -88.0144139, -83.2123638, -77.9886036, -70.9785561, -63.5267075, -55.8487674, -48.1604456, -40.8346679, -33.8670415, -27.4103897, -21.6175357, -17.5081853, -14.0215263, -10.963629, -8.14056397, -4.8928611, -1.67834712, 1.51069153, 4.68196843, 9.02829042, 12.8982405, 15.8254948, 17.3437298, 13.841766, 9.25607751, 4.37828269, 0., 0.430831866, 1.53721913, 2.70358662, 3.31435919, -0.0847916974, -4.13518649, -8.31689902, -12.1100031, -18.0667268, -21.3661281, -20.2594189, -12.9978114, -5.65089228, 14.4758508, 52.2585556, 112.57336, 127.359326, 203.604497, 375.359842, 676.676327, 1108.31316, 1804.19659, 2804.49135, 4284.70754, 6449.15256, 9251.59936, 12708.4545, 16836.1246}, 
{-217.723575, -189.98152, -168.069873, -151.165799, -138.446465, -129.089038, -122.270683, -117.168567, -112.959856, -106.589909, -100.360423, -94.3412868, -88.602391, -83.9929754, -79.4918381, -74.8571281, -69.8469946, -63.1281994, -55.9868336, -48.6176009, -41.2152054, -34.0498118, -27.2104788, -20.8617257, -15.168072, -11.1893742, -7.83668007, -4.91637443, -2.23484203, 0.921668114, 4.01858121, 7.04145819, 9.97586002, 13.921929, 17.3048121, 19.6642379, 20.5399346, 16.3430655, 10.9933504, 5.28194367, 0., -0.573489477, -0.530341535, -0.4745363, -1.0100539, -5.63864982, -10.9074187, -16.2612305, -21.1449552, -27.7808742, -31.7254815, -31.3126825, -24.8763827, -18.1852808, 0.835428433, 36.825757, 94.4257167, 113.608134, 189.54708, 352.749442, 633.722106, 1037.49814, 1671.00588, 2568.45021, 3871.51835, 5752.88672, 8153.79254, 11073.3328, 14510.6046}, 
{-197.850438, -175.509646, -157.670492, -143.67912, -132.881674, -124.624299, -118.253141, -113.114343, -108.55405, -102.012813, -95.5046078, -89.1378169, -83.0208222, -78.0605024, -73.2473444, -68.3703316, -63.2184475, -56.6677061, -49.7852481, -42.7252446, -35.6418668, -28.701748, -22.0416126, -15.8106468, -10.1580371, -5.99367853, -2.40176533, 0.772799518, 3.68511306, 7.00195942, 10.1620738, 13.1158783, 15.8137953, 19.0748347, 21.6333959, 23.0924661, 23.0550322, 18.2977904, 12.3805353, 6.03677068, 0., -1.83334156, -3.15785451, -4.50520788, -6.40707065, -12.1179021, -18.3574648, -24.5683118, -30.1929959, -36.8968677, -41.0105633, -41.0875167, -35.6811617, -31.1989674, -15.1987182, 16.907766, 69.7086655, 94.3265236, 170.201412, 323.307766, 579.620021, 944.956329, 1505.75997, 2294.23711, 3414.41715, 5000.44737, 6994.1804, 9385.54253, 12164.4601}, 
{-177.34246, -162.076478, -149.494462, -139.181246, -130.721668, -123.700563, -117.702768, -112.313118, -107.11645, -100.319447, -93.4363589, -86.6032825, -79.9563146, -74.3632984, -68.935886, -63.5174757, -57.9514658, -51.4277386, -44.704615, -37.8868996, -31.0793971, -24.3548212, -17.8629043, -11.7212873, -6.04761157, -1.48905451, 2.57809353, 6.24800591, 9.61485598, 13.2341966, 16.5542697, 19.4846971, 21.9351002, 24.3390085, 25.8725726, 26.2358508, 25.1289017, 19.9045551, 13.5489893, 6.70115434, 0., -3.27262161, -6.17782377, -9.13381751, -12.5588139, -19.2366111, -26.2735981, -33.1417511, -39.3130462, -45.8229357, -49.9545293, -50.5544126, -46.4691713, -45.1176209, -33.3452256, -6.56967931, 39.791324, 69.4419812, 144.193951, 284.980784, 512.73603, 831.948853, 1312.88596, 1987.69175, 2922.60026, 4211.28985, 5803.57192, 7685.22462, 9842.0261}, 
{-160.880546, -152.368455, -144.741251, -137.82473, -131.444689, -125.426926, -119.597237, -113.78142, -107.805271, -100.690886, -93.3892444, -86.0476237, -78.8133015, -72.4315259, -66.2124162, -60.0640621, -53.8945533, -47.2389595, -40.5275983, -33.8177671, -27.1667634, -20.5733762, -14.1768151, -8.05778129, -2.29697582, 2.75654163, 7.39777206, 11.6533579, 15.5499414, 19.5045633, 22.9973084, 25.8986598, 28.0791008, 29.5664513, 30.0109229, 29.2200637, 27.0014222, 21.3619739, 14.6300691, 7.33348529, 0., -4.81522689, -9.42275327, -14.1057546, -19.1474064, -26.7088395, -34.4440918, -41.8851566, -48.5640273, -54.9673069, -59.290535, -60.683861, -58.2974344, -60.3669104, -53.3227302, -32.6808411, 6.04281005, 38.8819919, 110.151155, 235.714466, 431.436092, 699.736845, 1096.81095, 1654.65389, 2405.26399, 3404.86948, 4612.77606, 6012.52008, 7587.63786}, 
{-117.464876, -123.918365, -127.492752, -128.518553, -127.326281, -124.246452, -119.609581, -113.746181, -106.986767, -100.040866, -92.7083752, -85.1682045, -77.5992633, -70.6372629, -63.8215906, -57.1484355, -50.6139865, -44.2722108, -38.038408, -31.8856556, -25.7870314, -19.5202853, -13.3319539, -7.27324596, -1.3953703, 4.18790751, 9.51295815, 14.5535955, 19.2836336, 24.2243411, 28.5830951, 32.1147276, 34.5740706, 35.1779179, 34.4343549, 32.3134286, 28.7851863, 22.7771943, 15.718973, 7.99756168, 0., -6.61155838, -13.2352339, -19.9940327, -27.010961, -35.9903708, -44.8413841, -53.054469, -60.1200932, -64.6781795, -67.4099591, -68.1461179, -66.7173418, -74.5060328, -75.1704745, -63.9206668, -35.9666093, 0.717528466, 68.7915847, 178.151228, 338.692126, 545.47663, 852.90036, 1304.26458, 1880.44809, 2572.04849, 3393.21069, 4352.19278, 5457.25286}, 
{-136.527966, -139.159539, -139.607421, -138.121857, -134.953088, -130.351358, -124.566911, -117.849989, -110.450837, -102.93602, -95.112929, -87.1052792, -79.0367843, -71.261443, -63.5805711, -56.0257688, -48.6286365, -41.5501343, -34.6407587, -27.880366, -21.2488126, -14.5876696, -8.07039264, -1.7321522, 4.3918813, 10.3326244, 15.962385, 21.2195577, 26.0425374, 30.7830152, 34.8007711, 37.8688812, 39.7604217, 39.674461, 38.1876865, 35.3027776, 31.0224139, 24.5660901, 17.0329442, 8.73892973, 0., -8.04750415, -16.2361718, -24.5782045, -33.0858038, -42.6581978, -52.0657511, -60.9658549, -69.0159002, -75.50637, -80.6083267, -84.1259248, -85.8633187, -94.7693226, -97.845567, -91.2383424, -71.0939391, -40.7141918, 13.7723711, 99.0816769, 221.929653, 380.179382, 607.105325, 934.522007, 1336.34321, 1800.19873, 2316.75508, 2875.91963, 3467.59973}, 
{-194.664249, -181.230941, -169.302868, -158.648225, -149.035212, -140.232023, -132.006856, -124.127907, -116.363375, -107.812862, -99.1805959, -90.502209, -81.8133354, -73.2086007, -64.6410495, -56.1227188, -47.6656453, -39.1592333, -30.7872051, -22.6106505, -14.6906592, -7.17655082, -0.00589349249, 6.79551487, 13.2018764, 19.3079585, 24.9191718, 29.9614923, 34.3608958, 38.1894679, 41.1686312, 43.1659177, 44.0488595, 43.4935761, 41.6355773, 38.4189602, 33.7878219, 26.7947351, 18.6319311, 9.6001171, 0., -9.24386883, -18.7041645, -28.3297179, -38.0693597, -47.9114983, -57.7495558, -67.516532, -77.1454268, -88.0404396, -98.0748906, -106.593299, -112.940186, -120.019097, -122.191915, -117.379548, -103.502907, -87.2876613, -54.3280535, 0.976913021, 84.2282345, 202.753802, 356.973927, 552.751261, 786.612743, 1083.52895, 1364.03162, 1574.7622, 1662.36213}, 
{-175.36164, -166.624748, -158.159801, -149.925768, -141.881618, -133.986319, -126.198841, -118.478152, -110.783222, -103.030633, -95.2386948, -87.3833331, -79.440472, -71.2174429, -62.9261997, -54.6101037, -46.3125163, -38.121428, -30.017719, -22.0268988, -14.1744769, -6.30376836, 1.30464509, 8.55237609, 15.3410373, 21.2894899, 26.6951987, 31.5728771, 35.9372383, 40.4691856, 44.2507664, 47.0302182, 48.5557784, 48.0057096, 45.9262139, 42.2935185, 37.0838506, 29.5326914, 20.6533127, 10.7182403, 0., -10.8380537, -22.0705555, -33.5813111, -45.2541263, -57.0379543, -68.7253946, -80.1741938, -91.2420989, -101.801961, -111.69038, -120.759063, -128.859713, -138.865711, -146.398417, -150.100866, -148.616093, -144.359406, -130.692659, -104.749978, -63.6654913, -2.1174713, 75.3923952, 179.345261, 277.973065, 352.592929, 374.43045, 317.234105, 154.752373}, 
{-191.788584, -182.390037, -173.740918, -165.639879, -157.885566, -150.276631, -142.611723, -134.68949, -126.308584, -116.144987, -105.569081, -94.8285798, -84.1712001, -74.6575672, -65.3973217, -56.3130147, -47.3271972, -37.9729249, -28.7180425, -19.6408993, -10.8198446, -2.39999017, 5.63378214, 13.2298281, 20.3365033, 26.8791912, 32.8384087, 38.1717007, 42.8366119, 47.2220462, 50.6816452, 53.00041, 53.9633413, 52.8498506, 50.1527642, 45.859319, 39.9567519, 31.9102235, 22.4378773, 11.7357803, 0., -12.5765892, -25.7934503, -39.453239, -53.3586108, -66.8576235, -80.38937, -93.9383452, -107.489044, -121.262378, -134.911858, -148.327413, -161.398972, -174.938597, -187.545228, -198.739941, -208.043809, -213.510736, -216.715838, -217.767058, -216.77234, -208.338342, -209.076867, -219.416934, -264.909434, -355.793203, -538.433422, -850.163683, -1328.31758}, 
{-198.71236, -197.308699, -193.101836, -186.510063, -177.951671, -167.844951, -156.608195, -144.659693, -132.417739, -121.285996, -110.303233, -99.4935919, -88.8812147, -78.707498, -68.6924271, -58.7732421, -48.8871829, -38.459504, -28.1442253, -18.0833812, -8.41900594, 0.229792022, 8.38888204, 16.1070594, 23.4331194, 30.88596, 37.8562326, 44.2046916, 49.7920911, 54.7738453, 58.5981843, 61.0079983, 61.7461775, 59.6423142, 55.7179154, 50.0811906, 42.840349, 34.2000096, 24.1334081, 12.7101898, 0., -14.5670711, -30.0260009, -46.0513217, -62.3175661, -77.5333094, -92.7254243, -107.954826, -123.28243, -138.151262, -153.487283, -169.598565, -186.793178, -205.872273, -226.45361, -248.648032, -272.566377, -292.582224, -316.83858, -347.741193, -387.695808, -434.404278, -504.384028, -598.974796, -746.35905, -974.171705, -1285.39517, -1689.17496, -2194.65663}
};

  const double a5tab[50][69] = {
{4296.48547, 4145.14906, 4028.39347, 3938.25617, 3866.77463, 3805.98633, 3747.92872, 3684.63929, 3608.15549, 3470.8602, 3320.30733, 3164.39619, 3011.02609, 2895.53131, 2787.40219, 2683.56406, 2580.94224, 2467.08826, 2352.05074, 2236.50453, 2121.12447, 2008.42576, 1896.50673, 1785.30605, 1674.76241, 1564.89358, 1455.52754, 1346.57134, 1237.93205, 1128.57981, 1019.73336, 911.67453, 804.685139, 699.870801, 596.360042, 494.10517, 393.058497, 293.389367, 194.746242, 96.9946212, 0., -95.8821865, -191.472354, -287.100981, -383.098545, -476.415091, -572.113704, -671.877035, -777.387735, -883.098866, -1000.81451, -1135.10914, -1290.55726, -1445.5122, -1641.25806, -1892.85779, -2215.37433, -2550.46591, -3015.96209, -3656.28769, -4515.86756, -5562.5241, -7070.48939, -8999.14401, -12015.737, -15367.8837, -23273.8659, -38709.2987, -64649.7973}, 
{3989.50138, 3890.74144, 3791.28858, 3691.26149, 3590.77885, 3489.95935, 3388.92168, 3287.78452, 3186.66657, 3088.65256, 2989.7087, 2888.76729, 2784.76057, 2664.7319, 2544.25806, 2427.02688, 2316.72622, 2240.60956, 2169.37286, 2097.27769, 2018.58566, 1902.04814, 1777.64101, 1649.82996, 1523.08064, 1415.81637, 1312.96214, 1213.40058, 1116.01432, 1016.79356, 918.670313, 821.68419, 725.874783, 631.536102, 538.351566, 446.259007, 355.196258, 265.266634, 176.17629, 87.7968658, 0., -86.9254916, -173.692018, -260.594811, -347.929102, -432.365887, -519.274328, -610.399354, -707.48589, -805.717825, -916.025541, -1042.77838, -1190.34568, -1339.01104, -1526.86385, -1767.90776, -2076.14641, -2395.68728, -2838.38865, -3446.21264, -4261.12134, -5252.34273, -6680.04137, -8501.55957, -11361.8428, -14541.9955, -22090.8829, -36865.43, -61722.5616}, 
{3757.56237, 3688.98233, 3602.74083, 3502.53272, 3392.05283, 3274.99599, 3155.05703, 3035.9308, 2921.31211, 2833.32287, 2749.86003, 2667.24761, 2581.80961, 2470.83397, 2357.29521, 2245.13178, 2138.28212, 2062.17793, 1990.6671, 1919.09078, 1842.79012, 1736.93836, 1625.1117, 1510.71846, 1397.16693, 1298.16665, 1202.70421, 1110.06741, 1019.54408, 928.906343, 839.563964, 751.411028, 664.341619, 578.018704, 492.659931, 408.25183, 324.780933, 242.505048, 161.030918, 80.2365616, 0., -79.1939341, -158.316759, -237.733181, -317.807904, -395.541068, -476.007772, -560.918545, -651.983919, -744.762196, -849.577029, -970.599841, -1112.00206, -1254.41911, -1434.97281, -1667.24897, -1964.8334, -2274.39487, -2703.20307, -3291.61063, -4079.9702, -5039.50056, -6417.95592, -8173.42849, -10924.8947, -14014.4248, -21230.8298, -35251.2358, -58752.7686}, 
{3573.43954, 3518.13048, 3441.64113, 3347.99242, 3241.20528, 3125.30065, 3004.29946, 2882.22265, 2763.09114, 2665.78538, 2573.523, 2484.38113, 2396.4369, 2303.02573, 2208.86315, 2113.92298, 2018.17905, 1918.81705, 1819.71417, 1721.95946, 1626.642, 1540.32011, 1456.42589, 1373.86068, 1291.52585, 1204.16451, 1116.49954, 1029.09559, 942.517315, 859.174238, 777.048176, 695.965823, 615.753875, 535.637321, 456.285245, 377.765026, 300.144041, 223.94418, 148.596504, 73.9865859, 0., -72.5271837, -145.040088, -218.033339, -292.001563, -364.645482, -440.371185, -520.79086, -607.516694, -696.21608, -796.823916, -913.330305, -1049.72535, -1186.00008, -1359.74331, -1584.54477, -1873.9942, -2177.10138, -2597.86799, -3175.71572, -3950.0663, -4894.3482, -6245.96286, -7966.33323, -10645.6449, -13703.476, -20613.072, -33834.2617, -55826.8739}, 
{3409.90395, 3356.44466, 3286.88038, 3203.56312, 3108.8449, 3005.07774, 2894.61363, 2779.80461, 2663.00268, 2546.95431, 2433.45928, 2324.71182, 2222.90615, 2150.49542, 2081.31138, 2009.44471, 1928.98609, 1797.5506, 1660.29474, 1525.89938, 1403.04545, 1339.4171, 1289.09066, 1245.14575, 1200.66196, 1126.02999, 1046.0939, 963.008853, 878.930002, 801.853306, 725.756809, 650.459352, 575.779776, 500.710664, 426.227618, 352.479982, 279.6171, 208.423598, 138.159426, 68.7198163, 0., -66.7649103, -133.555516, -201.012536, -269.776693, -338.38397, -410.421719, -487.372554, -570.71909, -656.063581, -753.121147, -865.726549, -997.714547, -1128.01761, -1295.3337, -1513.4585, -1796.18769, -2094.31954, -2509.84608, -3081.76197, -3849.06182, -4787.23624, -6125.79199, -7831.85622, -10464.8454, -13527.4537, -20156.9751, -32582.0535, -53031.3329}, 
{3239.72667, 3182.18361, 3117.34947, 3045.16737, 2965.58039, 2878.53165, 2783.96423, 2681.82123, 2572.04576, 2437.74389, 2302.43055, 2172.78365, 2055.48109, 2002.43127, 1956.9894, 1907.74119, 1843.27233, 1685.40227, 1516.18949, 1350.92621, 1204.90465, 1161.45301, 1140.61312, 1130.46279, 1119.07984, 1055.98318, 983.233088, 904.330944, 822.778121, 751.19961, 680.323725, 610.002396, 540.087553, 469.557447, 399.487158, 330.078085, 261.53163, 194.782874, 129.006064, 64.1091301, 0., -61.7467835, -123.556551, -186.188021, -250.399908, -315.46138, -384.216524, -458.019878, -538.225981, -620.2888, -713.823669, -822.545353, -950.168616, -1074.73535, -1235.90234, -1447.6535, -1723.97275, -2016.56203, -2426.60002, -2992.98342, -3754.60894, -4688.5153, -6019.17312, -7721.57993, -10323.2482, -13404.6625, -19781.9045, -31462.1566, -50452.6013}, 
{3035.67878, 2973.60609, 2911.93932, 2848.72771, 2782.02045, 2709.86678, 2630.31591, 2541.41704, 2441.2194, 2299.06835, 2153.19849, 2013.14059, 1888.4254, 1848.02151, 1818.24671, 1784.85662, 1733.60687, 1569.39574, 1391.17913, 1217.05561, 1065.12374, 1033.65155, 1028.50034, 1035.70092, 1041.28404, 986.24413, 919.662878, 845.585609, 768.05765, 701.469213, 635.382783, 569.705732, 504.34543, 438.496382, 373.063973, 308.240723, 244.21915, 181.861577, 120.422796, 59.8274049, 0., -57.3124734, -114.736706, -173.077041, -233.137821, -294.582554, -359.81275, -430.089086, -506.672238, -584.87584, -674.286427, -778.543494, -901.286533, -1020.41696, -1175.60758, -1380.79311, -1649.90829, -1934.34159, -2335.59247, -2892.61414, -3644.35983, -4568.53596, -5887.83607, -7587.08678, -10161.6055, -13253.4069, -19407.2257, -30442.1168, -48177.1347}, 
{2770.53136, 2708.97085, 2649.54083, 2590.16668, 2528.77379, 2463.28755, 2391.63335, 2311.73658, 2221.52263, 2091.84191, 1958.52477, 1830.32661, 1716.0028, 1676.45437, 1647.43281, 1616.83524, 1572.5588, 1436.55469, 1289.04433, 1144.30323, 1016.60689, 983.236401, 970.259435, 966.74926, 961.779143, 909.032918, 847.129054, 779.296591, 708.764568, 646.91818, 585.567846, 524.68014, 464.221637, 403.846181, 343.958173, 284.649282, 226.011177, 168.499276, 111.696, 55.5475182, 0., -53.3016495, -106.789491, -161.196848, -217.257045, -274.452337, -335.267548, -400.93643, -472.692737, -545.808803, -629.864368, -728.477751, -845.267274, -959.326113, -1108.60779, -1306.54069, -1566.55318, -1838.17091, -2224.28609, -2763.88818, -3495.96666, -4397.64884, -5693.51063, -7379.95924, -9920.66925, -12991.9915, -18952.3044, -29489.4796, -46291.3886}, 
{2417.05548, 2366.53666, 2309.0449, 2245.40682, 2176.44909, 2102.99835, 2025.88123, 1945.9244, 1863.95449, 1776.9788, 1691.17107, 1608.88567, 1532.477, 1476.91808, 1426.8972, 1379.72126, 1332.69721, 1273.9028, 1213.56578, 1152.68474, 1092.25826, 1037.43127, 983.397481, 929.496943, 875.06971, 816.569606, 757.3774, 697.987635, 638.894853, 581.802574, 525.512772, 470.036397, 415.384401, 361.925557, 309.169865, 256.985149, 205.239231, 153.535543, 102.112055, 50.9423473, 0., -49.5539819, -99.4084146, -150.064689, -202.024195, -253.775575, -308.638068, -367.918164, -432.922351, -499.071792, -575.912436, -667.104903, -776.309814, -885.726459, -1029.06132, -1218.55955, -1466.46631, -1718.56271, -2080.14355, -2590.0396, -3287.08162, -4146.20453, -5397.92663, -7051.77975, -9541.19154, -12538.7207, -18336.5059, -28571.7906, -44881.8187}, 
{2185.86243, 2129.66864, 2072.54176, 2014.18321, 1954.2944, 1892.57674, 1828.73164, 1762.46052, 1693.46479, 1613.81027, 1533.88821, 1456.45424, 1384.26402, 1331.09789, 1284.27693, 1242.1469, 1203.05357, 1164.97587, 1126.77315, 1086.93792, 1043.96269, 989.074843, 930.938071, 870.950935, 810.511996, 754.725473, 699.802008, 645.657899, 592.209444, 538.934767, 486.363612, 434.587547, 383.69814, 334.548323, 286.163756, 238.331466, 190.838476, 142.966997, 95.210792, 47.5588098, 0., -46.4100171, -93.1756315, -140.734562, -189.524529, -237.936684, -289.27394, -344.79264, -405.74913, -467.474595, -539.520605, -625.513568, -729.079894, -833.298458, -970.562217, -1152.71659, -1391.60701, -1630.84359, -1977.80117, -2471.61931, -3151.43754, -3991.40953, -5225.63246, -6875.89452, -9345.32147, -12339.851, -17972.8354, -27755.7736, -43200.1643}, 
{1981.4486, 1917.74747, 1860.59169, 1808.1705, 1758.6731, 1710.28871, 1661.20656, 1609.61586, 1553.70582, 1477.58662, 1399.15813, 1322.2412, 1250.65666, 1201.62972, 1160.21509, 1124.87184, 1094.05904, 1069.0927, 1044.43221, 1017.39389, 985.294068, 932.67365, 874.734515, 813.90314, 752.606002, 700.319407, 649.600069, 600.054533, 551.289343, 501.502803, 452.272994, 403.769755, 356.162928, 310.714537, 266.065364, 221.948375, 178.096535, 133.578698, 89.0575891, 44.5318191, 0., -43.5656426, -87.5295569, -132.282577, -178.215537, -223.723063, -271.99068, -324.207705, -381.563456, -439.300701, -506.933925, -588.031063, -686.16005, -784.990737, -915.94838, -1090.56015, -1320.35322, -1546.9207, -1879.69744, -2358.18424, -3021.88189, -3844.23415, -5062.91291, -6711.52012, -9164.32682, -12162.1323, -17638.9538, -26980.504, -41572.4958}, 
{1794.4152, 1724.40808, 1667.3188, 1620.06974, 1579.58327, 1542.78174, 1506.58753, 1467.923, 1423.71051, 1351.59344, 1275.48475, 1200.01839, 1129.82833, 1084.63602, 1047.95293, 1018.37801, 994.510209, 980.512175, 967.19371, 950.928297, 928.089424, 877.880279, 820.71276, 759.828472, 698.469016, 649.593536, 602.839079, 557.560233, 513.111583, 466.653192, 420.611983, 375.220355, 330.710706, 288.654201, 247.408964, 206.671883, 166.13985, 124.746758, 83.2576926, 41.6747436, 0., -40.8980775, -82.2303199, -124.341187, -167.575137, -210.340302, -255.692002, -304.749227, -358.630969, -412.54525, -475.886418, -552.13785, -644.782925, -738.021909, -862.334537, -1028.91743, -1248.96722, -1462.16565, -1779.83017, -2241.76334, -2887.76773, -3690.53634, -4891.20049, -6534.3606, -8967.55178, -11963.4481, -17291.233, -26220.9625, -40022.6926}, 
{1628.55148, 1554.65838, 1498.11445, 1454.99224, 1421.36432, 1393.30328, 1366.88168, 1338.17208, 1303.24707, 1235.66334, 1163.01568, 1090.38301, 1022.84422, 981.368108, 948.787774, 923.826191, 905.206329, 898.844712, 893.393335, 884.697746, 868.603494, 821.350866, 766.232776, 706.936878, 647.150825, 601.764906, 558.783084, 517.41196, 476.858131, 433.670357, 390.776215, 348.44544, 306.947768, 268.010996, 229.863574, 192.192013, 154.682824, 116.244139, 77.6521993, 38.9048655, 0., -38.3123138, -77.0873241, -116.628058, -157.237541, -197.364371, -239.907775, -285.912554, -336.423505, -386.690844, -445.871786, -517.328963, -604.425008, -691.796217, -809.022092, -966.955801, -1176.45051, -1375.54392, -1677.03085, -2120.89065, -2747.10266, -3527.644, -4706.50075, -6337.96158, -8745.4887, -11730.9564, -16913.4172, -25457.6291, -38528.3499}, 
{1487.64668, 1413.5063, 1358.36997, 1318.04925, 1288.35567, 1265.10079, 1244.09613, 1221.15325, 1192.0837, 1129.6289, 1061.89856, 993.932245, 930.769552, 893.077305, 864.016948, 842.37717, 826.946657, 823.700716, 821.366769, 815.858856, 803.091017, 759.741547, 708.65453, 653.438306, 597.701215, 556.050561, 516.696131, 478.846677, 441.710953, 401.838723, 362.161324, 322.951103, 284.480408, 248.428604, 213.098217, 178.198789, 143.439861, 107.843807, 72.0822061, 36.1394671, 0., -35.7133438, -71.9099732, -108.860855, -146.836957, -184.371235, -224.167871, -267.193037, -314.412902, -361.220084, -416.38373, -483.099433, -564.562788, -645.717902, -755.312447, -903.842611, -1101.80458, -1286.02099, -1570.131, -1994.10022, -2597.89427, -3352.88507, -4504.81926, -6115.86863, -8488.62993, -11451.815, -16489.2504, -24670.9838, -37067.063}, 
{1375.49005, 1305.95973, 1253.47674, 1214.35207, 1184.89673, 1161.42172, 1140.23803, 1117.65665, 1089.9886, 1033.32271, 972.280995, 911.263313, 854.669523, 821.014935, 794.937775, 775.19172, 760.530448, 754.690595, 749.449696, 741.568244, 727.806733, 689.708458, 645.337991, 597.542705, 549.16997, 511.667545, 475.842264, 441.101348, 406.852016, 370.442713, 334.162944, 298.243438, 262.914923, 229.550707, 196.781911, 164.382231, 132.125366, 99.3187258, 66.3888102, 33.2958309, 0., -33.0061597, -66.5076708, -100.757245, -136.007596, -170.936864, -208.00216, -248.086024, -292.070995, -335.615571, -386.915951, -448.944294, -524.672755, -599.191206, -700.507005, -838.745224, -1024.03094, -1192.56234, -1457.96212, -1859.92612, -2438.15015, -3163.58747, -4282.16159, -5861.62737, -8187.4678, -11113.1819, -16002.4765, -23841.5066, -35616.4271}, 
{1295.87085, 1237.02659, 1188.82609, 1149.01198, 1115.32691, 1085.51353, 1057.31449, 1028.47242, 996.729976, 946.577353, 894.310623, 842.97341, 795.609339, 766.432322, 742.847581, 723.430622, 706.756955, 691.424756, 675.977799, 658.982526, 639.005378, 607.907736, 573.643128, 537.46002, 500.606877, 467.832902, 435.48553, 403.412934, 371.463285, 338.766752, 306.176712, 273.828538, 241.857605, 211.020988, 180.583677, 150.432364, 120.453741, 90.4418588, 60.4131083, 30.291239, 0., -30.0957535, -60.6898206, -92.0348943, -124.383668, -156.637225, -190.940513, -228.086866, -268.869621, -309.359907, -356.962153, -414.358576, -484.231398, -551.62037, -643.907167, -770.830997, -942.131067, -1094.13344, -1339.35572, -1716.90239, -2265.87789, -2957.07913, -4034.53329, -5568.78338, -7832.49467, -10702.2149, -15436.8397, -22949.6777, -34154.0375}, 
{1204.65587, 1157.00006, 1114.47906, 1076.03051, 1040.59206, 1007.10134, 974.496, 941.713669, 907.69199, 863.90708, 819.742711, 777.121132, 737.964591, 712.778866, 691.469265, 672.524625, 654.433783, 631.935693, 608.769028, 584.922579, 560.385137, 535.126853, 509.162613, 482.488661, 455.101244, 426.520612, 397.409404, 367.954264, 338.341835, 309.081426, 279.907949, 250.878982, 222.052103, 193.704483, 165.586269, 137.667202, 109.917023, 82.4541578, 55.0401876, 27.5853794, 0., -27.3906496, -55.2573507, -83.8558505, -113.441896, -143.126009, -174.767252, -209.079461, -246.776473, -284.355627, -328.433857, -381.411599, -445.689291, -506.351102, -590.040242, -706.083651, -863.808272, -999.956678, -1225.47393, -1578.72072, -2098.05773, -2753.3297, -3788.44518, -5278.67022, -7480.52598, -10279.4618, -14863.9593, -22096.742, -32840.5336}, 
{1074.63988, 1033.0459, 996.321555, 963.439721, 933.373293, 905.095166, 877.57823, 849.795378, 820.719503, 782.391091, 743.488403, 705.757293, 670.943616, 648.479622, 629.350211, 612.226682, 595.780332, 575.305119, 554.200613, 532.489048, 510.192655, 487.464298, 464.143325, 440.199715, 415.603449, 389.833046, 363.546529, 336.910461, 310.091403, 283.54333, 257.03043, 230.6043, 204.316538, 178.421378, 152.686729, 127.083133, 101.581135, 76.2399296, 50.9059493, 25.5142782, 0., -25.256491, -50.9435458, -77.3042049, -104.581509, -131.993253, -161.217822, -192.908352, -227.717982, -262.453961, -303.153673, -352.008612, -411.210273, -466.199332, -542.618428, -649.359384, -795.314019, -917.313637, -1125.13479, -1456.4935, -1949.1058, -2570.59931, -3568.95528, -5028.68341, -7183.5431, -9910.74615, -14364.4451, -21379.5571, -31790.9996}, 
{952.045299, 916.813429, 886.000482, 858.645986, 833.789466, 810.47045, 787.728463, 764.603033, 740.133686, 707.0761, 673.26719, 640.260019, 609.607655, 589.582905, 572.331196, 556.717696, 541.607574, 522.902498, 503.616539, 483.800265, 463.504248, 443.03075, 422.07797, 400.595797, 378.534124, 355.357107, 331.694667, 307.69099, 283.49026, 259.485015, 235.471748, 211.495303, 187.600525, 164.018334, 140.533068, 117.115138, 93.7349577, 70.3914253, 47.0150752, 23.564927, 0., -23.2517332, -46.8877878, -71.1367257, -96.2271088, -121.468367, -148.375847, -177.545764, -209.574332, -241.585969, -279.037406, -323.913578, -378.199419, -427.713909, -497.074318, -594.73196, -729.138153, -837.32558, -1027.73164, -1337.3751, -1803.27474, -2390.79442, -3351.91758, -4779.57717, -6886.99293, -9541.19849, -13873.6851, -20699.8296, -30835.0089}, 
{837.812697, 808.936927, 783.919108, 761.884766, 741.959432, 723.268633, 704.937899, 686.092758, 665.858739, 637.79203, 608.815237, 580.281624, 553.544456, 535.682276, 520.03296, 505.659659, 491.625529, 474.469091, 456.787981, 438.655203, 420.143763, 401.668489, 382.823831, 363.546065, 343.771465, 322.977172, 301.742249, 280.186624, 258.430227, 236.798139, 215.123074, 193.4429, 171.795483, 150.388204, 129.021612, 107.665771, 86.2907419, 64.8428209, 43.325344, 21.7178806, 0., -21.3729548, -43.0874978, -65.3543693, -88.3843094, -111.556974, -136.246622, -162.996426, -192.349561, -221.745105, -256.071965, -297.114955, -346.658886, -390.927315, -453.490813, -542.358694, -665.540274, -760.3361, -933.747761, -1222.06808, -1661.58988, -2215.46107, -3139.40919, -4533.64928, -6593.803, -9176.58984, -13398.9314, -20062.6987, -29969.7628}, 
{732.882625, 710.050697, 690.4807, 673.391519, 658.002041, 643.53115, 629.197732, 614.220674, 597.81886, 574.368803, 549.868711, 525.47442, 502.341765, 486.371293, 472.076243, 458.714566, 445.544212, 429.746158, 413.486118, 396.852833, 379.935045, 363.21979, 346.238194, 328.919679, 311.19367, 292.577619, 273.577705, 254.288137, 234.803126, 215.37436, 195.875577, 176.337996, 156.792834, 137.423841, 118.04869, 98.6375877, 79.1607395, 59.5282925, 39.7945341, 19.9536936, 0., -19.6167344, -39.5400971, -59.9580924, -81.0587245, -102.2647, -124.83544, -149.265068, -176.047706, -202.924822, -234.244259, -271.6012, -316.590832, -355.872032, -411.950817, -492.396897, -604.779981, -686.688787, -843.66641, -1111.27496, -1525.07654, -2046.14531, -2933.5072, -4293.19749, -6306.90082, -8822.69125, -12947.4358, -19473.3033, -29192.4625}, 
{638.19564, 620.789034, 606.088527, 593.401702, 582.036144, 571.299435, 560.499161, 548.942905, 535.938251, 516.636341, 496.163778, 475.490721, 455.58733, 441.243514, 428.081786, 415.544409, 403.073644, 388.474958, 373.482125, 358.192127, 342.701942, 327.526931, 312.178341, 296.585799, 280.678933, 264.042825, 247.089466, 229.886301, 212.500777, 195.105336, 177.620429, 160.071501, 142.483999, 125.018098, 107.510627, 89.933143, 72.2572021, 54.3820161, 36.3804242, 18.2529207, 0., -17.9796505, -36.2430068, -54.9488516, -74.2559679, -93.5971696, -114.147596, -136.356417, -160.672804, -185.118576, -213.541195, -247.360771, -287.997414, -322.580542, -372.537235, -445.003884, -547.116876, -616.727231, -757.970853, -1005.69828, -1394.76004, -1884.39319, -2736.28873, -4060.51959, -6029.21391, -8485.27376, -12526.4503, -18936.7826, -28500.3094}, 
{554.692301, 541.786237, 531.145858, 522.150773, 514.180592, 506.614925, 498.833382, 490.215574, 480.14111, 464.424568, 447.436604, 429.98284, 412.8689, 399.892497, 387.670329, 375.811181, 363.923841, 350.396751, 336.547181, 322.472055, 308.268298, 294.432189, 280.501557, 266.413583, 252.105452, 237.257168, 222.165962, 206.871889, 191.415003, 175.882728, 160.248799, 144.53432, 128.760398, 113.063829, 97.3037508, 81.4549918, 65.4923811, 49.3381679, 33.0407928, 16.5961166, 0., -16.4582818, -33.193648, -50.3276035, -67.9816533, -85.560007, -104.188383, -124.275204, -146.228893, -168.319821, -193.949684, -224.382126, -260.880791, -291.085326, -335.332971, -400.336967, -492.810558, -550.795024, -677.144355, -906.040578, -1271.66572, -1731.75074, -2549.83088, -3837.91334, -5763.66979, -8170.1084, -12143.2269, -18458.2756, -27890.5049}, 
{483.313163, 473.676603, 466.055961, 459.874189, 454.554236, 449.519054, 444.191593, 437.994804, 430.351638, 417.563406, 403.423355, 388.60309, 373.77422, 361.911802, 350.462611, 339.176876, 327.804822, 315.252799, 302.452463, 289.491591, 276.457961, 263.777843, 251.065126, 238.272192, 225.351423, 212.105025, 198.695625, 185.135674, 171.437623, 157.598192, 143.651858, 129.617363, 115.513453, 101.453888, 87.3243877, 73.1056889, 58.7785281, 44.3309239, 29.7334185, 14.9638358, 0., -15.0492067, -30.3894418, -46.0953047, -62.2413945, -78.1588373, -94.9630951, -113.026156, -132.72001, -152.52201, -175.456634, -202.653723, -235.243119, -261.418867, -300.420927, -358.55346, -442.120629, -489.235757, -601.67018, -813.004394, -1156.81889, -1589.764, -2376.21074, -3627.6765, -5513.19596, -7882.96623, -11805.0174, -18042.9214, -27360.2502}, 
{430.484418, 421.113456, 414.158086, 408.98146, 404.94673, 401.417046, 397.755561, 393.325426, 387.489794, 376.483753, 364.049742, 350.802139, 337.355319, 326.256844, 315.414632, 304.669786, 293.863408, 282.296992, 270.567093, 258.73066, 246.844638, 235.147337, 223.441797, 211.712417, 199.9436, 188.128389, 176.239085, 164.256634, 152.161977, 139.83341, 127.395587, 114.870511, 102.280187, 89.7188522, 77.1073842, 64.4388941, 51.7064931, 39.0573716, 26.2689294, 13.2726458, 0., -13.7764105, -27.9026719, -42.3837523, -57.2246196, -71.5863824, -86.6554121, -102.77422, -120.285319, -137.798375, -158.081883, -182.171493, -211.102857, -233.791427, -268.24113, -320.335691, -395.958842, -432.859693, -532.310436, -727.448645, -1051.41189, -1460.81344, -2218.36381, -3429.83, -5277.77804, -7636.36413, -11536.5847, -17720.2261, -26929.0746}, 
{392.308466, 381.186439, 373.266059, 367.790599, 364.003329, 361.147518, 358.466437, 355.203357, 350.601548, 340.48549, 328.88476, 316.410145, 303.672432, 293.084004, 282.733413, 272.510809, 262.30634, 251.700434, 241.016849, 230.269621, 219.472788, 208.611231, 197.739805, 186.884208, 176.07014, 165.51361, 154.973883, 144.400534, 133.743138, 122.753431, 111.657963, 100.485445, 89.2645885, 78.0710907, 66.8678814, 55.664877, 44.4719937, 33.6607917, 22.7348862, 11.559536, 0., -12.6123, -25.6647331, -39.0785057, -52.774824, -65.6832532, -79.1132972, -93.3828187, -108.809681, -124.072008, -141.783297, -162.917303, -188.447785, -208.078692, -238.561512, -285.377927, -354.009618, -381.598833, -469.302458, -649.937946, -956.32275, -1345.70917, -2077.61012, -3248.67972, -5062.9382, -7430.12351, -11329.7366, -17473.8377, -26574.487}, 
{360.66773, 347.89344, 338.935119, 332.947164, 329.083975, 326.499951, 324.349491, 321.786994, 317.966858, 308.406182, 297.351585, 285.41239, 273.197914, 263.040569, 253.137348, 243.408333, 233.77361, 224.009152, 214.236795, 204.434269, 194.579301, 184.439194, 174.28627, 164.182425, 154.189557, 144.710428, 135.329724, 125.972995, 116.565793, 106.761765, 96.8671266, 86.916191, 76.9432708, 67.0016617, 57.0990999, 47.2623045, 37.5179949, 28.466852, 19.3320488, 9.91071997, 0., -11.5081985, -23.54943, -35.9644712, -48.5940986, -60.1461305, -72.0474852, -84.5921225, -98.0740026, -111.204635, -126.493411, -144.867268, -167.253148, -184.019583, -210.875282, -252.97055, -315.455689, -335.026177, -412.669073, -580.916609, -872.30102, -1244.6194, -1954.6094, -3090.28309, -4876.46181, -7258.87633, -11162.8104, -17268.6597, -26256.8197}, 
{332.848826, 319.339406, 309.84285, 303.496527, 299.437804, 296.804049, 294.732628, 292.360909, 288.82626, 279.6239, 268.990204, 257.519397, 245.805708, 236.102767, 226.681634, 217.472775, 208.406655, 199.343792, 190.312576, 181.27145, 172.178856, 162.676833, 153.16679, 143.733731, 134.462662, 125.865105, 117.42894, 109.068565, 100.698377, 91.9270899, 83.0970604, 74.2449607, 65.4074631, 56.6069349, 47.9000748, 39.3292765, 30.936934, 23.542257, 16.1000972, 8.34212235, 0., -10.4410908, -21.5043985, -32.9596599, -44.576612, -54.88014, -65.3827734, -76.3521899, -88.0560672, -99.2001149, -112.238766, -128.064485, -147.569738, -161.556655, -185.044169, -222.960877, -280.235378, -293.441953, -363.205243, -521.795574, -801.483269, -1160.10306, -1853.23206, -2962.50825, -4729.38986, -7132.2067, -11041.982, -17106.5658, -25973.8078}, 
{304.952152, 292.952218, 284.33093, 278.349426, 274.268842, 271.350315, 268.854983, 266.043982, 262.17845, 253.319746, 243.208695, 232.386348, 221.393753, 212.329699, 203.554397, 194.985801, 186.541861, 178.001582, 169.477443, 160.942978, 152.371717, 143.465933, 134.578921, 125.792715, 117.189349, 109.241916, 101.484971, 93.8441289, 86.2450012, 78.3468177, 70.4481283, 62.5810991, 54.7778966, 47.0124833, 39.3985105, 31.9914262, 24.8466778, 18.9719631, 13.08958, 6.87407634, 0., -9.38029515, -19.4611898, -29.9592584, -40.5910756, -49.7813155, -59.0552134, -68.6461041, -78.7873225, -88.1276716, -99.1188308, -112.627947, -129.522169, -140.631068, -160.874397, -195.134334, -248.293057, -257.393244, -322.292372, -475.008418, -747.559361, -1096.57811, -1780.23784, -2878.9862, -4641.40046, -7069.0396, -10983.4062, -17001.0171, -25738.3888}, 
{279.78566, 269.98824, 262.640456, 257.177906, 253.036189, 249.650905, 246.457651, 242.892026, 238.389629, 229.788609, 220.160994, 209.981364, 199.724297, 191.311845, 183.192125, 175.260728, 167.413243, 159.281274, 151.129992, 142.960581, 134.774227, 126.404438, 118.087145, 109.890602, 101.883063, 94.4378202, 87.1960756, 80.1040686, 73.1080386, 65.9558819, 58.8715185, 51.8805254, 45.0084794, 38.1787473, 31.56, 25.2186982, 19.2213027, 14.7377201, 10.2895871, 5.50198627, 0., -8.3384809, -17.4423055, -26.9875147, -36.6501492, -44.7971021, -52.9372211, -61.2702062, -69.9957574, -77.6869199, -86.8207102, -98.2474902, -112.817622, -121.191888, -138.48606, -169.626333, -219.538899, -226.031281, -287.995812, -437.206154, -705.435969, -1047.95037, -1726.04867, -2820.76157, -4583.33267, -7035.48163, -10946.811, -16899.9525, -25477.5381}, 
{257.302813, 250.074573, 244.204199, 239.330028, 235.090398, 231.123647, 227.068111, 222.562129, 217.244038, 208.851392, 199.683625, 190.139389, 180.617334, 172.839114, 165.351176, 158.02297, 150.723945, 142.906536, 135.024016, 127.112641, 119.208668, 111.321731, 103.525358, 95.8664567, 88.3919317, 81.3337627, 74.4797527, 67.8027781, 61.2757154, 54.7583607, 48.3819031, 42.1644511, 36.1241133, 30.1364421, 24.4191248, 19.0472927, 14.0960769, 10.8651988, 7.71536327, 4.23186531, 0., -7.30986874, -15.4315422, -24.013753, -32.7052334, -39.8530662, -46.9282936, -54.100308, -61.5385019, -67.7347666, -75.206996, -84.795583, -97.3409206, -103.184816, -117.865682, -146.423345, -193.897633, -199.101319, -259.792105, -407.500641, -673.757574, -1012.59923, -1688.03923, -2782.74657, -4546.80967, -7020.11695, -10915.9374, -16779.1949, -25154.8133}, 
{237.457073, 232.838319, 228.454933, 224.153853, 219.782019, 215.186369, 210.213843, 204.71138, 198.525917, 190.328997, 181.613115, 172.69537, 163.892858, 156.701414, 149.787905, 142.997937, 136.177112, 128.601035, 120.91331, 113.187537, 105.497321, 98.0471923, 90.727455, 83.5593415, 76.5640841, 69.8106881, 63.253503, 56.8946518, 50.7362571, 44.7583326, 38.993954, 33.4540878, 28.1497004, 22.9162824, 18.0104665, 13.5134097, 9.5062687, 7.38006997, 5.38215293, 3.0697266, 0., -6.28867938, -13.4126964, -21.0072975, -28.7077289, -34.8747742, -40.9279278, -47.0122213, -53.2726862, -58.1281187, -64.1402798, -72.1446952, -82.9768906, -86.555554, -98.9997838, -125.511841, -171.293985, -176.348612, -237.157795, -385.00374, -651.168655, -988.904083, -1663.58422, -2759.85339, -4523.45463, -7011.52972, -10874.5267, -16614.567, -24733.7721}, 
{221.112345, 217.273009, 213.245771, 208.971399, 204.390662, 199.444326, 194.07316, 188.21793, 181.819406, 174.07591, 165.967632, 157.732316, 149.607708, 142.759168, 136.12578, 129.574241, 122.971252, 115.636065, 108.201803, 100.754142, 93.3787588, 86.3649399, 79.5133073, 72.8280932, 66.31353, 59.9400399, 53.7591891, 47.7887339, 42.0464305, 36.5675385, 31.3453094, 26.3904981, 21.7138595, 17.1544906, 12.9634673, 9.22020788, 6.00413027, 4.68026363, 3.52817058, 2.11302482, 0., -5.20156086, -11.1883212, -17.612775, -24.1274164, -29.1515342, -34.0641099, -39.0109194, -44.1377391, -47.8775291, -52.7740078, -59.6580777, -69.3606412, -71.3464121, -82.3589571, -107.775654, -152.973881, -158.692157, -220.802263, -370.537118, -639.129644, -981.466477, -1657.81939, -2751.51569, -4503.96735, -6999.36377, -10792.7001, -16321.2699, -24022.367}, 
{207.178203, 203.732178, 199.820678, 195.45373, 190.64136, 185.393593, 179.720456, 173.631975, 167.138176, 159.874527, 152.375435, 144.800749, 137.310317, 130.721435, 124.273526, 117.863461, 111.388109, 104.253972, 97.0444365, 89.8525206, 82.7712422, 76.1562566, 69.7328892, 63.4891028, 57.4128605, 51.3843422, 45.5424068, 39.9181301, 34.5425882, 29.4933277, 24.7363656, 20.2841896, 16.1492874, 12.1875237, 8.63065835, 5.55382834, 3.03217056, 2.36462843, 1.91300975, 1.26292884, 0., -4.1202108, -8.95002157, -14.1717987, -19.4679084, -23.3384285, -27.1209295, -30.9706931, -35.0430011, -37.8103107, -41.7838582, -47.7920555, -56.6633142, -57.3923226, -67.3747053, -92.1723635, -137.347198, -144.531429, -208.788512, -361.25222, -633.056326, -981.531677, -1659.22083, -2747.46148, -4481.02972, -6970.9883, -10669.8149, -15944.6335, -23162.5681}, 
{195.375656, 192.004692, 188.015746, 183.466195, 178.413415, 172.914785, 167.027681, 160.809479, 154.317557, 147.537176, 140.626675, 133.672277, 126.760207, 120.360119, 114.021432, 107.676999, 101.25967, 94.2939214, 87.28433, 80.3270978, 73.5184264, 67.2620822, 61.2236764, 55.3763849, 49.6933835, 43.9759018, 38.4378405, 33.1211539, 28.0677965, 23.3846307, 19.0227396, 14.9981141, 11.3267455, 7.89149902, 4.89474193, 2.40571575, 0.493661944, 0.363790214, 0.494986523, 0.502105057, 0., -3.05610252, -6.72786903, -10.7365941, -14.8035724, -17.5234724, -20.1948652, -22.9896958, -26.0799092, -27.9994607, -31.213481, -36.549111, -44.8334916, -44.6143205, -53.9099591, -78.4593255, -124.001338, -133.312004, -200.278317, -355.824359, -630.8742119999999, -985.890452, -1663.18168, -2741.44858, -4446.50068, -6914.71693, -10493.0898, -15475.9586, -22157.6622}, 
{185.425713, 181.879415, 177.667063, 172.87414, 167.58613, 161.888515, 155.86678, 149.606408, 143.192881, 136.876184, 130.511499, 124.11851, 117.716901, 111.447124, 105.159785, 98.8262595, 92.4179225, 85.5950808, 78.7646057, 72.0222999, 65.4639662, 59.5233565, 53.8231446, 48.3239537, 42.9864068, 37.5470253, 32.2801747, 27.2361189, 22.4651217, 18.090378, 14.0600481, 10.3952234, 7.11699527, 4.1425339, 1.63842043, -0.332685234, -1.7081232, -1.39162518, -0.767583015, -0.186780138, 0., -2.0207093, -4.55193494, -7.35938705, -10.2087758, -11.7946811, -13.3823954, -15.166081, -17.3399, -18.5179766, -21.1065264, -25.9317269, -33.8197554, -32.9334411, -41.8276492, -66.3938965, -112.5237, -124.479458, -194.433453, -352.92885, -630.508812, -991.333571, -1665.09509, -2727.23475, -4392.23916, -6818.86329, -10249.7435, -14906.5457, -21010.9359}, 
{177.049383, 173.145211, 168.610719, 163.542914, 158.038804, 152.195397, 146.109701, 139.878725, 133.599477, 127.703878, 121.820057, 115.911057, 109.93992, 103.754355, 97.4788721, 91.1226459, 84.6948522, 77.9966166, 71.3283849, 64.7825529, 58.4515163, 52.7810191, 47.3687694, 42.1658237, 37.1232385, 31.9300195, 26.9040941, 22.1013388, 17.5776303, 13.4595, 9.70390802, 6.33846907, 3.39079796, 0.816745601, -1.25560381, -2.76992995, -3.66991249, -2.97099193, -1.91638277, -0.821060353, 0., -1.02550445, -2.45229068, -4.09240305, -5.75788596, -6.24006989, -6.77999864, -7.59800213, -8.91441027, -9.4388557, -11.5066445, -15.9423856, -23.5706878, -22.2707196, -30.9907063, -55.7334329, -102.501684, -117.479367, -190.415696, -351.241008, -629.885639, -994.651801, -1660.35421, -2698.57779, -4310.10409, -6671.741, -9926.99437, -14227.6955, -19725.6759}, 
{187.82453, 175.380638, 164.522254, 155.048751, 146.759502, 139.453881, 132.931261, 126.991014, 121.432515, 115.682444, 110.061944, 104.519465, 99.0034562, 93.4094581, 87.7599952, 82.0246817, 76.1731321, 69.8497126, 63.4793845, 57.1618613, 50.9968561, 45.3113449, 39.8868738, 34.732251, 29.8562852, 25.3107344, 21.0442775, 17.0485431, 13.3151599, 9.78749463, 6.52474233, 3.53783623, 0.837709548, -1.92123783, -4.22892593, -5.93180815, -6.87633789, -5.53361613, -3.67558966, -1.69885283, 0., -0.623134421, -1.65833802, -2.8432016, -3.91531598, -3.23888229, -2.47423689, -1.90832646, -1.82809767, -0.571915493, -1.15474102, -4.64295363, -12.1029327, -13.5154792, -25.4667823, -53.4574527, -102.988101, -117.986443, -190.155143, -349.623968, -626.522688, -999.557189, -1653.12889, -2646.10523, -4169.23759, -6411.62763, -9437.09199, -13340.7686, -18217.7952}, 
{202.594033, 179.751196, 161.257808, 146.529047, 134.980092, 126.026122, 119.082316, 113.563853, 108.885911, 103.059218, 97.4651852, 92.0807714, 86.8829361, 82.0677256, 77.3053773, 72.4852159, 67.4965658, 61.6813602, 55.6952715, 49.6465805, 43.6435684, 37.8332563, 32.2696889, 27.0456511, 22.2539279, 18.4601648, 15.0951418, 12.0624993, 9.26587801, 6.39380025, 3.65107215, 1.02738158, -1.48758357, -4.64361367, -7.41575128, -9.51851721, -10.6664323, -8.53514926, -5.69360421, -2.67186512, 0., -0.508669667, -1.50711222, -2.60515855, -3.41263955, -1.63513838, 0.451567238, 2.47594729, 4.06647178, 7.34644229, 8.45156461, 6.01237612, -1.34058582, -6.53591649, -22.7602928, -54.7595243, -107.279421, -121.97668, -191.921867, -347.096438, -617.481845, -993.994698, -1623.81099, -2550.28408, -3950.79094, -6012.26715, -8749.39028, -12219.9024, -16481.5456}, 
{202.282333, 175.871048, 154.829792, 138.423698, 125.9179, 116.577528, 109.667717, 104.453599, 100.200306, 94.3003516, 88.640535, 83.2350361, 78.0980346, 73.6586754, 69.3501876, 65.0207652, 60.5186022, 55.0509575, 49.3633342, 43.5603002, 37.7464236, 31.9534243, 26.3878575, 21.1834305, 16.4738504, 13.064521, 10.1487743, 7.59163891, 5.25814349, 2.701185, 0.222776282, -2.18720142, -4.53886688, -7.71996386, -10.5119361, -12.5738524, -13.5647814, -10.8155658, -7.24479091, -3.44281592, 0., -0.08641868, -0.680828277, -1.34170145, -1.62751084, 1.08043447, 4.17563281, 7.22874608, 9.81043621, 14.0900627, 16.0001108, 14.0717634, 6.83620339, -0.0235436439, -17.9888746, -51.3893436, -104.554505, -121.858854, -191.569027, -341.996602, -601.453156, -968.32838, -1560.69951, -2412.60149, -3678.80089, -5530.33367, -7952.42665, -10976.7589, -14635.0096}, 
{191.434042, 166.516391, 146.685199, 131.231755, 119.447344, 110.623254, 104.050772, 99.021184, 94.8257772, 88.8485881, 83.051054, 77.4873619, 72.211699, 67.8168808, 63.6030143, 59.4088351, 55.0730788, 49.8012396, 44.318591, 38.7171651, 33.088994, 27.4028755, 21.9233699, 16.7918031, 12.1495011, 8.82957003, 6.00484342, 3.539935, 1.29945857, -1.20513783, -3.61480841, -5.9236731, -8.12585181, -11.0434715, -13.5114422, -15.1926812, -15.7501056, -12.5223676, -8.42635531, -4.05469173, 0., 0.584774086, 0.691459186, 0.751561269, 1.1965863, 4.6894944, 8.53775573, 12.2802946, 15.4560354, 19.9833366, 22.0699146, 20.30292, 13.2695034, 6.39295062, -11.3101771, -43.9871831, -95.785371, -117.497179, -187.966722, -332.683249, -577.13601, -923.480797, -1467.20723, -2237.65049, -3360.76265, -4981.95802, -7072.18514, -9646.00898, -12717.9945}, 
{174.593769, 154.463423, 138.271025, 125.452267, 115.442846, 107.678454, 101.594787, 96.6275385, 92.2124029, 86.1466711, 80.1598024, 74.3428524, 68.7868768, 64.1769151, 59.7724458, 55.426931, 50.993833, 45.7749414, 40.3960599, 34.9313197, 29.4548519, 23.9126366, 18.5582165, 13.5169832, 8.91432833, 5.46107933, 2.46301766, -0.188639368, -2.60267448, -5.23995495, -7.71634509, -9.99979375, -12.0582498, -14.5073195, -16.4082314, -17.4698714, -17.4011254, -13.8030567, -9.33550295, -4.55047905, 0., 1.44606418, 2.48069555, 3.47902116, 4.81616805, 8.97369958, 13.3777319, 17.5608176, 21.055509, 25.3506551, 27.2399934, 25.4735579, 18.8013827, 13.0848776, -2.88184972, -33.1933156, -81.9440363, -108.755869, -179.985053, -317.515168, -543.229795, -860.374511, -1346.74691, -2030.02411, -3004.17146, -4383.27103, -6134.6498, -8262.3235, -10770.3079}, 
{156.306122, 142.488343, 131.034262, 121.584286, 113.778824, 107.258283, 101.66307, 96.6335931, 91.81026, 85.6373442, 79.4298407, 73.3066107, 67.3865153, 62.3733517, 57.5670705, 52.8525584, 48.1147021, 42.814798, 37.4307592, 32.0169088, 26.6275694, 21.2137341, 15.9743875, 11.0051847, 6.40178052, 2.66481611, -0.677034553, -3.69011116, -6.44075338, -9.31805307, -11.9364969, -14.2333237, -16.1457723, -18.0046911, -19.1962656, -19.5002906, -18.6965613, -14.805135, -10.0694394, -4.97316436, 0., 2.43860714, 4.5578262, 6.6450698, 8.98775053, 13.7147082, 18.5353573, 23.0005397, 26.6610969, 30.516409, 32.0893643, 30.3513891, 24.2739098, 20.4235486, 7.13845799, -19.6480136, -64.0025181, -95.4991373, -166.494121, -294.851148, -498.433899, -779.932084, -1202.73129, -1794.31537, -2616.52254, -3750.40353, -5165.80466, -6860.37338, -8831.75711}, 
{141.115714, 133.367348, 126.421904, 120.126861, 114.329698, 108.877894, 103.618928, 98.4002787, 93.0694261, 86.7633508, 80.3242294, 73.8837403, 67.5735619, 62.0407637, 56.6954766, 51.4632226, 46.2695237, 40.7635444, 35.2577072, 29.788077, 24.3907188, 19.0371946, 13.8538735, 8.90262193, 4.24530605, 0.146547664, -3.61564479, -7.06050731, -10.2072759, -13.3542189, -16.129927, -18.4420232, -20.1981307, -21.4287693, -21.8695065, -21.3788065, -19.8151337, -15.6761045, -10.7253702, -5.36573418, 0., 3.50355851, 6.79379651, 10.0540987, 13.4678499, 18.6941784, 23.8504278, 28.5296855, 32.3250388, 35.8049895, 37.1970445, 35.7041252, 30.5291532, 28.7802748, 18.5930966, -3.99155006, -42.9328337, -77.5911981, -146.364027, -263.049978, -441.44771, -683.076077, -1038.57315, -1535.11731, -2205.31111, -3099.48633, -4191.63378, -5474.82953, -6942.14962}, 
{103.940274, 109.246159, 112.100836, 112.786499, 111.585343, 108.779563, 104.651353, 99.482908, 93.556423, 87.4462621, 81.0255828, 74.4597118, 67.9139761, 62.0289608, 56.3046312, 50.7162112, 45.2389241, 39.780798, 34.4111305, 29.1320235, 23.945579, 18.7904806, 13.7576159, 8.87445415, 4.16846479, -0.307552597, -4.56159118, -8.57631379, -12.3343833, -16.2187059, -19.6516038, -22.4556421, -24.4533865, -25.0813324, -24.702543, -23.294012, -20.8327326, -16.4901945, -11.3710965, -5.77463332, 0., 4.75478168, 9.44887446, 14.1426142, 18.8963367, 25.0263267, 30.8345917, 35.8790878, 39.7177713, 41.4854709, 41.3325214, 38.9861302, 34.1735046, 36.3321689, 31.5948865, 15.8047379, -15.1951963, -53.3876082, -119.973336, -223.978989, -374.43118, -567.576215, -850.781612, -1259.91447, -1784.82538, -2422.06506, -3187.40028, -4091.7938, -5146.20839}, 
{119.696459, 121.904076, 122.249633, 120.94698, 118.209966, 114.25244, 109.288254, 103.531255, 97.1952948, 90.7420274, 84.0383753, 77.1990662, 70.3388277, 63.8365474, 57.4371289, 51.1496357, 44.9831315, 38.9566103, 33.0652331, 27.3140912, 21.708276, 16.2226188, 10.9045753, 5.77134083, 0.84011102, -3.93325859, -8.44569664, -12.6554717, -16.5208523, -20.3073567, -23.5431039, -26.0634625, -27.7038011, -27.8850365, -27.02277, -25.1181511, -22.1723292, -17.5689495, -12.1736673, -6.23463379, 0., 5.72320321, 11.4692564, 17.2135602, 22.9315156, 29.3276159, 35.3565322, 40.7020283, 45.0478681, 48.0136362, 49.3729471, 48.8352361, 46.1099386, 48.6623724, 45.3437374, 32.7611163, 7.52159136, -26.3360594, -82.6091273, -167.663208, -287.863898, -441.843924, -659.167485, -966.431705, -1346.42536, -1792.66723, -2299.56444, -2861.30205, -3472.06512}, 
{169.150239, 157.676836, 147.513519, 138.460033, 130.316123, 122.881534, 115.956011, 109.3393, 102.831146, 95.6550193, 88.4174499, 81.1486922, 73.8790011, 66.6745961, 59.5153815, 52.4172262, 45.3959992, 38.4031129, 31.5446755, 24.8623387, 18.3977542, 12.2411791, 6.36621736, 0.795078651, -4.45002757, -9.41572308, -13.9834346, -18.10342, -21.7259373, -24.9162069, -27.4635394, -29.2722079, -30.2464854, -30.1545278, -29.0911722, -27.0151386, -23.8851471, -18.9590206, -13.1767348, -6.77736834, 0., 6.50069279, 13.0674696, 19.6274911, 26.1079182, 32.5190318, 38.6716244, 44.4596088, 49.7768977, 55.8502074, 60.7075256, 63.7096435, 64.2173524, 64.9421617, 60.5538573, 49.0729433, 28.5199236, 5.13411772, -34.5703119, -95.8603875, -184.003131, -305.004472, -461.914711, -660.672514, -900.983227, -1205.90322, -1511.31808, -1774.90607, -1954.34547}, 
{157.588269, 149.87012, 142.397846, 135.138514, 128.05919, 121.126941, 114.308835, 107.571936, 100.883312, 94.184065, 87.4776118, 80.7414055, 73.952899, 66.9294522, 59.8726476, 52.8239751, 45.8249243, 38.9794848, 32.2416468, 25.6279002, 19.1547349, 12.7006072, 6.47525399, 0.550378712, -5.00231531, -9.89727766, -14.3621909, -18.4108905, -22.0572119, -25.8017927, -28.9769454, -31.4017847, -32.8954252, -32.8713233, -31.7165151, -29.4123786, -25.9402917, -20.6885535, -14.4688524, -7.49979804, 0., 7.57072875, 15.3314634, 23.1600759, 30.9344379, 38.6372245, 45.9995829, 52.8574638, 59.0468175, 64.5920658, 69.0652995, 72.2270807, 73.8379713, 76.5085846, 76.0094109, 70.9609919, 59.9838691, 45.5410236, 20.8735817, -16.934891, -70.8008287, -145.218014, -238.370836, -359.648158, -484.977861, -600.901629, -684.966921, -716.969752, -676.706137}, 
{170.100056, 162.077358, 154.716959, 147.843675, 141.282325, 134.857726, 128.394696, 121.718051, 114.652611, 106.038997, 97.0798998, 87.993815, 78.9992377, 71.0318317, 63.3060567, 55.7535409, 48.3059124, 40.5721135, 32.9355325, 25.4568718, 18.1968339, 11.2698263, 4.66136436, -1.58933151, -7.44304092, -12.8267762, -17.748591, -22.1827721, -26.1036058, -29.8062716, -32.815806, -34.9781382, -36.1391978, -35.7891867, -34.2720525, -31.5760155, -27.6892961, -22.1624404, -15.5964124, -8.15450233, 0., 8.75252273, 17.8715705, 27.174366, 36.478132, 45.3117409, 53.896106, 62.1637903, 70.0473565, 77.8136669, 84.9272654, 91.1869954, 96.3917, 101.4355, 104.58385, 105.197482, 102.63713, 95.6793504, 84.50272, 68.7016409, 47.8705152, 17.3576345, -10.504268, -36.2030578, -39.5162004, -13.1099957, 81.1941944, 273.863774, 595.366146}, 
{177.272638, 176.171898, 172.67048, 167.127495, 159.902055, 151.353269, 141.840249, 131.722107, 121.357952, 111.932907, 102.649669, 93.5369434, 84.6234388, 76.1690857, 67.8788777, 59.6890322, 51.5357667, 42.9013283, 34.3574928, 26.0220658, 18.0128529, 10.8353733, 4.06463329, -2.336647, -8.40574756, -14.5091726, -20.223288, -25.4536843, -30.1059515, -34.3111544, -37.6592191, -39.965546, -41.0455356, -40.0655795, -37.7496906, -34.1728732, -29.4101313, -23.5468726, -16.6435361, -8.77096437, 0., 10.1187145, 20.7860564, 31.7231031, 42.6509319, 52.6140477, 62.2807294, 71.6426832, 80.6916156, 89.0636486, 97.2483064, 105.379529, 113.591258, 122.654121, 131.810694, 140.94024, 149.922024, 154.37936, 160.149841, 168.815111, 181.956813, 197.339586, 227.99609, 274.768336, 359.778501, 504.715443, 711.721124, 987.824231, 1340.05345}
};
  
  /* Stuff for interpolating the NQC data */
  gsl_spline    *spline = NULL;
  gsl_interp_accel *acc = NULL;
  /* Interpolate the spin NQC data in 2-D parameter space -- spin and mass ratio */
  /* First, interpolate in spin dimension for all mass ratios */
  spline = gsl_spline_alloc( gsl_interp_cspline, adim );
  acc    = gsl_interp_accel_alloc();
  for (i = 0; i < qdim; i++)
  {
    gsl_spline_init( spline, alist, a3stab[i], adim );
    gsl_interp_accel_reset( acc );
    a3slist[i] = gsl_spline_eval( spline, a/(1.0-2.0*eta), acc );
    gsl_spline_init( spline, alist, a4tab[i], adim );
    gsl_interp_accel_reset( acc );
    a4list[i] = gsl_spline_eval( spline, a/(1.0-2.0*eta), acc );
    gsl_spline_init( spline, alist, a5tab[i], adim );
    gsl_interp_accel_reset( acc );
    a5list[i] = gsl_spline_eval( spline, a/(1.0-2.0*eta), acc );
//printf("%.15f\n",a3slist[i]);
  }
//printf("%.15f\n",a);
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
  /* Second, interpolate in mass ratio dimension */
  spline = gsl_spline_alloc( gsl_interp_cspline, qdim );
  acc    = gsl_interp_accel_alloc();
  gsl_spline_init( spline, etalist, a3slist, qdim );
  gsl_interp_accel_reset( acc );
  coeffs->a3S = gsl_spline_eval( spline, eta, acc );
  gsl_spline_init( spline, etalist, a4list, qdim );
  gsl_interp_accel_reset( acc );
  coeffs->a4 = gsl_spline_eval( spline, eta, acc );
  gsl_spline_init( spline, etalist, a5list, qdim );
  gsl_interp_accel_reset( acc );
  coeffs->a5 = gsl_spline_eval( spline, eta, acc );
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
 
  /* Interpolate nonspin NQC data in the mass ratio dimension */
  spline = gsl_spline_alloc( gsl_interp_cspline, nsqdim );
  acc    = gsl_interp_accel_alloc();
  gsl_spline_init( spline, nsetalist, a1list, nsqdim );
  gsl_interp_accel_reset( acc );
  coeffs->a1 = gsl_spline_eval( spline, eta, acc );
  gsl_spline_init( spline, nsetalist, a2list, nsqdim );
  gsl_interp_accel_reset( acc );
  coeffs->a2 = gsl_spline_eval( spline, eta, acc );
  gsl_spline_init( spline, nsetalist, a3list, nsqdim );
  gsl_interp_accel_reset( acc );
  coeffs->a3 = gsl_spline_eval( spline, eta, acc );
  gsl_spline_init( spline, nsetalist, b1list, nsqdim );
  gsl_interp_accel_reset( acc );
  coeffs->b1 = gsl_spline_eval( spline, eta, acc );
  gsl_spline_init( spline, nsetalist, b2list, nsqdim );
  gsl_interp_accel_reset( acc );
  coeffs->b2 = gsl_spline_eval( spline, eta, acc );
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
  /* Andrea and I have different sign conventions, so I need to put a minus sign in front */
  coeffs->b1 = - coeffs->b1;
  coeffs->b2 = - coeffs->b2;
#if 0
  coeffs->a1 = -8.02798637014;
  coeffs->a2 = 48.7446843797;
  coeffs->a3 = -45.7900277224;
  coeffs->a3S= 0.;
  coeffs->a4 = 0.;
  coeffs->a5 = 0.;
  coeffs->b1 = 0.834742923041;
  coeffs->b2 = -2.33512320852; // q=1
#endif
#if 1
  coeffs->a1 = -7.79667;
  coeffs->a2 = 47.182;
  coeffs->a3 = 2238.85334023;
  coeffs->a3S= 0.;
  coeffs->a4 = -7143.16738899;
  coeffs->a5 = 5596.0086893;
  coeffs->b1 = 0.85069;
  coeffs->b2 = -2.47071; // q=1, chi1=chi2=0.98
#endif
#if 0
  coeffs->a1 = -6.82562365707;
  coeffs->a2 = 41.5711482044;
  coeffs->a3 = -39.4329799959;
  coeffs->a3S= 0.;
  coeffs->a4 = 0.;
  coeffs->a5 = 0.;
  coeffs->b1 = 0.461925688819;
  coeffs->b2 = -1.38733263299; // q=8
#endif
#if 0
  coeffs->a1 = -7.5758;
  coeffs->a2 = 46.9323;
  coeffs->a3 = -118.368375152;
  //coeffs->a3 = -45.0036; // NS part 
  coeffs->a3S= 0.;
  coeffs->a4 = 125.555824111;
  coeffs->a5 = -22.0751068073;
  //coeffs->a4 = 0.;
  //coeffs->a5 = 0.;
  coeffs->b1 = 0.51305;
  coeffs->b2 = -1.55133; // q=8, chi1=0.5
#endif
   
  /* Obsolete polynomial fitting of nonspin NQC coefficients a1, a2, a3, b1 and b2 */
  /*
  coeffs->a1 = -12.67955358602124 + 75.41927959573084 * eta - 106.15933052937714 * eta2;
  coeffs->a2 = 101.45522216901628 - 757.3158549733314 * eta + 1473.314771676588 * eta2;
  coeffs->a3 = -107.6647834845902 + 857.6219519536213 * eta - 1776.2776804623143 * eta2;
  // Andrea and I have different sign conventions, so I need to put a minus sign in front 
  coeffs->b1 = - (-1.464129495621165 + 12.81732978488213 * eta - 60.09957767247623 * eta2);
  coeffs->b2 = - ( 7.477426352542122 - 85.26122117590637 * eta + 353.3251639728075 * eta2);
  */

  return XLAL_SUCCESS;

}

/**
 * Function to interpolate known amplitude NQC coeffcients of spin terms for SEOBNRv2,
 * namely a3s, a4 and a5.
 * The a3s, a4 and a5 values were calculated for
 * 11 mass ratios q=1,1.5,2,4,6,8,10,11,12,14 and 20, and
 * 19 spin (\f$\chi\f$ defined in Taracchini et al. PRD 86, 024011 (2012)) values
 * chi = -1, -0.9, -0.8, ......, 0.3, 0.4, 0.5, 0.55, 0.6, 0.65.
 * The calculation was done by Andrea Taracchini using a C++ code of the UMaryland group.
 * In principle, these numbers can be automatically calculated iteratively by the LAL code.
 * However, since such calcualtion increase the cost of each waveform generation by
 * about an order of magnitude, we prepare these numbers in advance reduce cost.
 * These number can be verified by confirming that
 * the peak amplitude and frequency agree well with the NR-fits predicted values,
 * and to get exact NR-fits predicted values, corrections on these numbers are ~1%.
 */
UNUSED static int XLALSimIMRGetEOBCalibratedSpinNQCv2( EOBNonQCCoeffs * restrict coeffs, 
                                    INT4 UNUSED l, 
                                    INT4 UNUSED m, 
                                    REAL8 eta, 
                                    REAL8 a )
{
  const unsigned int nsqdim = 101;
  const unsigned int   qdim = 21;
  const unsigned int   adim = 81;
  UINT4 i;
  /* REAL8 eta2 = eta*eta;*/
  REAL8 a3slist[qdim], a4list[qdim], a5list[qdim];
  REAL8 b3list[qdim], b4list[qdim];

  memset( coeffs, 0, sizeof( *coeffs ) );
  const double nsetalist[101] = {0.021, 0.0233, 0.0256, 0.0279, 0.0302, 0.0324, 0.0347, 0.037, 0.0393, 0.0416, 0.0439, 0.0462, 0.0485, 0.0508, 0.0531, 0.0553, 0.0576, 0.0599, 0.0622, 0.0645, 0.0668, 0.0691, 0.0714, 0.0737, 0.076, 0.0782, 0.0805, 0.0828, 0.0851, 0.0874, 0.0897, 0.092, 0.0943, 0.0966, 0.0989, 0.1012, 0.1034, 0.1057, 0.108, 0.1103, 0.1126, 0.1149, 0.1172, 0.1195, 0.1218, 0.124, 0.1263, 0.1286, 0.1309, 0.1332, 0.1355, 0.1378, 0.1401, 0.1424, 0.1447, 0.147, 0.1492, 0.1515, 0.1538, 0.1561, 0.1584, 0.1607, 0.163, 0.1653, 0.1676, 0.1698, 0.1721, 0.1744, 0.1767, 0.179, 0.1813, 0.1836, 0.1859, 0.1882, 0.1905, 0.1927, 0.195, 0.1973, 0.1996, 0.2019, 0.2042, 0.2065, 0.2088, 0.2111, 0.2134, 0.2156, 0.2179, 0.2202, 0.2225, 0.2248, 0.2271, 0.2294, 0.2317, 0.234, 0.2363, 0.2385, 0.2408, 0.2431, 0.2454, 0.2477, 0.25};
  const double etalist[21] = {0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25};
  const double alist[81]   = {-1., -0.975, -0.95, -0.925, -0.9, -0.875, -0.85, -0.825, -0.8, -0.775, -0.75, -0.725, -0.7, -0.675, -0.65, -0.625, -0.6, -0.575, -0.55, -0.525, -0.5, -0.475, -0.45, -0.425, -0.4, -0.375, -0.35, -0.325, -0.3, -0.275, -0.25, -0.225, -0.2, -0.175, -0.15, -0.125, -0.1, -0.075, -0.05, -0.025, 0., 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475, 0.5, 0.525, 0.55, 0.575, 0.6, 0.625, 0.65, 0.675, 0.7, 0.725, 0.75, 0.775, 0.8, 0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 1.};

  const double a1list[101] = {-9.70365622284, -9.5487796595, -9.40705218072, -9.28242447497, -9.16287060127, -9.03811669547, -8.92394749379, -8.81661027414, -8.7095851328, -8.60723315388, -8.50836995352, -8.40985391052, -8.314500912, -8.22035943611, -8.12770592382, -8.03757606119, -7.94997595927, -7.86421160183, -7.77953168332, -7.69584399308, -7.61324344682, -7.53192249536, -7.45216063215, -7.37374612092, -7.29648777301, -7.22046897906, -7.14581213839, -7.07265232481, -7.00073270258, -6.92991647363, -6.85992514074, -6.79104613524, -6.72328769905, -6.65664446919, -6.59158891131, -6.52783781218, -6.4667874811, -6.41088676175, -6.35876275395, -6.3073217784, -6.25613493225, -6.20587481079, -6.15851113649, -6.1141742151, -6.07327317284, -6.03621329966, -6.0001893395, -5.96655531812, -5.93844433597, -5.91027864135, -5.88037825648, -5.85358373242, -5.83336708498, -5.82744285781, -5.82818819089, -5.8231589417, -5.81812839946, -5.81887431779, -5.82568678472, -5.83621354379, -5.84931204972, -5.86442954488, -5.88109070102, -5.89755497618, -5.91392008194, -5.93177097875, -5.95199881766, -5.97697408745, -6.00401909417, -6.032531535, -6.06148086596, -6.09023163793, -6.11920032326, -6.14820900838, -6.17641365359, -6.20357599269, -6.23098735047, -6.25903214508, -6.28832401322, -6.32033354968, -6.35344345473, -6.38745961465, -6.42231976572, -6.46800466293, -6.50219865357, -6.52990506622, -6.54965342702, -6.58271099648, -6.626558215, -6.67664802777, -6.72810147477, -6.77989789496, -6.83928403282, -6.91446932169, -6.99994425628, -7.08201636107, -7.21472596367, -7.36827339046, -7.55475356656, -7.81665470799, -8.24936799702};

  const double a2list[101] = {71.8355837876, 70.2888163877, 68.858392382, 67.578368761, 66.3483263783, 65.0816770892, 63.9139561883, 62.8042133759, 61.7074812267, 60.6519669001, 59.6324792699, 58.6222094978, 57.6434687679, 56.6791933949, 55.7311388293, 54.8071777364, 53.9062860992, 53.0218701557, 52.1495389357, 51.2895492368, 50.4445863927, 49.6138840037, 48.797074544, 47.9943062521, 47.2056828988, 46.4310888273, 45.6708900349, 44.9258152616, 44.1929201091, 43.4702452949, 42.7537119379, 42.0477065689, 41.3531483593, 40.6709264131, 40.0029858455, 39.3453826346, 38.7105183605, 38.119829408, 37.5632095631, 37.0166865581, 36.4716876499, 35.9345595307, 35.4276417334, 34.9454114602, 34.4897359229, 34.0666726398, 33.6585411373, 33.268152249, 32.9062444867, 32.5494651377, 32.1932422911, 31.8674681969, 31.598287958, 31.4450951286, 31.3504698588, 31.2224679676, 31.1007751184, 31.0231242897, 30.991643584, 30.991534397, 31.0161048282, 31.0616750366, 31.1249606145, 31.1939358163, 31.2696363303, 31.3636652679, 31.4828878696, 31.6456708233, 31.8311325739, 32.034527166, 32.2477681231, 32.4658509298, 32.6919766483, 32.9247097971, 33.1580140019, 33.3900369294, 33.6294910387, 33.8788248397, 34.1421252583, 34.4296046145, 34.7294034713, 35.0398350888, 35.3602204133, 35.7665471711, 36.0881399675, 36.3632831298, 36.5789495976, 36.9023316975, 37.3120873104, 37.7708043813, 38.2382659746, 38.7083475975, 39.2358324951, 39.8798584792, 40.5947669891, 41.2766466634, 42.3250573144, 43.5282809602, 44.9693438982, 46.9606371658, 50.0808316558};

  const double a3list[101] = {-75.3870348213, -73.5563054743, -71.8586610238, -70.3326010764, -68.8669937325, -67.3652081804, -65.9797869993, -64.6591310991, -63.3595036445, -62.1068621881, -60.8985027736, -59.7042086628, -58.5479728539, -57.4103808438, -56.2930281651, -55.2041741738, -54.1420295685, -53.0987898958, -52.0706872299, -51.0586173069, -50.0667002046, -49.0926497207, -48.1344431455, -47.1933388413, -46.2704158358, -45.3650449971, -44.4772042283, -43.6074103996, -42.7520134632, -41.9084025308, -41.0711187682, -40.2460776811, -39.4347613879, -38.6386667672, -37.8587722914, -37.0898368809, -36.3458764907, -35.6510020946, -34.9947343733, -34.3515708492, -33.7091654129, -33.0750000728, -32.4772646605, -31.9060849394, -31.3624600752, -30.8542478078, -30.364610911, -29.8930785724, -29.4449207659, -29.0043606667, -28.5703691139, -28.171109369, -27.8341554058, -27.6229720264, -27.4766951091, -27.2985750163, -27.1291617273, -27.0061494789, -26.9317625784, -26.8914625991, -26.8788122164, -26.8902006455, -26.9223881792, -26.9631792328, -27.013848345, -27.0865681299, -27.1887211003, -27.3400465839, -27.5181167553, -27.7178394005, -27.9305187037, -28.1507835284, -28.3821145422, -28.6229907203, -28.8672132002, -29.1128371338, -29.36865915, -29.6371071707, -29.922311561, -30.2347674515, -30.5620356831, -30.9022473284, -31.254609498, -31.7001465622, -32.0565827322, -32.3648712883, -32.6103616483, -32.9736989432, -33.4312219324, -33.9418551023, -34.4614297518, -34.9838283376, -35.5674613231, -36.2737268484, -37.0519888454, -37.7923554764, -38.9140097, -40.1981922876, -41.7294871636, -43.8332602049, -47.0600574352};

  const double b1list[101] = {-0.454231789035, -0.453740895324, -0.45114951071, -0.445455559729, -0.44030174389, -0.438246415293, -0.435617386371, -0.430628416274, -0.428155263233, -0.424563996869, -0.421496742404, -0.419542314684, -0.417786543725, -0.416452769889, -0.415181542102, -0.413649756645, -0.411562161412, -0.408719678496, -0.40593475615, -0.403713197659, -0.402947057134, -0.402692696086, -0.401882254635, -0.401361543256, -0.401857592043, -0.403135870926, -0.405023096371, -0.407445456424, -0.410163626113, -0.412943367674, -0.415232607237, -0.417748668538, -0.420844132902, -0.424905043844, -0.429142196414, -0.432959586521, -0.436617400517, -0.440401062979, -0.445048680123, -0.451651705843, -0.457760235328, -0.463733230566, -0.472390233324, -0.480892678638, -0.488468389259, -0.496054924285, -0.505287254412, -0.516712248163, -0.530702320051, -0.543402030632, -0.553376187444, -0.564523083803, -0.576716924164, -0.588673429838, -0.599807622576, -0.6083625871, -0.61629220663, -0.626333930393, -0.638492384934, -0.651370663496, -0.663578126675, -0.675399039036, -0.687337552218, -0.699811934302, -0.712809946825, -0.725971047683, -0.73913751946, -0.752013364118, -0.764774460887, -0.777402692471, -0.789886995534, -0.802209798764, -0.814336203683, -0.826241102491, -0.837954892052, -0.849457493421, -0.860601375162, -0.871314573433, -0.881462999336, -0.890680302818, -0.899437045192, -0.907818061982, -0.915897388877, -0.92305068548, -0.931064125926, -0.939248450894, -0.947470287824, -0.953643840407, -0.957725818058, -0.959847025364, -0.960115447844, -0.959051577361, -0.957073291984, -0.952866371063, -0.945294863841, -0.937774627549, -0.921959850725, -0.901024675728, -0.87472187553, -0.834748962688, -0.733155709571};

  const double b2list[101] = {1.24250927283, 1.23820882351, 1.22814744033, 1.20899501569, 1.19187482987, 1.18459068292, 1.17544439367, 1.16014645984, 1.15217431136, 1.14122006693, 1.13199616785, 1.12640121441, 1.121556275, 1.11826795578, 1.11552355075, 1.11238433978, 1.10807035787, 1.10209548232, 1.09668053255, 1.09318287767, 1.09397420715, 1.09653673342, 1.09802015935, 1.10068222122, 1.10643573897, 1.11470254998, 1.12510511039, 1.13755406941, 1.15120517214, 1.16532741744, 1.17824955526, 1.19205557182, 1.20770493667, 1.22624497606, 1.24563734056, 1.26435121321, 1.28310293005, 1.30265970794, 1.32505082997, 1.35327207812, 1.38087484583, 1.40864740792, 1.44367567593, 1.47847338398, 1.51107248182, 1.54412932876, 1.58210883511, 1.62865620041, 1.68826303807, 1.74233039828, 1.78373725125, 1.82868267503, 1.87681370736, 1.92312517747, 1.9656186804, 1.99785596132, 2.02740592217, 2.06495795508, 2.11055684374, 2.15878076242, 2.20440433139, 2.24836632819, 2.29241478185, 2.33795222667, 2.38484419984, 2.43171272272, 2.47786827853, 2.52194361567, 2.56488422632, 2.60668479559, 2.6474159293, 2.6870781889, 2.72548564593, 2.76255774767, 2.79845498187, 2.83311863463, 2.86593650322, 2.89662455284, 2.9246508029, 2.94855630091, 2.97022227107, 2.98999047671, 3.00814380105, 3.02093394891, 3.03817745817, 3.05662290422, 3.07635559396, 3.08574133623, 3.0854061584, 3.07702110074, 3.06226203175, 3.04321477098, 3.01960302068, 2.98621442345, 2.94157099234, 2.89879486315, 2.82337825218, 2.72933462053, 2.61436316871, 2.4455062038, 2.1171325295};

  const double a3stab[21][81] = {
{194.034062, 193.315403, 192.340696, 191.108131, 189.615898, 187.862187, 185.845189, 183.563092, 181.014087, 178.196365, 175.108115, 171.747527, 168.112791, 164.202098, 159.88123, 155.232138, 150.346797, 145.250479, 139.968456, 134.526002, 128.948389, 123.451375, 117.843441, 112.123552, 106.290673, 100.343771, 94.2818106, 88.0767764, 81.7047066, 75.2290824, 68.6598248, 62.0068545, 55.2800926, 48.4894598, 41.7952745, 35.1123156, 28.3439776, 21.4694187, 14.4677972, 7.31827158, 0, -7.06413602, -14.3999825, -22.0896625, -30.2152987, -38.8590141, -48.1029316, -57.9430959, -68.3884858, -79.7270473, -92.0758537, -105.551978, -120.272495, -136.354477, -161.225662, -190.378371, -218.751539, -244.966902, -267.646199, -285.411167, -296.883542, -310.354742, -328.37036, -319.492152, -252.281869, 29.4084482, 249.400555, 588.783155, 1228.64495, 2363.34086, 4134.16095, 7287.44826, 10172.0754, 11510.5972, 12299.3266, 12652.4652, 12433.2947, 11460.0626, 9982.07373, 8101.6811, 5921.23755}, 
{160.112548, 160.667985, 160.855096, 160.683878, 160.164332, 159.306454, 158.120244, 156.615701, 154.802822, 152.691607, 150.292055, 147.614163, 144.66793, 141.463355, 137.974809, 134.234829, 130.278648, 126.123553, 121.786828, 117.285761, 112.637636, 107.962912, 103.161454, 98.2362991, 93.1904845, 88.0270475, 82.7490251, 77.3472032, 71.8142083, 66.1823723, 60.4597069, 54.6542236, 48.7739343, 42.8268505, 36.9526253, 31.0759936, 25.1117255, 19.0409065, 12.8446223, 6.50395828, 0, -6.31861818, -12.8891206, -19.7811827, -27.06448, -34.8086881, -43.0834823, -52.0954881, -62.0307527, -72.6314885, -83.9117648, -95.8856513, -108.567217, -121.970532, -144.794469, -171.559059, -196.126931, -216.73575, -231.62318, -239.026886, -237.184533, -230.598764, -228.655499, -197.207571, -102.107815, 154.463951, 449.595197, 893.840096, 1597.75282, 2688.58018, 4226.79843, 6614.56539, 8915.92748, 10410.9574, 11588.5513, 12442.7012, 12813.2517, 12518.0513, 11617.0807, 10087.775, 7907.56954}, 
{135.75872, 136.138893, 136.299307, 136.238612, 135.955459, 135.448498, 134.71638, 133.757755, 132.571273, 131.155584, 129.50934, 127.63119, 125.519785, 123.173775, 120.539911, 117.649675, 114.538477, 111.215586, 107.690265, 103.971782, 100.069403, 95.9633653, 91.6959723, 87.2804996, 82.7302224, 78.058416, 73.2783556, 68.4234731, 63.5241716, 58.5455294, 53.4926375, 48.3705869, 43.1844686, 37.9393735, 32.7564707, 27.5674199, 22.2950938, 16.9208406, 11.4260085, 5.79194558, 0, -5.65721456, -11.5526034, -17.7478063, -24.304463, -31.2842132, -38.7486967, -47.177716, -56.9882396, -67.2420324, -77.8309468, -88.646835, -99.5815494, -110.526942, -130.898858, -154.564236, -174.669094, -189.157228, -195.972438, -193.058522, -178.359277, -95.540229, -64.7497577, -36.3204737, 39.4150121, 182.864142, 515.868654, 1052.98428, 1808.76676, 2783.21579, 4034.5552, 5797.96238, 7575.84812, 8998.68735, 10279.6823, 11358.5028, 12072.546, 12243.6693, 11858.1609, 10846.2528, 9138.17686}, 
{117.594822, 118.216324, 118.587481, 118.713031, 118.597712, 118.24626, 117.663414, 116.853912, 115.822491, 114.573888, 113.112843, 111.444091, 109.572371, 107.502421, 105.228449, 102.761854, 100.114831, 97.2942708, 94.3070652, 91.1601054, 87.8602826, 84.4097412, 80.820775, 77.1009309, 73.2577557, 69.2987964, 65.2315998, 61.0732898, 56.8395514, 52.5150315, 48.1033884, 43.6082803, 39.0333656, 34.3823025, 29.7493628, 25.0808824, 20.3163379, 15.4408534, 10.439553, 5.29756047, 0, -5.35071899, -10.9029559, -16.6877846, -22.736279, -29.0795131, -35.7485608, -43.2877245, -52.1641879, -61.1818369, -70.1633582, -78.9314385, -87.3087643, -95.1180224, -113.214171, -134.440841, -150.806571, -159.877486, -159.219713, -146.399379, -118.982611, -70.3022292, -22.0562145, 53.8038232, 185.326274, 404.156778, 729.470528, 1193.63258, 1829.00799, 2680.16117, 3742.85913, 5113.35478, 6506.86997, 7696.20069, 8895.67902, 10040.8571, 10932.5801, 11339.6312, 11293.7482, 10737.025, 9611.55528}, 
{104.536193, 105.151992, 105.52303, 105.656338, 105.558946, 105.237884, 104.700185, 103.952877, 103.002992, 101.85756, 100.523612, 99.0081792, 97.3182909, 95.4609783, 93.4570201, 91.3047497, 89.00146, 86.5513695, 83.9586968, 81.2276605, 78.362479, 75.3512364, 72.216514, 68.9647586, 65.6024172, 62.1359365, 58.5717633, 54.9331837, 51.2409534, 47.4612548, 43.5936975, 39.637891, 35.593445, 31.4599689, 27.2943493, 23.0599623, 18.715848, 14.2499006, 9.65001414, 4.90408258, 0, -5.02324517, -10.2359098, -15.6571562, -21.306147, -27.2020446, -33.3640114, -40.5534775, -49.4203385, -58.20891, -66.6369699, -74.4222961, -81.2826666, -86.9358592, -101.298089, -117.635564, -128.44252, -131.350726, -123.991956, -103.997981, -69.0005725, -33.6912932, 21.7514718, 116.165345, 268.387948, 496.490784, 821.201233, 1260.43758, 1832.11811, 2569.02123, 3444.48487, 4457.74625, 5529.91911, 6593.29966, 7712.26975, 8816.62971, 9754.9376, 10360.7413, 10608.8717, 10425.9949, 9738.77691}, 
{96.9661604, 97.1499366, 97.1688322, 97.0244892, 96.7185494, 96.2526546, 95.6284466, 94.8475673, 93.9116586, 92.8223622, 91.58132, 90.190174, 88.6505658, 86.9641373, 85.1097322, 83.1034146, 80.9629742, 78.6947161, 76.3049452, 73.7999668, 71.1860856, 68.4994692, 65.712436, 62.8271668, 59.8458427, 56.7706443, 53.6037525, 50.3545028, 47.0311555, 43.6187841, 40.1166645, 36.5240723, 32.8402837, 29.0645742, 25.2349668, 21.3262257, 17.3101827, 13.178188, 8.92159236, 4.53174618, 0, -4.6801893, -9.51986849, -14.5279779, -19.713458, -25.0852492, -30.6522918, -37.2783443, -45.6987185, -53.8783872, -61.4792071, -68.1630351, -73.591728, -77.4271424, -90.8801985, -106.306773, -115.188531, -114.825063, -102.515961, -75.5608158, -31.2592188, 33.1149764, 110.137717, 217.631759, 373.41986, 572.166076, 888.81781, 1313.4076, 1835.96799, 2443.46734, 3135.13071, 3899.77877, 4762.59709, 5753.91388, 6767.53656, 7738.27797, 8606.45078, 9341.27386, 9790.84147, 9844.54054, 9391.75801}, 
{91.9984129, 92.0724899, 91.9723426, 91.7032712, 91.2705759, 90.6795572, 89.9355152, 89.0437503, 88.0095627, 86.8382527, 85.5351205, 84.1054666, 82.554591, 80.8877943, 79.1289793, 77.2716783, 75.3080151, 73.2394849, 71.067583, 68.7938044, 66.4196446, 63.9380846, 61.3603095, 58.6889905, 55.9267986, 53.0764049, 50.1404805, 47.134602, 44.0724069, 40.9257076, 37.6919349, 34.3685201, 30.9528939, 27.4424876, 23.8547107, 20.1743557, 16.3847026, 12.479096, 8.45088037, 4.29340019, 0, -4.52502743, -9.18698719, -13.9802361, -18.899131, -23.9380288, -29.0912862, -35.37357, -43.6462329, -51.4639546, -58.4068133, -64.0548872, -67.9882546, -69.7869938, -79.8850741, -90.9963626, -95.0131348, -89.2953953, -71.2031485, -38.096399, 12.6648486, 60.4616544, 143.801577, 269.670606, 445.05473, 673.029945, 970.22686, 1338.93951, 1781.46194, 2307.71228, 2897.11226, 3505.21395, 4196.11738, 5020.27394, 5911.38909, 6815.19861, 7663.7143, 8405.30481, 8936.63545, 9171.99775, 9025.68322}, 
{88.7697502, 88.7243429, 88.5229502, 88.1697385, 87.668874, 87.0245232, 86.2408524, 85.322028, 84.2722164, 83.0955839, 81.796297, 80.378522, 78.8464252, 77.2041731, 75.4608637, 73.6175437, 71.6748862, 69.6360488, 67.504189, 65.2824646, 62.9740331, 60.5905035, 58.1254149, 55.5807576, 52.9585221, 50.2606987, 47.4892778, 44.6635083, 41.8000455, 38.8596132, 35.8371944, 32.7277721, 29.5263291, 26.2278486, 22.8396849, 19.3489948, 15.7419984, 12.0111483, 8.14889689, 4.14769663, 0, -4.57373064, -9.27148902, -14.063259, -18.9190244, -23.8087691, -28.7024769, -34.6285481, -42.4563432, -49.6250515, -55.6749058, -60.1461388, -62.5789832, -62.5136718, -69.857389, -77.5921962, -77.915422, -68.2468249, -46.0061634, -8.61319596, 46.512319, 98.1593023, 187.472996, 318.664443, 495.944685, 718.350267, 1002.85596, 1347.46346, 1750.17445, 2217.72359, 2721.91372, 3194.0735, 3742.39868, 4446.91421, 5256.03088, 6117.02483, 6949.94367, 7676.88355, 8238.97487, 8570.85908, 8607.17759}, 
{86.991682, 86.9009299, 86.6584244, 86.2679651, 85.7333516, 85.0583835, 84.2468603, 83.3025818, 82.2293474, 81.0309567, 79.7112095, 78.2739051, 76.7228433, 75.0618236, 73.2860615, 71.4047868, 69.4278794, 67.3608948, 65.2093882, 62.9789152, 60.6750311, 58.3371277, 55.932251, 53.4612835, 50.9251075, 48.3246055, 45.6606599, 42.952423, 40.2163017, 37.4094932, 34.5254615, 31.5576709, 28.4995856, 25.3446698, 22.0706172, 18.6808685, 15.1800581, 11.5648762, 7.83201255, 3.97815717, 0, -4.44925282, -8.98199005, -13.5540843, -18.121408, -22.639834, -27.0652346, -32.4829047, -39.8084298, -46.2971051, -51.4462215, -54.7530699, -55.7149412, -53.8291262, -58.6640354, -63.3459328, -60.2387428, -46.7797921, -20.4064072, 21.4440851, 81.3343584, 146.856546, 244.47561, 378.327141, 552.546728, 759.30391, 1032.25046, 1361.16283, 1735.81745, 2148.15637, 2581.45923, 2984.80058, 3447.78767, 4035.66729, 4736.79845, 5514.39717, 6296.19014, 7018.13887, 7618.75482, 8035.68015, 8206.55699}, 
{87.4052431, 87.2488223, 86.9362164, 86.4718191, 85.8600238, 85.1052242, 84.2118138, 83.184186, 82.0267346, 80.7438529, 79.3399346, 77.8193732, 76.1865622, 74.4458952, 72.5913588, 70.6339301, 68.5853738, 66.4522119, 64.2409669, 61.9581608, 59.6103158, 57.2444895, 54.8210706, 52.3409829, 49.8051505, 47.2144972, 44.5699471, 41.8924106, 39.199795, 36.4452343, 33.6215372, 30.7215125, 27.737969, 24.6637155, 21.4396148, 18.0913369, 14.648484, 11.1144899, 7.49278858, 3.78681398, 0, -4.20404561, -8.43513042, -12.6428882, -16.7769527, -20.7869576, -24.6225368, -29.4472857, -36.2423878, -42.0547553, -46.3411145, -48.5581916, -48.162713, -44.6114048, -47.2420788, -49.2606505, -43.1250808, -26.2710014, 3.86595568, 49.8501587, 114.245976, 188.587294, 293.492819, 432.021615, 607.232747, 806.962399, 1071.81931, 1386.59523, 1736.08191, 2100.7761, 2478.35461, 2845.05278, 3255.50918, 3749.4546, 4345.50524, 5025.01448, 5735.30498, 6443.8395, 7071.55114, 7555.23546, 7831.68799}, 
{90.6686478, 90.338945, 89.8498427, 89.2065534, 88.4142895, 87.4782634, 86.4036875, 85.1957743, 83.8597362, 82.4007857, 80.8241351, 79.1349969, 77.3385836, 75.4401075, 73.4416099, 71.3503094, 69.1736645, 66.9175363, 64.5877858, 62.1902741, 59.7308624, 57.2362237, 54.6885327, 52.0907763, 49.445941, 46.7570139, 44.0269815, 41.2851226, 38.5567652, 35.7820291, 32.9532257, 30.0626665, 27.1026627, 24.0655259, 20.8727666, 17.5614852, 14.1741411, 10.7175274, 7.19843724, 3.62366374, 0, -4.06090718, -8.10254597, -12.0635508, -15.8825563, -19.4981968, -22.8491068, -27.132479, -33.3563927, -38.4501259, -41.8412979, -42.9575283, -41.2264364, -36.0756416, -36.4842853, -35.8376593, -26.798047, -6.83938301, 26.5643983, 75.9393624, 143.811575, 212.08357, 319.941909, 466.480743, 650.794218, 862.613943, 1124.12827, 1423.1964, 1747.67751, 2084.00816, 2424.31545, 2731.24107, 3082.92956, 3537.96935, 4077.76369, 4684.63584, 5326.20908, 5987.9228, 6598.65193, 7106.5548, 7459.78976}, 
{95.4194996, 94.9744199, 94.3458563, 93.5410786, 92.5673567, 91.4319604, 90.1421596, 88.7052241, 87.1284237, 85.4190283, 83.5843078, 81.6315319, 79.5679705, 77.4008935, 75.1326394, 72.7735975, 70.3345307, 67.8237172, 65.2494357, 62.6199646, 59.9435824, 57.2613442, 54.5442252, 51.7959772, 49.020352, 46.2211015, 43.4019773, 40.5995521, 37.8454667, 35.0649948, 32.2485619, 29.3865936, 26.4695154, 23.4877529, 20.3632018, 17.1296398, 13.8260322, 10.4568216, 7.02645092, 3.53936279, 0, -4.2462037, -8.44478288, -12.5002812, -16.3172421, -19.8002093, -22.8537263, -26.6855879, -32.3077596, -36.5085628, -38.6633793, -38.1475908, -34.3365791, -26.605726, -22.9025218, -17.179601, -2.74072752, 22.792071, 61.7967668, 116.651332, 189.733739, 263.596099, 376.57103, 526.669187, 711.901228, 926.647159, 1177.87322, 1459.23331, 1764.38136, 2097.49779, 2420.65693, 2632.95065, 2903.44544, 3359.22954, 3869.45585, 4404.85502, 4965.83074, 5566.44504, 6136.45351, 6634.35268, 7018.6391}, 
{101.846145, 101.147731, 100.262496, 99.1975513, 97.9600035, 96.5569615, 94.9955335, 93.2828279, 91.4259529, 89.432017, 87.3081284, 85.0613955, 82.6989265, 80.2278299, 77.6380801, 74.9476243, 72.1757068, 69.3329406, 66.4299384, 63.4773132, 60.485678, 57.4779833, 54.4508003, 51.4130383, 48.3736061, 45.3414129, 42.3253675, 39.3907795, 36.5944834, 33.8105285, 31.0249232, 28.2236762, 25.3927962, 22.5182917, 19.4346176, 16.2236559, 12.9787608, 9.7169399, 6.45520101, 3.21055181, 0, -3.68378724, -7.22704665, -10.5403555, -13.5342911, -16.1194306, -18.2063515, -21.0997824, -25.8949643, -29.1688988, -30.2660929, -28.5310535, -23.3082876, -13.9423019, -8.20238269, -0.103490751, 16.9771316, 45.4181963, 87.5984154, 145.896501, 222.691165, 302.797812, 421.032832, 575.026689, 762.409844, 975.885659, 1225.23808, 1502.1851, 1798.44471, 2114.12912, 2415.77363, 2618.3873, 2857.47376, 3230.62205, 3675.97438, 4169.32865, 4684.71809, 5189.26697, 5679.50007, 6141.91923, 6563.02631}, 
{109.957867, 108.923771, 107.696307, 106.283107, 104.691802, 102.930024, 101.005405, 98.9255759, 96.6981698, 94.330818, 91.8311523, 89.2068045, 86.4654063, 83.6145897, 80.6454548, 77.5760913, 74.4258405, 71.2057155, 67.9267295, 64.5998959, 61.2362277, 57.7997181, 54.3548943, 50.9192633, 47.510332, 44.1456075, 40.8425968, 37.7092791, 34.8400388, 32.036054, 29.278097, 26.5469404, 23.8233566, 21.088118, 18.1070639, 14.9969343, 11.8907385, 8.81321199, 5.78908988, 2.84310746, 0, -3.28293647, -6.33516022, -9.05356875, -11.3350596, -13.0765301, -14.174878, -16.0101729, -19.7396202, -21.7136831, -21.2270435, -17.5743834, -10.0503846, 2.05027093, 12.7196716, 26.9099491, 50.0813972, 84.3124696, 131.68162, 194.267302, 274.14797, 350.227635, 466.768801, 620.509952, 808.189573, 1022.38726, 1270.10006, 1543.07583, 1833.06243, 2139.56356, 2431.05959, 2635.46524, 2856.25309, 3163.22021, 3542.38549, 3979.04626, 4441.16822, 4893.91506, 5330.99006, 5736.36707, 6094.01992}, 
{119.074672, 117.759823, 116.217882, 114.459582, 112.495658, 110.336844, 107.993873, 105.477481, 102.7984, 99.9673644, 96.9951093, 93.8923682, 90.6698751, 87.3383642, 83.9063765, 80.3860334, 76.7896223, 73.1283257, 69.4133263, 65.6558067, 61.8669495, 57.9512264, 54.0412684, 50.1629958, 46.3423287, 42.6051873, 38.9774916, 35.5978285, 32.587855, 29.7040929, 26.9267158, 24.2358979, 21.6118129, 19.0346348, 16.2844216, 13.4679411, 10.7071094, 8.02303235, 5.43681598, 2.96956625, 0.642389138, -2.31001716, -4.9715306, -7.21243689, -8.90302176, -9.91357094, -10.1143701, -10.8729696, -13.3319374, -13.7814292, -11.4837931, -5.70137737, 4.30346968, 19.2683999, 34.8047905, 54.8931933, 83.9021962, 123.617988, 175.826758, 242.314694, 324.867985, 399.313533, 514.196004, 666.062089, 851.458481, 1063.45811, 1307.17626, 1574.99112, 1859.28093, 2159.87331, 2446.79825, 2659.00754, 2872.53587, 3134.43182, 3455.95783, 3831.93732, 4236.88674, 4657.31296, 5052.06814, 5389.4322, 5637.68509}, 
{128.495326, 127.074107, 125.349104, 123.337657, 121.057102, 118.524777, 115.758021, 112.77417, 109.590562, 106.224536, 102.693428, 99.0145773, 95.2053203, 91.2829951, 87.2836963, 83.2128958, 79.074646, 74.8824479, 70.6498026, 66.3902114, 62.1171752, 57.720481, 53.3544299, 49.049609, 44.836605, 40.746005, 36.8083958, 33.1657559, 29.9433259, 26.9053427, 24.0371649, 21.3241507, 18.7516587, 16.305047, 13.8190975, 11.3744239, 9.06303794, 6.90109728, 4.90475952, 3.09018229, 1.4735232, -0.793048375, -2.71006388, -4.14204299, -4.95350539, -5.00897077, -4.17295882, -3.71878542, -4.70807799, -3.63676357, 0.202654687, 7.51767362, 19.0157901, 35.4045009, 51.4404702, 71.5957135, 100.792698, 140.956113, 194.010651, 261.881001, 346.491854, 436.641211, 558.397569, 710.286632, 890.834104, 1091.36933, 1328.16899, 1591.12322, 1870.12216, 2157.70833, 2435.81478, 2673.27241, 2895.86752, 3111.89571, 3385.10739, 3718.32051, 4076.70075, 4436.6436, 4767.2776, 5040.20866, 5227.04271}, 
{138.900398, 137.337422, 135.391602, 133.087098, 130.44807, 127.498678, 124.263082, 120.765443, 117.02992, 113.080674, 108.941865, 104.637653, 100.192198, 95.6296595, 91.0129641, 86.3417484, 81.6127149, 76.8420946, 72.0461186, 67.2410178, 62.4430231, 57.528747, 52.6733214, 47.9122593, 43.2810742, 38.8152792, 34.5503876, 30.751528, 27.6493265, 24.73026, 21.9366103, 19.2106597, 16.49469, 13.7309833, 9.83315317, 5.39422106, 1.08507419, -2.94159973, -6.533113, -9.53677793, -11.7999068, -12.1167774, -11.5331679, -10.0418217, -7.63548268, -4.30689424, -0.0488000308, 3.44288832, 4.72818055, 7.88680239, 13.6175512, 22.6192243, 35.5906192, 53.2305331, 70.1560327, 90.9132406, 120.508639, 160.884993, 213.985069, 281.75163, 366.127442, 463.170325, 585.968126, 734.458352, 908.578505, 1100.01084, 1329.05576, 1585.74455, 1860.10848, 2144.0068, 2421.9869, 2675.82075, 2901.51559, 3083.49661, 3307.1074, 3583.65042, 3885.38171, 4198.36692, 4492.47968, 4743.11362, 4925.66238}, 
{150.492141, 148.740377, 146.534433, 143.904424, 140.880463, 137.492664, 133.771141, 129.746008, 125.447379, 120.905368, 116.150089, 111.211655, 106.120181, 100.905781, 95.6494103, 90.3490205, 84.9987135, 79.6182038, 74.2272058, 68.8454341, 63.4926031, 58.0444393, 52.6845308, 47.4524779, 42.3878808, 37.5303397, 32.9194547, 28.9855324, 26.100171, 23.3687472, 20.672222, 17.8915564, 14.9077111, 11.601647, 5.59600934, -1.79962227, -9.10472355, -15.9764138, -22.0718122, -27.0480381, -30.5622108, -28.5501535, -24.904218, -19.7954601, -13.3949353, -5.87369959, 2.59719142, 9.71708841, 13.505339, 18.8829889, 26.5436657, 37.1809971, 51.4886108, 70.1601344, 88.6285408, 110.915367, 141.440359, 181.973167, 234.283443, 300.140838, 381.315003, 474.099497, 590.161596, 730.114744, 894.572384, 1077.6267, 1295.97687, 1542.41088, 1809.71674, 2095.08978, 2378.09595, 2639.604, 2862.82661, 3020.70131, 3195.47472, 3406.67752, 3645.53151, 3922.39098, 4199.38914, 4455.09829, 4668.09072}, 
{161.612584, 159.667115, 157.206125, 154.264113, 150.875577, 147.075014, 142.896922, 138.3758, 133.546144, 128.442453, 123.099225, 117.550958, 111.832149, 105.977296, 100.052322, 94.0718457, 88.0481062, 82.0091739, 75.9831191, 69.9980119, 64.0819226, 58.2016273, 52.4549554, 46.878442, 41.5086225, 36.3820321, 31.5352061, 27.2761023, 23.8718937, 20.7101138, 17.7170917, 14.8191564, 11.942637, 9.01386256, 4.4868292, -0.780730804, -5.81963604, -10.4024039, -14.3015518, -17.289597, -19.139057, -17.5768668, -14.7036351, -10.5743881, -5.24415237, 1.23204592, 8.79918041, 15.7981135, 20.8107444, 27.6176602, 36.8151553, 48.9995239, 64.7670605, 84.7140596, 104.239416, 127.227319, 157.955375, 198.082965, 249.269468, 313.174262, 391.456729, 487.744948, 599.674335, 730.018462, 881.550901, 1048.86479, 1254.91206, 1492.64982, 1755.03523, 2037.67667, 2325.57748, 2619.73529, 2861.71536, 2991.38112, 3129.50183, 3302.22954, 3496.20653, 3710.8285, 3933.98728, 4157.49651, 4373.16983}, 
{177.48227, 175.503869, 172.943216, 169.838793, 166.229083, 162.152567, 157.647727, 152.753046, 147.507005, 141.948087, 136.114773, 130.045546, 123.778887, 117.353279, 110.83074, 104.234863, 97.5874593, 90.9221961, 84.2727413, 77.6727629, 71.1559286, 64.7702854, 58.533136, 52.4761623, 46.6310463, 41.0294699, 35.7031149, 30.7482842, 26.2515706, 22.09014, 18.2694361, 14.7949029, 11.671984, 8.90612311, 6.4000957, 4.2242939, 2.45688145, 1.1243021, 0.252999572, -0.130582374, 0, 0.107665064, 0.86020865, 2.36190121, 4.71701321, 8.02981511, 12.4045774, 17.1481302, 21.687129, 28.0326146, 36.6126438, 47.8552734, 62.1885602, 80.0405611, 97.444392, 117.608365, 144.073471, 178.166715, 221.215102, 274.545638, 339.48533, 421.429031, 513.437617, 619.067665, 741.875754, 879.019633, 1049.83962, 1250.21376, 1476.0201, 1723.93173, 1987.44155, 2289.09571, 2536.59024, 2645.32769, 2757.98307, 2901.89402, 3057.58613, 3212.0205, 3373.22828, 3542.05763, 3719.35668}, 
{203.880451, 202.208214, 199.864409, 196.892449, 193.33575, 189.237727, 184.641796, 179.59137, 174.129866, 168.300697, 162.14728, 155.71303, 149.041361, 142.175688, 135.235082, 128.215097, 121.105564, 113.934421, 106.72961, 99.5190709, 92.3307445, 85.1683129, 78.087325, 71.1190715, 64.294843, 57.6459303, 51.2036241, 45.0011085, 39.0712832, 33.4409113, 28.1405147, 23.2006152, 18.6517348, 14.5243953, 10.9819919, 7.97099029, 5.42779827, 3.35575953, 1.75821776, 0.638516679, 0, -0.204729396, 0.0854206845, 0.880801607, 2.19176474, 4.02866144, 6.40184309, 8.80951187, 10.8268265, 13.6887456, 17.6135694, 22.8195984, 29.525133, 37.9484735, 47.1793729, 58.1500548, 71.8781689, 88.8128507, 109.403236, 134.098459, 163.347656, 184.145996, 219.852907, 270.187002, 334.866894, 420.243061, 509.675555, 610.841165, 731.41668, 883.008518, 1061.50458, 1288.88822, 1489.49727, 1597.87824, 1694.42859, 1798.16764, 1908.18382, 2037.43503, 2163.3321, 2275.40424, 2363.18065}
};

  const double a4tab[21][81] = {
{-631.969628, -629.931559, -626.93986, -622.993704, -618.092259, -612.234699, -605.420193, -597.647913, -588.91703, -579.226715, -568.576139, -556.964473, -544.390888, -530.854555, -515.846825, -499.688117, -482.737292, -465.097392, -446.871459, -428.162534, -409.073659, -390.431317, -371.515197, -352.328428, -332.874139, -313.155461, -293.175523, -272.835853, -252.053245, -231.073772, -209.941815, -188.701758, -167.397985, -146.074879, -125.285446, -104.752312, -84.159603, -63.4476689, -42.5568579, -21.4275187, 0, 20.492378, 41.5808063, 63.5035042, 86.4986913, 110.804587, 136.659411, 164.137905, 193.339377, 224.892938, 259.103185, 296.274714, 336.71212, 380.720002, 451.856144, 535.71509, 616.131221, 688.652886, 748.828432, 792.206209, 814.334565, 847.454547, 888.473551, 839.16961, 601.32076, -334.087968, -1051.99176, -2143.56025, -4199.96305, -7855.60243, -13571.9501, -23786.9485, -33106.7677, -37360.6537, -39761.0933, -40691.7902, -39750.1263, -36377.5393, -31425.5653, -25265.023, -18266.7309}, 
{-505.215416, -508.196276, -509.772603, -509.990387, -508.895619, -506.534288, -502.952385, -498.195899, -492.310822, -485.343142, -477.33885, -468.343937, -458.404391, -447.566204, -435.751584, -423.084826, -409.699594, -395.667196, -381.058941, -365.946138, -350.400093, -334.883419, -319.02208, -302.833344, -286.334476, -269.542745, -252.475417, -235.102086, -217.39951, -199.498947, -181.437021, -163.250356, -144.975574, -126.649301, -108.741298, -91.0141828, -73.197789, -55.240145, -37.0892796, -18.6932216, 0, 17.996917, 36.5393224, 55.8235696, 76.0460123, 97.4030037, 120.090898, 144.841731, 172.307048, 201.402322, 232.1064, 264.398133, 298.256369, 333.659957, 398.01185, 473.942308, 542.005708, 596.571523, 632.009229, 642.688299, 622.97821, 589.700061, 568.112915, 450.1228, 127.635743, -717.548639, -1675.26318, -3105.72753, -5369.16134, -8880.19603, -13835.8159, -21538.445, -28962.0227, -33785.1132, -37519.4657, -40137.7205, -41147.9007, -39965.1597, -36840.7475, -31740.8799, -24631.7726}, 
{-416.815775, -419.057598, -420.469505, -421.053329, -420.810899, -419.744046, -417.854601, -415.144394, -411.615257, -407.26902, -402.107514, -396.13257, -389.346018, -381.749689, -373.154365, -363.682735, -353.471949, -342.562916, -330.996543, -318.813739, -306.055411, -292.666119, -278.796426, -264.500547, -249.832696, -234.847087, -219.597935, -204.212585, -188.80739, -173.261705, -157.600052, -141.846951, -126.026926, -110.164496, -94.647616, -79.2708986, -63.8014475, -48.1894474, -32.385083, -16.3385389, 0, 15.8033661, 32.1194751, 49.1192595, 66.9736521, 85.8535854, 105.929992, 128.799096, 155.842954, 183.824465, 212.335847, 240.969319, 269.317099, 296.971406, 353.48309, 419.49243, 473.371969, 508.586135, 518.599357, 496.876061, 436.880677, 151.559476, 37.8684007, -64.67, -316.533178, -780.861526, -1858.22995, -3592.31237, -6026.7827, -9156.64679, -13171.5827, -18819.4885, -24531.9546, -29153.6948, -33275.8654, -36682.7175, -38862.2881, -39230.2155, -37787.9462, -34337.2185, -28679.7708}, 
{-352.810903, -355.83479, -357.93123, -359.12342, -359.43456, -358.887848, -357.506483, -355.313665, -352.332593, -348.586464, -344.098478, -338.891834, -332.989731, -326.415368, -319.156081, -311.257756, -302.768993, -293.720326, -284.14229, -274.065418, -263.520245, -252.52262, -241.11979, -229.344318, -217.228767, -204.805698, -192.107676, -179.198016, -166.135417, -152.8789, -139.44854, -125.864412, -112.14659, -98.3151495, -84.6689645, -71.0517401, -57.2860767, -43.3350226, -29.1616263, -14.728936, 0, 14.8074113, 30.0201505, 45.7103479, 61.9501343, 78.81164, 96.3669955, 116.389741, 140.39776, 164.394917, 187.762516, 209.881862, 230.134257, 247.901005, 297.350154, 355.856786, 398.16278, 416.534098, 403.236704, 350.53656, 250.699628, 81.0333436, -86.8683061, -342.335077, -774.696726, -1485.36227, -2533.86603, -4024.03276, -6059.68721, -8785.088, -12182.7583, -16539.1067, -20994.3126, -24867.3396, -28756.8221, -32428.176, -35236.2064, -36409.007, -36092.1278, -34118.3483, -30320.4477}, 
{-308.350086, -311.217498, -313.195697, -314.31487, -314.605207, -314.096896, -312.820128, -310.805091, -308.081975, -304.680968, -300.632259, -295.966039, -290.712495, -284.901817, -278.614509, -271.848932, -264.599635, -256.886517, -248.729473, -240.148403, -231.163203, -221.738348, -211.956813, -201.846149, -191.433908, -180.747641, -169.8149, -158.71565, -147.52198, -136.136116, -124.564327, -112.812884, -100.888056, -88.7961136, -76.7068069, -64.522988, -52.135134, -39.5160767, -26.6386475, -13.4756781, 0, 13.7532862, 27.8822177, 42.4225625, 57.4100889, 72.8805651, 88.8697592, 107.833319, 131.863276, 155.20919, 176.924281, 196.061768, 211.674871, 222.816808, 260.477512, 303.506922, 328.337507, 327.490095, 293.485516, 218.844597, 96.0881681, -27.8598655, -216.405701, -529.145665, -1025.67608, -1762.42935, -2806.80613, -4214.60606, -6041.62879, -8393.47272, -11176.5412, -14348.9799, -17747.3316, -21230.7092, -24884.6691, -28450.6124, -31439.9222, -33295.7233, -33964.1906, -33220.0041, -30837.8436}, 
{-284.293575, -285.455423, -286.02854, -286.023086, -285.449223, -284.317112, -282.636914, -280.418789, -277.6729, -274.409407, -270.638472, -266.370255, -261.614918, -256.382622, -250.60771, -244.348306, -237.668273, -230.593279, -223.148995, -215.361088, -207.255227, -198.95188, -190.368826, -181.518641, -172.413902, -163.067186, -153.49107, -143.717927, -133.777153, -123.633993, -113.292987, -102.758672, -92.03559, -81.128279, -70.1406469, -59.0143722, -47.6836118, -36.1325801, -24.3454913, -12.3065598, 0, 12.6663563, 25.6213633, 38.8702577, 52.4182761, 66.2706552, 80.4326318, 97.6599695, 120.295132, 141.766533, 160.962606, 176.771785, 188.082505, 193.783199, 229.149613, 270.051435, 289.603997, 279.253008, 230.444179, 134.623223, -16.7641508, -233.206538, -489.797272, -842.946308, -1349.0636, -1986.97268, -3004.46287, -4364.84115, -6031.4145, -7957.99019, -10136.3744, -12487.2395, -15199.5083, -18463.6459, -21789.2919, -24932.21, -27708.8938, -30015.7269, -31368.8032, -31413.8615, -29796.6404}, 
{-269.261202, -269.893529, -269.922537, -269.369733, -268.256623, -266.604716, -264.435518, -261.770537, -258.63128, -255.039254, -251.015968, -246.582927, -241.76164, -236.573613, -231.101996, -225.329301, -219.233375, -212.823116, -206.107424, -199.095198, -191.795336, -184.183849, -176.307067, -168.178431, -159.811381, -151.21936, -142.415809, -133.451427, -124.371315, -115.099826, -105.635271, -95.975965, -86.1202194, -76.0663476, -65.8536243, -55.4544499, -44.8381242, -33.9945817, -22.9137572, -11.5855851, 0, 12.2126797, 24.6532377, 37.2820741, 50.0595888, 62.9461818, 75.902253, 92.1509699, 114.42523, 134.883784, 152.162243, 164.896218, 171.721322, 171.273165, 196.398655, 224.091212, 228.754614, 202.026835, 135.545848, 20.9496279, -150.123852, -310.660599, -586.222934, -997.812488, -1566.43089, -2298.62462, -3251.05128, -4427.36643, -5831.22566, -7492.21992, -9336.19869, -11173.0819, -13320.9669, -16036.5623, -18968.5552, -21908.6681, -24641.4471, -26991.1174, -28630.3656, -29283.5396, -28674.987}, 
{-260.19861, -260.329225, -259.924836, -259.002871, -257.580756, -255.67592, -253.305787, -250.487786, -247.239344, -243.577888, -239.520843, -235.085639, -230.2897, -225.150455, -219.702272, -213.951861, -207.904649, -201.574597, -194.975667, -188.121822, -181.027023, -173.729491, -166.215579, -158.495896, -150.581055, -142.481666, -134.20834, -125.822067, -117.376265, -108.761087, -99.9666886, -90.9832271, -81.8008583, -72.4097385, -62.81914, -53.0071261, -42.9503135, -32.6349483, -22.0472768, -11.1735453, 0, 12.4164605, 25.0318929, 37.7317017, 50.4012913, 62.9260662, 75.1914308, 90.4528171, 111.45927, 130.038073, 144.706292, 153.980991, 156.379236, 150.418091, 167.309568, 184.88779, 178.512015, 140.011827, 61.2168106, -66.0434513, -249.939375, -421.336883, -714.638337, -1141.92747, -1715.28801, -2427.92313, -3338.48847, -4436.41128, -5711.11882, -7182.28864, -8748.59704, -10135.8843, -11813.1906, -14138.1958, -16810.2724, -19624.9526, -22323.0395, -24636.7153, -26387.3746, -27365.9987, -27363.5686}, 
{-255.977882, -255.881463, -255.267847, -254.15314, -252.553448, -250.484879, -247.963539, -245.005536, -241.626976, -237.843965, -233.672612, -229.129022, -224.229303, -218.989561, -213.39858, -207.489751, -201.298535, -194.846627, -188.155723, -181.247519, -174.14371, -166.976683, -159.642156, -152.146537, -144.496235, -136.697658, -128.757215, -120.733546, -112.677443, -104.470422, -96.0976845, -87.54443, -78.7958597, -69.837174, -60.5898415, -51.07938, -41.3361316, -31.3583329, -21.1442203, -10.6920305, 0, 12.0834578, 24.2514848, 36.3470461, 48.2131066, 59.6926314, 70.6285854, 84.449911, 104.046714, 120.68748, 132.759149, 138.648667, 136.742975, 125.429017, 134.94212, 143.521823, 126.996883, 77.2399175, -13.8764533, -154.479611, -352.696935, -567.399046, -886.363643, -1321.38453, -1884.25553, -2545.61083, -3418.77101, -4466.13874, -5650.11671, -6941.4336, -8277.51408, -9446.43994, -10847.9371, -12784.0016, -15099.507, -17645.0049, -20184.1708, -22494.7141, -24385.0854, -25653.2552, -26097.1941}, 
{-259.301517, -258.908377, -257.986116, -256.552805, -254.626518, -252.225327, -249.367304, -246.070523, -242.353056, -238.232975, -233.728354, -228.857265, -223.63778, -218.087972, -212.193491, -205.992919, -199.527296, -192.821326, -185.899714, -178.787165, -171.508383, -164.22357, -156.80322, -149.253326, -141.57988, -133.788873, -125.886297, -117.935423, -109.990914, -101.921802, -93.710822, -85.3407098, -76.7941999, -68.0540274, -58.9282556, -49.5101175, -39.9060675, -30.1345682, -20.2140821, -10.1630719, 0, 11.3798946, 22.6798009, 33.7261314, 44.3452987, 54.3637154, 63.6077941, 75.75625, 93.9089438, 108.680992, 118.33465, 121.13217, 115.335806, 99.2078126, 102.39748, 103.311507, 77.9807358, 18.24746, -84.0460249, -237.057424, -448.944442, -690.660789, -1032.33718, -1482.18179, -2048.4028, -2687.14723, -3535.04044, -4537.81771, -5641.21434, -6777.12204, -7932.80681, -8992.50289, -10228.8405, -11850.2394, -13815.456, -16038.211, -18346.2706, -20623.88, -22615.3228, -24112.5969, -24907.7002}, 
{-272.576592, -271.509326, -269.905439, -267.78581, -265.171318, -262.08284, -258.541257, -254.567445, -250.182285, -245.406655, -240.261433, -234.767498, -228.945729, -222.817004, -216.396127, -209.707821, -202.777266, -195.626582, -188.277893, -180.753317, -173.074977, -165.336801, -157.479185, -149.514334, -141.454452, -133.311744, -125.098412, -116.905882, -108.813673, -100.644566, -92.3785988, -83.9958093, -75.4762357, -66.7999162, -57.7079606, -48.3315539, -38.8199673, -29.2021095, -19.5068898, -9.76321692, 0, 11.0216661, 21.8299644, 32.2208922, 41.9904468, 50.9346252, 58.8494248, 69.5324439, 86.1799924, 99.0197821, 106.223032, 105.960959, 96.4047834, 75.725723, 72.5747565, 65.8414949, 32.1083658, -36.6877954, -148.610153, -311.721872, -534.086117, -754.965378, -1105.54732, -1581.21004, -2177.33164, -2855.74813, -3694.57455, -4648.93882, -5673.96884, -6720.15437, -7756.5378, -8621.99595, -9668.11695, -11166.9082, -12949.962, -14934.2741, -17017.0813, -19145.3929, -21086.8863, -22670.375, -23724.6727}, 
{-290.803667, -289.274472, -287.128883, -284.394981, -281.100843, -277.274549, -272.944177, -268.137807, -262.883517, -257.209386, -251.143493, -244.713917, -237.948737, -230.876032, -223.51659, -215.903102, -208.068811, -200.043286, -191.856098, -183.536816, -175.115011, -166.733726, -158.293387, -149.807891, -141.291137, -132.757024, -124.219449, -115.792997, -107.577122, -99.3449735, -91.0695665, -82.7239183, -74.2810457, -65.7139654, -56.7644909, -47.551899, -38.2129678, -28.7680048, -19.2373174, -9.64121332, 0, 11.7686104, 23.2513321, 34.1374748, 44.1163481, 52.8772616, 60.109525, 69.661843, 84.7579215, 95.1414814, 98.8129842, 93.7728915, 78.021665, 49.5597661, 33.8139514, 11.434708, -39.0035261, -125.110111, -254.494408, -434.765777, -673.533578, -911.692133, -1278.54529, -1765.94828, -2365.75631, -3055.88021, -3862.5713, -4760.95166, -5726.14334, -6770.29687, -7757.44895, -8304.76651, -9080.37076, -10594.8809, -12284.7376, -14027.3107, -15837.8677, -17770.7606, -19586.105, -21144.3248, -22305.8438}, 
{-314.860794, -312.399881, -309.315105, -305.634028, -301.38421, -296.593211, -291.288593, -285.497916, -279.248741, -272.568628, -265.485138, -258.025831, -250.218269, -242.090012, -233.621957, -224.871184, -215.898305, -206.740427, -197.434653, -188.018091, -178.527844, -169.058221, -159.581225, -150.12606, -140.721932, -131.398046, -122.183608, -113.282612, -104.872792, -96.565407, -88.3186921, -80.0908819, -71.8402108, -63.5249134, -54.6077698, -45.360441, -36.0920914, -26.8622969, -17.730633, -8.75667548, 0, 10.1361592, 19.7071329, 28.4245933, 36.0002126, 42.145663, 46.5726166, 53.4597832, 66.3146467, 74.1676824, 74.9168018, 66.4599163, 46.6949373, 13.5197764, -8.04792363, -37.2662924, -95.4463037, -190.237163, -329.288076, -520.248248, -770.766884, -1027.86707, -1411.75897, -1912.81564, -2521.41019, -3210.41698, -4013.37242, -4899.65133, -5838.62854, -6829.70559, -7752.17711, -8281.29834, -8959.95466, -10188.4093, -11658.3146, -13267.0828, -14935.01, -16552.2472, -18106.642, -19552.7363, -20845.072}, 
{-344.679856, -341.022052, -336.721252, -331.8069, -326.308441, -320.255316, -313.67697, -306.602846, -299.062388, -291.085039, -282.700244, -273.937445, -264.826086, -255.395611, -245.632885, -235.594288, -225.33942, -214.906433, -204.33348, -193.658713, -182.920285, -172.031953, -161.173444, -150.40009, -139.767223, -129.330176, -119.144279, -109.546635, -100.832006, -92.3819813, -84.1374849, -76.0394417, -68.0287763, -60.0464131, -51.3358712, -42.2792565, -33.3114654, -24.5160711, -15.9766467, -7.77676527, 0, 9.05698653, 17.2769261, 24.3294609, 29.8842333, 33.6108856, 35.1790599, 39.0222633, 48.8581746, 52.9654922, 49.0795762, 34.9357869, 8.26948448, -33.1839709, -69.8652883, -117.845006, -194.845637, -307.595725, -462.823816, -667.258454, -927.628182, -1172.06384, -1551.03265, -2052.14941, -2663.0289, -3356.3428, -4156.56543, -5033.37999, -5956.46971, -6923.29205, -7820.20756, -8370.28056, -8997.37166, -9998.44616, -11243.2334, -12663.9936, -14157.567, -15604.3705, -16985.934, -18249.317, -19341.579}, 
{-377.73975, -373.104753, -367.713572, -361.606407, -354.823458, -347.404927, -339.391013, -330.821917, -321.73784, -312.178982, -302.185543, -291.797725, -281.055728, -269.999751, -258.680042, -247.130445, -235.384047, -223.478992, -211.453428, -199.3455, -187.193353, -174.725993, -162.333401, -150.096417, -138.095882, -126.412635, -115.127518, -104.688851, -95.4897345, -86.7323237, -78.3482499, -70.2691445, -62.4266393, -54.752366, -46.39279, -37.7762419, -29.3904859, -21.3277529, -13.6802736, -6.54027905, 0, 8.0344215, 14.9827519, 20.4508465, 24.0445606, 25.3697495, 24.0322684, 24.5292753, 30.6229514, 30.2234886, 20.9507164, 0.424463843, -33.7354397, -83.909165, -135.640282, -201.959857, -297.173412, -427.104905, -597.578294, -814.417538, -1083.44659, -1323.86307, -1698.16163, -2193.72605, -2797.94011, -3485.15224, -4274.9, -5138.92486, -6048.9683, -7003.2791, -7894.07684, -8485.75838, -9099.45099, -9946.79918, -10994.753, -12210.4725, -13512.3499, -14852.487, -16096.9162, -17141.0119, -17880.1487}, 
{-411.395283, -406.424621, -400.436968, -393.495203, -385.662204, -377.00085, -367.574018, -357.444586, -346.675432, -335.329434, -323.469471, -311.158421, -298.459161, -285.43457, -272.235267, -258.868624, -245.335373, -231.680445, -217.94877, -204.185281, -190.434908, -176.384163, -162.485896, -148.83454, -135.524528, -122.650289, -110.306257, -98.9695798, -89.0598987, -79.7565246, -70.9984936, -62.7248418, -54.8746053, -47.3868203, -39.3994688, -31.3583368, -23.7698461, -16.736882, -10.3623296, -4.74907387, 0, 5.82942382, 10.3061333, 13.044481, 13.6588196, 11.7635017, 6.97287981, 3.62121897, 5.33356082, 0.438408409, -13.3663218, -38.3827132, -76.9128493, -131.258814, -184.276109, -250.568796, -346.212957, -477.488314, -650.674588, -872.051499, -1147.89877, -1439.87244, -1836.72835, -2332.59123, -2921.58585, -3572.41332, -4341.90989, -5193.6922, -6091.37687, -7008.53331, -7878.91978, -8574.4393, -9231.47707, -9918.11911, -10798.793, -11871.3261, -13019.5422, -14158.9815, -15192.3213, -16026.6241, -16568.9523}, 
{-447.704261, -442.311139, -435.636192, -427.765098, -418.783533, -408.777175, -397.831702, -386.03279, -373.466118, -360.217362, -346.372199, -332.016307, -317.235363, -302.115045, -286.898914, -271.57277, -256.110459, -240.565366, -224.990874, -209.440365, -193.967224, -178.224329, -162.720881, -147.565575, -132.867107, -118.734172, -105.275468, -93.0029713, -82.3680642, -72.5151463, -63.3891673, -54.9350772, -47.097826, -39.8223634, -32.2517941, -24.8383183, -18.0948322, -12.1302966, -7.05367222, -2.97391977, 0, 3.9144123, 6.20741074, 6.47237482, 4.30268407, -0.708281985, -8.96714383, -16.6217659, -20.4599418, -31.0714639, -50.5921446, -81.1577962, -124.904231, -183.967261, -239.909578, -307.881626, -404.591142, -536.382, -709.598075, -930.583241, -1205.68137, -1520.42345, -1920.72148, -2405.59936, -2974.08102, -3596.46265, -4342.62963, -5177.1329, -6064.52341, -6976.43908, -7856.16996, -8624.20503, -9304.38621, -9874.12416, -10584.9548, -11466.2867, -12424.6604, -13408.4535, -14321.5634, -15083.0508, -15611.977}, 
{-486.985302, -481.063369, -473.621951, -464.766211, -454.601306, -443.232398, -430.764646, -417.30321, -402.95325, -387.819926, -372.008398, -355.623825, -338.771369, -321.556188, -304.281376, -286.926879, -269.457662, -251.938399, -234.433764, -217.008433, -199.727079, -182.251363, -165.104632, -148.40722, -132.27946, -116.841686, -102.214231, -88.9322037, -77.46839, -66.9513455, -57.3329911, -48.5652476, -40.6000359, -33.3892771, -26.1172825, -19.2215655, -13.1977443, -8.15474813, -4.20150593, -1.44694686, 0, 2.46866343, 3.07511633, 1.37368948, -3.08128632, -10.7354803, -22.0345617, -33.804155, -43.4138398, -59.9652257, -85.3738362, -121.555195, -170.424826, -233.898252, -294.951034, -367.480185, -466.816392, -598.74919, -769.068117, -983.56271, -1248.02251, -1549.17979, -1927.42098, -2384.28343, -2921.30444, -3517.40583, -4229.90974, -5033.21508, -5901.72077, -6825.70807, -7731.92885, -8541.16776, -9225.3043, -9715.70976, -10261.4576, -10921.8731, -11668.0294, -12532.7126, -13390.1094, -14168.4748, -14796.0635}, 
{-523.530987, -517.18179, -509.099572, -499.404296, -488.215926, -475.654428, -461.839764, -446.8919, -430.930799, -414.076427, -396.448746, -378.167723, -359.353319, -340.125501, -320.743696, -301.239643, -281.634523, -262.019777, -242.486841, -223.127154, -204.032154, -185.1605, -166.754747, -148.924671, -131.780048, -115.430653, -99.9862631, -85.8913394, -73.5400536, -62.2419086, -51.970787, -42.7005711, -34.4051436, -27.0583868, -19.9670615, -13.5270739, -8.18482821, -4.05066127, -1.23490977, 0.152089566, 0, 0.63444541, -0.827294392, -4.83197899, -11.826368, -22.2572209, -36.5712973, -52.069956, -66.5271897, -87.9107707, -117.944593, -158.352549, -210.858533, -277.186437, -340.946978, -415.322511, -514.866481, -645.007681, -811.174906, -1018.79695, -1273.30261, -1585.92326, -1950.09153, -2374.72401, -2868.73727, -3413.26677, -4085.75571, -4861.78357, -5716.92985, -6636.49509, -7566.89564, -8502.24447, -9262.53211, -9658.4135, -10079.7924, -10610.8568, -11207.0811, -11865.721, -12545.62, -13218.7166, -13856.9493}, 
{-573.544328, -567.3913, -559.264478, -549.297849, -537.625402, -524.381123, -509.699, -493.713021, -476.557174, -458.365445, -439.271823, -419.410296, -398.91485, -377.919473, -356.673279, -335.237426, -313.664358, -292.064516, -270.548338, -249.226266, -228.208738, -207.724455, -187.749263, -168.377269, -149.702582, -131.819308, -114.821554, -98.9907231, -84.5800736, -71.2358692, -58.9761693, -47.8190334, -37.782521, -28.8846914, -20.80112, -13.7665228, -8.04154082, -3.71428576, -0.872869407, 0.394596441, 0, -0.352641544, -2.79107424, -7.65091515, -15.2677814, -25.9772899, -40.115058, -55.4388162, -70.0936536, -90.5792079, -118.277801, -154.571754, -200.84339, -258.47503, -314.439123, -379.233729, -464.447671, -574.410692, -713.452534, -885.902939, -1096.09165, -1361.90342, -1660.1381, -2002.53274, -2400.82438, -2844.97238, -3400.43194, -4052.80709, -4787.70187, -5593.56535, -6447.46639, -7421.88042, -8212.53356, -8539.75889, -8877.49863, -9315.97155, -9789.15075, -10251.9446, -10731.7661, -11231.3108, -11753.2741}, 
{-655.005025, -650.22828, -643.153049, -633.930991, -622.713764, -609.653026, -594.900437, -578.607655, -560.926338, -542.008144, -522.004733, -501.067762, -479.348891, -456.999777, -434.44623, -411.666479, -388.618003, -365.396384, -342.097206, -318.816053, -295.648506, -272.707411, -250.068705, -227.82559, -206.071264, -184.898927, -164.401779, -144.67741, -125.822753, -107.920505, -91.0620844, -75.338907, -60.8423896, -47.6639487, -36.2708463, -26.5167378, -18.2268289, -11.4156602, -6.09777208, -2.28770517, 0, 0.899018952, 0.226125436, -2.05369074, -5.97543976, -11.5741318, -18.8847771, -26.3284314, -32.5686664, -41.4996417, -53.8116848, -70.1951234, -91.3402849, -117.937497, -147.032956, -181.622281, -224.976933, -278.532618, -343.72504, -421.989903, -514.762911, -581.305808, -694.756475, -854.432254, -1059.65049, -1330.79573, -1615.21967, -1937.52009, -2322.2948, -2806.55382, -3377.65823, -4109.87456, -4751.12504, -5083.61292, -5376.10302, -5691.28314, -6024.29668, -6414.84462, -6791.94206, -7122.33809, -7372.7818}
};

  const double a5tab[21][81] = {
{519.523893, 518.131239, 515.8465, 512.673834, 508.617401, 503.681359, 497.869868, 491.187084, 483.637169, 475.224279, 465.952575, 455.826214, 444.849356, 433.026159, 419.876634, 405.711214, 390.876979, 375.477116, 359.614813, 343.393256, 326.915633, 310.97921, 294.897239, 278.677049, 262.325971, 245.851335, 229.260471, 212.468357, 195.403847, 178.295097, 161.186937, 144.124194, 127.151698, 110.314275, 94.0951208, 78.2617512, 62.5485019, 46.9105375, 31.3030226, 15.6811219, 0, -14.8448983, -29.9685477, -45.5456425, -61.7508773, -78.7589465, -96.7445447, -115.834203, -136.161689, -158.016862, -181.593973, -207.087271, -234.691007, -264.599431, -315.508826, -375.908989, -432.889195, -482.859257, -522.228987, -547.408199, -554.806705, -574.962906, -596.142264, -541.649647, -334.789923, 439.372044, 1023.20821, 1899.29847, 3550.22271, 6493.90189, 11106.8927, 19378.4641, 26907.6181, 30293.5484, 32120.0118, 32707.9753, 31765.1304, 28865.1246, 24731.0574, 19694.1679, 14085.6953}, 
{401.023924, 404.584455, 406.799776, 407.720454, 407.397055, 405.880146, 403.220293, 399.468062, 394.674018, 388.888729, 382.162761, 374.54668, 366.091051, 356.846442, 346.759978, 335.947662, 324.533328, 312.588699, 300.185499, 287.395453, 274.290286, 261.318544, 248.123086, 234.723595, 221.139753, 207.391243, 193.497745, 179.435614, 165.187713, 150.877329, 136.541738, 122.218216, 107.944036, 93.7564757, 80.056375, 66.6510166, 53.320162, 40.0267222, 26.733608, 13.4037303, 0, -12.7748555, -25.7999156, -39.2144432, -53.157701, -67.768952, -83.1874591, -100.053584, -118.932391, -138.764961, -159.487092, -181.034587, -203.343245, -226.348868, -271.660391, -325.50311, -372.461762, -408.039085, -427.737821, -427.060708, -401.510488, -366.642495, -338.58577, -231.772114, 39.3666742, 733.614349, 1509.27518, 2659.77438, 4478.53714, 7303.46344, 11294.5542, 17504.0752, 23491.4724, 27390.0091, 30357.5262, 32364.0939, 33034.5068, 31900.5993, 29211.3213, 24971.1506, 19184.5649}, 
{320.853928, 323.645506, 325.647365, 326.8671, 327.312305, 326.990576, 325.909508, 324.076696, 321.499735, 318.186219, 314.143744, 309.379905, 303.902297, 297.718514, 290.663545, 282.854178, 274.420263, 265.404701, 255.850392, 245.800236, 235.297134, 224.310145, 212.966208, 201.31842, 189.419881, 177.323689, 165.082942, 152.819472, 140.644784, 128.447624, 116.253182, 104.086648, 91.9732114, 79.938062, 68.2956295, 56.8881296, 45.5358834, 34.2049193, 22.8612652, 11.4709494, 0, -10.9706148, -22.1747344, -33.7312581, -45.7590852, -58.377115, -71.7042469, -87.0660335, -105.606713, -124.55994, -143.554672, -162.21987, -180.184491, -197.077497, -236.11891, -282.013897, -317.681883, -337.92647, -337.551259, -311.359852, -254.155851, -10.8720224, 90.9067418, 180.871687, 388.714058, 763.065479, 1633.56374, 3032.62932, 4992.68271, 7505.62609, 10725.4348, 15243.7119, 19832.0917, 23594.4221, 26919.4491, 29614.2054, 31277.3477, 31430.3293, 30111.7931, 27184.0668, 22509.478}, 
{264.604952, 268.01221, 270.555054, 272.259659, 273.152196, 273.258839, 272.605761, 271.219136, 269.125135, 266.349933, 262.919703, 258.860616, 254.198848, 248.96057, 243.144059, 236.793135, 229.953732, 222.657727, 214.937001, 206.823432, 198.348899, 189.538801, 180.432392, 171.062448, 161.461742, 151.663049, 141.699143, 131.629487, 121.509533, 111.308241, 101.047548, 90.7493939, 80.4357158, 70.1284518, 60.0684983, 50.1392797, 40.2076459, 30.2507491, 20.2457412, 10.1697742, 0, -10.1716849, -20.49872, -31.0198002, -41.7736206, -52.7988762, -64.1342619, -77.2252061, -93.3057587, -109.050954, -123.92831, -137.40534, -148.949562, -158.028492, -191.573447, -231.678243, -258.357804, -265.462155, -246.841316, -196.345311, -107.824161, 37.7516773, 181.758735, 395.308243, 749.511432, 1325.47786, 2169.65629, 3365.15583, 4995.08563, 7176.1872, 9890.67249, 13346.7352, 16906.2956, 20068.6758, 23231.2194, 26180.9312, 28396.171, 29233.1402, 28845.3183, 27114.3841, 23922.0167}, 
{226.936794, 230.072047, 232.394903, 233.936573, 234.728264, 234.801186, 234.186548, 232.915558, 231.019425, 228.529359, 225.476568, 221.89226, 217.807645, 213.253932, 208.310342, 202.977713, 197.253243, 191.158322, 184.714336, 177.942674, 170.864725, 163.458199, 155.794193, 147.900128, 139.803422, 131.531496, 123.11177, 114.613821, 106.100893, 97.4996009, 88.8202487, 80.0731387, 71.2685737, 62.4168564, 53.6474839, 44.895356, 36.0863504, 27.2063895, 18.2413957, 9.17729174, 0, -9.32918894, -18.7963188, -28.4120651, -38.1871037, -48.13211, -58.2577598, -70.546693, -86.685238, -101.968879, -115.607611, -126.811428, -134.790326, -138.754299, -162.950269, -190.74978, -203.62919, -195.677397, -160.983299, -93.6357938, 12.2762194, 119.095052, 277.884525, 535.748638, 939.791391, 1533.80335, 2373.06163, 3500.69416, 4959.82882, 6835.43808, 9045.11621, 11519.1146, 14208.7872, 17071.2823, 20064.301, 22950.6412, 25337.2994, 26758.9674, 27195.3575, 26474.4695, 24424.3036}, 
{208.120791, 209.580113, 210.509836, 210.922517, 210.830712, 210.24698, 209.183878, 207.653963, 205.669793, 203.243925, 200.388916, 197.117324, 193.441706, 189.374619, 184.867629, 179.971877, 174.743121, 169.206393, 163.386727, 157.309156, 150.998712, 144.5579, 137.923582, 131.11009, 124.131758, 117.00292, 109.737907, 102.365653, 94.9128929, 87.3590567, 79.7125504, 71.9817802, 64.1751523, 56.3010729, 48.4310367, 40.53354, 32.5723035, 24.5428292, 16.4406192, 8.26117543, 0, -8.46945348, -17.0127789, -25.6176185, -34.2716146, -42.9624094, -51.6776452, -62.6192298, -77.6563512, -91.4820897, -103.185025, -111.853735, -116.576801, -116.442801, -139.244831, -165.912887, -175.204164, -160.335968, -114.525607, -30.9903895, 97.0523769, 277.669417, 490.299248, 779.530743, 1189.95278, 1700.99006, 2517.96931, 3607.28304, 4935.32375, 6461.28195, 8173.15611, 9967.4643, 12094.3003, 14791.2017, 17529.173, 20080.8976, 22307.089, 24120.9478, 25137.6962, 25074.1402, 23647.0861}, 
{196.99112, 197.852928, 198.189355, 198.02151, 197.370502, 196.257439, 194.703432, 192.729589, 190.35702, 187.606833, 184.500137, 181.058042, 177.301656, 173.252089, 168.982587, 164.481277, 159.732339, 154.746218, 149.53336, 144.104209, 138.469211, 132.610235, 126.570248, 120.363641, 114.004805, 107.50813, 100.888009, 94.1862905, 87.4406999, 80.5999704, 73.6673442, 66.6460632, 59.5393694, 52.3505049, 45.099863, 37.7798357, 30.3815167, 22.9046398, 15.3489389, 7.71414766, 0, -8.13710569, -16.3066832, -24.4615817, -32.5546505, -40.5387386, -48.3666954, -58.602855, -73.4191447, -86.5240549, -96.8100879, -103.169746, -104.495532, -99.6799483, -114.611839, -131.105722, -128.932611, -101.463129, -42.0679019, 55.8824482, 199.017297, 332.674538, 559.635715, 895.530662, 1355.98921, 1943.25865, 2705.97899, 3643.72114, 4756.05602, 6064.72725, 7502.78762, 8874.38743, 10537.2961, 12783.959, 15205.7867, 17604.9312, 19811.9288, 21677.5237, 22943.2431, 23387.8477, 22790.0979}, 
{190.832536, 191.182498, 191.070976, 190.515396, 189.533185, 188.141768, 186.358572, 184.201023, 181.686546, 178.832569, 175.656516, 172.175814, 168.40789, 164.370168, 160.095692, 155.592007, 150.865481, 145.930345, 140.800832, 135.491174, 130.015603, 124.4074, 118.659118, 112.782357, 106.788719, 100.689805, 94.4972168, 88.259642, 82.0201968, 75.7018023, 69.3010013, 62.8143362, 56.2383498, 49.5695846, 42.8059539, 35.9431333, 28.9766944, 21.9028994, 14.7180105, 7.41828999, 0, -8.32900583, -16.675172, -24.9333518, -32.9983983, -40.7651645, -48.1285038, -57.6696175, -71.566051, -83.289292, -91.6434499, -95.4326343, -93.4609546, -84.5325202, -93.2601919, -102.121369, -91.6398416, -55.3407527, 13.2507551, 120.609539, 273.210457, 414.532412, 654.724254, 1002.23974, 1465.53261, 2035.81164, 2764.06811, 3638.06196, 4645.55311, 5800.59664, 7012.06701, 8011.00495, 9284.86209, 11212.7157, 13428.9288, 15737.5879, 17929.6157, 19775.3061, 21140.6326, 21858.9509, 21763.6162}, 
{188.572122, 188.667871, 188.320343, 187.545776, 186.360407, 184.780474, 182.822216, 180.50187, 177.835674, 174.839867, 171.530685, 167.924367, 164.03715, 159.885273, 155.464333, 150.803623, 145.934004, 140.875935, 135.649876, 130.276286, 124.775626, 119.260585, 113.646656, 107.941562, 102.153024, 96.2887649, 90.3565074, 84.4012916, 78.4625502, 72.4587743, 66.3825339, 60.226399, 53.9829398, 47.6447264, 41.1445932, 34.5128992, 27.7845253, 20.9642599, 14.0568915, 7.06720873, 0, -8.0991995, -16.1338685, -23.9667389, -31.460543, -38.4780125, -44.8818797, -53.3854886, -66.2738445, -76.5935407, -83.0498692, -84.3481219, -79.1935907, -66.2915677, -69.5515217, -71.7344903, -53.6936662, -8.97903516, 68.8594169, 186.271704, 349.70784, 524.871315, 784.98631, 1138.2079, 1592.69116, 2121.10369, 2819.1361, 3652.35886, 4586.34238, 5594.64432, 6618.87351, 7444.53231, 8495.24079, 10098.3722, 12019.8806, 14111.2512, 16180.113, 18033.9887, 19524.6826, 20488.903, 20763.3579}, 
{192.640596, 192.425679, 191.760187, 190.661977, 189.148906, 187.238831, 184.949609, 182.299097, 179.305152, 175.985632, 172.358392, 168.44129, 164.252184, 159.80893, 155.105186, 150.174117, 145.050721, 139.757805, 134.318175, 128.754639, 123.090004, 117.461761, 111.762194, 105.99827, 100.176957, 94.3052222, 88.3900346, 82.4792313, 76.6145088, 70.7051111, 64.7414116, 58.7137839, 52.6126016, 46.4282381, 40.0060334, 33.4281099, 26.7875683, 20.1044475, 13.3987865, 6.69062437, 0, -7.58044865, -14.9752011, -22.0361381, -28.61514, -34.5640876, -39.7348614, -47.0396638, -58.9308472, -67.9427133, -72.6845537, -71.76566, -63.7953237, -47.3828366, -46.0843521, -42.7276859, -18.2522863, 33.8352296, 120.028244, 246.820141, 420.704303, 616.581004, 894.493792, 1259.59357, 1717.03125, 2227.35983, 2905.73808, 3703.79987, 4573.17903, 5454.46437, 6332.42487, 7077.50902, 7998.51352, 9337.6381, 10966.8088, 12792.4299, 14674.4371, 16510.3169, 18093.8205, 19254.0633, 19820.1611}, 
{205.201042, 204.348048, 203.041196, 201.300742, 199.146942, 196.600052, 193.680328, 190.408026, 186.803402, 182.886711, 178.67821, 174.198155, 169.466801, 164.504404, 159.33076, 153.966415, 148.431952, 142.747722, 136.934076, 131.011362, 124.999932, 118.983059, 112.90948, 106.790855, 100.638845, 94.4651098, 88.2813103, 82.1587966, 76.15995, 70.1537049, 64.1274852, 58.0687151, 51.9648188, 45.8032204, 39.3720536, 32.7848145, 26.1700843, 19.5560505, 12.9709002, 6.44282098, 0, -7.3490753, -14.4156901, -21.0308295, -27.0254786, -32.2306226, -36.4772465, -42.7807182, -53.6779141, -61.3855962, -64.4417882, -61.3845132, -50.7517947, -31.0816562, -25.2471131, -16.391784, 14.1827643, 72.9160298, 166.247511, 300.616705, 482.463112, 660.681848, 945.390279, 1331.37205, 1813.41082, 2356.24576, 3028.77002, 3789.71416, 4597.80879, 5408.08129, 6190.37179, 6777.22915, 7544.01912, 8786.78629, 10268.4158, 11898.0197, 13595.0589, 15311.8408, 16858.9378, 18095.2411, 18879.6421}, 
{221.690131, 220.398247, 218.58656, 216.281616, 213.509957, 210.298129, 206.672674, 202.660137, 198.287063, 193.579994, 188.565475, 183.270051, 177.720264, 171.94266, 165.966042, 159.815525, 153.516051, 147.093703, 140.574562, 133.984708, 127.350225, 120.795949, 114.235566, 107.681522, 101.146258, 94.6422182, 88.1818448, 81.8547985, 75.7391366, 69.6626667, 63.6064793, 57.5516645, 51.4793127, 45.3705141, 39.013682, 32.5117949, 25.9824161, 19.4460464, 12.9231866, 6.43433754, 0, -8.04107527, -15.7543762, -22.8871409, -29.1866072, -34.4000132, -38.2745968, -43.8795544, -53.7849183, -59.7947447, -60.3074499, -53.7214501, -38.4351615, -12.8470002, 2.69274262, 23.672487, 67.242313, 139.493879, 246.518844, 394.408865, 589.255602, 780.905091, 1078.63779, 1474.37621, 1960.04287, 2514.26833, 3161.76182, 3878.49565, 4640.44212, 5456.18109, 6203.86222, 6522.65532, 7063.54625, 8330.475, 9738.76737, 11163.6008, 12630.3005, 14190.7865, 15640.4513, 16862.1399, 17738.6977}, 
{243.041167, 240.916442, 238.268739, 235.124196, 231.508948, 227.44913, 222.97088, 218.100331, 212.863622, 207.286886, 201.396261, 195.217881, 188.777884, 182.102404, 175.188008, 168.079538, 160.824073, 153.453798, 146.000897, 138.497553, 130.975952, 123.529313, 116.120354, 108.772831, 101.510497, 94.3571074, 87.336416, 80.6075209, 74.3091827, 68.1415349, 62.0733782, 56.0735135, 50.1107416, 44.1538631, 37.7686538, 31.1788696, 24.6387738, 18.1996027, 11.9125925, 5.82897951, 0, -6.83643511, -13.134385, -18.661234, -23.1843668, -26.4711679, -28.2890216, -31.9863919, -40.3736424, -44.6553882, -43.1449803, -34.15577, -16.0011085, 13.005653, 32.7119582, 58.6200168, 107.88517, 186.661095, 301.10147, 457.359973, 661.59028, 867.621722, 1179.24059, 1586.96251, 2081.30312, 2637.14399, 3283.56464, 3992.32, 4735.16499, 5510.68417, 6210.14351, 6525.55868, 6991.93862, 8010.82803, 9232.1407, 10550.4892, 11905.9158, 13206.7859, 14442.7419, 15575.7181, 16567.6489}, 
{269.088009, 265.924239, 262.224134, 258.015465, 253.326007, 248.183532, 242.615814, 236.650625, 230.315739, 223.638929, 216.647968, 209.370629, 201.834685, 194.067909, 186.073201, 177.89407, 169.575906, 161.151569, 152.653921, 144.115822, 135.570133, 126.971984, 118.442702, 110.025882, 101.765121, 93.7040138, 85.886157, 78.5745132, 71.9990828, 65.6789304, 59.5685817, 53.6225628, 47.7953994, 42.0416175, 35.7514821, 29.2364724, 22.8507801, 16.6643459, 10.7471106, 5.16901506, 0, -6.10211882, -11.4602514, -15.8094327, -18.884698, -20.4210821, -20.1536201, -21.6460347, -27.8867423, -29.4559488, -24.5341191, -11.3017181, 12.0607894, 47.3729384, 78.7156359, 119.129938, 182.999637, 275.772557, 402.896523, 569.81936, 781.988894, 978.133378, 1286.266, 1694.78546, 2192.09049, 2753.06761, 3399.44548, 4103.40834, 4837.14042, 5597.6513, 6282.64914, 6628.25706, 7059.49464, 7882.44393, 8911.69889, 10073.7181, 11286.1547, 12446.5037, 13541.6988, 14528.1828, 15362.399}, 
{297.538563, 293.556225, 288.943238, 283.736585, 277.973249, 271.690214, 264.924463, 257.712979, 250.092744, 242.100743, 233.773957, 225.149371, 216.263967, 207.154729, 197.88078, 188.465097, 178.928981, 169.304887, 159.625269, 149.922583, 140.229283, 130.357193, 120.589868, 110.990234, 101.621217, 92.545741, 83.8267321, 75.8132842, 68.8114909, 62.2000155, 55.9255896, 49.9349446, 44.1748124, 38.5919243, 32.4983618, 26.2423386, 20.2201084, 14.5082154, 9.18320368, 4.32161718, 0, -5.43618466, -9.94023776, -13.1965411, -14.8894766, -14.7034259, -12.322771, -11.3742724, -14.8923013, -13.1345613, -4.18470192, 13.8736276, 42.956778, 84.9810999, 128.019953, 182.746814, 260.883015, 367.176375, 506.374714, 683.225853, 902.477611, 1096.46704, 1401.4294, 1805.90828, 2298.44724, 2855.37901, 3495.36723, 4192.30257, 4920.07567, 5676.13651, 6363.69775, 6757.39355, 7186.55806, 7877.68362, 8738.18483, 9726.41403, 10777.6599, 11849.6135, 12833.4199, 13642.9101, 14191.9155}, 
{325.967927, 321.759463, 316.699275, 310.843686, 304.249017, 296.971589, 289.067725, 280.593744, 271.605969, 262.160722, 252.314323, 242.123094, 231.643357, 220.931433, 210.135985, 199.25492, 188.279152, 177.246116, 166.193245, 155.157973, 144.177733, 133.034222, 122.055929, 111.315608, 100.886012, 90.8398934, 81.2500054, 72.4853239, 64.8703136, 57.7694243, 51.135133, 44.9199163, 39.076251, 33.5566141, 27.6663412, 21.767293, 16.2703157, 11.2602536, 6.82195087, 3.0402516, 0, -3.83572388, -6.53633503, -7.79301253, -7.29693548, -4.73928297, 0.188765894, 3.97288046, 3.67320435, 8.74015049, 21.0348622, 42.4184827, 74.7521554, 119.897023, 163.757313, 218.288572, 296.654215, 403.979242, 545.38865, 726.007439, 950.960608, 1187.36036, 1510.91855, 1915.91188, 2396.61706, 2924.84282, 3550.2867, 4240.26406, 4962.09027, 5692.37005, 6370.55074, 6867.86278, 7344.32832, 7893.98241, 8608.93485, 9476.692, 10400.6984, 11305.6616, 12115.2992, 12753.6661, 13144.8175}, 
{355.924553, 351.463197, 345.92638, 339.389563, 331.928211, 323.617786, 314.533751, 304.75157, 294.346705, 283.394619, 271.970776, 260.150639, 248.009671, 235.623334, 223.220949, 210.780647, 198.268917, 185.72975, 173.207141, 160.745081, 148.387565, 135.896053, 123.636089, 111.690687, 100.142859, 89.0756176, 78.5719749, 69.0240383, 60.7774699, 53.1762019, 46.1777443, 39.7396075, 33.8193017, 28.3743371, 22.7181984, 17.2158088, 12.2808414, 8.00253643, 4.47013431, 1.77287539, 0, -2.47767474, -3.61515922, -3.08588707, -0.563291916, 4.27919262, 11.7681329, 18.7689569, 22.6671699, 32.0653158, 48.6952942, 74.2890049, 110.578348, 159.295223, 205.2998, 260.997097, 340.069247, 447.699122, 589.069595, 769.363538, 993.763824, 1249.0487, 1575.51435, 1971.68525, 2436.08585, 2942.24165, 3550.34111, 4228.91015, 4946.47472, 5678.36192, 6372.69412, 6946.10806, 7452.8416, 7901.32133, 8470.01808, 9175.75662, 9940.06053, 10715.5055, 11425.7601, 12004.3475, 12384.791}, 
{387.430025, 382.686339, 376.665071, 369.457889, 361.156459, 351.852449, 341.637525, 330.603354, 318.841603, 306.44394, 293.50203, 280.107541, 266.352139, 252.327492, 238.313641, 224.283086, 210.194068, 196.099725, 182.053193, 168.107609, 154.316108, 140.455871, 126.894102, 113.72205, 101.030964, 88.9120906, 77.4566788, 67.070536, 58.1122035, 49.9207821, 42.4597989, 35.6927807, 29.5832546, 24.0947477, 18.578586, 13.385577, 8.91286911, 5.24920986, 2.48334686, 0.704027705, 0, -1.49811338, -1.47300968, 0.434489138, 4.58356111, 11.3333843, 21.0431367, 31.1274123, 39.4432642, 53.390707, 74.5245146, 104.399461, 144.570319, 196.591864, 246.451664, 305.552364, 386.474976, 493.958418, 632.741605, 807.563455, 1023.16289, 1267.6883, 1576.12367, 1949.31661, 2388.11474, 2873.75188, 3455.4562, 4110.53896, 4816.31143, 5564.37504, 6289.17058, 6910.67669, 7431.21712, 7813.68352, 8241.97583, 8760.16318, 9344.65662, 10021.3633, 10685.864, 11278.1692, 11738.2894}, 
{415.644597, 410.80384, 404.496223, 396.826118, 387.897897, 377.815933, 366.684597, 354.608262, 341.6913, 328.038083, 313.752984, 298.940374, 283.704625, 268.15011, 252.526896, 236.847188, 221.112163, 205.396391, 189.774444, 174.320894, 159.110312, 144.161103, 129.611761, 115.544616, 102.041996, 89.1862286, 77.0596429, 65.993561, 56.2818907, 47.411588, 39.3638818, 32.1200012, 25.6611749, 19.9686321, 14.4900332, 9.54414462, 5.49012131, 2.41832923, 0.419134351, -0.417097379, 0, -0.19189763, 1.31994565, 4.89557566, 10.8950383, 19.6783793, 31.6056445, 44.4793832, 56.4864388, 74.1021239, 98.7249103, 131.75327, 174.585674, 228.620594, 280.369375, 340.648141, 421.399902, 527.068169, 662.096454, 830.928265, 1038.00711, 1291.98541, 1588.48038, 1934.65271, 2337.66309, 2781.05892, 3330.24684, 3964.0518, 4661.29879, 5409.67817, 6161.41866, 6903.2878, 7499.09516, 7801.54674, 8123.45042, 8532.17509, 8990.76641, 9496.34206, 10014.0548, 10519.9545, 10990.0912}, 
{453.17581, 448.837173, 442.811713, 435.216307, 426.167836, 415.783177, 404.179211, 391.472816, 377.78087, 363.220254, 347.907847, 331.960526, 315.495171, 298.628662, 281.605983, 264.462973, 247.225773, 229.98506, 212.831509, 195.855796, 179.148598, 162.951874, 147.184124, 131.915128, 117.214671, 103.152535, 89.7985016, 77.3576024, 66.0145451, 55.5157189, 45.8759913, 37.1102299, 29.2333022, 22.2600757, 15.9192601, 10.4067475, 5.94009081, 2.59268896, 0.43794067, -0.450755246, 0, 0.43664773, 2.55661682, 6.63037865, 12.9284046, 21.7211662, 33.2791347, 45.7879692, 57.7465969, 74.4105134, 96.8966909, 126.322102, 163.803718, 210.458512, 255.584784, 307.776041, 376.520412, 465.35228, 577.806029, 717.416043, 887.716706, 1103.5243, 1345.46792, 1623.23514, 1946.51352, 2306.45731, 2758.46946, 3289.99774, 3888.48991, 4543.93457, 5236.15697, 6023.08053, 6654.59158, 6900.07118, 7152.52932, 7486.30119, 7845.51724, 8191.2698, 8547.00849, 8914.92465, 9297.20964}, 
{514.044111, 511.244973, 506.460749, 499.824339, 491.468637, 481.526543, 470.130952, 457.414762, 443.51087, 428.552173, 412.671569, 396.001954, 378.676226, 360.827281, 342.838575, 324.6845, 306.320485, 287.828178, 269.289225, 250.785274, 232.397974, 214.298224, 196.466094, 178.970904, 161.881976, 145.268632, 129.200193, 113.748283, 98.9841805, 84.9717004, 71.779229, 59.4751529, 48.1278582, 37.8057315, 28.8385816, 21.1294181, 14.5614618, 9.14962737, 4.90882923, 1.85398193, 0, -0.744376298, -0.243394776, 1.53252249, 4.61295341, 9.02747592, 14.8056679, 20.7043685, 25.671661, 32.7803872, 42.576898, 55.6075441, 72.4186764, 93.5566457, 116.626628, 144.035572, 178.41306, 220.907037, 272.665445, 334.836229, 408.567332, 461.903817, 552.18437, 679.028025, 842.053812, 1057.64344, 1284.06812, 1541.06193, 1848.35894, 2235.51207, 2692.79891, 3282.76076, 3795.71593, 4050.6964, 4271.96945, 4511.18603, 4762.95519, 5057.70467, 5339.62931, 5582.32165, 5759.37424}
};

  const double b3tab[21][81] = {
{43.0068882, 42.0198432, 41.0202058, 40.0082996, 38.9844479, 37.9489742, 36.902202, 35.8444546, 34.7760554, 33.697328, 32.6085957, 31.510182, 30.4024104, 29.2856042, 28.1351348, 26.9671105, 25.7995281, 24.6378149, 23.487398, 22.3537045, 21.2421617, 20.2361252, 19.2523312, 18.2854444, 17.3301295, 16.3810513, 15.4328744, 14.4741566, 13.494373, 12.5027911, 11.4965554, 10.4728101, 9.42869978, 8.36136873, 7.43944692, 6.55159635, 5.5734983, 4.46722121, 3.19483351, 1.71840363, 0, -1.82864689, -3.99856215, -6.57110885, -9.60765004, -13.1695488, -17.3181681, -23.1599594, -31.6443369, -40.3337383, -48.8651873, -56.8757074, -64.0023224, -69.8820558, -135.302016, -221.215459, -283.946851, -310.625471, -288.380599, -204.341513, -45.6374928, 1094.39585, 1322.02807, 1275.54747, 1593.24234, 5319.21786, 7157.42945, 10633.1226, 19271.5426, 33460.8266, 66137.5454, 158302.586, 233321.855, 232304.72, 214245.066, 191781.421, 164045.903, 133506.452, 102794.137, 74484.4135, 51152.7378}, 
{34.8479771, 34.0828232, 33.2979234, 32.4946199, 31.6742548, 30.8381703, 29.9877087, 29.124212, 28.2490225, 27.3634824, 26.4689337, 25.5667189, 24.6581798, 23.7446589, 22.8204984, 21.8914687, 20.9638698, 20.0404758, 19.1240605, 18.2173979, 17.323262, 16.48684, 15.6626349, 14.8475632, 14.0385411, 13.2324851, 12.4263115, 11.6119553, 10.7820998, 9.94557265, 9.10131288, 8.24825948, 7.38535144, 6.51152776, 5.87655893, 5.32070703, 4.66524731, 3.85781326, 2.84603836, 1.57755611, 0, -1.69213164, -3.83216388, -6.50655697, -9.80177113, -13.8042666, -18.6005036, -25.4550677, -35.4555167, -45.8712791, -56.3104583, -66.381158, -75.6914817, -83.8495329, -163.679532, -268.472799, -345.978518, -380.829042, -357.656722, -261.093911, -75.7729602, 1382.32395, 1669.1229, 1567.52846, 1860.4452, 5398.05478, 7863.99151, 12521.8726, 22635.3153, 38866.4163, 72283.3546, 159204.733, 230487.863, 232792.773, 220036.698, 203067.024, 179727.633, 149951.728, 118318.723, 87824.0906, 61463.3019}, 
{27.3016278, 26.3201609, 25.4316053, 24.6262982, 23.8945769, 23.2267788, 22.613241, 22.044301, 21.510296, 21.0015632, 20.5084399, 20.0212635, 19.5303712, 19.0261003, 18.4401786, 17.80002, 17.1374748, 16.4548684, 15.7545261, 15.0387734, 14.3099356, 13.5774014, 12.8354574, 12.0854534, 11.3287394, 10.5666652, 9.80058079, 9.02693396, 8.2429092, 7.46157766, 6.68627963, 5.92035536, 5.16714512, 4.42998917, 4.08113488, 3.89055052, 3.60027976, 3.13820596, 2.43221249, 1.41018271, 0, -1.63570072, -3.83620844, -6.7060606, -10.3497946, -14.8719479, -20.377058, -28.3310151, -39.9951502, -52.2188502, -64.5538995, -76.5520826, -87.7651839, -97.7449876, -189.784157, -310.45969, -400.009772, -440.857681, -415.426698, -306.140102, -95.4211725, 1828.26796, 2199.66409, 1975.70095, 2113.31228, 4971.62859, 8048.84049, 13984.5044, 25418.1767, 42730.9598, 75337.7724, 154867.656, 220221.658, 224088.91, 215027.204, 202272.431, 182777.578, 154477.671, 123662.487, 93773.5817, 68252.5085}, 
{19.3140409, 18.8127435, 18.3268192, 17.8544265, 17.3937237, 16.9428693, 16.5000216, 16.0633391, 15.6309802, 15.2011032, 14.7718666, 14.3414288, 13.9079482, 13.4695833, 13.0132296, 12.5441705, 12.0685417, 11.5868053, 11.0994235, 10.6068582, 10.1095716, 9.60921725, 9.10490122, 8.59692112, 8.08557452, 7.57115897, 7.05397204, 6.53172619, 6.00252231, 5.47283926, 4.94402423, 4.41742444, 3.89438708, 3.37625936, 3.26677856, 3.31273845, 3.23047363, 2.93902574, 2.35743642, 1.40474729, 0, -1.64066316, -3.93633285, -7.00899908, -10.9806518, -15.9732812, -22.108877, -31.1013368, -44.4253548, -58.3964885, -72.4903626, -86.1826014, -98.9488294, -110.264671, -213.913405, -349.711289, -450.335908, -495.973035, -466.808444, -343.027907, -104.817199, 1622.64393, 2013.60032, 1968.70522, 2388.61189, 5741.17106, 9061.29224, 15130.2506, 26729.3214, 45277.9892, 77642.9008, 149564.093, 208953.624, 215315.125, 207096.358, 192638.252, 171985.692, 147259.988, 120572.795, 94039.0621, 69773.735}, 
{12.7881556, 12.5733124, 12.333449, 12.0708547, 11.7878189, 11.4866307, 11.1695797, 10.838955, 10.497046, 10.146142, 9.78853229, 9.42650626, 9.06235319, 8.69836239, 8.35563103, 8.02455051, 7.69408857, 7.36268754, 7.02878975, 6.69083752, 6.34727319, 5.96573091, 5.57971601, 5.19192562, 4.8050569, 4.42180699, 4.04487306, 3.68477709, 3.35086534, 3.0271248, 2.71307548, 2.40823739, 2.11213051, 1.82427487, 2.01413104, 2.38393289, 2.60034188, 2.56675554, 2.18657142, 1.36318706, 0, -1.7361858, -4.24175745, -7.64969563, -12.092981, -17.7045944, -24.6175163, -34.7775333, -49.8580361, -65.6573797, -81.5724874, -97.0002825, -111.337688, -123.981628, -235.651602, -381.647502, -489.599547, -538.179962, -506.060977, -371.914819, -114.413714, 1427.64716, 1823.98063, 1923.64426, 2575.69561, 6217.63839, 9780.37322, 16019.0778, 27688.9297, 46640.1621, 78342.7864, 145553.488, 200330.975, 204985.626, 196733.976, 183520.692, 164854.785, 142255.175, 117636.798, 92779.5835, 69463.4617}, 
{8.67660924, 8.32370542, 8.00786559, 7.72441883, 7.46869424, 7.23602091, 7.02172791, 6.82114435, 6.6295993, 6.44242186, 6.25494111, 6.06248614, 5.86038604, 5.64396989, 5.36030655, 5.0352547, 4.69832613, 4.35472116, 4.0096401, 3.66828328, 3.33585102, 3.07039037, 2.81695643, 2.57345101, 2.33777595, 2.10783307, 1.88152419, 1.63946016, 1.36484999, 1.09694022, 0.840653377, 0.600911977, 0.382638539, 0.190755583, 0.61184877, 1.28287839, 1.79677541, 2.03948795, 1.89696414, 1.25515211, 0, -1.78730835, -4.44302392, -8.10816197, -12.9237377, -19.0307665, -26.5702635, -37.9782779, -55.3460025, -73.3267645, -91.1297216, -107.964032, -123.038852, -135.563341, -253.40089, -407.025591, -518.686891, -565.369623, -524.05862, -371.738712, -85.3947335, 1642.39205, 2114.71184, 2223.45834, 2860.52525, 5983.2463, 9855.43663, 16647.5077, 28529.8711, 46469.7646, 76247.1207, 140366.291, 192371.438, 195823.168, 187334.435, 174699.255, 157267.568, 136599.006, 113858.512, 90345.5841, 67359.7194}, 
{6.43126977, 6.27531737, 6.08715296, 5.87043508, 5.62882228, 5.3659731, 5.08554609, 4.7911998, 4.48659277, 4.17538354, 3.86123066, 3.54779269, 3.23872816, 2.93769562, 2.6778981, 2.44430419, 2.21964611, 2.00153933, 1.78759931, 1.57544153, 1.36268145, 1.11185147, 0.860495329, 0.611073717, 0.366047309, 0.127876785, -0.100977178, -0.325939738, -0.551251114, -0.755594685, -0.933307868, -1.07872808, -1.18619273, -1.25003925, -0.612569755, 0.319398433, 1.08969179, 1.57060448, 1.63443064, 1.15346443, 0, -1.83248891, -4.60962345, -8.47584522, -13.5755958, -20.0533167, -28.0534496, -40.4072577, -59.5422772, -79.1784521, -98.3692881, -116.168291, -131.628966, -143.804819, -261.626797, -414.639536, -524.070431, -566.498467, -518.502627, -356.661892, -57.5552467, 1284.42393, 1733.93985, 2047.08474, 2979.95078, 6611.84084, 10434.936, 16793.1685, 28030.4706, 45936.7797, 74518.013, 131486.776, 178779.249, 186637.764, 180981.383, 168083.622, 149637.232, 129364.855, 107450.345, 85283.9669, 64255.9838}, 
{5.29448937, 5.04323248, 4.77966155, 4.50506682, 4.22073848, 3.92796675, 3.62804184, 3.32225395, 3.0118933, 2.69825009, 2.38261454, 2.06627686, 1.75052725, 1.43665593, 1.13009179, 0.829506878, 0.53326, 0.241794828, -0.0444449594, -0.325015685, -0.599473671, -0.889773659, -1.16998018, -1.43655617, -1.68596459, -1.91466838, -2.11913049, -2.29959983, -2.45575642, -2.57501055, -2.65228792, -2.68251428, -2.66061534, -2.58151683, -1.71284295, -0.509613024, 0.518100112, 1.22660722, 1.47221906, 1.1112464, 0, -2.14884821, -5.31551114, -9.6238405, -15.197688, -22.1609055, -30.6373446, -43.7251104, -64.075391, -84.7002585, -104.515922, -122.438591, -137.384475, -148.269783, -263.962577, -413.823124, -518.888817, -555.586124, -500.34151, -329.581443, -19.7323883, 1253.13033, 1722.01745, 2095.46885, 3082.02438, 6471.78476, 10305.4467, 16589.4126, 27330.0851, 43794.8627, 70207.1597, 123503.677, 168501.148, 178068.16, 173279.753, 159892.662, 141050.456, 122301.905, 102188.135, 81639.254, 61585.3689}, 
{4.84351675, 4.51636659, 4.18454238, 3.84828129, 3.50782051, 3.16339724, 2.81524864, 2.46361192, 2.10872426, 1.75082284, 1.39014485, 1.02692748, 0.661407912, 0.293823333, -0.0856447143, -0.470397806, -0.853076352, -1.23138637, -1.60303387, -1.96572487, -2.31716538, -2.67235893, -3.00932511, -3.32338105, -3.60984386, -3.86403063, -4.0812585, -4.26199265, -4.40592472, -4.49606216, -4.52563179, -4.48786044, -4.37597493, -4.18320209, -3.08069225, -1.58172083, -0.261790799, 0.717722075, 1.19544202, 1.00999325, 0, -2.44350491, -5.97249041, -10.6865168, -16.6851445, -24.0679338, -32.934445, -46.6456857, -68.0725919, -89.5162313, -109.751912, -127.554941, -141.700627, -150.964279, -263.324215, -408.47352, -507.839091, -537.859663, -474.973969, -295.620744, 23.7612774, 1419.05114, 1935.38659, 2283.16973, 3172.80264, 6234.3338, 9909.7101, 16012.9004, 26357.8738, 42320.9255, 67029.0451, 114093.299, 154438.174, 165687.996, 163355.177, 151948.899, 134753.464, 116800.69, 97357.2159, 77655.0931, 58926.3714}, 
{4.82260682, 4.44699608, 4.06271369, 3.67030988, 3.27033492, 2.86333905, 2.44987252, 2.0304856, 1.60572854, 1.17615158, 0.742304977, 0.304738993, -0.135996125, -0.579350123, -1.0338924, -1.49330357, -1.95057395, -2.40328793, -2.84902991, -3.2853843, -3.70993548, -4.15038129, -4.57003383, -4.96231861, -5.32066118, -5.63848706, -5.90922177, -6.13521168, -6.31746272, -6.42806964, -6.45683582, -6.39356464, -6.22805947, -5.9501237, -4.63367847, -2.8479176, -1.23136316, 0.0388453389, 0.785568377, 0.831666437, 0, -2.65228347, -6.43086025, -11.4071197, -17.652451, -25.2382436, -34.2358866, -48.1992661, -70.1589809, -91.8593735, -111.957828, -129.111727, -141.978456, -149.215398, -255.278804, -391.897309, -482.791038, -504.977092, -435.47257, -251.294572, 70.5398011, 1522.76762, 2063.48009, 2388.61223, 3194.09905, 6026.91803, 9483.77153, 15281.8373, 25138.2933, 40758.8599, 63895.0867, 104532.188, 139701.079, 151597.138, 151500.514, 142684.992, 127836.263, 110418.605, 91370.037, 72489.8072, 55577.1629}, 
{5.12537602, 4.67830998, 4.22600324, 3.7684454, 3.30562605, 2.8375348, 2.36416125, 1.88549499, 1.40152563, 0.91224276, 0.417635984, -0.0823050996, -0.587590892, -1.09823179, -1.62810779, -2.16845532, -2.70946099, -3.24829829, -3.78214069, -4.30816171, -4.82353482, -5.36894761, -5.89204989, -6.38400556, -6.83597852, -7.23913267, -7.58463191, -7.87544036, -8.11274905, -8.25950615, -8.30208429, -8.22685609, -8.02019418, -7.66847118, -6.16973237, -4.13556985, -2.25238776, -0.708712917, 0.306927866, 0.606007777, 0, -2.81622665, -6.7603839, -11.8667879, -18.1697546, -25.7036003, -34.5026409, -48.1725756, -69.7824605, -90.8270269, -109.890495, -125.557084, -136.411015, -141.036506, -236.673829, -359.497024, -438.212308, -451.224619, -376.938899, -193.760087, 119.906876, 1347.48589, 1851.26156, 2237.22252, 3111.35744, 5837.1687, 9152.11593, 14571.197, 23609.4097, 37622.2754, 58603.2218, 96635.7552, 129169.689, 138671.701, 138033.642, 130473.028, 117553.059, 101287.315, 83542.1718, 66234.1913, 51279.9336}, 
{6.07046285, 5.585625, 5.08344396, 4.56462091, 4.02985707, 3.47985361, 2.91531175, 2.33693266, 1.74541756, 1.14146763, 0.52578407, -0.100931924, -0.737979157, -1.38465643, -2.06000494, -2.75083438, -3.44246004, -4.13014258, -4.80914267, -5.47472097, -6.12213813, -6.78464611, -7.41426743, -8.00101589, -8.53490532, -9.0059495, -9.40416227, -9.73401699, -9.99781433, -10.1509936, -10.1776975, -10.0620689, -9.78825047, -9.34038497, -7.63854184, -5.34000051, -3.18258604, -1.3680886, -0.0982983435, 0.424994576, 0, -3.08931915, -7.32137822, -12.6888395, -19.1843653, -26.8006179, -35.5302595, -49.0355546, -70.427366, -90.9239181, -109.027897, -123.241988, -132.068878, -134.011252, -220.744836, -331.830711, -399.777169, -404.029121, -324.031477, -139.229147, 170.932957, 1189.77782, 1682.44204, 2145.68871, 3076.28091, 5592.05626, 8657.798, 13515.5527, 21407.3669, 32906.4806, 51261.3599, 89151.949, 120385.793, 124416.287, 122483.423, 118297.419, 109902.396, 95063.6149, 78327.8777, 61930.0868, 48105.1446}, 
{7.30722818, 6.78165834, 6.2222456, 5.63092945, 5.00964936, 4.36034482, 3.6849553, 2.98542031, 2.26367931, 1.52167178, 0.761337221, -0.0153848974, -0.80655509, -1.61023387, -2.44139153, -3.2873914, -4.13431683, -4.97676958, -5.80935141, -6.62666409, -7.42330939, -8.25801351, -9.05239776, -9.79220786, -10.4631895, -11.0510886, -11.5416507, -11.9223401, -12.1803631, -12.2973559, -12.2583665, -12.0484429, -11.652633, -11.0559848, -9.17566822, -6.67227475, -4.28722941, -2.22400572, -0.686077233, 0.123082519, 0, -2.92215445, -6.95377737, -12.0686214, -18.2404393, -25.4429835, -33.6500069, -46.3563417, -66.5277358, -85.7186384, -102.473131, -115.335294, -122.849209, -123.558957, -201.015913, -299.988234, -358.809399, -358.63595, -280.62443, -105.93138, 184.286659, 1171.02488, 1639.01505, 2052.88883, 2877.27786, 5016.97965, 7850.88622, 12371.824, 19572.6195, 29795.4396, 45985.0894, 79044.5536, 106928.722, 112598.551, 111337.789, 106287.883, 97341.2319, 84708.7183, 70576.344, 56456.772, 43862.6649}, 
{9.20514842, 8.63722588, 8.01374915, 7.3382654, 6.6143218, 5.8454655, 5.03524369, 4.18720353, 3.30489219, 2.39185683, 1.45164463, 0.487802757, -0.496121627, -1.49658135, -2.52715272, -3.5734562, -4.61981607, -5.65918269, -6.68450644, -7.68873768, -8.66482677, -9.67964643, -10.6420155, -11.5346751, -12.3403666, -13.0418309, -13.6218094, -14.0573956, -14.3265315, -14.4254627, -14.3392234, -14.0528478, -13.5513704, -12.8198253, -10.7519076, -8.02303678, -5.39124208, -3.06478227, -1.25191614, -0.160902459, 0, -2.89490659, -6.87163049, -11.873619, -17.8443194, -24.727179, -32.4656451, -44.1354082, -62.3415014, -79.5378213, -94.3960238, -105.587765, -111.784701, -111.658487, -179.78201, -266.811531, -317.657921, -315.467867, -243.388055, -84.5651731, 177.854093, 1055.49709, 1480.67478, 1860.11891, 2600.56127, 4472.5578, 6985.40964, 10982.4341, 17306.9484, 26267.2918, 40311.7159, 68486.5005, 92528.7394, 98412.2179, 97620.7655, 92766.4715, 84466.186, 73852.4758, 61990.28, 49967.5925, 38872.4072}, 
{11.9539984, 11.396919, 10.7419251, 9.9961151, 9.16658737, 8.26044025, 7.28477208, 6.24668123, 5.15326603, 4.01162485, 2.82885603, 1.61205793, 0.368328902, -0.895232707, -2.18239559, -3.47918686, -4.77081105, -6.04794703, -7.3012737, -8.52146993, -9.69921461, -10.888297, -12.0075697, -13.0389955, -13.9645373, -14.7661581, -15.4258207, -15.9226181, -16.2360743, -16.3550149, -16.262568, -15.9418619, -15.3760248, -14.548185, -12.3130815, -9.36766573, -6.49430298, -3.90692415, -1.81946011, -0.445841767, 0, -2.83192104, -6.72447656, -11.5965936, -17.367199, -23.95522, -31.2795835, -42.0007638, -58.3672847, -73.7427208, -86.9328428, -96.7434215, -101.980228, -101.449032, -160.189052, -235.106538, -278.094161, -274.410187, -209.312884, -68.0605165, 164.088647, 884.061076, 1259.43462, 1620.23501, 2296.48794, 3967.08945, 6101.52072, 9448.44848, 14756.5395, 22360.0286, 34250.8782, 58000.5016, 78152.0784, 82774.3965, 82064.3976, 78248.386, 71601.9274, 62799.6613, 52870.6706, 42723.2044, 33265.5121}, 
{15.6671109, 15.2502665, 14.6526298, 13.8888203, 12.9734576, 11.921161, 10.7465502, 9.46424457, 8.0888636, 6.63502681, 5.11735366, 3.55046366, 1.94897629, 0.327511038, -1.28246908, -2.87699986, -4.45339206, -6.0004714, -7.5070636, -8.96199437, -10.3540894, -11.7138662, -12.9827007, -14.1436609, -15.1798145, -16.0742294, -16.8099733, -17.3781596, -17.7686925, -17.9454023, -17.8880901, -17.576557, -16.9906041, -16.1100324, -13.7676545, -10.6688606, -7.60586307, -4.79307008, -2.4448898, -0.775730382, 0, -2.42321699, -5.8798827, -10.2956087, -15.5960064, -21.7066874, -28.5532632, -38.631874, -54.0524075, -68.5940473, -81.1386877, -90.568223, -95.7645474, -95.6095552, -144.935686, -207.230314, -241.745606, -235.919234, -177.18887, -52.9921887, 149.233139, 762.447176, 1093.2023, 1412.08142, 1989.66745, 3399.62048, 5164.93501, 7919.88365, 12298.739, 18650.8127, 28465.2603, 47623.0852, 63949.4113, 68060.681, 67420.7648, 63840.1625, 57976.2463, 51008.8496, 43181.0877, 35055.4132, 27194.2787}, 
{20.7028796, 20.5861279, 20.1539088, 19.4333935, 18.4517532, 17.2361592, 15.8137827, 14.2117952, 12.4573677, 10.5776715, 8.59987803, 6.55115842, 4.45868398, 2.34962596, 0.323808252, -1.63755819, -3.5587619, -5.42749212, -7.23143807, -8.95828898, -10.5957341, -12.1542692, -13.5956273, -14.9043477, -16.06497, -17.0620336, -17.8800781, -18.5192514, -18.9773562, -19.2016102, -19.1702154, -18.8613735, -18.2532866, -17.3241566, -14.8846365, -11.6515239, -8.42999473, -5.43706336, -2.88974419, -1.00505161, 0, -2.21275942, -5.44624247, -9.62451741, -14.6716525, -20.511716, -27.0687761, -36.3999129, -50.2416955, -63.4179179, -74.9865771, -84.00567, -89.5331937, -90.6271451, -131.051248, -181.582493, -209.613846, -205.059116, -157.832113, -57.8466478, 104.98347, 594.229003, 862.120608, 1121.81634, 1586.47426, 2713.72533, 4113.69602, 6292.90951, 9757.88901, 14790.8997, 22571.2387, 37853.8867, 50706.6617, 53429.424, 52491.7453, 49449.0084, 44759.2656, 39356.2846, 33368.3027, 27198.171, 21248.7405}, 
{26.2786764, 26.8149455, 26.8307674, 26.3716265, 25.4830075, 24.2103949, 22.5992733, 20.6951271, 18.5434409, 16.1896994, 13.679387, 11.0579882, 8.37098771, 5.66386998, 3.13542106, 0.734146418, -1.60305272, -3.86204835, -6.02871249, -8.08891713, -10.0285343, -11.8312131, -13.4853554, -14.9771402, -16.2927466, -17.4183534, -18.3401396, -19.061775, -19.584301, -19.8520745, -19.8441726, -19.5396726, -18.9176517, -17.9571871, -15.4599763, -12.1499111, -8.84000616, -5.75016203, -3.10027925, -1.11025839, 0, -2.18440409, -5.3852273, -9.51922572, -14.5031554, -20.2537726, -26.6878332, -35.1634238, -46.8220536, -58.1340934, -68.4310726, -77.0445209, -83.305968, -86.5469435, -119.067102, -159.342218, -183.352521, -183.686185, -152.931381, -83.6762811, 31.4909424, 378.75362, 568.315242, 756.471386, 1099.51763, 1936.39386, 2972.8742, 4584.42565, 7146.5152, 10821.9368, 16624.1765, 28601.6773, 38298.484, 39062.8219, 37650.5301, 35460.6162, 32288.2609, 28204.785, 23780.7761, 19416.6161, 15512.6874}, 
{30.7209326, 32.3181309, 33.1233752, 33.2056257, 32.6338428, 31.4769868, 29.804018, 27.6838966, 25.1855831, 22.3780376, 19.3302205, 16.1110921, 12.7896127, 9.43474254, 6.36402724, 3.48917116, 0.703061833, -1.97618648, -4.53045955, -6.94164312, -9.19162296, -11.192986, -13.0064874, -14.6235837, -16.0357312, -17.2343863, -18.2110054, -18.9885104, -19.5850951, -19.9169783, -19.9628406, -19.7013622, -19.1112237, -18.1711055, -15.6711163, -12.3418315, -9.00379556, -5.87880157, -3.18864237, -1.15511088, 0, -2.05100706, -5.1331807, -9.17747401, -14.1148401, -19.8762319, -26.3926027, -34.3899995, -44.4749973, -54.6773882, -64.6052914, -73.8668264, -82.0701122, -88.8232684, -115.820738, -148.698752, -171.421702, -179.080132, -166.764588, -129.565612, -62.5737494, 153.460398, 267.996353, 383.268615, 601.511685, 1178.28366, 1837.62182, 2865.74751, 4548.88207, 7038.95973, 11025.0632, 19420.5602, 26011.4796, 25852.1266, 24394.1761, 22732.7513, 20557.3096, 17767.5076, 14827.8707, 12010.3547, 9586.91506}, 
{33.2118217, 37.2674517, 40.093284, 41.7904447, 42.4600599, 42.2032557, 41.1211581, 39.3148933, 36.8855873, 33.9343662, 30.5623561, 26.870683, 22.9604731, 18.9328525, 15.2375346, 11.7551279, 8.33985421, 5.02153894, 1.83000748, -1.20491475, -4.05340236, -6.52948658, -8.78104985, -10.7998312, -12.5775698, -14.1060045, -15.3768746, -16.3943759, -17.1608327, -17.6381981, -17.8131533, -17.6723794, -17.2025577, -16.3903692, -14.1693483, -11.1924007, -8.19215066, -5.37069192, -2.93011813, -1.07252293, 0, -1.59714616, -4.1511867, -7.63185005, -12.0088647, -17.251959, -23.3308615, -30.4691984, -38.8524451, -47.8432309, -57.3081937, -67.1139712, -77.1272012, -87.2145214, -108.962979, -134.824838, -156.365179, -171.053327, -176.358606, -169.750341, -148.697858, -43.8984837, 2.44346809, 43.0311626, 130.567765, 428.429069, 716.326467, 1179.76921, 2004.26657, 3330.98416, 5478.46212, 9981.34654, 13433.8445, 13098.0778, 12052.9428, 10942.9317, 9628.42492, 8131.29836, 6611.55752, 5182.11805, 3955.89562}, 
{32.9833959, 41.9638812, 49.0019939, 54.2435725, 57.834455, 59.9204799, 60.6474854, 60.1613096, 58.6077911, 56.1327678, 52.8820783, 49.0015607, 44.6370532, 39.9343943, 35.4826326, 31.14723, 26.7600988, 22.3764224, 18.0513838, 13.8401663, 9.79795311, 6.2160221, 2.8808555, -0.184969778, -2.95887684, -5.41828877, -7.54062865, -9.21369604, -10.3387575, -11.1075362, -11.5338454, -11.6314981, -11.4143075, -10.8960868, -9.48107211, -7.56869726, -5.60453816, -3.72709184, -2.07485534, -0.786325716, 0, -0.605736375, -1.88690262, -3.87822766, -6.6144404, -10.1302698, -14.4604447, -19.5041537, -25.1809518, -31.8496605, -39.6000426, -48.5218608, -58.7048776, -70.2388559, -84.6194615, -101.047078, -118.615667, -137.127426, -156.384551, -176.18924, -196.343689, -174.409178, -186.787085, -223.983128, -276.503022, -332.503198, -388.284267, -431.532825, -449.935465, -441.679614, -362.949371, -45.2891724, 163.197081, -19.4728738, -202.899705, -322.640033, -440.529833, -632.727459, -791.822389, -872.984671, -831.384353}
};

  const double b4tab[21][81] = {
{-76.9329587, -74.9120768, -72.8831959, -70.8465369, -68.8023206, -66.750768, -64.6921, -62.6265373, -60.554301, -58.4756119, -56.3906909, -54.2997588, -52.2030367, -50.1007452, -47.9461158, -45.769095, -43.6031866, -41.4582228, -39.3440357, -37.2704577, -35.247321, -33.4197451, -31.6435908, -29.9100064, -28.2101399, -26.5351394, -24.8761532, -23.2108204, -21.5188101, -19.8235718, -18.1217387, -16.4099441, -14.6848212, -12.943003, -11.4512563, -10.0353269, -8.50051257, -6.78819309, -4.83974807, -2.59655716, 0, 2.76386468, 6.03212818, 9.89720271, 14.4515005, 19.7874337, 25.9974146, 34.7489177, 47.4727453, 60.495156, 73.2690363, 85.2472724, 95.8827507, 104.628358, 205.529346, 338.199126, 435.090883, 476.309414, 441.959521, 312.146002, 66.9736571, -1739.1593, -2089.86491, -1995.19665, -2465.208, -8373.90658, -11200.2837, -16591.1008, -30193.1191, -52541.1941, -104617.805, -252880.445, -373123.584, -369854.392, -339140.269, -301672.616, -256155.659, -206854.116, -157777.913, -113023.84, -76688.683}, 
{-61.0665593, -59.5405326, -57.9862446, -56.4062639, -54.803159, -53.1794984, -51.5378507, -49.8807845, -48.2108684, -46.5306708, -44.8427603, -43.1497056, -41.4540751, -39.7584374, -38.0534102, -36.3491224, -34.6566071, -32.9808773, -31.326946, -29.6998262, -28.1045309, -26.6183587, -25.1640537, -23.736646, -22.3311654, -20.9426417, -19.5661048, -18.1863601, -16.7897493, -15.3957502, -14.0035442, -12.6123126, -11.2212368, -9.82949807, -8.82586569, -7.96306628, -6.96433642, -5.74917049, -4.23706286, -2.3475079, 0, 2.5232123, 5.71548762, 9.70743045, 14.6296452, 20.6127364, 27.7873085, 38.0516951, 53.0386074, 68.6518075, 84.304145, 99.4084695, 113.377631, 125.624478, 248.967192, 411.077997, 531.044831, 585.084421, 549.413493, 400.248775, 113.806991, -2210.93178, -2654.26793, -2462.11744, -2880.39627, -8492.70369, -12312.0249, -19589.4638, -35576.1242, -61261.328, -114681.525, -254898.409, -369567.039, -372005.489, -350069.387, -321533.332, -283016.134, -234619.989, -183669.342, -135026.499, -93553.7637}, 
{-46.6827849, -44.8401848, -43.166788, -41.64566, -40.2598663, -38.9924727, -37.8265445, -36.7451475, -35.7313471, -34.768209, -33.8387986, -32.9261816, -32.0134235, -31.0835898, -30.014852, -28.8566321, -27.666292, -26.4483527, -25.207335, -23.9477595, -22.6741473, -21.4036985, -20.1265034, -18.8453316, -17.5629529, -16.2821369, -15.0056532, -13.72913, -12.4492681, -11.1859141, -9.94473733, -8.7314073, -7.55159344, -6.41096523, -5.88317934, -5.61459444, -5.20857505, -4.55461355, -3.5422023, -2.06083366, 0, 2.40324205, 5.65496541, 9.91367872, 15.3378906, 22.0861097, 30.3168448, 42.2279784, 59.7109529, 78.0479007, 96.569279, 114.605545, 131.487156, 146.54457, 288.992582, 475.996345, 614.889334, 678.435238, 639.397747, 470.54055, 144.627335, -2932.79777, -3509.49445, -3113.86505, -3274.31192, -7801.70177, -12594.3728, -21919.6626, -40044.909, -67533.3141, -119764.622, -248260.16, -353625.192, -358898.889, -343194.683, -321641.463, -289396.778, -243276.737, -193434.545, -145476.652, -105009.505}, 
{-31.8683333, -30.9374697, -30.0361184, -29.1612377, -28.3097862, -27.4787223, -26.6650044, -25.865591, -25.0774406, -24.2975115, -23.5227623, -22.7501514, -21.9766372, -21.1991781, -20.3947106, -19.5728592, -18.7447638, -17.9114784, -17.0740568, -16.2335527, -15.3910199, -14.5491341, -13.7071032, -12.8657571, -12.0259258, -11.188439, -10.3541265, -9.51942752, -8.68144093, -7.85149525, -7.03220309, -6.22617708, -5.43602981, -4.66437391, -4.52962731, -4.64484193, -4.57645592, -4.201124, -3.3955009, -2.03624133, 0, 2.38547115, 5.75691422, 10.2999741, 16.2002957, 23.643524, 32.8153038, 46.2870767, 66.2717898, 87.2503757, 108.439769, 129.056905, 148.318717, 165.442141, 326.167472, 537.020431, 693.43557, 764.659679, 719.939552, 528.521982, 159.653759, -2608.91273, -3214.82827, -3099.76368, -3705.3898, -9026.72777, -14207.2007, -23752.4752, -42168.218, -71709.0263, -123633.775, -239923.275, -335811.755, -345417.838, -331331.384, -307187.809, -273154.191, -232810.415, -189553.475, -146842.712, -108137.47}, 
{-20.0302572, -19.6381658, -19.2053708, -18.7358815, -18.2337072, -17.7028571, -17.1473405, -16.5711667, -15.978345, -15.3728846, -14.7587948, -14.1400848, -13.520764, -12.9048416, -12.3275698, -11.7731935, -11.2235926, -10.6763859, -10.1291922, -9.57963022, -9.02531876, -8.41389128, -7.79985521, -7.18773265, -6.58204571, -5.98731649, -5.40806709, -4.86092485, -4.36069817, -3.88296412, -3.42732968, -2.99340183, -2.58078756, -2.18909385, -2.53214569, -3.15773313, -3.5595831, -3.59121517, -3.10614896, -1.95790403, 0, 2.50229889, 6.16408247, 11.1886959, 17.7794844, 26.139793, 36.472967, 51.7037552, 74.3479831, 98.1018087, 122.0636, 145.331724, 167.00455, 186.180444, 359.580345, 586.590006, 754.66385, 830.659851, 781.435983, 573.850219, 174.760534, -2303.86234, -2917.28492, -3028.24439, -3999.47794, -9783.43757, -15354.9138, -25184.2768, -43741.8967, -73989.5527, -124923.387, -233695.183, -322248.095, -329244.738, -315251.626, -293243.744, -262504.614, -225610.988, -185658.156, -145564.631, -108248.927}, 
{-12.6902541, -12.0962501, -11.5645809, -11.0876955, -10.6580429, -10.2680723, -9.91023253, -9.5769728, -9.26074213, -8.95398956, -8.64916416, -8.33871498, -8.01509107, -7.67074149, -7.21941031, -6.70333566, -6.17071324, -5.63009056, -5.09001513, -4.55903445, -4.04569602, -3.6434085, -3.26413833, -2.90471312, -2.56196047, -2.23270796, -1.91378321, -1.57505508, -1.19044332, -0.821236911, -0.475209569, -0.160134995, 0.116213103, 0.34606102, -0.359519246, -1.45087998, -2.31167396, -2.76944208, -2.65172521, -1.78606422, 0, 2.55805603, 6.42574797, 11.8198493, 18.9571337, 28.0543746, 39.3283455, 56.4541114, 82.5870872, 109.674877, 136.530078, 161.965285, 184.793094, 203.826103, 387.349059, 626.965857, 801.452096, 874.956242, 811.626764, 575.612131, 131.060809, -2646.10171, -3379.52842, -3502.81122, -4449.542, -9401.98794, -15480.3531, -26216.6228, -45142.7826, -73801.6903, -121692.714, -225579.161, -309758.77, -314876.602, -300584.785, -279611.063, -250946.073, -217184.806, -180241.945, -142264.098, -105397.871}, 
{-8.67512051, -8.41984217, -8.11450304, -7.76504885, -7.37742535, -6.95757828, -6.51145338, -6.04499641, -5.56415309, -5.07486918, -4.58309042, -4.09476254, -3.6158313, -3.15224244, -2.75719541, -2.4067431, -2.07336145, -1.75333087, -1.44293172, -1.1384444, -0.836149293, -0.476906124, -0.120069905, 0.230425007, 0.570644256, 0.896653486, 1.20451834, 1.50199322, 1.79507616, 2.05188364, 2.26373529, 2.42195074, 2.51784962, 2.54275157, 1.50268045, 0.0122575862, -1.23931314, -2.059178, -2.25448326, -1.63237517, 0, 2.61669024, 6.65957361, 12.3517224, 19.916209, 29.5761055, 41.5544845, 60.1388366, 89.0059252, 118.660334, 147.67485, 174.622262, 198.075357, 216.606923, 400.555509, 639.834256, 811.354392, 878.555439, 804.876918, 553.758348, 88.6392516, -2073.4964, -2770.75087, -3221.23308, -4643.05199, -10416.5871, -16426.3604, -26485.1848, -44405.8732, -73070.8069, -119084.093, -211386.946, -287977.184, -300412.52, -290833.052, -269475.859, -239172.615, -206071.305, -170460.078, -134624.763, -100851.19}, 
{-6.56467003, -6.16903661, -5.75562568, -5.32659038, -4.88408385, -4.43025922, -3.96726963, -3.49726822, -3.02240813, -2.54484249, -2.06672443, -1.5902071, -1.11744363, -0.650587158, -0.199108732, 0.23946785, 0.668172719, 1.08634956, 1.49334205, 1.88849387, 2.2711487, 2.67469139, 3.05972314, 3.42088631, 3.75282326, 4.05017636, 4.30758797, 4.52566205, 4.70410679, 4.82330965, 4.87549239, 4.85287677, 4.74768455, 4.55213748, 3.1606537, 1.25992955, -0.380239742, -1.54308576, -2.01184009, -1.5697343, 0, 3.08044718, 7.69909034, 14.0497285, 22.3261608, 32.7221862, 45.4316038, 65.1716839, 95.9799915, 127.234189, 157.2914, 184.508749, 207.243361, 223.85236, 404.978382, 639.990109, 805.214965, 863.770499, 778.774256, 513.343783, 30.596626, -2021.49672, -2750.28202, -3299.07912, -4811.2079, -10199.2299, -16242.3468, -26206.9516, -43359.4371, -69736.4892, -112293.622, -198712.494, -271662.186, -286963.163, -278845.454, -256704.337, -225749.505, -195124.667, -162423.052, -129175.57, -96913.1306}, 
{-5.56124504, -5.05499527, -4.54410271, -4.02898833, -3.51007309, -2.98777796, -2.46252392, -1.93473193, -1.40482298, -0.873218021, -0.340338033, 0.193396014, 0.727563151, 1.26174241, 1.81049577, 2.36392399, 2.91099369, 3.44821924, 3.97211504, 4.47919548, 4.96597494, 5.4534668, 5.91030297, 6.32961438, 6.70453194, 7.02818655, 7.29370913, 7.50304539, 7.65681649, 7.72707617, 7.70337626, 7.57526857, 7.33230493, 6.96403715, 5.21864841, 2.87098488, 0.789858834, -0.781266194, -1.59892669, -1.41965913, 0, 3.51847295, 8.6793345, 15.6411182, 24.5623574, 35.6015858, 48.9173368, 69.6540832, 102.207101, 134.812973, 165.605781, 192.719607, 214.288534, 228.446643, 404.925127, 633.30511, 790.178568, 838.580815, 741.547159, 462.112912, -36.6866138, -2284.05963, -3089.80009, -3601.01998, -4964.8313, -9839.3587, -15644.5593, -25340.7457, -41888.2307, -67517.4983, -107378.347, -183709.23, -249142.889, -267254.155, -263181.338, -244282.62, -216004.911, -186655.952, -155011.43, -123096.606, -92936.7449}, 
{-5.27365126, -4.70375745, -4.12317741, -3.53283143, -2.93363981, -2.32652285, -1.71240084, -1.09219407, -0.466822831, 0.162792571, 0.795731847, 1.4310747, 2.06790084, 2.70528998, 3.35530245, 4.00880637, 4.65568726, 5.29236974, 5.91527844, 6.52083796, 7.10547294, 7.70873069, 8.27795759, 8.80362274, 9.27619522, 9.6861441, 10.0239385, 10.2953378, 10.5038039, 10.603245, 10.5779216, 10.4120945, 10.0900245, 9.59597211, 7.53027008, 4.75265703, 2.22762887, 0.222517333, -0.995345823, -1.15862885, 0, 3.83543963, 9.37725417, 16.7425749, 26.048533, 37.4122599, 50.9508867, 72.1330057, 105.623088, 138.74231, 169.43494, 195.645242, 215.317483, 226.395929, 393.567701, 609.294928, 753.418037, 789.775422, 682.205474, 394.546589, -109.362842, -2448.65246, -3295.332, -3774.37078, -5010.73815, -9531.90674, -15003.3467, -24233.018, -40028.8803, -65171.9598, -102551.017, -168472.081, -225525.938, -244736.343, -244336.067, -229668.948, -205209.544, -176713.735, -145681.396, -115063.321, -87810.301}, 
{-5.56216828, -4.89659609, -4.22510258, -3.54780906, -2.86483687, -2.17630733, -1.48234175, -0.783061461, -0.0785877887, 0.630957947, 1.34545442, 2.06478032, 2.78881431, 3.51743507, 4.26997829, 5.03401407, 5.7956397, 6.55075409, 7.29525615, 8.02504478, 8.7360189, 9.48666324, 10.2016473, 10.8682266, 11.4736563, 12.005192, 12.4500888, 12.8159221, 13.1072135, 13.2626317, 13.2611816, 13.081868, 12.7036959, 12.1056701, 9.77307445, 6.62984725, 3.71299499, 1.30705091, -0.303451771, -0.833979815, 0, 4.0926515, 9.89657631, 17.4740066, 26.8871747, 38.1983126, 51.4696528, 72.2842673, 105.395656, 137.665093, 166.91317, 190.960483, 207.627624, 214.735188, 366.046031, 560.731102, 686.157209, 708.248243, 592.928098, 306.120667, -186.250159, -2169.11327, -2959.58379, -3540.04735, -4892.88955, -9245.41646, -14509.2179, -23160.5721, -37675.7573, -60261.3122, -94202.7332, -155991.43, -208832.763, -224134.522, -222840.103, -210216.716, -188897.934, -162270.56, -133342.115, -105249.421, -81129.3001}, 
{-6.81819222, -6.11133366, -5.37916772, -4.62300956, -3.84417433, -3.04397719, -2.22373331, -1.38475784, -0.528365933, 0.34412724, 1.23140652, 2.13215676, 3.04506279, 3.96880946, 4.92891572, 5.90709105, 6.88301372, 7.84987985, 8.8008856, 9.72922711, 10.6281005, 11.5449591, 12.4112486, 13.2126719, 13.9349318, 14.5637312, 15.0847729, 15.5091816, 15.8442621, 16.0149308, 15.9965684, 15.7645555, 15.2942731, 14.5611019, 11.9296181, 8.39420326, 5.07116982, 2.26537526, 0.281677096, -0.575067161, 0, 4.52295874, 10.7820589, 18.7767731, 28.5065737, 39.9709333, 53.1693243, 73.8052084, 106.724954, 138.289145, 166.181252, 188.084742, 201.683088, 204.659757, 342.470831, 519.318287, 628.218366, 636.62887, 512.007604, 221.812371, -266.499026, -1911.56176, -2686.45663, -3397.68549, -4851.75019, -8873.56267, -13757.553, -21532.3054, -34226.4041, -52759.5347, -82486.9772, -144214.754, -195045.015, -201335.676, -197880.8, -190780.071, -176865.566, -152524.57, -125191.137, -98533.3628, -76219.3402}, 
{-8.50563467, -7.75712742, -6.95769135, -6.11052505, -5.21882713, -4.28579619, -3.31463085, -2.30852971, -1.27069137, -0.204314439, 0.887402479, 2.00126078, 3.13406185, 4.2826071, 5.46637267, 6.66781582, 7.86767744, 9.05812097, 10.2313099, 11.3794076, 12.4945777, 13.6589526, 14.7623014, 15.7843622, 16.7048731, 17.5035721, 18.1601975, 18.6635354, 19.0010128, 19.1307322, 19.0287581, 18.6711545, 18.0339858, 17.0933162, 14.2022937, 10.3636453, 6.70294558, 3.52821283, 1.14746542, -0.131278303, 0, 4.29175776, 10.2722567, 17.9121975, 27.182281, 38.0532079, 50.495679, 69.9784769, 101.144231, 130.817102, 156.735391, 176.637403, 188.261439, 189.345803, 312.838843, 470.944493, 565.631885, 566.983223, 445.080713, 170.006561, -288.157028, -1883.52972, -2621.96049, -3258.96568, -4550.06164, -7974.90776, -12502.8055, -19758.2359, -31365.6801, -47873.2169, -74134.5337, -128055.042, -173488.059, -182501.555, -180156.964, -171629.173, -156782.909, -136018.94, -112904.031, -89925.6229, -69571.1584}, 
{-11.1687064, -10.3817552, -9.50984147, -8.55865164, -7.53387206, -6.44118915, -5.28628928, -4.07485884, -2.81258425, -1.50515187, -0.158248107, 1.22244065, 2.631228, 4.06242757, 5.53237265, 7.02144711, 8.50836803, 9.98294507, 11.4349879, 12.8543061, 14.2307095, 15.6582985, 17.0081887, 18.2557864, 19.3764979, 20.3457297, 21.138888, 21.7312377, 22.0980647, 22.2151141, 22.0578497, 21.6017354, 20.8222351, 19.6948126, 16.5310094, 12.3588927, 8.33248143, 4.76758106, 1.97999707, 0.28553489, 0, 4.27317723, 10.1955007, 17.6913839, 26.6852403, 37.1014835, 48.8645268, 66.8206199, 95.0764473, 121.787751, 144.88052, 162.280745, 171.914416, 171.707522, 280.53568, 419.811704, 501.793056, 499.666481, 386.61872, 135.836519, -279.493379, -1700.41279, -2372.90976, -2959.46396, -4122.55508, -7123.68991, -11147.7518, -17576.0471, -27789.8824, -42285.6309, -65099.3987, -111099.244, -150308.439, -159708.128, -158151.438, -149949.174, -136152.92, -118673.915, -99244.0018, -79656.1808, -61703.4528}, 
{-15.0034096, -14.2724615, -13.3902868, -12.3680307, -11.2168388, -9.94785629, -8.57222864, -7.10110119, -5.54561934, -3.91692846, -2.22617391, -0.484501089, 1.29694464, 3.10701791, 4.94583587, 6.79512844, 8.63577297, 10.4543204, 12.2373218, 13.971328, 15.6428901, 17.328053, 18.911514, 20.3674642, 21.6700948, 22.793597, 23.712162, 24.4049746, 24.8504692, 25.0108968, 24.8584212, 24.3652061, 23.5034152, 22.2452121, 18.8383676, 14.3459026, 9.96136703, 6.00963498, 2.81558045, 0.704077447, 0, 4.20045209, 10.0199711, 17.3453243, 26.0632791, 36.0606028, 47.2240627, 63.77577, 89.2863988, 113.276414, 133.872294, 149.200516, 157.387561, 156.559905, 250.571061, 370.6642, 440.079422, 435.317615, 332.879666, 109.266463, -259.021108, -1425.77493, -2020.73533, -2582.30767, -3648.89731, -6332.28436, -9756.55346, -15148.1542, -23733.536, -36053.0643, -55391.4414, -94200.9791, -127087.449, -134450.543, -133052.038, -126571.412, -115494.953, -100976.185, -84692.8819, -68144.0964, -52828.88}, 
{-20.0462326, -19.591834, -18.8574763, -17.8657248, -16.639145, -15.200302, -13.5717613, -11.7760882, -9.8358481, -7.77360624, -5.611928, -3.37337873, -1.08052375, 1.24407159, 3.54735161, 5.82603928, 8.07916554, 10.2904016, 12.4434186, 14.5218878, 16.5094805, 18.4470624, 20.2532111, 21.9036989, 23.3742981, 24.6407807, 25.6789191, 26.486672, 27.0586643, 27.3176179, 27.2302965, 26.7634639, 25.8838841, 24.5583208, 20.9964813, 16.2766499, 11.6101121, 7.32302301, 3.74153784, 1.19181178, 0, 3.60418307, 8.78881258, 15.4502655, 23.4849188, 32.7891494, 43.2593342, 58.8619025, 82.9517093, 105.693172, 125.330096, 140.106285, 148.265546, 148.051682, 227.311307, 327.605162, 383.711452, 375.546468, 283.026502, 86.0678438, -235.413214, -1228.99871, -1754.62657, -2253.85336, -3168.23562, -5438.45292, -8274.89496, -12718.0611, -19808.4507, -30115.6316, -46092.8977, -77411.6011, -104063.578, -110627.87, -109379.491, -103312.972, -93537.7931, -82031.5276, -69184.4152, -55927.8351, -43193.1661}, 
{-26.6172289, -26.7224939, -26.3390842, -25.5085305, -24.272363, -22.6721124, -20.7493091, -18.5454836, -16.1021664, -13.460888, -10.6631789, -7.7505696, -4.76459053, -1.74677221, 1.14811707, 3.95018146, 6.69809673, 9.37349424, 11.9580054, 14.4332615, 16.7808939, 19.0100211, 21.0709912, 22.9416395, 24.5998012, 26.0233114, 27.1900054, 28.1121784, 28.7969475, 29.1397503, 29.1044301, 28.6548302, 27.754794, 26.3681648, 22.6649347, 17.7437439, 12.8390194, 8.2815219, 4.40201229, 1.53125138, 0, 3.30355828, 8.17110317, 14.4963509, 22.1730178, 31.0948201, 41.155474, 55.6343922, 77.3022458, 97.9571645, 116.118163, 130.304256, 139.034459, 140.827786, 205.93302, 287.492664, 333.21959, 326.961083, 252.564425, 93.876899, -165.25421, -957.708558, -1384.02209, -1791.95534, -2529.26882, -4346.7656, -6598.03722, -10114.4914, -15727.5357, -23898.7166, -36569.0259, -61556.1781, -82536.9717, -86855.432, -85156.2204, -80011.5735, -72196.0449, -63274.5668, -53448.24, -43383.1097, -33745.2212}, 
{-33.3425809, -34.5797467, -35.0142483, -34.7149696, -33.7507949, -32.1906082, -30.1032937, -27.5577354, -24.6228174, -21.3674239, -17.8604389, -14.1707467, -10.3672312, -6.51877659, -2.92507814, 0.49099202, 3.82403315, 7.05237166, 10.1543339, 13.1082464, 15.8924355, 18.472886, 20.8419704, 22.9797195, 24.8661641, 26.4813351, 27.8052633, 28.8560246, 29.6459781, 30.0641847, 30.0752278, 29.6436907, 28.7341567, 27.3112092, 23.5246129, 18.4870122, 13.4484249, 8.74464055, 4.71144862, 1.68463861, 0, 3.27487408, 8.11029377, 14.3888438, 21.9931088, 30.8056735, 40.7091226, 53.8891774, 72.1854853, 89.9735684, 106.200864, 119.81481, 129.762842, 134.992399, 187.316708, 252.240115, 291.265256, 292.52631, 244.157451, 134.292858, -48.9332942, -611.195059, -912.713132, -1208.30045, -1752.76996, -3101.57027, -4767.94823, -7367.47672, -11515.7286, -17477.6065, -26920.694, -46499.002, -62319.7035, -63452.8615, -61014.9319, -57316.892, -52034.8045, -45306.6837, -38057.6528, -30944.5013, -24624.0188}, 
{-37.4574626, -40.5138186, -42.3515766, -43.0745254, -42.7864534, -41.5911492, -39.5924013, -36.8939983, -33.5997288, -29.8133813, -25.6387443, -21.1796065, -16.5397563, -11.8229823, -7.50072163, -3.444187, 0.503238223, 4.31296495, 7.9564041, 11.4049666, 14.6300634, 17.4882773, 20.0817059, 22.3976185, 24.4232848, 26.1459741, 27.5529558, 28.6894543, 29.5919855, 30.1092422, 30.2049623, 29.8428838, 28.9867447, 27.6002828, 23.8125351, 18.747878, 13.6695669, 8.91661274, 4.82802626, 1.74281839, 0, 3.08152086, 7.74939156, 13.9085616, 21.4639804, 30.3205975, 40.3833624, 52.8391957, 68.6823865, 84.7525449, 100.434099, 115.111476, 128.169105, 138.991414, 182.3011, 235.127472, 271.82489, 284.54963, 265.45797, 206.706186, 100.450556, -248.727491, -430.444526, -610.947882, -956.48489, -1884.44788, -2941.59479, -4595.54455, -7313.91607, -11342.7906, -17814.4, -31506.4558, -42234.9855, -41893.8713, -39432.2948, -36652.8, -33051.8461, -28475.7232, -23677.4838, -19101.0555, -15190.3661}, 
{-37.6999664, -44.6073811, -49.6529908, -52.9865426, -54.757784, -55.1164621, -54.2123245, -52.1951182, -49.2145908, -45.4204894, -40.9625614, -35.9905541, -30.6542148, -25.1032908, -19.9894335, -15.1448838, -10.3638909, -5.69354162, -1.18092301, 3.12687795, 7.18277429, 10.6916503, 13.8887022, 16.7610974, 19.2960035, 21.480588, 23.3020182, 24.7764267, 25.9155933, 26.6374274, 26.9173356, 26.7307249, 26.053002, 24.8595739, 21.493007, 16.9616489, 12.3974458, 8.10978864, 4.40806809, 1.60167494, 0, 2.41668115, 6.31100838, 11.6465188, 18.3867495, 26.4952375, 35.9355201, 47.0895956, 60.2765839, 74.4594328, 89.4317708, 104.987226, 120.919427, 137.022002, 171.755359, 213.104453, 247.641248, 271.341234, 280.179905, 270.132754, 237.175274, 68.9326989, -3.94063173, -66.6118504, -204.248089, -682.503503, -1141.34506, -1882.52258, -3207.78589, -5345.72218, -8817.56921, -16127.4291, -21722.5291, -21138.8662, -19403.4913, -17570.8586, -15415.6481, -12978.8405, -10515.9876, -8209.49146, -6241.75435}, 
{-32.9009557, -47.1404618, -58.5117179, -67.2253254, -73.4918858, -77.5220004, -79.5262706, -79.7152978, -78.2996834, -75.4900288, -71.4969354, -66.5310046, -60.8028378, -54.5230364, -48.5130929, -42.5971569, -36.5531358, -30.4666784, -24.4234333, -18.5090492, -12.8091746, -7.78310138, -3.0912324, 1.23238631, 5.15370873, 8.63868887, 11.6532807, 14.0420896, 15.667955, 16.7869894, 17.4144182, 17.5654669, 17.255361, 16.4993258, 14.343992, 11.4173216, 8.42059515, 5.56715624, 3.07034838, 1.14351512, 0, 0.972835271, 3.00103908, 6.14331829, 10.4583798, 16.0049303, 22.8416769, 30.8168166, 39.8101783, 50.3838214, 62.6819266, 76.8486743, 93.0282452, 111.36482, 134.169525, 160.215722, 188.112747, 217.561549, 248.263077, 279.918281, 312.228111, 279.105702, 299.57109, 358.81613, 442.032678, 530.434992, 619.026309, 688.225405, 718.451057, 706.837009, 583.657133, 81.7638928, -247.651353, 41.7735923, 332.221774, 521.477762, 706.957138, 1007.80801, 1256.61738, 1383.66396, 1319.22649}
};
  
  /* Stuff for interpolating the NQC data */
  gsl_spline    *spline = NULL;
  gsl_interp_accel *acc = NULL;
  /* Interpolate the spin NQC data in 2-D parameter space -- spin and mass ratio */
  /* First, interpolate in spin dimension for all mass ratios */
  spline = gsl_spline_alloc( gsl_interp_cspline, adim );
  acc    = gsl_interp_accel_alloc();
  for (i = 0; i < qdim; i++)
  {
    gsl_spline_init( spline, alist, a3stab[i], adim );
    gsl_interp_accel_reset( acc );
    a3slist[i] = gsl_spline_eval( spline, a/(1.0-2.0*eta), acc );
    gsl_spline_init( spline, alist, a4tab[i], adim );
    gsl_interp_accel_reset( acc );
    a4list[i] = gsl_spline_eval( spline, a/(1.0-2.0*eta), acc );
    gsl_spline_init( spline, alist, a5tab[i], adim );
    gsl_interp_accel_reset( acc );
    a5list[i] = gsl_spline_eval( spline, a/(1.0-2.0*eta), acc );
    gsl_spline_init( spline, alist, b3tab[i], adim );
    gsl_interp_accel_reset( acc );
    b3list[i] = gsl_spline_eval( spline, a/(1.0-2.0*eta), acc );
    gsl_spline_init( spline, alist, b4tab[i], adim );
    gsl_interp_accel_reset( acc );
    b4list[i] = gsl_spline_eval( spline, a/(1.0-2.0*eta), acc );
//  printf("%.15f\n",a3slist[i]);
  }
//printf("%.15f\n",a);
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
  /* Second, interpolate in mass ratio dimension */
  spline = gsl_spline_alloc( gsl_interp_cspline, qdim );
  acc    = gsl_interp_accel_alloc();
  gsl_spline_init( spline, etalist, a3slist, qdim );
  gsl_interp_accel_reset( acc );
  coeffs->a3S = gsl_spline_eval( spline, eta, acc );
  gsl_spline_init( spline, etalist, a4list, qdim );
  gsl_interp_accel_reset( acc );
  coeffs->a4 = gsl_spline_eval( spline, eta, acc );
  gsl_spline_init( spline, etalist, a5list, qdim );
  gsl_interp_accel_reset( acc );
  coeffs->a5 = gsl_spline_eval( spline, eta, acc );
  gsl_spline_init( spline, etalist, b3list, qdim );
  gsl_interp_accel_reset( acc );
  coeffs->b3 = gsl_spline_eval( spline, eta, acc );
  gsl_spline_init( spline, etalist, b4list, qdim );
  gsl_interp_accel_reset( acc );
  coeffs->b4 = gsl_spline_eval( spline, eta, acc );
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
 
  /* Interpolate nonspin NQC data in the mass ratio dimension */
  spline = gsl_spline_alloc( gsl_interp_cspline, nsqdim );
  acc    = gsl_interp_accel_alloc();
  gsl_spline_init( spline, nsetalist, a1list, nsqdim );
  gsl_interp_accel_reset( acc );
  coeffs->a1 = gsl_spline_eval( spline, eta, acc );
  gsl_spline_init( spline, nsetalist, a2list, nsqdim );
  gsl_interp_accel_reset( acc );
  coeffs->a2 = gsl_spline_eval( spline, eta, acc );
  gsl_spline_init( spline, nsetalist, a3list, nsqdim );
  gsl_interp_accel_reset( acc );
  coeffs->a3 = gsl_spline_eval( spline, eta, acc );
  gsl_spline_init( spline, nsetalist, b1list, nsqdim );
  gsl_interp_accel_reset( acc );
  coeffs->b1 = gsl_spline_eval( spline, eta, acc );
  gsl_spline_init( spline, nsetalist, b2list, nsqdim );
  gsl_interp_accel_reset( acc );
  coeffs->b2 = gsl_spline_eval( spline, eta, acc );
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
  /* Andrea and I have different sign conventions, so I need to put a minus sign in front */
  coeffs->b1 = - coeffs->b1;
  coeffs->b2 = - coeffs->b2;
#if 0
  coeffs->a1 = -8.02798637014;
  coeffs->a2 = 48.7446843797;
  coeffs->a3 = -45.7900277224;
  coeffs->a3S= 0.;
  coeffs->a4 = 0.;
  coeffs->a5 = 0.;
  coeffs->b1 = 0.834742923041;
  coeffs->b2 = -2.33512320852; // q=1
#endif
#if 0
  coeffs->a1 = -7.79667;
  coeffs->a2 = 47.182;
  coeffs->a3 = 2238.85334023;
  coeffs->a3S= 0.;
  coeffs->a4 = -7143.16738899;
  coeffs->a5 = 5596.0086893;
  coeffs->b1 = 0.85069;
  coeffs->b2 = -2.47071; // q=1, chi1=chi2=0.98
#endif
#if 0
  coeffs->a1 = -6.82562365707;
  coeffs->a2 = 41.5711482044;
  coeffs->a3 = -39.4329799959;
  coeffs->a3S= 0.;
  coeffs->a4 = 0.;
  coeffs->a5 = 0.;
  coeffs->b1 = 0.461925688819;
  coeffs->b2 = -1.38733263299; // q=8
#endif
#if 0
  coeffs->a1 = -7.5758;
  coeffs->a2 = 46.9323;
  coeffs->a3 = -118.368375152;
  //coeffs->a3 = -45.0036; // NS part 
  coeffs->a3S= 0.;
  coeffs->a4 = 125.555824111;
  coeffs->a5 = -22.0751068073;
  //coeffs->a4 = 0.;
  //coeffs->a5 = 0.;
  coeffs->b1 = 0.51305;
  coeffs->b2 = -1.55133; // q=8, chi1=0.5
#endif
#if 0
  coeffs->a1 = -5.79216;
  coeffs->a2 = 30.3902;
  coeffs->a3 = 3329.40898295;
  coeffs->a3S= 0.0;
  coeffs->a4 = -11082.7670039;
  coeffs->a5 = 9152.28860893;
  coeffs->b1 = 0.974779;
  coeffs->b2 = -3.34626;
#endif
   
  /* Obsolete polynomial fitting of nonspin NQC coefficients a1, a2, a3, b1 and b2 */
  /*
  coeffs->a1 = -12.67955358602124 + 75.41927959573084 * eta - 106.15933052937714 * eta2;
  coeffs->a2 = 101.45522216901628 - 757.3158549733314 * eta + 1473.314771676588 * eta2;
  coeffs->a3 = -107.6647834845902 + 857.6219519536213 * eta - 1776.2776804623143 * eta2;
  // Andrea and I have different sign conventions, so I need to put a minus sign in front 
  coeffs->b1 = - (-1.464129495621165 + 12.81732978488213 * eta - 60.09957767247623 * eta2);
  coeffs->b2 = - ( 7.477426352542122 - 85.26122117590637 * eta + 353.3251639728075 * eta2);
  */

  return XLAL_SUCCESS;

}

/**
 * This function computes the coefficients a3s, a4, etc. used in the
 * non-quasicircular correction. The details of the calculation of these
 * coefficients are found in the DCC document T1100433.
 * In brief, this function populates and solves the linear equations
 * Eq. 18 (for amplitude) and Eq. 19 (for phase) of the DCC document T1100433v2.
 */
UNUSED static int XLALSimIMRSpinEOBCalculateNQCCoefficients(
                 REAL8Vector    * restrict amplitude,   /**<< Waveform amplitude, func of time */
                 REAL8Vector    * restrict phase,       /**<< Waveform phase(rad), func of time */
                 REAL8Vector    * restrict rVec,        /**<< Position-vector, function of time */
                 REAL8Vector    * restrict prVec,       /**<< Momentum vector, function of time */
                 REAL8Vector    * restrict orbOmegaVec, /**<< Orbital frequency, func of time */
                 INT4                      l,           /**<< Mode index l */
                 INT4                      m,           /**<< Mode index m */
                 REAL8                     timePeak,    /**<< Time of peak orbital frequency */
                 REAL8                     deltaT,      /**<< Sampling interval */
                 REAL8                     eta,         /**<< Symmetric mass ratio */
                 REAL8                     a,           /**<< Normalized spin of deformed-Kerr */
                 EOBNonQCCoeffs * restrict coeffs,      /**<< OUTPUT, NQC coefficients */
                 UINT4                     SpinAlignedEOBversion  /**<< 1 for SEOBNRv1, 2 for SEOBNRv2 */
)
{

  /* For gsl permutation stuff */

  int signum;

  REAL8Vector * restrict timeVec = NULL;

  /* Vectors which are used in the computation of the NQC coefficients */
  REAL8Vector *q3 = NULL, *q4 = NULL, *q5 = NULL;
  REAL8Vector *p3 = NULL, *p4 = NULL;

  REAL8Vector *qNS = NULL, *pNS = NULL;

  /* Since the vectors we actually want are q etc * A, we will have to generate them here */
  REAL8Vector *q3LM  = NULL;
  REAL8Vector *q4LM  = NULL;
  REAL8Vector *q5LM  = NULL; 
  REAL8Vector *qNSLM = NULL;

  REAL8 amp, aDot, aDDot;
  REAL8 omega, omegaDot;

  REAL8 qNSLMPeak, qNSLMDot, qNSLMDDot;
  REAL8 pNSLMDot, pNSLMDDot;

  REAL8 nra, nraDDot;
  REAL8 nromega, nromegaDot;

  REAL8 nrDeltaT, nrTimePeak;

  /* Stuff for finding numerical derivatives */
  gsl_spline    *spline = NULL;
  gsl_interp_accel *acc = NULL;

  /* Matrix stuff for calculating coefficients */
  gsl_matrix *qMatrix = NULL, *pMatrix = NULL;
  gsl_vector *aCoeff  = NULL, *bCoeff  = NULL;

  gsl_vector *amps = NULL, *omegaVec = NULL;

  gsl_permutation *perm1 = NULL, *perm2 = NULL;

  memset( coeffs, 0, sizeof( EOBNonQCCoeffs ) );

  /* Populate the time vector */
  /* It is okay to assume initial t = 0 */
  timeVec = XLALCreateREAL8Vector( rVec->length );
  q3    = XLALCreateREAL8Vector( rVec->length );
  q4    = XLALCreateREAL8Vector( rVec->length );
  q5    = XLALCreateREAL8Vector( rVec->length );
  p3    = XLALCreateREAL8Vector( rVec->length );
  p4    = XLALCreateREAL8Vector( rVec->length );
  qNS   = XLALCreateREAL8Vector( rVec->length );
  pNS   = XLALCreateREAL8Vector( rVec->length );
  q3LM  = XLALCreateREAL8Vector( rVec->length );
  q4LM  = XLALCreateREAL8Vector( rVec->length );
  q5LM  = XLALCreateREAL8Vector( rVec->length );
  qNSLM = XLALCreateREAL8Vector( rVec->length );

  if ( !timeVec || !q3 || !q4 || !q5 || !p3 || !p4 || !qNS || !pNS || !q3LM
          || !q4LM || !q5LM || !qNSLM )
  {
    XLALDestroyREAL8Vector( timeVec );
    XLALDestroyREAL8Vector( q3 );
    XLALDestroyREAL8Vector( q4 );
    XLALDestroyREAL8Vector( q5 );
    XLALDestroyREAL8Vector( p3 );
    XLALDestroyREAL8Vector( p4 );
    XLALDestroyREAL8Vector( qNS );
    XLALDestroyREAL8Vector( pNS );
    XLALDestroyREAL8Vector( q3LM );
    XLALDestroyREAL8Vector( q4LM );
    XLALDestroyREAL8Vector( q5LM );
    XLALDestroyREAL8Vector( qNSLM );
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* We need the calibrated non-spinning NQC coefficients */
  switch ( SpinAlignedEOBversion)
  {
   case 1:
     if ( XLALSimIMRGetEOBCalibratedSpinNQC( coeffs, l, m, eta, a ) == XLAL_FAILURE )
     {
       XLALDestroyREAL8Vector( timeVec );
       XLALDestroyREAL8Vector( q3 );
       XLALDestroyREAL8Vector( q4 );
       XLALDestroyREAL8Vector( q5 );
       XLALDestroyREAL8Vector( p3 );
       XLALDestroyREAL8Vector( p4 );
       XLALDestroyREAL8Vector( qNS );
       XLALDestroyREAL8Vector( pNS );
       XLALDestroyREAL8Vector( q3LM );
       XLALDestroyREAL8Vector( q4LM );
       XLALDestroyREAL8Vector( q5LM );
       XLALDestroyREAL8Vector( qNSLM );
       XLAL_ERROR( XLAL_EFUNC );
     }
     break;
   case 2:
     if ( XLALSimIMRGetEOBCalibratedSpinNQCv2( coeffs, l, m, eta, a ) == XLAL_FAILURE )
     {
       XLALDestroyREAL8Vector( timeVec );
       XLALDestroyREAL8Vector( q3 );
       XLALDestroyREAL8Vector( q4 );
       XLALDestroyREAL8Vector( q5 );
       XLALDestroyREAL8Vector( p3 );
       XLALDestroyREAL8Vector( p4 );
       XLALDestroyREAL8Vector( qNS );
       XLALDestroyREAL8Vector( pNS );
       XLALDestroyREAL8Vector( q3LM );
       XLALDestroyREAL8Vector( q4LM );
       XLALDestroyREAL8Vector( q5LM );
       XLALDestroyREAL8Vector( qNSLM );
       XLAL_ERROR( XLAL_EFUNC );
     }
     break;
   default:
     XLALPrintError( "XLAL Error - %s: Unknown SEOBNR version!\nAt present only v1 and v2 are available.\n", __func__);
     XLAL_ERROR( XLAL_EINVAL );
     break;
  }

  /* Populate vectors as necessary. Eqs. 14 - 17 of the LIGO DCC document T1100433v2 */
  for ( unsigned int i = 0; i < timeVec->length; i++ )
  {
    
    REAL8 rootR  = sqrt(rVec->data[i]);
    REAL8 rOmega = rVec->data[i] * orbOmegaVec->data[i];

    /* We don't need these as vectors as their coefficients are calibrated */
    REAL8 q1, q2, p1, p2;

    timeVec->data[i] = i * deltaT;
    q1            = prVec->data[i]*prVec->data[i] / (rOmega*rOmega);
    q2            = q1 / rVec->data[i];
    q3->data[i]   = q2 / rootR;
    q4->data[i]   = q2 / rVec->data[i];
    q5->data[i]   = q3->data[i] / rVec->data[i];

    p1          = prVec->data[i] / rOmega;
    p2          = p1 * prVec->data[i] * prVec->data[i];
    p3->data[i] = p2 / rootR;
    p4->data[i] = p2 / rVec->data[i];

    qNS->data[i]  = coeffs->a1 * q1 + coeffs->a2 * q2 + coeffs->a3 * q3->data[i];
    pNS->data[i]  = coeffs->b1 * p1 + coeffs->b2 * p2;
    q3LM->data[i] = q3->data[i] * amplitude->data[i];
    q4LM->data[i] = q4->data[i] * amplitude->data[i];
    q5LM->data[i] = q5->data[i] * amplitude->data[i];

    qNSLM->data[i] = qNS->data[i] * amplitude->data[i];
  }
  /* Allocate all the memory we need */
  XLAL_CALLGSL(
    /* a stuff */
    qMatrix = gsl_matrix_alloc( 3, 3 );
    aCoeff  = gsl_vector_alloc( 3 );
    amps    = gsl_vector_alloc( 3 );
    perm1   = gsl_permutation_alloc( 3 );

    /* b stuff */
    pMatrix  = gsl_matrix_alloc( 2, 2 );
    bCoeff   = gsl_vector_alloc( 2 );
    omegaVec = gsl_vector_alloc( 2 );
    perm2    = gsl_permutation_alloc( 2 );
  );

  if ( !qMatrix || !aCoeff || !amps || !pMatrix || !bCoeff || !omegaVec )
  {
    XLALDestroyREAL8Vector( timeVec );
    XLALDestroyREAL8Vector( q3 );
    XLALDestroyREAL8Vector( q4 );
    XLALDestroyREAL8Vector( q5 );
    XLALDestroyREAL8Vector( p3 );
    XLALDestroyREAL8Vector( p4 );
    XLALDestroyREAL8Vector( qNS );
    XLALDestroyREAL8Vector( pNS );
    XLALDestroyREAL8Vector( q3LM );
    XLALDestroyREAL8Vector( q4LM );
    XLALDestroyREAL8Vector( q5LM );
    XLALDestroyREAL8Vector( qNSLM );
    XLAL_ERROR( XLAL_ENOMEM );
  }

  /* The time we want to take as the peak time depends on l and m */
  /* Calculate the adjustment we need to make here */
  nrDeltaT   = XLALSimIMREOBGetNRSpinPeakDeltaT( l, m, eta, a );
  if ( XLAL_IS_REAL8_FAIL_NAN( nrDeltaT ) )
  {
    XLALDestroyREAL8Vector( timeVec );
    XLALDestroyREAL8Vector( q3 );
    XLALDestroyREAL8Vector( q4 );
    XLALDestroyREAL8Vector( q5 );
    XLALDestroyREAL8Vector( p3 );
    XLALDestroyREAL8Vector( p4 );
    XLALDestroyREAL8Vector( qNS );
    XLALDestroyREAL8Vector( pNS );
    XLALDestroyREAL8Vector( q3LM );
    XLALDestroyREAL8Vector( q4LM );
    XLALDestroyREAL8Vector( q5LM );
    XLALDestroyREAL8Vector( qNSLM );
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* nrDeltaT defined in XLALSimIMREOBGetNRSpinPeakDeltaT is a minus sign different from Eq. (33) of Taracchini et al.
   * Therefore, the plus sign in Eq. (21) of Taracchini et al and Eq. (18) of DCC document T1100433v2 is 
   * changed to a minus sign here.
   */
  nrTimePeak = timePeak - nrDeltaT;

  /* We are now in a position to use the interp stuff to calculate the derivatives we need */
  /* We will start with the quantities used in the calculation of the a coefficients */
  spline = gsl_spline_alloc( gsl_interp_cspline, amplitude->length );
  acc    = gsl_interp_accel_alloc();

  /* Populate the Q matrix in Eq. 18 of the LIGO DCC document T1100433v2 */
  /* Q3 */
  gsl_spline_init( spline, timeVec->data, q3LM->data, q3LM->length );
  gsl_matrix_set( qMatrix, 0, 0, gsl_spline_eval( spline, nrTimePeak, acc ) );
  gsl_matrix_set( qMatrix, 1, 0, gsl_spline_eval_deriv( spline, nrTimePeak, acc ) );
  gsl_matrix_set( qMatrix, 2, 0, gsl_spline_eval_deriv2( spline, nrTimePeak, acc ) );

  /* Q4 */
  gsl_spline_init( spline, timeVec->data, q4LM->data, q4LM->length );
  gsl_interp_accel_reset( acc );
  gsl_matrix_set( qMatrix, 0, 1, gsl_spline_eval( spline, nrTimePeak, acc ) );
  gsl_matrix_set( qMatrix, 1, 1, gsl_spline_eval_deriv( spline, nrTimePeak, acc ) );
  gsl_matrix_set( qMatrix, 2, 1, gsl_spline_eval_deriv2( spline, nrTimePeak, acc ) );

  /* Q5 */
  gsl_spline_init( spline, timeVec->data, q5LM->data, q5LM->length );
  gsl_interp_accel_reset( acc );
  gsl_matrix_set( qMatrix, 0, 2, gsl_spline_eval( spline, nrTimePeak, acc ) );
  gsl_matrix_set( qMatrix, 1, 2, gsl_spline_eval_deriv( spline, nrTimePeak, acc ) );
  gsl_matrix_set( qMatrix, 2, 2, gsl_spline_eval_deriv2( spline, nrTimePeak, acc ) );

  /* Populate the r.h.s vector of Eq. 18 of the LIGO DCC document T1100433v2 */
  /* Amplitude */
  gsl_spline_init( spline, timeVec->data, amplitude->data, amplitude->length );
  gsl_interp_accel_reset( acc );
  amp   = gsl_spline_eval( spline, nrTimePeak, acc );
  aDot  = gsl_spline_eval_deriv( spline, nrTimePeak, acc );
  aDDot = gsl_spline_eval_deriv2( spline, nrTimePeak, acc );

  /* qNSLM */
  gsl_spline_init( spline, timeVec->data, qNSLM->data, qNSLM->length );
  gsl_interp_accel_reset( acc );
  qNSLMPeak = gsl_spline_eval( spline, nrTimePeak, acc );
  qNSLMDot  = gsl_spline_eval_deriv( spline, nrTimePeak, acc );
  qNSLMDDot = gsl_spline_eval_deriv2( spline, nrTimePeak, acc );

  nra = GetNRSpinPeakAmplitude( l, m, eta, a );
  nraDDot = - GetNRSpinPeakADDot( l, m, eta, a );

  if ( XLAL_IS_REAL8_FAIL_NAN( nra ) || XLAL_IS_REAL8_FAIL_NAN( nraDDot ) )
  {
    XLALDestroyREAL8Vector( timeVec );
    XLALDestroyREAL8Vector( q3 );
    XLALDestroyREAL8Vector( q4 );
    XLALDestroyREAL8Vector( q5 );
    XLALDestroyREAL8Vector( p3 );
    XLALDestroyREAL8Vector( p4 );
    XLALDestroyREAL8Vector( qNS );
    XLALDestroyREAL8Vector( pNS );
    XLALDestroyREAL8Vector( q3LM );
    XLALDestroyREAL8Vector( q4LM );
    XLALDestroyREAL8Vector( q5LM );
    XLALDestroyREAL8Vector( qNSLM );
    XLAL_ERROR( XLAL_EFUNC );
  }

  gsl_vector_set( amps, 0, nra - amp - qNSLMPeak );
  gsl_vector_set( amps, 1, - aDot - qNSLMDot );
  gsl_vector_set( amps, 2, nraDDot - aDDot - qNSLMDDot );

  /* We have now set up all the stuff to calculate the a coefficients */
  /* So let us do it! */
  gsl_linalg_LU_decomp( qMatrix, perm1, &signum );
  gsl_linalg_LU_solve( qMatrix, perm1, amps, aCoeff );

  /* Now we (should) have calculated the a values. Now we can do the b values */

  /* Populate the P matrix in Eq. 18 of the LIGO DCC document T1100433v2 */
  /* P3 */
  gsl_spline_init( spline, timeVec->data, p3->data, p3->length );
  gsl_interp_accel_reset( acc );
  gsl_matrix_set( pMatrix, 0, 0, - gsl_spline_eval_deriv( spline, nrTimePeak, acc ) );
  gsl_matrix_set( pMatrix, 1, 0, - gsl_spline_eval_deriv2( spline, nrTimePeak, acc ) );

  /* P4 */
  gsl_spline_init( spline, timeVec->data, p4->data, p4->length );
  gsl_interp_accel_reset( acc );
  gsl_matrix_set( pMatrix, 0, 1, - gsl_spline_eval_deriv( spline, nrTimePeak, acc ) );
  gsl_matrix_set( pMatrix, 1, 1, - gsl_spline_eval_deriv2( spline, nrTimePeak, acc ) );

  /* Populate the r.h.s vector of Eq. 18 of the LIGO DCC document T1100433v2 */
  /* Phase */
  gsl_spline_init( spline, timeVec->data, phase->data, phase->length );
  gsl_interp_accel_reset( acc );
  omega    = gsl_spline_eval_deriv( spline, nrTimePeak, acc );
  omegaDot = gsl_spline_eval_deriv2( spline, nrTimePeak, acc );

  /* pNSLM */
  gsl_spline_init( spline, timeVec->data, pNS->data, pNS->length );
  gsl_interp_accel_reset( acc );
  pNSLMDot  = gsl_spline_eval_deriv( spline, nrTimePeak, acc );
  pNSLMDDot = gsl_spline_eval_deriv2( spline, nrTimePeak, acc );

  /* Since the phase can be decreasing, we need to take care not to have a -ve frequency */
  if ( omega * omegaDot > 0.0 )
  {
    omega    = fabs( omega );
    omegaDot = fabs( omegaDot );
  }
  else
  {
    omega    = fabs( omega );
    omegaDot = - fabs( omegaDot );
  }

  //nromega = GetNRPeakOmega( l, m, eta );
  //nromegaDot = GetNRPeakOmegaDot( l, m, eta );
  switch ( SpinAlignedEOBversion )
  {
   case 1:
     nromega = GetNRSpinPeakOmega( l, m, eta, a );
     nromegaDot = GetNRSpinPeakOmegaDot( l, m, eta, a );
     break;
   case 2:
     nromega = GetNRSpinPeakOmegav2( l, m, eta, a );
     nromegaDot = GetNRSpinPeakOmegaDotv2( l, m, eta, a );
     break;
   default:
     XLALPrintError( "XLAL Error - %s: Unknown SEOBNR version!\nAt present only v1 and v2 are available.\n", __func__);
     XLAL_ERROR( XLAL_EINVAL );
     break;
  }

  /*printf("NR inputs: %.16e, %.16e, %.16e, %.16e\n",nra,nraDDot,nromega,nromegaDot);
  printf("NR inputs: %.16e, %.16e, %.16e, %.16e\n",pNSLMDot, pNSLMDDot,omega,omegaDot);*/

  if ( XLAL_IS_REAL8_FAIL_NAN( nromega ) || XLAL_IS_REAL8_FAIL_NAN( nromegaDot ) )
  {
    XLALDestroyREAL8Vector( timeVec );
    XLALDestroyREAL8Vector( q3 );
    XLALDestroyREAL8Vector( q4 );
    XLALDestroyREAL8Vector( q5 );
    XLALDestroyREAL8Vector( p3 );
    XLALDestroyREAL8Vector( p4 );
    XLALDestroyREAL8Vector( qNS );
    XLALDestroyREAL8Vector( pNS );
    XLALDestroyREAL8Vector( q3LM );
    XLALDestroyREAL8Vector( q4LM );
    XLALDestroyREAL8Vector( q5LM );
    XLALDestroyREAL8Vector( qNSLM );
    XLAL_ERROR( XLAL_EFUNC );
  }

  gsl_vector_set( omegaVec, 0, nromega - omega + pNSLMDot );
  gsl_vector_set( omegaVec, 1, nromegaDot - omegaDot + pNSLMDDot );

  /*printf( "P MATRIX\n" );
  for (unsigned int i = 0; i < 2; i++ )
  {
    for (unsigned int j = 0; j < 2; j++ )
    {
      printf( "%.12e\t", gsl_matrix_get( pMatrix, i, j ));
    }
    printf( "= %.12e\n", gsl_vector_get( omegaVec, i ) );
  }*/

  /* And now solve for the b coefficients */
  gsl_linalg_LU_decomp( pMatrix, perm2, &signum );
  gsl_linalg_LU_solve( pMatrix, perm2, omegaVec, bCoeff );

  /* We can now populate the coefficients structure */
/*  coeffs->a3S = gsl_vector_get( aCoeff, 0 );
  coeffs->a4  = gsl_vector_get( aCoeff, 1 );
  coeffs->a5  = gsl_vector_get( aCoeff, 2 );*/
  switch ( SpinAlignedEOBversion )
  {
   case 1:
     coeffs->b3  = gsl_vector_get( bCoeff, 0 );
     coeffs->b4  = gsl_vector_get( bCoeff, 1 );
     break;
   case 2:
     break;
   default:
     XLALPrintError( "XLAL Error - %s: Unknown SEOBNR version!\nAt present only v1 and v2 are available.\n", __func__);
     XLAL_ERROR( XLAL_EINVAL );
     break;
  }
  coeffs->b3  *= 1.0;
  coeffs->b4  *= 1.0;
//  coeffs->b3  = -876.669217307; 
//  coeffs->b4  = 1386.13223658;
//  coeffs->b3 = 41583.9402122;
//  coeffs->b4 = 68359.70064;

  /*printf( "NQC coefficients:\n" );
  printf( "a1 = %.16e, a2 = %.16e, a3 = %.16e, a3s = %.16e, a4 = %.16e, a5 = %.16e\n",
    coeffs->a1, coeffs->a2, coeffs->a3, coeffs->a3S, coeffs->a4, coeffs->a5 );

  printf( "b1 = %.16e, b2 = %.16e, b3 = %.16e, b4 = %.16e\n",
    coeffs->b1, coeffs->b2, coeffs->b3, coeffs->b4 );*/

  /* Free memory and exit */
  gsl_matrix_free( qMatrix );
  gsl_vector_free( amps );
  gsl_vector_free( aCoeff );
  gsl_permutation_free( perm1 );

  gsl_matrix_free( pMatrix );
  gsl_vector_free( omegaVec );
  gsl_vector_free( bCoeff );
  gsl_permutation_free( perm2 );

  gsl_spline_free( spline );
  gsl_interp_accel_free( acc );

  XLALDestroyREAL8Vector( timeVec );
  XLALDestroyREAL8Vector( q3 );
  XLALDestroyREAL8Vector( q4 );
  XLALDestroyREAL8Vector( q5 );
  XLALDestroyREAL8Vector( p3 );
  XLALDestroyREAL8Vector( p4 );
  XLALDestroyREAL8Vector( qNS );
  XLALDestroyREAL8Vector( pNS );
  XLALDestroyREAL8Vector( q3LM );
  XLALDestroyREAL8Vector( q4LM );
  XLALDestroyREAL8Vector( q5LM );
  XLALDestroyREAL8Vector( qNSLM );

  return XLAL_SUCCESS;
}

#endif /*_LALSIMIMRNQCCORRECTION_C*/
