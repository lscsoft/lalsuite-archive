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
{194.03406207, 193.315403054, 192.340696042, 191.108131058, 189.615898125, 187.862187267, 185.845188507, 183.563091869, 181.014087376, 178.196365052, 175.108114921, 171.747527005, 168.112791328, 164.202097915, 159.881229844, 155.232138256, 150.346796911, 145.250478579, 139.96845603, 134.526002032, 128.948389356, 123.45137525, 117.843440792, 112.123551541, 106.290673055, 100.34377089, 94.2818106058, 88.0767763505, 81.7047065536, 75.2290824016, 68.6598247712, 62.006854539, 55.2800925814, 48.4894597752, 41.7952744946, 35.1123156102, 28.3439775648, 21.4694186698, 14.4677972364, 7.31827157602, 0., -7.06413601656, -14.3999825353, -22.089662454, -30.2152986705, -38.8590140827, -48.1029315883, -57.9430958504, -68.3884858166, -79.7270472756, -92.0758536549, -105.551978382, -120.272494883, -136.354476586, -161.225661949, -190.378371019, -218.751538717, -244.966902356, -267.646199249, -285.411166709, -296.883542052, -310.354741537, -328.370360469, -319.492151805, -252.281868503, 29.4084482122, 249.400554575, 588.783154503, 1228.64495192, 2363.3408566, 4134.16095488, 7287.44826007, 10172.075365, 11510.5971716, 12299.3266332, 12652.4651786, 12433.2946947, 11460.0625525, 9982.07372846, 8101.68110164, 5921.2375509}, 
{160.112547782, 160.667985093, 160.855095784, 160.683878451, 160.164331691, 159.3064541, 158.120244275, 156.615700811, 154.802822306, 152.691607356, 150.292054557, 147.614162506, 144.6679298, 141.463355034, 137.974809288, 134.234829257, 130.278648479, 126.123552866, 121.786828327, 117.285760773, 112.637636115, 107.962912061, 103.161453992, 98.2362990855, 93.1904845213, 88.0270474773, 82.7490251321, 77.3472031958, 71.8142083082, 66.1823723079, 60.4597068632, 54.6542236423, 48.7739343137, 42.8268505456, 36.9526253384, 31.0759935733, 25.1117254627, 19.0409065284, 12.844622292, 6.50395827526, 0., -6.3186181821, -12.889120591, -19.7811827168, -27.0644800495, -34.8086880791, -43.0834822956, -52.0954881413, -62.0307526891, -72.6314884523, -83.9117648395, -95.8856512595, -108.567217121, -121.970531833, -144.794468742, -171.559058589, -196.126930728, -216.735749697, -231.623180032, -239.026886272, -237.184532955, -230.598764417, -228.655499265, -197.207571435, -102.107814864, 154.463950585, 449.595197101, 893.840096377, 1597.75282011, 2688.58018497, 4226.79842773, 6614.56538705, 8915.92747576, 10410.9574341, 11588.5513029, 12442.701159, 12813.2516669, 12518.0513156, 11617.0807085, 10087.7750498, 7907.5695438}, 
{135.758719821, 136.138892686, 136.299306582, 136.238611887, 135.955458979, 135.448498234, 134.716380031, 133.757754747, 132.571272759, 131.155584444, 129.509340181, 127.631190347, 125.519785318, 123.173775473, 120.539911026, 117.649674588, 114.53847738, 111.215585531, 107.690265168, 103.971782416, 100.069403403, 95.9633652912, 91.6959722706, 87.2804995674, 82.7302224074, 78.0584160164, 73.2783556203, 68.4234731059, 63.5241715804, 58.5455293912, 53.4926375121, 48.3705869167, 43.1844685788, 37.9393734723, 32.7564707424, 27.5674198879, 22.2950938015, 16.9208406265, 11.4260085058, 5.79194558261, 0., -5.65721456314, -11.5526034027, -17.7478062787, -24.3044629514, -31.2842131808, -38.7486967269, -47.1777160283, -56.9882395771, -67.2420324103, -77.8309468035, -88.6468350325, -99.5815493728, -110.5269421, -130.898857812, -154.564236483, -174.669093835, -189.157228225, -195.972438007, -193.058521539, -178.359277176, -95.5402290233, -64.7497576531, -36.320473725, 39.4150121017, 182.864141589, 515.868654177, 1052.984283, 1808.76676121, 2783.21579022, 4034.55519832, 5797.96237879, 7575.84812321, 8998.68734808, 10279.682266, 11358.5028254, 12072.5459869, 12243.669337, 11858.1609461, 10846.2527933, 9138.17685761}, 
{117.899489604, 118.477822781, 118.808138344, 118.895098789, 118.743366611, 118.357604308, 117.742474375, 116.902639308, 115.842761604, 114.567503758, 113.081528268, 111.389497629, 109.496074338, 107.40592089, 105.108255316, 102.617510339, 99.9492877561, 97.1114091028, 94.1116959124, 90.9579697189, 87.6580520563, 84.2349186636, 80.67914397, 76.99645661, 73.1925852179, 69.2732584282, 65.2442048752, 61.1155102537, 56.8966055582, 52.5828011938, 48.1780566883, 43.6863315697, 39.1115853657, 34.4577776041, 29.8218389375, 25.1489150208, 20.3771140048, 15.4913789421, 10.4766528853, 5.31787888703, 0., -5.38554486058, -10.975016313, -16.7981791129, -22.8847980159, -29.2646377777, -35.967463154, -43.5769005201, -52.5933519154, -61.7222357053, -70.7684310302, -79.5368170306, -87.832272847, -95.4596776199, -113.071515913, -133.610432087, -149.197954571, -157.420174465, -155.86318287, -142.113070886, -113.755929611, -71.0505990986, -22.5536702847, 57.2863878641, 194.021106381, 426.448385193, 747.445463509, 1198.45938787, 1820.93720481, 2676.07658973, 3746.07270261, 5107.50170243, 6493.27281216, 7687.08762929, 8881.51936706, 10014.984599, 10904.3576612, 11351.2889681, 11339.8668589, 10791.9846542, 9629.5356747}, 
{102.855578984, 103.832750167, 104.518854214, 104.924456378, 105.06012191, 104.936416061, 104.563904084, 103.95315123, 103.114722751, 102.059183899, 100.797099926, 99.3390360827, 97.6955576216, 95.8772297943, 93.944343313, 91.8760069217, 89.6475653663, 87.2594129665, 84.7119440421, 82.005552913, 79.1406338989, 76.0270270215, 72.7681870678, 69.3770145265, 65.8664098865, 62.2492736366, 58.5385062657, 54.7800969793, 51.0150630074, 47.1771859379, 43.2659311773, 39.2807641321, 35.2211502089, 31.0865548142, 26.9222572373, 22.6987408553, 18.3830211364, 13.9651926359, 9.43534990934, 4.78358751223, 0., -4.86835538571, -9.89621371286, -15.1013473632, -20.5015287186, -26.1145301608, -31.9581240716, -38.4425623387, -45.9191221076, -53.4671124355, -60.944946413, -68.2110371307, -75.1237976793, -81.5416411492, -99.9335545198, -122.180455817, -139.208698189, -148.297306059, -146.725303844, -131.771715967, -100.715566847, -46.6446937922, 11.0230872987, 97.7389714413, 238.954153652, 429.375684183, 770.290118373, 1250.25597266, 1857.8317635, 2563.49814514, 3410.04722055, 4523.85864428, 5665.97230621, 6652.4368958, 7770.19181506, 8941.26001047, 9889.92364452, 10211.1185657, 10214.9066624, 9967.55579086, 9535.33380749}, 
{97.0841696579, 97.222503046, 97.2054101217, 97.0335763996, 96.7076873941, 96.2284286198, 95.5964855913, 94.812543823, 93.8772888295, 92.7914061254, 91.5555812251, 90.1704996432, 88.6368468943, 86.9553084928, 85.0950471287, 83.0766897705, 80.9232495254, 78.6418596412, 76.2396533656, 73.7237639465, 71.1013246315, 68.419613301, 65.6400743214, 62.7642966916, 59.7938694107, 56.7303814776, 53.5754218913, 50.3363258633, 47.0195651682, 43.6129889541, 40.1158530729, 36.5274133765, 32.846925717, 29.0736459462, 25.2464606045, 21.3395550258, 17.3241146753, 13.1912892754, 8.93222854807, 4.53808221555, 0., -4.68705986094, -9.53628214267, -14.5570431051, -19.7587190082, -25.1506861118, -30.7423206759, -37.4179572678, -45.9304574993, -54.1970723864, -61.8719170007, -68.6091064137, -74.062755697, -77.886979922, -91.1696391353, -106.331819609, -114.929123503, -114.276987846, -101.690849669, -74.4861460011, -29.9783138722, 36.6595606103, 114.011063675, 220.556346098, 374.775558658, 570.116370593, 886.805560268, 1313.28454092, 1837.9947258, 2446.12911056, 3135.87436119, 3889.51572572, 4747.37559906, 5753.13369345, 6769.55354536, 7731.16959208, 8593.74822762, 9353.31992828, 9823.62688448, 9880.55529164, 9399.99134517}, 
{93.0636841069, 92.7091563051, 92.2735422747, 91.7525701679, 91.1419681368, 90.4374643335, 89.6347869103, 88.7296640192, 87.7178238124, 86.5949944421, 85.3569040604, 83.9992808194, 82.5178528714, 80.9083483685, 79.0986192073, 77.1273323069, 75.03829252, 72.8411115256, 70.5454010031, 68.1607726316, 65.6968380906, 63.2533509026, 60.7373316973, 58.1459429476, 55.4763471267, 52.7257067074, 49.891184163, 46.9683728859, 43.9531020416, 40.8462859541, 37.6457241934, 34.3492163298, 30.9545619333, 27.4595605743, 23.8903102512, 20.2267088452, 16.4465121708, 12.5417315887, 8.50437845915, 4.32646414277, 0., -4.5549191721, -9.26443327947, -14.1265987915, -19.1394721775, -24.3011099068, -29.6095684489, -36.2079604554, -45.0673401849, -53.4478039953, -60.8824787511, -66.9044913168, -71.0469685569, -72.8430373358, -82.3414005385, -92.4229932616, -95.1727678332, -87.9729769478, -68.2058732995, -33.2537095827, 19.5012615083, 96.6976371383, 186.141747245, 302.616340538, 460.904165729, 648.025541468, 947.243699599, 1340.02673928, 1807.84275967, 2331.19414387, 2894.44613919, 3416.56633498, 4076.45749588, 5032.08703851, 5947.18225381, 6764.8583376, 7552.88380905, 8491.86308813, 9186.40366826, 9437.4764883, 9046.05248712}, 
{89.4001902854, 89.1496619218, 88.7835548368, 88.3020782703, 87.705441462, 86.9938536517, 86.1675240792, 85.2266619842, 84.1714766065, 83.0021771859, 81.7189729622, 80.322073175, 78.8116870643, 77.1880238697, 75.4216733999, 73.5315822548, 71.538939095, 69.45001157, 67.2710673297, 65.0083740235, 62.6681993012, 60.3056880161, 57.8714803373, 55.3650936379, 52.7860452907, 50.1338526687, 47.4080331447, 44.6184565644, 41.7734371906, 38.8477384514, 35.8366742833, 32.7355586225, 29.5397054054, 26.2444285683, 22.8594504828, 19.3709702564, 15.7640963418, 12.0311955489, 8.16463468763, 4.15678056801, 0., -4.56089583683, -9.25179537336, -14.0461426709, -18.9173817907, -23.838956794, -28.7843117423, -34.86285636, -43.013307938, -50.4928868952, -56.813798731, -61.4882489444, -64.0284430348, -63.9465865013, -70.9147046313, -78.0178638469, -77.5720769814, -67.011442361, -43.7700583119, -5.28202316033, 51.0185647674, 114.47605353, 206.695068996, 334.737310585, 505.664477715, 712.827984173, 997.108146195, 1349.11445206, 1759.45639003, 2224.54206207, 2717.58511545, 3153.94914556, 3690.33809503, 4452.54518111, 5266.51176508, 6080.42098861, 6880.27632335, 7701.51410096, 8344.01888305, 8687.97732142, 8613.57606789}, 
{87.883203151, 87.776768013, 87.5080234844, 87.0823440319, 86.5051041224, 85.7816782224, 84.9174407988, 83.9177663182, 82.7880292474, 81.533604053, 80.1598652019, 78.6721871606, 77.075944396, 75.3765113748, 73.590367283, 71.7158617014, 69.7505036324, 67.6973961641, 65.5596423845, 63.340345382, 61.0426082446, 58.6728919088, 56.2304778726, 53.7180054822, 51.1381140838, 48.4934430239, 45.7866316485, 43.0439637996, 40.2881704503, 37.4653542579, 34.5685541438, 31.5908090292, 28.5251578355, 25.3646394841, 22.1023299588, 18.7312447351, 15.244396483, 11.6348165429, 7.89553625538, 4.01958696091, 0., -4.59541160783, -9.28967263747, -14.0310261848, -18.7677153458, -23.4479832162, -28.0200728919, -33.5316268369, -40.8650896042, -47.3399148364, -52.4579541583, -55.7210591946, -56.6310815703, -54.68987291, -59.1530233295, -63.3521338873, -59.8804906174, -46.244898545, -19.9521626952, 21.4909119068, 80.5775202359, 130.212789464, 223.762993343, 361.381878189, 543.22319032, 771.245599638, 1041.15071917, 1355.25818621, 1715.88763805, 2138.46571177, 2585.99104527, 2980.5656999, 3430.34641526, 4007.51220707, 4725.52399373, 5538.23238851, 6340.0604201, 7004.59713009, 7548.19509298, 7938.53728464, 8143.30668098}, 
{89.001695046, 88.839773766, 88.4993360361, 87.9877118473, 87.3122311907, 86.4802240571, 85.4990204378, 84.3759503236, 83.1183437056, 81.7335305748, 80.2288409223, 78.611604739, 76.8891520161, 75.0688127444, 73.1796299352, 71.2151978438, 69.1674671453, 67.039326612, 64.8336650159, 62.5533711293, 60.2013337243, 57.7738067425, 55.2812301022, 52.7274088912, 50.1161481973, 47.4512531081, 44.7365287114, 42.0059461789, 39.2889438671, 36.5171963118, 33.6822602014, 30.7756922242, 27.7890490685, 24.7138874228, 21.5260091658, 18.2269375358, 14.8193883409, 11.2981407854, 7.65797407386, 3.89366741057, 0., -4.58172329066, -9.22681032155, -13.86404329, -18.4222043934, -22.8300758291, -27.0164397944, -32.0552608273, -38.8484261158, -44.586456087, -48.7331521092, -50.7523155504, -50.1077477789, -46.2632501627, -47.8172245231, -48.4548906443, -41.1700193074, -23.5580048496, 6.78575839161, 52.2658760789, 115.286953875, 162.985467163, 260.493937289, 404.295422089, 590.872979399, 826.606677009, 1083.56701154, 1370.11335777, 1694.60509046, 2083.23117329, 2490.86221431, 2834.34343264, 3213.18556477, 3686.49304282, 4318.65002814, 5070.76556647, 5822.39042939, 6414.30286591, 6923.92948943, 7348.02098256, 7683.32802791}, 
{91.624829329, 91.3631054632, 90.9128956613, 90.2824355449, 89.4799607355, 88.5137068545, 87.3919095235, 86.1228043639, 84.7146269973, 83.1756130452, 81.513998129, 79.7380178703, 77.8559078906, 75.8759038114, 73.8251347088, 71.6998841361, 69.495005495, 67.2148699069, 64.8638484935, 62.4463123761, 59.9666326763, 57.4290420225, 54.838069156, 52.1981043252, 49.5135377783, 46.7887597637, 44.0281605296, 41.2700196099, 38.5475242676, 35.784031599, 32.9701717134, 30.0965747204, 27.1538707292, 24.1326898494, 20.9860001562, 17.7282568992, 14.3767661782, 10.9298615605, 7.38587661341, 3.74314490415, 0., -4.48398332399, -8.98349574252, -13.4119867225, -17.6829057308, -21.7097022344, -25.4058257002, -29.8911460872, -36.1042542394, -41.0739104392, -44.2237186669, -44.9772829031, -42.758207128, -36.9900953222, -35.7967022478, -33.0978847489, -22.154936592, -0.611917816509, 33.8871115382, 83.6980914328, 151.176961828, 205.541667424, 308.969452584, 456.929288815, 644.890147625, 878.012805911, 1127.85984391, 1401.53030656, 1706.12323879, 2066.58021576, 2436.47269176, 2726.72576278, 3049.43610401, 3474.63445714, 4035.89844866, 4700.29090829, 5374.37656133, 5932.44637646, 6438.38445323, 6891.44058726, 7290.86457417}, 
{95.4194996158, 94.9744198833, 94.3458562556, 93.5410785727, 92.5673566743, 91.4319604004, 90.142159591, 88.7052240859, 87.128423725, 85.4190283482, 83.5843077955, 81.6315319067, 79.5679705217, 77.4008934805, 75.1326393639, 72.7735975442, 70.3345306684, 67.8237172249, 65.2494357016, 62.619964587, 59.9435823691, 57.2613441887, 54.5442252012, 51.7959772145, 49.0203520365, 46.221101475, 43.4019773381, 40.5995521447, 37.8454667093, 35.064994771, 32.2485618743, 29.3865935638, 26.4695153841, 23.4877528797, 20.3632017694, 17.1296398156, 13.8260321866, 10.4568216366, 7.02645091937, 3.53936278917, 0., -4.2462036952, -8.44478288285, -12.5002811505, -16.3172420856, -19.8002092758, -22.8537263086, -26.6855878777, -32.3077596035, -36.5085628312, -38.6633793212, -38.1475908336, -34.3365791285, -26.6057259661, -22.9025217571, -17.1796010448, -2.74072752032, 22.7920710373, 61.7967668492, 116.651332136, 189.73373912, 263.59609922, 376.571029829, 526.669187413, 711.901228439, 926.647158506, 1177.87321632, 1459.23331215, 1764.38135629, 2097.49779158, 2420.6569306, 2632.95065195, 2903.44543575, 3359.22954237, 3869.45585076, 4404.85502454, 4965.83074221, 5566.44504354, 6136.45350836, 6634.35268051, 7018.63910385}, 
{118.318050658, 118.129966604, 117.694885592, 117.025064998, 116.132762195, 115.03023456, 113.729739468, 112.243534294, 110.583876412, 108.763023198, 106.793232028, 104.686760276, 102.455865318, 100.112804529, 97.8024253653, 95.4531082331, 92.9831970528, 90.377828994, 87.6221412264, 84.7012709198, 81.6003552438, 77.9938652565, 74.2205093578, 70.308329836, 66.2853689796, 62.1796690766, 58.0192724156, 54.4141936875, 51.886999266, 49.0741689561, 45.7674455497, 41.7585718386, 36.8392906148, 30.8013446702, 19.4803395649, 5.17068494195, -9.53374740007, -24.0320194936, -37.723193371, -50.0063310647, -60.2804946072, -59.358433914, -56.411351669, -52.024138439, -46.7816847909, -41.2688812917, -36.0706185083, -34.9813991339, -41.313443233, -47.9770922022, -54.2540309486, -59.4259443793, -62.7745174013, -63.5814349217, -64.6966546077, -63.1445589957, -55.6794280976, -40.8530861168, -17.2173572571, 16.675934278, 62.2749642849, 66.6397646057, 129.983543891, 241.577288542, 390.691984959, 581.571158002, 766.553628363, 952.877284857, 1147.7800163, 1377.66409637, 1592.27425931, 1657.61932947, 1807.50387552, 2222.87277396, 2666.89982447, 3106.43833251, 3587.72527696, 4203.25424439, 4781.10871714, 5239.51237915, 5496.68891437}, 
{132.206193534, 131.96780308, 131.43231472, 130.616258442, 129.536164233, 128.20856208, 126.64998197, 124.876953891, 122.906007831, 120.753673775, 118.436481713, 115.97096163, 113.373643514, 110.661057353, 108.044769849, 105.417930051, 102.658923589, 99.7443873024, 96.6509580293, 93.3552726091, 89.8339678807, 85.5861996387, 81.1320291703, 76.5140367184, 71.7748025259, 66.9569068354, 62.1029298899, 58.1251047587, 55.8049887676, 53.1057229264, 49.716780485, 45.3276346936, 39.6277588023, 32.3066260612, 17.1116669304, -2.50868738502, -22.6562192243, -42.426059592, -60.9133394927, -77.2131899307, -90.4207419108, -86.4996608835, -79.4900887735, -70.3007019518, -59.8401767892, -49.0171896566, -38.7404169249, -34.0867251897, -39.5066605881, -45.9422806062, -52.6098424442, -58.7256033027, -63.5058203821, -66.1667508828, -67.1873083737, -64.9849890285, -57.8817114975, -44.835467056, -24.8042469795, 3.25395745659, 40.381154977, 20.8615141556, 67.0261347286, 164.87798928, 300.420050394, 481.526959302, 640.251712569, 788.843075189, 939.549812158, 1126.77061224, 1296.30446911, 1308.00033794, 1403.47732397, 1769.77747748, 2150.25815959, 2515.93763129, 2931.21921158, 3521.22104418, 4074.81691601, 4495.50987217, 4686.80295778}, 
{139.783317364, 139.27687506, 138.431892401, 137.268659381, 135.807465992, 134.068602227, 132.072358081, 129.839023546, 127.388888616, 124.742243284, 121.919377543, 118.940581386, 115.826144807, 112.5963578, 109.472318802, 106.347285652, 103.099315873, 99.7076257497, 96.1514315682, 92.4099496133, 88.4623961703, 83.7807401897, 78.9214996241, 73.9339450911, 68.8673472085, 63.7709765941, 58.6941038654, 54.6347148753, 52.4482396599, 49.9154623554, 46.7004447085, 42.4672484658, 36.8799353739, 29.6025671795, 13.8768710732, -6.57029891943, -27.5259463467, -48.0123740051, -67.0518846912, -83.6667812017, -96.879366333, -91.7044628262, -83.1063825047, -72.0419571362, -59.4680184886, -46.3413983298, -33.6189284274, -26.6289477851, -30.0429241342, -34.3649161629, -38.776783229, -42.46038469, -44.5975799036, -44.3702282273, -41.906580159, -35.7898045672, -24.5314326811, -7.11974737192, 17.4569684891, 50.2104320307, 92.1523603816, 81.6273702355, 133.554593855, 234.752376756, 372.039064455, 552.444741479, 712.932259801, 864.573857721, 1018.44177354, 1206.50979614, 1377.14551207, 1404.43029117, 1495.42037109, 1807.90135078, 2128.51966525, 2438.04443489, 2798.30925876, 3332.46489074, 3837.3289868, 4220.96960809, 4391.45481578}, 
{143.748811267, 142.846128294, 141.567367398, 139.936268289, 137.976570677, 135.712014272, 133.166338785, 130.363283925, 127.326589403, 124.079994928, 120.647240212, 117.052064965, 113.318208896, 109.469411716, 105.697718213, 101.920397689, 98.0419794977, 94.0517779357, 89.9391073009, 85.6932818908, 81.3036160032, 76.3198797641, 71.2316357045, 66.0889021836, 60.9416975607, 55.8400401953, 50.8339484466, 46.8454535544, 44.6655563886, 42.259192917, 39.3223162486, 35.5508794923, 30.640835757, 24.2881381517, 10.3956392005, -7.67597141542, -26.1274918276, -44.0780380788, -60.646726212, -74.9526722702, -86.114992296, -81.007417965, -72.6856211707, -61.9598894395, -49.6405102973, -36.5377712704, -23.4619598849, -15.2354486603, -16.0777463366, -17.2057869925, -17.8008227241, -17.0441056274, -14.1168877982, -8.20042133264, -0.386567067304, 11.3520809641, 29.1875533905, 54.3293973438, 87.9871599559, 131.370388359, 185.688629684, 204.303355211, 275.424712784, 389.3815381, 536.502666852, 720.254683608, 898.560060925, 1077.51280702, 1263.20693013, 1478.57410503, 1679.19534021, 1767.24580294, 1887.3876069, 2145.20209177, 2410.67315057, 2677.45025735, 2990.99213118, 3449.10523187, 3888.39581133, 4234.51603543, 4413.11807007}, 
{146.802064361, 145.464508531, 143.712488474, 141.573085643, 139.073381495, 136.240457483, 133.101395065, 129.683275693, 126.013180825, 122.118191914, 118.025390416, 113.761857786, 109.354675479, 104.83092495, 100.333614068, 95.8164888168, 91.2245200588, 86.5610774599, 81.8295306853, 77.0332494002, 72.1756032698, 66.9460112164, 61.7051523969, 56.4997552253, 51.3765481156, 46.3822594818, 41.5636177379, 37.6597503131, 35.3057433995, 32.8927202848, 30.1862731336, 26.9519941102, 22.955475379, 17.9623091044, 7.28765851952, -6.48752662711, -20.4454187275, -33.9101271594, -46.2057613007, -56.6564305289, -64.586244222, -60.4431045226, -53.6531930801, -44.7664843088, -34.3329526233, -22.9025724377, -11.0253181664, -2.53360955529, -0.766639660697, 1.57431878513, 5.27206964932, 11.109416799, 19.8691641013, 32.3341154234, 45.8406337977, 63.3517541413, 88.0913917562, 121.547291499, 165.207198228, 220.558856799, 289.090012071, 344.255491447, 438.492283029, 566.946560439, 724.764497296, 910.886964763, 1111.09990681, 1325.10309767, 1552.5963116, 1804.65599587, 2046.8519056, 2216.78348705, 2383.43362139, 2589.63739831, 2805.7074247, 3038.8466127, 3311.26454152, 3683.26151538, 4047.76827142, 4354.77360272, 4554.26630235}, 
{151.642465766, 149.920961521, 147.74100439, 145.133111917, 142.127801649, 138.75559113, 135.046997905, 131.032539519, 126.742733517, 122.208097444, 117.459148846, 112.526405267, 107.440384253, 102.231603348, 96.9926523574, 91.7147816883, 86.3845431508, 81.019757922, 75.6382471793, 70.2578320998, 64.8963338609, 59.4015274009, 53.9847646866, 48.6873514458, 43.5505934064, 38.6157962964, 33.9242658436, 29.9800346684, 27.2176051383, 24.5718501326, 21.8961933918, 19.0440586565, 15.8688696672, 12.2240501643, 5.17261623763, -3.66678630863, -12.464290107, -20.7957165931, -28.2368872025, -34.3636233707, -38.7517465332, -36.0461007217, -31.4344865409, -25.1737271917, -17.5206458747, -8.73206579091, 0.935189858974, 8.84918778995, 12.7348834281, 18.0146130501, 25.3959244703, 35.586365503, 49.2934839625, 67.2248276633, 85.2429253322, 107.12030154, 136.996227455, 176.569259503, 227.537954113, 291.600867711, 370.456556726, 456.849801308, 568.613096105, 705.628530905, 867.778195499, 1050.27176402, 1264.51658831, 1504.78790423, 1765.36094763, 2046.44792559, 2324.51316028, 2573.37995731, 2787.61300456, 2949.16496828, 3122.61129677, 3326.92501497, 3561.12320245, 3847.05318907, 4135.1972489, 4400.36675847, 4617.37309433}, 
{160.969404601, 159.004433012, 156.526663908, 153.570347587, 150.169734346, 146.35907448, 142.172618289, 137.644616068, 132.809318114, 127.700974725, 122.353836198, 116.80215283, 111.080174917, 105.222152757, 99.287479067, 93.2944989575, 87.2596543682, 81.2120529217, 75.1808022407, 69.1950099482, 63.2837836667, 57.4288211722, 51.7131875589, 46.1725380744, 40.8425279662, 35.7588124817, 30.9570468683, 26.7087361373, 23.2499460507, 20.052388134, 17.0559550516, 14.2005394684, 11.426034049, 8.672331458, 4.67019956207, 0.124427785877, -4.16866902665, -8.02188172609, -11.248001163, -13.6598181881, -15.0701236518, -13.8509847851, -11.4548898616, -7.89360353514, -3.17889045974, 2.67748471071, 9.66375732229, 16.2855616354, 21.2713104644, 28.1543076824, 37.5247723176, 49.9729233983, 66.0889799528, 86.4631610092, 106.288210432, 129.568809737, 160.718205526, 201.430625764, 253.400298417, 318.321451452, 397.888312835, 497.452307162, 611.642943523, 743.608536624, 896.497401171, 1064.33926045, 1272.7748963, 1514.01040128, 1780.25186791, 2065.64235117, 2356.57705631, 2657.37182749, 2903.98034641, 3031.74249957, 3170.37357591, 3346.3769782, 3542.56482667, 3752.59970075, 3970.43362559, 4189.91995122, 4404.91202768}, 
{177.482269984, 175.503868754, 172.943215793, 169.838793128, 166.22908279, 162.152566804, 157.647727201, 152.753046007, 147.507005251, 141.948086961, 136.114773165, 130.045545892, 123.778887169, 117.353279024, 110.830740185, 104.234863278, 97.5874593057, 90.9221960584, 84.2727413275, 77.6727629037, 71.155928578, 64.7702853846, 58.5331359992, 52.4761623408, 46.6310463282, 41.0294698804, 35.7031149163, 30.748284237, 26.2515705827, 22.0901399627, 18.2694361413, 14.7949028829, 11.6719839518, 8.90612311218, 6.4000957001, 4.22429390229, 2.45688145316, 1.12430209553, 0.252999572246, -0.130582373875, 0., 0.107665064499, 0.860208649684, 2.36190121341, 4.71701321354, 8.02981510794, 12.4045773545, 17.1481302409, 21.6871289828, 28.0326145621, 36.6126437703, 47.8552733989, 62.1885602395, 80.0405610837, 97.4443919945, 117.608365308, 144.073471007, 178.166714689, 221.215101948, 274.545638381, 339.485329585, 421.429031372, 513.437616799, 619.067664728, 741.875754023, 879.019633126, 1049.83962165, 1250.21376339, 1476.02010212, 1723.93172955, 1987.44154572, 2289.09571141, 2536.59023694, 2645.32769005, 2757.98307125, 2901.89401641, 3057.58612684, 3212.02049823, 3373.22828331, 3542.05762948, 3719.35668413}, 
{203.880451035, 202.208214496, 199.864408804, 196.892449016, 193.335750186, 189.237727371, 184.641795624, 179.591370003, 174.129865562, 168.300697356, 162.147280441, 155.713029873, 149.041360706, 142.175687997, 135.235081699, 128.215097304, 121.105563558, 113.934420932, 106.729609897, 99.5190709247, 92.3307444852, 85.1683128927, 78.0873249929, 71.1190714743, 64.2948430255, 57.6459303353, 51.2036240921, 45.0011084845, 39.0712831798, 33.4409112924, 28.1405146892, 23.2006152371, 18.6517348029, 14.5243952535, 10.981991859, 7.97099028647, 5.42779827188, 3.35575952556, 1.7582177579, 0.63851667926, 0., -0.204729395662, 0.0854206845436, 0.880801606736, 2.19176473703, 4.02866144156, 6.40184308642, 8.80951186649, 10.8268265181, 13.6887455693, 17.6135694072, 22.8195984185, 29.5251329901, 37.948473509, 47.1793729149, 58.1500548284, 71.8781689391, 88.8128506871, 109.403235513, 134.098458857, 163.34765616, 184.145996305, 219.852907445, 270.187002345, 334.866893766, 420.243061128, 509.67555522, 610.841165126, 731.41667993, 883.008517692, 1061.50458057, 1288.88822286, 1489.49726615, 1597.87823761, 1694.42859193, 1798.16764364, 1908.18381565, 2037.43502931, 2163.33210389, 2275.40424177, 2363.18064536}
};

  const double a4tab[21][81] = {
{-631.969627765, -629.931558839, -626.939860408, -622.99370354, 
-618.0922593, -612.234698757, -605.420192976, -597.647913023, -588.917029966, -579.226714872, -568.576138807, -556.964472837, -544.39088803, -530.854555451, -515.846824755, -499.688116687, -482.737291838, -465.097391946, -446.871458752, -428.162533995, -409.073659417, -390.431317392, -371.515196916, -352.32842762, -332.874139137, -313.155461098, -293.175523133, -272.835853015, -252.053245413, -231.073771605, -209.941814844, -188.70175838, -167.397985465, -146.074879349, -125.285445789, -104.752311583, -84.1596030453, -63.4476688976, -42.5568578602, -21.4275186541, 0, 20.4923780371, 41.5808062822, 63.5035042161, 86.4986913197, 110.804587074, 136.659410959, 164.137905434, 193.339377336, 224.892938377, 259.103184964, 296.2747135, 336.712120392, 380.720002045, 451.856144032, 535.715090425, 616.131221344, 688.652885521, 748.828431691, 792.206208588, 814.334564947, 847.454547018, 888.473550546, 839.169609957, 601.320759679, -334.087967929, -1051.9917624, -2143.5602482, -4199.9630498, -7855.60243017, -13571.9500983, -23786.9485008, -33106.7676873, -37360.653705, -39761.0933291, -40691.7902113, -39750.1262669, -36377.5392781, -31425.5652961, -25265.0229552, -18266.7308898}, 
{-505.215415872, -508.196275589, -509.772602555, -509.990386863, -508.895618606, -506.534287876, -502.952384768, -498.195899373, -492.310821785, -485.343142096, -477.3388504, -468.343936788, -458.404391355, -447.566204193, -435.751584423, -423.084826434, -409.699594285, -395.667196444, -381.058941378, -365.946137556, -350.400093444, -334.883419369, -319.022080471, -302.833343751, -286.334476209, -269.542744842, -252.475416653, -235.102085808, -217.39950989, -199.498947083, -181.437021088, -163.250355605, -144.975574336, -126.649300981, -108.741297522, -91.0141828036, -73.19778901, -55.2401449874, -37.0892795811, -18.6932216368, 0, 17.9969170254, 36.5393224046, 55.8235696442, 76.046012251, 97.4030037317, 120.090897593, 144.841730923, 172.307047943, 201.402321533, 232.106400453, 264.398133462, 298.256369319, 333.659956784, 398.011849858, 473.942308238, 542.005707696, 596.571522905, 632.009228542, 642.68829928, 622.978209795, 589.70006138, 568.112915425, 450.122799986, 127.635743119, -717.548638845, -1675.26318027, -3105.72753016, -5369.16133753, -8880.19602788, -13835.8159207, -21538.4449876, -28962.0227248, -33785.1131727, -37519.4657086, -40137.7204646, -41147.9006977, -39965.1597083, -36840.7474987, -31740.8798606, -24631.7725858}, 
{-416.815774728, -419.05759752, -420.469505287, -421.053328761, -420.810898675, -419.744045761, -417.854600751, -415.144394377, -411.615257372, -407.269020467, -402.107514396, -396.13256989, -389.346017681, -381.749688503, -373.154364742, -363.682735012, -353.471949442, -342.562916123, -330.996543149, -318.81373861, -306.055410599, -292.66611878, -278.796426051, -264.500546883, -249.832695745, -234.847087107, -219.597935441, -204.212584583, -188.807389811, -173.261704922, -157.600051513, -141.846951179, -126.026925517, -110.164496124, -94.647615966, -79.2708986292, -63.801447545, -48.1894474268, -32.3850829878, -16.3385389411, 0, 15.8033661418, 32.11947511, 49.1192595499, 66.9736521067, 85.8535854257, 105.929992152, 128.799096163, 155.842954272, 183.824465004, 212.335846965, 240.96931876, 269.317098996, 296.971406276, 353.483090379, 419.492430494, 473.371969443, 508.586135416, 518.599356602, 496.87606119, 436.880677369, 151.559475822, 37.868400673, -64.6700000254, -316.533178219, -780.861525805, -1858.22995323, -3592.31237482, -6026.78270487, -9156.64679336, -13171.5827477, -18819.4884661, -24531.9545763, -29153.6948003, -33275.8654062, -36682.7175389, -38862.2880824, -39230.2154833, -37787.9461988, -34337.2185386, -28679.7708119}, 
{-353.99543721, -356.846771122, -358.781066799, -359.821143027, -359.989818592, -359.30991228, -357.804242876, -355.495629166, -352.406889936, -348.560843972, -343.98031006, -338.688106985, -332.707053534, -326.059968491, -318.716481992, -310.733060156, -302.170195253, -293.061585367, -283.440928582, -273.341922985, -262.798266661, -251.899022138, -240.614876848, -228.971882667, -216.996091472, -204.713555139, -192.150325543, -179.343936697, -166.330197286, -153.107704085, -139.697846853, -126.122015352, -112.401599339, -98.5579885754, -84.9010555343, -71.2696947575, -57.4809639164, -43.4972460157, -29.2809240596, -14.7943810529, 0, 14.9211025452, 30.2554857294, 46.0709751501, 62.4353964044, 79.4165750898, 97.0823368036, 117.335048534, 141.801415471, 166.162665223, 189.74174063, 211.861584536, 231.845139783, 249.015349213, 296.887652213, 353.154693521, 392.920896561, 408.515411585, 392.267388849, 336.505978608, 233.560331114, 83.6241198527, -85.1602406681, -353.844804647, -803.481626284, -1559.26357874, -2593.36231194, -4039.79843909, -6032.59257339, -8771.42456159, -12193.3373163, -16518.2937717, -20947.1316025, -24836.534726, -28711.176873, -32345.1710442, -35145.8060321, -36448.1361303, -36242.9778352, -34297.5913714, -30379.2369634}, 
{-301.81291907, -306.111578654, -309.333143977, -311.522135868, -312.723075159, -312.98048268, -312.338879262, -310.842785736, -308.536722934, -305.465211685, -301.67277282, -297.20392717, -292.103195567, -286.41509884, -280.368561535, -273.89145002, -266.897672114, -259.394030394, -251.387327432, -242.884365806, -233.891948088, -224.09953253, -213.875093457, -203.26926087, -192.33266477, -181.115935157, -169.669702033, -158.158283636, -146.728915186, -135.160385123, -123.457162371, -111.623715852, -99.6645144897, -87.584027207, -75.5054138548, -63.3580590421, -51.0623629979, -38.5985173968, -25.9467139136, -13.0871442231, 0, 13.2534853556, 26.7853543118, 40.6266076113, 54.8082459971, 69.3612702119, 84.3166809985, 100.976510232, 120.451802149, 139.73437689, 158.339155715, 175.781059884, 191.575010658, 205.235929297, 256.091490621, 318.470922689, 363.689803549, 383.11968475, 368.132117838, 310.098654363, 200.390845871, 15.8388731288, -179.533402536, -466.97946729, -927.752807298, -1540.59831927, -2639.62335398, -4183.07207432, -6129.18864318, -8375.48674927, -11062.401978, -14579.6216436, -18215.7203306, -21430.6833921, -25067.8557787, -28848.6760519, -31873.7742811, -32804.6184259, -32680.1333768, -31728.7438069, -30178.8743889}, 
{-284.781252339, -285.771104018, -286.206934048, -286.095461564, -285.443405697, -284.257485578, -282.544420342, -280.310929119, -277.563731042, -274.309545243, -270.555090856, -266.307087011, -261.572252841, -256.357307478, -250.562720525, -244.262424953, -237.538396358, -230.419086336, -222.932946481, -215.10842839, -206.973983659, -198.687787305, -190.130651808, -181.313113071, -172.245706996, -162.938969485, -153.40343644, -143.664460555, -133.745168124, -123.620666407, -113.295475215, -102.774114355, -92.0611036353, -81.1609628657, -70.1804633526, -59.0594402339, -47.7299886594, -36.175713739, -24.3802205824, -12.3271142995, 0, 12.6871401941, 25.6714627735, 38.9597472169, 52.5587730027, 66.4753196096, 80.716166516, 98.1049404153, 121.042003, 142.79873034, 162.23817143, 178.223375263, 189.617390832, 195.283267133, 230.089209048, 270.121272553, 288.740081975, 277.444129231, 227.731906238, 131.101904912, -20.9473828286, -244.85163344, -502.495548545, -852.470424987, -1353.36755961, -1979.98633117, -2997.60423152, -4364.26306021, -6038.00461676, -7966.88568998, -10138.9031116, -12451.972398, -15147.2177527, -18461.0347505, -21797.171786, -24909.918568, -27667.4602203, -30054.9921093, -31475.7732159, -31531.4507908, -29823.6720845}, 
{-273.657349424, -272.684103392, -271.445420469, -269.928714809, -268.121400565, -266.01089189, -263.584602939, -260.829947864, -257.734340819, -254.285195958, -250.469927433, -246.275949399, -241.690676009, -236.701521416, -231.057822328, -224.897601558, -218.376903337, -211.531838591, -204.398518248, -197.013053236, -189.411554483, -181.930394098, -174.2639537, -166.40687609, -158.353804069, -150.099380438, -141.638247999, -132.952397953, -124.025722555, -114.883116027, -105.524358124, -95.9492286033, -86.1575072218, -76.1489737363, -65.9907178979, -55.6399389155, -45.0487409277, -34.2031359885, -23.0891361517, -11.6927534709, 0, 12.2961137417, 24.8761086464, 37.7135071409, 50.781831652, 64.0546046063, 77.5053484307, 94.7805372135, 118.974740105, 141.279032613, 160.175600765, 174.146630588, 181.674308108, 181.240819351, 204.397425441, 228.706694809, 229.228012624, 197.679558897, 125.779513639, 5.24605686281, -172.20263142, -429.141372332, -724.391694961, -1104.69335281, -1616.78609938, -2214.08339703, -3173.52925348, -4429.87276878, -5917.86304298, -7570.66086314, -9327.78026823, -10869.4590696, -12910.3172602, -16076.8026498, -19097.7219117, -21751.8649469, -24278.4468593, -27272.369335, -29445.3625284, -30151.1496903, -28743.4540715}, 
{-262.758924574, -262.113472023, -261.081805142, -259.667040125, -257.872293172, -255.70068048, -253.155318244, -250.239322664, -246.955809936, -243.307896257, -239.298697824, -234.931330836, -230.208911489, -225.13455598, -219.6072599, -213.696006522, -207.47766059, -200.976635311, -194.217343893, -187.224199541, -180.021615464, -172.795254687, -165.386010893, -157.796027583, -150.027448261, -142.082416427, -133.963075583, -125.698086784, -117.312126502, -108.743931745, -99.9848790391, -91.0263449067, -81.8597058721, -72.4763384594, -62.8934817729, -53.086151574, -43.0274059405, -32.7033314195, -22.1000145577, -11.2035419023, 0, 12.371283182, 24.9605193553, 37.6646788112, 50.3807318412, 63.0056487368, 75.4363997892, 91.193625081, 113.253464947, 132.848281342, 148.40371818, 158.345419375, 161.099028843, 155.090190498, 170.755036726, 186.269259609, 177.385498307, 135.981915863, 53.9366753214, -76.8720602719, -264.566127873, -474.727112665, -777.463354667, -1194.26756567, -1746.63245746, -2408.93019051, -3318.88386829, -4441.44199337, -5741.5530683, -7205.02519288, -8734.22807794, -9998.84669636, -11635.0630563, -14157.6474741, -16849.4723403, -19507.4035558, -22093.1103792, -24715.230032, -26729.7245257, -27749.1508132, -27386.0658472}, 
{-259.275990836, -259.123882778, -258.418165934, -257.180325843, -255.431848043, -253.194218075, -250.488921476, -247.337443786, -243.761270544, -239.781887289, -235.42077956, -230.699432895, -225.639332835, -220.261964917, -214.628001493, -208.734138382, -202.574105133, -196.161371948, -189.509409025, -182.631686567, -175.541674773, -168.261426091, -160.794643208, -153.153611056, -145.350614569, -137.39793868, -129.307868322, -121.163341855, -113.036681117, -104.771230627, -96.3505877579, -87.7583498815, -78.9781143694, -69.9934785935, -60.7719554264, -51.3073173527, -41.5945543826, -31.6205538392, -21.3722030456, -10.8363893248, 0, 12.5572793642, 25.2470154816, 37.8879768246, 50.2989318654, 62.2986490763, 73.7058969297, 87.8304028973, 107.457336088, 124.058183296, 136.034271464, 141.786927532, 139.717478442, 128.227251134, 136.537729166, 143.553295923, 125.845384743, 75.5113142419, -15.3515969646, -154.646030261, -350.274667033, -513.504732378, -819.32883499, -1266.603406, -1854.18487656, -2584.95999694, -3447.84377019, -4446.52897172, -5584.70837696, -6911.11538595, -8294.32090024, -9428.86076201, -10782.620903, -12685.5643737, -15060.758387, -17726.13186, -20334.1358205, -22453.175157, -24153.3172041, -25333.6072619, -25893.0906303}, 
{-265.219436067, -264.812067215, -263.797412645, -262.203671835, -260.059044262, -257.391729402, -254.229926734, -250.601835734, -246.535655879, -242.059586646, -237.201827514, -231.990577958, -226.454037457, -220.620405487, -214.594700787, -208.35652669, -201.879671646, -195.176622389, -188.259865651, -181.141888162, -173.835176657, -166.329316483, -158.662858593, -150.851452557, -142.910747946, -134.856394327, -126.704041273, -118.561145264, -110.521367698, -102.380837329, -94.1179271695, -85.7110102336, -77.1384595339, -68.3786480833, -59.3472390162, -50.0622758436, -40.5465489384, -30.791258101, -20.7876031318, -10.5267838313, 0, 12.6060558568, 25.246461822, 37.6808032737, 49.66866559, 60.969634149, 71.3432943288, 84.1830230763, 102.336176673, 116.873529832, 126.079231789, 128.237431777, 121.63227903, 104.547922783, 104.223051879, 100.625068123, 71.5296422594, 9.29851173747, -93.7065859944, -245.123913488, -452.591733293, -608.034342415, -925.808021943, -1392.75892049, -1995.73318667, -2752.26407127, -3573.50306728, -4483.12055131, -5504.7869, -6723.05204564, -7978.94769789, -8949.69780629, -10072.1187352, -11630.2425955, -13722.4852483, -16192.1788187, -18641.4825582, -20530.5584654, -22128.7304242, -23429.4704719, -24426.250646}, 
{-276.149614772, -275.332737385, -273.87684957, -271.813326606, -269.173543774, -265.988876353, -262.29069962, -258.110388857, -253.479319341, -248.428866353, -242.990405172, -237.195311077, -231.074959346, -224.66072526, -218.05385791, -211.241530146, -204.205625679, -196.963227679, -189.531419314, -181.927283752, -174.167904162, -166.271413666, -158.253700473, -150.131702746, -141.922358648, -133.642606341, -125.309383989, -117.043658597, -108.95076559, -100.798898289, -92.5627558036, -84.2170372438, -75.7364417194, -67.0956683406, -58.1385388872, -48.9225459458, -39.5150889403, -29.9176368567, -20.1316586813, -10.1586234004, 0, 12.4032352855, 24.7045525509, 36.6179146432, 47.8572844091, 58.1366246954, 67.169898349, 78.5055660332, 95.1159092306, 107.546168562, 113.953377545, 112.494569696, 101.326778532, 78.607037569, 70.1879876032, 56.7063199097, 16.7022903322, -57.3319751938, -172.904350733, -337.522710351, -558.694928111, -734.519253892, -1070.84192717, -1551.18733437, -2159.07986188, -2908.05280354, -3707.60876617, -4577.282537, -5536.60890327, -6666.21730929, -7802.35857108, -8602.48396025, -9543.97303295, -10946.0445439, -12803.8491331, -14981.7793292, -17174.8197468, -18960.2886513, -20555.5776181, -21959.8704604, -23172.3509912}, 
{-290.803666917, -289.274471634, -287.12888331, -284.394980828, -281.100843072, -277.274548925, -272.944177271, -268.137806994, -262.883516977, -257.209386103, -251.143493257, -244.713917322, -237.948737182, -230.87603172, -223.516589752, -215.903101885, -208.068810552, -200.04328576, -191.856097516, -183.536815826, -175.115010698, -166.733726383, -158.293387073, -149.807891202, -141.291137207, -132.757023524, -124.219448589, -115.792996576, -107.577122422, -99.3449734635, -91.0695664915, -82.7239182966, -74.2810456697, -65.7139654016, -56.764490888, -47.5518990314, -38.2129678122, -28.768004765, -19.2373174241, -9.64121332433, 0, 11.7686103835, 23.2513320955, 34.1374747746, 44.1163480595, 52.8772615886, 60.1095250006, 69.6618430243, 84.7579215045, 95.1414813656, 98.812984183, 93.7728915319, 78.0216649879, 49.559766126, 33.8139513502, 11.434708033, -39.0035261187, -125.110111475, -254.494408406, -434.765777282, -673.533578474, -911.692133432, -1278.54529359, -1765.94828002, -2365.7563138, -3055.88020564, -3862.57129962, -4760.95165824, -5726.14334399, -6770.29687025, -7757.44894691, -8304.76650745, -9080.37075902, -10594.8808504, -12284.7375827, -14027.3107432, -15837.8677096, -17770.7606029, -19586.1049804, -21144.3247566, -22305.8438458}, 
{-366.917658252, -366.04082161, -364.3353285, -361.84652871, -358.619772027, -354.700408243, -350.133787144, -344.96525852, -339.240172159, -333.003877851, -326.301725384, -319.179064547, -311.681245129, -303.853616919, -296.165810556, -288.394774215, -280.285339879, -271.796074193, -262.885543806, -253.512315365, -243.634955515, -232.259700719, -220.428971134, -208.232856729, -195.761447476, -183.104833344, -170.353104306, -157.864767808, -145.957998358, -134.081069086, -122.215083842, -110.341146476, -98.4403608387, -86.4938307804, -74.5241682595, -62.4862189455, -50.3316865268, -38.0331847132, -25.5633272144, -12.8947277401, 0, 13.9806334723, 28.1272088661, 42.352153547, 56.5678948806, 70.6868602322, 84.6214769676, 102.120048684, 126.354492674, 148.065218129, 165.607163406, 177.335266857, 181.604466839, 176.769701706, 175.705105693, 167.580720852, 140.467550085, 89.7507502222, 10.8154780919, -100.953109476, -250.169855651, -262.408271935, -468.330676854, -831.084029719, -1313.81528984, -1926.70854343, -2516.88613781, -3103.93913301, -3707.45858908, -4412.81050955, -5042.26112399, -5101.50126061, -5475.2765331, -6853.49878319, -8326.85021145, -9760.51686418, -11319.4203541, -13314.5514006, -15172.8887303, -16622.7801097, -17392.5733055}, 
{-414.959592893, -413.794235729, -411.63508889, -408.542018192, -404.574889456, -399.793568497, -394.257921136, -388.027813189, -381.163110476, -373.723678814, -365.769384022, -357.360091917, -348.555668318, -339.415979043, -330.624653239, -321.846961766, -312.700959348, -303.11892617, -293.033142416, -282.37588827, -271.079443916, -257.608089298, -243.564845731, -229.084734294, -214.302776062, -199.353992111, -184.373403517, -169.886016083, -156.358235741, -142.992585544, -129.765740236, -116.654374562, -103.635163268, -90.6847810983, -77.8287716484, -65.0128950505, -52.1792122807, -39.2944023649, -26.325144329, -13.2381171987, 0, 13.8679306068, 27.8920801119, 42.0442563708, 56.2962672394, 70.6199205733, 84.9870242281, 103.072365783, 127.994315803, 150.8704347, 170.169001572, 184.358295518, 191.906595634, 191.282181019, 189.728565601, 180.162782898, 154.837619819, 110.426453698, 43.6026618723, -48.9603783248, -170.589289558, -105.383549457, -255.919705834, -574.717356221, -1014.29609815, -1596.90983573, -2103.06732323, -2568.96865665, -3030.81393194, -3600.50829153, -4087.13669206, -3978.29577273, -4177.31107248, -5387.57191481, -6647.53588253, -7838.6230605, -9183.61260824, -11099.3954574, -12885.0563765, -14220.7520797, -14786.6392811}, 
{-443.141012235, -441.040932434, -437.811578902, -433.525380004, -428.254764105, -422.072159568, -415.049994758, -407.260698039, -398.776697777, -389.670422335, -380.014300078, -369.88075937, -359.342228576, -348.47113606, -337.987644767, -327.554423801, -316.785111615, -305.619647818, -293.997972017, -281.860023824, -269.145742846, -254.238808814, -238.850351368, -223.135240271, -207.248345284, -191.344536168, -175.578682685, -160.571709249, -146.874509986, -133.527564754, -120.496510151, -107.746982774, -95.2446192185, -82.9550560824, -70.7265382265, -58.5989646882, -46.6211205418, -34.7826538936, -23.0732128497, -11.4824455165, 0, 12.0604105385, 23.9398518006, 35.5553244327, 46.8238290816, 57.6623663936, 67.9879370154, 81.4009117357, 100.948190904, 117.739416471, 130.196022445, 136.739442836, 135.791111652, 125.772462901, 113.329208489, 91.6800742116, 53.421244269, -4.70805540473, -85.9685988759, -193.621160211, -330.926513476, -295.282721036, -464.969660124, -795.376229309, -1241.89132716, -1824.34163156, -2339.16993116, -2819.0898411, -3296.81497655, -3876.56884781, -4376.53538463, -4336.96366225, -4532.42636247, -5557.15719197, -6612.0770037, -7614.59022016, -8776.60434366, -10507.1028866, -12132.2622412, -13347.3992015, -13847.8305613}, 
{-459.673457673, -456.287130169, -451.648212962, -445.840544874, -438.947964724, -431.054311332, -422.24342352, -412.599140106, -402.205299912, -391.145741757, -379.504304462, -367.364826846, -354.811147732, -341.927105937, -329.349312106, -316.811919587, -304.007239335, -290.906045263, -277.479111284, -263.697211308, -249.531119249, -233.611775966, -217.434828548, -201.156091033, -184.931377458, -168.91650186, -153.267278276, -138.636815154, -125.603496505, -113.134048633, -101.182368067, -89.702351331, -78.6478949531, -67.9728954596, -57.2657051657, -46.7114652001, -36.5229867777, -26.7289354724, -17.3579768584, -8.43877650958, 0, 9.18798201886, 17.6644299001, 25.2268989195, 31.6729443526, 36.8001214752, 40.4059855628, 45.9878039567, 56.4869179886, 62.8544488343, 63.3857497856, 56.3761741343, 40.1210751723, 12.9158061915, -16.0480882252, -55.3682251631, -114.467015375, -197.277831446, -307.734045964, -449.769031515, -627.316160685, -686.770541705, -919.137819319, -1291.58861877, -1771.29356529, -2366.86147998, -2943.24903683, -3517.35229832, -4106.06732689, -4784.17081002, -5394.91693521, -5579.31854758, -5886.5743884, -6722.31156131, -7585.75598269, -8441.37923115, -9444.55543213, -10919.4838013, -12324.1606465, -13411.8800102, -13935.9359347}, 
{-472.7684706, -468.039047376, -461.928405493, -454.531443529, -445.943060062, -436.25815367, -425.571622931, -413.978366425, -401.573282729, -388.451270422, -374.707228081, -360.436054285, -345.732647613, -330.691906642, -315.804182225, -300.914208387, -285.836785163, -270.585924635, -255.175638887, -239.619940004, -223.93284007, -207.186907449, -190.467617774, -173.919002954, -157.685094903, -141.909925533, -126.737526756, -112.796301644, -100.641870711, -89.2600790969, -78.598288459, -68.6038604526, -59.2241567334, -50.4065389571, -41.4945096378, -32.817433928, -24.7503864561, -17.3642432748, -10.7298804367, -4.91817399422, 0, 5.88055379932, 10.4597203784, 13.4005210178, 14.3659769981, 13.0191095999, 9.02294010361, 5.7151598599, 5.88129706964, 0.397817184416, -12.5642926442, -34.8340452646, -68.2404535253, -114.612530275, -160.958447126, -218.482935184, -299.512597924, -408.907928785, -551.529421204, -732.23756862, -955.892864469, -1134.5117665, -1442.08146301, -1861.88249438, -2377.19540096, -2982.32693, -3633.35971552, -4326.80564026, -5059.17658697, -5866.49280978, -6626.74107731, -7107.17404714, -7585.70713558, -8243.09196948, -8933.85522722, -9671.95098146, -10533.6257454, -11718.3483146, -12870.4059146, -13823.3530411, -14410.7441898}, 
{-490.63759241, -484.8029025, -477.435570919, -468.642006697, -458.528618867, -447.201816459, -434.768008504, -421.333604034, -407.005012079, -391.888641673, -376.090901844, -359.718201625, -342.876950048, -325.673556142, -308.446782089, -291.156049468, -273.743191753, -256.267092059, -238.786633501, -221.360699193, -204.048172252, -186.424119961, -169.098059548, -152.195692408, -135.842719937, -120.164843531, -105.287764587, -91.7651365679, -80.0863080166, -69.3536980609, -59.5192458062, -50.5348903579, -42.3525708214, -34.924226302, -27.4611888148, -20.3839082134, -14.1688950448, -8.91957347409, -4.73936766634, -1.73170178662, 0, 2.76803463136, 3.71962920327, 2.41773191446, -1.57470903633, -8.69474545037, -19.3794291289, -30.5349031403, -39.5978718414, -55.4481930862, -79.9565810825, -114.993750038, -162.430414161, -224.13728766, -283.956990798, -355.164875808, -452.400942191, -581.22340178, -747.190466407, -955.860347906, -1212.79125811, -1493.17115046, -1857.4578708, -2304.78582592, -2834.28942261, -3428.59553066, -4127.55704246, -4910.49947888, -5756.74836077, -6666.71347872, -7556.46754446, -8322.34377935, -8975.77658934, -9479.55536314, -10021.657145, -10659.2663591, -11389.9751552, -12285.5065397, -13180.6523674, -13990.976829, -14632.0441153}, 
{-521.492364498, -515.084913984, -506.953123664, -497.216165108, -485.993209889, -473.403429578, -459.565995747, -444.600079967, -428.624853812, -411.759488852, -394.12315666, -375.835028807, -357.014276864, -337.780072405, -318.371638664, -298.832202094, -279.195901762, -259.557353664, -240.011173797, -220.651978157, -201.57438274, -182.783330198, -164.475494374, -146.757875767, -129.737474876, -113.521292199, -98.2163282368, -84.2582877718, -72.0334838359, -60.8629474416, -50.720214586, -41.5788212661, -33.4123034788, -26.1941972213, -19.2139798684, -12.8779253981, -7.64408801148, -3.62592224355, -0.936882629099, 0.309576297031, 0, 0.480333266428, -1.16193765743, -5.37992720387, -12.6267498052, -23.3555198936, -38.0193519016, -53.8802676298, -68.6797887324, -90.5012965843, -121.093591768, -162.205474864, -215.585746457, -282.983207127, -347.598841825, -422.914866992, -523.817486986, -655.849304789, -824.552923384, -1035.47094576, -1294.14597489, -1617.41344861, -1988.92432228, -2418.82658319, -2917.26821866, -3463.52483099, -4143.89609292, -4931.48342612, -5799.38825228, -6728.01144847, -7668.55607016, -8626.64136262, -9402.734735, -9791.75868893, -10214.4441438, -10756.2862519, -11359.7635332, -12002.7685898, -12664.5543272, -13323.9099093, -13959.6244998}, 
{-573.544328257, -567.391300272, -559.264478152, -549.297849489, -537.625401876, -524.381122905, -509.69900017, -493.713021263, -476.557173777, -458.365445304, -439.271823436, -419.410295768, -398.914849891, -377.919473397, -356.673278919, -335.23742553, -313.664357843, -292.064515578, -270.54833845, -249.226266178, -228.208738478, -207.724454857, -187.749262756, -168.377269406, -149.702582038, -131.819307882, -114.821554169, -98.9907231035, -84.5800735815, -71.2358691551, -58.9761692761, -47.8190333964, -37.7825209678, -28.8846914421, -20.8011199706, -13.7665228236, -8.04154082388, -3.71428575654, -0.872869406709, 0.394596440511, 0, -0.352641544002, -2.79107423583, -7.65091515046, -15.2677813629, -25.977289948, -40.1150579809, -55.4388161944, -70.0936535914, -90.5792079169, -118.277800938, -154.57175442, -200.843390132, -258.475029839, -314.439122789, -379.233728695, -464.447671122, -574.410692172, -713.452533947, -885.902938549, -1096.09164808, -1361.903416, -1660.13809704, -2002.53273596, -2400.82437753, -2844.97238002, -3400.43194215, -4052.80709395, -4787.70186548, -5593.56535065, -6447.46638793, -7421.88041538, -8212.53355789, -8539.75889351, -8877.49863129, -9315.97154808, -9789.1507513, -10251.944578, -10731.7661161, -11231.3108169, -11753.2741318}, 
{-655.005025082, -650.228279807, -643.153048806, -633.930990567, -622.713763576, -609.65302632, -594.900437285, -578.607654957, -560.926337823, -542.008144369, -522.004733082, -501.067762449, -479.348890955, -456.999777087, -434.44622982, -411.666479041, -388.618002653, -365.396383927, -342.097206131, -318.816052536, -295.648506411, -272.707410633, -250.068705196, -227.825589699, -206.071263742, -184.898926925, -164.401778848, -144.67741041, -125.822752667, -107.920505117, -91.0620843543, -75.3389069681, -60.8423895501, -47.6639486914, -36.2708462931, -26.5167378317, -18.2268289497, -11.4156601864, -6.09777208091, -2.28770517239, 0, 0.899018951533, 0.226125435952, -2.05369073856, -5.97543976381, -11.5741318316, -18.8847771338, -26.3284314197, -32.5686664064, -41.4996416908, -53.8116848306, -70.1951233836, -91.3402849074, -117.93749696, -147.032956277, -181.622280873, -224.97693341, -278.532618287, -343.725039907, -421.989902669, -514.762910975, -581.305807653, -694.75647468, -854.432254023, -1059.65048765, -1330.79572679, -1615.2196654, -1937.52009431, -2322.29480437, -2806.55381688, -3377.65823128, -4109.87455604, -4751.12504332, -5083.61292352, -5376.10301522, -5691.28313547, -6024.29668111, -6414.84461751, -6791.94205615, -7122.338087, -7372.78180002}
};

  const double a5tab[21][81] = {
{519.52389259, 518.131238614, 515.84649954, 512.673834129, 508.617401144, 503.681359348, 497.869867501, 491.187084367, 483.637168707, 475.224279284, 465.95257486, 455.826214197, 444.849356057, 433.026159203, 419.876634404, 405.71121421, 390.876979038, 375.477116122, 359.614812695, 343.393255989, 326.915633237, 310.979210278, 294.897238728, 278.677048807, 262.325970739, 245.851334744, 229.260471043, 212.468356931, 195.403846833, 178.295097371, 161.186937022, 144.124194265, 127.151697578, 110.314275439, 94.0951208275, 78.2617512383, 62.5485019061, 46.9105374736, 31.3030225833, 15.6811218779, 0, -14.8448983333, -29.9685476749, -45.545642503, -61.7508772958, -78.7589465317, -96.7445446888, -115.834202855, -136.161689243, -158.0168625, -181.593972988, -207.087271069, -234.691007106, -264.599431462, -315.508826205, -375.908988916, -432.889194904, -482.859256853, -522.228987449, -547.408199377, -554.806705321, -574.962905685, -596.142263911, -541.649647251, -334.789922961, 439.372043702, 1023.20821195, 1899.29846984, 3550.22270543, 6493.90189401, 11106.892662, 19378.4641475, 26907.6180856, 30293.5483768, 32120.0118328, 32707.9753079, 31765.1304306, 28865.1245932, 24731.0573588, 19694.1678859, 14085.6953332}, 
{401.023924002, 404.584454518, 406.799775776, 407.720454015, 407.39705547, 405.88014638, 403.22029298, 399.468061509, 394.674018203, 388.888729298, 382.162761033, 374.546679643, 366.091051367, 356.84644244, 346.759978418, 335.94766249, 324.53332785, 312.588698674, 300.185499134, 287.395453404, 274.290285658, 261.318543775, 248.123086287, 234.723595435, 221.139753455, 207.391242587, 193.497745068, 179.435613896, 165.187712805, 150.877329198, 136.541738377, 122.218215645, 107.944036305, 93.7564756591, 80.0563750461, 66.6510165674, 53.3201619776, 40.0267221513, 26.7336079632, 13.4037302879, 0, -12.7748555162, -25.7999156246, -39.2144431797, -53.1577010357, -67.768952047, -83.187459068, -100.053583656, -118.932391291, -138.764960599, -159.487092257, -181.034586946, -203.343245343, -226.348868128, -271.660391067, -325.503110484, -372.461761881, -408.03908527, -427.737820663, -427.060708071, -401.510487508, -366.642494705, -338.585770163, -231.772113862, 39.3666742213, 733.614348879, 1509.27518361, 2659.77438129, 4478.53714483, 7303.46344497, 11294.554181, 17504.0752039, 23491.4724031, 27390.009096, 30357.5261523, 32364.0939479, 33034.5068354, 31900.5993, 29211.321319, 24971.1506125, 19184.5649005}, 
{320.853928082, 323.645506111, 325.647364959, 326.867099617, 327.312305078, 326.990576332, 325.909508371, 324.076696185, 321.499734768, 318.186219109, 314.143744201, 309.379905034, 303.902296601, 297.718513892, 290.663545466, 282.854177775, 274.420262814, 265.404700795, 255.850391928, 245.800236423, 235.297134492, 224.310145034, 212.966207561, 201.318420276, 189.41988138, 177.323689075, 165.082941563, 152.819471579, 140.644783673, 128.447623925, 116.253182108, 104.086647993, 91.9732113538, 79.9380619615, 68.2956295161, 56.888129556, 45.5358834317, 34.2049192785, 22.8612652319, 11.4709494273, 0, -10.9706147872, -22.1747343712, -33.7312580616, -45.759085168, -58.3771149998, -71.7042468668, -87.066033482, -105.606713446, -124.559939965, -143.554672294, -162.219869686, -180.184491394, -197.077496672, -236.118910276, -282.013897325, -317.681883233, -337.926469772, -337.551258712, -311.359851825, -254.155850881, -10.8720223816, 90.9067418446, 180.87168694, 388.714058047, 763.065479413, 1633.56374228, 3032.62932487, 4992.68270543, 7505.62609173, 10725.4347734, 15243.7118833, 19832.0916707, 23594.4220564, 26919.4490701, 29614.2054123, 31277.3477455, 31430.3292581, 30111.7930734, 27184.0667948, 22509.4780255}, 
{265.757057868, 268.992014755, 271.374000094, 272.928727155, 273.681909204, 273.659259508, 272.886491334, 271.38931795, 269.193452623, 266.32460862, 262.808499208, 258.670837655, 253.937337227, 248.633711192, 242.742431083, 236.316565045, 229.412454434, 222.064657256, 214.307731519, 206.176235228, 197.70472639, 188.982828262, 179.982428712, 170.730480859, 161.253937821, 151.579752716, 141.734878663, 131.755229673, 121.675373277, 111.500836074, 101.254932771, 90.9609780769, 80.642286699, 70.3221733452, 60.2526836524, 50.3124363079, 40.3627360757, 30.3801126421, 20.341095694, 10.2222149177, 0, -10.2644562319, -20.6908659813, -31.3143783106, -42.1701422825, -53.2933069594, -64.719021404, -77.9982493207, -94.4541851258, -110.497491948, -125.547849398, -139.024937083, -150.348434614, -158.938021599, -191.199015796, -229.479883636, -254.086397652, -258.918548191, -237.876325597, -184.859720216, -93.7687223946, 35.5117657919, 180.291277618, 404.817816634, 773.339386386, 1386.7392502, 2218.89867197, 3378.02761529, 4972.33604374, 7164.73153685, 9899.33121033, 13328.2449008, 16865.3756879, 20042.6392326, 23194.4795398, 26114.4381133, 28323.8597239, 29265.975914, 28968.6678726, 27260.5035394, 23970.0508542}, 
{220.583032734, 225.133541135, 228.681820835, 231.273563387, 232.954460346, 233.770203268, 233.766483707, 232.988993218, 231.483423357, 229.295465679, 226.470811737, 223.055153088, 219.094181286, 214.633587886, 209.892492842, 204.80687627, 199.299590523, 193.38085381, 187.060884339, 180.34990032, 173.258119962, 165.52221022, 157.46371986, 149.130646394, 140.570987334, 131.832740192, 122.963902479, 114.111941205, 105.409376897, 96.6663641382, 87.8905126268, 79.0894320625, 70.2707321446, 61.4420225726, 52.6863692597, 43.9636480289, 35.2280230006, 26.4716699291, 17.6867645689, 8.86548267436, 0, -8.92613737653, -17.9109322065, -26.9610169181, -36.0830239393, -45.2835856983, -54.569334623, -64.9752994756, -77.3819801629, -89.3369960233, -100.429417527, -110.248315145, -118.382759346, -124.421820602, -159.42287808, -203.068035322, -232.657383567, -241.343410808, -222.278605036, -168.615454246, -73.5064464281, 82.2729673551, 246.276045279, 483.34908777, 858.338395253, 1350.4448778, 2235.77291725, 3476.3951678, 5034.38428363, 6820.92925718, 8950.75372792, 11720.2504834, 14611.909575, 17240.3546917, 20208.6474594, 23267.9402507, 25685.5804587, 26355.8905134, 26149.0976754, 25261.8223463, 23890.6849277}, 
{208.61921005, 209.915347094, 210.713859519, 211.024196795, 210.855808392, 210.218143781, 209.120652433, 207.572783816, 205.583987403, 203.163712662, 200.321409065, 197.066526081, 193.408513182, 189.356819836, 184.833331841, 179.902891767, 174.636970368, 169.06292754, 163.20812318, 157.099917184, 150.765669448, 144.339762451, 137.727752954, 130.942220304, 123.995743843, 116.900902917, 109.670276871, 102.326656576, 94.8912985012, 87.3523650613, 79.7182893921, 71.9975046296, 64.1984439096, 56.3295403682, 48.4647341476, 40.5710184124, 32.6104277358, 24.5779963371, 16.4687584354, 8.27774824993, 0, -8.48508618694, -17.0508721395, -25.6863537782, -34.3805270234, -43.1223877955, -51.9009320148, -62.9740331917, -78.2586388286, -92.3185993958, -104.221673605, -113.035620167, -117.828197794, -117.667165198, -140.008456596, -165.96136727, -174.486332707, -158.844544931, -112.297195964, -28.1054778299, 100.46941745, 287.246173884, 500.723222361, 787.30024641, 1193.37692956, 1695.04967351, 2512.1324737, 3606.66175036, 4940.67392375, 6468.70711687, 8175.29264168, 9937.16077603, 12049.3921511, 14789.0148884, 17536.7183287, 20063.4979192, 22273.316662, 24152.955124, 25224.9560583, 25170.1257998, 23669.2706833}, 
{201.476866176, 200.830761019, 199.966682398, 198.875303392, 197.547297078, 195.973336536, 194.144094843, 192.050245077, 189.682460317, 187.03141364, 184.087778125, 180.84222685, 177.285432893, 173.408069333, 168.993149371, 164.162712481, 159.054517075, 153.701711341, 148.13744347, 142.39486165, 136.507114072, 130.758902677, 124.897080605, 118.92005475, 112.826232006, 106.614019265, 100.28182342, 93.8125967298, 87.1916136977, 80.4542350115, 73.6051426422, 66.6490185607, 59.5905447382, 52.4344031455, 45.2237038034, 37.9388189585, 30.5572120232, 23.0757048287, 15.4911192059, 7.80027698603, 0, -8.19377298947, -16.4645452067, -24.7767029443, -33.094632495, -41.3827201514, -49.6053522063, -60.6757125119, -77.0638928014, -91.6833676142, -103.301216678, -110.684519721, -112.600356471, -107.815806656, -121.134044267, -134.851620781, -129.294270065, -97.9036039359, -34.1212342132, 68.611227285, 216.85216874, 429.644449025, 672.520768152, 982.366869827, 1396.06849776, 1871.91357258, 2640.66196524, 3644.88278628, 4827.14514625, 6130.18230973, 7496.07092516, 8614.35078945, 10184.990503, 12818.2125978, 15321.765809, 17483.4945034, 19514.7942972, 21906.0560746, 23608.1383984, 24096.7426856, 22847.5703534}, 
{193.412758673, 193.027177063, 192.317482484, 191.28812499, 189.943554636, 188.288221475, 186.326575562, 184.063066952, 181.502145697, 178.648261854, 175.505865475, 172.079406616, 168.37333533, 164.392101672, 160.04959652, 155.407557993, 150.534579142, 145.453633152, 140.187693209, 134.759732498, 129.192724205, 123.644311847, 117.984199382, 112.2167611, 106.346371289, 100.37740424, 94.3142342432, 88.177761346, 81.9864024045, 75.7050167314, 69.3312686272, 62.8628223926, 56.297342328, 49.6324927341, 42.8730634985, 36.0122112495, 29.0425532394, 21.9602962888, 14.7616472183, 7.44281284856, 0, -8.28948217435, -16.6112444927, -24.869795441, -32.969643505, -40.8152971705, -48.3112649234, -58.2559074825, -73.0127496402, -85.5667686808, -94.6474760933, -98.984383367, -97.3070019913, -88.3448434552, -96.0719154922, -103.247969136, -90.7218582927, -52.0595413758, 19.1730232021, 129.409877028, 285.08506169, 458.251570389, 706.120238051, 1044.91253722, 1490.84994043, 2019.54915205, 2747.38981293, 3641.86806361, 4670.48004459, 5819.51712876, 7000.08975955, 7893.9715777, 9132.50272843, 11229.5039323, 13465.2015364, 15643.3770615, 17739.9249003, 19837.8530978, 21419.602218, 22172.3526321, 21783.2847115}, 
{191.601359398, 191.648579923, 191.221035197, 190.339554318, 189.024966385, 187.298100499, 185.179785759, 182.690851264, 179.852126114, 176.684439408, 173.208620246, 169.445497728, 165.415900952, 161.140659018, 156.675988439, 152.020332108, 147.169454246, 142.13694575, 136.936397516, 131.581400439, 126.085545418, 120.469396301, 114.738608019, 108.905808456, 102.983625497, 96.9846870237, 90.9216209216, 84.8600608487, 78.8576757115, 72.8003505332, 66.6791911499, 60.4853033978, 54.209793113, 47.8437661318, 41.3546690887, 34.7485747172, 28.0333466479, 21.2049300113, 14.2592699376, 7.19231155702, 0, -8.48424425866, -16.9413988129, -25.2149659123, -33.1484478061, -40.585346744, -47.3691649753, -56.1187592943, -69.0359562647, -79.3280684403, -85.7117149429, -86.9035148947, -81.6200874176, -68.5780516336, -70.8642320785, -71.7785669886, -52.7762774367, -7.58925843224, 70.0505950155, 186.411387897, 347.761225203, 481.252448158, 730.754880425, 1093.9234584, 1568.41311849, 2153.44676898, 2842.81169228, 3636.04437602, 4532.68130783, 5571.06200878, 6634.31386621, 7426.76738499, 8434.45780506, 10012.3901728, 11986.6611226, 14180.3275732, 16308.2407799, 18002.4165024, 19334.4181456, 20225.9626724, 20598.7670457}, 
{198.085243192, 197.864701648, 197.124898376, 195.892428711, 194.193887989, 192.055871547, 189.504974719, 186.567792842, 183.270921252, 179.640955283, 175.704490273, 171.488121556, 167.018444469, 162.322054347, 157.494003788, 152.517581815, 147.370895535, 142.066537918, 136.617101932, 131.035180546, 125.333366728, 119.506471894, 113.58732632, 107.590978729, 101.532477844, 95.4268723871, 89.289211082, 83.2045589259, 77.2474601333, 71.2655465179, 65.2454379453, 59.1737542811, 53.0371153909, 46.8221411402, 40.4574148969, 33.9662706269, 27.376435403, 20.6863999769, 13.8946551, 6.99969152385, 0, -8.57763657726, -17.0597666603, -25.2446465588, -32.9305325824, -39.9156810407, -45.9983482433, -53.8635192418, -65.7624135512, -74.59078936, -78.9747983476, -77.5405921932, -68.9143225761, -51.7221411755, -47.5526023711, -40.5057580908, -12.9435075615, 41.2047672398, 128.009684336, 253.54186175, 423.871917504, 549.93267103, 808.527753602, 1187.47275984, 1674.58328438, 2281.20823494, 2937.11451778, 3658.35747732, 4460.992458, 5413.1625911, 6374.65986045, 7035.4556386, 7854.09403056, 9145.57429634, 10886.517868, 12922.0963926, 14924.6585679, 16436.9310183, 17692.9516097, 18691.0828929, 19429.6874184}, 
{208.515273621, 207.895394079, 206.73117931, 205.052008192, 202.887259608, 200.266312438, 197.218545563, 193.773337865, 189.960068223, 185.80811552, 181.346858636, 176.605676451, 171.613947848, 166.401051707, 161.061008648, 155.582304953, 149.948533806, 144.175852155, 138.280416948, 132.278385132, 126.185913654, 120.022885944, 113.801217814, 107.536551558, 101.244529472, 94.9407938489, 88.6409869843, 82.4406798579, 76.4234332082, 70.4137707142, 64.3948810557, 58.3499529126, 52.2621749648, 46.1147358922, 39.7796836381, 33.3105149178, 26.7691393111, 20.1614784169, 13.4934538342, 6.77098716227, 0, -8.4782748928, -16.7633221431, -24.619315217, -31.8104275809, -38.1008327008, -43.2547040432, -50.0869586774, -60.9541029408, -68.3256339127, -70.7270251248, -66.683750109, -54.721282397, -33.3650955206, -23.1959971881, -8.78717603874, 26.9576708996, 90.0243064272, 186.398493344, 322.065994451, 503.012572547, 644.710459703, 917.932133632, 1307.76711159, 1799.30491083, 2400.61028018, 3040.10064607, 3730.43558307, 4484.27466575, 5366.68562445, 6233.1035665, 6756.94413795, 7429.54030234, 8594.31013232, 10141.1541134, 11933.7843154, 13724.2665973, 15157.5506057, 16418.6970474, 17508.6370555, 18428.3017627}, 
{221.690130768, 220.398246625, 218.586560333, 216.281615915, 213.509957393, 210.298128789, 206.672674127, 202.660137429, 198.287062719, 193.579994018, 188.565475349, 183.270050736, 177.7202642, 171.942659766, 165.966042024, 159.815524955, 153.516051419, 147.093703058, 140.574561515, 133.98470843, 127.350225447, 120.795948539, 114.235566364, 107.681521916, 101.146258186, 94.6422181653, 88.181844846, 81.854798501, 75.7391365828, 69.6626667447, 63.6064792811, 57.5516644863, 51.4793126547, 45.3705140807, 39.0136819752, 32.5117948587, 25.9824160579, 19.4460463723, 12.9231866011, 6.43433754388, 0, -8.04107527006, -15.7543762247, -22.8871408612, -29.1866071769, -34.400013169, -38.274596835, -43.8795544398, -53.7849182675, -59.7947446899, -60.3074499009, -53.7214500945, -38.4351614648, -12.8470002057, 2.69274262492, 23.672486987, 67.2423130375, 139.493879167, 246.518843766, 394.408865224, 589.255601934, 780.905090651, 1078.63778946, 1474.3762113, 1960.04286909, 2514.26833182, 3161.76182244, 3878.49564874, 4640.44211851, 5456.18109363, 6203.86221957, 6522.65531823, 7063.54625334, 8330.47500224, 9738.76736699, 11163.6008056, 12630.3004549, 14190.7865456, 15640.4512534, 16862.1398814, 17738.697733}, 
{284.283845877, 283.384498341, 281.788474232, 279.537192289, 276.672071249, 273.234529852, 269.265986836, 264.807860939, 259.9015709, 254.588535456, 248.910173348, 242.907903312, 236.623144088, 230.097314414, 223.712012113, 217.293457393, 210.642120404, 203.729839064, 196.528451291, 189.009795003, 181.145708119, 172.183119805, 162.918891595, 153.424976274, 143.773326623, 134.035895427, 124.284635468, 114.812673441, 105.879902013, 97.0294218855, 88.2433821283, 79.5039318099, 70.7932199986, 62.0933957631, 53.3972560408, 44.6802140197, 35.9168768908, 27.0872157894, 18.1712018504, 9.14880620896, 0, -9.9799139608, -20.0317388329, -30.0809460952, -40.0530072265, -49.8733937058, -59.4675770119, -71.817505352, -89.4458547306, -104.969711683, -117.07351992, -124.441723153, -125.758765095, -119.709089455, -116.532694541, -107.603938055, -83.4224731341, -40.3091502973, 25.4151799369, 117.42966705, 239.413460523, 247.334438304, 414.645739477, 709.643707954, 1100.62468764, 1592.68020266, 2063.3453879, 2525.07045029, 2990.30559674, 3528.13496043, 3981.10696954, 3894.09243459, 4099.75603874, 5251.57161774, 6482.59405139, 7658.4287082, 8926.61725621, 10549.1951049, 12048.2863985, 13198.2588585, 13773.4802068}, 
{324.726459799, 323.512901482, 321.466687699, 318.641532775, 315.091151033, 310.869256799, 306.029564395, 300.625788146, 294.711642376, 288.340841409, 281.567099569, 274.44413118, 267.025650567, 259.365372053, 252.016226892, 244.716122955, 237.165179595, 229.315000439, 221.117189114, 212.523349246, 203.485084464, 192.831974261, 181.792605618, 170.473541382, 158.981344401, 147.422577522, 135.903803594, 124.851731116, 114.644962875, 104.624556897, 94.7670866063, 85.0491254281, 75.4472467874, 65.938024109, 56.5112876625, 47.1352245446, 37.7770183671, 28.410530977, 19.0096242212, 9.54815994655, 0, -10.0460190514, -20.1778354517, -30.3684127246, -40.5907143941, -50.8177039839, -61.0223450179, -74.1335870573, -92.6362073893, -109.435069313, -123.302911035, -133.012470762, -137.336486699, -135.047697052, -131.911383725, -122.276776261, -99.9637994581, -62.314924498, -6.67262256244, 69.6206351666, 169.222377507, 114.370792549, 237.099483708, 497.050219728, 853.864769351, 1322.70419852, 1726.26303154, 2090.80566081, 2442.59647874, 2872.82843337, 3212.98025014, 2994.89972292, 3060.26065094, 4067.81540286, 5119.08086806, 6095.25355033, 7189.80406456, 8750.7005116, 10195.6469372, 11259.4842445, 11677.0533363}, 
{349.285236663, 347.293610919, 344.358226682, 340.543292462, 335.913016769, 330.531608113, 324.463275004, 317.772225953, 310.522669469, 302.778814063, 294.604868245, 286.065040526, 277.223539414, 268.144573421, 259.41513866, 250.768726297, 241.899254775, 232.764000722, 223.320240763, 213.525251524, 203.336309632, 191.523430235, 179.395120382, 167.072625643, 154.677191589, 142.330063793, 130.152487823, 118.642805237, 108.242694403, 98.1717213272, 88.3980184938, 78.889718387, 69.6149534907, 60.5418562891, 51.5233513761, 42.6004521678, 33.8228928097, 25.1823706078, 16.6705828683, 8.279226897, 0, -8.74739247548, -17.320265256, -25.6479250271, -33.6596784739, -41.2848322819, -48.4526931364, -58.0463925353, -72.5052144524, -84.6955330841, -93.3473073247, -97.1904960684, -94.9550582096, -85.3709526427, -73.8185643003, -54.8207682937, -22.3970317877, 26.0829760078, 93.2495858829, 181.733128628, 294.163935033, 263.69608624, 402.396955947, 672.431934703, 1035.96641306, 1506.47879281, 1919.56269848, 2298.95864152, 2668.40713339, 3113.03675472, 3472.42380943, 3326.98010234, 3399.00851661, 4245.93318796, 5120.8825451, 5938.14966868, 6880.08091883, 8286.84746309, 9599.26795989, 10564.7598518, 10930.7405815}, 
{364.2274406, 361.236781526, 357.200117129, 352.191126438, 346.283488478, 339.550882278, 332.066986864, 323.905481263, 315.140044502, 305.844355609, 296.09209361, 285.956937533, 275.512566404, 264.832659251, 254.445199718, 244.136472073, 233.658371726, 222.991653451, 212.117072019, 201.015382202, 189.667338773, 177.038406055, 164.264847836, 151.467637457, 138.767748257, 126.286153578, 114.143826758, 102.856729512, 92.881471827, 83.3945647238, 74.3566035699, 65.7281837322, 57.469900578, 49.5423494744, 41.6010217175, 33.7995226233, 26.3145525412, 19.1691133545, 12.3862069466, 5.98883520064, 0, -6.59203616645, -12.5866274209, -17.8178672799, -22.1198492601, -25.3266668781, -27.2724136504, -30.7745147654, -38.2021141289, -42.2558184978, -41.5583928644, -34.7326022212, -20.4012115609, 2.81301412421, 27.3705243565, 60.2893451838, 109.262464024, 177.490970663, 268.175954888, 384.518506486, 529.719715244, 576.99225223, 766.93608999, 1071.62504097, 1463.1329176, 1946.14245895, 2412.29233384, 2872.78675361, 3338.82992959, 3871.53516672, 4332.37939548, 4392.23649196, 4569.94806241, 5252.65180337, 5960.65381066, 6651.19139989, 7457.66785792, 8649.18065681, 9775.88455657, 10633.0294931, 11015.8654022}, 
{375.820335738, 371.852568172, 366.729384988, 360.533689791, 353.348386186, 345.256377779, 336.340568175, 326.683860979, 316.369159798, 305.479368235, 294.097389897, 282.306128388, 270.188487315, 257.827370282, 245.642862364, 233.504564937, 221.256556229, 208.912772163, 196.48714866, 183.993621644, 171.446127036, 158.157820047, 144.94019993, 131.903985226, 119.159894475, 106.818646219, 94.9909589975, 84.1643376538, 74.7696703799, 66.0167366348, 57.8632676133, 50.2669945098, 43.1856485191, 36.5769608356, 29.9118732225, 23.4560616454, 17.5120498846, 12.1371378899, 7.38862561102, 3.32381299781, 0, -4.08795205786, -7.10452112122, -8.77662376076, -8.83117654716, -6.99509605107, -2.99529884318, 0.46345327266, 1.1238553722, 6.37935894314, 17.7121482685, 36.6044076314, 64.5383213147, 102.996073602, 141.280642871, 188.578823509, 254.999322123, 344.515478911, 461.100634069, 608.728127794, 791.371300284, 935.941223372, 1187.11481963, 1530.4657266, 1951.56758183, 2443.83367034, 2973.49988274, 3535.54735828, 4124.95723622, 4771.09891165, 5365.78875633, 5692.57181091, 6027.027715, 6554.69807946, 7111.04939284, 7698.4530806, 8382.78492071, 9329.24479017, 10242.2318175, 10983.2369808, 11413.7512584}, 
{390.331186208, 385.651125731, 379.683056206, 372.519637609, 364.253529915, 354.977393099, 344.783887137, 333.765672005, 322.015407677, 309.625754129, 296.689371337, 283.298919276, 269.547057922, 255.526447249, 241.544578897, 227.558209545, 213.507834067, 199.442170397, 185.409936467, 171.459850214, 157.640629571, 143.66259054, 129.959588613, 116.62707735, 103.760510311, 91.4553410556, 79.8070231433, 69.2364633713, 60.1156652933, 51.761886608, 44.138436403, 37.2086237656, 30.9357577833, 25.2831475435, 19.623480427, 14.2936949682, 9.67543716262, 5.85282288658, 2.90996801651, 0.930988428844, 0, -1.74314208322, -2.00154553197, -0.422578747037, 3.3463898708, 9.65799192076, 18.8648590021, 28.4489185994, 36.3234558421, 49.705283736, 70.112631997, 99.0637303411, 138.076808484, 188.670096142, 237.536551867, 295.572926018, 374.798176657, 479.762920194, 615.017773039, 785.113351604, 994.600272298, 1222.22493252, 1519.33107866, 1884.7901797, 2317.47370459, 2801.69090041, 3372.23329029, 4010.49781672, 4697.88142217, 5434.50323175, 6145.59364005, 6729.88897832, 7224.19590103, 7618.79884656, 8044.72401978, 8544.00904745, 9115.6521461, 9818.58456057, 10515.0448327, 11134.3261275, 11605.7216099}, 
{414.027256137, 409.142609073, 402.798156731, 395.097624981, 386.144739689, 376.043226723, 364.896811952, 352.809221243, 339.884180463, 326.225415482, 311.936652166, 297.121616383, 281.884034001, 266.327630888, 250.686801615, 234.982610552, 219.226231022, 203.49466169, 187.86490122, 172.413948277, 157.218801526, 142.333635861, 127.861425835, 113.882322232, 100.476475835, 87.7240374288, 75.7051577958, 64.7439403749, 55.1278317991, 46.3536641913, 38.4025357179, 31.2555445455, 24.8937888404, 19.2983667691, 13.9034178667, 9.036048326, 5.06476669816, 2.08254701741, 0.182363317966, -0.542810365915, 0, -0.0656081760354, 1.5947001719, 5.34588348375, 11.5529001994, 20.5807087589, 32.7942676021, 45.9632882353, 58.2474490718, 76.2172403782, 101.291374244, 134.88856276, 178.427518014, 233.326952098, 285.76301197, 346.79691205, 428.643661773, 535.839713956, 672.921521416, 844.425536969, 1054.88821343, 1317.52531253, 1619.98280088, 1970.43458833, 2377.05458472, 2821.85262255, 3377.54050162, 4020.89549014, 4728.69485634, 5484.52336928, 6244.7357947, 7006.09091329, 7615.40104718, 7911.68093504, 8234.33241959, 8651.93363709, 9116.48957297, 9608.74466543, 10111.0586925, 10605.2407457, 11073.0999166}, 
{453.175809657, 448.83717307, 442.811712512, 435.216306994, 426.167835531, 415.783177135, 404.179210819, 391.472815596, 377.780870479, 363.22025448, 347.907846614, 331.960525893, 315.495171329, 298.628661936, 281.60598282, 264.462972611, 247.225772875, 229.98505958, 212.831508698, 195.855796198, 179.148598052, 162.951874337, 147.184123545, 131.915128272, 117.214671117, 103.152534679, 89.7985015557, 77.357602375, 66.0145451294, 55.5157189325, 45.875991337, 37.1102298953, 29.2333021598, 22.2600756831, 15.9192600775, 10.4067474529, 5.94009081402, 2.59268895508, 0.437940670296, -0.450755246132, 0, 0.436647730172, 2.55661681543, 6.6303786541, 12.9284046445, 21.721166185, 33.2791346738, 45.7879692008, 57.7465968524, 74.4105133671, 96.8966909333, 126.322101739, 163.803717973, 210.458511822, 255.584783805, 307.77604094, 376.520411617, 465.35227964, 577.806028813, 717.416042941, 887.716705828, 1103.52429625, 1345.46792008, 1623.2351406, 1946.51352106, 2306.45731019, 2758.46946184, 3289.99773973, 3888.48990761, 4543.93456651, 5236.15696833, 6023.08053496, 6654.59158009, 6900.07117523, 7152.52932041, 7486.30118616, 7845.51724021, 8191.26980215, 8547.00848695, 8914.92464792, 9297.20963838}, 
{514.044110896, 511.244972594, 506.460749494, 499.824338737, 491.468637465, 481.526542818, 470.130951938, 457.414761967, 443.510870045, 428.552173315, 412.671568916, 396.001953991, 378.676225682, 360.827281128, 342.838574808, 324.684500378, 306.320485408, 287.828177606, 269.289224681, 250.785274342, 232.397974298, 214.298224298, 196.466093692, 178.970903871, 161.881976225, 145.268632146, 129.200193024, 113.748283082, 98.9841805159, 84.9717003796, 71.7792290391, 59.4751528609, 48.1278582112, 37.8057314563, 28.8385815953, 21.1294180831, 14.561461833, 9.14962737234, 4.90882922841, 1.85398192853, 0, -0.744376298108, -0.243394776343, 1.53252248652, 4.6129534117, 9.02747592042, 14.8056679339, 20.7043685163, 25.6716609749, 32.7803872002, 42.5768979871, 55.6075441306, 72.4186764256, 93.5566456671, 116.626627996, 144.035572027, 178.413060337, 220.907036689, 272.665444847, 334.836228573, 408.567331632, 461.903816535, 552.184370064, 679.028024583, 842.05381246, 1057.64343674, 1284.06811606, 1541.06192668, 1848.35894487, 2235.51206568, 2692.798909, 3282.76076244, 3795.71592643, 4050.69639749, 4271.96945035, 4511.1860313, 4762.95518671, 5057.70466815, 5339.62930629, 5582.32164675, 5759.37423514}
};

const double b3tab[21][81] = {
{43.006888187, 42.0198431635, 41.0202058237, 40.0082996009, 38.9844479288, 37.9489742407, 36.90220197, 35.8444545503, 34.7760554148, 33.6973279971, 32.6085957307, 31.5101820489, 30.4024103852, 29.285604173, 28.1351348106, 26.96711048, 25.7995281222, 24.6378149032, 23.4873979892, 22.353704546, 21.2421617399, 20.2361251944, 19.2523311641, 18.2854443608, 17.3301294967, 16.3810512839, 15.4328744343, 14.4741565671, 13.4943729651, 12.5027910972, 11.4965553517, 10.4728101165, 9.42869977973, 8.36136872951, 7.43944692215, 6.55159634942, 5.57349830121, 4.46722120984, 3.19483350763, 1.71840362691, 0, -1.82864689055, -3.99856215427, -6.57110885044, -9.60765003835, -13.1695487773, -17.3181681266, -23.1599594087, -31.6443369287, -40.3337383472, -48.8651873061, -56.8757074473, -64.0023224128, -69.8820558445, -135.302015635, -221.215458956, -283.946851061, -310.625471226, -288.380598727, -204.341512841, -45.6374928453, 1094.39585398, 1322.02807435, 1275.54747001, 1593.24234271, 5319.217858, 7157.42945067, 10633.1225634, 19271.5426388, 33460.8266451, 66137.5454483, 158302.586314, 233321.855053, 232304.720441, 214245.06634, 191781.420917, 164045.903301, 133506.452376, 102794.137042, 74484.4134509, 51152.7377553}, 
{34.8479771104, 34.0828232354, 33.2979234396, 32.4946199049, 31.6742548131, 30.8381703461, 29.9877086858, 29.1242120141, 28.2490225128, 27.3634823638, 26.468933749, 25.5667188503, 24.6581798495, 23.7446589285, 22.8204984398, 21.8914686829, 20.9638698135, 20.0404757708, 19.1240604938, 18.2173979217, 17.3232619936, 16.4868399804, 15.6626349178, 14.8475631731, 14.0385411136, 13.2324851069, 12.4263115202, 11.6119552641, 10.7820997724, 9.94557264507, 9.1013128811, 8.24825947944, 7.38535143907, 6.51152775896, 5.87655892888, 5.32070703129, 4.66524731014, 3.85781325592, 2.84603835915, 1.57755611034, 0, -1.69213163609, -3.83216387926, -6.5065569656, -9.80177113118, -13.8042666121, -18.6005036443, -25.4550676649, -35.4555167015, -45.8712790792, -56.3104583448, -66.3811580455, -75.6914817284, -83.8495329406, -163.679531509, -268.472798615, -345.978518061, -380.829041893, -357.656722162, -261.093910914, -75.7729601991, 1382.32394614, 1669.12289502, 1567.52845544, 1860.44519637, 5398.05478386, 7863.99151125, 12521.8726183, 22635.3153448, 38866.4163056, 72283.354615, 159204.733468, 230487.863001, 232792.772626, 220036.698221, 203067.02445, 179727.632588, 149951.727682, 118318.723077, 87824.0905572, 61463.3019091}, 
{27.3016278376, 26.3201609262, 25.4316052792, 24.6262981821, 23.8945769209, 23.2267787813, 22.6132410491, 22.0443010101, 21.5102959501, 21.0015631548, 20.5084399101, 20.0212635017, 19.5303712154, 19.026100337, 18.4401786469, 17.8000200179, 17.1374748041, 16.4548683812, 15.7545261249, 15.038773411, 14.3099356152, 13.5774013959, 12.8354573584, 12.0854533908, 11.3287393811, 10.5666652173, 9.80058078738, 9.02693396133, 8.24290919596, 7.46157766436, 6.6862796314, 5.92035536198, 5.16714512095, 4.42998917319, 4.08113487694, 3.89055052347, 3.60027976262, 3.13820596214, 2.43221248972, 1.41018271311, 0, -1.63570072198, -3.83620844389, -6.70606059688, -10.3497946121, -14.8719479208, -20.377057954, -28.3310150951, -39.9951501561, -52.2188501768, -64.5538995445, -76.5520826465, -87.7651838701, -97.7449876027, -189.784156758, -310.45969035, -400.009771887, -440.85768112, -415.426697802, -306.140101684, -95.421172519, 1828.26795622, 2199.66408917, 1975.70095431, 2113.31227959, 4971.62858844, 8048.84048732, 13984.5043978, 25418.1767413, 42730.9598168, 75337.772414, 154867.655701, 220221.658437, 224088.909956, 215027.203614, 202272.430773, 182777.57841, 154477.671148, 123662.487418, 93773.5816982, 68252.5084653}, 
{19.2586377425, 18.8003073152, 18.3462394291, 17.8955830901, 17.4474873038, 17.0011010763, 16.5555734132, 16.1100533203, 15.6636898036, 15.2156318687, 14.7650285216, 14.3110287681, 13.8527816138, 13.3894360648, 12.9095360277, 12.4189393328, 11.9243065695, 11.4269559288, 10.9282056019, 10.42937378, 9.93177865431, 9.45331350472, 8.97643230028, 8.50016409882, 8.02353795816, 7.54558293613, 7.06532809055, 6.57811646971, 6.07984498892, 5.57835637584, 5.07417632661, 4.5678305374, 4.05984470435, 3.55074452361, 3.43289261098, 3.45893719721, 3.35186618484, 3.03205826641, 2.41989213443, 1.43574648145, 0, -1.66411453275, -3.98089873307, -7.07180013261, -11.058266263, -16.0617446559, -22.2036828429, -31.1912447578, -44.49332147, -58.4397314858, -72.5080705226, -86.1759342979, -98.9209185293, -110.220618934, -213.907692349, -349.770378471, -450.454654645, -496.138568622, -467.000168155, -343.217500997, -104.9686149, 1568.99792525, 1952.72649168, 1931.13169949, 2389.12816383, 5878.89254073, 9172.77041475, 15156.3748445, 26715.3188887, 45517.3206734, 78101.6780551, 149670.498386, 208686.902733, 215022.654843, 206287.0998, 190983.590201, 169878.757991, 146471.669021, 121068.374164, 95063.4568778, 69851.5006174}, 
{12.9663745266, 12.5937778428, 12.2387102972, 11.8997403971, 11.5754366499, 11.2643675629, 10.9651016434, 10.6762073988, 10.3962533363, 10.1238079634, 9.85743978732, 9.59571731539, 9.33720905496, 9.08048351336, 8.84361875511, 8.61284147046, 8.37290156214, 8.11837701879, 7.84384582904, 7.54388598155, 7.21307546493, 6.75455496135, 6.26696788446, 5.75752034143, 5.23341843943, 4.70186828563, 4.17007598721, 3.65135689791, 3.15810838407, 2.68292864655, 2.23054324209, 1.80567772744, 1.41305765934, 1.05740859454, 1.26135689085, 1.70200216816, 2.01724123779, 2.10586732777, 1.86667366611, 1.19845348085, 0, -1.56438303578, -3.89808173212, -7.13897168199, -11.4249284784, -16.8938277143, -23.6835449826, -33.8384723961, -49.1165251466, -65.0968800814, -81.1433058232, -96.6195709948, -110.889444219, -123.316694118, -235.579017799, -382.316019323, -490.422589208, -538.335027566, -504.489634512, -367.322710159, -105.270554619, 1722.27467819, 2170.01091914, 2162.19015615, 2623.06437714, 5529.95243652, 9209.54812445, 15849.7735593, 27638.5508593, 45055.6863125, 75413.4495285, 144231.595381, 201271.839517, 207597.542022, 202292.196166, 192313.587573, 174523.689315, 144497.159081, 113053.728647, 85815.9227619, 68406.266175}, 
{8.62274930115, 8.24927815446, 7.92075497275, 7.63150801547, 7.37586554207, 7.14815581201, 6.94270708473, 6.75384761969, 6.57590567635, 6.40320951416, 6.23008739258, 6.05086757105, 5.85987830903, 5.65144786597, 5.36009859326, 5.01799218606, 4.66181738325, 4.29813524851, 3.93350684555, 3.57449323805, 3.22765548971, 2.96807328859, 2.72432618252, 2.49351234367, 2.27272994425, 2.05907715644, 1.84965215244, 1.62273936823, 1.35945022879, 1.10186869834, 0.854731930536, 0.622777079014, 0.410741297426, 0.223361739418, 0.646835748046, 1.31806663905, 1.82994375403, 2.0682713703, 1.91885376518, 1.26749521597, 0, -1.81162671834, -4.49349716712, -8.18552268732, -13.0276146199, -19.1596843059, -26.7216430863, -38.1410697971, -55.5017938588, -73.4736523523, -91.2676901046, -108.094951942, -123.166482693, -135.693327183, -253.675451087, -407.50363087, -519.333801065, -566.125131811, -524.83679325, -372.427955523, -85.8577887716, 1663.82753952, 2138.85646874, 2238.70590781, 2862.85276565, 5951.07744117, 9836.78155925, 16667.7962299, 28591.9525633, 46520.1623884, 76311.0146582, 140728.893922, 192794.641992, 195621.318958, 186677.526934, 173944.796445, 156686.736867, 136607.022523, 114366.989097, 90983.3165417, 67472.6848087}, 
{6.00950006493, 5.59048440889, 5.23387262956, 4.93174208649, 4.67617013922, 4.45923414729, 4.27301147024, 4.10957946762, 3.96101549897, 3.81939692383, 3.67680110174, 3.52530539225, 3.35698715491, 3.16392374924, 2.84325618466, 2.44711890229, 2.03483260661, 1.617893103, 1.20779619683, 0.816037693466, 0.454113398287, 0.259251642306, 0.0998511819772, -0.0299567006062, -0.136040723349, -0.224269604156, -0.300512060932, -0.39886329836, -0.549177153723, -0.689829591235, -0.815228468148, -0.919781641714, -0.997896969189, -1.04398230782, -0.390790200988, 0.548705927741, 1.31144826309, 1.76769282152, 1.7876956195, 1.24171267351, 0, -2.02829191107, -5.01717310999, -9.10175917358, -14.4171656787, -21.098508202, -29.2809023205, -41.761055684, -60.9357398845, -80.5758998661, -99.7359592973, -117.470341846, -132.833471182, -144.879770973, -264.468407471, -419.925716147, -531.114906121, -574.221664962, -525.431680241, -360.930639529, -56.9042303956, 1594.94356216, 2089.56677847, 2282.86588267, 3030.74133893, 6119.83515652, 10110.2246116, 16974.6910876, 28686.0159675, 46142.5696035, 74540.3664726, 134749.094752, 183101.963485, 185223.623915, 175537.567216, 162020.602198, 145252.492382, 129664.17681, 111678.937722, 90449.9522673, 65130.3975988}, 
{5.11010590324, 4.76722634736, 4.44513271877, 4.14103371318, 3.85213802628, 3.57565435378, 3.30879139136, 3.04875783472, 2.79276237957, 2.5380137216, 2.28172055651, 2.02109157999, 1.75333548774, 1.47566097546, 1.14386920032, 0.781363418014, 0.414681252428, 0.0495009690587, -0.308499166598, -0.653640889047, -0.980245932794, -1.24175377417, -1.47901452583, -1.69199604224, -1.88066617787, -2.04499278719, -2.18494372468, -2.31351038997, -2.44172723919, -2.53842134288, -2.59827257771, -2.61596082031, -2.58616594733, -2.50356783541, -1.63095778755, -0.426011063169, 0.59819633833, 1.29728278117, 1.52686662958, 1.14256624778, 0, -2.20997862996, -5.44423123854, -9.8243843026, -15.472064299, -22.5088977046, -31.0565109963, -44.2119298624, -64.6250906557, -85.3031005118, -105.159478218, -123.107742559, -138.061412324, -148.934006298, -265.390018676, -416.281497008, -522.077128024, -559.03723378, -503.422136332, -331.492157736, -19.5076200467, 1364.58774607, 1849.72175887, 2181.06040736, 3103.76968052, 6329.96509383, 10219.6766594, 16678.4005155, 27611.6328, 44066.9767406, 70523.6072071, 124790.004577, 169886.804059, 177123.065366, 170747.696554, 157165.011731, 139091.885023, 122421.285321, 104029.24244, 83865.7561161, 61880.8260805}, 
{5.15399406081, 4.86732126432, 4.55749710807, 4.22709360111, 3.87868275245, 3.51483657113, 3.13812706617, 2.75112624661, 2.35640612147, 1.95653869979, 1.55409599059, 1.1516500029, 0.751772745752, 0.357036228175, -0.0150230558662, -0.371299685127, -0.719820983932, -1.06107580497, -1.39555300094, -1.72374142453, -2.04612992843, -2.41422047005, -2.77044353808, -3.10824272591, -3.42106162696, -3.70234383463, -3.94553294231, -4.14140282694, -4.2811285223, -4.36453721879, -4.38615650048, -4.34051395144, -4.22213715575, -4.02555369749, -2.94277721813, -1.47600855067, -0.188195297326, 0.761987363164, 1.21586425205, 1.01476019058, 0, -2.41110971808, -5.89466610136, -10.5507845072, -16.479580293, -23.781168816, -32.5556654336, -46.1493394636, -67.420686747, -88.7078977722, -108.793045078, -126.458201203, -140.485438686, -149.656830065, -260.331960641, -403.238975769, -500.96287712, -530.28164606, -467.973263956, -290.815712172, 24.4130279265, 1140.68681744, 1611.06884479, 2064.63976776, 3130.48024415, 6454.8293666, 10157.7834967, 16089.0036495, 26098.1508399, 41371.1182663, 65748.8703933, 114172.513928, 156044.737934, 168099.786951, 165181.991585, 152082.531554, 133172.11335, 115035.487808, 95821.7016099, 76555.6222204, 58262.117104}, 
{5.56760782684, 5.25695301437, 4.90469201207, 4.51547130935, 4.09393739563, 3.6447367603, 3.17251589278, 2.68192128247, 2.1775994188, 1.66419679115, 1.14635988895, 0.628735201597, 0.115969218505, -0.387291570918, -0.843684749148, -1.2692600423, -1.68254362809, -2.08558078983, -2.48041681083, -2.8690969744, -3.25366656385, -3.72203882403, -4.17853212288, -4.61333278987, -5.01662715449, -5.37860154622, -5.68944229452, -5.93249088567, -6.09211732554, -6.17487474549, -6.17372870127, -6.08164474867, -5.89158844345, -5.59652534138, -4.30944755491, -2.57999505157, -1.02441791103, 0.184327015565, 0.873282877043, 0.869492822243, 0, -2.64531908372, -6.40044687295, -11.3365324546, -17.5247249156, -25.0361733429, -33.9420268233, -47.7833967126, -69.5799910136, -91.1058827169, -111.023304726, -127.994489943, -140.681671273, -147.747081618, -251.275371018, -384.503373847, -472.839303231, -493.79119137, -424.867070465, -243.574972718, 72.5770696725, 1021.76580801, 1486.20386496, 2009.59026725, 3135.62404164, 6405.41436789, 9902.38827869, 15367.1274093, 24540.2133952, 38615.4318301, 60973.7524736, 104656.703708, 143126.24692, 156403.89239, 155188.277647, 143359.6563, 125478.951261, 107766.339134, 89124.3124556, 70950.0396685, 54640.6892174}, 
{6.02277142244, 5.64900739765, 5.23240635567, 4.77740005463, 4.28842025264, 3.76989870784, 3.22626717835, 2.66195742231, 2.08140119784, 1.48903026307, 0.889276376125, 0.286571295135, -0.314653221772, -0.909965416469, -1.46931037936, -2.00446565144, -2.52914833339, -3.04416767095, -3.55033290986, -4.04845329586, -4.53933807471, -5.1117399235, -5.66637906618, -6.19191915808, -6.6770238545, -7.11035681077, -7.48058168218, -7.7737772803, -7.97641082038, -8.08332661438, -8.08423784725, -7.96885770391, -7.72689936932, -7.34807602841, -5.85626405229, -3.85216851922, -2.00960409835, -0.515837912833, 0.44186291417, 0.6762312595, 0, -2.87543854517, -6.88390826439, -12.0605736037, -18.4405990091, -26.0591489267, -34.9513878025, -48.7693058796, -70.6214220455, -91.8946571806, -111.155628235, -126.970952161, -137.907245908, -142.531126426, -237.821836757, -360.089807422, -438.1960668, -450.577742495, -375.671962111, -191.915853252, 122.253456477, 1019.91083791, 1489.14034784, 2022.94308597, 3114.32015202, 6139.56704367, 9412.89766537, 14487.2579148, 22915.5936898, 35716.0685946, 56045.9754069, 96058.6866338, 131006.517651, 142116.301651, 141440.889095, 132417.846527, 117757.769152, 101000.950853, 83095.8621459, 65874.5535525, 51169.0755925}, 
{6.07046284803, 5.58562499995, 5.08344395691, 4.56462091454, 4.02985706848, 3.47985361435, 2.91531174779, 2.33693266445, 1.74541755995, 1.14146762992, 0.525784070003, -0.100931924167, -0.737979156956, -1.38465643273, -2.06000494023, -2.7508343794, -3.44246003872, -4.13014258099, -4.80914266899, -5.47472096555, -6.12213813345, -6.78464610876, -7.41426742588, -8.00101589246, -8.53490531617, -9.00594950466, -9.4041622656, -9.73401698649, -9.99781432998, -10.1509935968, -10.1776975405, -10.0620689148, -9.78825047319, -9.34038496929, -7.63854184167, -5.34000050643, -3.18258604081, -1.36808860109, -0.0982983435339, 0.424994575586, 0, -3.08931914507, -7.32137822272, -12.6888395146, -19.1843653022, -26.8006178672, -35.5302594912, -49.0355546094, -70.4273660222, -90.9239181495, -109.027896881, -123.241988106, -132.068877713, -134.011251594, -220.744835824, -331.830711028, -399.777169254, -404.029121038, -324.031476918, -139.229147432, 170.932956884, 1189.77782153, 1682.44204265, 2145.6887096, 3076.28091169, 5592.05626014, 8657.79799852, 13515.5526801, 21407.3668583, 32906.4806078, 51261.3599175, 89151.948951, 120385.79309, 124416.28671, 122483.423152, 118297.418691, 109902.395721, 95063.6148664, 78327.8777031, 61930.0867585, 48105.1445597}, 
{4.62204935772, 6.90952299937, 8.90489381107, 10.6166154287, 12.053141488, 13.222925625, 14.1344214755, 14.7960826753, 15.2163628603, 15.4037156663, 15.3665947292, 15.1134536849, 14.6527461692, 13.992925818, 13.0697255559, 11.9376024061, 10.6565180292, 9.24980048145, 7.74077781926, 6.15277809897, 4.50912937691, 2.7150431699, 0.928276777232, -0.811529041351, -2.46473352609, -3.99169591723, -5.35277545502, -6.38191504171, -6.93205314705, -7.2658249137, -7.39491867735, -7.33102277367, -7.08582553837, -6.67101530713, -5.46458777394, -3.87910755233, -2.37510699648, -1.09389103485, -0.176764595915, 0.234967391857, 0, -2.16592564661, -5.10331075854, -8.79561049299, -13.2262800072, -18.3787744583, -24.2365490037, -34.5712865407, -52.5854434187, -69.2043896062, -82.8734380684, -92.0379017703, -95.1430936768, -90.6343267531, -160.802107284, -251.05097294, -300.438749721, -290.260909227, -201.812923059, -16.3902628179, 284.711599896, 1250.30744487, 1730.83342155, 2172.12099829, 3020.00164346, 5114.1378625, 7928.91375283, 12382.7542465, 19394.0842758, 29180.1366911, 44762.9126691, 77227.9291106, 104204.637663, 108254.0882, 107079.823212, 103776.894566, 96782.3955451, 84281.453845, 70125.0180851, 56232.4985964, 44523.3057096}, 
{6.19769347791, 9.42888597242, 12.2483204165, 14.6692433757, 16.7049014157, 18.3685411019, 19.673409, 20.6327516756, 21.2598156941, 21.5678476211, 21.5700940223, 21.2798014632, 20.7102165094, 19.8745857264, 18.7018694304, 17.258633976, 15.6178258154, 13.8099315698, 11.8654378606, 9.81483130911, 7.68859853661, 5.36189965534, 3.04199945084, 0.780836199559, -1.36965182203, -3.35752633747, -5.13084907029, -6.45767396485, -7.13310328834, -7.53561792336, -7.6863685361, -7.60650579275, -7.31718035951, -6.83954290256, -5.70322397316, -4.24031165605, -2.82009994279, -1.56427587698, -0.594526502181, -0.0325388619938, 0, -1.74511611364, -4.10747507888, -7.05318392555, -10.5483496835, -14.5590793825, -19.0514800524, -27.6194636209, -43.3118209187, -57.4201631262, -68.4375919099, -74.8572089366, -75.172115873, -67.8754143858, -125.572841188, -199.873574535, -236.776792788, -219.616470609, -131.726582662, 43.5588963896, 322.905991882, 1199.51768019, 1646.56765537, 2054.61806207, 2814.23104499, 4571.69844227, 7086.78412465, 11056.923414, 17179.5516322, 25481.1173254, 38672.0161427, 66339.4749419, 89413.1721394, 93082.7299842, 92218.1236281, 89448.0205471, 83544.3903503, 73106.6845961, 61281.1249116, 49663.8909095, 39851.1622024}, 
{9.97694973869, 12.9226000518, 15.4319897049, 17.522364248, 19.2109692312, 20.5150502047, 21.4518527187, 22.0386223232, 22.2926045683, 22.2310450043, 21.8711891812, 21.2302826492, 20.3255709583, 19.1742996588, 17.7327145584, 16.0566498668, 14.2065572055, 12.2121591179, 10.1031781472, 7.9093368367, 5.6603577298, 3.22979769586, 0.82511250165, -1.50240775995, -3.70147299607, -5.72079311383, -7.50907802036, -8.84517000965, -9.53343601724, -9.92875881416, -10.0488200169, -9.91130124205, -9.53388410605, -8.93425022546, -7.57534171638, -5.82577042786, -4.09613954466, -2.51759806353, -1.22129498122, -0.338379294459, 0, -1.68086597422, -3.98301150913, -6.85203077637, -10.2335179476, -14.0730671945, -18.3162726886, -26.1586452108, -40.3073566825, -52.9370794831, -62.6738372091, -68.1436534573, -67.9725518243, -60.7865559067, -109.267448662, -171.519346538, -200.79762726, -182.62623823, -102.529126848, 53.9697594852, 301.346473367, 1062.39119843, 1459.06028797, 1823.82401596, 2489.15265637, 3973.06493914, 6150.34613441, 9576.12288204, 14805.5218219, 21799.9612154, 32895.6928378, 56258.4742576, 75718.2889753, 78727.852258, 77855.9156365, 75372.2877597, 70331.061764, 61723.557156, 52007.9034739, 42460.3580422, 34357.1781857}, 
{15.1393726701, 17.1695513703, 18.7741676078, 19.9758375375, 20.7971773143, 21.260803093, 21.3893310286, 21.205377276, 20.73155799, 19.9904893256, 19.0047874377, 17.7970684811, 16.3899486107, 14.8060439815, 13.0585488147, 11.1760396148, 9.18780608518, 7.11839155938, 4.99233937099, 2.83419285356, 0.668495340652, -1.60738836023, -3.82417216361, -5.93974851011, -7.91200984031, -9.69884859481, -11.2581572142, -12.4282794298, -13.0655225972, -13.4136337745, -13.4790460199, -13.2681923915, -12.7875059475, -12.0434197457, -10.3018322807, -8.03764147801, -5.78180046696, -3.69221253057, -1.92678095191, -0.643409014032, 0, -1.82715065645, -4.37906037469, -7.58262200962, -11.3647284161, -15.6522724491, -20.3721469634, -28.1673906718, -41.2729088703, -53.1209794508, -62.5356527702, -68.3409791854, -69.3610090533, -64.4197927308, -106.096340832, -159.209119671, -184.507581944, -169.820645137, -102.977226738, 28.1937557633, 235.86338488, 863.910670532, 1197.72686323, 1510.38297496, 2074.95001773, 3326.5642928, 5138.53680251, 7978.41535016, 12313.7477391, 18127.2070657, 27340.9652539, 46756.8148703, 62826.8806268, 65015.0952167, 63950.7904722, 61611.1873303, 57285.0914138, 50316.3215612, 42517.0590632, 34857.9943387, 28309.817807}, 
{20.8645168022, 21.9486260606, 22.5931200571, 22.8295227367, 22.6893580448, 22.2041499265, 21.4054223272, 20.3246991919, 18.9935044661, 17.4433620949, 15.7057960236, 13.8123301973, 11.7944885613, 9.68379506092, 7.57566007456, 5.46119275664, 3.32666634008, 1.19053732817, -0.928737775791, -3.01270246851, -5.04290024669, -7.07578416461, -9.00764263824, -10.8096736411, -12.4530751468, -13.9090451289, -15.1487815608, -16.0908784792, -16.6618342915, -16.9586289927, -16.9738194452, -16.6999625112, -16.1296150532, -15.2553339333, -13.103586943, -10.2780824168, -7.45565737458, -4.82647421416, -2.58069533314, -0.908483129204, 0, -2.0379455884, -4.94476200095, -8.63542858941, -13.0249247056, -18.0282297011, -23.5603229279, -31.624259365, -43.9093356423, -55.3377038032, -64.9765173971, -71.8929299737, -75.1540950824, -73.8271662728, -110.269928824, -156.16372466, -179.912985645, -171.730124376, -121.827553455, -20.4176854817, 142.287066945, 629.058767427, 891.982925002, 1144.93905411, 1601.80666918, 2640.52344295, 4070.29314937, 6301.86351789, 9745.98227794, 14453.393581, 21914.8558904, 37606.3845925, 50445.8395499, 51770.099056, 50460.3393703, 48226.2103849, 44549.1609273, 39069.2278481, 33020.2969708, 27092.8941434, 21977.5452138}, 
{26.3319366651, 27.0387102555, 27.2071129842, 26.8832793376, 26.1133438023, 24.9434408647, 23.4197050114, 21.5882707289, 19.4952725035, 17.186844822, 14.7091221707, 12.1082390361, 9.43032990481, 6.72152926327, 4.18033621289, 1.75649882871, -0.611768144157, -2.909495142, -5.12171260108, -7.23345095769, -9.22974064808, -11.101515369, -12.8270870155, -14.3906707432, -15.7764817077, -16.9687350645, -17.9516459691, -18.7168434114, -19.2548423637, -19.5321306573, -19.5299131929, -19.2293948712, -18.6117805926, -17.6582752581, -15.2014969806, -11.9492508544, -8.69628493242, -5.65873805034, -3.0527490438, -1.09445674847, 0, -2.16722619815, -5.3292567133, -9.40092147991, -14.2970504325, -19.9324735054, -26.2220206333, -34.5078106514, -45.9174951588, -56.9530933141, -66.949909894, -75.243249675, -81.168417434, -84.0607179475, -115.998623764, -155.603992228, -179.020167169, -178.885108996, -147.836778121, -78.5131349555, 36.4478600894, 382.818160061, 571.244017161, 758.136368412, 1099.90615084, 1923.26932926, 2964.55219542, 4584.53008476, 7143.97833277, 10769.0594661, 16524.3872467, 28579.0712367, 38282.0582009, 38818.5039711, 37342.1535664, 35278.8480497, 32265.9519319, 28166.5260531, 23729.3224881, 19401.1518005, 15628.8245537}, 
{30.7211867887, 32.2186900876, 32.9344123209, 32.9369668326, 32.2949669667, 31.0770260674, 29.3517574787, 27.1877745446, 24.6536906092, 21.8181190166, 18.7496731108, 15.5169662359, 12.1886117359, 8.83322295492, 5.76886510474, 2.9063473675, 0.13759651811, -2.51979741734, -5.04824441279, -7.43015444215, -9.64793747936, -11.6107076251, -13.3842933887, -14.9612274067, -16.3340423151, -17.4952707501, -18.4374453481, -19.1900504802, -19.7770180771, -20.1025249566, -20.1441001632, -19.8792727413, -19.2855717353, -18.3405261899, -15.8164536703, -12.4533044011, -9.08225780536, -5.92735897518, -3.21265300278, -1.16218498032, 0, -2.06896791379, -5.18168483713, -9.26957164526, -14.2640492134, -20.0965384168, -26.6984601306, -34.7966038923, -44.99824558, -55.3329887574, -65.4093090647, -74.8356821423, -83.2205836306, -90.1724891699, -117.492836778, -150.7507531, -173.835455322, -181.816032042, -169.761571857, -132.741163364, -65.82389516, 150.171519373, 264.925683569, 380.619032894, 599.432002815, 1183.12889144, 1840.25096105, 2864.47775029, 4549.4887978, 7064.74342554, 11076.5818222, 19446.7626156, 26042.4290357, 25985.9501575, 24553.8242957, 22830.5914509, 20578.1460551, 17792.4662126, 14855.8409064, 12018.8616543, 9532.11997425}, 
{33.2118217033, 37.2674516898, 40.0932839988, 41.7904447138, 42.4600599178, 42.2032556942, 41.1211581262, 39.3148932971, 36.8855872902, 33.9343661887, 30.5623560759, 26.870683035, 22.9604731494, 18.9328525022, 15.2375346252, 11.7551279095, 8.33985421252, 5.02153893592, 1.83000748119, -1.20491475009, -4.05340235637, -6.52948658449, -8.7810498512, -10.7998312216, -12.577569761, -14.1060045343, -15.3768746068, -16.3943759394, -17.1608326952, -17.638198079, -17.8131532561, -17.6723793916, -17.2025576507, -16.3903691987, -14.1693482892, -11.1924006672, -8.19215065829, -5.37069192474, -2.93011812895, -1.07252293325, 0, -1.5971461634, -4.15118669784, -7.6318500496, -12.0088646649, -17.25195899, -23.3308614712, -30.4691984487, -38.8524450661, -47.8432309068, -57.3081937134, -67.1139712284, -77.1272011946, -87.2145213545, -108.962978992, -134.824837999, -156.365178911, -171.053326562, -176.358605785, -169.750341412, -148.697858277, -43.8984836961, 2.44346809246, 43.0311625721, 130.567765226, 428.429069182, 716.326466699, 1179.76921399, 2004.26656725, 3330.98416411, 5478.46211654, 9981.34654168, 13433.8445105, 13098.0778108, 12052.9427936, 10942.9317147, 9628.42492444, 8131.29836296, 6611.55751688, 5182.1180492, 3955.89562293}, 
{32.9833959388, 41.9638811947, 49.0019939497, 54.2435724734, 57.8344550352, 59.9204799046, 60.6474853511, 60.1613096442, 58.6077910533, 56.132767848, 52.8820782976, 49.0015606718, 44.6370532399, 39.9343942714, 35.4826326492, 31.1472299911, 26.7600988247, 22.3764223516, 18.051383773, 13.8401662903, 9.79795310504, 6.216022101, 2.88085550384, -0.184969778351, -2.95887683747, -5.41828876543, -7.54062865413, -9.21369604263, -10.3387574813, -11.1075362129, -11.5338453717, -11.631498092, -11.4143075082, -10.8960867544, -9.48107211448, -7.56869726299, -5.6045381561, -3.72709183508, -2.0748553412, -0.786325715747, 0, -0.605736375047, -1.88690262082, -3.87822765706, -6.61444040352, -10.1302697799, -14.460444706, -19.5041536819, -25.1809517774, -31.8496605362, -39.600042644, -48.5218607862, -58.7048776485, -70.2388559164, -84.6194615314, -101.047077648, -118.615666741, -137.127425602, -156.384551026, -176.189239806, -196.343688735, -174.409178207, -186.787085404, -223.983127535, -276.503021814, -332.503197832, -388.284267226, -431.532824617, -449.935464628, -441.679613538, -362.949370992, -45.2891724322, 163.197081379, -19.4728737909, -202.899704633, -322.640032926, -440.529832534, -632.72745934, -791.822389068, -872.984670572, -831.384352706}
};

const double b4tab[21][81] = {
{-76.9329587423, -74.9120768374, -72.8831959254, -70.8465368937, -68.8023206297, -66.750768021, -64.692099955, -62.6265373192, -60.5543010011, -58.4756118881, -56.3906908678, -54.2997588276, -52.203036655, -50.1007452374, -47.9461157855, -45.7690950277, -43.6031866043, -41.4582227526, -39.3440357104, -37.2704577151, -35.2473210043, -33.4197450905, -31.6435908367, -29.9100063808, -28.2101398607, -26.5351394144, -24.8761531797, -23.2108203634, -21.5188100506, -19.8235717769, -18.1217387403, -16.4099441389, -14.6848211709, -12.9430030343, -11.4512563215, -10.0353268596, -8.50051256789, -6.78819309128, -4.83974807456, -2.59655716253, 0, 2.76386467791, 6.03212817555, 9.89720270697, 14.4515004862, 19.7874337273, 25.9974146443, 34.748917662, 47.4727452607, 60.4951560173, 73.2690362758, 85.2472723802, 95.8827506744, 104.628357502, 205.529346109, 338.199126315, 435.090882584, 476.309414359, 441.959521081, 312.146002193, 66.9736571378, -1739.15930399, -2089.86491288, -1995.19665002, -2465.20799591, -8373.90657653, -11200.2837093, -16591.1007556, -30193.1190768, -52541.1941081, -104617.80499, -252880.445399, -373123.584022, -369854.391849, -339140.269409, -301672.615926, -256155.658628, -206854.115601, -157777.913493, -113023.840048, -76688.6830109}, 
{-61.0665592504, -59.5405325996, -57.9862446317, -56.4062639038, -54.8031589728, -53.1794983957, -51.5378507295, -49.8807845312, -48.2108683576, -46.5306707659, -44.8427603129, -43.1497055557, -41.4540750511, -39.7584373562, -38.0534102497, -36.3491223946, -34.6566070748, -32.9808772805, -31.3269460014, -29.6998262275, -28.1045309487, -26.6183586579, -25.1640537188, -23.7366459979, -22.3311653619, -20.9426416776, -19.5661048115, -18.1863601135, -16.7897492896, -15.3957501937, -14.0035441748, -12.6123125819, -11.221236764, -9.82949807006, -8.8258656896, -7.96306628189, -6.96433642338, -5.74917049001, -4.2370628577, -2.34750790239, 0, 2.52321229601, 5.71548762464, 9.70743044737, 14.6296452257, 20.6127364211, 27.787308495, 38.051695097, 53.0386074135, 68.6518075209, 84.3041450268, 99.4084695388, 113.377630664, 125.624478011, 248.96719231, 411.077996819, 531.044830946, 585.084421041, 549.413493457, 400.248774545, 113.806990657, -2210.93177734, -2654.26793031, -2462.11743736, -2880.39626762, -8492.70369226, -12312.0248963, -19589.4637699, -35576.124203, -61261.3280215, -114681.525308, -254898.409019, -369567.039, -372005.488924, -350069.38676, -321533.331733, -283016.134207, -234619.988664, -183669.341978, -135026.498778, -93553.7636964}, 
{-46.6827849325, -44.8401848162, -43.1667879676, -41.6456599515, -40.2598663327, -38.9924726759, -37.8265445458, -36.7451475072, -35.7313471249, -34.7682089636, -33.8387985881, -32.9261815631, -32.0134234533, -31.0835898236, -30.0148520157, -28.8566320627, -27.6662920279, -26.4483527227, -25.207334958, -23.9477595451, -22.6741472952, -21.4036985249, -20.1265034133, -18.845331645, -17.5629529045, -16.2821368763, -15.0056532452, -13.7291299653, -12.4492681216, -11.1859140912, -9.94473733052, -8.73140729623, -7.5515934448, -6.41096523278, -5.88317934205, -5.61459443748, -5.20857505374, -4.55461355448, -3.54220230334, -2.06083366397, 0, 2.40324205281, 5.65496541032, 9.9136787163, 15.3378906145, 22.0861097486, 30.3168447625, 42.227978353, 59.7109528622, 78.0479006515, 96.5692789666, 114.605545053, 131.487156157, 146.544569523, 288.992582446, 475.996344633, 614.889333583, 678.435238199, 639.397747385, 470.540550044, 144.627335079, -2932.79776771, -3509.494448, -3113.86505125, -3274.31192294, -7801.70177205, -12594.3727884, -21919.6626421, -40044.9090035, -67533.3141439, -119764.62193, -248260.160271, -353625.192219, -358898.889463, -343194.682946, -321641.463116, -289396.777872, -243276.73674, -193434.545286, -145476.651948, -105009.505165}, 
{-31.7685999272, -30.9140803884, -30.069319311, -29.2330547729, -28.404024852, -27.5809676261, -26.762621173, -25.9477235708, -25.1350128971, -24.3232272298, -23.5111046469, -22.6973832261, -21.8808010453, -21.0600961825, -20.2160548816, -19.3587716161, -18.4996997342, -17.6412492127, -16.785830029, -15.9358521599, -15.0937255826, -14.2888854945, -13.4929842842, -12.7046995608, -11.9227089332, -11.1456900105, -10.3723204016, -9.59452646825, -8.80524902778, -8.01930869257, -7.23812429495, -6.46311466726, -5.69569864183, -4.937295051, -4.78900175182, -4.87287843703, -4.76559173962, -4.34590110341, -3.4925659722, -2.08434578979, 0, 2.42085545575, 5.82396507125, 10.3942168427, 16.3164987663, 23.7756988382, 32.9567050545, 46.4216152109, 66.3753228843, 87.3185359349, 108.470918792, 129.052135884, 148.281851642, 165.379730493, 326.163992173, 537.115520004, 693.618332648, 764.906571595, 720.214378334, 528.775894353, 159.825261141, -2523.47384838, -3117.91790692, -3040.05068028, -3706.41593424, -9248.08662758, -14385.7852229, -23793.4644584, -42145.077072, -72096.6026604, -124375.913386, -240089.396646, -335370.114979, -344937.678074, -330024.108378, -304523.678207, -269763.84717, -231537.355072, -190342.633636, -148483.000214, -108261.77216}, 
{-20.3254024719, -19.6641221415, -19.0350186708, -18.4356736921, -17.8636688373, -17.3165857385, -16.7920060277, -16.2875113371, -15.8006832986, -15.3291035444, -14.8703537065, -14.4220154169, -13.9816703078, -13.5469000111, -13.1457268425, -12.7564755465, -12.3551666421, -11.9331553712, -11.4817969754, -10.9924466966, -10.4564597767, -9.71331813501, -8.92722508221, -8.11051060658, -7.27550469635, -6.43453733979, -5.59993852513, -4.7964153973, -4.04681528532, -3.33387298221, -2.6648929529, -2.0471796623, -1.48803757532, -0.994771156883, -1.36259949059, -2.10038937803, -2.65736595606, -2.87967226598, -2.61345134908, -1.70484624665, 0, 2.24618514113, 5.6527777913, 10.4300873566, 16.7884232431, 24.938094857, 35.0894116043, 50.3112441763, 73.2429302252, 97.2596326421, 121.410993095, 144.746653253, 166.316254784, 185.169439357, 359.575466591, 587.861066754, 756.300746632, 831.377090892, 779.572684205, 567.370111238, 161.251956662, -2770.28073486, -3464.93510446, -3404.99012875, -4072.72478437, -8676.76636825, -14441.0790487, -24919.543214, -43666.0392523, -71423.6780248, -120178.648501, -231571.832103, -323807.043601, -333516.01987, -324255.681672, -307422.677717, -278071.652523, -229272.413294, -178376.120688, -134456.16511, -106585.936964}, 
{-12.6005017941, -11.9728803772, -11.4204448092, -10.9340172239, -10.5044197553, -10.1224745372, -9.77900370347, -9.46482938794, -9.1707737245, -8.887658847, -8.6063068893, -8.31753998527, -8.01218026877, -7.68104987367, -7.21755378068, -6.67414057027, -6.11063287214, -5.53777869374, -4.9663260425, -4.40702292587, -3.8706173513, -3.4781989565, -3.11493518347, -2.77633510447, -2.45790779178, -2.15516231767, -1.8636077544, -1.54917832458, -1.18225222586, -0.829055373876, -0.49710518022, -0.193919056472, 0.0729855857832, 0.296091334962, -0.412906419956, -1.50433589369, -2.36185359113, -2.81281918137, -2.68459233352, -1.8045327167, 0, 2.5936542577, 6.49959806732, 11.93300755, 19.1090588269, 28.2429280192, 39.5497912479, 56.6924011244, 82.8154706378, 109.890626767, 136.733218795, 162.158596005, 184.982107679, 204.019103101, 387.764131388, 627.69286075, 802.436991657, 876.105395361, 812.806943112, 576.650506163, 131.744955765, -2680.68058813, -3418.44304241, -3527.34547536, -4453.19095528, -9350.06518548, -15450.0183341, -26248.7762322, -45242.0647105, -73881.4550752, -121795.136732, -226173.231242, -310453.986792, -314547.267876, -299514.210114, -278387.625833, -250008.944779, -217192.420389, -181052.39391, -143282.72464, -105577.271875}, 
{-7.98239760804, -7.29399786559, -6.71050790948, -6.21896699896, -5.80641439332, -5.45988935183, -5.16643113375, -4.91307899835, -4.68687220491, -4.47485001269, -4.26405168096, -4.041516469, -3.79428363607, -3.50939244145, -3.01894399289, -2.40799198977, -1.7733182668, -1.13365340167, -0.507727972078, 0.0857274442849, 0.627982269727, 0.895130072613, 1.10195232043, 1.25805462674, 1.37304260507, 1.45652186899, 1.51809803203, 1.61275669825, 1.7886645851, 1.9469185659, 2.07869851983, 2.17518432605, 2.22755586375, 2.22699301211, 1.1650208586, -0.334795448054, -1.57322129634, -2.35460568096, -2.48329759664, -1.76364603808, 0, 2.90233959077, 7.25391782516, 13.264327862, 21.14316286, 31.1000159781, 43.344480375, 62.1158490198, 91.0479042651, 120.714766783, 149.689842121, 176.546535825, 199.858253442, 218.19840052, 404.917386665, 648.024046794, 822.297366765, 890.561560369, 815.640841397, 560.359423642, 87.5415208934, -2571.3685586, -3340.47974173, -3598.11833983, -4722.6106642, -9623.43596565, -15902.0864374, -26774.2572192, -45455.6434508, -73384.7446094, -119108.842822, -216742.708147, -295084.64224, -298106.040214, -281952.09004, -259635.596463, -232099.085313, -206514.179655, -177206.392525, -142880.611759, -102241.725194}, 
{-6.27186590723, -5.72375064278, -5.21238901175, -4.73323817132, -4.28175527869, -3.85339749102, -3.44362196552, -3.04788585935, -2.66164632971, -2.28036053378, -1.89948562874, -1.51447877177, -1.12079712006, -0.713897830798, -0.222281750843, 0.31623715635, 0.858776629682, 1.39618415205, 1.91930720634, 2.41899327545, 2.88608984228, 3.24266505273, 3.55784462439, 3.83197493784, 4.06540237367, 4.25847331247, 4.41153413484, 4.54610471123, 4.68052334188, 4.76450817485, 4.78980847345, 4.748173501, 4.63135252083, 4.43109479625, 3.0345115368, 1.13215007009, -0.501826043145, -1.64972295484, -2.09384681693, -1.61650378134, 0, 3.16932106442, 7.88612678066, 14.3410472067, 22.7247124005, 33.2277524202, 46.0407973236, 65.8806646373, 96.7840580597, 128.118976466, 158.23827236, 185.494798244, 208.24140662, 224.830949992, 407.156048185, 643.787663546, 810.160562735, 869.131473049, 783.557121786, 516.294236242, 30.1995437148, -2200.44302838, -2955.13410144, -3436.04892543, -4845.36275028, -9970.80989339, -16103.5300525, -26348.3541374, -43810.1130575, -70170.954719, -112803.759042, -200820.136095, -273937.973824, -285413.960155, -274710.647808, -252276.857897, -222592.029026, -195298.859233, -165356.391326, -132729.481512, -97382.9859978}, 
{-6.18844134475, -5.74473376836, -5.26598316987, -4.7564466173, -4.2203811787, -3.66204392209, -3.0856919155, -2.49558222697, -1.89597192452, -1.2911180762, -0.685277750021, -0.0827080140262, 0.512334063756, 1.09559141529, 1.63700573228, 2.14864190545, 2.64451786573, 3.12565397554, 3.59307059733, 4.04778809353, 4.49082682657, 5.00211767511, 5.49287240239, 5.9532132879, 6.3732626111, 6.74314265146, 7.05297568846, 7.2891782772, 7.43872380224, 7.50059535321, 7.46641985275, 7.32782422351, 7.07643538811, 6.7038802692, 4.99074622516, 2.69476822676, 0.665351603365, -0.858168588012, -1.63645729037, -1.4301794467, 0, 3.47197851865, 8.5666253211, 15.4433720311, 24.2616502724, 35.1808916689, 48.3605278441, 68.9225268074, 101.243176093, 133.615916047, 164.185229055, 191.095597505, 212.491503783, 226.517430276, 400.331557899, 625.159686523, 779.434023791, 826.72869685, 730.617832848, 454.675558934, -37.5239977463, -1839.01188812, -2571.49027083, -3251.95764485, -4897.41250916, -10186.8098627, -16037.666141, -25461.3116883, -41469.0768489, -65967.5576099, -105282.287387, -183837.39382, -251768.862259, -271191.524383, -266154.799741, -244481.867666, -213401.875653, -183777.627122, -152531.333642, -121337.702012, -91871.4390343}, 
{-6.72461798546, -6.25227089398, -5.71528399759, -5.1213013641, -4.47796706136, -3.79292515719, -3.07381971941, -2.32829481586, -1.56399451437, -0.78856288276, -0.00964398886505, 0.765118099488, 1.52807931447, 2.27159558825, 2.93249873787, 3.53826945337, 4.12059135565, 4.68317737575, 5.22974044473, 5.76399349363, 6.28964945349, 6.94417833777, 7.57906322557, 8.1795442784, 8.73086165775, 9.21825552513, 9.62696604203, 9.93243877378, 10.1115910409, 10.1730830032, 10.1061317429, 9.89995434217, 9.54376788347, 9.02678944903, 7.00970977686, 4.32217656321, 1.89440109448, -0.0126604627676, -1.13805194197, -1.22081717657, 0, 3.82889354967, 9.33916212403, 16.6476421701, 25.8711701347, 37.126582465, 50.5307156079, 71.5325146943, 104.779711721, 137.639462279, 168.063595352, 194.003939919, 213.412324962, 224.240579462, 387.405524105, 597.766739367, 737.848073951, 772.26821459, 665.64584802, 382.599660973, -112.251659816, -1646.0492215, -2370.91230837, -3167.71579244, -4917.33454572, -10127.3570029, -15666.1083256, -24366.5040447, -39061.4596911, -61671.3990537, -97766.7128899, -168676.893352, -231130.809087, -252580.391652, -250330.670612, -230726.243569, -201313.895601, -172378.917284, -142050.389723, -112601.574425, -86305.732896}, 
{-7.21285730645, -6.65573378737, -6.03213333636, -5.34938179523, -4.61480500578, -3.8357288098, -3.0194790491, -2.17338156546, -1.30476220069, -0.420946796584, 0.470738805055, 1.36296876243, 2.24841723375, 3.11975837721, 3.92444756023, 4.68376451309, 5.42241182928, 6.14231279695, 6.84539070424, 7.53356883931, 8.2087704903, 9.00877105673, 9.78087960798, 10.5082573248, 11.1740653879, 11.7614649781, 12.2536172762, 12.631298252, 12.8756422822, 12.9755138587, 12.9150426345, 12.6783582624, 12.2495903956, 11.6128686868, 9.29257638772, 6.19628405219, 3.34205860444, 1.01260875991, -0.509356765975, -0.941129257757, 0, 4.18541530255, 10.0897428028, 17.7762842333, 27.3083413269, 38.7492158163, 52.162209434, 73.1988669265, 106.671032657, 139.279881069, 168.819707661, 193.08480793, 209.869477375, 216.968011496, 367.740304703, 561.569286874, 686.062257543, 707.193854606, 590.938715958, 303.271479494, -189.83321689, -1640.86790786, -2375.69055353, -3194.40968758, -4897.13384372, -9724.7186013, -14923.4390021, -23022.2862036, -36550.2513631, -57144.7903552, -90009.5001861, -155045.718963, -211835.526591, -229737.746728, -228360.256732, -213338.598864, -189180.339259, -161776.855883, -132612.255506, -104673.130267, -80946.0723106}, 
{-6.8181922214, -6.11133366465, -5.37916772436, -4.62300955953, -3.84417432918, -3.04397719233, -2.22373330798, -1.38475783516, -0.528365932866, 0.344127239877, 1.23140652406, 2.13215676067, 3.04506279069, 3.96880945511, 4.92891571647, 5.90709105253, 6.88301371628, 7.84987985093, 8.80088559965, 9.72922710565, 10.6281005121, 11.5449591268, 12.4112486422, 13.2126719153, 13.9349318029, 14.5637311622, 15.0847728502, 15.5091815855, 15.8442621479, 16.0149308196, 15.9965683649, 15.764555548, 15.2942731332, 14.5611018847, 11.9296180707, 8.39420325835, 5.07116981743, 2.26537525951, 0.281677096219, -0.575067160826, 0, 4.52295874412, 10.7820589301, 18.7767730705, 28.5065736776, 39.9709332639, 53.1693243418, 73.8052083742, 106.72495423, 138.289145375, 166.18125152, 188.084742372, 201.683087639, 204.659757032, 342.470831378, 519.318287284, 628.218365849, 636.628870427, 512.007604372, 221.812371039, -266.499026219, -1911.56175663, -2686.45663171, -3397.68549434, -4851.75018742, -8873.56267204, -13757.5530173, -21532.3054313, -34226.4041221, -52759.5346508, -82486.9771663, -144214.753529, -195045.01482, -201335.675597, -197880.799806, -190780.070985, -176865.565514, -152524.569567, -125191.137308, -98533.3628351, -76219.340246}, 
{-4.86582622572, -7.87275909216, -10.5167004997, -12.8052894916, -14.7461651115, -16.3469664027, -17.6153324085, -18.5589021725, -19.185314738, -19.5022091484, -19.5172244471, -19.2379996775, -18.6721738831, -17.8273861073, -16.5859843899, -15.0348673205, -13.2704174654, -11.3259011072, -9.23458452878, -7.02973401276, -4.74461584187, -2.2424789485, 0.249912473321, 2.67581157988, 4.97847152747, 7.10114547237, 8.98708657089, 10.4029571991, 11.141954611, 11.5695806871, 11.7007901148, 11.5505375819, 11.1337777758, 10.4654653841, 8.58619594868, 6.12730730863, 3.79389199692, 1.80020181188, 0.360488551831, -0.310995984904, 0, 3.21081626952, 7.5974085403, 13.1388204936, 19.8140958108, 27.6022781731, 36.4824112619, 52.3619331941, 80.2574698828, 105.972596819, 127.079239043, 141.149321594, 145.754769511, 138.467507833, 250.45247075, 394.770158156, 474.368667824, 459.494578836, 320.394470275, 27.3149212221, -449.497489239, -2010.4964622, -2770.26591163, -3454.56281844, -4789.14416353, -8150.23269078, -12659.1073286, -19822.0777292, -31145.4535445, -46976.3436207, -72298.6600277, -125347.503982, -169366.086866, -175703.315771, -173464.87181, -167766.405333, -156083.446987, -135500.615836, -112306.296988, -89651.2143427, -70686.0918014}, 
{-7.26295953104, -11.5392956561, -15.2984976198, -18.5535934061, -21.3176109988, -23.603578382, -25.4245235395, -26.7934744552, -27.7234591131, -28.2275054971, -28.3186415912, -28.0098953793, -27.3142948452, -26.244867973, -24.6660202991, -22.6847989101, -20.4195009569, -17.9135539169, -15.2103852673, -12.3534224854, -9.38609304844, -6.12818177765, -2.87764538954, 0.291202055461, 3.30404649694, 6.08657387446, 8.56447012759, 10.4111715995, 11.3380181537, 11.8738533246, 12.0467844938, 11.8849190426, 11.4163643525, 10.6692278049, 8.90293646315, 6.63186768331, 4.42858453416, 2.47842160346, 0.966713479013, 0.0787947485979, 0, 2.61986287953, 6.18415734974, 10.6428564317, 15.9459331465, 22.0433605153, 28.885111559, 42.1261243707, 66.5640977749, 88.50778173, 105.590750116, 115.446576815, 115.708835706, 104.011100669, 196.481306804, 315.793062631, 375.650475665, 349.450077485, 210.588399667, -67.5380262126, -511.532668576, -1929.95144417, -2639.21283854, -3276.31667361, -4478.26277134, -7308.13585195, -11351.9290486, -17757.9399707, -27674.4662278, -41142.3569321, -62632.2551849, -107941.467664, -145662.228694, -151400.841708, -149673.084282, -144840.860869, -134925.328803, -117690.936301, -98269.7151758, -79282.0100537, -63348.1655615}, 
{-12.6652121681, -16.5941850036, -19.9590976418, -22.7796466972, -25.0755287846, -26.8664405185, -28.1720785134, -29.0121393841, -29.406319745, -29.3743162107, -28.9358253959, -28.1105439151, -26.918168383, -25.3783954141, -23.4000374728, -21.0729368917, -18.4953294453, -15.7095921845, -12.7581021606, -9.68323642457, -6.52737202763, -3.10862063081, 0.275402737707, 3.55134843982, 6.64586683742, 9.4856082924, 11.9972231667, 13.8715525391, 14.8348706551, 15.3776748044, 15.5223614954, 15.2913272362, 14.706968535, 13.7916819001, 11.6981648293, 8.99898859889, 6.3325568677, 3.89915588899, 1.89907191599, 0.532591201912, 0, 2.53934085691, 6.03171562615, 10.3959824547, 15.5509994894, 21.4156248773, 27.9087167653, 40.0803092376, 62.2094157938, 81.9456668721, 97.1244983057, 105.581345928, 105.15164557, 93.6708330674, 171.63073071, 271.925969606, 319.680544758, 291.732376541, 164.919385331, -83.9205084983, -477.949384571, -1710.16706257, -2340.85376322, -2912.90000436, -3969.1963038, -6363.71952508, -9872.01585756, -15410.2764475, -23894.6924412, -35262.4673388, -53366.7552247, -91674.0856748, -123519.279303, -128207.27605, -126494.261908, -122154.021871, -113665.050626, -99428.4992352, -83450.1835612, -67825.229948, -54648.7647392}, 
{-19.7282041677, -22.5206687815, -24.7330391222, -26.3951747594, -27.5369352623, -28.1881802004, -28.3787691432, -28.1385616599, -27.4974173202, -26.4851956933, -25.1317563486, -23.4669588557, -21.5206627839, -19.3227277026, -16.8768813728, -14.2295144412, -11.4289956175, -8.51052951663, -5.5093207534, -2.46057394272, 0.600506300526, 3.81973322185, 6.95588453916, 9.94875583077, 12.738142675, 15.2638406501, 17.4656453345, 19.1218277704, 20.0349299944, 20.5309711986, 20.6153311133, 20.2933894687, 19.5705259951, 18.4521204229, 15.7902062624, 12.3197742718, 8.86311843613, 5.66152592322, 2.95628390098, 0.988679537273, 0, 2.75849248442, 6.62949363721, 11.5053001324, 17.278208644, 23.8405158461, 31.0845184128, 43.2070151284, 63.758001827, 82.3372190107, 97.0891971773, 106.158466825, 107.689558451, 99.8270025535, 166.974133639, 252.707847982, 294.075630489, 271.578486173, 165.717420049, -43.0065628681, -374.092457566, -1391.38367739, -1922.74503651, -2414.26575523, -3312.03505387, -5333.43107969, -8255.36543617, -12849.4714516, -19887.3824542, -29341.5675946, -44381.1527342, -76222.7991113, -102523.077694, -105901.641436, -103917.229371, -99856.4726177, -92576.4521195, -81046.2729114, -68216.4938334, -55680.3540055, -45031.0925479}, 
{-27.1075555605, -28.8019886369, -29.8548606182, -30.3119029868, -30.2188472253, -29.621424816, -28.5653672412, -27.0964059836, -25.2602725254, -23.102698349, -20.669414937, -18.0061537716, -15.1586463354, -12.1726241107, -9.18539746089, -6.18476473404, -3.15159216058, -0.112879519672, 2.90437340951, 5.87316684779, 8.76650101599, 11.6604085101, 14.4106276989, 16.9759293261, 19.3150841357, 21.3868628715, 23.1500362774, 24.4997250457, 25.3406140509, 25.7836685794, 25.8135033411, 25.4147330462, 24.5719724046, 23.2698361265, 19.9973859775, 15.6853289185, 11.377578678, 7.36465296089, 3.9370694722, 1.38534591692, 0, 3.06656004483, 7.46690165057, 13.0779110348, 19.776474415, 27.4394780088, 35.9438080337, 48.4887693765, 67.7744337625, 85.7334049107, 100.893560298, 111.7827774, 116.928933693, 114.859906654, 173.584906762, 247.677666658, 286.452488242, 274.225416547, 195.312496604, 34.0297734445, -225.3067079, -1013.84164865, -1432.44300919, -1830.36687073, -2556.8693145, -4233.71788531, -6537.97546491, -10145.9092747, -15733.7865362, -23384.5504533, -35554.4403001, -61265.0490715, -82259.462866, -84262.9605047, -81930.8113565, -78098.7973863, -71933.3729466, -62877.2256033, -52937.4376817, -43246.862206, -34938.3522004}, 
{-33.4588863774, -34.9213862167, -35.5591006864, -35.4415567739, -34.6382814671, -33.2188017533, -31.2526446204, -28.8093370558, -25.9584060472, -22.7693785822, -19.3117816485, -15.6551422335, -11.868987325, -8.02284391048, -4.41443119891, -0.968920945903, 2.40578823857, 5.68684419988, 8.85139478342, 11.8765878346, 14.7395711987, 17.4169339639, 19.886459901, 22.1253740237, 24.1109013457, 25.8202668805, 27.2306956417, 28.3429721177, 29.1543407038, 29.5856930188, 29.6046881724, 29.1789852747, 28.2762434352, 26.864121764, 23.1380291901, 18.1867567556, 13.2332470319, 8.60765825673, 4.64014866789, 1.66087650307, 0, 3.25278582091, 8.03334993393, 14.2209167317, 21.6947106069, 30.3339559521, 40.0178771601, 52.9080993156, 70.8232894881, 88.1851913372, 103.946301233, 117.059115545, 126.476130642, 131.149842896, 182.536441249, 246.374394532, 284.427873405, 284.910177832, 236.034607778, 126.014463206, -56.936955917, -617.781336354, -917.504032084, -1211.15629541, -1753.78937864, -3081.02731146, -4755.84362429, -7369.97420881, -11515.1549567, -17396.3086685, -26765.6105093, -46478.2766534, -62314.2738182, -63070.2558967, -60523.8325501, -57031.5804547, -52009.6527709, -45254.3255842, -37981.8067955, -30924.2345294, -24813.7469101}, 
{-37.4378166491, -40.3621031681, -42.0802978837, -42.695861515, -42.312254781, -41.0329384008, -38.9613730932, -36.2010195774, -32.8553385724, -29.0277907971, -24.8218369706, -20.3409378119, -15.68855404, -10.968146374, -6.65282804863, -2.61221625241, 1.31205289293, 5.09212803554, 8.7001578236, 12.1082909053, 15.2886759287, 18.0928383129, 20.6302088296, 22.8895950214, 24.859804431, 26.529644601, 27.8879230739, 28.9892967387, 29.8785278322, 30.3869705889, 30.4766956009, 30.1097734602, 29.2482747588, 27.8542700886, 24.0304611152, 18.9151619993, 13.7874329363, 8.98966306549, 4.86424152626, 1.75355745798, 0, 3.1064120954, 7.81824875495, 14.0414187931, 21.6818310244, 30.6453942632, 40.8380173241, 53.4475322792, 69.4691468915, 85.7435450552, 101.656133549, 116.592319153, 129.937508645, 141.077108805, 184.90212827, 238.337000505, 275.618541362, 288.869780195, 270.213746352, 211.773469184, 105.671978043, -243.443100504, -425.484455994, -606.586973787, -952.885539247, -1891.80672766, -2944.96759483, -4592.05054564, -7312.737985, -11381.7349941, -17893.6559489, -31539.9229548, -42273.3495504, -42102.5502514, -39685.1176363, -36805.406101, -33079.1312559, -28510.5411273, -23718.3928643, -19111.9509556, -15100.4798899}, 
{-37.6999664064, -44.6073811379, -49.652990767, -52.9865426044, -54.7577839606, -55.1164621464, -54.2123244724, -52.1951182492, -49.2145907875, -45.4204893979, -40.962561391, -35.9905540776, -30.6542147681, -25.1032907733, -19.9894334718, -15.1448838292, -10.3638908845, -5.69354161915, -1.18092301468, 3.1268779475, 7.18277428592, 10.691650287, 13.8887021687, 16.761097417, 19.2960035177, 21.4805879568, 23.3020182201, 24.7764266611, 25.9155933152, 26.637427362, 26.9173356202, 26.730724909, 26.0530020471, 24.8595738535, 21.4930069681, 16.9616488663, 12.3974458299, 8.1097886419, 4.40806808553, 1.60167494387, 0, 2.41668115109, 6.31100838131, 11.6465187889, 18.386749472, 26.4952375289, 35.9355200578, 47.0895956007, 60.2765838605, 74.4594328297, 89.4317708137, 104.987226118, 120.919427047, 137.022001908, 171.755358998, 213.104453475, 247.6412475, 271.341233803, 280.179905109, 270.132754145, 237.175273639, 68.9326988884, -3.9406317289, -66.6118503895, -204.24808927, -682.503503456, -1141.34505706, -1882.52257707, -3207.78589052, -5345.72218375, -8817.56920571, -16127.4290734, -21722.5290621, -21138.8662085, -19403.4912999, -17570.8586029, -15415.6480652, -12978.840506, -10515.9875773, -8209.49146456, -6241.75435307}, {-32.90095568, -47.1404617732, -58.5117178931, -67.2253254363, -73.4918857993, -77.5220003786, -79.5262705707, -79.7152977721, -78.2996833793, -75.4900287888, -71.4969353972, -66.5310046009, -60.8028377966, -54.5230363805, -48.5130929303, -42.5971568519, -36.5531357806, -30.4666783707, -24.4234332762, -18.5090491511, -12.8091746497, -7.78310138406, -3.09123239759, 1.23238630822, 5.1537087319, 8.63868887197, 11.6532807269, 14.0420896375, 15.6679550323, 16.7869894101, 17.414418224, 17.5654669269, 17.2553609718, 16.4993258117, 14.343991964, 11.4173215729, 8.4205951511, 5.56715624071, 3.07034838394, 1.14351512297, 0, 0.972835270734, 3.00103908068, 6.14331828888, 10.4583797544, 16.0049303362, 22.8416768934, 30.8168166139, 39.8101782829, 50.3838214256, 62.6819265918, 76.8486743315, 93.0282451945, 111.364819731, 134.169524603, 160.215722341, 188.112747205, 217.561548824, 248.263076831, 279.918280855, 312.228110529, 279.105701817, 299.571089903, 358.816130254, 442.032678337, 530.434991644, 619.026308522, 688.225405054, 718.451057321, 706.837008684, 583.657133393, 81.7638927573, -247.651352966, 41.7735923421, 332.22177424, 521.477761603, 706.957137841, 1007.80800644, 1256.61737605, 1383.66396395, 1319.2264874}};
  
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
//  coeffs->b3  = 0.;
//  coeffs->b4  = 0.;
//  coeffs->b3  = -876.669217307; 
//  coeffs->b4  = 1386.13223658;

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
