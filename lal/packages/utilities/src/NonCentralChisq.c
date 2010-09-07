/**
* @file
* @author Craig Robinson
* 
*  Copyright (C) 2010 Craig Robinson
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
*
* @section DESCRIPTION
*
* Functions for calculating the CDF and PDF of a non-central chi-squared
* distribution.
*/

#include <math.h>

#include <lal/NonCentralChisq.h>
#include <lal/LALConstants.h>
#include <lal/LALError.h>

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_errno.h>

#define LAL_NCCHISQ_MAXITER 1000


REAL8 XLALNonCentralChisqPDF( REAL8 x,
                              INT4  dof, 
                              REAL8 lambda )
{

  static const char func[] = "XLALNonCentralChisqPDF";

  INT4 i;
  INT4 initialIndex;
  INT4 count = 0;

  REAL8 besselArg, besselApprox;

  REAL8 logMultiplier;

  /* Stuff for calculating the bessel function */
  REAL8 besselNum, besselDenom, besselContrib;
  REAL8 besselOrder;

  REAL8 finalValue;

  if ( x < 0.0 || lambda < 0 || dof < 0 )
  {
    XLAL_ERROR_REAL8( func, XLAL_EINVAL );
  }

  logMultiplier  = - ( x + lambda ) / 2.0;
  logMultiplier += (0.25 * (REAL8)dof - 0.5) * log(x / lambda);
  logMultiplier -= LAL_LN2;

  besselArg     = sqrt( lambda * x);
  besselOrder   = (REAL8)dof / 2.0 - 1.0;
  besselApprox  = 0.0;
  besselDenom   = 0.0;
  besselContrib = 1.0;

  /* Find the starting index */
  /* We want to choose the point which contributed most to the sum */
  initialIndex = i = floor( besselArg / 2.0 );

  while (besselContrib > 1.0e-6 * besselApprox || i > initialIndex )
  {

    if ( i <= initialIndex )
    {
      i = initialIndex + count;
    }
    else
    {
      i = initialIndex - count;
    }

    if ( i >= 0 )
    {
      besselNum = log(besselArg) - LAL_LN2;
      besselNum *= (2.0 * i + besselOrder);
      
      besselDenom  = i ? gsl_sf_lnfact((REAL8)i) : 0.0;
      besselDenom += gsl_sf_lngamma((REAL8)i + besselOrder + 1.0);

      besselContrib = exp(besselNum - besselDenom + logMultiplier);

      besselApprox  += besselContrib;
    }

    if ( i <= initialIndex )
    {
      count++;
    }
  }
  finalValue = besselApprox;

  if ( isnan( finalValue ) )
  {
    XLAL_ERROR_REAL8( func, XLAL_EFPOVRFLW );
  }

  return finalValue;
}

REAL8 XLALNonCentralChisqCDF( REAL8 x,
                              INT4  dof,
                              REAL8 lambda )
{

  static const char func[] = "XLALNonCentralChisqCDF";

  REAL8 finalValue = 0.0;
  REAL8 lastValue  = 1.0;
  REAL8 multiplier;
  REAL8 q;

  INT4 i = 0;
  INT4 count = 0;
  INT4 initialIndex;

  INT4 maxIter;

  if ( x < 0.0 || lambda < 0 || dof < 0 )
  {
    XLAL_ERROR_REAL8( func, XLAL_EINVAL );
  }

  /* We start from the index with the greatest Poisson weighting */
  /* Integrate backwards and forwards from there */
  initialIndex = floor( lambda );

  /* Ensure that the number of iterations will allow us to include all */
  /* contributions which may prove to be significant. */
  maxIter = initialIndex > LAL_NCCHISQ_MAXITER ? initialIndex : LAL_NCCHISQ_MAXITER;

  i = initialIndex;

  while ( fabs( lastValue - finalValue ) > 1.0e-8 * finalValue || i > initialIndex 
            || ( finalValue == 0 && count < initialIndex) )
  {
    if ( i <= initialIndex )
    {
      i = initialIndex + count;
    }
    else
    {
      i = initialIndex - count;
    }

    if ( i >= 0 )
    {
      lastValue = finalValue;
      q = gsl_cdf_chisq_P( x, dof + 2.0*i );

      if ( q )
      {
        multiplier = - 0.5 * lambda;
        multiplier += (REAL8)i * (log(lambda) - LAL_LN2);
        multiplier -= i > 0 ? gsl_sf_lnfact((REAL8)i) : 0;

        finalValue += exp(multiplier + log(q));
      }
    }

    if ( i <= initialIndex )
    {
      count++;
    }
 
    if ( count > maxIter)
    {
      XLAL_ERROR_REAL8( func, XLAL_EMAXITER );
    }
  }

  return finalValue;
}
