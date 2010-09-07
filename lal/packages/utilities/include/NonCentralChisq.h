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
* Prototypes for functions for calculating the CDF and PDF of a 
* non-central chi-squared distribution.
*/

#ifndef _LALNCCHISQ_H
#define _LALNCCHISQ_H

#include <lal/LALDatatypes.h>
#include <lal/LALError.h>

#ifdef  __cplusplus
extern "C" {
#endif

#define LAL_NCCHISQ_MAXITER 1000

/**
 * Returns the value of the non-central chi-squared PDF
 * @f$ p(x; k, \lambda) @f$, where k is the number of degrees
 * of freedom, and @f$ \lambda @f$ is the non-centrality
 * parameter.
 *
 * @param x      The x value for which one required the PDF value
 * @param dof    The number of degrees of freedom
 * @param lambda The non-centrality parameter
 *
 * @return The value of the PDF at x.
 */
REAL8 XLALNonCentralChisqPDF( REAL8 x,
                              INT4  dof,
                              REAL8 lambda );


/**
 * Returns the value of the non-central chi-squared CDF
 * @f$ P(x; k, \lambda) @f$, where k is the number of degrees
 * of freedom, and @f$ \lambda @f$ is the non-centrality
 * parameter.
 *
 * @param x      The x value for which one required the CDF value
 * @param dof    The number of degrees of freedom
 * @param lambda The non-centrality parameter
 *
 * @return The value of the CDF at x.
 */
REAL8 XLALNonCentralChisqCDF( REAL8 x,
                              INT4  dof,
                              REAL8 lambda );

#ifdef  __cplusplus
}
#endif

#endif /* _LALNCCHISQ_H */
