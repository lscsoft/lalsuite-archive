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
 * \brief In newer versions of the EOBNR approximant, we
 * do not have an analytic expression for the derivative of the waveform.
 * As such, it is necessary to calculate the derivatives numerically. This
 * function provides the means to do that for the SEOBNRv1 Hamiltonian
 * in the spin-aligned case, i.e. equatorial orbits with four dynamical variables:
 * r, phi, pr and pphi. It then combines energy flux and numerical derivatives of
 * the Hamiltonian to calculate the right hand side of the ODEs given by
 * Eqs. 10a - 10d of Pan et al. PRD 84, 124052 (2011).
 * Since SEOBNRv1 is a spin-aligned model, the ODEs are written
 * in the same format as ODEs of a nonspinning model.
 *
 */

#ifndef FINITE_DIFFERENCE_C
#define FINITE_DIFFERENCE_C

#include <unistd.h>

#include <gsl/gsl_deriv.h>

/*------------------------------------------------------------------------------------------
 *
 *          Prototypes of functions defined in this code.
 *
 *------------------------------------------------------------------------------------------
 */


 static int gsl_deriv_fixed( const gsl_function * f, 
                            double x, 
                            double h, 
                            double * result, 
                            double * abserr ); 
/*------------------------------------------------------------------------------------------
 *
 *          Defintions of functions.
 *
 *------------------------------------------------------------------------------------------
 */

 static int gsl_deriv_fixed( 
      const gsl_function * f, 
      double x, 
      double h, 
      double * result, 
      double * abserr )
{
   // Function to be called as -- Func(x) = (*(f->function))(double x, (void*) f->params );
   // *result = Func(x+h) - Func(x-h) / 2h 

   // Treat h as a fractional step size instead of absolute
   //h = h * x;
   *result = (*(f->function))( x + 0.5*h, (void*) f->params ) - (*(f->function))( x - 0.5*h, (void*) f->params );
   *result /= h;

   double ddF = (*(f->function))( x + 0.5*h, (void*) f->params ) + (*(f->function))( x - 0.5*h, (void*) f->params ) - 2*(*(f->function))( x , (void*) f->params );
   ddF /= (0.25*h*h);
   *abserr = abs( 0.25 * h * ddF );

   return GSL_SUCCESS;

}
#endif /* FINITE_DIFFERENCE_C */
