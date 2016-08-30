/*
*  Copyright (C) 2010 Michele Vallisneri, Will Farr, Evan Ochsner, 2014 A. Klein
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
#ifndef _LALSIMIMRSPINALIGNEDEOBGSLOPTIMIZEDINTERPOLATION_C
#define _LALSIMIMRSPINALIGNEDEOBGSLOPTIMIZEDINTERPOLATION_C

#include "LALSimIMRSpinAlignedEOBOptimizedInterpolatorGeneral.c"

SEOBNRv2OptimizedInterpolatorNoAmpPhase (REAL8Array * yin, REAL8 tinit,
					 REAL8 deltat, UINT4 num_input_times,
					 REAL8Array ** yout)
{
  int errnum = 0;

  /* needed for the final interpolation */
  gsl_spline *interp = NULL;
  gsl_interp_accel *accel = NULL;
  int outputlen = 0;
  REAL8Array *output = NULL;
  REAL8 *times, *vector;	/* aliases */

  /* needed for the integration */
  size_t dim = 4;
  
  interp = gsl_spline_alloc (gsl_interp_cspline, num_input_times);
  accel = gsl_interp_accel_alloc ();

  outputlen = (int) (yin->data[num_input_times - 1] / deltat) + 1;
  output = XLALCreateREAL8ArrayL (2, dim + 1, outputlen);	/* Only dim + 1 rather than dim + 3 since we're not adding amp & phase */

  if (!interp || !accel || !output)
    {
      errnum = XLAL_ENOMEM;	/* ouch again, ran out of memory */
      if (output)
	XLALDestroyREAL8Array (output);
      outputlen = 0;
      goto bail_out;
    }

  /* make an array of times */
  times = output->data;
  for (int j = 0; j < outputlen; j++)
    times[j] = tinit + deltat * j;

  /* interpolate! */
  for (unsigned int i = 1; i <= dim; i++)
    {				/* only up to dim (4) because we are not interpolating amplitude and phase */
      //gsl_spline_init(interp, &yin->data[0], &yin->data[num_input_times * i], num_input_times + 1);
      gsl_spline_init (interp, &yin->data[0], &yin->data[num_input_times * i],
		       num_input_times);

      vector = output->data + outputlen * i;
      unsigned int index_old = 0;
      double x_lo_old = 0, y_lo_old = 0, b_i_old = 0, c_i_old = 0, d_i_old =
	0;
      for (int j = 0; j < outputlen; j++)
	{
	  optimized_gsl_spline_eval_e (interp, times[j], accel, &(vector[j]),
				       &index_old, &x_lo_old, &y_lo_old,
				       &b_i_old, &c_i_old, &d_i_old);
	}
    }

  /* deallocate stuff and return */
bail_out:

  if (interp)
    XLAL_CALLGSL (gsl_spline_free (interp));
  if (accel)
    XLAL_CALLGSL (gsl_interp_accel_free (accel));
  if (errnum)
    XLAL_ERROR (errnum);

  *yout = output;
  return outputlen;
}

UNUSED static int
SEOBNRv2OptimizedInterpolatorOnlyAmpPhase (REAL8Array * yin, REAL8 tinit,
					   REAL8 deltat,
					   UINT4 num_input_times,
					   REAL8Array ** yout)
{
  int errnum = 0;

  /* needed for the integration */
  size_t dim = 4;

  /* needed for the final interpolation */
  gsl_spline *interp = NULL;
  gsl_interp_accel *accel = NULL;
  int outputlen = 0;
  REAL8Array *output = NULL;
  REAL8 *times, *vector;	/* aliases */

  /* note: for speed, this replaces the single CALLGSL wrapper applied before each GSL call */
  interp = gsl_spline_alloc (gsl_interp_cspline, num_input_times);
  accel = gsl_interp_accel_alloc ();


  outputlen = (int) (yin->data[num_input_times - 1] / deltat) + 1;
  output = XLALCreateREAL8ArrayL (2, dim + 3, outputlen);	/* Original (dim+1), Optimized (dim+3), since we're adding amp & phase */

  if (!interp || !accel || !output)
    {
      errnum = XLAL_ENOMEM;	/* ouch again, ran out of memory */
      if (output)
	XLALDestroyREAL8Array (output);
      outputlen = 0;
      goto bail_out;
    }

  /* make an array of times */
  times = output->data;
  for (int j = 0; j < outputlen; j++)
    times[j] = tinit + deltat * j;

  /* interpolate! */
  for (unsigned int i = dim + 1; i <= dim + 2; i++)
    {				/* Original (dim), Optimized (dim+2), since we're also interpolating amp & phase */
      //gsl_spline_init(interp, &yin->data[0], &yin->data[num_input_times * i], num_input_times + 1);
      gsl_spline_init (interp, &yin->data[0], &yin->data[num_input_times * i],
		       num_input_times);
      vector = output->data + outputlen * i;
      unsigned int index_old = 0;
      double x_lo_old = 0, y_lo_old = 0, b_i_old = 0, c_i_old = 0, d_i_old =
	0;
      for (int j = 0; j < outputlen; j++)
	{
	  optimized_gsl_spline_eval_e (interp, times[j], accel, &(vector[j]),
				       &index_old, &x_lo_old, &y_lo_old,
				       &b_i_old, &c_i_old, &d_i_old);
	}
    }
  /* deallocate stuff and return */
bail_out:

  if (interp)
    XLAL_CALLGSL (gsl_spline_free (interp));
  if (accel)
    XLAL_CALLGSL (gsl_interp_accel_free (accel));

  if (errnum)
    XLAL_ERROR (errnum);

  *yout = output;
  return outputlen;
}

UNUSED static int
SEOBNRv2OptimizedInterpolatorIncludeAmpPhase (REAL8Array * yin, REAL8 tinit, REAL8 deltat, UINT4 num_input_times, REAL8Array ** yout)	// Davids: maybe this function can be removed since it is split into noampphase and onlyampphase
{
  int errnum = 0;

  /* needed for the integration */
  size_t dim = 4;

  /* needed for the final interpolation */
  gsl_spline *interp = NULL;
  gsl_interp_accel *accel = NULL;
  int outputlen = 0;
  REAL8Array *output = NULL;
  REAL8 *times, *vector;	/* aliases */

  /* note: for speed, this replaces the single CALLGSL wrapper applied before each GSL call */
  interp = gsl_spline_alloc (gsl_interp_cspline, num_input_times);
  accel = gsl_interp_accel_alloc ();


  outputlen = (int) (yin->data[num_input_times - 1] / deltat) + 1;
  output = XLALCreateREAL8ArrayL (2, dim + 3, outputlen);	/* Original (dim+1), Optimized (dim+3), since we're adding amp & phase */

  if (!interp || !accel || !output)
    {
      errnum = XLAL_ENOMEM;	/* ouch again, ran out of memory */
      if (output)
	XLALDestroyREAL8Array (output);
      outputlen = 0;
      goto bail_out;
    }

  /* make an array of times */
  times = output->data;
  for (int j = 0; j < outputlen; j++)
    times[j] = tinit + deltat * j;
  /* interpolate! */
  for (unsigned int i = 1; i <= dim + 2; i++)
    {				/* Original (dim), Optimized (dim+2), since we're also interpolating amp & phase */
      //gsl_spline_init(interp, &yin->data[0], &yin->data[num_input_times * i], num_input_times + 1);
      gsl_spline_init (interp, &yin->data[0], &yin->data[num_input_times * i],
		       num_input_times);
      vector = output->data + outputlen * i;
      unsigned int index_old = 0;
      double x_lo_old = 0, y_lo_old = 0, b_i_old = 0, c_i_old = 0, d_i_old =
	0;
      for (int j = 0; j < outputlen; j++)
	{
	  optimized_gsl_spline_eval_e (interp, times[j], accel, &(vector[j]),
				       &index_old, &x_lo_old, &y_lo_old,
				       &b_i_old, &c_i_old, &d_i_old);
	}
    }
  /* deallocate stuff and return */
bail_out:

  if (interp)
    XLAL_CALLGSL (gsl_spline_free (interp));
  if (accel)
    XLAL_CALLGSL (gsl_interp_accel_free (accel));

  if (errnum)
    XLAL_ERROR (errnum);

  *yout = output;
  return outputlen;
}


#endif /* _LALSIMIMRSPINALIGNEDEOBGSLOPTIMIZEDINTERPOLATION_C */
