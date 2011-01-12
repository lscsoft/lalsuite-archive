/*
*  Copyright (C) 2007 Xavier Siemens
*  Copyright (C) 2010 Andrew Mergl
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
/*            Cosmic string burst rate computation code for small loops          */
/*                                                                               */
/*                  Xavier Siemens, Jolien Creighton, Irit Maor                  */
/*                                                                               */
/*                         UWM/Caltech - September 2006                          */
/*********************************************************************************/
/*Modified June 2010 by Andrew Mergl for use with Python*/

#include <Python.h>
#include <numpy/arrayobject.h>
#include <math.h>
#include <stdlib.h>
#include "gsl/gsl_interp.h"
#include <gsl/gsl_errno.h>
#include <lal/cs_cosmo.h>
#include <lal/cs_lambda_cosmo.h>

#define CUSPS_PER_LOOP 1.0		/* c */
#define LOOP_RAD_POWER 50.0		/* Gamma */

double H0 = LAMBDA_H_0;

static PyObject *cs_gamma_findzofA(PyObject *self, PyObject *args);
static PyObject *cs_gamma_finddRdz(PyObject *self, PyObject *args);

/*****************************************************************************/
/*int finddRdz(double Gmu, double alpha, double f, double Gamma, int Namp, double *zofA, double *dRdz)
 * Find dR/dz given Gmu, alpha, f, and Gamma. This is a C function that was taken from cs_gamma.c
 * and modified so that it could be called from Python. The global variables
 * that the function used to use are now all passed through arguments.
 * 
 * Arguments:
 * Gmu, alpha, gamma: parameters calculated by the main program. See technical
 *  document for what they represent cosmologically.
 * f: the frequency that is passed to the main program with --frequency opt
 * Namp: The size of the data arrays. This is set when opening the data file
 * *zofA, *dRdz: 1D arrays of length Namp. See technical document for what 
 *  they represent cosmologically.
 */
/*****************************************************************************/
static PyObject *cs_gamma_finddRdz(PyObject *self, PyObject *args)
{
  //Declare incoming variables
  double Gmu=0.0, alpha=0.0, f=0.0, Gamma=0.0, *zofA=NULL, *dRdz=NULL;
  int Namp=0;
  PyObject *Numpy_zofA=NULL, *Numpy_dRdz=NULL; 
  //Extract the variables from *args
  if (!PyArg_ParseTuple(args, "ddddiOO", &Gmu, &alpha, &f,
			  &Gamma, &Namp, &Numpy_zofA, &Numpy_dRdz))
    return NULL;
  //Get the pointers to the actual and ensure the data is in a format that 
  // plays nicely with C: NPY_INOUT_ARRAY flag
  zofA = PyArray_DATA(
	       PyArray_FROM_OTF(Numpy_zofA, NPY_DOUBLE, NPY_INOUT_ARRAY));
  dRdz = PyArray_DATA(
	       PyArray_FROM_OTF(Numpy_dRdz, NPY_DOUBLE, NPY_INOUT_ARRAY));

  cs_cosmo_functions_t cosmofns;
  int j;

  cosmofns = XLALCSCosmoFunctions( zofA, (size_t) Namp);
  
  for ( j = 0; j < Namp; j++ )
    {


      double theta = pow((1+cosmofns.z[j]) * f * alpha * cosmofns.phit[j] / H0, -1.0/3.0);
      
      if (theta > 1.0)
	{
	  dRdz[j] = 0.0;
	}
      else
	{
	
	  dRdz[j] = 0.5 * H0 * pow(f/H0,-2.0/3.0) * pow(alpha, -5.0/3.0) / (Gamma*Gmu) *
	    pow(cosmofns.phit[j],-14.0/3.0) * cosmofns.phiV[j] * pow(1+cosmofns.z[j],-5.0/3.0);
	}

    }


  XLALCSCosmoFunctionsFree( cosmofns );

      return Py_BuildValue("i",0);
}
/*****************************************************************************/
/*int findzofA(double Gmu, double alpha, int Namp, double *zofA, double *amp)
 *Find z(A) given Gmu and alpha. This function was taken from cs_gamma.c and
 * modified so that it could be called from Python. The global variables it
 * used to use are now passed to it.
 *
 * Arugments:
 * Gmu, alpha: values calculated by the main program. See S4 technical
 *  documentation for their cosmological meanings.
 * Namp: The length of the data arrays. This is set when reading in the data
 * zofA, amp: 1D data arrays of length Namp. See S4 technical documentation
 *  for their cosmological meanings.
 */
/*****************************************************************************/
static PyObject *cs_gamma_findzofA(PyObject *self, PyObject *args)
{
  //Declare incoming variables
  double Gmu=0, alpha=0, *zofA=NULL, *amp=NULL;
  int Namp=0;
  PyObject *Numpy_zofA=NULL, *Numpy_amp=NULL;

  //Extract variables from *args
  if (!PyArg_ParseTuple(args, "ddiOO", &Gmu, &alpha, &Namp, 
			&Numpy_zofA, &Numpy_amp))
    return NULL;
 
  //Get the pointers to the actual and ensure the data is in a format that 
  // plays nicely with C: NPY_INOUT_ARRAY flag
  amp = PyArray_DATA(
	    PyArray_FROM_OTF(Numpy_amp, NPY_DOUBLE, NPY_IN_ARRAY));
  zofA  = PyArray_DATA(
            PyArray_FROM_OTF(Numpy_zofA,  NPY_DOUBLE, NPY_INOUT_ARRAY)); 
  double lnz_min = log(1e-20), lnz_max = log(1e10), dlnz =0.05;
  size_t numz       = floor( (lnz_max - lnz_min)/dlnz );
  int i,j;
  cs_cosmo_functions_t cosmofns;
  double *fz,*z;
  double a;
  gsl_interp *zofa_interp; 
  gsl_interp_accel *acc_zofa = gsl_interp_accel_alloc(); 

  cosmofns = XLALCSCosmoFunctionsAlloc( exp( lnz_min ), dlnz, numz );

  zofa_interp = gsl_interp_alloc (gsl_interp_linear, cosmofns.n);

  fz   = calloc( cosmofns.n, sizeof( *fz ) ); 
  z   = calloc( cosmofns.n, sizeof( *z ) ); 
  
  /* first compute the function that relates A and z */
  /* invert order; b/c fz is a monotonically decreasing func of z */
  j=0;
  for ( i = cosmofns.n-1 ; i >=  0; i-- )
    {
      z[j]=cosmofns.z[i];
      fz[j] = pow(cosmofns.phit[i],2.0/3.0) * pow(1+z[j],-1.0/3.0) / cosmofns.phiA[i];
      j=j+1;
    }

  gsl_interp_init (zofa_interp, fz, z, cosmofns.n);

  /* now compute the amplitudes (suitably multiplied) that are equal to fz for some z*/
  for ( j = 0; j < Namp; j++ )
    {
      a = amp[j] * pow(H0,-1.0/3.0) * pow(alpha,-2.0/3.0) / Gmu;
      zofA[j] = gsl_interp_eval (zofa_interp, fz, z, a, acc_zofa );
    }

  XLALCSCosmoFunctionsFree( cosmofns );
  free(fz);
  free(z);
  gsl_interp_free (zofa_interp);
  gsl_interp_accel_free(acc_zofa);

  return Py_BuildValue("i",0);
  
}

/*******************************************************************************/

//List of functions available to the Python module.
static PyMethodDef cs_gammaMethods[] = {
  {"findzofA", cs_gamma_findzofA, METH_VARARGS,
  "Function to find z(A). From cs_gamma.c; modified to be called from Python."},
  {"finddRdz", cs_gamma_finddRdz, METH_VARARGS,
  "Function to find dR/dz. From cs_gamma.c; modified to be called from Python."},
  {NULL, NULL, 0, NULL}
};

//They Python module initialization function.
PyMODINIT_FUNC
initcs_gamma(void)
{
  (void) Py_InitModule("cs_gamma", cs_gammaMethods);
  import_array();
}
