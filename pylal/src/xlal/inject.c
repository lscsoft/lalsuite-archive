/*
 * Copyright (C) 2006  Kipp Cannon
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
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */


/*
 * ============================================================================
 *
 *           Python Wrapper For LAL's Inject Package (and Friends)
 *
 * ============================================================================
 */


#include <Python.h>
#include <structmember.h>
#include <numpy/arrayobject.h>
#include <lal/DetResponse.h>


#define MODULE_NAME "pylal.xlal.inject"


/*
 * ============================================================================
 *
 *                                 Functions
 *
 * ============================================================================
 */


/*
 * Amplitude response.
 */


static PyObject *pylal_XLALComputeDetAMResponse(PyObject *self, PyObject *args)
{
	double fplus, fcross;
	double ra;
	double dec;
	double psi;
	double gmst;
	PyObject *response;

	/* 3x3 array, double, double, double, double */
	if(!PyArg_ParseTuple(args, "Odddd:XLALComputeDetAMResponse", &response, &ra, &dec, &psi, &gmst))
		return NULL;
	response = PyArray_FromAny(response, PyArray_DescrFromType(NPY_FLOAT32), 2, 2, NPY_CONTIGUOUS, NULL);
	if(!response || (PyArray_DIM(response, 0) != 3) || (PyArray_DIM(response, 1) != 3)) {
		Py_XDECREF(response);
		return NULL;
	}

	XLALComputeDetAMResponse(&fplus, &fcross, PyArray_DATA(response), ra, dec, psi, gmst);
	Py_DECREF(response);

	/* (double, double) */
	return Py_BuildValue("(dd)", fplus, fcross);
}


/*
 * ============================================================================
 *
 *                            Module Registration
 *
 * ============================================================================
 */


static struct PyMethodDef methods[] = {
	{"XLALComputeDetAMResponse", pylal_XLALComputeDetAMResponse, METH_VARARGS, "Compute the F_{+} and F_{\\times} amplitude response factors for a gravitational wave antenna.\n\nExample:\n\n>>> from xlal.tools import cached_detector\n>>> resp = cached_detector[\"LHO_4k\"].response\n>>> ra, dec = 0.0, 0.0\n>>> psi = 0.0\n>>> gmst = 0.0\n>>> fp, fc = XLALComputeDetAMResponse(resp, ra, dec, psi, gmst)"},
	{NULL,}
};


void initinject(void)
{
	Py_InitModule3(MODULE_NAME, methods, "Wrapper for LAL's inject package.");

	import_array();
}
