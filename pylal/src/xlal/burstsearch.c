/*
 * Copyright (C) 2007  Kipp Cannon
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
 *               Python Wrapper For LAL's Burst Search Package
 *
 * ============================================================================
 */


#include <Python.h>
#include <math.h>


#include <lal/LALDatatypes.h>
#include <lal/TFTransform.h>
#include <lal/Date.h>


#define MODULE_NAME "pylal.xlal.burstsearch"


/*
 * ============================================================================
 *
 *                                 Functions
 *
 * ============================================================================
 */


static PyObject *pylal_XLALEPGetTimingParameters(PyObject *self, PyObject *args)
{
	int window_length;
	int max_tile_length;
	double fractional_tile_stride;
	int psd_length;
	int psd_shift;
	int window_shift;
	int window_pad;
	int tiling_length;

	psd_length = -1;
	if(!PyArg_ParseTuple(args, "iid|i", &window_length, &max_tile_length, &fractional_tile_stride, &psd_length))
		return NULL;

	if(XLALEPGetTimingParameters(window_length, max_tile_length, fractional_tile_stride, psd_length < 0 ? NULL : &psd_length, psd_length < 0 ? NULL : &psd_shift, &window_shift, &window_pad, &tiling_length) < 0) {
	}

	if(psd_length < 0)
		return Py_BuildValue("{s:i,s:i,s:i}", "window_shift", window_shift, "window_pad", window_pad, "tiling_length", tiling_length);
	return Py_BuildValue("{s:i,s:is:i,s:i,s:i}", "psd_length", psd_length, "psd_shift", psd_shift, "window_shift", window_shift, "window_pad", window_pad, "tiling_length", tiling_length);
}


/*
 * ============================================================================
 *
 *                            Module Registration
 *
 * ============================================================================
 */


static struct PyMethodDef methods[] = {
	{"XLALEPGetTimingParameters", pylal_XLALEPGetTimingParameters, METH_VARARGS, NULL},
	{NULL,}
};


void initburstsearch(void)
{
	Py_InitModule3(MODULE_NAME, methods, "Wrapper for LAL's burst search package.");
}
