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
 * Returns Py_False if the time-frequency tiles of the sngl_burst rows a
 * and b intersect after allowing for the light travel time.  Returns
 * Py_True if they do not intersect.  If the return values seem backwards,
 * it's the convention of comparison operators behaving like subtraction:
 * 0 == arguments are "equal".
 */


static PyObject *pylal_ExcessPowerCoincCompare(PyObject *self, PyObject *args)
{
	PyObject *a, *b;
	double a_central_freq, b_central_freq;
	double a_bandwidth, b_bandwidth;
	LIGOTimeGPS a_start, b_start;
	double a_duration, b_duration;
	double light_travel_time;
	double deltat;
	PyObject *attribute;

	if(!PyArg_ParseTuple(args, "OOd", &a, &b, &light_travel_time))
		return NULL;

	if((attribute = PyObject_GetAttrString(a, "central_freq"))) {
		a_central_freq = PyFloat_AsDouble(attribute);
		Py_DECREF(attribute);
	} else
		return NULL;
	if((attribute = PyObject_GetAttrString(a, "bandwidth"))) {
		a_bandwidth = PyFloat_AsDouble(attribute);
		Py_DECREF(attribute);
	} else
		return NULL;

	if((attribute = PyObject_GetAttrString(b, "central_freq"))) {
		b_central_freq = PyFloat_AsDouble(attribute);
		Py_DECREF(attribute);
	} else
		return NULL;
	if((attribute = PyObject_GetAttrString(b, "bandwidth"))) {
		b_bandwidth = PyFloat_AsDouble(attribute);
		Py_DECREF(attribute);
	} else
		return NULL;

	if(PyErr_Occurred())
		return NULL;

	if(fabs(a_central_freq - b_central_freq) > (a_bandwidth + b_bandwidth) / 2) {
		Py_INCREF(Py_True);
		return Py_True;
	}

	if((attribute = PyObject_GetAttrString(a, "start_time"))) {
		a_start.gpsSeconds = PyInt_AsLong(attribute);
		Py_DECREF(attribute);
	} else
		return NULL;
	if((attribute = PyObject_GetAttrString(a, "start_time_ns"))) {
		a_start.gpsNanoSeconds = PyInt_AsLong(attribute);
		Py_DECREF(attribute);
	} else
		return NULL;
	if((attribute = PyObject_GetAttrString(a, "duration"))) {
		a_duration = PyFloat_AsDouble(attribute);
		Py_DECREF(attribute);
	} else
		return NULL;

	if((attribute = PyObject_GetAttrString(b, "start_time"))) {
		b_start.gpsSeconds = PyInt_AsLong(attribute);
		Py_DECREF(attribute);
	} else
		return NULL;
	if((attribute = PyObject_GetAttrString(b, "start_time_ns"))) {
		b_start.gpsNanoSeconds = PyInt_AsLong(attribute);
		Py_DECREF(attribute);
	} else
		return NULL;
	if((attribute = PyObject_GetAttrString(b, "duration"))) {
		b_duration = PyFloat_AsDouble(attribute);
		Py_DECREF(attribute);
	} else
		return NULL;

	if(PyErr_Occurred())
		return NULL;

	deltat = XLALGPSDiff(&b_start, &a_start);
	if(deltat >= 0) {
		/* b starts at the same time as or after a */
		if(deltat > a_duration + light_travel_time) {
			Py_INCREF(Py_True);
			return Py_True;
		}
	} else {
		/* b starts before a */
		if(-deltat > b_duration + light_travel_time) {
			Py_INCREF(Py_True);
			return Py_True;
		}
	}

	/* time-frequency tiles intersect */

	Py_INCREF(Py_False);
	return Py_False;
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
	{"ExcessPowerCoincCompare", pylal_ExcessPowerCoincCompare, METH_VARARGS, NULL},
	{NULL,}
};


void initburstsearch(void)
{
	Py_InitModule3(MODULE_NAME, methods, "Wrapper for LAL's burst search package.");
}
