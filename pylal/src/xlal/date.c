/*
 * Copyright (C) 2006  Kipp C. Cannon
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
 *                   Python Wrapper For LAL's Date Package
 *
 * ============================================================================
 */


#include <Python.h>
#include <structmember.h>
#include <stdlib.h>
#include <time.h>
#include <numpy/arrayobject.h>
#include <lal/Date.h>
#include <lal/TimeDelay.h>
#include <datatypes/ligotimegps.h>


#define MODULE_NAME "pylal.xlal.date"


/*
 * ============================================================================
 *
 *                       LIGOTimeGPS Function Wrappers
 *
 * ============================================================================
 */


static PyObject *pylal_XLALGPSToINT8NS(PyObject *self, PyObject *args)
{
	pylal_LIGOTimeGPS *s;

	/* LIGOTimeGPS */
	if(!PyArg_ParseTuple(args, "O!:XLALGPSToINT8NS", pylal_LIGOTimeGPS_Type, &s))
		return NULL;

	/* long */
	return PyLong_FromLongLong(XLALGPSToINT8NS(&s->gps));
}


static PyObject *pylal_XLALINT8NSToGPS(PyObject *self, PyObject *args)
{
	long long ns;
	pylal_LIGOTimeGPS *new;

	/* long */
	if(!PyArg_ParseTuple(args, "L:XLALINT8NSToGPS", &ns))
		return NULL;

	new = (pylal_LIGOTimeGPS *) _PyObject_New(pylal_LIGOTimeGPS_Type);
	XLALINT8NSToGPS(&new->gps, ns);

	/* LIGOTimeGPS */
	return (PyObject *) new;
}


/*
 * ============================================================================
 *
 *                              Time Conversion
 *
 * ============================================================================
 */


/*
 * Convert a struct tm of the kind Python uses to/from a struct tm of the
 * kind the C library uses.
 *
 * Python:
 *	1900 = year 1900
 *	January = month 1
 *	Monday = week day 0
 *	January 1st = year day 1
 *
 * C:
 *	1900 = year 0
 *	January = month 0
 *	Monday = week day 1
 *	January 1st = year day 0
 */


static struct tm *struct_tm_python_to_c(struct tm *tm)
{
	tm->tm_year -= 1900;
	tm->tm_mon -= 1;
	tm->tm_wday = (tm->tm_wday + 8) % 7;
	tm->tm_yday -= 1;

	return tm;
}


static struct tm *struct_tm_c_to_python(struct tm *tm)
{
	tm->tm_year += 1900;
	tm->tm_mon += 1;
	tm->tm_wday = (tm->tm_wday + 6) % 7;
	tm->tm_yday += 1;

	return tm;
}


/*
 * Leap seconds.
 */


static PyObject *pylal_XLALLeapSeconds(PyObject *self, PyObject *args)
{
	int gpssec;

	/* int */
	if(!PyArg_ParseTuple(args, "i:XLALLeapSeconds", &gpssec))
		return NULL;

	/* int */
	return PyInt_FromLong(XLALLeapSeconds(gpssec));
}


static PyObject *pylal_XLALLeapSecondsUTC(PyObject *self, PyObject *args)
{
	struct tm utc;

	/* time.struct_time */
	if(!PyArg_ParseTuple(args, "(iiiiiiiii):XLALLeapSecondsUTC", &utc.tm_year, &utc.tm_mon, &utc.tm_mday, &utc.tm_hour, &utc.tm_min, &utc.tm_sec, &utc.tm_wday, &utc.tm_yday, &utc.tm_isdst))
		return NULL;

	/* int */
	return PyInt_FromLong(XLALLeapSecondsUTC(struct_tm_python_to_c(&utc)));
}


/*
 * GPS to/from UTC
 */


static PyObject *pylal_XLALGPSToUTC(PyObject *self, PyObject *args)
{
	pylal_LIGOTimeGPS *gps;
	struct tm utc;

	/* LIGOTimeGPS */
	if(!PyArg_ParseTuple(args, "O!:XLALGPSToUTC", pylal_LIGOTimeGPS_Type, &gps))
		return NULL;
	if(gps->gps.gpsNanoSeconds) {
		PyErr_SetString(PyExc_TypeError, "cannot convert non-integer seconds");
		return NULL;
	}

	XLALGPSToUTC(&utc, gps->gps.gpsSeconds);
	struct_tm_c_to_python(&utc);

	/* time.struct_time */
	return Py_BuildValue("(iiiiiiiii)", utc.tm_year, utc.tm_mon, utc.tm_mday, utc.tm_hour, utc.tm_min, utc.tm_sec, utc.tm_wday, utc.tm_yday, utc.tm_isdst);
}


static PyObject *pylal_XLALUTCToGPS(PyObject *self, PyObject *args)
{
	LIGOTimeGPS gps;
	struct tm utc;

	/* time.struct_time */
	if(!PyArg_ParseTuple(args, "(iiiiiiiii):XLALUTCToGPS", &utc.tm_year, &utc.tm_mon, &utc.tm_mday, &utc.tm_hour, &utc.tm_min, &utc.tm_sec, &utc.tm_wday, &utc.tm_yday, &utc.tm_isdst))
		return NULL;

	XLALGPSSet(&gps, XLALUTCToGPS(struct_tm_python_to_c(&utc)), 0);

	/* LIGOTimeGPS */
	return pylal_LIGOTimeGPS_new(gps);
}


/*
 * Julian day
 */


static PyObject *pylal_XLALJulianDay(PyObject *self, PyObject *args)
{
	struct tm utc;

	/* time.struct_time */
	if(!PyArg_ParseTuple(args, "(iiiiiiiii):XLALJulianDay", &utc.tm_year, &utc.tm_mon, &utc.tm_mday, &utc.tm_hour, &utc.tm_min, &utc.tm_sec, &utc.tm_wday, &utc.tm_yday, &utc.tm_isdst))
		return NULL;

	/* float */
	return PyFloat_FromDouble(XLALJulianDay(struct_tm_python_to_c(&utc)));
}


static PyObject *pylal_XLALModifiedJulianDay(PyObject *self, PyObject *args)
{
	struct tm utc;

	/* time.struct_time */
	if(!PyArg_ParseTuple(args, "(iiiiiiiii):XLALModifiedJulianDay", &utc.tm_year, &utc.tm_mon, &utc.tm_mday, &utc.tm_hour, &utc.tm_min, &utc.tm_sec, &utc.tm_wday, &utc.tm_yday, &utc.tm_isdst))
		return NULL;

	/* int */
	return PyInt_FromLong(XLALModifiedJulianDay(struct_tm_python_to_c(&utc)));
}


/*
 * Sidereal time
 */


static PyObject *pylal_XLALGreenwichSiderealTime(PyObject *self, PyObject *args)
{
	pylal_LIGOTimeGPS *gps;
	double equation_of_equinoxes;

	/* LIGOTimeGPS, float */
	if(!PyArg_ParseTuple(args, "O!d:XLALGreenwichSiderealTime", pylal_LIGOTimeGPS_Type, &gps, &equation_of_equinoxes))
		return NULL;

	/* float */
	return PyFloat_FromDouble(XLALGreenwichSiderealTime(&gps->gps, equation_of_equinoxes));
}


static PyObject *pylal_XLALGreenwichMeanSiderealTimeToGPS(PyObject *self, PyObject *args)
{
	LIGOTimeGPS gps;
	double gmst;

	/* float */
	if(!PyArg_ParseTuple(args, "d:XLALGreenwichMeanSiderealTimeToGPS", &gmst))
		return NULL;

	XLALGreenwichMeanSiderealTimeToGPS(gmst, &gps);

	/* LIGOTimeGPS */
	return pylal_LIGOTimeGPS_new(gps);
}


/*
 * ============================================================================
 *
 *                             Propogation Delay
 *
 * ============================================================================
 */


static PyObject *pylal_XLALArrivalTimeDiff(PyObject *self, PyObject *args)
{
	PyObject *pos1, *pos2;
	double ra, dec;
	pylal_LIGOTimeGPS *gps;
	PyObject *result;

	/* 3-element list, 3-element list, float, float, LIGOTimeGPS */
	if(!PyArg_ParseTuple(args, "OOddO!:XLALArrivalTimeDiff", &pos1, &pos2, &ra, &dec, pylal_LIGOTimeGPS_Type, &gps))
		return NULL;
	pos1 = PyArray_FromAny(pos1, PyArray_DescrFromType(NPY_FLOAT64), 1, 1, NPY_CONTIGUOUS, NULL);
	if(!pos1) {
		PyErr_SetString(PyExc_TypeError, "XLALArrivalTimeDiff() unable to convert argument 1 to array");
		return NULL;
	}
	pos2 = PyArray_FromAny(pos2, PyArray_DescrFromType(NPY_FLOAT64), 1, 1, NPY_CONTIGUOUS, NULL);
	if(!pos2) {
		PyErr_SetString(PyExc_TypeError, "XLALArrivalTimeDiff() unable to convert argument 2 to array");
		Py_DECREF(pos1);
		return NULL;
	}
	if(PyArray_DIM(pos1, 0) != 3) {
		PyErr_SetString(PyExc_TypeError, "XLALArrivalTimeDiff() argument 1 must have length 3");
		Py_DECREF(pos1);
		Py_DECREF(pos2);
	}
	if(PyArray_DIM(pos2, 0) != 3) {
		PyErr_SetString(PyExc_TypeError, "XLALArrivalTimeDiff() argument 2 must have length 3");
		Py_DECREF(pos1);
		Py_DECREF(pos2);
	}

	result = PyFloat_FromDouble(XLALArrivalTimeDiff(PyArray_DATA(pos1), PyArray_DATA(pos2), ra, dec, &gps->gps));
	Py_DECREF(pos1);
	Py_DECREF(pos2);

	/* float */
	return result;
}


/*
 * ============================================================================
 *
 *                            Module Registration
 *
 * ============================================================================
 */


static struct PyMethodDef module_methods[] = {
	{"XLALArrivalTimeDiff", pylal_XLALArrivalTimeDiff, METH_VARARGS, NULL},
	{"XLALGPSToINT8NS", pylal_XLALGPSToINT8NS, METH_VARARGS, NULL},
	{"XLALGPSToUTC", pylal_XLALGPSToUTC, METH_VARARGS, NULL},
	{"XLALGreenwichSiderealTime", pylal_XLALGreenwichSiderealTime, METH_VARARGS, NULL},
	{"XLALGreenwichMeanSiderealTimeToGPS", pylal_XLALGreenwichMeanSiderealTimeToGPS, METH_VARARGS, NULL},
	{"XLALINT8NSToGPS", pylal_XLALINT8NSToGPS, METH_VARARGS, NULL},
	{"XLALJulianDay", pylal_XLALJulianDay, METH_VARARGS, NULL},
	{"XLALLeapSeconds", pylal_XLALLeapSeconds, METH_VARARGS, NULL},
	{"XLALLeapSecondsUTC", pylal_XLALLeapSecondsUTC, METH_VARARGS, NULL},
	{"XLALModifiedJulianDay", pylal_XLALModifiedJulianDay, METH_VARARGS, NULL},
	{"XLALUTCToGPS", pylal_XLALUTCToGPS, METH_VARARGS, NULL},
	{NULL,}
};


void initdate(void)
{
	/* commented out to silence warning */
	/*PyObject *module = */Py_InitModule3(MODULE_NAME, module_methods, "Wrapper for LAL's date package.");

	import_array();
	pylal_ligotimegps_import();
}
