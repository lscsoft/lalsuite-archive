/*
 * $Id$
 *
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
#include <lal/Date.h>
#include <lal/TimeDelay.h>


/*
 * ============================================================================
 *
 *                              LIGOTimeGPS Type
 *
 * ============================================================================
 */

/*
 * forward references
 */

static PyObject *pylal_LIGOTimeGPS_New(LIGOTimeGPS);
static int pylal_LIGOTimeGPS_Check(PyObject *);


/*
 * Structure
 */

typedef struct {
	PyObject_HEAD
	LIGOTimeGPS gps;
} pylal_LIGOTimeGPS;


/*
 * Member access
 */

static struct PyMemberDef pylal_LIGOTimeGPS_members[] = {
	{"seconds", T_INT, offsetof(pylal_LIGOTimeGPS, gps.gpsSeconds), 0, "integer seconds"},
	{"nanoseconds", T_INT, offsetof(pylal_LIGOTimeGPS, gps.gpsNanoSeconds), 0, "integer nanoseconds"},
	{NULL,}
};


/*
 * Type conversion.
 */

static LIGOTimeGPS *pyobject_to_ligotimegps(LIGOTimeGPS *gps, PyObject *obj)
{
	if(pylal_LIGOTimeGPS_Check(obj)) {
		*gps = ((pylal_LIGOTimeGPS *) obj)->gps;
	} else if(PyInt_Check(obj)) {
		XLALGPSSet(gps, PyInt_AsLong(obj), 0);
	} else if(PyLong_Check(obj)) {
		XLALGPSSet(gps, PyLong_AsLongLong(obj), 0);
	} else if(PyFloat_Check(obj)) {
		XLALGPSSetREAL8(gps, PyFloat_AsDouble(obj));
	} else if(PyComplex_Check(obj)) {
		if(PyComplex_ImagAsDouble(obj) != 0.0) {
			XLALGPSSet(gps, 0, 0);
			PyErr_SetObject(PyExc_ValueError, obj);
			return NULL;
		}
		XLALGPSSetREAL8(gps, PyComplex_RealAsDouble(obj));
	} else {
		PyObject *s_attr = PyObject_GetAttrString(obj, "seconds");
		PyObject *n_attr = PyObject_GetAttrString(obj, "nanoseconds");
		XLALGPSSet(gps, PyInt_AsLong(s_attr), PyInt_AsLong(n_attr));
		Py_XDECREF(s_attr);
		Py_XDECREF(n_attr);
		if(PyErr_Occurred()) {
			PyErr_SetObject(PyExc_TypeError, obj);
			return NULL;
		}
	}
	return gps;
}


/*
 * Methods
 */

static PyObject *pylal_LIGOTimeGPS___abs__(PyObject *self)
{
	LIGOTimeGPS gps;

	XLALINT8NSToGPS(&gps, llabs(XLALGPSToINT8NS(&((pylal_LIGOTimeGPS *) self)->gps)));

	return pylal_LIGOTimeGPS_New(gps);
}


static PyObject *pylal_LIGOTimeGPS___add__(PyObject *self, PyObject *other)
{
	LIGOTimeGPS gps;

	if(!pyobject_to_ligotimegps(&gps, other))
		return NULL;

	XLALINT8NSToGPS(&gps, XLALGPSToINT8NS(&((pylal_LIGOTimeGPS *) self)->gps) + XLALGPSToINT8NS(&gps));

	return pylal_LIGOTimeGPS_New(gps);
}


static PyObject *pylal_LIGOTimeGPS___float__(PyObject *self)
{
	LIGOTimeGPS *gps = &((pylal_LIGOTimeGPS *) self)->gps;

	return PyFloat_FromDouble(XLALGPSGetREAL8(gps));
}


static int pylal_LIGOTimeGPS___init__(PyObject *self, PyObject *args, PyObject *kwds)
{
	LIGOTimeGPS *gps = &((pylal_LIGOTimeGPS *) self)->gps;
	int seconds;
	long long nanoseconds;

	if(!PyArg_ParseTuple(args, "iL", &seconds, &nanoseconds))
		return -1;

	XLALINT8NSToGPS(gps, nanoseconds);
	XLALGPSSet(gps, seconds + gps->gpsSeconds, gps->gpsNanoSeconds);

	return 0;
}


static PyObject *pylal_LIGOTimeGPS___int__(PyObject *self)
{
	LIGOTimeGPS *gps = &((pylal_LIGOTimeGPS *) self)->gps;

	return PyInt_FromLong(gps->gpsSeconds);
}


static PyObject *pylal_LIGOTimeGPS___long__(PyObject *self)
{
	LIGOTimeGPS *gps = &((pylal_LIGOTimeGPS *) self)->gps;

	return PyLong_FromLong(gps->gpsSeconds);
}


static PyObject *pylal_LIGOTimeGPS___neg__(PyObject *self)
{
	LIGOTimeGPS gps;

	XLALINT8NSToGPS(&gps, -XLALGPSToINT8NS(&((pylal_LIGOTimeGPS *) self)->gps));

	return pylal_LIGOTimeGPS_New(gps);
}


static int pylal_LIGOTimeGPS___nonzero__(PyObject *self)
{
	LIGOTimeGPS *gps = &((pylal_LIGOTimeGPS *) self)->gps;

	return gps->gpsSeconds || gps->gpsNanoSeconds;
}


static PyObject *pylal_LIGOTimeGPS___pos__(PyObject *self)
{
	return pylal_LIGOTimeGPS_New(((pylal_LIGOTimeGPS *) self)->gps);
}


static PyObject *pylal_LIGOTimeGPS__richcompare__(PyObject *self, PyObject *other, int op_id)
{
	LIGOTimeGPS gps;
	int d;
	PyObject *result;

	if(!pyobject_to_ligotimegps(&gps, other))
		return NULL;

	d = XLALGPSCmp(&((pylal_LIGOTimeGPS *) self)->gps, &gps);
	switch(op_id) {
	case Py_LT:
		result = (d < 0) ? Py_True : Py_False;
		break;

	case Py_LE:
		result = (d <= 0) ? Py_True : Py_False;
		break;

	case Py_EQ:
		result = (d == 0) ? Py_True : Py_False;
		break;

	case Py_NE:
		result = (d != 0) ? Py_True : Py_False;
		break;

	case Py_GT:
		result = (d > 0) ? Py_True : Py_False;
		break;

	case Py_GE:
		result = (d >= 0) ? Py_True : Py_False;
		break;

	default:
		PyErr_BadInternalCall();
		return NULL;
	}

	Py_INCREF(result);
	return result;
}


static PyObject *pylal_LIGOTimeGPS___sub__(PyObject *self, PyObject *other)
{
	LIGOTimeGPS gps;

	if(!pyobject_to_ligotimegps(&gps, other))
		return NULL;

	XLALINT8NSToGPS(&gps, XLALGPSToINT8NS(&((pylal_LIGOTimeGPS *) self)->gps) - XLALGPSToINT8NS(&gps));

	return pylal_LIGOTimeGPS_New(gps);
}


/*
 * Type information
 */

static PyNumberMethods pylal_LIGOTimeGPS_as_number = {
	.nb_absolute = pylal_LIGOTimeGPS___abs__,
	.nb_add = pylal_LIGOTimeGPS___add__,
	.nb_float =  pylal_LIGOTimeGPS___float__,
	.nb_int = pylal_LIGOTimeGPS___int__,
	.nb_long = pylal_LIGOTimeGPS___long__,
	.nb_negative = pylal_LIGOTimeGPS___neg__,
	.nb_nonzero = pylal_LIGOTimeGPS___nonzero__,
	.nb_positive = pylal_LIGOTimeGPS___pos__,
	.nb_subtract = pylal_LIGOTimeGPS___sub__,
};


static PyTypeObject pylal_LIGOTimeGPS_Type = {
	PyObject_HEAD_INIT(NULL)
	.tp_as_number = &pylal_LIGOTimeGPS_as_number,
	.tp_basicsize = sizeof(pylal_LIGOTimeGPS),
	.tp_doc = "A GPS time with nanosecond precision",
	.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_CHECKTYPES,
	.tp_init = pylal_LIGOTimeGPS___init__,
	.tp_members = pylal_LIGOTimeGPS_members,
	.tp_name = "LIGOTimeGPS",
	.tp_richcompare = pylal_LIGOTimeGPS__richcompare__,
};


/*
 * Utilities
 */

static PyObject *pylal_LIGOTimeGPS_New(LIGOTimeGPS gps)
{
	pylal_LIGOTimeGPS *new = (pylal_LIGOTimeGPS *) _PyObject_New(&pylal_LIGOTimeGPS_Type);

	XLALGPSSet(&new->gps, gps.gpsSeconds, gps.gpsNanoSeconds);

	return (PyObject *) new;
}


static int pylal_LIGOTimeGPS_Check(PyObject *obj)
{
	return obj ? PyObject_TypeCheck(obj, &pylal_LIGOTimeGPS_Type) : 0;
}


/*
 * ============================================================================
 *
 *                       LIGOTimeGPS Function Wrappers
 *
 * ============================================================================
 */

static PyObject *pylal_XLALREAL8ToGPS(PyObject *self, PyObject *args)
{
	pylal_LIGOTimeGPS *new;
	double t;

	/* float */
	if(!PyArg_ParseTuple(args, "d:XLALREAL8ToGPS", &t))
		return NULL;

	new = (pylal_LIGOTimeGPS *) _PyObject_New(&pylal_LIGOTimeGPS_Type);
	XLALGPSSetREAL8(&new->gps, t);

	/* LIGOTimeGPS */
	return (PyObject *) new;
}


static PyObject *pylal_XLALGPSToINT8NS(PyObject *self, PyObject *args)
{
	pylal_LIGOTimeGPS *s;

	/* LIGOTimeGPS */
	if(!PyArg_ParseTuple(args, "O!:XLALGPSToINT8NS", &pylal_LIGOTimeGPS_Type, &s))
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

	new = (pylal_LIGOTimeGPS *) _PyObject_New(&pylal_LIGOTimeGPS_Type);
	XLALINT8NSToGPS(&new->gps, ns);

	/* LIGOTimeGPS */
	return (PyObject *) new;
}


static PyObject *pylal_XLALStrToGPS(PyObject *self, PyObject *args)
{
	LIGOTimeGPS gps;
	int result;
	char *str, *end;

	/* string */
	if(!PyArg_ParseTuple(args, "s:XLALStrToGPS", &str))
		return NULL;

	result = XLALStrToGPS(&gps, str, &end);
	if((result < 0) || (end == str)) {
		PyErr_SetString(PyExc_ValueError, str);
		return NULL;
	}

	/* LIGOTimeGPS */
	return pylal_LIGOTimeGPS_New(gps);
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

static void struct_tm_python_to_c(struct tm *tm)
{
	tm->tm_year -= 1900;
	tm->tm_mon -= 1;
	tm->tm_wday = (tm->tm_wday + 8) % 7;
	tm->tm_yday -= 1;
}


static void struct_tm_c_to_python(struct tm *tm)
{
	tm->tm_year += 1900;
	tm->tm_mon += 1;
	tm->tm_wday = (tm->tm_wday + 6) % 7;
	tm->tm_yday += 1;
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

	struct_tm_python_to_c(&utc);

	/* int */
	return PyInt_FromLong(XLALLeapSecondsUTC(&utc));
}


/*
 * GPS to/from UTC
 */

static PyObject *pylal_XLALGPSToUTC(PyObject *self, PyObject *args)
{
	pylal_LIGOTimeGPS *gps;
	struct tm utc;

	/* LIGOTimeGPS */
	if(!PyArg_ParseTuple(args, "O!:XLALGPSToUTC", &pylal_LIGOTimeGPS_Type, &gps))
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

	struct_tm_python_to_c(&utc);
	XLALGPSSet(&gps, XLALUTCToGPS(&utc), 0);

	/* LIGOTimeGPS */
	return pylal_LIGOTimeGPS_New(gps);
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

	struct_tm_python_to_c(&utc);

	/* float */
	return PyFloat_FromDouble(XLALJulianDay(&utc));
}


static PyObject *pylal_XLALModifiedJulianDay(PyObject *self, PyObject *args)
{
	struct tm utc;

	/* time.struct_time */
	if(!PyArg_ParseTuple(args, "(iiiiiiiii):XLALModifiedJulianDay", &utc.tm_year, &utc.tm_mon, &utc.tm_mday, &utc.tm_hour, &utc.tm_min, &utc.tm_sec, &utc.tm_wday, &utc.tm_yday, &utc.tm_isdst))
		return NULL;

	struct_tm_python_to_c(&utc);

	/* int */
	return PyInt_FromLong(XLALModifiedJulianDay(&utc));
}


/*
 * Sidereal time
 */

static PyObject *pylal_XLALGreenwichSiderealTime(PyObject *self, PyObject *args)
{
	pylal_LIGOTimeGPS *gps;
	double equation_of_equinoxes;

	/* LIGOTimeGPS, float */
	if(!PyArg_ParseTuple(args, "O!d:XLALGreenwichSiderealTime", &pylal_LIGOTimeGPS_Type, &gps, &equation_of_equinoxes))
		return NULL;

	/* float */
	return PyFloat_FromDouble(XLALGreenwichSiderealTime(&gps->gps, equation_of_equinoxes));
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
	PyObject *pos1_list, *pos2_list;
	double pos1[3], pos2[3];
	double ra, dec;
	pylal_LIGOTimeGPS *gps;
	int i;

	/* 3-element list, 3-element list, float, float, LIGOTimeGPS */
	if(!PyArg_ParseTuple(args, "O!O!ddO!:XLALArrivalTimeDiff", &PyList_Type, &pos1_list, &PyList_Type, &pos2_list, &ra, &dec, &pylal_LIGOTimeGPS_Type, &gps))
		return NULL;
	if(PyList_Size(pos1_list) != 3) {
		PyErr_Format(PyExc_TypeError, "XLALArrivalTimeDiff() argument 1 must have length 3, not %d", PyList_Size(pos1_list));
		return NULL;
	}
	if(PyList_Size(pos2_list) != 3) {
		PyErr_Format(PyExc_TypeError, "XLALArrivalTimeDiff() argument 2 must have length 3, not %d", PyList_Size(pos1_list));
		return NULL;
	}

	for(i = 0; i < 3; i++) {
		pos1[i] = PyFloat_AsDouble(PyList_GetItem(pos1_list, i));
		if(PyErr_Occurred())
			return NULL;
		pos2[i] = PyFloat_AsDouble(PyList_GetItem(pos2_list, i));
		if(PyErr_Occurred())
			return NULL;
	}

	/* float */
	return PyFloat_FromDouble(XLALArrivalTimeDiff(pos1, pos2, ra, dec, &gps->gps));
}


/*
 * ============================================================================
 *
 *                            Module Registration
 *
 * ============================================================================
 */

static struct PyMethodDef methods[] = {
	{"XLALArrivalTimeDiff", pylal_XLALArrivalTimeDiff, 1},
	{"XLALGPSToINT8NS", pylal_XLALGPSToINT8NS, 1},
	{"XLALGPSToUTC", pylal_XLALGPSToUTC, 1},
	{"XLALGreenwichSiderealTime", pylal_XLALGreenwichSiderealTime, 1},
	{"XLALINT8NSToGPS", pylal_XLALINT8NSToGPS, 1},
	{"XLALJulianDay", pylal_XLALJulianDay, 1},
	{"XLALLeapSeconds", pylal_XLALLeapSeconds, 1},
	{"XLALLeapSecondsUTC", pylal_XLALLeapSecondsUTC, 1},
	{"XLALModifiedJulianDay", pylal_XLALModifiedJulianDay, 1},
	{"XLALREAL8ToGPS", pylal_XLALREAL8ToGPS, 1},
	{"XLALStrToGPS", pylal_XLALStrToGPS, 1},
	{"XLALUTCToGPS", pylal_XLALUTCToGPS, 1},
	{NULL,}
};


void initdate(void)
{
	PyObject *m = Py_InitModule3("pylal.xlal.date", methods, "Wrapper for LAL's date package.");

	/* LIGOTimeGPS */
	pylal_LIGOTimeGPS_Type.tp_new = PyType_GenericNew;
	if(PyType_Ready(&pylal_LIGOTimeGPS_Type) < 0)
		return;
	Py_INCREF(&pylal_LIGOTimeGPS_Type);
	PyModule_AddObject(m, "LIGOTimeGPS", (PyObject *) &pylal_LIGOTimeGPS_Type);
}
