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
#include <numpy/arrayobject.h>
#include <lal/Date.h>
#include <lal/TimeDelay.h>


#define MODULE_NAME "pylal.xlal.date"


/*
 * ============================================================================
 *
 *                              LIGOTimeGPS Type
 *
 * ============================================================================
 */


/*
 * Forward references
 */


static PyTypeObject pylal_LIGOTimeGPS_Type;


/*
 * Structure
 */


typedef struct {
	PyObject_HEAD
	LIGOTimeGPS gps;
} pylal_LIGOTimeGPS;


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
 * Converter function
 */


static int pyobject_to_ligotimegps(PyObject *obj, LIGOTimeGPS *gps)
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
			return 0;
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
			return 0;
		}
	}
	return 1;
}


/*
 * Methods
 */


static PyObject *pylal_LIGOTimeGPS___abs__(PyObject *self)
{
	LIGOTimeGPS gps;

	if(!pyobject_to_ligotimegps(self, &gps))
		return NULL;

	XLALINT8NSToGPS(&gps, llabs(XLALGPSToINT8NS(&gps)));

	return pylal_LIGOTimeGPS_New(gps);
}


static PyObject *pylal_LIGOTimeGPS___add__(PyObject *self, PyObject *other)
{
	LIGOTimeGPS self_gps;
	LIGOTimeGPS other_gps;

	if(!pyobject_to_ligotimegps(self, &self_gps))
		return NULL;
	if(!pyobject_to_ligotimegps(other, &other_gps))
		return NULL;

	XLALGPSAddGPS(&self_gps, &other_gps);

	return pylal_LIGOTimeGPS_New(self_gps);
}


static PyObject *pylal_LIGOTimeGPS___div__(PyObject *self, PyObject *other)
{
	LIGOTimeGPS self_gps;
	/* FIXME:  what about type(other) == LIGOTimeGPS */
	double other_double = PyFloat_AsDouble(other);

	if(PyErr_Occurred())
		return NULL;
	if(!pyobject_to_ligotimegps(self, &self_gps))
		return NULL;

	XLALGPSDivide(&self_gps, other_double);

	return pylal_LIGOTimeGPS_New(self_gps);
}


static PyObject *pylal_LIGOTimeGPS___float__(PyObject *self)
{
	LIGOTimeGPS gps;

	if(!pyobject_to_ligotimegps(self, &gps))
		return NULL;

	return PyFloat_FromDouble(XLALGPSGetREAL8(&gps));
}


static int pylal_LIGOTimeGPS___init__(PyObject *self, PyObject *args, PyObject *kwds)
{
	LIGOTimeGPS *gps = &((pylal_LIGOTimeGPS *) self)->gps;
	PyObject *seconds;
	long long nanoseconds = 0;

	if(!PyArg_ParseTuple(args, "O|L", &seconds, &nanoseconds))
		return -1;

	if(PyUnicode_Check(seconds)) {
		/* convert to ascii string */
		PyObject *str = PyUnicode_AsASCIIString(seconds);
		if(!str)
			return -1;
		Py_DECREF(seconds);
		seconds = str;
	}
	if(PyString_Check(seconds)) {
		char *end, *str = PyString_AsString(seconds);
		int result = XLALStrToGPS(gps, str, &end);
		if((result < 0) || (end == str)) {
			PyErr_SetObject(PyExc_ValueError, seconds);
			return -1;
		}
	} else if(!pyobject_to_ligotimegps(seconds, gps)) {
		PyErr_SetObject(PyExc_ValueError, seconds);
		return -1;
	}

	XLALINT8NSToGPS(gps, XLALGPSToINT8NS(gps) + nanoseconds);

	return 0;
}


static PyObject *pylal_LIGOTimeGPS___int__(PyObject *self)
{
	LIGOTimeGPS gps;

	if(!pyobject_to_ligotimegps(self, &gps))
		return NULL;

	return PyInt_FromLong(gps.gpsSeconds);
}


static PyObject *pylal_LIGOTimeGPS___long__(PyObject *self)
{
	LIGOTimeGPS gps;

	if(!pyobject_to_ligotimegps(self, &gps))
		return NULL;

	return PyLong_FromLong(gps.gpsSeconds);
}


static PyObject *pylal_LIGOTimeGPS___mod__(PyObject *self, PyObject *other)
{
	LIGOTimeGPS gps;
	const double other_double = PyFloat_AsDouble(other);

	if(PyErr_Occurred())
		return NULL;
	if(!pyobject_to_ligotimegps(self, &gps))
		return NULL;

	/* FIXME: loss of precision */
	XLALINT8NSToGPS(&gps, XLALGPSToINT8NS(&gps) % (long long) (other_double * 1e9));

	return pylal_LIGOTimeGPS_New(gps);
}


static PyObject *pylal_LIGOTimeGPS___mul__(PyObject *self, PyObject *other)
{
	LIGOTimeGPS gps;
	double factor;

	if(pylal_LIGOTimeGPS_Check(self) && !pylal_LIGOTimeGPS_Check(other)) {
		gps = ((pylal_LIGOTimeGPS *) self)->gps;
		factor = PyFloat_AsDouble(other);
	} else if(!pylal_LIGOTimeGPS_Check(self) && pylal_LIGOTimeGPS_Check(other)) {
		gps = ((pylal_LIGOTimeGPS *) other)->gps;
		factor = PyFloat_AsDouble(self);
	} else {
		Py_INCREF(Py_NotImplemented);
		return Py_NotImplemented;
	}
	if(PyErr_Occurred())
		return NULL;

	XLALGPSMultiply(&gps, factor);

	return pylal_LIGOTimeGPS_New(gps);
}


static PyObject *pylal_LIGOTimeGPS___neg__(PyObject *self)
{
	LIGOTimeGPS gps;

	if(!pyobject_to_ligotimegps(self, &gps))
		return NULL;

	XLALINT8NSToGPS(&gps, -XLALGPSToINT8NS(&gps));

	return pylal_LIGOTimeGPS_New(gps);
}


static int pylal_LIGOTimeGPS___nonzero__(PyObject *self)
{
	LIGOTimeGPS gps;

	if(!pyobject_to_ligotimegps(self, &gps))
		return -1;

	return gps.gpsSeconds || gps.gpsNanoSeconds;
}


static PyObject *pylal_LIGOTimeGPS___pos__(PyObject *self)
{
	LIGOTimeGPS gps;

	if(!pyobject_to_ligotimegps(self, &gps))
		return NULL;

	return pylal_LIGOTimeGPS_New(gps);
}


static PyObject *pylal_LIGOTimeGPS___reduce__(PyObject *self, PyObject *args)
{
	LIGOTimeGPS gps;

	if(!pyobject_to_ligotimegps(self, &gps))
		return NULL;

	Py_INCREF(&pylal_LIGOTimeGPS_Type);
	return Py_BuildValue("(O,(i,i))", &pylal_LIGOTimeGPS_Type, gps.gpsSeconds, gps.gpsNanoSeconds);
}


static PyObject *pylal_LIGOTimeGPS___repr__(PyObject *self)
{
	LIGOTimeGPS gps;

	if(!pyobject_to_ligotimegps(self, &gps))
		return NULL;

	return PyString_FromFormat("LIGOTimeGPS(%d,%d)", gps.gpsSeconds, gps.gpsNanoSeconds);
}


static PyObject *pylal_LIGOTimeGPS_richcompare(PyObject *self, PyObject *other, int op_id)
{
	LIGOTimeGPS self_gps;
	LIGOTimeGPS other_gps;
	int d;
	PyObject *result;

	if(!pyobject_to_ligotimegps(self, &self_gps))
		return NULL;
	if(!pyobject_to_ligotimegps(other, &other_gps)) {
		PyErr_Clear();
		Py_INCREF(Py_NotImplemented);
		return Py_NotImplemented;
	}

	d = XLALGPSCmp(&self_gps, &other_gps);
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


static long pylal_LIGOTimeGPS_hash(PyObject *self)
{
	LIGOTimeGPS *gps = &((pylal_LIGOTimeGPS *) self)->gps;
	long hash = (long) gps->gpsSeconds ^ (long) gps->gpsNanoSeconds;
	return hash == -1 ? -2 : hash;
}


static PyObject *pylal_LIGOTimeGPS___str__(PyObject *self)
{
	LIGOTimeGPS gps;
	char str[22];
	int i;

	if(!pyobject_to_ligotimegps(self, &gps))
		return NULL;

	if(gps.gpsNanoSeconds) {
		if((gps.gpsSeconds == 0) && (gps.gpsNanoSeconds < 0))
			snprintf(str, 21, "-0.%09u", abs(gps.gpsNanoSeconds));
		else
			snprintf(str, 21, "%d.%09u", gps.gpsSeconds, abs(gps.gpsNanoSeconds));
		for(i = strlen(str); str[--i] == '0'; str[i] = '\0');
	} else
		snprintf(str, 21, "%d", gps.gpsSeconds);

	return PyString_FromString(str);
}


static PyObject *pylal_LIGOTimeGPS___sub__(PyObject *self, PyObject *other)
{
	LIGOTimeGPS self_gps;
	LIGOTimeGPS other_gps;

	if(!pyobject_to_ligotimegps(self, &self_gps))
		return NULL;
	if(!pyobject_to_ligotimegps(other, &other_gps))
		return NULL;

	XLALINT8NSToGPS(&self_gps, XLALGPSToINT8NS(&self_gps) - XLALGPSToINT8NS(&other_gps));

	return pylal_LIGOTimeGPS_New(self_gps);
}


/*
 * Type information
 */


static struct PyMemberDef pylal_LIGOTimeGPS_members[] = {
	{"seconds", T_INT, offsetof(pylal_LIGOTimeGPS, gps.gpsSeconds), 0, "integer seconds"},
	{"nanoseconds", T_INT, offsetof(pylal_LIGOTimeGPS, gps.gpsNanoSeconds), 0, "integer nanoseconds"},
	{NULL,}
};


static PyNumberMethods pylal_LIGOTimeGPS_as_number = {
	.nb_absolute = pylal_LIGOTimeGPS___abs__,
	.nb_add = pylal_LIGOTimeGPS___add__,
	.nb_divide = pylal_LIGOTimeGPS___div__,
	.nb_float =  pylal_LIGOTimeGPS___float__,
	.nb_int = pylal_LIGOTimeGPS___int__,
	.nb_long = pylal_LIGOTimeGPS___long__,
	.nb_remainder = pylal_LIGOTimeGPS___mod__,
	.nb_multiply = pylal_LIGOTimeGPS___mul__,
	.nb_negative = pylal_LIGOTimeGPS___neg__,
	.nb_nonzero = pylal_LIGOTimeGPS___nonzero__,
	.nb_positive = pylal_LIGOTimeGPS___pos__,
	.nb_subtract = pylal_LIGOTimeGPS___sub__,
};


static struct PyMethodDef pylal_LIGOTimeGPS_methods[] = {
	{"__reduce__", pylal_LIGOTimeGPS___reduce__, METH_NOARGS, NULL},
	{NULL,}
};


static PyTypeObject pylal_LIGOTimeGPS_Type = {
	PyObject_HEAD_INIT(NULL)
	.tp_as_number = &pylal_LIGOTimeGPS_as_number,
	.tp_basicsize = sizeof(pylal_LIGOTimeGPS),
	.tp_doc = "A GPS time with nanosecond precision",
	.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_CHECKTYPES,
	.tp_init = pylal_LIGOTimeGPS___init__,
	.tp_members = pylal_LIGOTimeGPS_members,
	.tp_methods = pylal_LIGOTimeGPS_methods,
	.tp_name = MODULE_NAME ".LIGOTimeGPS",
	.tp_new = PyType_GenericNew,
	.tp_repr = pylal_LIGOTimeGPS___repr__,
	.tp_richcompare = pylal_LIGOTimeGPS_richcompare,
	.tp_hash = pylal_LIGOTimeGPS_hash,
	.tp_str = pylal_LIGOTimeGPS___str__,
};


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

	XLALGPSSet(&gps, XLALUTCToGPS(struct_tm_python_to_c(&utc)), 0);

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
	if(!PyArg_ParseTuple(args, "O!d:XLALGreenwichSiderealTime", &pylal_LIGOTimeGPS_Type, &gps, &equation_of_equinoxes))
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
	return pylal_LIGOTimeGPS_New(gps);
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
	if(!PyArg_ParseTuple(args, "OOddO!:XLALArrivalTimeDiff", &pos1, &pos2, &ra, &dec, &pylal_LIGOTimeGPS_Type, &gps))
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
	PyObject *module = Py_InitModule3(MODULE_NAME, module_methods, "Wrapper for LAL's date package.");

	import_array();

	/* LIGOTimeGPS */
	if(PyType_Ready(&pylal_LIGOTimeGPS_Type) < 0)
		return;
	Py_INCREF(&pylal_LIGOTimeGPS_Type);
	PyModule_AddObject(module, "LIGOTimeGPS", (PyObject *) &pylal_LIGOTimeGPS_Type);
}
