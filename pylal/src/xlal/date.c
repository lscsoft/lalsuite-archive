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


/*
 * ============================================================================
 *
 *                              LIGOTimeGPS Type
 *
 * ============================================================================
 */

static PyObject *pylal_LIGOTimeGPS_New(LIGOTimeGPS);


typedef struct {
	PyObject_HEAD
	LIGOTimeGPS gps;
} pylal_LIGOTimeGPS;


static struct PyMemberDef pylal_LIGOTimeGPS_members[] = {
	{"gpsSeconds", T_INT, offsetof(pylal_LIGOTimeGPS, gps.gpsSeconds), 0, "integer seconds"},
	{"gpsNanoSeconds", T_INT, offsetof(pylal_LIGOTimeGPS, gps.gpsNanoSeconds), 0, "integer nanoseconds"},
	{NULL,}
};


static PyObject *pylal_LIGOTimeGPS___abs__(PyObject *self)
{
	LIGOTimeGPS gps;

	XLALINT8NSToGPS(&gps, llabs(XLALGPSToINT8NS(&((pylal_LIGOTimeGPS *) self)->gps)));

	return pylal_LIGOTimeGPS_New(gps);
}


static PyObject *pylal_LIGOTimeGPS___add__(PyObject *self, PyObject *other)
{
	LIGOTimeGPS result;
	LIGOTimeGPS *gps1 = &((pylal_LIGOTimeGPS *) self)->gps;
	LIGOTimeGPS *gps2 = &((pylal_LIGOTimeGPS *) other)->gps;

	XLALINT8NSToGPS(&result, XLALGPSToINT8NS(gps1) + XLALGPSToINT8NS(gps2));

	return pylal_LIGOTimeGPS_New(result);
}


static int pylal_LIGOTimeGPS_compare(PyObject *self, PyObject *other)
{
	LIGOTimeGPS *gps1 = &((pylal_LIGOTimeGPS *) self)->gps;
	LIGOTimeGPS *gps2 = &((pylal_LIGOTimeGPS *) other)->gps;
	int result = gps1->gpsSeconds - gps2->gpsSeconds;
	if(!result)
		result = gps2->gpsNanoSeconds - gps2->gpsNanoSeconds;
	return result;
}


static PyObject *pylal_LIGOTimeGPS___div__(PyObject *self, PyObject *other)
{
	LIGOTimeGPS result;
	PyObject *nanoseconds = PyLong_FromLongLong(XLALGPSToINT8NS(&((pylal_LIGOTimeGPS *) self)->gps));

	PyNumber_InPlaceDivide(nanoseconds, other);
	XLALINT8NSToGPS(&result, PyLong_AsLongLong(nanoseconds));

	Py_DECREF(nanoseconds);

	return pylal_LIGOTimeGPS_New(result);
}


static PyObject *pylal_LIGOTimeGPS___float__(PyObject *self)
{
	LIGOTimeGPS *gps = &((pylal_LIGOTimeGPS *) self)->gps;

	return PyFloat_FromDouble(XLALGPSGetREAL8(gps));
}


static int pylal_LIGOTimeGPS___init__(PyObject *self, PyObject *args, PyObject *kwds)
{
	LIGOTimeGPS *gps = &((pylal_LIGOTimeGPS *) self)->gps;

	if(!PyArg_ParseTuple(args, "ii", &gps->gpsSeconds, &gps->gpsNanoSeconds))
		return -1;

	/* wrap nanoseconds if needed */
	XLALGPSSet(gps, gps->gpsSeconds, gps->gpsNanoSeconds);

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


static PyObject *pylal_LIGOTimeGPS___mod__(PyObject *self, PyObject *other)
{
	LIGOTimeGPS result;
	PyObject *val1 = PyLong_FromLongLong(XLALGPSToINT8NS(&((pylal_LIGOTimeGPS *) self)->gps));
	PyObject *val2 = PyNumber_InPlaceMultiply(PyLong_FromLong(1000000000), other);

	/* nanoseconds %= other * 1000000000 */
	PyNumber_InPlaceRemainder(val1, val2);
	XLALINT8NSToGPS(&result, PyLong_AsLongLong(val1));

	Py_DECREF(val1);
	Py_DECREF(val2);

	return pylal_LIGOTimeGPS_New(result);
}


static PyObject *pylal_LIGOTimeGPS___mul__(PyObject *self, PyObject *other)
{
	LIGOTimeGPS result;
	PyObject *nanoseconds = PyLong_FromLongLong(XLALGPSToINT8NS(&((pylal_LIGOTimeGPS *) self)->gps));

	PyNumber_InPlaceMultiply(nanoseconds, other);
	XLALINT8NSToGPS(&result, PyLong_AsLongLong(nanoseconds));

	Py_DECREF(nanoseconds);

	return pylal_LIGOTimeGPS_New(result);
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


static PyObject *pylal_LIGOTimeGPS___repr__(PyObject *self)
{
	LIGOTimeGPS *gps = &((pylal_LIGOTimeGPS *) self)->gps;

	return PyString_FromFormat("LIGOTimeGPS(%d,%d)", gps->gpsSeconds, gps->gpsNanoSeconds);
}


static PyObject *pylal_LIGOTimeGPS___str__(PyObject *self)
{
	LIGOTimeGPS *gps = &((pylal_LIGOTimeGPS *) self)->gps;
	char str[40];

	/* can't use PyString_FromFormat() because it can't 0-pad the
	 * fractional part */
	sprintf(str, "%d.%09d", gps->gpsSeconds, abs(gps->gpsNanoSeconds));
	return PyString_FromString(str);
}


static PyObject *pylal_LIGOTimeGPS___sub__(PyObject *self, PyObject *other)
{
	LIGOTimeGPS result;
	LIGOTimeGPS *gps1 = &((pylal_LIGOTimeGPS *) self)->gps;
	LIGOTimeGPS *gps2 = &((pylal_LIGOTimeGPS *) other)->gps;

	XLALINT8NSToGPS(&result, XLALGPSToINT8NS(gps1) - XLALGPSToINT8NS(gps2));

	return pylal_LIGOTimeGPS_New(result);
}


static PyNumberMethods pylal_LIGOTimeGPS_as_number = {
	.nb_absolute = pylal_LIGOTimeGPS___abs__,
	.nb_add = pylal_LIGOTimeGPS___add__,
	.nb_divide = pylal_LIGOTimeGPS___div__,
	.nb_float =  pylal_LIGOTimeGPS___float__,
	.nb_int = pylal_LIGOTimeGPS___int__,
	.nb_long = pylal_LIGOTimeGPS___long__,
	.nb_multiply = pylal_LIGOTimeGPS___mul__,
	.nb_negative = pylal_LIGOTimeGPS___neg__,
	.nb_nonzero = pylal_LIGOTimeGPS___nonzero__,
	.nb_positive = pylal_LIGOTimeGPS___pos__,
	.nb_remainder = pylal_LIGOTimeGPS___mod__,
	.nb_subtract = pylal_LIGOTimeGPS___sub__,
};


static PyTypeObject pylal_LIGOTimeGPS_Type = {
	PyObject_HEAD_INIT(NULL)
	.tp_as_number = &pylal_LIGOTimeGPS_as_number,
	.tp_basicsize = sizeof(pylal_LIGOTimeGPS),
	.tp_compare = pylal_LIGOTimeGPS_compare,
	.tp_doc = "A GPS time with nanosecond precision",
	.tp_flags = Py_TPFLAGS_DEFAULT,
	.tp_init = pylal_LIGOTimeGPS___init__,
	.tp_members = pylal_LIGOTimeGPS_members,
	.tp_name = "LIGOTimeGPS",
	.tp_repr = pylal_LIGOTimeGPS___repr__,
	.tp_str = pylal_LIGOTimeGPS___str__,
};


static PyObject *pylal_LIGOTimeGPS_New(LIGOTimeGPS gps)
{
	pylal_LIGOTimeGPS *new = (pylal_LIGOTimeGPS *) _PyObject_New(&pylal_LIGOTimeGPS_Type);

	XLALGPSSet(&new->gps, gps.gpsSeconds, gps.gpsNanoSeconds);

	return (PyObject *) new;
}


/*
 * ============================================================================
 *
 *                       LIGOTimeGPS Function Wrappers
 *
 * ============================================================================
 */

static PyObject *pylal_XLALGPSSetREAL8(PyObject *self, PyObject *args)
{
	pylal_LIGOTimeGPS *s;
	double t;

	if(!PyArg_ParseTuple(args, "Od:XLALGPSSetREAL8", &s, &t))
		return NULL;

	XLALGPSSetREAL8(&s->gps, t);

	return (PyObject *) s;
}


static PyObject *pylal_XLALGPSToINT8NS(PyObject *self, PyObject *args)
{
	pylal_LIGOTimeGPS *s;

	if(!PyArg_ParseTuple(args, "O:XLALGPSToINT8NS", &s))
		return NULL;

	return PyLong_FromLongLong(XLALGPSToINT8NS(&s->gps));
}


static PyObject *pylal_XLALINT8NSToGPS(PyObject *self, PyObject *args)
{
	long long ns;
	pylal_LIGOTimeGPS *new;

	if(!PyArg_ParseTuple(args, "L:XLALINT8NSToGPS", &ns))
		return NULL;

	new = (pylal_LIGOTimeGPS *) _PyObject_New(&pylal_LIGOTimeGPS_Type);
	XLALINT8NSToGPS(&new->gps, ns);

	return (PyObject *) new;
}


/*
 * ============================================================================
 *
 *                              Time Conversion
 *
 * ============================================================================
 */

static PyObject *pylal_XLALGPSToUTC(PyObject *self, PyObject *args)
{
	struct tm tm;
	int gpssec;

	/* arguments:  GPS integer seconds */

	if(!PyArg_ParseTuple(args, "i:XLALGPSToUTC", &gpssec))
		return NULL;

	XLALGPSToUTC(&tm, gpssec);

	/* return a tuple compatible with Python's time.struct_time type:
	 *    - 1900 = year 1900
	 *    - January = month 1
	 *    - Monday = week day 0
	 *    - January 1st = year day 1
	 */

	return Py_BuildValue("(iiiiiiiii)", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec, (tm.tm_wday + 6) % 7, tm.tm_yday + 1, tm.tm_isdst);
}


static PyObject *pylal_XLALUTCToGPS(PyObject *self, PyObject *args)
{
	struct tm tm;

	/* arguments:  9 element tuple compatible with Python's
	 * time.struct_time type. */

	if(!PyArg_ParseTuple(args, "(iiiiiiiii):XLALUTCToGPS", &tm.tm_year, &tm.tm_mon, &tm.tm_mday, &tm.tm_hour, &tm.tm_min, &tm.tm_sec, &tm.tm_wday, &tm.tm_yday, &tm.tm_isdst))
		return NULL;

	/* convert from Python's conventions to C library's */

	tm.tm_year -= 1900;
	tm.tm_mon -= 1;
	tm.tm_wday = (tm.tm_wday + 8) % 7;
	tm.tm_yday -= 1;

	return PyInt_FromLong(XLALUTCToGPS(&tm));
}


/*
 * ============================================================================
 *
 *                            Module Registration
 *
 * ============================================================================
 */

static struct PyMethodDef methods[] = {
	{"XLALGPSSetREAL8", pylal_XLALGPSSetREAL8, 1},
	{"XLALGPSToINT8NS", pylal_XLALGPSToINT8NS, 1},
	{"XLALINT8NSToGPS", pylal_XLALINT8NSToGPS, 1},
	{"XLALGPSToUTC", pylal_XLALGPSToUTC, 1},
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
