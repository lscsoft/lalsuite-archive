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
	XLALGPSSet(gps, gps->gpsSeconds + seconds, gps->gpsNanoSeconds);

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
	.tp_compare = pylal_LIGOTimeGPS_compare,
	.tp_doc = "A GPS time with nanosecond precision",
	.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
	/*.tp_getset = pylal_LIGOTimeGPS_getset,*/
	.tp_init = pylal_LIGOTimeGPS___init__,
	.tp_members = pylal_LIGOTimeGPS_members,
	.tp_name = "LIGOTimeGPS",
	.tp_repr = pylal_LIGOTimeGPS___repr__,
	.tp_str = pylal_LIGOTimeGPS___str__,
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

	/* LIGOTimeGPS, float */
	if(!PyArg_ParseTuple(args, "Od:XLALGPSSetREAL8", &s, &t))
		return NULL;

	XLALGPSSetREAL8(&s->gps, t);

	/* LIGOTimeGPS */
	return (PyObject *) s;
}


static PyObject *pylal_XLALGPSToINT8NS(PyObject *self, PyObject *args)
{
	pylal_LIGOTimeGPS *s;

	/* LIGOTimeGPS */
	if(!PyArg_ParseTuple(args, "O:XLALGPSToINT8NS", &s))
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
	char *str, *end;

	/* string */
	if(!PyArg_ParseTuple(args, "s:XLALStrToGPS", &str))
		return NULL;

	XLALStrToGPS(&gps, str, &end);

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
	struct tm utc;
	int gpssec;

	/* int */
	if(!PyArg_ParseTuple(args, "i:XLALGPSToUTC", &gpssec))
		return NULL;

	XLALGPSToUTC(&utc, gpssec);

	struct_tm_c_to_python(&utc);

	/* time.struct_time */
	return Py_BuildValue("(iiiiiiiii)", utc.tm_year, utc.tm_mon, utc.tm_mday, utc.tm_hour, utc.tm_min, utc.tm_sec, utc.tm_wday, utc.tm_yday, utc.tm_isdst);
}


static PyObject *pylal_XLALUTCToGPS(PyObject *self, PyObject *args)
{
	struct tm utc;

	/* time.struct_time */
	if(!PyArg_ParseTuple(args, "(iiiiiiiii):XLALUTCToGPS", &utc.tm_year, &utc.tm_mon, &utc.tm_mday, &utc.tm_hour, &utc.tm_min, &utc.tm_sec, &utc.tm_wday, &utc.tm_yday, &utc.tm_isdst))
		return NULL;

	struct_tm_python_to_c(&utc);

	/* int */
	return PyInt_FromLong(XLALUTCToGPS(&utc));
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
	if(!PyArg_ParseTuple(args, "Od:XLALGreenwichSiderealTime", &gps, &equation_of_equinoxes))
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
	if(!PyArg_ParseTuple(args, "OOddO:XLALArrivalTimeDiff", &pos1_list, &pos2_list, &ra, &dec, &gps))
		return NULL;

	for(i = 0; i < 3; i++) {
		pos1[i] = PyFloat_AsDouble(PyList_GetItem(pos1_list, i));
		pos2[i] = PyFloat_AsDouble(PyList_GetItem(pos2_list, i));
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
	{"XLALGPSSetREAL8", pylal_XLALGPSSetREAL8, 1},
	{"XLALGPSToINT8NS", pylal_XLALGPSToINT8NS, 1},
	{"XLALGPSToUTC", pylal_XLALGPSToUTC, 1},
	{"XLALGreenwichSiderealTime", pylal_XLALGreenwichSiderealTime, 1},
	{"XLALINT8NSToGPS", pylal_XLALINT8NSToGPS, 1},
	{"XLALJulianDay", pylal_XLALJulianDay, 1},
	{"XLALLeapSeconds", pylal_XLALLeapSeconds, 1},
	{"XLALLeapSecondsUTC", pylal_XLALLeapSecondsUTC, 1},
	{"XLALModifiedJulianDay", pylal_XLALModifiedJulianDay, 1},
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
