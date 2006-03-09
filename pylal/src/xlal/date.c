/*
 * ============================================================================
 *
 *                   Python Wrapper For LAL's Date Package
 *
 * ============================================================================
 */

#include <Python.h>
#include <time.h>
#include <lal/Date.h>


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

	return Py_BuildValue("i", XLALUTCToGPS(&tm));
}


/*
 * registration table
 */

static struct PyMethodDef methods[] = {
	{"XLALGPSToUTC", pylal_XLALGPSToUTC, 1},
	{"XLALUTCToGPS", pylal_XLALUTCToGPS, 1},
	{NULL, NULL}
};

/* module initializer */
void initdate(void)
{
	Py_InitModule("pylal.xlal.date", methods);
}
