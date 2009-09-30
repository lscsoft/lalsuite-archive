/*
 * Copyright (C) 2009  Kipp Cannon
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
 *                 Python Wrapper For LAL's LIGOTimeGPS Type
 *
 * ============================================================================
 */


#include <Python.h>
#include <structmember.h>
#include <ligotimegps.h>
#include <lal/LALDatatypes.h>
#include <lal/Date.h>


#define MODULE_NAME PYLAL_LIGOTIMEGPS_MODULE_NAME


/*
 * ============================================================================
 *
 *                              LIGOTimeGPS Type
 *
 * ============================================================================
 */


/*
 * Utilities
 */


static int pylal_LIGOTimeGPS_Check(PyObject *obj)
{
	return obj ? PyObject_TypeCheck(obj, pylal_LIGOTimeGPS_Type) : 0;
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

	return pylal_LIGOTimeGPS_new(gps);
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

	return pylal_LIGOTimeGPS_new(self_gps);
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

	return pylal_LIGOTimeGPS_new(self_gps);
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

	return pylal_LIGOTimeGPS_new(gps);
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

	return pylal_LIGOTimeGPS_new(gps);
}


static PyObject *pylal_LIGOTimeGPS___neg__(PyObject *self)
{
	LIGOTimeGPS gps;

	if(!pyobject_to_ligotimegps(self, &gps))
		return NULL;

	XLALINT8NSToGPS(&gps, -XLALGPSToINT8NS(&gps));

	return pylal_LIGOTimeGPS_new(gps);
}


static int pylal_LIGOTimeGPS___nonzero__(PyObject *self)
{
	LIGOTimeGPS gps;

	if(!pyobject_to_ligotimegps(self, &gps))
		return -1;

	return gps.gpsSeconds || gps.gpsNanoSeconds;
}


static PyObject *pylal_LIGOTimeGPS_ns(PyObject *self, PyObject *args)
{
	LIGOTimeGPS gps;

	if(!pyobject_to_ligotimegps(self, &gps))
		return NULL;

	return PyLong_FromLong(XLALGPSToINT8NS(&gps));
}


static PyObject *pylal_LIGOTimeGPS___pos__(PyObject *self)
{
	LIGOTimeGPS gps;

	if(!pyobject_to_ligotimegps(self, &gps))
		return NULL;

	return pylal_LIGOTimeGPS_new(gps);
}


static PyObject *pylal_LIGOTimeGPS___reduce__(PyObject *self, PyObject *args)
{
	LIGOTimeGPS gps;

	if(!pyobject_to_ligotimegps(self, &gps))
		return NULL;

	Py_INCREF(pylal_LIGOTimeGPS_Type);
	return Py_BuildValue("(O,(i,i))", pylal_LIGOTimeGPS_Type, gps.gpsSeconds, gps.gpsNanoSeconds);
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

	return pylal_LIGOTimeGPS_new(self_gps);
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
	{"ns", pylal_LIGOTimeGPS_ns, METH_NOARGS, NULL},
	{"__reduce__", pylal_LIGOTimeGPS___reduce__, METH_NOARGS, NULL},
	{NULL,}
};


static PyTypeObject _pylal_LIGOTimeGPS_Type = {
	PyObject_HEAD_INIT(NULL)
	.tp_as_number = &pylal_LIGOTimeGPS_as_number,
	.tp_basicsize = sizeof(pylal_LIGOTimeGPS),
	.tp_doc =
"A GPS time with nanosecond precision.  This is an immutable type.\n" \
"\n" \
"Example:\n" \
"\n" \
">>> LIGOTimeGPS(100.5)\n" \
"LIGOTimeGPS(100, 500000000)\n" \
">>> LIGOTimeGPS(\"100.5\")\n" \
"LIGOTimeGPS(100, 500000000)\n" \
">>> LIGOTimeGPS(100, 500000000)\n" \
"LIGOTimeGPS(100, 500000000)\n" \
">>> LIGOTimeGPS(0, 100500000000L)\n" \
"LIGOTimeGPS(100, 500000000)\n" \
">>> LIGOTimeGPS(100.2, 300000000)\n" \
"LIGOTimeGPS(100, 500000000)\n" \
">>> x = LIGOTimeGPS(100.5)\n" \
">>> x.seconds\n" \
"100\n" \
">>> x.nanoseconds\n" \
"500000000\n" \
">>> str(x)\n" \
"'100.5'\n" \
">>> x + 100\n" \
"LIGOTimeGPS(200,500000000)",
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
 *                            Module Registration
 *
 * ============================================================================
 */


void initligotimegps(void)
{
	PyObject *module = Py_InitModule3(MODULE_NAME, NULL, "Wrapper for LAL's LIGOTimeGPS type.");

	/* LIGOTimeGPS */
	pylal_LIGOTimeGPS_Type = &_pylal_LIGOTimeGPS_Type;
	if(PyType_Ready(pylal_LIGOTimeGPS_Type) < 0)
		return;
	Py_INCREF(pylal_LIGOTimeGPS_Type);
	PyModule_AddObject(module, "LIGOTimeGPS", (PyObject *) pylal_LIGOTimeGPS_Type);
}
