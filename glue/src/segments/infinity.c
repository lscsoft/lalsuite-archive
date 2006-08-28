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
 *                Segments Module Component --- Infinity Class
 *
 * ============================================================================
 */

#include <Python.h>
#include <structmember.h>
#include <stdlib.h>

#include <segments.h>


/*
 * ============================================================================
 *
 *                               Infinity Class
 *
 * ============================================================================
 */


/*
 * Preallocated instances
 */

segments_Infinity *segments_PosInfinity;
segments_Infinity *segments_NegInfinity;


/*
 * Utilities
 */


static int segments_Infinity_Check(PyObject *obj)
{
	return obj ? PyObject_TypeCheck(obj, &segments_Infinity_Type) : 0;
}


/*
 * Methods
 */

static PyObject *__new__(PyTypeObject *subtype, PyObject *args, PyObject *kwds)
{
	Py_INCREF(segments_PosInfinity);
	return (PyObject *) segments_PosInfinity;
}


static PyObject *__add__(PyObject *self, PyObject *other)
{
	if(segments_Infinity_Check(self)) {
		/* __add__ case */
		Py_INCREF(self);
		return self;
	} else if(segments_Infinity_Check(other)) {
		/* __radd__ case */
		Py_INCREF(other);
		return other;
	}
	PyErr_SetObject(PyExc_TypeError, other);
	return NULL;
}


static PyObject *__neg__(PyObject *self)
{
	PyObject *result;

	if(!segments_Infinity_Check(self)) {
		PyErr_SetObject(PyExc_TypeError, self);
		return NULL;
	}

	if(self == (PyObject *) segments_PosInfinity)
		result = (PyObject *) segments_NegInfinity;
	else
		result = (PyObject *) segments_PosInfinity;
	Py_INCREF(result);
	return result;
}


static int __nonzero__(PyObject *self)
{
	if(segments_Infinity_Check(self))
		return 1;
	PyErr_SetObject(PyExc_TypeError, self);
	return 0;
}


static PyObject *__pos__(PyObject *self)
{
	if(segments_Infinity_Check(self)) {
		Py_INCREF(self);
		return self;
	}
	PyErr_SetObject(PyExc_TypeError, self);
	return NULL;
}


static PyObject *__repr__(PyObject *self)
{
	return PyString_FromString(self == (PyObject *) segments_PosInfinity ? "infinity" : "-infinity");
}


static PyObject *__reduce__(PyObject *self, PyObject *args)
{
	if(!segments_Infinity_Check(self)) {
		PyErr_SetObject(PyExc_TypeError, self);
		return NULL;
	}

	Py_INCREF(&segments_Infinity_Type);
	/* FIXME */
	return Py_BuildValue("(O,())", &segments_Infinity_Type);
}


static PyObject *richcompare(PyObject *self, PyObject *other, int op_id)
{
	int s = segments_Infinity_Check(self) ? self == (PyObject *) segments_PosInfinity ? +1 : -1 : 0;
	int o = segments_Infinity_Check(other) ? other == (PyObject *) segments_PosInfinity ? +1 : -1 : 0;
	int d = s - o;
	PyObject *result;

	if(!s && !o) {
		PyErr_SetObject(PyExc_TypeError, other);
		return NULL;
	}

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


static PyObject *__sub__(PyObject *self, PyObject *other)
{
	PyObject *result;

	if(segments_Infinity_Check(self)) {
		/* __sub__ case */
		if(self == other) {
			PyErr_SetObject(PyExc_ValueError, other);
			return NULL;
		}
		result = self;
	} else {
		/* __rsub__ case */
		if(!segments_Infinity_Check(other)) {
			PyErr_SetObject(PyExc_TypeError, self);
			return NULL;
		}
		if(other == (PyObject *) segments_PosInfinity)
			result = (PyObject *) segments_NegInfinity;
		else
			result = (PyObject *) segments_PosInfinity;
	}
	Py_INCREF(result);
	return result;
}


/*
 * Type information
 */

static PyNumberMethods as_number = {
	.nb_add = __add__,
	.nb_negative = __neg__,
	.nb_nonzero = __nonzero__,
	.nb_positive = __pos__,
	.nb_subtract = __sub__,
};


static struct PyMethodDef methods[] = {
	{"__reduce__", __reduce__, METH_NOARGS, "pickle helper"},
	{NULL,}
};


PyTypeObject segments_Infinity_Type = {
	PyObject_HEAD_INIT(NULL)
	.tp_as_number = &as_number,
	.tp_basicsize = sizeof(segments_Infinity),
	.tp_doc =
	"The infinity object possess the algebraic properties necessary for\n" \
	"use as a bound on semi-infinite and infinite segments.\n" \
	"\n" \
	"Example:\n" \
	"\n" \
	">>> x = infinity()\n" \
	">>> x > 0\n" \
	"True\n" \
	">>> x + 10\n" \
	"infinity\n" \
	">>> segment(-10, 10) - segment(-x, 0)\n" \
	"segment(0, 10)",
	.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_CHECKTYPES,
	.tp_methods = methods,
	.tp_name = MODULE_NAME ".infinity",
	.tp_new = __new__,
	.tp_repr = __repr__,
	.tp_richcompare = richcompare,
	.tp_str = __repr__,
};
