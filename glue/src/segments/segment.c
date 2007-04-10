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
 *                Segments Module Component --- segment Class
 *
 * ============================================================================
 */

#include <Python.h>
#include <stdlib.h>

#include <segments.h>


/*
 * ============================================================================
 *
 *                               segment Class
 *
 * ============================================================================
 */


/*
 * Utilities
 */


static int segments_Segment_Check(PyObject *obj)
{
	return obj ? PyObject_TypeCheck(obj, &segments_Segment_Type) : 0;
}


/*
 * Basic methods
 */


static PyObject *segments_Segment_New(PyTypeObject *type, PyObject *a, PyObject *b)
{
	PyObject *new = type->tp_alloc(type, 2);
	PyTuple_SET_ITEM(new, 0, a);
	PyTuple_SET_ITEM(new, 1, b);
	return new;
}


static PyObject *__new__(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
	PyObject *a, *b;
	int delta;

	if(!PyArg_ParseTuple(args, "OO", &a, &b))
		if(!PyArg_ParseTuple(args, "(OO)", &a, &b)) {
			PyErr_SetString(PyExc_TypeError, "__new__() takes 2 arguments, or 1 arguments when it is a sequence of length 2");
			return NULL;
		}

	Py_INCREF(a);
	Py_INCREF(b);

	PyErr_Clear();
	delta = PyObject_Compare(a, b);
	if(PyErr_Occurred()) {
		Py_DECREF(a);
		Py_DECREF(b);
		return NULL;
	}

	if(delta <= 0)
		return segments_Segment_New(type, a, b);
	return segments_Segment_New(type, b, a);
}


static PyObject *__repr__(PyObject *self)
{
	PyObject *a = PyObject_Repr(PyTuple_GET_ITEM(self, 0));
	PyObject *b = PyObject_Repr(PyTuple_GET_ITEM(self, 1));
	PyObject *result = PyString_FromFormat("segment(%s, %s)", PyString_AsString(a), PyString_AsString(b));
	Py_DECREF(a);
	Py_DECREF(b);
	return result;
}


static PyObject *__str__(PyObject *self)
{
	PyObject *a = PyObject_Str(PyTuple_GET_ITEM(self, 0));
	PyObject *b = PyObject_Str(PyTuple_GET_ITEM(self, 1));
	PyObject *result = PyString_FromFormat("[%s ... %s)", PyString_AsString(a), PyString_AsString(b));
	Py_DECREF(a);
	Py_DECREF(b);
	return result;
}


/*
 * Accessors
 */


static PyObject *__abs__(PyObject *self)
{
	return PyNumber_Subtract(PyTuple_GET_ITEM(self, 1), PyTuple_GET_ITEM(self, 0));
}


/*
 * Comparisons
 */


static int __nonzero__(PyObject *self)
{
	return PyObject_Compare(PyTuple_GET_ITEM(self, 0), PyTuple_GET_ITEM(self, 1)) != 0;
}


static PyObject *intersects(PyObject *self, PyObject *other)
{
	PyObject *sa = PyTuple_GET_ITEM(self, 0);
	PyObject *sb = PyTuple_GET_ITEM(self, 1);
	PyObject *oa = PyTuple_GET_ITEM(other, 0);
	PyObject *ob = PyTuple_GET_ITEM(other, 1);
	PyObject *result = (PyObject_Compare(sb, oa) > 0) && (PyObject_Compare(sa, ob) < 0) ? Py_True : Py_False;
	Py_INCREF(result);
	return result;
}


static int __contains__(PyObject *self, PyObject *other)
{
	PyObject *sa = PyTuple_GET_ITEM(self, 0);
	PyObject *sb = PyTuple_GET_ITEM(self, 1);
	if(segments_Segment_Check(other)) {
		PyObject *oa = PyTuple_GET_ITEM(other, 0);
		PyObject *ob = PyTuple_GET_ITEM(other, 1);
		return (PyObject_Compare(sa, oa) <= 0) && (PyObject_Compare(sb, ob) >= 0);
	} else
		return (PyObject_Compare(sa, other) <= 0) && (PyObject_Compare(other, sb) < 0);
}


static PyObject *not_continuous(PyObject *self, PyObject *other)
{
	PyObject *sa = PyTuple_GET_ITEM(self, 0);
	PyObject *sb = PyTuple_GET_ITEM(self, 1);
	PyObject *oa = PyTuple_GET_ITEM(other, 0);
	PyObject *ob = PyTuple_GET_ITEM(other, 1);
	PyObject *result;
	if(PyObject_Compare(sa, ob) > 0)
		return PyInt_FromLong(1);
	if(PyObject_Compare(sb, oa) < 0)
		return PyInt_FromLong(-1);
	return PyInt_FromLong(0);
}


static PyObject *continuous(PyObject *self, PyObject *other)
{
	PyObject *sa = PyTuple_GET_ITEM(self, 0);
	PyObject *sb = PyTuple_GET_ITEM(self, 1);
	PyObject *oa = PyTuple_GET_ITEM(other, 0);
	PyObject *ob = PyTuple_GET_ITEM(other, 1);
	PyObject *result = (PyObject_Compare(sb, oa) >= 0) && (PyObject_Compare(sa, ob) <= 0) ? Py_True : Py_False;
	Py_INCREF(result);
	return result;
}


static PyObject *order(PyObject *self, PyObject *other)
{
	if(__contains__(self, other))
		return PyInt_FromLong(0);
	return PyInt_FromLong(PyObject_Compare(self, other));
}


/*
 * Arithmetic
 */


static PyObject *__and__(PyObject *self, PyObject *other)
{
	PyObject *sa = PyTuple_GET_ITEM(self, 0);
	PyObject *sb = PyTuple_GET_ITEM(self, 1);
	PyObject *oa = PyTuple_GET_ITEM(other, 0);
	PyObject *ob = PyTuple_GET_ITEM(other, 1);
	PyObject *a, *b;

	if((PyObject_Compare(sb, oa) <= 0) || (PyObject_Compare(sa, ob) >= 0)) {
		/* self and other don't intersect */
		PyErr_SetObject(PyExc_ValueError, other);
		return NULL;
	}
	a = (PyObject_Compare(sa, oa) >= 0) ? sa : oa;
	b = (PyObject_Compare(sb, ob) <= 0) ? sb : ob;
	if((a == sa) && (b == sb)) {
		/* re-use self */
		Py_INCREF(self);
		return self;
	}
	if((a == oa) && (b == ob)) {
		/* re-use other */
		Py_INCREF(other);
		return other;
	}
	Py_INCREF(a);
	Py_INCREF(b);
	return segments_Segment_New(&segments_Segment_Type, a, b);
}


static PyObject *__or__(PyObject *self, PyObject *other)
{
	PyObject *sa = PyTuple_GET_ITEM(self, 0);
	PyObject *sb = PyTuple_GET_ITEM(self, 1);
	PyObject *oa = PyTuple_GET_ITEM(other, 0);
	PyObject *ob = PyTuple_GET_ITEM(other, 1);
	PyObject *a, *b;

	if((PyObject_Compare(sb, oa) < 0) || (PyObject_Compare(sa, ob) > 0)) {
		/* self and other are disjoint */
		PyErr_SetObject(PyExc_ValueError, other);
		return NULL;
	}
	a = (PyObject_Compare(sa, oa) <= 0) ? sa : oa;
	b = (PyObject_Compare(sb, ob) >= 0) ? sb : ob;
	if((a == sa) && (b == sb)) {
		/* re-use self */
		Py_INCREF(self);
		return self;
	}
	if((a == oa) && (b == ob)) {
		/* re-use other */
		Py_INCREF(other);
		return other;
	}
	Py_INCREF(a);
	Py_INCREF(b);
	return segments_Segment_New(&segments_Segment_Type, a, b);
}


static PyObject *__sub__(PyObject *self, PyObject *other)
{
	PyObject *sa = PyTuple_GET_ITEM(self, 0);
	PyObject *sb = PyTuple_GET_ITEM(self, 1);
	PyObject *oa = PyTuple_GET_ITEM(other, 0);
	PyObject *ob = PyTuple_GET_ITEM(other, 1);
	PyObject *a, *b;

	if((PyObject_Compare(sb, oa) <= 0) || (PyObject_Compare(sa, ob) >= 0)) {
		/* self and other do not intersect */
		Py_INCREF(self);
		return self;
	}
	if(__contains__(other, self) || ((PyObject_Compare(sa, oa) < 0) && (PyObject_Compare(sb, ob) > 0))) {
		/* result is not exactly 1 segment */
		PyErr_SetObject(PyExc_ValueError, other);
		return NULL;
	}
	if(PyObject_Compare(sa, oa) < 0) {
		a = sa;
		b = oa;
	} else {
		a = ob;
		b = sb;
	}
	Py_INCREF(a);
	Py_INCREF(b);
	return segments_Segment_New(&segments_Segment_Type, a, b);
}


/*
 * Protraction and contraction and shifting
 */


static PyObject *protract(PyObject *self, PyObject *delta)
{
	PyObject *a = PyNumber_Subtract(PyTuple_GET_ITEM(self, 0), delta);
	PyObject *b = PyNumber_Add(PyTuple_GET_ITEM(self, 1), delta);
	if(PyErr_Occurred()) {
		Py_DECREF(a);
		Py_DECREF(b);
		return NULL;
	}
	if(PyObject_Compare(a, b) <= 0)
		return segments_Segment_New(&segments_Segment_Type, a, b);
	return segments_Segment_New(&segments_Segment_Type, b, a);
}


static PyObject *contract(PyObject *self, PyObject *delta)
{
	PyObject *a = PyNumber_Add(PyTuple_GET_ITEM(self, 0), delta);
	PyObject *b = PyNumber_Subtract(PyTuple_GET_ITEM(self, 1), delta);
	if(PyErr_Occurred()) {
		Py_DECREF(a);
		Py_DECREF(b);
		return NULL;
	}
	if(PyObject_Compare(a, b) <= 0)
		return segments_Segment_New(&segments_Segment_Type, a, b);
	return segments_Segment_New(&segments_Segment_Type, b, a);
}


static PyObject *shift(PyObject *self, PyObject *delta)
{
	PyObject *a = PyNumber_Add(PyTuple_GET_ITEM(self, 0), delta);
	PyObject *b = PyNumber_Add(PyTuple_GET_ITEM(self, 1), delta);
	if(PyErr_Occurred()) {
		Py_DECREF(a);
		Py_DECREF(b);
		return NULL;
	}
	return segments_Segment_New(&segments_Segment_Type, a, b);
}


/*
 * Type information
 */


static PyNumberMethods as_number = {
	.nb_add = __or__,
	.nb_and = __and__,
	.nb_absolute = __abs__,
	.nb_nonzero = __nonzero__,
	.nb_or = __or__,
	.nb_subtract = __sub__,
};


static PySequenceMethods as_sequence = {
	.sq_contains = __contains__,
};


static struct PyMethodDef methods[] = {
	{"not_continuous", not_continuous, METH_O, ""},
	{"order", order, METH_O, ""},
	{"intersects", intersects, METH_O, ""},
	{"continuous", continuous, METH_O, ""},
	{"protract", protract, METH_O, ""},
	{"contract", contract, METH_O, ""},
	{"shift", shift, METH_O, ""},
	{NULL,}
};


PyTypeObject segments_Segment_Type = {
	PyObject_HEAD_INIT(NULL)
	.tp_base = &PyTuple_Type,
	.tp_as_number = &as_number,
	.tp_as_sequence = &as_sequence,
	.tp_doc =
	"",
	.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_CHECKTYPES | Py_TPFLAGS_BASETYPE,
	.tp_methods = methods,
	.tp_name = MODULE_NAME ".segment",
	.tp_new = __new__,
	.tp_repr = __repr__,
	.tp_str = __str__,
};
