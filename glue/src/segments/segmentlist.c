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
 *              Segments Module Component --- segmentlist Class
 *
 * ============================================================================
 */


#include <Python.h>
#include <stdlib.h>


#include <segments.h>


/*
 * ============================================================================
 *
 *                             segmentlist Class
 *
 * ============================================================================
 */


/*
 * Utilities
 */


static int segments_SegmentList_Check(PyObject *obj)
{
	return obj ? PyObject_TypeCheck(obj, &segments_SegmentList_Type) : 0;
}


static int bisect_left(PyObject *seglist, PyObject *seg, int lo, int hi)
{
	if(lo < 0)
		lo = 0;
	if(hi < 0) {
		hi = PyList_GET_SIZE(seglist);
		if(hi < 0)
			return -1;
	}
	while(lo < hi) {
		int mid = (lo + hi) / 2;
		PyObject *item = PyList_GET_ITEM(seglist, mid);
		int result;
		if(!item)
			return -1;
		Py_INCREF(item);
		result = PyObject_RichCompareBool(item, seg, Py_LT);
		Py_DECREF(item);
		if(result > 0)
			/* item < seg */
			lo = mid + 1;
		else if(result == 0)
			/* item >= seg */
			hi = mid;
		else
			/* error */
			return -1;
	}

	return lo;
}


static PyListObject *segments_SegmentList_New(PyTypeObject *type, PyObject *sequence)
{
	PyListObject *new;

	new = (PyListObject *) type->tp_alloc(type, 0);
	if(!new)
		return NULL;

	if(sequence)
		if(!_PyList_Extend(new, sequence)) {
			Py_DECREF(new);
			return NULL;
		}

	return new;
}


/*
 * Accessors
 */


static PyObject *__abs__(PyObject *self)
{
	int n = PyList_GET_SIZE(self);
	int i;
	PyObject *abs;

	if(n < 0)
		return NULL;

	abs = PyInt_FromLong(0);
	if(!abs)
		return NULL;

	for(i = 0; i < n; i++) {
		PyObject *seg, *segsize, *newabs;
		seg = PyList_GET_ITEM(self, i);
		Py_INCREF(seg);
		segsize = PyNumber_Absolute(seg);
		Py_DECREF(seg);
		if(!segsize) {
			Py_DECREF(abs);
			return NULL;
		}
		newabs = PyNumber_InPlaceAdd(abs, segsize);
		Py_DECREF(segsize);
		Py_DECREF(abs);
		abs = newabs;
		if(!abs)
			return NULL;
	}

	return abs;
}


static PyObject *extent(PyObject *self, PyObject *nul)
{
	int n = PyList_GET_SIZE(self);
	int i;
	PyObject *seg, *min, *max;

	if(n < 0)
		return NULL;
	if(n < 1) {
		PyErr_SetString(PyExc_ValueError, "empty list");
		return NULL;
	}

	seg = PyList_GET_ITEM(self, 0);
	min = PyTuple_GetItem(seg, 0);
	max = PyTuple_GetItem(seg, 1);
	if(!(min && max))
		return NULL;
	Py_INCREF(min);
	Py_INCREF(max);

	for(i = 1; i < n; i++) {
		PyObject *item_min, *item_max;
		int result;

		seg = PyList_GET_ITEM(self, i);
		item_min = PyTuple_GetItem(seg, 0);
		item_max = PyTuple_GetItem(seg, 1);
		if(!(item_min && item_max)) {
			Py_DECREF(min);
			Py_DECREF(max);
			return NULL;
		}
		Py_INCREF(item_min);
		Py_INCREF(item_max);

		result = PyObject_RichCompareBool(min, item_min, Py_GT);
		if(result > 0) {
			Py_DECREF(min);
			min = item_min;
		} else if(result < 0) {
			Py_DECREF(min);
			Py_DECREF(max);
			Py_DECREF(item_min);
			Py_DECREF(item_max);
			return NULL;
		} else
			Py_DECREF(item_min);

		result = PyObject_RichCompareBool(max, item_max, Py_LT);
		if(result > 0) {
			Py_DECREF(max);
			max = item_max;
		} else if(result < 0) {
			Py_DECREF(min);
			Py_DECREF(max);
			Py_DECREF(item_min);
			Py_DECREF(item_max);
			return NULL;
		} else
			Py_DECREF(item_max);
	}

	/* this consumes the references to min and max */
	return segments_Segment_New(&segments_Segment_Type, min, max);
}


static PyObject *find(PyObject *self, PyObject *item)
{
	int n = PyList_GET_SIZE(self);
	int i;

	if(n < 0)
		return NULL;
	Py_INCREF(item);
	for(i = 0; i < n; i++) {
		int result;
		PyObject *seg = PyList_GET_ITEM(self, i);
		Py_INCREF(seg);
		result = PySequence_Contains(seg, item);
		Py_DECREF(seg);
		if(!result)
			/* not a match */
			continue;
		Py_DECREF(item);
		if(result > 0)
			/* match found */
			return PyInt_FromLong(i);
		/* result < 0 == error */
		return NULL;
	}

	PyErr_SetObject(PyExc_ValueError, item);
	return NULL;
}


/*
 * Comparisons
 */


static PyObject *intersects(PyObject *self, PyObject *other)
{
	int n_self = PyList_GET_SIZE(self);
	int n_other = PyList_Size(other);
	PyObject *seg;
	PyObject *seg_lo;
	PyObject *seg_hi;
	PyObject *oth_lo;
	PyObject *oth_hi;
	int i = 0;
	int j = 0;

	if((n_self < 0) || (n_other < 0))
		return NULL;
	if((n_self < 1) || (n_other < 1)) {
		Py_INCREF(Py_False);
		return Py_False;
	}

	seg = PyList_GET_ITEM(self, 0);
	seg_lo = PyTuple_GetItem(seg, 0);
	seg_hi = PyTuple_GetItem(seg, 1);
	seg = PyList_GetItem(other, 0);
	oth_lo = PyTuple_GetItem(seg, 0);
	oth_hi = PyTuple_GetItem(seg, 1);
	if(!(seg_lo && seg_hi && oth_lo && oth_hi))
		return NULL;
	Py_INCREF(seg_lo);
	Py_INCREF(seg_hi);
	Py_INCREF(oth_lo);
	Py_INCREF(oth_hi);
	while(1) {
		int result;

		result = PyObject_RichCompareBool(seg_hi, oth_lo, Py_LE);
		if(result < 0) {
			Py_DECREF(seg_lo);
			Py_DECREF(seg_hi);
			Py_DECREF(oth_lo);
			Py_DECREF(oth_hi);
			return NULL;
		}
		if(result > 1) {
			Py_DECREF(seg_lo);
			Py_DECREF(seg_hi);
			if(++i >= n_self) {
				Py_DECREF(oth_lo);
				Py_DECREF(oth_hi);
				Py_INCREF(Py_False);
				return Py_False;
			}
			seg = PyList_GET_ITEM(self, i);
			seg_lo = PyTuple_GetItem(seg, 0);
			seg_hi = PyTuple_GetItem(seg, 1);
			if(!(seg_lo && seg_hi))
				return NULL;
			Py_INCREF(seg_lo);
			Py_INCREF(seg_hi);
			continue;
		}

		result = PyObject_RichCompareBool(oth_hi, seg_lo, Py_LE);
		if(result < 0) {
			Py_DECREF(seg_lo);
			Py_DECREF(seg_hi);
			Py_DECREF(oth_lo);
			Py_DECREF(oth_hi);
			return NULL;
		}
		if(result > 1) {
			Py_DECREF(oth_lo);
			Py_DECREF(oth_hi);
			if(++j >= n_other) {
				Py_DECREF(seg_lo);
				Py_DECREF(seg_hi);
				Py_INCREF(Py_False);
				return Py_False;
			}
			seg = PyList_GetItem(other, j);
			oth_lo = PyTuple_GetItem(seg, 0);
			oth_hi = PyTuple_GetItem(seg, 1);
			if(!(oth_lo && oth_hi))
				return NULL;
			Py_INCREF(oth_lo);
			Py_INCREF(oth_hi);
			continue;
		}

		/* self[i] and other[j] intersect */
		Py_INCREF(Py_True);
		return Py_True;
	}
}


static PyObject *intersects_segment(PyObject *self, PyObject *other)
{
	int i = bisect_left(self, other, -1, -1);
	PyObject *a, *b;
	int result;

	if(i < 0)
		/* error */
		return NULL;

	if(i != 0) {
		a = PyTuple_GetItem(other, 0);
		b = PyTuple_GetItem(PyList_GET_ITEM(self, i - 1), 1);
		Py_INCREF(a);
		Py_INCREF(b);
		result = PyObject_RichCompareBool(a, b, Py_LT);
		Py_DECREF(a);
		Py_DECREF(b);
		if(result > 0) {
			Py_INCREF(Py_True);
			return Py_True;
		}
		if(result < 0)
			/* error */
			return NULL;
	}

	if(i != PyList_GET_SIZE(self)) {
		a = PyTuple_GetItem(other, 1);
		b = PyTuple_GetItem(PyList_GET_ITEM(self, i), 0);
		Py_INCREF(a);
		Py_INCREF(b);
		result = PyObject_RichCompareBool(a, b, Py_GT);
		Py_DECREF(a);
		Py_DECREF(b);
		if(result > 0) {
			Py_INCREF(Py_True);
			return Py_True;
		}
		if(result < 0)
			/* error */
			return NULL;
	}

	Py_INCREF(Py_False);
	return Py_False;
}


static int __contains__(PyObject *self, PyObject *other)
{
	int n = PyList_GET_SIZE(self);
	int i;

	for(i = 0; i < n; i++) {
		PyObject *seg = PyList_GET_ITEM(self, i);
		int result;
		Py_INCREF(seg);
		result = PySequence_Contains(seg, other);
		Py_DECREF(seg);
		if(!result)
			continue;
		return result > 0 ? 1 : result;
	}

	return 0;
}


/*
 * Arithmetic
 */


static PyObject *__iand__(PyObject *self, PyObject *other)
{
	/* FIXME */
	return NULL;
}


static PyObject *__and__(PyObject *self, PyObject *other)
{
	PyObject *result;

	self = (PyObject *) segments_SegmentList_New(&segments_SegmentList_Type, self);
	if(!self)
		return NULL;
	result = PyNumber_InPlaceAnd(self, other);
	Py_DECREF(self);
	return result;
}


static PyObject *__ior__(PyObject *self, PyObject *other)
{
	/* FIXME */
	return NULL;
}


static PyObject *__or__(PyObject *self, PyObject *other)
{
	PyObject *result;

	self = (PyObject *) segments_SegmentList_New(&segments_SegmentList_Type, self);
	if(!self)
		return NULL;
	result = PyNumber_InPlaceOr(self, other);
	Py_DECREF(self);
	return result;
}


static PyObject *__xor__(PyObject *self, PyObject *other)
{
	/* FIXME */
	return NULL;
}


static PyObject *__isub__(PyObject *self, PyObject *other)
{
	/* FIXME */
	return NULL;
}


static PyObject *__sub__(PyObject *self, PyObject *other)
{
	PyObject *result;

	self = (PyObject *) segments_SegmentList_New(&segments_SegmentList_Type, self);
	if(!self)
		return NULL;
	result = PyNumber_InPlaceSubtract(self, other);
	Py_DECREF(self);
	return result;
}


static PyObject *__invert__(PyObject *self)
{
	/* FIXME */
	return NULL;
}


static PyObject *coalesce(PyObject *self, PyObject *nul)
{
	/* FIXME */
	return NULL;
}


/*
 * Protraction and contraction and shifting
 */


static PyObject *protract(PyObject *self, PyObject *delta)
{
	/* FIXME */
	return NULL;
}


static PyObject *contract(PyObject *self, PyObject *delta)
{
	/* FIXME */
	return NULL;
}


static PyObject *shift(PyObject *self, PyObject *delta)
{
	/* FIXME */
	return NULL;
}


/*
 * Type information
 */


static PyNumberMethods as_number = {
	.nb_inplace_and = __iand__,
	.nb_and = __and__,
	.nb_inplace_or = __ior__,
	.nb_or = __or__,
	.nb_xor = __xor__,
	.nb_inplace_add = __ior__,
	.nb_add = __or__,
	.nb_inplace_subtract = __isub__,
	.nb_subtract = __sub__,
	.nb_invert = __invert__,
	.nb_absolute = __abs__,
};


static PySequenceMethods as_sequence = {
	.sq_contains = __contains__,
};


static struct PyMethodDef methods[] = {
	{"extent", extent, METH_NOARGS, ""},
	{"find", find, METH_O, ""},
	{"intersects", intersects, METH_O, ""},
	{"intersects_segment", intersects_segment, METH_O, ""},
	{"coalesce", coalesce, METH_NOARGS, ""},
	{"protract", protract, METH_O, ""},
	{"contract", contract, METH_O, ""},
	{"shift", shift, METH_O, ""},
	{NULL,}
};


PyTypeObject segments_SegmentList_Type = {
	PyObject_HEAD_INIT(NULL)
	.tp_base = &PyList_Type,
	.tp_as_number = &as_number,
	.tp_as_sequence = &as_sequence,
	.tp_doc =
	"",
	.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_CHECKTYPES | Py_TPFLAGS_BASETYPE,
	.tp_methods = methods,
	.tp_name = MODULE_NAME ".segmentlist",
};
