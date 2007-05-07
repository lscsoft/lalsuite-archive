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


/* copied from bisect.py */

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
		if(result < 0)
			return -1;
		else if(result > 0)
			/* item < seg */
			lo = mid + 1;
		else
			/* item >= seg */
			hi = mid;
	}

	return lo;
}


/* copied from bisect.py */

static int bisect_right(PyObject *seglist, PyObject *seg, int lo, int hi)
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
		result = PyObject_RichCompareBool(seg, item, Py_LT);
		Py_DECREF(item);
		if(result < 0)
			return -1;
		else if(result > 0)
			/* seg < item */
			hi = mid;
		else
			/* seg >= item */
			lo = mid + 1;
	}

	return lo;
}


static int segs_are_disjoint(PyObject *a, PyObject *b)
{
	PyObject *lo, *hi;
	int result;

	if(!(a && b))
		return -1;

	lo = PyTuple_GetItem(a, 0);
	hi = PyTuple_GetItem(b, 1);
	if(!(lo && hi))
		return -1;
	Py_INCREF(lo);
	Py_INCREF(hi);
	result = PyObject_RichCompareBool(hi, lo, Py_LT);
	Py_DECREF(lo);
	Py_DECREF(hi);
	if(result)
		return result;

	lo = PyTuple_GetItem(b, 0);
	hi = PyTuple_GetItem(a, 1);
	if(!(lo && hi))
		return -1;
	Py_INCREF(lo);
	Py_INCREF(hi);
	result = PyObject_RichCompareBool(hi, lo, Py_LT);
	Py_DECREF(lo);
	Py_DECREF(hi);

	return result;
}


static PyObject *make_segment(PyObject *lo, PyObject *hi)
{
	PyObject *seg = segments_Segment_New(&segments_Segment_Type, lo, hi);
	if(!seg) {
		Py_DECREF(lo);
		Py_DECREF(hi);
	}
	return seg;
}


static PyListObject *segments_SegmentList_New(PyTypeObject *type, PyObject *sequence)
{
	PyListObject *new = (PyListObject *) type->tp_alloc(type, 0);
	if(new && sequence)
		if(!_PyList_Extend(new, sequence)) {
			Py_DECREF(new);
			new = NULL;
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

		if((result = PyObject_RichCompareBool(min, item_min, Py_GT)) < 0) {
			Py_DECREF(min);
			Py_DECREF(max);
			Py_DECREF(item_min);
			Py_DECREF(item_max);
			return NULL;
		} else if(result > 0) {
			Py_DECREF(min);
			min = item_min;
		} else
			Py_DECREF(item_min);

		if((result = PyObject_RichCompareBool(max, item_max, Py_LT)) < 0) {
			Py_DECREF(min);
			Py_DECREF(max);
			Py_DECREF(item_min);
			Py_DECREF(item_max);
			return NULL;
		} else if(result > 0) {
			Py_DECREF(max);
			max = item_max;
		} else
			Py_DECREF(item_max);
	}

	return make_segment(min, max);
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
		if(result < 0)
			return NULL;
		else if(result > 0)
			/* match found */
			return PyInt_FromLong(i);
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
	int n_other = PySequence_Size(other);
	PyObject *seg;
	PyObject *lo;
	PyObject *hi;
	PyObject *olo;
	PyObject *ohi;
	int result;
	int i, j;

	if((n_self < 0) || (n_other < 0))
		return NULL;
	if((n_self < 1) || (n_other < 1)) {
		Py_INCREF(Py_False);
		return Py_False;
	}

	i = j = 0;

	seg = PyList_GET_ITEM(self, 0);
	lo = PyTuple_GetItem(seg, 0);
	hi = PyTuple_GetItem(seg, 1);
	seg = PySequence_GetItem(other, 0);
	olo = PyTuple_GetItem(seg, 0);
	ohi = PyTuple_GetItem(seg, 1);
	if(!(lo && hi && olo && ohi)) {
		Py_DECREF(seg);
		return NULL;
	}
	Py_INCREF(lo);
	Py_INCREF(hi);
	Py_INCREF(olo);
	Py_INCREF(ohi);
	Py_DECREF(seg);

	while(1) {
		if((result = PyObject_RichCompareBool(hi, olo, Py_LE)) < 0) {
			Py_DECREF(lo);
			Py_DECREF(hi);
			Py_DECREF(olo);
			Py_DECREF(ohi);
			return NULL;
		} else if(result > 0) {
			Py_DECREF(lo);
			Py_DECREF(hi);
			if(++i >= n_self) {
				Py_DECREF(olo);
				Py_DECREF(ohi);
				Py_INCREF(Py_False);
				return Py_False;
			}
			seg = PyList_GET_ITEM(self, i);
			lo = PyTuple_GetItem(seg, 0);
			hi = PyTuple_GetItem(seg, 1);
			if(!(lo && hi)) {
				Py_DECREF(olo);
				Py_DECREF(ohi);
				return NULL;
			}
			Py_INCREF(lo);
			Py_INCREF(hi);
		} else if((result = PyObject_RichCompareBool(ohi, lo, Py_LE)) < 0) {
			Py_DECREF(lo);
			Py_DECREF(hi);
			Py_DECREF(olo);
			Py_DECREF(ohi);
			return NULL;
		} else if(result > 0) {
			Py_DECREF(olo);
			Py_DECREF(ohi);
			if(++j >= n_other) {
				Py_DECREF(lo);
				Py_DECREF(hi);
				Py_INCREF(Py_False);
				return Py_False;
			}
			seg = PySequence_GetItem(other, j);
			olo = PyTuple_GetItem(seg, 0);
			ohi = PyTuple_GetItem(seg, 1);
			if(!(olo && ohi)) {
				Py_DECREF(lo);
				Py_DECREF(hi);
				Py_DECREF(seg);
				return NULL;
			}
			Py_INCREF(olo);
			Py_INCREF(ohi);
			Py_DECREF(seg);
		} else {
			/* self[i] and other[j] intersect */
			Py_DECREF(lo);
			Py_DECREF(hi);
			Py_DECREF(olo);
			Py_DECREF(ohi);
			Py_INCREF(Py_True);
			return Py_True;
		}
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
		if(!(a && b))
			return NULL;
		Py_INCREF(a);
		Py_INCREF(b);
		result = PyObject_RichCompareBool(a, b, Py_LT);
		Py_DECREF(a);
		Py_DECREF(b);
		if(result < 0)
			return NULL;
		else if(result > 0) {
			Py_INCREF(Py_True);
			return Py_True;
		}
	}

	if(i != PyList_GET_SIZE(self)) {
		a = PyTuple_GetItem(other, 1);
		b = PyTuple_GetItem(PyList_GET_ITEM(self, i), 0);
		if(!(a && b))
			return NULL;
		Py_INCREF(a);
		Py_INCREF(b);
		result = PyObject_RichCompareBool(a, b, Py_GT);
		Py_DECREF(a);
		Py_DECREF(b);
		if(result < 0)
			return NULL;
		else if(result > 0) {
			Py_INCREF(Py_True);
			return Py_True;
		}
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
		if(result)
			return result > 0 ? 1 : result;
	}

	return 0;
}


/*
 * Arithmetic
 */


static PyObject *__iand__(PyObject *self, PyObject *other)
{
	PyObject *new = NULL;
	other = PyNumber_Invert(other);
	if(other) {
		new = PyNumber_InPlaceSubtract(self, other);
		Py_DECREF(other);
	}
	return new;
}


static PyObject *__and__(PyObject *self, PyObject *other)
{
	PyObject *new = NULL;
	self = (PyObject *) segments_SegmentList_New(&segments_SegmentList_Type, self);
	if(self) {
		new = PyNumber_InPlaceAnd(self, other);
		Py_DECREF(self);
	}
	return new;
}


static PyObject *__ior__(PyObject *self, PyObject *other)
{
	PyObject *seg;
	int result;
	int i, j, n;

	i = 0;
	other = PyObject_GetIter(other);
	if(!other)
		return NULL;
	while((seg = PyIter_Next(other))) {
		i = j = bisect_right(self, seg, i, -1);
		if(i < 0) {
			Py_DECREF(seg);
			Py_DECREF(other);
			return NULL;
		} else if(i > 0) {
			if((result = segs_are_disjoint(PyList_GET_ITEM(self, i - 1), seg)) < 0) {
				Py_DECREF(seg);
				Py_DECREF(other);
				return NULL;
			} else if(result == 0) {
				PyObject *new = PyNumber_Or(seg, PyList_GET_ITEM(self, --i));
				Py_DECREF(seg);
				seg = new;
				if(!seg) {
					Py_DECREF(other);
					return NULL;
				}
			}
		}
		n = PyList_GET_SIZE(self);
		result = 0;
		while((j < n) && !(result = segs_are_disjoint(seg, PyList_GET_ITEM(self, j)))) {
			j++;
		}
		if(result < 0) {
			Py_DECREF(seg);
			Py_DECREF(other);
			return NULL;
		}
		if(j > i) {
			PyObject *new = PyNumber_Or(seg, PyList_GET_ITEM(self, j - 1));
			Py_DECREF(seg);
			seg = new;
			if(!seg) {
				Py_DECREF(other);
				return NULL;
			}
			if(PyList_SetSlice(self, i + 1, j, NULL) < 0) {
				Py_DECREF(seg);
				Py_DECREF(other);
				return NULL;
			}
			/* _SetItem consumes a ref count */
			if(PyList_SetItem(self, i, seg) < 0) {
				Py_DECREF(seg);
				Py_DECREF(other);
				return NULL;
			}
		} else {
			/* _Insert increments seg's ref count */
			if(PyList_Insert(self, i, seg) < 0) {
				Py_DECREF(seg);
				Py_DECREF(other);
				return NULL;
			}
			Py_DECREF(seg);
		}
		i++;
	}
	Py_DECREF(other);
	if(PyErr_Occurred())
		return NULL;

	Py_INCREF(self);
	return self;
}


static PyObject *__or__(PyObject *self, PyObject *other)
{
	PyObject *new = NULL;
	int nself = PyList_GET_SIZE(self);
	int nother = PySequence_Size(other);

	if(nself < 0 || nother < 0)
		return NULL;

	if(nself >= nother) {
		self = (PyObject *) segments_SegmentList_New(&segments_SegmentList_Type, self);
		if(self) {
			new = PyNumber_InPlaceOr(self, other);
			Py_DECREF(self);
		}
	} else {
		other = (PyObject *) segments_SegmentList_New(&segments_SegmentList_Type, other);
		if(other) {
			new = PyNumber_InPlaceOr(other, self);
			Py_DECREF(other);
		}
	}
	return new;
}


static PyObject *__xor__(PyObject *self, PyObject *other)
{
	PyObject *a, *b, *new;

	a = PyNumber_Subtract(self, other);
	if(!a)
		return NULL;
	b = PyNumber_Subtract(other, self);
	if(!b) {
		Py_DECREF(a);
		return NULL;
	}
	new = PyNumber_Or(a, b);
	Py_DECREF(a);
	Py_DECREF(b);

	return new;
}


static PyObject *__isub__(PyObject *self, PyObject *other)
{
	PyObject *seg;
	PyObject *olo, *ohi;
	PyObject *lo, *hi;
	int result;
	int i, j;
	int n;
	
	n = PySequence_Size(other);
	if(n < 0)
		return NULL;
	if(n < 1) {
		Py_INCREF(self);
		return self;
	}

	i = j = 0;

	seg = PySequence_GetItem(other, j);
	if(!seg)
		return NULL;
	olo = PyTuple_GetItem(seg, 0);
	ohi = PyTuple_GetItem(seg, 1);
	if(!(olo && ohi)) {
		Py_DECREF(seg);
		return NULL;
	}
	Py_INCREF(olo);
	Py_INCREF(ohi);
	Py_DECREF(seg);

	while(i < PyList_GET_SIZE(self)) {
		seg = PyList_GET_ITEM(self, i);
		if(!seg) {
			Py_DECREF(olo);
			Py_DECREF(ohi);
			return NULL;
		}
		lo = PyTuple_GetItem(seg, 0);
		hi = PyTuple_GetItem(seg, 1);
		if(!(lo && hi)) {
			Py_DECREF(olo);
			Py_DECREF(ohi);
			return NULL;
		}
		Py_INCREF(lo);
		Py_INCREF(hi);

		while((result = PyObject_RichCompareBool(ohi, lo, Py_LE))) {
			if(result < 0) {
				Py_DECREF(olo);
				Py_DECREF(ohi);
				Py_DECREF(lo);
				Py_DECREF(hi);
				return NULL;
			}
			if(++j >= n) {
				Py_DECREF(olo);
				Py_DECREF(ohi);
				Py_DECREF(lo);
				Py_DECREF(hi);
				Py_INCREF(self);
				return self;
			}
			Py_DECREF(olo);
			Py_DECREF(ohi);
			seg = PySequence_GetItem(other, j);
			if(!seg) {
				Py_DECREF(lo);
				Py_DECREF(hi);
				return NULL;
			}
			olo = PyTuple_GetItem(seg, 0);
			ohi = PyTuple_GetItem(seg, 1);
			if(!(olo && ohi)) {
				Py_DECREF(lo);
				Py_DECREF(hi);
				Py_DECREF(seg);
				return NULL;
			}
			Py_INCREF(olo);
			Py_INCREF(ohi);
			Py_DECREF(seg);
		}

		if((result = PyObject_RichCompareBool(hi, olo, Py_LE)) < 0) {
			Py_DECREF(olo);
			Py_DECREF(ohi);
			Py_DECREF(lo);
			Py_DECREF(hi);
			return NULL;
		} else if(result > 0) {
			/* seg[1] <= otherseg[0] */
			i++;
		} else if((result = PyObject_RichCompareBool(olo, lo, Py_LE)) < 0) {
			Py_DECREF(olo);
			Py_DECREF(ohi);
			Py_DECREF(lo);
			Py_DECREF(hi);
			return NULL;
		} else if(result > 0) {
			/* otherseg[0] <= seg[0] */
			if((result = PyObject_RichCompareBool(ohi, hi, Py_GE)) < 0) {
				Py_DECREF(olo);
				Py_DECREF(ohi);
				Py_DECREF(lo);
				Py_DECREF(hi);
				return NULL;
			} else if(result > 0) {
				/* otherseg[1] >= seg[1] */
				if(PySequence_DelItem(self, i) < 0) {
					Py_DECREF(olo);
					Py_DECREF(ohi);
					Py_DECREF(lo);
					Py_DECREF(hi);
					return NULL;
				}
			} else {
				/* else */
				PyObject *newseg = make_segment(ohi, hi);
				if(!newseg) {
					Py_DECREF(olo);
					Py_DECREF(lo);
					return NULL;
				}
				/* _SetItem consumes a ref count */
				if(PyList_SetItem(self, i, newseg) < 0) {
					Py_DECREF(olo);
					Py_DECREF(lo);
					Py_DECREF(newseg);
					return NULL;
				}
				/* make_segment() consumed references,
				 * which we need */
				Py_INCREF(ohi);
				Py_INCREF(hi);
			}
		} else {
			/* else */
			PyObject *newseg = make_segment(lo, olo);
			if(!newseg) {
				Py_DECREF(ohi);
				Py_DECREF(hi);
				return NULL;
			}
			/* _SetItem consumes a ref count */
			if(PyList_SetItem(self, i++, newseg) < 0) {
				Py_DECREF(ohi);
				Py_DECREF(hi);
				Py_DECREF(newseg);
				return NULL;
			}
			/* make_segment() consumed references, which we
			 * need */
			Py_INCREF(lo);
			Py_INCREF(olo);
			if((result = PyObject_RichCompareBool(ohi, hi, Py_LT)) < 0) {
				Py_DECREF(olo);
				Py_DECREF(ohi);
				Py_DECREF(lo);
				Py_DECREF(hi);
				return NULL;
			} else if(result > 0) {
				/* otherseg[1] < seg[1] */
				newseg = make_segment(ohi, hi);
				if(!newseg) {
					Py_DECREF(olo);
					Py_DECREF(lo);
					return NULL;
				}
				/* _Insert increments the ref count */
				if(PyList_Insert(self, i, newseg) < 0) {
					Py_DECREF(olo);
					Py_DECREF(lo);
					Py_DECREF(newseg);
					return NULL;
				}
				Py_DECREF(newseg);
				/* make_segment() consumed references,
				 * which we need */
				Py_INCREF(ohi);
				Py_INCREF(hi);
			}
		}
		Py_DECREF(lo);
		Py_DECREF(hi);
	}
	Py_DECREF(olo);
	Py_DECREF(ohi);

	Py_INCREF(self);
	return self;
}


static PyObject *__sub__(PyObject *self, PyObject *other)
{
	PyObject *new = NULL;
	self = (PyObject *) segments_SegmentList_New(&segments_SegmentList_Type, self);
	if(self) {
		new = PyNumber_InPlaceSubtract(self, other);
		Py_DECREF(self);
	}
	return new;
}


static PyObject *__invert__(PyObject *self)
{
	PyObject *seg, *newseg;
	PyObject *a, *last;
	PyObject *new;
	int result;
	int n;
	int i;

	n = PyList_GET_SIZE(self);
	if(n < 0)
		return NULL;

	new = (PyObject *) segments_SegmentList_New(&segments_SegmentList_Type, NULL);
	if(!new)
		return NULL;

	if(n < 1) {
		Py_INCREF(segments_NegInfinity);
		Py_INCREF(segments_PosInfinity);
		newseg = make_segment((PyObject *) segments_NegInfinity, (PyObject *) segments_PosInfinity);
		if(!newseg) {
			Py_DECREF(new);
			return NULL;
		}
		/* _Append increments newseg's ref count */
		if(PyList_Append(new, newseg) < 0) {
			Py_DECREF(newseg);
			Py_DECREF(new);
			return NULL;
		}
		Py_DECREF(newseg);
		return new;
	}

	seg = PyList_GET_ITEM(self, 0);
	if(!seg) {
		Py_DECREF(new);
		return NULL;
	}
	a = PyTuple_GetItem(seg, 0);
	if(!a) {
		Py_DECREF(new);
		return NULL;
	}
	Py_INCREF(segments_NegInfinity);
	Py_INCREF(a);
	if((result = PyObject_RichCompareBool(a, (PyObject *) segments_NegInfinity, Py_GT)) < 0) {
		Py_DECREF(segments_NegInfinity);
		Py_DECREF(a);
		Py_DECREF(new);
		return NULL;
	} else if(result > 0) {
		newseg = make_segment((PyObject *) segments_NegInfinity, a);
		if(!newseg) {
			Py_DECREF(new);
			return NULL;
		}
		/* _Append increments newseg's ref count */
		if(PyList_Append(new, newseg) < 0) {
			Py_DECREF(newseg);
			Py_DECREF(new);
			return NULL;
		}
		Py_DECREF(newseg);
	} else {
		Py_DECREF(segments_NegInfinity);
		Py_DECREF(a);
	}

	last = PyTuple_GetItem(seg, 1);
	if(!last) {
		Py_DECREF(new);
		return NULL;
	}
	for(i = 1; i < n; i++) {
		seg = PyList_GET_ITEM(self, i);
		if(!seg) {
			Py_DECREF(new);
			return NULL;
		}
		a = PyTuple_GetItem(seg, 0);
		if(!a) {
			Py_DECREF(new);
			return NULL;
		}
		Py_INCREF(last);
		Py_INCREF(a);
		newseg = make_segment(last, a);
		if(!newseg) {
			Py_DECREF(new);
			return NULL;
		}
		/* _Append increments newseg's ref count */
		if(PyList_Append(new, newseg) < 0) {
			Py_DECREF(newseg);
			Py_DECREF(new);
			return NULL;
		}
		Py_DECREF(newseg);
		seg = PyList_GET_ITEM(self, i);
		if(!seg) {
			Py_DECREF(new);
			return NULL;
		}
		last = PyTuple_GetItem(seg, 1);
		if(!last) {
			Py_DECREF(new);
			return NULL;
		}
	}

	Py_INCREF(last);
	Py_INCREF(segments_PosInfinity);
	if((result = PyObject_RichCompareBool(last, (PyObject *) segments_PosInfinity, Py_LT)) < 0) {
		Py_DECREF(last);
		Py_DECREF(segments_PosInfinity);
		Py_DECREF(new);
		return NULL;
	} else if(result > 0) {
		newseg = make_segment(last, (PyObject *) segments_PosInfinity);
		if(!newseg) {
			Py_DECREF(new);
			return NULL;
		}
		/* _Append increments newseg's ref count */
		if(PyList_Append(new, newseg) < 0) {
			Py_DECREF(newseg);
			Py_DECREF(new);
			return NULL;
		}
		Py_DECREF(newseg);
	} else {
		Py_DECREF(last);
		Py_DECREF(segments_PosInfinity);
	}

	return new;
}


static PyObject *coalesce(PyObject *self, PyObject *nul)
{
	PyObject *seg;
	int result;
	int i, j;
	int n;

	if(PyList_Sort(self) < 0)
		return NULL;

	n = PyList_GET_SIZE(self);
	if(n < 0)
		return NULL;

	i = j = 0;
	while(j < n) {
		seg = PyList_GET_ITEM(self, j++);
		if(!seg)
			return NULL;
		Py_INCREF(seg);
		result = 0;
		while((j < n) && !(result = segs_are_disjoint(seg, PyList_GET_ITEM(self, j)))) {
			PyObject *new = PyNumber_Or(seg, PyList_GET_ITEM(self, j++));
			Py_DECREF(seg);
			seg = new;
			if(!seg)
				return NULL;
		}
		if(result < 0) {
			Py_DECREF(seg);
			return NULL;
		}
		/* _SetItem consumes a ref count */
		if(PyList_SetItem(self, i, seg) < 0) {
			Py_DECREF(seg);
			return NULL;
		}
		i++;
	}
	if(PyList_SetSlice(self, i, n, NULL) < 0)
		return NULL;

	Py_INCREF(self);
	return self;
}


/*
 * Protraction and contraction and shifting
 */


static PyObject *protract(PyObject *self, PyObject *delta)
{
	PyObject *protract;
	PyObject *seg, *new;
	int i;
	int n;

	n = PyList_GET_SIZE(self);
	if(n < 0)
		return NULL;

	protract = PyString_FromString("protract");
	if(!protract)
		return NULL;

	for(i = 0; i < n; i++) {
		seg = PyList_GET_ITEM(self, i);
		if(!seg) {
			Py_DECREF(protract);
			return NULL;
		}
		new = PyObject_CallMethodObjArgs(seg, protract, delta, NULL);
		if(!new) {
			Py_DECREF(protract);
			return NULL;
		}
		/* _SetItem consumes a ref count */
		if(PyList_SetItem(self, i, new) < 0) {
			Py_DECREF(protract);
			return NULL;
		}
	}

	Py_DECREF(protract);

	return coalesce(self, NULL);
}


static PyObject *contract(PyObject *self, PyObject *delta)
{
	PyObject *contract;
	PyObject *seg, *new;
	int i;
	int n;

	n = PyList_GET_SIZE(self);
	if(n < 0)
		return NULL;

	contract = PyString_FromString("contract");
	if(!contract)
		return NULL;

	for(i = 0; i < n; i++) {
		seg = PyList_GET_ITEM(self, i);
		if(!seg) {
			Py_DECREF(contract);
			return NULL;
		}
		new = PyObject_CallMethodObjArgs(seg, contract, delta, NULL);
		if(!new) {
			Py_DECREF(contract);
			return NULL;
		}
		/* _SetItem consumes a ref count */
		if(PyList_SetItem(self, i, new) < 0) {
			Py_DECREF(contract);
			return NULL;
		}
	}

	Py_DECREF(contract);

	return coalesce(self, NULL);
}


static PyObject *shift(PyObject *self, PyObject *delta)
{
	PyObject *shift;
	PyObject *seg, *new;
	int i;
	int n;

	n = PyList_GET_SIZE(self);
	if(n < 0)
		return NULL;

	shift = PyString_FromString("shift");
	if(!shift)
		return NULL;

	for(i = 0; i < n; i++) {
		seg = PyList_GET_ITEM(self, i);
		if(!seg) {
			Py_DECREF(shift);
			return NULL;
		}
		new = PyObject_CallMethodObjArgs(seg, shift, delta, NULL);
		if(!new) {
			Py_DECREF(shift);
			return NULL;
		}
		/* _SetItem consumes a ref count */
		if(PyList_SetItem(self, i, new) < 0) {
			Py_DECREF(shift);
			return NULL;
		}
	}

	Py_DECREF(shift);

	Py_INCREF(self);
	return self;
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
