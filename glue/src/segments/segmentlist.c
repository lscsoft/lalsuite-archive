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


/*
 * Basic methods
 */


/*
 * Accessors
 */


static PyObject *duration(PyObject *self, PyObject *nul)
{
	/* FIXME */
	return NULL;
}


static PyObject *extent(PyObject *self, PyObject *nul)
{
	/* FIXME */
	return NULL;
}


static PyObject *find(PyObject *self, PyObject *item)
{
	/* FIXME */
	return NULL;
}


/*
 * Comparisons
 */


static PyObject *intersects(PyObject *self, PyObject *other)
{
	/* FIXME */
	return NULL;
}


static PyObject *intersects_segment(PyObject *self, PyObject *other)
{
	/* FIXME */
	return NULL;
}


static int __contains__(PyObject *self, PyObject *other)
{
	/* FIXME */
	return -1;
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
	/* FIXME */
	return NULL;
}


static PyObject *__ior__(PyObject *self, PyObject *other)
{
	/* FIXME */
	return NULL;
}


static PyObject *__or__(PyObject *self, PyObject *other)
{
	/* FIXME */
	return NULL;
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
	/* FIXME */
	return NULL;
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
};


static PySequenceMethods as_sequence = {
	.sq_contains = __contains__,
};


static struct PyMethodDef methods[] = {
	{"find", find, METH_O, ""},
	{"duration", duration, METH_NOARGS, ""},
	{"extent", extent, METH_NOARGS, ""},
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
	.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_CHECKTYPES,
	.tp_methods = methods,
	.tp_name = MODULE_NAME ".segmentlist",
};
