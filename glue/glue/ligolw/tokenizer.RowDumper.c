/*
 * $Id$
 *
 * Copyright (C) 2007  Kipp C. Cannon
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
 *                         tokenizer.RowDumper Class
 *
 * ============================================================================
 */


#include <Python.h>
#include <structmember.h>
#include <stdlib.h>
#include <tokenizer.h>


/*
 * ============================================================================
 *
 *                              Row Dumper Type
 *
 * ============================================================================
 */


/*
 * Structure
 */


typedef struct {
	PyObject_HEAD
	/* delimiter */
	PyObject *delimiter;
	/* tuple of attribute names */
	PyObject *attributes;
	/* tuple of format strings */
	PyObject *formats;
	/* the source of row objects to be turned to strings */
	PyObject *iter;
	/* number of rows converted so far */
	int rows_converted;
	/* tuple of unicode representations of values in most recently
	 * converted row */
	PyObject *representations;
} ligolw_RowDumper;


/*
 * __del__() method
 */


static void __del__(PyObject *self)
{
	ligolw_RowDumper *rowdumper = (ligolw_RowDumper *) self;

	Py_XDECREF(rowdumper->delimiter);
	Py_XDECREF(rowdumper->attributes);
	Py_XDECREF(rowdumper->formats);
	Py_XDECREF(rowdumper->iter);
	Py_XDECREF(rowdumper->representations);

	self->ob_type->tp_free(self);
}


/*
 * __init__() method
 */


static int __init__(PyObject *self, PyObject *args, PyObject *kwds)
{
	ligolw_RowDumper *rowdumper = (ligolw_RowDumper *) self;

	if(!PyArg_ParseTuple(args, "OOOO", &rowdumper->attributes, &rowdumper->formats, &rowdumper->delimiter, &rowdumper->iter))
		return -1;

	rowdumper->delimiter = PyUnicode_FromObject(rowdumper->delimiter);
	rowdumper->attributes = _build_attributes(rowdumper->attributes);
	rowdumper->formats = _build_formats(rowdumper->formats);
	rowdumper->iter = PyObject_GetIter(rowdumper->iter);
	if(!rowdumper->delimiter || !rowdumper->attributes || !rowdumper->formats || !rowdumper->iter) {
		Py_XDECREF(rowdumper->delimiter);
		Py_XDECREF(rowdumper->attributes);
		Py_XDECREF(rowdumper->formats);
		Py_XDECREF(rowdumper->iter);
		return -1;
	}

	if(PyTuple_GET_SIZE(rowdumper->attributes) != PyTuple_GET_SIZE(rowdumper->formats)) {
		Py_DECREF(rowdumper->delimiter);
		Py_DECREF(rowdumper->attributes);
		Py_DECREF(rowdumper->formats);
		Py_DECREF(rowdumper->iter);
		PyErr_SetString(PyExc_ValueError, "len(attributes) != len(formats)");
		return -1;
	}

	rowdumper->rows_converted = 0;
	rowdumper->representations = Py_None;
	Py_INCREF(rowdumper->representations);

	return 0;
}


/*
 * __iter__() method
 */


static PyObject *__iter__(PyObject *self)
{
	Py_INCREF(self);
	return self;
}


/*
 * next() method
 */


static PyObject *next(PyObject *self)
{
	ligolw_RowDumper *rowdumper = (ligolw_RowDumper *) self;
	PyObject *row;

	row = PyIter_Next(rowdumper->iter);
	if(row) {
		const int n = PyTuple_GET_SIZE(rowdumper->attributes);
		PyObject *result;
		int i;

		Py_DECREF(rowdumper->representations);
		rowdumper->representations = PyTuple_New(n);
		if(!rowdumper->representations) {
			rowdumper->representations = Py_None;
			Py_INCREF(rowdumper->representations);
			Py_DECREF(row);
			return NULL;
		}

		for(i = 0; i < n; i++) {
			PyObject *val = PyObject_GetAttr(row, PyTuple_GET_ITEM(rowdumper->attributes, i));
			PyObject *r;

			if(!val) {
				Py_DECREF(rowdumper->representations);
				rowdumper->representations = Py_None;
				Py_INCREF(rowdumper->representations);
				Py_DECREF(row);
				return NULL;
			}

			/* the commented out bits would enable support for
			 * writing None values as empty entries in the
			 * table. */
#if 0
			if(val == Py_None)
				r = PyUnicode_FromUnicode(NULL, 0); /* u"" */
			else
#endif
				r = PyNumber_Remainder(PyTuple_GET_ITEM(rowdumper->formats, i), val);
			Py_DECREF(val);

			if(!r) {
				Py_DECREF(rowdumper->representations);
				rowdumper->representations = Py_None;
				Py_INCREF(rowdumper->representations);
				Py_DECREF(row);
				return NULL;
			}

			PyTuple_SET_ITEM(rowdumper->representations, i, r);
		}

		result = PyUnicode_Join(rowdumper->delimiter, rowdumper->representations);

		Py_DECREF(row);

		rowdumper->rows_converted += result != NULL;

		return result;
	}

	if(!PyErr_Occurred())
		PyErr_SetNone(PyExc_StopIteration);

	return NULL;
}


/*
 * Type information
 */


static struct PyMemberDef members[] = {
	{"delimiter", T_OBJECT, offsetof(ligolw_RowDumper, delimiter), READONLY, "The delimiter character."},
	{"attributes", T_OBJECT, offsetof(ligolw_RowDumper, attributes), READONLY, "In-order tuple of attribute names."},
	{"formats", T_OBJECT, offsetof(ligolw_RowDumper, formats), READONLY, "In-order tuple of format strings."},
	{"rows_converted", T_INT, offsetof(ligolw_RowDumper, rows_converted), READONLY, "Number of rows converted."},
	{"representations", T_OBJECT, offsetof(ligolw_RowDumper, representations), READONLY, "In-order tuple of unicode representations of values in most recently converted row."},
	{NULL,}
};


PyTypeObject ligolw_RowDumper_Type = {
	PyObject_HEAD_INIT(NULL)
	.tp_basicsize = sizeof(ligolw_RowDumper),
	.tp_dealloc = __del__,
	.tp_doc =
"An iterator for converting row objects into string tokens.\n" \
"\n" \
"Example:\n" \
"\n" \
">>> class Row(object):\n" \
"...     pass\n" \
"... \n" \
">>> rows = [Row(), Row(), Row()]\n" \
">>> rows[0].snr = 10.1\n" \
">>> rows[1].snr = 15.2\n" \
">>> rows[2].snr = 20.3\n" \
">>> rows[0].status = \"bad\"\n" \
">>> rows[1].status = \"bad\"\n" \
">>> rows[2].status = \"good\"\n" \
">>> rowdumper = RowDumper((\"snr\", \"status\"), (\"%.16g\", \"\\\"%s\\\"\"), \",\", rows)\n" \
">>> print \",\\n\".join(rowdumper)\n" \
"10.1,\"bad\",\n" \
"15.2,\"bad\",\n" \
"20.3,\"good\"\n" \
"\n" \
"An instance of RowDumper is initialized with four arguments.  The first\n" \
"argument is a sequence of attribute names.  The second argument is a\n" \
"sequence of Python format strings.  The third argument is a delimiter\n" \
"character.  And the final argument is an iterable object that provides row\n" \
"objects one-by-one.",
	.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
	.tp_init = __init__,
	.tp_iter = __iter__,
	.tp_iternext = next,
	.tp_members = members,
	.tp_name = MODULE_NAME ".RowDumper",
	.tp_new = PyType_GenericNew,
};
