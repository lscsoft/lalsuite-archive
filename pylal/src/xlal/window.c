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
 *                  Python Wrapper For LAL's Window Package
 *
 * ============================================================================
 */


#include <Python.h>
#include <numpy/arrayobject.h>
#include <lal/LALDatatypes.h>
#include <lal/Window.h>
#include <misc.h>
#include <datatypes/real8window.h>


#define MODULE_NAME "pylal.xlal.window"


/*
 * ============================================================================
 *
 *                                 Utilities
 *
 * ============================================================================
 */


/*
 * Add input error checking to pylal_REAL8Window_new()
 */


static PyObject *new(REAL8Window *w, PyObject *owner)
{
	if(!w) {
		pylal_set_exception_from_xlalerrno();
		return NULL;
	}

	return pylal_REAL8Window_new(w, owner);
}


/*
 * ============================================================================
 *
 *                                 Functions
 *
 * ============================================================================
 */


static PyObject *pylal_XLALUnitaryWindowREAL8Sequence(PyObject *self, PyObject *args)
{
	PyArrayObject *array;
	pylal_REAL8Window *window;
	REAL8Sequence sequence = {
		.length = 0,
		.data = NULL
	};

	if(!PyArg_ParseTuple(args, "O!O!", &PyArray_Type, &array, pylal_REAL8Window_Type, &window))
		return NULL;
	/* requier double precision floats */
	if(PyArray_TYPE(array) != NPY_DOUBLE) {
		PyErr_SetObject(PyExc_TypeError, (PyObject *) array);
		return NULL;
	}
	/* require exactly 1 dimension */
	if(array->nd != 1) {
		PyErr_SetObject(PyExc_ValueError, (PyObject *) array);
		return NULL;
	}

	sequence.length = PyArray_DIM(array, 0);
	sequence.data = PyArray_GETPTR1(array, 0);
	if(!XLALUnitaryWindowREAL8Sequence(&sequence, window->window)) {
		pylal_set_exception_from_xlalerrno();
		return NULL;
	}

	Py_INCREF(Py_None);
	return Py_None;
}


static PyObject *pylal_XLALCreateRectangularREAL8Window(PyObject *self, PyObject *args)
{
	int length;

	if(!PyArg_ParseTuple(args, "i", &length))
		return NULL;

	return new(XLALCreateRectangularREAL8Window(length), NULL);
}


static PyObject *pylal_XLALCreateHannREAL8Window(PyObject *self, PyObject *args)
{
	int length;

	if(!PyArg_ParseTuple(args, "i", &length))
		return NULL;

	return new(XLALCreateHannREAL8Window(length), NULL);
}


static PyObject *pylal_XLALCreateWelchREAL8Window(PyObject *self, PyObject *args)
{
	int length;

	if(!PyArg_ParseTuple(args, "i", &length))
		return NULL;

	return new(XLALCreateWelchREAL8Window(length), NULL);
}


static PyObject *pylal_XLALCreateBartlettREAL8Window(PyObject *self, PyObject *args)
{
	int length;

	if(!PyArg_ParseTuple(args, "i", &length))
		return NULL;

	return new(XLALCreateBartlettREAL8Window(length), NULL);
}


static PyObject *pylal_XLALCreateParzenREAL8Window(PyObject *self, PyObject *args)
{
	int length;

	if(!PyArg_ParseTuple(args, "i", &length))
		return NULL;

	return new(XLALCreateParzenREAL8Window(length), NULL);
}


static PyObject *pylal_XLALCreatePapoulisREAL8Window(PyObject *self, PyObject *args)
{
	int length;

	if(!PyArg_ParseTuple(args, "i", &length))
		return NULL;

	return new(XLALCreatePapoulisREAL8Window(length), NULL);
}


static PyObject *pylal_XLALCreateHammingREAL8Window(PyObject *self, PyObject *args)
{
	int length;

	if(!PyArg_ParseTuple(args, "i", &length))
		return NULL;

	return new(XLALCreateHammingREAL8Window(length), NULL);
}


static PyObject *pylal_XLALCreateKaiserREAL8Window(PyObject *self, PyObject *args)
{
	int length;
	double beta;

	if(!PyArg_ParseTuple(args, "id", &length, &beta))
		return NULL;

	return new(XLALCreateKaiserREAL8Window(length, beta), NULL);
}


static PyObject *pylal_XLALCreateCreightonREAL8Window(PyObject *self, PyObject *args)
{
	int length;
	double beta;

	if(!PyArg_ParseTuple(args, "id", &length, &beta))
		return NULL;

	return new(XLALCreateCreightonREAL8Window(length, beta), NULL);
}


static PyObject *pylal_XLALCreateTukeyREAL8Window(PyObject *self, PyObject *args)
{
	int length;
	double beta;

	if(!PyArg_ParseTuple(args, "id", &length, &beta))
		return NULL;

	return new(XLALCreateTukeyREAL8Window(length, beta), NULL);
}


static PyObject *pylal_XLALCreateGaussREAL8Window(PyObject *self, PyObject *args)
{
	int length;
	double beta;

	if(!PyArg_ParseTuple(args, "id", &length, &beta))
		return NULL;

	return new(XLALCreateGaussREAL8Window(length, beta), NULL);
}


/*
 * ============================================================================
 *
 *                            Module Registration
 *
 * ============================================================================
 */


static struct PyMethodDef methods[] = {
	{"XLALUnitaryWindowREAL8Sequence", pylal_XLALUnitaryWindowREAL8Sequence, METH_VARARGS, NULL},
	{"XLALCreateRectangularREAL8Window", pylal_XLALCreateRectangularREAL8Window, METH_VARARGS, NULL},
	{"XLALCreateHannREAL8Window", pylal_XLALCreateHannREAL8Window, METH_VARARGS, NULL},
	{"XLALCreateWelchREAL8Window", pylal_XLALCreateWelchREAL8Window, METH_VARARGS, NULL},
	{"XLALCreateBartlettREAL8Window", pylal_XLALCreateBartlettREAL8Window, METH_VARARGS, NULL},
	{"XLALCreateParzenREAL8Window", pylal_XLALCreateParzenREAL8Window, METH_VARARGS, NULL},
	{"XLALCreatePapoulisREAL8Window", pylal_XLALCreatePapoulisREAL8Window, METH_VARARGS, NULL},
	{"XLALCreateHammingREAL8Window", pylal_XLALCreateHammingREAL8Window, METH_VARARGS, NULL},
	{"XLALCreateKaiserREAL8Window", pylal_XLALCreateKaiserREAL8Window, METH_VARARGS, NULL},
	{"XLALCreateCreightonREAL8Window", pylal_XLALCreateCreightonREAL8Window, METH_VARARGS, NULL},
	{"XLALCreateTukeyREAL8Window", pylal_XLALCreateTukeyREAL8Window, METH_VARARGS, NULL},
	{"XLALCreateGaussREAL8Window", pylal_XLALCreateGaussREAL8Window, METH_VARARGS, NULL},
	{NULL,}
};


void initwindow(void)
{
	/* commented out to silence warning */
	/*PyObject *module = */Py_InitModule3(MODULE_NAME, methods, "Wrapper for LAL's window package.");

	import_array();
	pylal_real8window_import();
}
