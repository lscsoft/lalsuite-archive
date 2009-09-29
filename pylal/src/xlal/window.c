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
#include <structmember.h>
#include <numpy/arrayobject.h>
#include <lal/LALDatatypes.h>
#include <lal/Window.h>


#define MODULE_NAME "pylal.xlal.window"


/*
 * ============================================================================
 *
 *                              REAL8Window Type
 *
 * ============================================================================
 */


/*
 * Forward references
 */


static PyTypeObject pylal_REAL8Window_Type;


/*
 * Structure
 *
 * FIXME:  this is a bit of a mess.  really what should happen is that
 * first there should be a wrapping of the REAL8Sequence type, then the
 * REAL8Window wrapping would be built on top of that providing access to
 * the data object as a sequence.  that would ensure that data from inside
 * a window object could be passed to a LAL function that takes a sequence
 * as an argument, just like would be allowed in C.
 */


typedef struct {
	PyObject_HEAD
	REAL8Window *window;
	PyObject *owner;
} pylal_REAL8Window;


/*
 * Methods
 */


static void pylal_REAL8Window___del__(PyObject *self)
{
	pylal_REAL8Window *obj = (pylal_REAL8Window *) self;

	if(obj->owner)
		Py_DECREF(obj->owner);
	else
		/* we are the owner */
		XLALDestroyREAL8Window(obj->window);

	self->ob_type->tp_free(self);
}


static PyObject *pylal_REAL8Window___getattr__(PyObject *self, char *name)
{
	pylal_REAL8Window *window = (pylal_REAL8Window *) self;

	if(!strcmp(name, "sumofsquares"))
		return PyFloat_FromDouble(window->window->sumofsquares);
	if(!strcmp(name, "sum"))
		return PyFloat_FromDouble(window->window->sum);
	if(!strcmp(name, "data")) {
		npy_intp dims[] = {window->window->data->length};
		PyObject *array = PyArray_SimpleNewFromData(1, dims, NPY_FLOAT64, window->window->data->data);
		if(array) {
			/* incref self to prevent data from disappearing
			 * while array is still in use, and tell numpy to
			 * decref self when the array is deallocated */
			Py_INCREF(self);
			PyArray_BASE(array) = self;
		}
		return array;
	}
	PyErr_SetString(PyExc_AttributeError, name);
	return NULL;
}


/*
 * Type
 */


static PyTypeObject pylal_REAL8Window_Type = {
	PyObject_HEAD_INIT(NULL)
	.tp_basicsize = sizeof(pylal_REAL8Window),
	.tp_dealloc = pylal_REAL8Window___del__,
	.tp_doc = "REAL8Window structure",
	.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_CHECKTYPES,
	.tp_getattr = pylal_REAL8Window___getattr__,
	.tp_name = MODULE_NAME ".REAL8Window",
	.tp_new = PyType_GenericNew
};


/*
 * ============================================================================
 *
 *                                 Functions
 *
 * ============================================================================
 */


static PyObject *pylal_REAL8Window_new(REAL8Window *w, PyObject *owner)
{
	PyObject *new;

	if(!w)
		/* FIXME:  map XLAL error codes to Python exceptions */
		return PyErr_NoMemory();

	new = PyType_GenericNew(&pylal_REAL8Window_Type, NULL, NULL);
	if(!new) {
		XLALDestroyREAL8Window(w);
		return NULL;
	}

	((pylal_REAL8Window *) new)->window = w;
	if(owner)
		Py_INCREF(owner);
	((pylal_REAL8Window *) new)->owner = owner;

	return new;
}


static PyObject *pylal_XLALCreateRectangularREAL8Window(PyObject *self, PyObject *args)
{
	int length;

	if(!PyArg_ParseTuple(args, "i", &length))
		return NULL;

	return pylal_REAL8Window_new(XLALCreateRectangularREAL8Window(length), NULL);
}


static PyObject *pylal_XLALCreateHannREAL8Window(PyObject *self, PyObject *args)
{
	int length;

	if(!PyArg_ParseTuple(args, "i", &length))
		return NULL;

	return pylal_REAL8Window_new(XLALCreateHannREAL8Window(length), NULL);
}


static PyObject *pylal_XLALCreateWelchREAL8Window(PyObject *self, PyObject *args)
{
	int length;

	if(!PyArg_ParseTuple(args, "i", &length))
		return NULL;

	return pylal_REAL8Window_new(XLALCreateWelchREAL8Window(length), NULL);
}


static PyObject *pylal_XLALCreateBartlettREAL8Window(PyObject *self, PyObject *args)
{
	int length;

	if(!PyArg_ParseTuple(args, "i", &length))
		return NULL;

	return pylal_REAL8Window_new(XLALCreateBartlettREAL8Window(length), NULL);
}


static PyObject *pylal_XLALCreateParzenREAL8Window(PyObject *self, PyObject *args)
{
	int length;

	if(!PyArg_ParseTuple(args, "i", &length))
		return NULL;

	return pylal_REAL8Window_new(XLALCreateParzenREAL8Window(length), NULL);
}


static PyObject *pylal_XLALCreatePapoulisREAL8Window(PyObject *self, PyObject *args)
{
	int length;

	if(!PyArg_ParseTuple(args, "i", &length))
		return NULL;

	return pylal_REAL8Window_new(XLALCreatePapoulisREAL8Window(length), NULL);
}


static PyObject *pylal_XLALCreateHammingREAL8Window(PyObject *self, PyObject *args)
{
	int length;

	if(!PyArg_ParseTuple(args, "i", &length))
		return NULL;

	return pylal_REAL8Window_new(XLALCreateHammingREAL8Window(length), NULL);
}


static PyObject *pylal_XLALCreateKaiserREAL8Window(PyObject *self, PyObject *args)
{
	int length;
	double beta;

	if(!PyArg_ParseTuple(args, "id", &length, &beta))
		return NULL;

	return pylal_REAL8Window_new(XLALCreateKaiserREAL8Window(length, beta), NULL);
}


static PyObject *pylal_XLALCreateCreightonREAL8Window(PyObject *self, PyObject *args)
{
	int length;
	double beta;

	if(!PyArg_ParseTuple(args, "id", &length, &beta))
		return NULL;

	return pylal_REAL8Window_new(XLALCreateCreightonREAL8Window(length, beta), NULL);
}


static PyObject *pylal_XLALCreateTukeyREAL8Window(PyObject *self, PyObject *args)
{
	int length;
	double beta;

	if(!PyArg_ParseTuple(args, "id", &length, &beta))
		return NULL;

	return pylal_REAL8Window_new(XLALCreateTukeyREAL8Window(length, beta), NULL);
}


static PyObject *pylal_XLALCreateGaussREAL8Window(PyObject *self, PyObject *args)
{
	int length;
	double beta;

	if(!PyArg_ParseTuple(args, "id", &length, &beta))
		return NULL;

	return pylal_REAL8Window_new(XLALCreateGaussREAL8Window(length, beta), NULL);
}


/*
 * ============================================================================
 *
 *                            Module Registration
 *
 * ============================================================================
 */


static struct PyMethodDef methods[] = {
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
	PyObject *module = Py_InitModule3(MODULE_NAME, methods, "Wrapper for LAL's window package.");

	import_array();

	/* REAL8Window */
	if(PyType_Ready(&pylal_REAL8Window_Type) < 0)
		return;
	Py_INCREF(&pylal_REAL8Window_Type);
	PyModule_AddObject(module, "REAL8Window", (PyObject *) &pylal_REAL8Window_Type);
}
