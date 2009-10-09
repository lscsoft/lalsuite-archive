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
 *                 Python Wrapper For LAL's REAL8Window Type
 *
 * ============================================================================
 */


#include <Python.h>
#include <numpy/arrayobject.h>
#include <lal/LALDatatypes.h>
#include <lal/Window.h>
#include <real8window.h>


#define MODULE_NAME PYLAL_REAL8WINDOW_MODULE_NAME


/*
 * ============================================================================
 *
 *                              REAL8Window Type
 *
 * ============================================================================
 */


/*
 * Methods
 */


static void __del__(PyObject *self)
{
	pylal_REAL8Window *obj = (pylal_REAL8Window *) self;

	if(obj->owner)
		Py_DECREF(obj->owner);
	else
		/* we are the owner */
		XLALDestroyREAL8Window(obj->window);

	self->ob_type->tp_free(self);
}


static PyObject *__getattr__(PyObject *self, char *name)
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


static PyTypeObject _pylal_REAL8Window_Type = {
	PyObject_HEAD_INIT(NULL)
	.tp_basicsize = sizeof(pylal_REAL8Window),
	.tp_dealloc = __del__,
	.tp_doc = "REAL8Window structure",
	.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_CHECKTYPES,
	.tp_getattr = __getattr__,
	.tp_name = MODULE_NAME ".REAL8Window",
	.tp_new = PyType_GenericNew
};


/*
 * ============================================================================
 *
 *                            Module Registration
 *
 * ============================================================================
 */


void initreal8window(void)
{
	PyObject *module = Py_InitModule3(MODULE_NAME, NULL, "Wrapper for LAL's REAL8Window type.");

	import_array();

	/* REAL8Window */
	pylal_REAL8Window_Type = &_pylal_REAL8Window_Type;
	if(PyType_Ready(pylal_REAL8Window_Type) < 0)
		return;
	Py_INCREF(pylal_REAL8Window_Type);
	PyModule_AddObject(module, "REAL8Window", (PyObject *) pylal_REAL8Window_Type);
}
