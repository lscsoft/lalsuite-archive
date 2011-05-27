/*
 * Copyright (C) 2006-2011  Kipp Cannon
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
 *                                  Preamble
 *
 * ============================================================================
 */


#ifndef _PYLAL_XLAL_DATATYPES_LALDETECTOR_H_
#define _PYLAL_XLAL_DATATYPES_LALDETECTOR_H_


#include <Python.h>
#include <numpy/arrayobject.h>
#include <lal/LALDetectors.h>


#define PYLAL_LALDETECTOR_MODULE_NAME "pylal.xlal.datatypes.laldetector"


/*
 * ============================================================================
 *
 *                                    Type
 *
 * ============================================================================
 */


static PyTypeObject *_pylal_LALDetector_Type = NULL;
#define pylal_LALDetector_Type (*_pylal_LALDetector_Type)


typedef struct {
	PyObject_HEAD
	LALDetector detector;
	PyObject *location;
	PyObject *response;
} pylal_LALDetector;


static PyObject *pylal_laldetector_import(void)
{
	PyObject *name = PyString_FromString(PYLAL_LALDETECTOR_MODULE_NAME);
	PyObject *module = PyImport_Import(name);
	Py_DECREF(name);

	if(_import_array() < 0) {
		Py_DECREF(module);
		return NULL;
	}

	name = PyString_FromString("LALDetector");
	_pylal_LALDetector_Type = (PyTypeObject *) PyDict_GetItem(PyModule_GetDict(module), name);
	Py_INCREF(&pylal_LALDetector_Type);
	Py_DECREF(name);

	return module;
}


static PyObject *pylal_LALDetector_new(const LALDetector *detector)
{
	PyObject *empty_tuple = PyTuple_New(0);
	pylal_LALDetector *obj = (pylal_LALDetector *) PyType_GenericNew(&pylal_LALDetector_Type, empty_tuple, NULL);
	Py_DECREF(empty_tuple);
	if(!obj)
		return NULL;

	obj->detector = *detector;
	{
	npy_intp dims[] = {3};
	obj->location = PyArray_SimpleNewFromData(1, dims, NPY_FLOAT64, obj->detector.location);
	}
	{
	npy_intp dims[] = {3, 3};
	obj->response = PyArray_SimpleNewFromData(2, dims, NPY_FLOAT32, obj->detector.response);
	}

	return (PyObject *) obj;
}


#endif /* _PYLAL_XLAL_DATATYPES_LALDETECTOR_H_ */
