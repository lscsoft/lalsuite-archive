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
 *                   Python Wrapper For LAL's Tools Package
 *
 * ============================================================================
 */


#ifndef _PYLAL_XLAL_TOOLS_H_
#define _PYLAL_XLAL_TOOLS_H_


#include <Python.h>
#include <lal/LALDetectors.h>


/*
 * ============================================================================
 *
 *                              LALDetector Type
 *
 * ============================================================================
 */


/*
 * Type
 */
static PyTypeObject *_pylal_LALDetector_Type = NULL;
#define pylal_LALDetector_Type (*_pylal_LALDetector_Type)
#define PYLAL_LALDETECTOR_MODULE_NAME "pylal.xlal.tools"
/*
 * Structure
 */


typedef struct {
    PyObject_HEAD
    PyObject* owner;
    LALDetector detector;
    PyObject *location;
    PyObject *response;
} pylal_LALDetector;

static PyObject *pylal_laldetector_import(void)
{
    PyObject *name = PyString_FromString(PYLAL_LALDETECTOR_MODULE_NAME);
    PyObject *module = PyImport_Import(name);
    Py_DECREF(name);

    name = PyString_FromString("LALDetector");
    _pylal_LALDetector_Type = (PyTypeObject *) PyDict_GetItem(PyModule_GetDict(module), name);
    Py_INCREF(&pylal_LALDetector_Type);
    Py_DECREF(name);

    return module;
}


static PyObject *pylal_LALDetector_new(LALDetector detector, PyObject *owner)
{
    PyObject *empty_tuple = PyTuple_New(0);
    pylal_LALDetector *obj = (pylal_LALDetector *) PyType_GenericNew(&pylal_LALDetector_Type, empty_tuple, NULL);
    Py_DECREF(empty_tuple);
    if(!obj) {
        return NULL;
    }
    if(owner)
        Py_INCREF(owner);
    obj->owner = owner;
    obj->detector = detector;
    return (PyObject *) obj;
}

#endif /* _PYLAL_XLAL_TOOLS_H_ */
