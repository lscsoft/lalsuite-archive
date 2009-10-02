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
 * Use this file like this:
 */

#if 0
#define PYLAL_REAL8FREQUENCYSERIES_MODULE_NAME "pylal.xlal.datatypes.real8frequencyseries"

#include <frequencyseries.h>

MAKE_FREQUENCYSERIES_HEADER(PYLAL_REAL8FREQUENCYSERIES_MODULE_NAME, real8frequencyseries, REAL8FrequencySeries)

#undef MAKE_FREQUENCYSERIES_HEADER

#endif


/*
 * ============================================================================
 *
 *                                  Preamble
 *
 * ============================================================================
 */


#include <Python.h>
#include <lal/LALDatatypes.h>


#define MAKE_FREQUENCYSERIES_HEADER(__MODULE_NAME__, __PYLAL_NAME__, __LAL_NAME__) \
\
\
/*\
 * ============================================================================\
 *\
 *                                    Type\
 *\
 * ============================================================================\
 */\
\
\
static PyTypeObject *pylal_##__LAL_NAME__##_Type;\
\
\
typedef struct {\
	PyObject_HEAD\
	PyObject *owner;\
	__LAL_NAME__ *series;\
} pylal_##__LAL_NAME__;\
\
\
static PyObject *pylal_##__PYLAL_NAME__##_import(void)\
{\
	PyObject *name = PyString_FromString(__MODULE_NAME__);\
	PyObject *module = PyImport_Import(name);\
	Py_DECREF(name);\
\
	name = PyString_FromString("__LAL_NAME__");\
	pylal_##__LAL_NAME__##_Type = (PyTypeObject *) PyDict_GetItem(PyModule_GetDict(module), name);\
	Py_INCREF(pylal_##__LAL_NAME__##_Type);\
	Py_DECREF(name);\
\
	return module;\
}\
\
\
PyObject *pylal_##__LAL_NAME__##_new(__LAL_NAME__ *series, PyObject *owner)\
{\
	PyObject *empty_tuple = PyTuple_New(0);\
	pylal_##__LAL_NAME__ *obj = (pylal_##__LAL_NAME__ *) PyType_GenericNew(pylal_##__LAL_NAME__##_Type, empty_tuple, NULL);\
	Py_DECREF(empty_tuple);\
	if(owner)\
		Py_INCREF(owner);\
	obj->owner = owner;\
	obj->series = series;\
	return (PyObject *) obj;\
}
