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


#include <Python.h>
#include <structmember.h>
#include <numpy/arrayobject.h>
#include <lal/LALDetectors.h>
#include <laldetector.h>


#define MODULE_NAME PYLAL_LALDETECTOR_MODULE_NAME


/*
 * ============================================================================
 *
 *                                    Type
 *
 * ============================================================================
 */


/*
 * Member access
 */


static struct PyMemberDef members[] = {
	{"name", T_STRING_INPLACE, offsetof(pylal_LALDetector, detector.frDetector.name), READONLY, "name"},
	{"prefix", T_STRING_INPLACE, offsetof(pylal_LALDetector, detector.frDetector.prefix), READONLY, "prefix"},
	{"vertexLongitudeRadians", T_DOUBLE, offsetof(pylal_LALDetector, detector.frDetector.vertexLongitudeRadians), READONLY, "vertexLongitudeRadians"},
	{"vertexLatitudeRadians", T_DOUBLE, offsetof(pylal_LALDetector, detector.frDetector.vertexLatitudeRadians), READONLY, "vertexLatitudeRadians"},
	{"vertexElevation", T_FLOAT, offsetof(pylal_LALDetector, detector.frDetector.vertexElevation), READONLY, "vertexElevation"},
	{"xArmAltitudeRadians", T_FLOAT, offsetof(pylal_LALDetector, detector.frDetector.xArmAltitudeRadians), READONLY, "xArmAltitudeRadians"},
	{"xArmAzimuthRadians", T_FLOAT, offsetof(pylal_LALDetector, detector.frDetector.xArmAzimuthRadians), READONLY, "xArmAzimuthRadians"},
	{"yArmAltitudeRadians", T_FLOAT, offsetof(pylal_LALDetector, detector.frDetector.yArmAltitudeRadians), READONLY, "yArmAltitudeRadians"},
	{"yArmAzimuthRadians", T_FLOAT, offsetof(pylal_LALDetector, detector.frDetector.yArmAzimuthRadians), READONLY, "yArmAzimuthRadians"},
	{"xArmMidpoint", T_FLOAT, offsetof(pylal_LALDetector, detector.frDetector.xArmMidpoint), READONLY, "xArmMidpoint"},
	{"yArmMidpoint", T_FLOAT, offsetof(pylal_LALDetector, detector.frDetector.yArmMidpoint), READONLY, "yArmMidpoint"},
	{"location", T_OBJECT, offsetof(pylal_LALDetector, location), READONLY, "location"},
	{"response", T_OBJECT, offsetof(pylal_LALDetector, response), READONLY, "response"},
	{NULL,}
};


/*
 * Methods
 */


static PyObject *__new__(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
	/* call the generic __new__() */
	pylal_LALDetector *new = (pylal_LALDetector *) PyType_GenericNew(type, args, kwds);

	if(!new)
		return NULL;

	{
	npy_intp dims[] = {3};
	new->location = PyArray_SimpleNewFromData(1, dims, NPY_FLOAT64, new->detector.location);
	}
	{
	npy_intp dims[] = {3, 3};
	new->response = PyArray_SimpleNewFromData(2, dims, NPY_FLOAT32, new->detector.response);
	}

	/* done */
	return (PyObject *) new;
}


static void __del__(PyObject *obj)
{
	pylal_LALDetector *detector = (pylal_LALDetector *) obj;
	Py_DECREF(detector->location);
	detector->location = NULL;
	Py_DECREF(detector->response);
	detector->response = NULL;

	pylal_LALDetector_Type.tp_free(obj);
}


/*
 * Type
 */


PyTypeObject pylal_laldetector_type = {
	PyObject_HEAD_INIT(NULL)
	.tp_basicsize = sizeof(pylal_LALDetector),
	.tp_doc = "LAL's LALDetector structure",
	.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_CHECKTYPES,
	.tp_members = members,
	.tp_name = MODULE_NAME ".LALDetector",
	.tp_new = __new__,
	.tp_dealloc = __del__,
};


/*
 * ============================================================================
 *
 *                            Module Registration
 *
 * ============================================================================
 */


static struct PyMethodDef functions[] = {
	{NULL,}
};


void initlaldetector(void)
{
	PyObject *module = Py_InitModule3(MODULE_NAME, functions, "Wrapper for LAL's LALDetector type.");

	import_array();

	/* LALDetector */
	_pylal_LALDetector_Type = &pylal_laldetector_type;
	if(PyType_Ready(&pylal_LALDetector_Type) < 0)
		return;
	Py_INCREF(&pylal_LALDetector_Type);
	PyModule_AddObject(module, "LALDetector", (PyObject *) &pylal_LALDetector_Type);
}
