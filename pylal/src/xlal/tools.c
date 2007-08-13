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
 *                   Python Wrapper For LAL's Tools Package
 *
 * ============================================================================
 */


#include <Python.h>
#include <structmember.h>
#include <numpy/arrayobject.h>
#include <lal/DetectorSite.h>
#include <tools.h>


#define MODULE_NAME "pylal.xlal.tools"


/*
 * ============================================================================
 *
 *                              LALDetector Type
 *
 * ============================================================================
 */


/*
 * Member access
 */


static struct PyMemberDef pylal_LALDetector_members[] = {
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
 * Type
 */


PyTypeObject pylal_LALDetector_Type = {
	PyObject_HEAD_INIT(NULL)
	.tp_basicsize = sizeof(pylal_LALDetector),
	.tp_doc = "LALDetector structure",
	.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_CHECKTYPES,
	.tp_members = pylal_LALDetector_members,
	.tp_name = MODULE_NAME ".LALDetector",
	.tp_new = PyType_GenericNew
};


/*
 * ============================================================================
 *
 *                                 Functions
 *
 * ============================================================================
 */


/*
 * ============================================================================
 *
 *                            Module Registration
 *
 * ============================================================================
 */


static PyObject *make_cached_detectors(void)
{
	PyObject *cached_detector = PyDict_New();
	pylal_LALDetector *new;
	int i;

	for(i = 0; i < LALNumCachedDetectors; i++) {
		new = (pylal_LALDetector *) _PyObject_New(&pylal_LALDetector_Type);
		memcpy(&new->detector, &lalCachedDetectors[i], sizeof(new->detector));
		{
		npy_intp dims[] = {3};
		new->location = PyArray_SimpleNewFromData(1, dims, NPY_FLOAT64, new->detector.location);
		}
		{
		npy_intp dims[] = {3, 3};
		new->response = PyArray_SimpleNewFromData(2, dims, NPY_FLOAT32, new->detector.response);
		}

		PyDict_SetItemString(cached_detector, new->detector.frDetector.name, (PyObject *) new);
	}

	return cached_detector;
}


static struct PyMethodDef methods[] = {
	{NULL,}
};


void inittools(void)
{
	PyObject *module = Py_InitModule3(MODULE_NAME, methods, "Wrapper for LAL's tools package.");

	import_array();

	/* LALDetector */
	if(PyType_Ready(&pylal_LALDetector_Type) < 0)
		return;
	Py_INCREF(&pylal_LALDetector_Type);
	PyModule_AddObject(module, "LALDetector", (PyObject *) &pylal_LALDetector_Type);
	PyModule_AddObject(module, "cached_detector", make_cached_detectors());
}
