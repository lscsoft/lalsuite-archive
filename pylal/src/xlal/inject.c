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
 *           Python Wrapper For LAL's Inject Package (and Friends)
 *
 * ============================================================================
 */

#include <Python.h>
#include <structmember.h>
#include <numpy/arrayobject.h>
#include <lal/LALDetectors.h>
#include <lal/DetResponse.h>


#define MODULE_NAME "pylal.xlal.inject"


/*
 * ============================================================================
 *
 *                              LALDetector Type
 *
 * ============================================================================
 */

/*
 * Forward references
 */

static PyTypeObject pylal_LALDetector_Type;


/*
 * Structure
 */

typedef struct {
	PyObject_HEAD
	LALDetector detector;
	PyObject *location;
	PyObject *response;
} pylal_LALDetector;


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

static PyTypeObject pylal_LALDetector_Type = {
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
 * Amplitude response.
 */

static PyObject *pylal_XLALComputeDetAMResponse(PyObject *self, PyObject *args)
{
	double fplus, fcross;
	double ra;
	double dec;
	double psi;
	double gmst;
	PyObject *response;

	/* 3x3 array, double, double, double, double */
	if(!PyArg_ParseTuple(args, "Odddd:XLALComputeDetAMResponse", &response, &ra, &dec, &psi, &gmst))
		return NULL;
	response = PyArray_FromAny(response, PyArray_DescrFromType(NPY_FLOAT32), 2, 2, NPY_CONTIGUOUS, NULL);
	if(!response || (PyArray_DIM(response, 0) != 3) || (PyArray_DIM(response, 1) != 3)) {
		Py_XDECREF(response);
		return NULL;
	}

	XLALComputeDetAMResponse(&fplus, &fcross, PyArray_DATA(response), ra, dec, psi, gmst);
	Py_XDECREF(response);

	/* (double, double) */
	return Py_BuildValue("(dd)", fplus, fcross);
}


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
	{"XLALComputeDetAMResponse", pylal_XLALComputeDetAMResponse, METH_VARARGS, NULL},
	{NULL,}
};


void initinject(void)
{
	PyObject *module = Py_InitModule3(MODULE_NAME, methods, "Wrapper for LAL's inject package (and friends).");

	import_array();

	/* LALDetector */
	if(PyType_Ready(&pylal_LALDetector_Type) < 0)
		return;
	Py_INCREF(&pylal_LALDetector_Type);
	PyModule_AddObject(module, "LALDetector", (PyObject *) &pylal_LALDetector_Type);
	PyModule_AddObject(module, "cached_detector", make_cached_detectors());
}
