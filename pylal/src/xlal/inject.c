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
#include <numarray/libnumarray.h>
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
	{NULL,}
};


/*
 * Methods
 */

static PyObject *pylal_LALDetector_location(PyObject *self, PyObject *args)
{
	double *location = ((pylal_LALDetector *) self)->detector.location;

	return Py_BuildValue("(ddd)", location[0], location[1], location[2]);
}


static PyObject *pylal_LALDetector_response(PyObject *self, PyObject *args)
{
	float *response[3] = {((pylal_LALDetector *) self)->detector.response[0], ((pylal_LALDetector *) self)->detector.response[1], ((pylal_LALDetector *) self)->detector.response[2]};

	return Py_BuildValue("((ddd)(ddd)(ddd))", (double) response[0][0], (double) response[0][1], (double) response[0][2], (double) response[1][0], (double) response[1][1], (double) response[1][2], (double) response[2][0], (double) response[2][1], (double) response[2][2]);
}


static struct PyMethodDef pylal_LALDetector_methods[] = {
	{"location", pylal_LALDetector_location, METH_VARARGS, "get detector location"},
	{"response", pylal_LALDetector_response, METH_VARARGS, "get detector response"},
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
	.tp_methods = pylal_LALDetector_methods,
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
	PyObject *oresponse;
	PyArrayObject *response;

	/* 3x3 array, double, double, double, double */
	if(!PyArg_ParseTuple(args, "Odddd:XLALComputeDetAMResponse", &oresponse, &ra, &dec, &psi, &gmst))
		return NULL;
	response = NA_InputArray(oresponse, tFloat32, NUM_C_ARRAY);
	if(!response || (response->nd != 2) || (response->dimensions[0] != 3 && response->dimensions[1] != 3)) {
		Py_XDECREF(response);
		return NULL;
	}

	XLALComputeDetAMResponse(&fplus, &fcross, NA_OFFSETDATA(response), ra, dec, psi, gmst);
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

	import_libnumarray();

	/* LALDetector */
	if(PyType_Ready(&pylal_LALDetector_Type) < 0)
		return;
	Py_INCREF(&pylal_LALDetector_Type);
	PyModule_AddObject(module, "LALDetector", (PyObject *) &pylal_LALDetector_Type);
	PyModule_AddObject(module, "cached_detector", make_cached_detectors());
}
