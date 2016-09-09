/*
 * Copyright (C) 2006-2011,2013  Kipp Cannon
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
#include <string.h>
#include <lal/DetectorSite.h>
#include <misc.h>
#include <datatypes/snglringdowntable.h>
#include <lal/LIGOMetadataRingdownUtils.h>
#include <lal/RingUtils.h>


#define MODULE_NAME "pylal.xlal.tools"


/*
 * ============================================================================
 *
 *                                 Coinc Type
 *
 * ============================================================================
 */


/*
 * Cached ID types
 */


static PyObject *process_id_type = NULL;
static PyObject *coinc_event_id_type = NULL;
static PyObject *coinc_def_id_type = NULL;
static PyObject *time_slide_id_type = NULL;


/*
 * Structure
 */


typedef struct {
	PyObject_HEAD
	long process_id_i;
	long coinc_def_id_i;
	long coinc_event_id_i;
	long time_slide_id_i;
	PyObject *instruments;
	int nevents;
	double likelihood;
} ligolw_Coinc;


/*
 * Attributes
 */


static struct PyMemberDef ligolw_Coinc_members[] = {
	{"instruments", T_OBJECT, offsetof(ligolw_Coinc, instruments), 0, "instruments"},
	{"nevents", T_INT, offsetof(ligolw_Coinc, nevents), 0, "nevents"},
	{"likelihood", T_DOUBLE, offsetof(ligolw_Coinc, likelihood), 0, "likelihood"},
	{NULL,}
};


static struct PyGetSetDef ligolw_Coinc_getset[] = {
	{"process_id", pylal_ilwdchar_id_get, pylal_ilwdchar_id_set, "process_id", &(struct pylal_ilwdchar_id_description) {offsetof(ligolw_Coinc, process_id_i), &process_id_type}},
	{"coinc_def_id", pylal_ilwdchar_id_get, pylal_ilwdchar_id_set, "coinc_def_id", &(struct pylal_ilwdchar_id_description) {offsetof(ligolw_Coinc, coinc_def_id_i), &coinc_def_id_type}},
	{"coinc_event_id", pylal_ilwdchar_id_get, pylal_ilwdchar_id_set, "coinc_event_id", &(struct pylal_ilwdchar_id_description) {offsetof(ligolw_Coinc, coinc_event_id_i), &coinc_event_id_type}},
	{"time_slide_id", pylal_ilwdchar_id_get, pylal_ilwdchar_id_set, "time_slide_id", &(struct pylal_ilwdchar_id_description) {offsetof(ligolw_Coinc, time_slide_id_i), &time_slide_id_type}},
	{NULL,}
};


/*
 * Methods
 */


/* FIXME:  add get_instruments(), set_instruments() */


/*
 * Type
 */


PyTypeObject ligolw_Coinc_Type = {
	PyObject_HEAD_INIT(NULL)
	.tp_basicsize = sizeof(ligolw_Coinc),
	.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_CHECKTYPES,
	.tp_name = MODULE_NAME ".Coinc",
	.tp_new = PyType_GenericNew,
	.tp_members = ligolw_Coinc_members,
	.tp_getset = ligolw_Coinc_getset,
};


/*
 * ============================================================================
 *
 *                                 Functions
 *
 * ============================================================================
 */


/*
 * sngl_ringdown related coincidence stuff.
 */


static PyObject *pylal_XLALRingdownTimeError(PyObject *self, PyObject *args)
{
	pylal_SnglRingdownTable *row;
	double ds_sq_threshold;
	double delta_t;

	if(!PyArg_ParseTuple(args, "O!d", &pylal_SnglRingdownTable_Type, &row, &ds_sq_threshold))
		return NULL;

	delta_t = XLALRingdownTimeError(&row->sngl_ringdown, ds_sq_threshold);
	if(XLAL_IS_REAL8_FAIL_NAN(delta_t)) {
		pylal_set_exception_from_xlalerrno();
		return NULL;
	}

	return PyFloat_FromDouble(delta_t);
}


static PyObject *pylal_XLAL3DRinca(PyObject *self, PyObject *args)
{
	pylal_SnglRingdownTable *row1, *row2;
	double result;

	if(!PyArg_ParseTuple(args, "O!O!", &pylal_SnglRingdownTable_Type, &row1, &pylal_SnglRingdownTable_Type, &row2))
		return NULL;

	result = XLAL3DRinca(&row1->sngl_ringdown, &row2->sngl_ringdown);

	if(XLAL_IS_REAL8_FAIL_NAN(result)) {
		XLALClearErrno();
		PyErr_SetString(PyExc_ValueError, "not coincident");
		return NULL;
	}

	return PyFloat_FromDouble(result);
}


/*
 * ============================================================================
 *
 *                            Module Registration
 *
 * ============================================================================
 */


static struct PyMethodDef methods[] = {
	{"XLALRingdownTimeError", pylal_XLALRingdownTimeError, METH_VARARGS, "XLALRingdownTimeError(row, ds^2)\n\nFrom a sngl_ringdown event compute the \\Delta t interval corresponding to the given ds^2 threshold."},
	{"XLAL3DRinca", pylal_XLAL3DRinca, METH_VARARGS, "XLAL3DRinca(row1, row)\n\nTakes two SnglRingdown objects and\ncalculates the distance, ds^2, between them."},
	{NULL,}
};


PyMODINIT_FUNC inittools(void)
{
	PyObject *module = Py_InitModule3(MODULE_NAME, methods, "Wrapper for LAL's tools package.");
	if(!module)
		goto nomodule;

	pylal_snglringdowntable_import();

	/* Cached ID types */
	process_id_type = pylal_get_ilwdchar_class("process", "process_id");
	coinc_def_id_type = pylal_get_ilwdchar_class("coinc_definer", "coinc_def_id");
	coinc_event_id_type = pylal_get_ilwdchar_class("coinc_event", "coinc_event_id");
	time_slide_id_type = pylal_get_ilwdchar_class("time_slide", "time_slide_id");
	if(!process_id_type || !coinc_def_id_type || !coinc_event_id_type || !time_slide_id_type)
		goto noids;

	/* Coinc */
	if(PyType_Ready(&ligolw_Coinc_Type) < 0)
		goto nocoinctype;
	Py_INCREF(&ligolw_Coinc_Type);
	PyModule_AddObject(module, "Coinc", (PyObject *) &ligolw_Coinc_Type);

	return;

nocoinctype:
nocoincmaptype:
noids:
	Py_XDECREF(process_id_type);
	Py_XDECREF(coinc_def_id_type);
	Py_XDECREF(coinc_event_id_type);
	Py_XDECREF(time_slide_id_type);

nomodule:
	return;
}
