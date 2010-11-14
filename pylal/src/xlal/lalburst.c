/*
 * Copyright (C) 2010  Kipp Cannon
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
 *                        Python Wrapper For LALBurst
 *
 * ============================================================================
 */


#include <Python.h>
#include <lal/GenerateBurst.h>
#include <misc.h>
#include <datatypes/real8timeseries.h>
#include <datatypes/simburst.h>


#define MODULE_NAME "pylal.xlal.lalburst"


/*
 * ============================================================================
 *
 *                                 Functions
 *
 * ============================================================================
 */


/*
 * XLALGenerateSimBurst()
 */


static PyObject *pylal_XLALGenerateSimBurst(PyObject *self, PyObject *args)
{
	PyObject *hplus_obj, *hcross_obj;
	REAL8TimeSeries *hplus, *hcross;
	pylal_SimBurst *sim_burst;
	double delta_t;

	if(!PyArg_ParseTuple(args, "O!d", &pylal_SimBurst_Type, &sim_burst, &delta_t))
		return NULL;

	if(XLALGenerateSimBurst(&hplus, &hcross, &sim_burst->sim_burst, delta_t)) {
		pylal_set_exception_from_xlalerrno();
		return NULL;
	}

	hplus_obj = pylal_REAL8TimeSeries_new(hplus, NULL);
	hcross_obj = pylal_REAL8TimeSeries_new(hcross, NULL);
	if(!hplus_obj || !hcross_obj) {
		Py_XDECREF(hplus_obj);
		Py_XDECREF(hcross_obj);
		return NULL;
	}

	return Py_BuildValue("(OO)", hplus_obj, hcross_obj);
}


/*
 * ============================================================================
 *
 *                            Module Registration
 *
 * ============================================================================
 */


static struct PyMethodDef methods[] = {
	{"XLALGenerateSimBurst", pylal_XLALGenerateSimBurst, METH_VARARGS, "Compute the h+ and hx time series for a row in a LIGO Light Weight XML sim_burst table."},
	{NULL,}
};


void initlalburst(void)
{
	Py_InitModule3(MODULE_NAME, methods, "Wrapper for LALBurst package.");

	pylal_real8timeseries_import();
	pylal_simburst_import();
}
