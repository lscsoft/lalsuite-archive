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


#include <math.h>


#include <Python.h>
#include <numpy/arrayobject.h>


#include <lal/Date.h>
#include <lal/GenerateBurst.h>
#include <lal/LALSimBurst.h>
#include <lal/LALDatatypes.h>
#include <lal/Sequence.h>
#include <lal/TFTransform.h>
#include <lal/Thresholds.h>
#include <lal/Window.h>


#include <misc.h>
#include <datatypes/complex16frequencyseries.h>
#include <datatypes/real8fftplan.h>
#include <datatypes/real8frequencyseries.h>
#include <datatypes/real8timeseries.h>
#include <datatypes/real8window.h>
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

	return Py_BuildValue("(NN)", hplus_obj, hcross_obj);
}


/*
 * XLALMeasureHrss()
 */


static PyObject *pylal_XLALMeasureHrss(PyObject *self, PyObject *args)
{
	pylal_REAL8TimeSeries *hplus, *hcross;

	if(!PyArg_ParseTuple(args, "O!O!", &pylal_REAL8TimeSeries_Type, &hplus, &pylal_REAL8TimeSeries_Type, &hcross))
		return NULL;

	return PyFloat_FromDouble(XLALMeasureHrss(hplus->series, hcross->series));
}


/*
 * XLALMeasureEoverRsquared()
 */


static PyObject *pylal_XLALMeasureEoverRsquared(PyObject *self, PyObject *args)
{
	pylal_REAL8TimeSeries *hplus, *hcross;

	if(!PyArg_ParseTuple(args, "O!O!", &pylal_REAL8TimeSeries_Type, &hplus, &pylal_REAL8TimeSeries_Type, &hcross))
		return NULL;

	return PyFloat_FromDouble(XLALMeasureEoverRsquared(hplus->series, hcross->series));
}


static PyObject *pylal_XLALEPGetTimingParameters(PyObject *self, PyObject *args)
{
	int window_length;
	int max_tile_length;
	double fractional_tile_stride;
	int psd_length;
	int psd_shift;
	int window_shift;
	int window_pad;
	int tiling_length;

	psd_length = -1;
	if(!PyArg_ParseTuple(args, "iid|i", &window_length, &max_tile_length, &fractional_tile_stride, &psd_length))
		return NULL;

	if(XLALEPGetTimingParameters(window_length, max_tile_length, fractional_tile_stride, psd_length < 0 ? NULL : &psd_length, psd_length < 0 ? NULL : &psd_shift, &window_shift, &window_pad, &tiling_length) < 0) {
	}

	if(psd_length < 0)
		return Py_BuildValue("{s:i,s:i,s:i}", "window_shift", window_shift, "window_pad", window_pad, "tiling_length", tiling_length);
	return Py_BuildValue("{s:i,s:is:i,s:i,s:i}", "psd_length", psd_length, "psd_shift", psd_shift, "window_shift", window_shift, "window_pad", window_pad, "tiling_length", tiling_length);
}


/*
 * pylal_XLALREAL8WindowTwoPointSpectralCorrelation()
 */


static PyObject *pylal_XLALREAL8WindowTwoPointSpectralCorrelation(PyObject *self, PyObject *args)
{
	pylal_REAL8Window *window;
	pylal_REAL8FFTPlan *plan;
	REAL8Sequence *sequence;
	PyObject *result;

	if(!PyArg_ParseTuple(args, "O!O!", &pylal_REAL8Window_Type, &window, &pylal_REAL8FFTPlan_Type, &plan))
		return NULL;

	sequence = XLALREAL8WindowTwoPointSpectralCorrelation(window->window, plan->plan);
	if(!sequence) {
		pylal_set_exception_from_xlalerrno();
		return NULL;
	}

	{
	npy_intp dims[] = {sequence->length};
	result = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, sequence->data);
	}

	/* free sequence without freeing array data */
	if(result)
		sequence->data = NULL;
	XLALDestroyREAL8Sequence(sequence);

	return result;
}


/*
 * pylal_XLALExcessPowerFilterInnerProduct()
 */


static PyObject *pylal_XLALExcessPowerFilterInnerProduct(PyObject *self, PyObject *args)
{
	pylal_COMPLEX16FrequencySeries *filter1;
	pylal_COMPLEX16FrequencySeries *filter2;
	PyArrayObject *correlation;
	REAL8Sequence corr;
	pylal_REAL8FrequencySeries *psd = NULL;
	double result;

	if(!PyArg_ParseTuple(args, "O!O!O!|O", &pylal_COMPLEX16FrequencySeries_Type, &filter1, &pylal_COMPLEX16FrequencySeries_Type, &filter2, &PyArray_Type, &correlation, &psd))
		return NULL;
	correlation = PyArray_GETCONTIGUOUS(correlation);
	if(!correlation)
		return NULL;
	corr.length = PyArray_DIM(correlation, 0);
	corr.data = PyArray_DATA(correlation);
	if((PyObject *) psd == Py_None)
		psd = NULL;
	else if(psd && !PyObject_TypeCheck((PyObject *) psd, &pylal_REAL8FrequencySeries_Type)) {
		Py_DECREF(correlation);
		PyErr_SetObject(PyExc_TypeError, (PyObject *) psd);
		return NULL;
	}

	result = XLALExcessPowerFilterInnerProduct(filter1->series, filter2->series, &corr, psd ? psd->series : NULL);
	if(XLAL_IS_REAL8_FAIL_NAN(result)) {
		Py_DECREF(correlation);
		pylal_set_exception_from_xlalerrno();
		return NULL;
	}

	Py_DECREF(correlation);
	return PyFloat_FromDouble(result);
}


/*
 * pylal_XLALCreateExcessPowerFilter()
 */


static PyObject *pylal_XLALCreateExcessPowerFilter(PyObject *self, PyObject *args)
{
	double flow;
	double width;
	pylal_REAL8FrequencySeries *psd = NULL;
	PyArrayObject *correlation = NULL;
	REAL8Sequence corr;
	COMPLEX16FrequencySeries *filter;

	if(!PyArg_ParseTuple(args, "ddO!O!", &flow, &width, &pylal_REAL8FrequencySeries_Type, &psd, &PyArray_Type, &correlation))
		return NULL;
	correlation = PyArray_GETCONTIGUOUS(correlation);
	if(!correlation)
		return NULL;
	corr.length = PyArray_DIM(correlation, 0);
	corr.data = PyArray_DATA(correlation);

	filter = XLALCreateExcessPowerFilter(flow, width, psd->series, &corr);
	Py_DECREF(correlation);
	if(!filter) {
		pylal_set_exception_from_xlalerrno();
		return NULL;
	}

	return pylal_COMPLEX16FrequencySeries_new(filter, NULL);
}


/*
 * XLALChisqCdf()
 */


static PyObject *pylal_XLALChisqCdf(PyObject *self, PyObject *args)
{
	double chi2;
	double dof;
	double P;

	if(!PyArg_ParseTuple(args, "dd", &chi2, &dof))
		return NULL;

	P = XLALChisqCdf(chi2, dof);
	if(XLALIsREAL8FailNaN(P)) {
		pylal_set_exception_from_xlalerrno();
		return NULL;
	}

	return PyFloat_FromDouble(P);
}


/*
 * XLALOneMinusChisqCdf()
 */


static PyObject *pylal_XLALOneMinusChisqCdf(PyObject *self, PyObject *args)
{
	double chi2;
	double dof;
	double P;

	if(!PyArg_ParseTuple(args, "dd", &chi2, &dof))
		return NULL;

	P = XLALOneMinusChisqCdf(chi2, dof);
	if(XLALIsREAL8FailNaN(P)) {
		pylal_set_exception_from_xlalerrno();
		return NULL;
	}

	return PyFloat_FromDouble(P);
}


/*
 * XLALlnOneMinusChisqCdf()
 */


static PyObject *pylal_XLALlnOneMinusChisqCdf(PyObject *self, PyObject *args)
{
	double chi2;
	double dof;
	double P;

	if(!PyArg_ParseTuple(args, "dd", &chi2, &dof))
		return NULL;

	P = XLALlnOneMinusChisqCdf(chi2, dof);
	if(XLALIsREAL8FailNaN(P)) {
		pylal_set_exception_from_xlalerrno();
		return NULL;
	}

	return PyFloat_FromDouble(P);
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
	{"XLALMeasureHrss", pylal_XLALMeasureHrss, METH_VARARGS, "Measure h_{rss}"},
	{"XLALMeasureEoverRsquared", pylal_XLALMeasureEoverRsquared, METH_VARARGS, "Measure E_{GW}/r^{2}"},
	{"XLALEPGetTimingParameters", pylal_XLALEPGetTimingParameters, METH_VARARGS, NULL},
	{"XLALREAL8WindowTwoPointSpectralCorrelation", pylal_XLALREAL8WindowTwoPointSpectralCorrelation, METH_VARARGS, NULL},
	{"XLALExcessPowerFilterInnerProduct", pylal_XLALExcessPowerFilterInnerProduct, METH_VARARGS, NULL},
	{"XLALCreateExcessPowerFilter", pylal_XLALCreateExcessPowerFilter, METH_VARARGS, NULL},
	{"XLALChisqCdf", pylal_XLALChisqCdf, METH_VARARGS, NULL},
	{"XLALOneMinusChisqCdf", pylal_XLALOneMinusChisqCdf, METH_VARARGS, NULL},
	{"XLALlnOneMinusChisqCdf", pylal_XLALlnOneMinusChisqCdf, METH_VARARGS, NULL},
	{NULL,}
};


PyMODINIT_FUNC initlalburst(void)
{
	Py_InitModule3(MODULE_NAME, methods, "Wrapper for LALBurst package.");

	import_array()

	pylal_complex16frequencyseries_import();
	pylal_real8fftplan_import();
	pylal_real8frequencyseries_import();
	pylal_real8timeseries_import();
	pylal_real8window_import();
	pylal_simburst_import();
}
