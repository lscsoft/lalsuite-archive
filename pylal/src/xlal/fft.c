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
 * ============================================================================
 *
 *                    Python Wrapper For LAL's FFT Package
 *
 * ============================================================================
 */


#include <Python.h>
#include <lal/TimeFreqFFT.h>
#include <misc.h>
#include <datatypes/complex16frequencyseries.h>
#include <datatypes/real8frequencyseries.h>


#define MODULE_NAME "pylal.xlal.fft"


/*
 * ============================================================================
 *
 *                             Function Wrappers
 *
 * ============================================================================
 */


static PyObject *pylal_XLALWhitenCOMPLEX16FrequencySeries(PyObject *self, PyObject *args)
{
	pylal_COMPLEX16FrequencySeries *fseries;
	pylal_REAL8FrequencySeries *psd;

	if(!PyArg_ParseTuple(args, "O!O!:XLALWhitenCOMPLEX16FrequencySeries", &pylal_COMPLEX16FrequencySeries_Type, &fseries, &pylal_REAL8FrequencySeries_Type, &psd))
		return NULL;

	if(!XLALWhitenCOMPLEX16FrequencySeries(fseries->series, psd->series)) {
		pylal_set_exception_from_xlalerrno();
		XLALClearErrno();
		return NULL;
	}

	Py_INCREF(Py_None);
	return Py_None;
}


/*
 * ============================================================================
 *
 *                            Module Registration
 *
 * ============================================================================
 */


static struct PyMethodDef module_methods[] = {
	{"XLALWhitenCOMPLEX16FrequencySeries", pylal_XLALWhitenCOMPLEX16FrequencySeries, METH_VARARGS, NULL},
	{NULL,}
};


void initfft(void)
{
	/* commented out to silence warning */
	/*PyObject *module = */Py_InitModule3(MODULE_NAME, module_methods, "Wrapper for LAL's fft package.");

	pylal_complex16frequencyseries_import();
	pylal_real8frequencyseries_import();
}
