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
 *               Python Wrapper For LAL's Noise Models Package
 *
 * ============================================================================
 */


#include <Python.h>
#include <lal/LALDatatypes.h>
#include <lal/LALNoiseModels.h>


#define MODULE_NAME "pylal.xlal.noisemodels"


/*
 * ============================================================================
 *
 *                                 Functions
 *
 * ============================================================================
 */


static PyObject *pylal_XLALLIGOIPsd(PyObject *self, PyObject *args)
{
	double f;

	if(!PyArg_ParseTuple(args, "d", &f))
		return NULL;

	return PyFloat_FromDouble(XLALLIGOIPsd(f));
}


/*
 * ============================================================================
 *
 *                            Module Registration
 *
 * ============================================================================
 */


static struct PyMethodDef methods[] = {
	{"XLALLIGOIPsd", pylal_XLALLIGOIPsd, METH_VARARGS, NULL},
	{NULL,}
};


void initnoisemodels(void)
{
	/* commented out to silence warning;  uncomment when needed again */
	/*PyObject *module =*/ Py_InitModule3(MODULE_NAME, methods, "Wrapper for LAL's noisemodels package.");
}
