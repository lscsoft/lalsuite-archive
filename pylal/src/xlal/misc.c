/*
 * $Id$
 *
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
 *                              Helper Utilities
 *
 * ============================================================================
 */


#include <Python.h>
#include <lal/XLALError.h>
#include <misc.h>


/*
 * ============================================================================
 *
 *                               Error Handling
 *
 * ============================================================================
 */


PyObject *pylal_exception_from_errno(enum XLALErrorValue code, const char **msg)
{
	if(msg)
		*msg = XLALErrorString(code);

	switch(code) {
	case XLAL_SUCCESS:
		return NULL;

	case XLAL_EIO:
		return PyExc_IOError;

	case XLAL_ENOMEM:
	case XLAL_EFAULT:
		return PyExc_MemoryError;

	case XLAL_EINVAL:
	case XLAL_EDOM:
	case XLAL_ERANGE:
	case XLAL_EBADLEN:
	case XLAL_ESIZE:
	case XLAL_EDIMS:
		return PyExc_ValueError;

	case XLAL_ETYPE:
		return PyExc_TypeError;

	case XLAL_ETIME:
	case XLAL_EFREQ:
	case XLAL_EUNIT:
	case XLAL_ENAME:
	case XLAL_EDATA:
		return PyExc_ValueError;

	case XLAL_ESYS:
		return PyExc_SystemError;

	case XLAL_EERR:
		return PyExc_RuntimeError;

	case XLAL_EFPDIV0:
		return PyExc_ZeroDivisionError;

	case XLAL_EFPINVAL:
	case XLAL_EFPOVRFLW:
	case XLAL_EFPUNDFLW:
	case XLAL_EFPINEXCT:
	case XLAL_EMAXITER:
	case XLAL_EDIVERGE:
	case XLAL_ESING:
	case XLAL_ETOL:
	case XLAL_ELOSS:
		return PyExc_FloatingPointError;

	default:
		return PyExc_Exception;
	}
}


void pylal_set_exception_from_xlalerrno(void)
{
	const char *msg;

	PyErr_SetString(pylal_exception_from_errno(XLALGetBaseErrno(), &msg), msg);
	XLALClearErrno();
}
