/*
 * Copyright (C) 2015
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


#ifndef _PYLAL_XLAL_DATATYPES_POSTCOHINSPIRALTABLE_H_
#define _PYLAL_XLAL_DATATYPES_POSTCOHINSPIRALTABLE_H_


#include <Python.h>
#include <gstlal-ugly/postcohinspiral_table.h>
#include <lal/LIGOMetadataTables.h>


#define PYLAL_POSTCOHINSPIRALTABLE_MODULE_NAME "pylal.xlal.datatypes.postcohinspiraltable"


/*
 * ============================================================================
 *
 *                                    Type
 *
 * ============================================================================
 */


static PyTypeObject *_pylal_PostcohInspiralTable_Type = NULL;
#define pylal_PostcohInspiralTable_Type (*_pylal_PostcohInspiralTable_Type)


typedef struct {
	PyObject_HEAD
	PostcohInspiralTable postcoh_inspiral;
	/* FIXME:  these should be incorporated into the LAL structure */
	long process_id_i;
	long event_id_i;
} pylal_PostcohInspiralTable;


static PyObject *pylal_postcohinspiraltable_import(void)
{
	PyObject *name = PyString_FromString(PYLAL_POSTCOHINSPIRALTABLE_MODULE_NAME);
	PyObject *module = PyImport_Import(name);
	Py_DECREF(name);

	name = PyString_FromString("PostcohInspiralTable");
	_pylal_PostcohInspiralTable_Type = (PyTypeObject *) PyDict_GetItem(PyModule_GetDict(module), name);
	Py_INCREF(&pylal_PostcohInspiralTable_Type);
	Py_DECREF(name);

	return module;
}


static PyObject *pylal_PostcohInspiralTable_new(const PostcohInspiralTable *row)
{
	PyObject *empty_tuple = PyTuple_New(0);
	pylal_PostcohInspiralTable *obj = (pylal_PostcohInspiralTable *) PyType_GenericNew(&pylal_PostcohInspiralTable_Type, empty_tuple, NULL);
	Py_DECREF(empty_tuple);
	if(!obj)
		return NULL;

	obj->postcoh_inspiral = *row;
	obj->event_id_i = row->event_id;
	obj->process_id_i = 0;

	return (PyObject *) obj;
}


#endif
