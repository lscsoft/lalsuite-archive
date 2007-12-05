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
 *                   glue.ligolw.tokenizer Extension Module
 *
 * ============================================================================
 */


#include <Python.h>
#include <stdlib.h>
#include <tokenizer.h>


/*
 * ============================================================================
 *
 *                              Helper Functions
 *
 * ============================================================================
 */


/*
 * Convert a sequence of unicode and/or strings to a tuple of strings.
 */


PyObject *_build_attributes(PyObject *sequence)
{
	PyObject *attributes;

	sequence = PySequence_Tuple(sequence);
	if(!sequence)
		return NULL;

	attributes = PyTuple_New(PyTuple_GET_SIZE(sequence));
	if(attributes) {
		int i;
		for(i = 0; i < PyTuple_GET_SIZE(sequence); i++) {
			PyObject *item = PyTuple_GET_ITEM(sequence, i);
			if(item) {
				PyObject *str;
				if(PyString_Check(item)) {
					str = item;
					Py_INCREF(str);
				} else
					str = PyUnicode_AsEncodedString(item, NULL, "strict");
				if(str) {
					PyTuple_SET_ITEM(attributes, i, str);
					continue;
				}
			}
			Py_DECREF(attributes);
			attributes = NULL;
			break;
		}
	}

	Py_DECREF(sequence);

	return attributes;
}


/*
 * Convert a sequence of unicode and/or strings to a tuple of unicodes.
 */


PyObject *_build_formats(PyObject *sequence)
{
	PyObject *formats;

	sequence = PySequence_Tuple(sequence);
	if(!sequence)
		return NULL;

	formats = PyTuple_New(PyTuple_GET_SIZE(sequence));
	if(formats) {
		int i;
		for(i = 0; i < PyTuple_GET_SIZE(sequence); i++) {
			PyObject *item = PyTuple_GET_ITEM(sequence, i);
			if(item) {
				PyObject *unicd;
				if(PyUnicode_Check(item)) {
					unicd = item;
					Py_INCREF(unicd);
				} else
					unicd = PyUnicode_FromObject(item);
				if(unicd) {
					PyTuple_SET_ITEM(formats, i, unicd);
					continue;
				}
			}
			Py_DECREF(formats);
			formats = NULL;
			break;
		}
	}

	Py_DECREF(sequence);

	return formats;
}


/*
 * ============================================================================
 *
 *                            Module Registration
 *
 * ============================================================================
 */


void inittokenizer(void)
{
	/*
	 * Create the module.
	 */

	PyObject *module = Py_InitModule3(MODULE_NAME, NULL,
		"This module provides a tokenizer for LIGO Light Weight XML Stream and Array\n" \
		"elements, as well as other utilities to assist in packing parsed tokens into\n" \
		"various data storage units."
	);

	/*
	 * Add the Tokenizer class.
	 */

	if(PyType_Ready(&ligolw_Tokenizer_Type) < 0)
		return;
	Py_INCREF(&ligolw_Tokenizer_Type);
	PyModule_AddObject(module, "Tokenizer", (PyObject *) &ligolw_Tokenizer_Type);

	/*
	 * Add the RowBuilder class.
	 */

	if(PyType_Ready(&ligolw_RowBuilder_Type) < 0)
		return;
	Py_INCREF(&ligolw_RowBuilder_Type);
	PyModule_AddObject(module, "RowBuilder", (PyObject *) &ligolw_RowBuilder_Type);

	/*
	 * Add the RowDumper class.
	 */

	if(PyType_Ready(&ligolw_RowDumper_Type) < 0)
		return;
	Py_INCREF(&ligolw_RowDumper_Type);
	PyModule_AddObject(module, "RowDumper", (PyObject *) &ligolw_RowDumper_Type);
}
