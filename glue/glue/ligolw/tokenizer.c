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
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <wchar.h>
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
 * Creates a reference to a new object, does not decref its argument.
 */


PyObject *_build_attributes(PyObject *sequence)
{
	PyObject *result;
	int i;

	/* guaranteed to produce a new object */
	sequence = PySequence_List(sequence);
	if(!sequence)
		return NULL;

	for(i = 0; i < PyList_GET_SIZE(sequence); i++) {
		PyObject *item = PyList_GET_ITEM(sequence, i);
		if(!item) {
			Py_DECREF(sequence);
			return NULL;
		}
		if(!PyString_Check(item)) {
			PyObject *str = PyUnicode_AsEncodedString(item, NULL, "strict");
			if(!str) {
				Py_DECREF(sequence);
				return NULL;
			}
			Py_DECREF(item);
			PyList_SET_ITEM(sequence, i, str);
		}
	}

	result = PySequence_Tuple(sequence);
	Py_DECREF(sequence);

	return result;
}


/*
 * Convert a sequence of functions to a tuple of functions.  Creates a
 * reference to a new object, does not decref its argument.
 */


PyObject *_build_formats(PyObject *sequence)
{
	return PySequence_Tuple(sequence);
}


/*
 * ============================================================================
 *
 *                                   Blobs
 *
 * ============================================================================
 */


static PyObject *parse_blob(PyObject *self, PyObject *args)
{
	PyObject *unicode;
	PyObject *buffer;
	const Py_UNICODE *unicode_data;
	PyObject *string;
	char *pos;
	int i;
	int length;

	/*
	 * Get the input unicode object's data address and length.
	 */

	if(!PyArg_ParseTuple(args, "U", &unicode))
		return NULL;
	unicode_data = PyUnicode_AS_UNICODE(unicode);
	length = PyUnicode_GET_SIZE(unicode);

	/*
	 * Each character turns into (at most) one byte in the buffer.  We
	 * need to arrange for the data to be free()'ed when the buffer
	 * object is garbage collected.  The buffer object has a b_base
	 * element in its structure which points to something to decref
	 * when the buffer is garbage collected, so instead of allocating
	 * our own data we put it inside a string object (hidden from the
	 * rest of the interpreter) and then point the buffer object at the
	 * string object.  The string object's machinery should then take
	 * care of garbage collecting our data.
	 */

	string = PyString_FromStringAndSize(NULL, length);
	if(!string)
		return NULL;

	/*
	 * Translate unicode string.  Whitespace is skipped, a backslash
	 * initates a three-digit octal code, everything else goes in as
	 * itself and must be in the range [0, 128).
	 */

	pos = PyString_AS_STRING(string);
	for(i = 0; i < length; i++) {
		if(Py_UNICODE_ISSPACE(unicode_data[i]))
			continue;
		if(unicode_data[i] != '\\') {
			if(unicode_data[i] < 0 || unicode_data[i] > 127) {
				Py_DECREF(string);
				PyErr_SetObject(PyExc_ValueError, unicode);
				return NULL;
			}
			*(pos++) = unicode_data[i];
		} else if(length - i < 4) {
			Py_DECREF(string);
			PyErr_SetObject(PyExc_ValueError, unicode);
			return NULL;
		} else {
			int h = Py_UNICODE_TODECIMAL(unicode_data[++i]) * 64;
			int m = Py_UNICODE_TODECIMAL(unicode_data[++i]) * 8;
			int l = Py_UNICODE_TODECIMAL(unicode_data[++i]);

			if(h < 0 || m < 0 || l < 0 || h + m + l > 255) {
				Py_DECREF(string);
				PyErr_SetObject(PyExc_ValueError, unicode);
				return NULL;
			}

			*(pos++) = h + m + l;
		}
	}

	/*
	 * Adjust the strings's allocated size.
	 */

	length = pos - PyString_AS_STRING(string);
	if(_PyString_Resize(&string, length))
		return NULL;

	/*
	 * Construct buffer object.  PyBuffer_FromObject() increfs the
	 * string, so when we decref the string the buffer object is the
	 * last place holding a reference to it.
	 */

	buffer = PyBuffer_FromObject(string, 0, length);
	Py_DECREF(string);

	return buffer;
}


static PyObject *blob_format_func(PyObject *self, PyObject *args)
{
	PyObject *unicode;
	wchar_t *string;
	wchar_t *pos;
	const char *buffer;
	int i;
	int length;

	/*
	 * Get the input buffer's data address and length.
	 */

	if(!PyArg_ParseTuple(args, "t#", &buffer, &length))
		return NULL;

	/*
	 * Each char turns into (at most) a 4-character escaped octal code.
	 */

	string = malloc((4 * length + 1) * sizeof(*string));
	if(!string)
		return PyErr_NoMemory();

	/*
	 * Translate buffer.  Spaces, unprintable characters, backslashes,
	 * and things with the high bit set go in as octal codes,
	 * everything else goes in as itself.
	 */

	pos = string;
	for(i = 0; i < length; i++) {
		int count = swprintf(pos, 5, isspace(buffer[i]) || !isprint(buffer[i]) || (buffer[i] == '\\') || (buffer[i] < 0) ? L"\\%03hho" : L"%hhc", buffer[i]);
		if(count < 0) {
			free(string);
			PyErr_SetFromErrno(PyExc_ValueError);
			return NULL;
		}
		pos += count;
	}

	/*
	 * Build into unicode string.
	 */

	unicode = PyUnicode_FromWideChar(string, pos - string);
	free(string);

	return unicode;
}


/*
 * ============================================================================
 *
 *                            Module Registration
 *
 * ============================================================================
 */


static struct PyMethodDef module_methods[] = {
	{"parse_blob", parse_blob, METH_VARARGS, "Converts a unicode string of printable ASCII characters and escaped octal codes into a Python array object of unsigned char data."},
	{"blob_format_func", blob_format_func, METH_VARARGS, "Converts a Python array object of unsigned char data into a unicode string of printable ASCII characters and escaped octal codes."},
	{NULL,}
};


void inittokenizer(void)
{
	/*
	 * Create the module.
	 */

	PyObject *module = Py_InitModule3(MODULE_NAME,
		module_methods,
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
