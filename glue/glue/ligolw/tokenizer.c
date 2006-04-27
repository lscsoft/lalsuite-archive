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
 *              LIGO Light Weight XML Stream and Array Tokenizer
 *
 * ============================================================================
 */

#include <Python.h>
#include <structmember.h>
#include <ctype.h>
#include <stdlib.h>

#define MODULE_NAME "glue.ligolw.tokenizer"


/*
 * ============================================================================
 *
 *                               Tokenizer Type
 *
 * ============================================================================
 */

/*
 * Forward references
 */

static PyTypeObject ligolw_Tokenizer_Type;


/*
 * Structure
 */

typedef struct {
	PyObject_HEAD
	char delimiter;
	int allocation;
	char *data;
	char *length;
	char *pos;
} ligolw_Tokenizer;


/*
 * Member access
 */

static struct PyMemberDef ligolw_Tokenizer_members[] = {
	{"data", T_STRING, offsetof(ligolw_Tokenizer, data), 0, "text remaining after last iteration"},
	{NULL,}
};


/*
 * Utilities
 */


static void add_to_data(ligolw_Tokenizer *tokenizer, PyObject *string)
{
	int n = PyString_GET_SIZE(string);

	if(n) {
		if(tokenizer->length - tokenizer->data + n > tokenizer->allocation) {
			int pos = tokenizer->pos - tokenizer->data;
			int length = tokenizer->length - tokenizer->data;
			tokenizer->allocation += n;
			tokenizer->data = realloc(tokenizer->data, tokenizer->allocation);
			tokenizer->pos = &tokenizer->data[pos];
			tokenizer->length = &tokenizer->data[length];
		}
		memcpy(tokenizer->length, PyString_AS_STRING(string), n);
		tokenizer->length += n;
	}
}


static void shift(ligolw_Tokenizer *tokenizer, char *start)
{
	int n = start - tokenizer->data;

	if(n) {
		tokenizer->length -= n;
		tokenizer->pos -= n;
		memmove(tokenizer->data, start, tokenizer->length - tokenizer->data);
	}
}


/*
 * Methods
 */

static PyObject *ligolw_Tokenizer_add(PyObject *self, PyObject *data)
{
	if(PyString_Check(data)) {
		add_to_data((ligolw_Tokenizer *) self, data);
	} else if(PyUnicode_Check(data)) {
		PyObject *string = PyUnicode_AsASCIIString(data);
		add_to_data((ligolw_Tokenizer *) self, string);
		Py_DECREF(string);
	} else {
		PyErr_SetString(PyExc_TypeError, "Tokenzer.add(): argument must be a string or unicode string");
		return NULL;
	}

	Py_INCREF(self);
	return self;
}


static void ligolw_Tokenizer___del__(PyObject *self)
{
	ligolw_Tokenizer *tokenizer = (ligolw_Tokenizer *) self;

	free(tokenizer->data);
	tokenizer->data = NULL;
	tokenizer->allocation = 0;
	tokenizer->length = NULL;
	tokenizer->pos = NULL;

	self->ob_type->tp_free(self);
}


static int ligolw_Tokenizer___init__(PyObject *self, PyObject *args, PyObject *kwds)
{
	ligolw_Tokenizer *tokenizer = (ligolw_Tokenizer *) self;
	char *delimiter;
	int delimiter_length;

	if(!PyArg_ParseTuple(args, "s#", &delimiter, &delimiter_length))
		return -1;
	if(delimiter_length != 1) {
		PyErr_SetString(PyExc_TypeError, "argument must have length 1");
		return -1;
	}

	tokenizer->delimiter = *delimiter;
	tokenizer->allocation = 0;
	tokenizer->data = NULL;
	tokenizer->length = tokenizer->data;
	tokenizer->pos = tokenizer->data;

	return 0;
}


static PyObject *ligolw_Tokenizer___iter__(PyObject *self)
{
	Py_INCREF(self);
	return self;
}


static PyObject *ligolw_Tokenizer_next(PyObject *self)
{
	ligolw_Tokenizer *tokenizer = (ligolw_Tokenizer *) self;
	char *pos = tokenizer->pos;
	char *start, *end;

	/*
	 * a token is:
	 *
	 * any amount of white-space + " + non-quote characters + " + any
	 * amount of white-space + delimiter
	 *
	 * or
	 *
	 * any amount of white-space + non-quote, non-white-space,
	 * non-delimiter characters + any amount of white-space + delimiter
	 */

	if(pos >= tokenizer->length)
		goto stop_iteration;
	while(isspace(*pos))
		if(++pos >= tokenizer->length)
			goto stop_iteration;
	if(*pos == '"') {
		start = ++pos;
		if(pos >= tokenizer->length)
			goto stop_iteration;
		while(*pos != '"')
			if(++pos >= tokenizer->length)
				goto stop_iteration;
		end = pos;
		if(++pos >= tokenizer->length)
			goto stop_iteration;
	} else {
		start = pos;
		while(!isspace(*pos) && (*pos != tokenizer->delimiter) && (*pos != '"'))
			if(++pos >= tokenizer->length)
				goto stop_iteration;
		end = pos;
	}
	while(isspace(*pos))
		if(++pos >= tokenizer->length)
			goto stop_iteration;
	if(*pos != tokenizer->delimiter)
		goto parse_error;
	tokenizer->pos = ++pos;

	return PyString_FromStringAndSize(start, end - start);

stop_iteration:
	shift(tokenizer, tokenizer->pos);
	PyErr_SetNone(PyExc_StopIteration);
	return NULL;

parse_error:
	{
	char msg[tokenizer->length - start + 1];
	snprintf(msg, tokenizer->length - start + 1, "%s", start);
	PyErr_SetString(PyExc_ValueError, msg);
	}
	return NULL;
}


/*
 * Type information
 */

static struct PyMethodDef ligolw_Tokenizer_methods[] = {
	{"add", ligolw_Tokenizer_add, METH_O, "Append a string to the tokenizer's contents"},
	{NULL,}
};

static PyTypeObject ligolw_Tokenizer_Type = {
	PyObject_HEAD_INIT(NULL)
	.tp_basicsize = sizeof(ligolw_Tokenizer),
	.tp_dealloc = ligolw_Tokenizer___del__,
	.tp_doc = "A tokenizer for LIGO Light Weight XML Stream and Array elements",
	.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_CHECKTYPES,
	.tp_init = ligolw_Tokenizer___init__,
	.tp_iter = ligolw_Tokenizer___iter__,
	.tp_iternext = ligolw_Tokenizer_next,
	.tp_members = ligolw_Tokenizer_members,
	.tp_methods = ligolw_Tokenizer_methods,
	.tp_name = MODULE_NAME ".Tokenizer",
	.tp_new = PyType_GenericNew,
};


/*
 * ============================================================================
 *
 *                            Module Registration
 *
 * ============================================================================
 */


void inittokenizer(void)
{
	PyObject *module = Py_InitModule3(MODULE_NAME, NULL, "This module provides a tokenizer for LIGO Light Weight XML Stream and Array elements");

	if(PyType_Ready(&ligolw_Tokenizer_Type) < 0)
		return;
	Py_INCREF(&ligolw_Tokenizer_Type);
	PyModule_AddObject(module, "Tokenizer", (PyObject *) &ligolw_Tokenizer_Type);
}
