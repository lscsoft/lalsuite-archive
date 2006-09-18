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
	PyObject **types;
	PyObject **types_length;
	PyObject **type;
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
			tokenizer->data = realloc(tokenizer->data, tokenizer->allocation + 1);
			tokenizer->data[tokenizer->allocation] = '\0';
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
		*tokenizer->length = '\0';
	}
}


static void unref_types(ligolw_Tokenizer *tokenizer)
{
	for(tokenizer->type = tokenizer->types; tokenizer->type < tokenizer->types_length; tokenizer->type++)
		Py_DECREF(*tokenizer->type);

	free(tokenizer->types);
	tokenizer->types = NULL;
	tokenizer->types_length = NULL;
	tokenizer->type = NULL;
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

	unref_types(tokenizer);
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
	tokenizer->types = malloc(1 * sizeof(*tokenizer->types));
	tokenizer->types_length = &tokenizer->types[1];
	tokenizer->types[0] = (PyObject *) &PyString_Type;
	Py_INCREF(&PyString_Type);
	tokenizer->type = tokenizer->types;
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


static PyObject *tokenizer_next_string(ligolw_Tokenizer *tokenizer, char **start, char **end)
{
	char *pos = tokenizer->pos;
	PyObject *type = *tokenizer->type;

	/*
	 * The following code matches the pattern:
	 *
	 * any amount of white-space + " + non-quote characters + " + any
	 * amount of white-space + delimiter
	 *
	 * or
	 *
	 * any amount of white-space + non-quote, non-white-space,
	 * non-delimiter characters + any amount of white-space + delimiter
	 *
	 * The middle bit is returned as the token.
	 */

	/*
	 * start == a white-space to non-white-space transition outside of
	 * a quote, or a non-quoted to quoted transition.
	 *
	 * end == a non-white-space to white-space transition outside of a
	 * quote, or a delimiter outside of a quote, or a quoted to
	 * non-quoted transition.
	 */

	if(pos >= tokenizer->length)
		goto stop_iteration;
	while(isspace(*pos))
		if(++pos >= tokenizer->length)
			goto stop_iteration;
	if(*pos == '"') {
		*start = ++pos;
		if(pos >= tokenizer->length)
			goto stop_iteration;
		while(*pos != '"')
			if(++pos >= tokenizer->length)
				goto stop_iteration;
		*end = pos;
		if(++pos >= tokenizer->length)
			goto stop_iteration;
	} else {
		*start = pos;
		while(!isspace(*pos) && (*pos != tokenizer->delimiter) && (*pos != '"'))
			if(++pos >= tokenizer->length)
				goto stop_iteration;
		*end = pos;
	}
	while(*pos != tokenizer->delimiter) {
		if(!isspace(*pos))
			goto parse_error;
		if(++pos >= tokenizer->length)
			goto stop_iteration;
	}

	tokenizer->pos = ++pos;

	if(++tokenizer->type >= tokenizer->types_length)
		tokenizer->type = tokenizer->types;

	return type;

stop_iteration:
	shift(tokenizer, tokenizer->pos);
	PyErr_SetNone(PyExc_StopIteration);
	return NULL;

parse_error:
	PyErr_SetString(PyExc_ValueError, *start);
	return NULL;
}


static PyObject *ligolw_Tokenizer_next(PyObject *self)
{
	PyObject *type;
	PyObject *token;
	char *start, *end;

	do {
		type = tokenizer_next_string((ligolw_Tokenizer *) self, &start, &end);
		if(!type)
			return NULL;
	} while(type == Py_None);

	*end = '\0';
	if(type == (PyObject *) &PyString_Type) {
		token = PyString_FromStringAndSize(start, end - start);
	} else if(type == (PyObject *) &PyInt_Type) {
		token = PyInt_FromString(start, NULL, 0);
	} else if(type == (PyObject *) &PyFloat_Type) {
		double x = strtod(start, &end);
		if(*end == '\0')
			token = PyFloat_FromDouble(x);
		else {
			PyErr_Format(PyExc_ValueError, "invalid literal for float(): %s", start);
			token = NULL;
		}
	} else {
		PyErr_BadArgument();
		token = NULL;
	}

	return token;
}


static PyObject *ligolw_Tokenizer_set_types(PyObject *self, PyObject *list)
{
	ligolw_Tokenizer *tokenizer = (ligolw_Tokenizer *) self;
	int length, i;

	if(!PyList_Check(list))
		goto type_error;
	length = PyList_GET_SIZE(list);
	for(i = 0; i < length; i++) {
		PyObject *type = PyList_GET_ITEM(list, i);
		if((type != (PyObject *) &PyString_Type) && (type != (PyObject *) &PyInt_Type) && (type != (PyObject *) &PyFloat_Type) && (type != Py_None))
			goto type_error;
	}

	unref_types(tokenizer);

	tokenizer->types = malloc(length * sizeof(*tokenizer->types));
	tokenizer->types_length = &tokenizer->types[length];
	tokenizer->type = tokenizer->types;

	for(i = 0; i < length; i++) {
		tokenizer->types[i] = PyList_GET_ITEM(list, i);
		Py_INCREF(tokenizer->types[i]);
	}

	Py_INCREF(Py_None);
	return Py_None;

type_error:
	PyErr_SetString(PyExc_TypeError, "Tokenizer.set_types(): argument must be a list of str, int, or float types, or Nones");
	return NULL;
}


/*
 * Type information
 */

static struct PyMethodDef ligolw_Tokenizer_methods[] = {
	{"add", ligolw_Tokenizer_add, METH_O, "Append a string to the tokenizer's contents"},
	{"set_types", ligolw_Tokenizer_set_types, METH_O, "Set the list of Python types to be used cyclically for token parsing"},
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
