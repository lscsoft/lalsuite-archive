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
 *                         tokenizer.Tokenizer Class
 *
 * ============================================================================
 */


#include <Python.h>
#include <ctype.h>
#include <stdlib.h>
#include <tokenizer.h>


/*
 * ============================================================================
 *
 *                               Tokenizer Type
 *
 * ============================================================================
 */


/*
 * Structure
 */


typedef struct {
	PyObject_HEAD
	/* list of the types to which parsed tokens will be converted */
	PyObject **types;
	/* end of the types list */
	PyObject **types_length;
	/* the type to which the next parsed token will be converted */
	PyObject **type;
	/* delimiter character to be used in parsing */
	Py_UNICODE delimiter;
	/* size of internal buffer, minus null terminator */
	int allocation;
	/* internal buffer */
	Py_UNICODE *data;
	/* end of internal buffer's contents (null terminator) */
	Py_UNICODE *length;
	/* current offset in buffer */
	Py_UNICODE *pos;
} ligolw_Tokenizer;


/*
 * Append the contents of a unicode object to a tokenizer's internal
 * buffer, increasing the size of the buffer if needed.
 */


static int add_to_data(ligolw_Tokenizer *tokenizer, PyObject *unicode)
{
	int n = PyUnicode_GET_SIZE(unicode);

	if(n) {
		if(tokenizer->length - tokenizer->data + n > tokenizer->allocation) {
			/*
			 * convert pointers to integer offsets
			 */

			int pos = tokenizer->pos - tokenizer->data;
			int length = tokenizer->length - tokenizer->data;

			/*
			 * increase buffer size, adding 1 to leave room for
			 * the null terminator
			 */

			Py_UNICODE *old_data = tokenizer->data;

			tokenizer->data = realloc(tokenizer->data, (tokenizer->allocation + n + 1) * sizeof(*tokenizer->data));
			if(!tokenizer->data) {
				/*
				 * memory failure, restore pointer and exit
				 */

				tokenizer->data = old_data;
				return -1;
			}
			tokenizer->allocation += n;

			/*
			 * convert integer offsets back to pointers
			 */

			tokenizer->pos = &tokenizer->data[pos];
			tokenizer->length = &tokenizer->data[length];
		}

		/*
		 * copy data from unicode into buffer, appending null
		 * terminator
		 */

		memcpy(tokenizer->length, PyUnicode_AsUnicode(unicode), n * sizeof(*tokenizer->length));
		tokenizer->length += n;
		*tokenizer->length = 0;
	}

	/*
	 * success
	 */

	return 0;
}


/*
 * Shift the contents of the tokenizer's buffer so that the data starting
 * at pos is moved to the start of the buffer.  When moving data, add 1 to
 * the length to also move the null terminator.
 */


static void advance_to_pos(ligolw_Tokenizer *tokenizer)
{
	if(tokenizer->pos != tokenizer->data) {
		tokenizer->length -= tokenizer->pos - tokenizer->data;
		memmove(tokenizer->data, tokenizer->pos, (tokenizer->length - tokenizer->data + 1) * sizeof(*tokenizer->data));
		tokenizer->pos = tokenizer->data;
	}
}


/*
 * Free the tokenizer's types list.
 */


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
 * Identify the next token to extract from the tokenizer's internal buffer.
 * On success, start will be left pointing to the address of the start of
 * the string, and end will be pointing to the first character after the
 * string.  The return value is the Python type to which the text should be
 * converted, or NULL on error.  On error, the values of start and end are
 * undefined.  Raises StopIteration if the end of the tokenizer's internal
 * buffer is reached, or ValueError if a parse error occurs.
 */


static PyObject *next_token(ligolw_Tokenizer *tokenizer, Py_UNICODE **start, Py_UNICODE **end)
{
	Py_UNICODE *pos = tokenizer->pos;
	Py_UNICODE *bailout = tokenizer->length;
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

	if(pos >= bailout)
		goto stop_iteration;
	while(Py_UNICODE_ISSPACE(*pos))
		if(++pos >= bailout)
			goto stop_iteration;
	if(*pos == '"') {
		*start = ++pos;
		if(pos >= bailout)
			goto stop_iteration;
		while(*pos != '"')
			if(++pos >= bailout)
				goto stop_iteration;
		*end = pos;
		if(++pos >= bailout)
			goto stop_iteration;
	} else {
		*start = pos;
		while(!Py_UNICODE_ISSPACE(*pos) && (*pos != tokenizer->delimiter) && (*pos != '"'))
			if(++pos >= bailout)
				goto stop_iteration;
		*end = pos;
	}
	while(*pos != tokenizer->delimiter) {
		if(!Py_UNICODE_ISSPACE(*pos))
			goto parse_error;
		if(++pos >= bailout)
			goto stop_iteration;
	}

	tokenizer->pos = ++pos;

	/*
	 * Select the next type
	 */

	if(++tokenizer->type >= tokenizer->types_length)
		tokenizer->type = tokenizer->types;

	/*
	 * Done
	 */

	return type;

	/*
	 * Errors
	 */

stop_iteration:
	advance_to_pos(tokenizer);
	PyErr_SetNone(PyExc_StopIteration);
	return NULL;

parse_error:
	{
	PyObject *unicode = PyUnicode_FromUnicode(*start, tokenizer->length - *start);
	PyErr_SetObject(PyExc_ValueError, unicode);
	Py_DECREF(unicode);
	return NULL;
	}
}


/*
 * append() method
 */


static PyObject *append(PyObject *self, PyObject *data)
{
	int fail;

	if(PyUnicode_Check(data)) {
		fail = add_to_data((ligolw_Tokenizer *) self, data);
	} else if(PyString_Check(data)) {
		/* decode to unicode */
		if(!(data = PyUnicode_FromObject(data)))
			/* decode failure */
			return NULL;
		fail = add_to_data((ligolw_Tokenizer *) self, data);
		Py_DECREF(data);
	} else {
		PyErr_SetObject(PyExc_TypeError, data);
		return NULL;
	}

	if(fail < 0)
		return PyErr_NoMemory();

	Py_INCREF(self);
	return self;
}


/*
 * __del__() method
 */


static void __del__(PyObject *self)
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


/*
 * __init__() method
 */


static int __init__(PyObject *self, PyObject *args, PyObject *kwds)
{
	ligolw_Tokenizer *tokenizer = (ligolw_Tokenizer *) self;
	PyObject *arg;
	int delimiter_length;

	if(!PyArg_ParseTuple(args, "O", &arg))
		return -1;
	if(!(arg = PyUnicode_FromObject(arg)))
		return -1;
	delimiter_length = PyUnicode_GET_SIZE(arg);
	tokenizer->delimiter = *PyUnicode_AS_UNICODE(arg);
	Py_DECREF(arg);
	if(delimiter_length != 1) {
		PyErr_SetString(PyExc_ValueError, "delimiter must have length 1");
		return -1;
	}

	tokenizer->types = malloc(1 * sizeof(*tokenizer->types));
	tokenizer->types_length = &tokenizer->types[1];
	tokenizer->types[0] = (PyObject *) &PyUnicode_Type;
	Py_INCREF(tokenizer->types[0]);
	tokenizer->type = tokenizer->types;
	tokenizer->allocation = 0;
	tokenizer->data = NULL;
	tokenizer->length = tokenizer->data;
	tokenizer->pos = tokenizer->data;

	return 0;
}


/*
 * __iter__() method
 */


static PyObject *__iter__(PyObject *self)
{
	Py_INCREF(self);
	return self;
}


/*
 * next() method
 */


static PyObject *next(PyObject *self)
{
	PyObject *type;
	PyObject *token;
	Py_UNICODE *start, *end;

	/*
	 * Identify the start and end of the next token.
	 */

	do {
		type = next_token((ligolw_Tokenizer *) self, &start, &end);
		if(!type)
			return NULL;
	} while(type == Py_None);

	/*
	 * Null-terminate the token.
	 */

	*end = 0;

	/*
	 * Extract token as desired type.
	 */

	if(type == (PyObject *) &PyFloat_Type) {
		char ascii_buffer[end - start + 1];
		char *ascii_end;
		if(PyUnicode_EncodeDecimal(start, end - start, ascii_buffer, NULL))
			return NULL;
		token = PyFloat_FromDouble(strtod(ascii_buffer, &ascii_end));
		if(*ascii_end != 0) {
			Py_DECREF(token);
			PyErr_Format(PyExc_ValueError, "invalid literal for float(): '%s'", ascii_buffer);
			token = NULL;
		}
	} else if(type == (PyObject *) &PyUnicode_Type) {
		token = PyUnicode_FromUnicode(start, end - start);
	} else if(type == (PyObject *) &PyString_Type) {
		token = PyUnicode_Encode(start, end - start, NULL, NULL);
	} else if(type == (PyObject *) &PyInt_Type) {
		token = PyInt_FromUnicode(start, end - start, 0);
	} else {
		token = PyObject_CallFunction(type, "(u#)", start, end - start);
	}

	/*
	 * Done.
	 */

	return token;
}


/*
 * set_types() method
 */


static PyObject *set_types(PyObject *self, PyObject *list)
{
	ligolw_Tokenizer *tokenizer = (ligolw_Tokenizer *) self;
	int length, i;

	/*
	 * The argument must be a list.
	 */

	if(!PyList_Check(list))
		goto type_error;
	length = PyList_GET_SIZE(list);

	/*
	 * Clear the current internal type list.
	 */

	unref_types(tokenizer);

	/*
	 * Copy the new list's contents into the internal type list.
	 */

	tokenizer->types = malloc(length * sizeof(*tokenizer->types));
	if(!tokenizer->types)
		return PyErr_NoMemory();
	tokenizer->types_length = &tokenizer->types[length];
	tokenizer->type = tokenizer->types;

	for(i = 0; i < length; i++) {
		tokenizer->types[i] = PyList_GET_ITEM(list, i);
		Py_INCREF(tokenizer->types[i]);
	}

	/*
	 * Done.
	 */

	Py_INCREF(Py_None);
	return Py_None;

	/*
	 * Errors.
	 */

type_error:
	PyErr_SetString(PyExc_TypeError, "Tokenizer.set_types(): argument must be a list whose elements are chosen from the types unicode, str, int, and float, or None, or are callable objeccts");
	return NULL;
}


/*
 * Type information
 */


static struct PyMethodDef methods[] = {
	{"append", append, METH_O, "Append a string to the tokenizer's contents."},
	{"set_types", set_types, METH_O, "Set the list of Python types to be used cyclically for token parsing."},
	{NULL,}
};


PyTypeObject ligolw_Tokenizer_Type = {
	PyObject_HEAD_INIT(NULL)
	.tp_basicsize = sizeof(ligolw_Tokenizer),
	.tp_dealloc = __del__,
	.tp_doc =
		"A tokenizer for LIGO Light Weight XML Stream and Array elements.  Converts\n" \
		"(usually comma-) delimited text streams into sequences of Python objects.  An\n" \
		"instance is created by calling the class with the delimiter character as the\n" \
		"single argument.  Text is appended to the internal buffer by passing it to the\n" \
		"append() method.  Tokens are extracted by iterating over the instance.  The\n" \
		"Tokenizer is able to directly extract tokens as various Python types.  The\n" \
		"set_types() method is passed a list of the types to which tokens are to be\n" \
		"converted.  The types will be used in order, cyclically.  For example, passing\n" \
		"[int] to set_types() causes all tokens to be converted to ints, while\n" \
		"[str, int] causes the first token to be returned as a string, the second as an\n" \
		"int, then the third as a string again, and so on.  The default is to extract\n" \
		"all tokens as strings.  Note that the last token will not be extracted until\n" \
		"a delimiter character is seen to terminate it.\n" \
		"\n" \
		"Example:\n" \
		"\n" \
		">>> from glue.ligolw import tokenizer\n" \
		">>> t = tokenizer.Tokenizer(\",\")\n" \
		">>> t.set_types([str, int])\n" \
		">>> list(t.append(\"a,10,b,2\"))\n" \
		"['a', 10, 'b']\n" \
		">>> list(t.append(\"0,\"))\n" \
		"[20]\n",
	.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_CHECKTYPES,
	.tp_init = __init__,
	.tp_iter = __iter__,
	.tp_iternext = next,
	.tp_methods = methods,
	.tp_name = MODULE_NAME ".Tokenizer",
	.tp_new = PyType_GenericNew,
};
